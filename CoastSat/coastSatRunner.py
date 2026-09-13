#!/usr/bin/env python

import argparse
import pickle
import pandas as pd
import pytz
import numpy as np
import os
import sys
import ast
import logging
from logging.handlers import RotatingFileHandler
from datetime import datetime

from coastsat import SDS_download, SDS_shoreline, SDS_tools, SDS_transects

# ---------------- logging setup ----------------
LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO").upper()
LOG_DIR = os.getenv("LOG_DIR", "/run_data/logs")

logger = logging.getLogger("cs_runner")
logger.setLevel(LOG_LEVEL)

_formatter = logging.Formatter(
    fmt="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

class _FlushingStreamHandler(logging.StreamHandler):
    """Flushes after every record so logs appear live (e.g. `docker logs -f`)
    instead of sitting in a buffer until the process exits."""
    def emit(self, record):
        super().emit(record)
        self.flush()


class _FlushingFileHandler(RotatingFileHandler):
    def emit(self, record):
        super().emit(record)
        self.flush()


# stdout — always on, captured by `docker logs`
_console_handler = _FlushingStreamHandler(sys.stdout)
_console_handler.setFormatter(_formatter)
logger.addHandler(_console_handler)

# file — persists alongside output data, survives container exit
try:
    os.makedirs(LOG_DIR, exist_ok=True)
    _run_ts = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    _log_path = os.path.join(LOG_DIR, f"cs_runner_{_run_ts}.log")

    _file_handler = _FlushingFileHandler(
        _log_path, maxBytes=5 * 1024 * 1024, backupCount=3  # 5MB per file, keep last 3
    )
    _file_handler.setFormatter(_formatter)
    logger.addHandler(_file_handler)
    logger.info("Logging to file: %s", _log_path)
except OSError as e:
    # don't crash the whole run just because logs couldn't be written to disk
    logger.warning("Could not set up file logging at %s (%s) — stdout only", LOG_DIR, e)
# -------------------------------------------------


class CoastSatRunner():
    def __init__(
        self,
        startDate,
        endDate,
        savePath,
        coordinates,
        sitename,
        epsg,
        path_to_transects,
        path_to_tides,
        path_to_ref_shoreline
    ):
        self.startDate = startDate
        self.endDate = endDate
        self.savePath = savePath
        self.coordinates = coordinates
        self.sitename = sitename
        self.epsg = epsg
        self.path_to_transects = path_to_transects
        self.path_to_tides = path_to_tides
        self.path_to_ref_shoreline = path_to_ref_shoreline

        logger.debug(
            "Initialized CoastSatRunner: site=%s dates=%s->%s epsg=%s",
            self.sitename, self.startDate, self.endDate, self.epsg
        )

    def init_inputs(self):
        logger.info("Building inputs dict for site '%s'", self.sitename)
        polygon = SDS_tools.smallest_rectangle([self.coordinates])
        dates = [self.startDate, self.endDate]
        sat_list = ['L5', 'L7', 'L8', 'S2']
        collection = 'C02'

        inputs = {
            'polygon': polygon,
            'dates': dates,
            'sat_list': sat_list,
            'sitename': self.sitename,
            'filepath': f"/data/{self.sitename}",
            'landsat_collection': collection
        }
        logger.debug("Inputs: %s", inputs)

        return inputs

    def init_settings(self):
        logger.info("Loading reference shoreline from %s", self.path_to_ref_shoreline)
        try:
            with open(self.path_to_ref_shoreline, 'rb') as f:
                ref_shoreline_coords = pickle.load(f)
        except Exception as e:
            logger.error("Failed to load reference shoreline: %s", e)
            raise

        logger.debug("Building settings dict")
        settings = {
            # general parameters:
            'cloud_thresh': 0.1,        # threshold on maximum cloud cover
            'dist_clouds': 300,         # ditance around clouds where shoreline can't be mapped
            'output_epsg': self.epsg,       # epsg code of spatial reference system desired for the output

            # quality control:
            'check_detection': False,    # if True, shows each shoreline detection to the user for validation
            'adjust_detection': False,  # if True, allows user to adjust the postion of each shoreline by changing the threhold
            'save_figure': False,        # if True, saves a figure showing the mapped shoreline for each image

            # [ONLY FOR ADVANCED USERS] shoreline detection parameters:
            'min_beach_area': 1000,     # minimum area (in metres^2) for an object to be labelled as a beach
            'min_length_sl': 500,       # minimum length (in metres) of shoreline perimeter to be valid
            'cloud_mask_issue': False,  # switch this parameter to True if sand pixels are masked (in black) on many images
            'sand_color': 'default',    # 'default', 'latest', 'dark' (for grey/black sand beaches) or 'bright' (for white sand beaches)
            'pan_off': False,           # True to switch pansharpening off for Landsat 7/8/9 imagery
            's2cloudless_prob': 40,      # threshold to identify cloud pixels in the s2cloudless probability mask

            # add the inputs defined previously
            'inputs': self.inputs,

            # reference shoreline
            'reference_shoreline': ref_shoreline_coords,
            'max_dist_ref': 100
        }

        return settings

    def extract_shorelines(self, metadata, settings):
        logger.info("Extracting shorelines from imagery")
        output = SDS_shoreline.extract_shorelines(metadata, settings)
        logger.info("Extracted %d shoreline entries (pre-cleanup)", len(output.get('dates', [])))

        output = SDS_tools.remove_duplicates(output)
        output = SDS_tools.remove_inaccurate_georef(output, 10)
        logger.info("Shoreline cleanup done: %d entries remain", len(output.get('dates', [])))

        return output

    def load_transect_geojson(self):
        logger.info("Loading transects from %s", self.path_to_transects)
        transects = SDS_tools.transects_from_geojson(self.path_to_transects)
        logger.debug("Loaded %d transects", len(transects))

        return transects

    def compute_transect_shoreline_intersects(self, output, transects):
        logger.info("Computing transect/shoreline intersections")
        settings_transects = {  # parameters for computing intersections
            'along_dist':          25,        # along-shore distance to use for computing the intersection
            'min_points':          3,         # minimum number of shoreline points to calculate an intersection
            'max_std':             15,        # max std for points around transect
            'max_range':           30,        # max range for points around transect
            'min_chainage':        -100,      # largest negative value along transect (landwards of transect origin)
            'multiple_inter':      'auto',    # mode for removing outliers ('auto', 'nan', 'max')
            'auto_prc':            0.1,      # percentage to use in 'auto' mode to switch from 'nan' to 'max'
        }
        cross_distance = SDS_transects.compute_intersection_QC(output, transects, settings_transects)
        logger.debug("Computed intersections for %d transects", len(cross_distance))

        return cross_distance

    def tidal_correction(self, output, cross_distance):
        logger.info("Applying tidal correction using %s", self.path_to_tides)
        try:
            tide_data = pd.read_csv(self.path_to_tides, parse_dates=['dates'])
        except Exception as e:
            logger.error("Failed to read tide data: %s", e)
            raise

        dates_ts = [pd.to_datetime(_).to_pydatetime() for _ in tide_data['dates'].dt.tz_localize(pytz.timezone('UTC'))]
        tides_ts = np.array(tide_data['tide'])

        dates_sat = output['dates']
        tides_sat = SDS_tools.get_closest_datapoint(dates_sat, dates_ts, tides_ts)

        # tidal correction along each transect
        reference_elevation = 0.7  # elevation at which you would like the shoreline time-series to be
        beach_slope = 0.1
        cross_distance_tidally_corrected = {}
        for key in cross_distance.keys():
            correction = (tides_sat - reference_elevation) / beach_slope
            cross_distance_tidally_corrected[key] = cross_distance[key] + correction

        out_dict = dict([])
        out_dict['dates'] = dates_sat
        for key in cross_distance_tidally_corrected.keys():
            out_dict['Transect ' + key] = cross_distance_tidally_corrected[key]
        df = pd.DataFrame(out_dict)
        logger.debug("Tidal correction complete: %d rows", len(df))

        return df

    def run(self):
        logger.info("=== Starting CoastSat run for site '%s' ===", self.sitename)

        self.inputs = self.init_inputs()
        self.settings = self.init_settings()

        logger.info("Retrieving imagery for dates %s to %s", self.startDate, self.endDate)
        metadata = SDS_download.retrieve_images(self.inputs)
        settings = self.init_settings()

        output = self.extract_shorelines(metadata, settings)
        transects = self.load_transect_geojson()
        cross_distance = self.compute_transect_shoreline_intersects(output, transects)
        tidal_corrected_df = self.tidal_correction(output, cross_distance)

        try:
            tidal_corrected_df.to_csv(self.savePath, sep=',')
        except Exception as e:
            logger.error("Failed extracting data due to error: %s", e)
            sys.exit(1)
        else:
            logger.info("File saved to %s", self.savePath)
            logger.info("=== Run complete for site '%s' ===", self.sitename)


def assertfile_type_and_exists(file_path, expected_extension, assert_exist=True):
    if assert_exist:
        exists = os.path.isfile(file_path)
        if not exists:
            logger.error("Can't find file %s", file_path)
            sys.exit(1)

    _, extension = os.path.splitext(file_path)
    is_correct_extension = extension == expected_extension
    if is_correct_extension:
        return True
    else:
        logger.error("%s has wrong extension. expected %s", file_path, expected_extension)
        sys.exit(1)


def assert_dir_exists(dir_path):
    exists = os.path.isdir(dir_path)
    if not exists:
        logger.error("Can't find dir %s", dir_path)
        sys.exit(1)


def initializeCoastSatRunner(_args) -> CoastSatRunner:
    parser = argparse.ArgumentParser(
        prog="Coastsat",
        description="process shoreline data"
    )

    parser.add_argument("startDate", help="in YYYY-mm-dd format")
    parser.add_argument("endDate", help="in YYYY-mm-dd format")
    parser.add_argument("save_path", help="save file path in csv")
    parser.add_argument("coordinates", help="an array of coordinates of the area polygon")
    parser.add_argument("sitename", help="sitename")
    parser.add_argument("epsg")
    parser.add_argument("transects", help="path to transects geojson file")
    parser.add_argument("tides", help="path to tide data csv file")
    parser.add_argument("ref_shoreline", help="path to ref shoreline")

    args = parser.parse_args(_args)
    logger.debug("Parsed args: %s", vars(args))

    try:
        coordinates = ast.literal_eval(args.coordinates)
    except Exception as e:
        logger.error("Error parsing coordinates: %s", e)
        sys.exit(1)

    path_to_transects = args.transects
    path_to_tides = args.tides
    path_to_shoreline = args.ref_shoreline
    save_path = args.save_path

    assertfile_type_and_exists(path_to_transects, ".geojson")
    assertfile_type_and_exists(path_to_tides, ".csv")
    assertfile_type_and_exists(path_to_shoreline, ".pkl")
    assertfile_type_and_exists(save_path, ".csv", assert_exist=False)

    logger.info("All input files validated successfully")

    coastSatRunner = CoastSatRunner(
        args.startDate,
        args.endDate,
        args.save_path,
        coordinates,
        args.sitename,
        args.epsg,
        path_to_transects,
        path_to_tides,
        path_to_shoreline
    )

    return coastSatRunner


if __name__ == "__main__":
    coastSatRunner = initializeCoastSatRunner(sys.argv[1:])
    coastSatRunner.run()