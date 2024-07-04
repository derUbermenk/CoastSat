#!/usr/bin/env python

import argparse
import pickle
import pandas as pd
import numpy as np
import os
import sys
import ast

from coastsat import SDS_download, SDS_shoreline, SDS_tools, SDS_transects

class CoastSatRunner():
    def __init__(
        self,
        startDate,
        endDate,
        saveDir,
        coordinates,
        sitename,
        epsg,
        path_to_transects,
        path_to_tides,
        path_to_ref_shoreline
    ):
        self.startDate = startDate
        self.endDate = endDate
        self.saveDir = saveDir
        self.coordinates = coordinates
        self.sitename = sitename
        self.epsg = epsg
        self.path_to_transects = path_to_transects 
        self.path_to_tides = path_to_tides
        self.path_to_ref_shoreline = path_to_ref_shoreline

    def init_inputs(self):
        polygon = SDS_tools.smallest_rectangle([self.coordinates])
        dates = [self.startDate, self.endDate]
        sat_list = ['L5','L7','L8']
        collection = 'C02'

        inputs = {
            'polygon': polygon,
            'dates': dates,
            'sat_list': sat_list,
            'sitename': self.sitename,
            'filepath': self.saveDir,
            'landsat_collection': collection
        }

        return inputs

    def init_settings(self):
        # get reference shorline
        with open(self.path_to_ref_shoreline, 'rb') as f:
            ref_shoreline_coords = pickle.load(f)

        settings = {
            # general parameters:
            'cloud_thresh': 0.1,        # threshold on maximum cloud cover
            'dist_clouds': 300,         # ditance around clouds where shoreline can't be mapped
            'output_epsg': 6637,       # epsg code of spatial reference system desired for the output

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
        # extract shorelines from all images (also saves output.pkl and shorelines.kml)
        output = SDS_shoreline.extract_shorelines(metadata, settings)

        # remove duplicates (images taken on the same date by the same satellite)
        output = SDS_tools.remove_duplicates(output)
        # remove inaccurate georeferencing (set threshold to 10 m)
        output = SDS_tools.remove_inaccurate_georef(output, 10)

        return output

    def load_transect_geojson(self):
        transects = SDS_tools.transects_from_geojson(self.path_to_transects)

        return transects

    def compute_transect_shoreline_intersects(self, output, transects):
        settings_transects = { # parameters for computing intersections
                            'along_dist':          25,        # along-shore distance to use for computing the intersection
                            'min_points':          3,         # minimum number of shoreline points to calculate an intersection
                            'max_std':             15,        # max std for points around transect
                            'max_range':           30,        # max range for points around transect
                            'min_chainage':        -100,      # largest negative value along transect (landwards of transect origin)
                            'multiple_inter':      'auto',    # mode for removing outliers ('auto', 'nan', 'max')
                            'auto_prc':            0.1,      # percentage to use in 'auto' mode to switch from 'nan' to 'max'
                            }
        cross_distance = SDS_transects.compute_intersection_QC(output, transects, settings_transects)
        return cross_distance
    
    def tidal_correction(self, output, cross_distance):
        tide_data = pd.read_csv(self.path_to_tides , parse_dates=['dates'])
        dates_ts = [pd.to_datetime(_).to_pydatetime() for _ in tide_data['dates']]
        tides_ts = np.array(tide_data['tide'])

        dates_sat = output['dates']
        tides_sat = SDS_tools.get_closest_datapoint(dates_sat, dates_ts, tides_ts)

        # tidal correction along each transect
        reference_elevation = 0.7 # elevation at which you would like the shoreline time-series to be
        beach_slope = 0.1
        cross_distance_tidally_corrected = {}
        for key in cross_distance.keys():
            correction = (tides_sat-reference_elevation)/beach_slope
            cross_distance_tidally_corrected[key] = cross_distance[key] + correction
        
        out_dict = dict([])
        out_dict['dates'] = dates_sat
        for key in cross_distance_tidally_corrected.keys():
            out_dict['Transect '+ key] = cross_distance_tidally_corrected[key]
        df = pd.DataFrame(out_dict)
        return df

    def run(self):
        self.inputs = self.init_inputs()
        self.settings = self.init_settings()

        # download images
        metadata = SDS_download.retrieve_images(self.inputs)
        settings = self.init_settings()

        output = self.extract_shorelines(metadata, settings)
        transects = self.load_transect_geojson()
        cross_distance = self.compute_transect_shoreline_intersects(output, transects)
        tidal_corrected_df = self.tidal_correction(output, cross_distance)

        # save to csv
        save_path = os.path.join(
            self.saveDir, f"{self.startDate}_{self.endDate}_data.csv"
        )
        try:
            tidal_corrected_df.to_csv(save_path, sep=',')
        except Exception as e:
            print(f"failed extracting data due to error \n\t{e}")
            sys.exit(1)
        else:
            print(f"file saved in \n\t{save_path}")

def assertfile_type_and_exists(file_path, expected_extension):
    exists = os.path.isfile(file_path)
    if exists:
        _, extension = os.path.splitext(file_path)
        is_correct_extension = extension == expected_extension

        if is_correct_extension:
            return True
        else:
            message = f"{file_path} has wrong extension. expected {expected_extension}"
            sys.exit((1, message))
    else: 
        sys.exit((1, f"cant find file {file_path}"))

def assert_dir_exists(dir_path):
    exists = os.path.isdir(dir_path)
    if not exists:
        sys.exit((1,f"cant find dir {dir_path}"))

def initializeCoastSatRunner(_args) ->  CoastSatRunner:
    parser = argparse.ArgumentParser(
        prog="Coastsat",
        description="process shoreline data"
    )    

    parser.add_argument("startDate", help="in YYYY-mm-dd format")
    parser.add_argument("endDate", help="in YYYY-mm-dd format")
    parser.add_argument("saveDir", help="save dir of csv file")
    parser.add_argument("coordinates", help="an array of coordinates of the area polygon")
    parser.add_argument("sitename", help="sitename")
    parser.add_argument("epsg")
    parser.add_argument("transects", help="path to transects geojson file")
    parser.add_argument("tides", help="path to tide data csv file")
    parser.add_argument("ref_shoreline", help="path to ref shoreline")

    args = parser.parse_args(_args)

    try:
        coordinates = ast.literal_eval(args.coordinates)
    except Exception as e:
        print(f"Error encountered: {e} \n Exiting") 
        sys.exit(1)

    path_to_transects = args.transects
    path_to_tides = args.tides
    path_to_shoreline = args.ref_shoreline   

    assert_dir_exists(args.saveDir)
    assertfile_type_and_exists(path_to_transects, ".geojson")
    assertfile_type_and_exists(path_to_tides, ".csv")
    assertfile_type_and_exists(path_to_shoreline, ".pkl")
    
    coastSatRunner = CoastSatRunner(
    args.startDate,
    args.endDate,
    args.saveDir,
    coordinates,
    args.sitename,
    args.epsg,
    path_to_transects,
    path_to_tides,
    path_to_shoreline
    )

    return coastSatRunner