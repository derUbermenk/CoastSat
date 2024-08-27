#!/usr/bin/env python

import argparse
import pickle
import pandas as pd
import pytz
import numpy as np
import os
import sys
import ast
import json
import geopandas as gpd

from sqlalchemy import create_engine, text 
from sqlalchemy.orm import sessionmaker 

from coastsat import SDS_download, SDS_shoreline, SDS_tools, SDS_transects

class Baseline():
    def __init__(self, sitename: str, baseline_geom: dict, area_geom: dict) -> None:
        self.sitename = sitename
        self.baseline_geom = baseline_geom
        self.area_geom = area_geom

class CoastSatRunnerDB():
    def __init__(
        self,
        startDate,
        endDate,
        sitename,
        epsg,
        path_to_tides,
        connstring
    ):
        self.startDate = startDate
        self.endDate = endDate
        self.sitename = sitename
        self.epsg=epsg 
        self.path_to_tides = path_to_tides 
        self.connstring = connstring

        base_shoreline = self.retrieve_base_shoreline() 

        self.baseline_geom = base_shoreline.baseline_geom
        self.area_geom = base_shoreline.area_geom
    

    def init_inputs(self):
        # polygon = SDS_tools.smallest_rectangle(self.area_geom['coordinates'])
        coordinates = [[-125.895220405324,49.1237726477147],[-125.88841138016,49.1127817966321],[-125.899425059767,49.1098655680256],[-125.906924940215,49.1205546121385],[-125.895220405324,49.1237726477147]]
        polygon = SDS_tools.smallest_rectangle([coordinates])
        dates = [self.startDate, self.endDate]
        sat_list = ['L5','L7','L8', 'S2']
        collection = 'C02'

        print(
f"""
used the following coords for polygon:
    {self.area_geom['coordinates']}

resulting polygon:
    {polygon}
""")

        inputs = {
            'polygon': polygon,
            'dates': dates,
            'sat_list': sat_list,
            'sitename': self.sitename,
            'filepath': f"/data/{self.sitename}",
            'landsat_collection': collection
        }

        return inputs

    def init_settings(self):
        # get reference shorline
        # cast to np array, coastsat shoreline requires this
        ref_shoreline_coords = np.array(self.baseline_geom['coordinates'])

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
        dates_ts = [pd.to_datetime(_).to_pydatetime() for _ in tide_data['dates'].dt.tz_localize(pytz.timezone('UTC'))]
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
    
    def retrieve_base_shoreline(self):
        engine = create_engine(self.connstring)
        session = sessionmaker(bind=engine)()

        # shoreline = session.query(Shoreline).filter(Shoreline.sitename==self.sitename).one()
        sql_query = text(
            """
                SELECT
                    sitename,
                    ST_AsGeoJSON(baseline) as baseline_geom,
                    ST_AsGeoJSON(area) as area_geom
                FROM shorelines
                WHERE sitename = :sitename
            """
        )
        result = session.execute(sql_query, {"sitename": self.sitename})
        row = result.fetchone()

        return Baseline(
            row[0],
            ast.literal_eval(row[1]),
            ast.literal_eval(row[2])
        ) 

    def retrieve_transects(self):
        query = """
            SELECT transect_name, ST_AsGeoJSON(geom) as geom_json
                FROM transects;
        """

        transects_df = pd.read_sql(query, self.connstring)
        transects_df['geom_array'] = transects_df['geom_json'].apply(self.geojson_to_numpy)
        transects_dict = transects_df.set_index('transect_name')['geom_array'].to_dict()
        return transects_dict

    def geojson_to_numpy(self, geojson_str):
        geom_dict = json.loads(geojson_str)
        coordinates = geom_dict['coordinates']
        return np.array(coordinates)


    def save_profiles_to_db(self,gdf: gpd.GeoDataFrame):
        df = pd.DataFrame(gdf)

        # drop unnecessary columns
        df['shoreline_sitename'] = self.sitename
        try:
            df['date'] = pd.to_datetime(df['date'])
        except Exception as e:
            print("\n")
            print(f"here are df columns: {df.columns} \n")
            print(f"here are df columns: {df} \n")
            print(f"here are gdf columns: {gdf.columns} \n")
            print("\n")
            raise e


        df['record_date'] = df['date'].dt.strftime('%Y-%m-%d')
        df = df.drop(columns=['date', 'geometry'])       

        engine = create_engine(self.connstring)
        df.to_sql(
            name='profiles',
            con=engine,
            if_exists='append',
            index=False,
            index_label='record_date'
        )
    
    def save_intersects_to_db(self,intersects):
        engine = create_engine(self.connstring)

        # get df with id, transect name that have shoreline sitename as sitename
        sql_query = text("SELECT id as transect_id, transect_name FROM transects WHERE shoreline_sitename = :sitename")
        params = { 'sitename': self.sitename }
        sitename_transects = pd.read_sql(sql_query, engine, params=params)

        intersects['dates'] = pd.to_datetime(intersects['dates'])
        intersects['profile_record_date'] = intersects['dates'].dt.strftime('%Y-%m-%d')
        intersects = intersects.drop(columns=['dates'])       

        intersects_melted = pd.melt(intersects, id_vars=['profile_record_date'], var_name='transect_name', value_name='distance')

        # see schema
        transects_intersects_melted = pd.merge(
            sitename_transects,           # Left DataFrame
            intersects_melted,            # Right DataFrame
            on='transect_name',           # Column to join on
            how='right'                    # Join type
        )

        keep_columns = ['profile_record_date', 'transect_id', 'distance']
        transects_intersects_melted = transects_intersects_melted[keep_columns]
        transects_intersects_melted.to_sql(
            name='intersects',
            con=engine,
            if_exists='append',
            index=False
        )

    
    def run(self):
        self.inputs = self.init_inputs()
        self.settings = self.init_settings()

        # download images
        metadata = SDS_download.retrieve_images(self.inputs)
        settings = self.init_settings()

        output = self.extract_shorelines(metadata, settings)
        try:
            raise ValueError
        except:
            print(f"\nhere is output: {output} \n")
        transects = self.retrieve_transects()
        cross_distance = self.compute_transect_shoreline_intersects(output, transects)
        tidal_corrected_df = self.tidal_correction(output, cross_distance)

        # save to csv
        gdf = SDS_tools.output_to_gdf(output, 'lines')
        self.save_profiles_to_db(gdf)
        try:
            tidal_corrected_df.to_csv(self.savePath, sep=',')
        except Exception as e:
            print(f"failed extracting data due to error \n\t{e}")
            sys.exit(1)
        else:
            print(f"file saved in \n\t{self.savePath}")

def assertfile_type_and_exists(file_path, expected_extension, assert_exist = True):
    if assert_exist:
        exists = os.path.isfile(file_path)
        if not exists:
            sys.exit((1, f"cant find file {file_path}"))

    _, extension = os.path.splitext(file_path)
    is_correct_extension = extension == expected_extension
    if is_correct_extension:
        return True
    else:
        message = f"{file_path} has wrong extension. expected {expected_extension}"
        sys.exit((1, message))

def assert_dir_exists(dir_path):
    exists = os.path.isdir(dir_path)
    if not exists:
        sys.exit((1,f"cant find dir {dir_path}"))

def initializeCoastSatRunnerDB(_args) ->  CoastSatRunnerDB:
    parser = argparse.ArgumentParser(
        prog="Coastsat DB",
        description="process shoreline data"
    )    

    parser.add_argument("startDate", help="in YYYY-mm-dd format")
    parser.add_argument("endDate", help="in YYYY-mm-dd format")
    parser.add_argument("sitename", help="sitename")
    parser.add_argument("epsg", help="epsg")
    parser.add_argument("tides", help="path to tide data csv file")
    parser.add_argument("connstring", help="db connstring")

    args = parser.parse_args(_args)

    assertfile_type_and_exists(args.tides, ".csv")

    coastSatRunner = CoastSatRunnerDB(
        args.startDate,
        args.endDate,
        args.sitename,
        args.epsg,
        args.tides,
        args.connstring
    )

    return coastSatRunner

if __name__ == "__main__":
    coastSatRunner = initializeCoastSatRunnerDB(sys.argv[1:])
    coastSatRunner.run()