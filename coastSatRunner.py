#!/usr/bin/env python

import argparse

def CoastSatRunner():
    def __init__():
        pass

def initializeCoastSatRunner(_args) ->  CoastSatRunner:
    parser = argparse.ArgumentParser(
        prog="Coastsat",
        description="process shoreline data"
    )    

    parser.add_argument("startdate", help="in YYYY-mm-dd format")
    parser.add_argument("enddate", help="in YYYY-mm-dd format")
    parser.add_argument("saveDir", help="save dir of csv file")
    parser.add_argument("coordinates", help="an array of coordinates")
    parser.add_argument("sitename", help="sitename")
    parser.add_argument("epsg")
    parser.add_argument("transects", help="path to transects geojson file")
    parser.add_argument("tides", help="path to tide data csv file")

    args = parser.parse_args(_args)

    coastSatRunner = CoastSatRunner(
    args.startdate,
    args.enddate,
    args.saveDir,
    args.coordinates,
    args.sitename,
    args.epsg,
    args.transects,
    args.tides
    )

    return coastSatRunner
