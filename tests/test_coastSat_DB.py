from CoastSat import initializeCoastSatRunnerDB, CoastSatRunnerDB, Baseline 
from unittest.mock import Mock, patch
from sqlalchemy import create_engine, text
import pandas as pd


def test_initializeCoastSatRunner():
    startDate = "2024-01-01"
    endDate = "2024-02-01"
    sitename = "TEST1"
    tides = "/data/tides.csv"
    connstring = "postgresql://shoreline:shoreline@localhost:5436/shoreline_test"

    args = [
        startDate,
        endDate,
        sitename,
        tides,
        connstring
    ]

    with patch("os.path.isfile", return_value=True):
        with patch("os.path.isdir", return_value=True):
            coastSatRunner = initializeCoastSatRunnerDB(args)

            assert isinstance(coastSatRunner, CoastSatRunnerDB)
            assert coastSatRunner.startDate == startDate   
            assert coastSatRunner.endDate == endDate
            assert coastSatRunner.sitename == sitename
            assert coastSatRunner.tides == tides
            assert coastSatRunner.connstring == connstring 

def test_retrieve_base_shoreline():
    shoreline_baseline_coordinates = [
            [1007671.0201, 455030.6850],
            [1007738.1499, 455136.2692],
            [1007735.6044, 455475.0194],
            [1007659.8639, 455679.7725],
            [1007582.0699, 455898.7145],
            [1007518.5850, 456095.3603],
            [1007308.3600, 456271.5097]
    ]

    shoreline_area_coordinates=[[ 
        [144.7948, 13.4293], 
        [144.8004, 13.4286], 
        [144.7853, 13.4205],  
        [144.7948, 13.4293], 
    ]]

    expected_baseline = {"type":"LineString","coordinates": shoreline_baseline_coordinates}
    expected_area = {"type":"Polygon","coordinates": shoreline_area_coordinates}

    expected_base_shoreline = Baseline("TEST1", expected_baseline, expected_area)

    csRunner = CoastSatRunnerDB(
        "2024-01-01",
        "2024-02-01",
        "TEST1",
        "/data/tides.csv",
        "postgresql://shoreline:shoreline@localhost:5436/shoreline_test"
    )

    base_shoreline = csRunner.retrieve_base_shoreline()
    assert base_shoreline.sitename == expected_base_shoreline.sitename
    assert base_shoreline.baseline_geom == expected_base_shoreline.baseline_geom
    assert base_shoreline.area_geom == expected_base_shoreline.area_geom

def test_save_profiles_to_db():

    # Define the data
    data = {
        'geometry': [None] * 7,  # Assuming geometry is not provided here, hence set to None
        'date': [
            '2019-12-02 23:56:06',
            '2019-12-16 00:06:02',
            '2019-12-17 23:56:06',
            '2019-12-21 00:06:02',
            '2019-12-26 00:06:03',
            '2019-12-27 23:56:07',
            '2019-12-31 00:06:02'
        ],
        'satname': ['S2'] * 7,
        'geoaccuracy': ['PASSED'] * 7,
        'cloud_cover': [
            0,
            0.41537095271372,
            0,
            0.0404583042161659,
            0,
            0,
            0.0552061495457722
        ]
    }

    # Create the DataFrame
    gdf = pd.DataFrame(data)

    connstring = "postgresql://shoreline:shoreline@localhost:5436/shoreline_test"
    csRunner = CoastSatRunnerDB(
        "2024-01-01",
        "2024-02-01",
        "TEST1",
        "/data/tides.csv",
        connstring 
    )

    # cleanup
    engine = create_engine(connstring)
    with engine.connect() as connection:
        connection.execute(text("DELETE FROM profiles"))
        results = connection.execute(text("SELECT * FROM profiles where shoreline_sitename = 'TEST1'"))
        rows = results.fetchall()
        assert not rows

        csRunner.save_profiles_to_db(gdf) 
        results = connection.execute(text("SELECT * FROM profiles where shoreline_sitename = 'TEST1'"))
        rows = results.fetchall()
        
        assert rows
        connection.execute(text("DELETE FROM profiles"))

def test_save_profiles_to_db():
    pass