from CoastSat import initializeCoastSatRunnerDB, CoastSatRunnerDB, Baseline 
from unittest.mock import Mock, patch

def test_initializeCoastSatRunner():
    startDate = "2024-01-01"
    endDate = "2024-02-01"
    sitename = "TEST1"
    tides = "/data/tides.csv"
    connstring = "appropriate_db_connstring"

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

def test_retrieve_area_geometry_coords_from_db():
    csRunner = CoastSatRunnerDB(
        "2024-01-01",
        "2024-02-01",
        "TEST1",
        "/data/tides.csv",
        "postgresql://shoreline:shoreline@localhost:5436/shoreline"
    )

    expected_coordinates = [ 
        [144.7948, 13.4293], 
        [144.8004, 13.4286], 
        [144.7853, 13.4205],  
        [144.7948, 13.4293], 
    ]

    coordinates = csRunner.retrieve_area_geometry_coords_from_db()

    assert expected_coordinates, coordinates


def test_retrieve_baseline_geometry_from_db():
    return

def test_retrieve_transects_from_db():
    return
