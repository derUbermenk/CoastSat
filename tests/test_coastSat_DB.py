from CoastSat import initializeCoastSatRunnerDB, assertfile_type_and_exists, CoastSatRunnerDB
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

def test_retrieve_area_geometry_coords_from_db():
    # init csRunner
    csRunner = CoastSatRunner()

    expected_coordinates = []
    coordinates = csRunner.retrieve_area_geometry_coords_from_db()

    assert expected_coordinates, coordinates


def test_retrieve_baseline_geometry_from_db():
    return

def test_retrieve_transects_from_db():
    return
