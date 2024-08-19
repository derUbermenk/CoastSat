from CoastSat import initializeCoastSatRunnerDB, CoastSatRunnerDB, Shoreline
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
            [1007671.020112160709687, 455030.685090125480201],
            [1007738.149961497285403, 455136.269249796518125],
            [1007735.60449628578499, 455475.019404369115364],
            [1007659.86397591792047, 455679.772528785164468],
            [1007582.069948134245351, 455898.714581696374808],
            [1007518.585046245483682, 456095.360391221067403],
            [1007308.360074182623066, 456271.509743563947268]
    ]

    shoreline_area_coordinates=[[[ 
        [144.79485033965136, 13.429388175797682], 
        [144.80045079197197, 13.428688997897156], 
        [144.78536604893787, 13.42058047237599],  
        [144.78098868383339, 13.421999744999741]  
    ]]]

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
    assert base_shoreline.baseline == expected_base_shoreline.baseline
    assert base_shoreline.area == expected_base_shoreline.area

def test_retrieve_area_geometry_coords_from_db():
    csRunner = CoastSatRunnerDB(
        "2024-01-01",
        "2024-02-01",
        "TEST1",
        "/data/tides.csv",
        "postgresql://shoreline:shoreline@localhost:5436/shoreline"
    )

    expected_coordinates = [ 
        [144.79485033965136, 13.429388175797682], 
        [144.80045079197197, 13.428688997897156], 
        [144.78536604893787, 13.42058047237599],  
        [144.78098868383339, 13.421999744999741]  
    ]

    coordinates = csRunner.retrieve_area_geometry_coords_from_db()

    assert expected_coordinates, coordinates


def test_retrieve_baseline_geometry_from_db():
    return

def test_retrieve_transects_from_db():
    return
