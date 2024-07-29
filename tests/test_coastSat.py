from CoastSat import initializeCoastSatRunner, assertfile_type_and_exists, CoastSatRunner
from unittest.mock import Mock, patch

def test_assertfile_type_and_exists():
    with patch("os.path.isfile", return_value=True):
        # assuming file exists and is of correct type
        existing_file_of_correct_type = "/existing-dir/existing.csv"
        expected_extension = ".csv"
        try:
            assertfile_type_and_exists(existing_file_of_correct_type, expected_extension)
        except SystemExit as e:
            assert False 
        else:
            assert True 
        
        # assuming file exists and but not of correct type
        existing_file_of_wrong_type = "/existing-dir/existing.csv"
        expected_extension = ".geojson"
        try:
            assertfile_type_and_exists(existing_file_of_wrong_type, expected_extension)
        except SystemExit as e:
            assert True
            code, message= e.code
            assert code == 1
            assert message == "/existing-dir/existing.csv has wrong extension. expected .geojson"
        else:
            assert False 
    
    # fails if file dne 
    non_existent_file_path = "/absolutely_non_existent_file_path-7018ab22-59e5-499d-a7f1-a1afc35c310f/file.csv"
    expected_extension = ".csv"
    try:
        assertfile_type_and_exists(non_existent_file_path, expected_extension)
    except SystemExit as e:
        assert True
        code, message= e.code
        assert code == 1
        assert message == f"cant find file {non_existent_file_path}"
    else:
        assert False 

def test_initializeCoastSatRunner():

    startDate = "2024-01-01"
    endDate = "2024-02-01"
    saveDir = "/data"
    coordinates = "[ \
        [144.79485033965136, 13.429388175797682], \
        [144.80045079197197, 13.428688997897156], \
        [144.78536604893787, 13.42058047237599],  \
        [144.78098868383339, 13.421999744999741]  \
    ]"
    sitename = "PAGO",
    epsg = "6637"
    transects = "/data/transects.geojson"
    tides = "/data/tides.csv"
    path_to_ref_shoreline = "/data/path_to_ref_shoreline.pkl"

    args = [
        startDate,
        endDate,
        saveDir,
        coordinates,
        sitename,
        epsg,
        transects,
        tides,
        path_to_ref_shoreline
    ]

    with patch("os.path.isfile", return_value=True):
        with patch("os.path.isdir", return_value=True):
            coastSatRunner = initializeCoastSatRunner(args)

            assert isinstance(coastSatRunner, CoastSatRunner)
            assert coastSatRunner.startDate == startDate   
            assert coastSatRunner.endDate == endDate
            assert coastSatRunner.saveDir == saveDir
            assert coastSatRunner.coordinates == [
                [144.79485033965136, 13.429388175797682],
                [144.80045079197197, 13.428688997897156],
                [144.78536604893787, 13.42058047237599],
                [144.78098868383339, 13.421999744999741] 
            ]
            assert coastSatRunner.sitename == sitename
            assert coastSatRunner.epsg == epsg
            assert coastSatRunner.path_to_transects == transects
            assert coastSatRunner.path_to_tides == tides


