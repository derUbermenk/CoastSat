CREATE EXTENSION IF NOT EXISTS postgis;

CREATE TABLE Shorelines (
    sitename VARCHAR(20) PRIMARY KEY,
    loc VARCHAR(50),
    NSM DOUBLE PRECISION,
    SCE DOUBLE PRECISION,
    LRR DOUBLE PRECISION,
    WLR DOUBLE PRECISION,
    -- crs to use for calculating detected shoreline positions
    --  also used as crs for baseline
    output_epsg VARCHAR(50) NOT NULL DEFAULT 'EPSG:4326',
    -- used as reference for determining shoreline. must use
    --  output crs
    baseline geometry(LINESTRING) NOT NULL,
    -- uses EPSG:4326
    area  geometry(POLYGON) NOT NULL
);

CREATE TABLE Profiles (
    record_date DATE PRIMARY KEY,
    shoreline_sitename VARCHAR(20) NOT NULL REFERENCES Shorelines(sitename) ON DELETE CASCADE,
    satname VARCHAR(5),
    geoaccuracy VARCHAR(10),
    cloud_cover DOUBLE PRECISION,
    geom geometry(MULTILINESTRING)
); 

CREATE TABLE Transects (
    id SERIAL PRIMARY KEY,
    transect_name VARCHAR(5) NOT NULL,
    shoreline_sitename VARCHAR(20) NOT NULL REFERENCES Shorelines(sitename) ON DELETE CASCADE,
    geom geometry(LINESTRING),
    CONSTRAINT unique_transect_shoreline UNIQUE(transect_name, shoreline_sitename)
);

CREATE TABLE Intersects (
    id SERIAL PRIMARY KEY,
    profile_record_date DATE NOT NULL REFERENCES Profiles(record_date) ON DELETE CASCADE,
    transect_id INT NOT NULL REFERENCES Transects(id) ON DELETE CASCADE,
    distance DOUBLE PRECISION NOT NULL,
    geom geometry(POINT),

    CONSTRAINT unique_id_transect_shoreline UNIQUE(id, profile_record_date, transect_id)
);

INSERT INTO Shorelines (sitename, baseline, area)
VALUES (
    'TEST1',
    ST_GeomFromText('LINESTRING(1007671.0201 455030.6850,
                                1007738.1499 455136.2692,
                                1007735.6044 455475.0194,
                                1007659.8639 455679.7725,
                                1007582.0699 455898.7145,
                                1007518.5850 456095.3603,
                                1007308.3600 456271.5097)'),
    ST_GeomFromText('POLYGON((144.7948 13.4293, 
                              144.8004 13.4286, 
                              144.7853 13.4205,  
                              144.7948 13.4293))')  
), (
    'TEST2',
    ST_Transform(
        ST_SetSRID(
            ST_GeomFromGeoJSON('{
            "type": "LineString", 
            "coordinates": [ 
                [ -125.89586421741501, 49.110203027846936 ], 
                [ -125.893666071682333, 49.11309892870819 ], 
                [ -125.894848314904664, 49.118233447220312 ], 
                [ -125.897599572756462, 49.12182124087591 ], 
                [ -125.900544509274937, 49.123062024048529 ] 
            ]
            }'), 4326
        ),
        3005
    ),
    ST_SetSRID(
        ST_GeomFromGeoJSON('{
            "type": "Polygon", 
            "coordinates": [ [ 
                [ -125.896629083633542, 49.109333285130973 ], 
                [ -125.89163391797301, 49.112641084424745 ], 
                [ -125.898995069493452, 49.125868077976676 ], 
                [ -125.903407551879965, 49.122147641927576 ], 
                [ -125.896629083633542, 49.109333285130973 ] 
            ] ]
        }'), 4326
    )
);

INSERT INTO Transects (id, transect_name, shoreline_sitename, geom)
VALUES
(1, 'T1', 'TEST2', ST_SetSRID(
        ST_GeomFromGeoJSON('{ "type": "LineString", "coordinates": [ [ 1007639.722662, 456334.414515000011306 ], [ 1006898.373719, 456033.745684999972582 ] ] }'),
        3005
    )
),
(2, 'T2', 'TEST2', ST_SetSRID(
        ST_GeomFromGeoJSON('{ "type": "LineString", "coordinates": [ [ 1007756.983505, 456045.288427 ], [ 1007015.634563, 455744.619597000011709 ] ] }'),
        3005
    )
),
(3, 'T3', 'TEST2', ST_SetSRID(
    ST_GeomFromGeoJSON('{ "type": "LineString", "coordinates": [ [ 1007918.593002, 455646.813371 ], [ 1007177.244059, 455346.14454 ] ] }'),
    3005
    )
),
(4, 'T4', 'TEST2', ST_SetSRID(
    ST_GeomFromGeoJSON('{ "type": "LineString", "coordinates": [ [ 1008108.766037, 455177.910164 ], [ 1007367.417094, 454877.241334000020288 ] ] }'),
    3005
    )
);

INSERT INTO Profiles (record_date, shoreline_sitename, satname, geoaccuracy, cloud_cover, geom)
VALUES
    ('2024-01-06', 'TEST2', NULL, NULL, NULL, NULL),
    ('2024-01-11', 'TEST2', NULL, NULL, NULL, NULL),
    ('2024-01-12', 'TEST2', NULL, NULL, NULL, NULL);