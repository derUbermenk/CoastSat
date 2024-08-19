CREATE EXTENSION IF NOT EXISTS postgis;

CREATE TABLE Shorelines (
    sitename VARCHAR(20) PRIMARY KEY,
    loc VARCHAR(50),
    NSM DOUBLE PRECISION,
    SCE DOUBLE PRECISION,
    LRR DOUBLE PRECISION,
    WLR DOUBLE PRECISION,
    baseline geometry(LINESTRING) NOT NULL,
    area  geometry(POLYGON) NOT NULL
);

CREATE TABLE Profiles (
    id SERIAL PRIMARY KEY,
    shoreline_sitename VARCHAR(20) NOT NULL REFERENCES Shorelines(sitename) ON DELETE CASCADE,
    record_date DATE NOT NULL,
    geom geometry(MULTILINESTRING)
); 

CREATE TABLE Transects (
    id SERIAL PRIMARY KEY,
    shoreline_sitename VARCHAR(20) NOT NULL REFERENCES Shorelines(sitename) ON DELETE CASCADE,
    geom geometry(LINESTRING)
);

CREATE TABLE Intersects (
    id SERIAL PRIMARY KEY,
    profile_id INT NOT NULL REFERENCES Profiles(id) ON DELETE CASCADE,
    transect_id INT NOT NULL REFERENCES Transects(id) ON DELETE CASCADE,
    distance DOUBLE PRECISION NOT NULL,
    geom geometry(POINT)
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
);