CREATE EXTENSION IF NOT EXISTS postgis;

CREATE TABLE Shorelines (
    sitename VARCHAR(20) PRIMARY KEY,
    loc VARCHAR(50) NOT NULL,
    NSM DOUBLE PRECISION,
    SCE DOUBLE PRECISION,
    LRR DOUBLE PRECISION,
    WLR DOUBLE PRECISION
    baseline geometry(LINESTRING) NOT NULL,
    area  geometry(POLYGON) NOT NULL
);

CREATE TABLE Profiles (
    id SERIAL PRIMARY KEY,
    shoreline_id INT NOT NULL REFERENCES Shorelines(id) ON DELETE CASCADE,
    record_date DATE NOT NULL,
    geom geometry(MULTILINESTRING)
); 

CREATE TABLE Transects (
    id SERIAL PRIMARY KEY,
    shoreline_id INT NOT NULL REFERENCES Shorelines(id) ON DELETE CASCADE,
    geom geometry(LINESTRING)
);

CREATE TABLE Intersects (
    id SERIAL PRIMARY KEY,
    profile_id INT NOT NULL REFERENCES Profiles(id) ON DELETE CASCADE,
    transect_id INT NOT NULL REFERENCES Transects(id) ON DELETE CASCADE,
    distance DOUBLE PRECISION NOT NULL,
    geom geometry(POINT)
);