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
    ST_GeomFromText('LINESTRING(1007671.020112160709687 455030.685090125480201, 
                                 1007738.149961497285403 455136.269249796518125, 
                                 1007735.60449628578499 455475.019404369115364, 
                                 1007659.86397591792047 455679.772528785164468, 
                                 1007582.069948134245351 455898.714581696374808, 
                                 1007518.585046245483682 456095.360391221067403, 
                                 1007308.360074182623066 456271.509743563947268)'),
    ST_GeomFromText('POLYGON((144.79485033965136 13.429388175797682, 
                              144.80045079197197 13.428688997897156, 
                              144.78536604893787 13.42058047237599,  
                              144.78098868383339 13.421999744999741, 
                              144.79485033965136 13.429388175797682))')
);