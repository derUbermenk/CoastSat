"""
This module contains all the functions needed to download the satellite images 
from the Google Earth Engine server
    
Author: Kilian Vos, Water Research Laboratory, University of New South Wales
"""

# load modules
import os
import numpy as np
import matplotlib.pyplot as plt
import pdb

# earth engine modules
import ee
from urllib.request import urlretrieve
import zipfile
import copy

# additional modules
from datetime import datetime, timedelta
import pytz
import pickle
from skimage import morphology, transform
from scipy import ndimage

# CoastSat modules
from coastsat import SDS_preprocess, SDS_tools, gdal_merge

np.seterr(all='ignore') # raise/ignore divisions by 0 and nans


def download_tif(image, polygon, bandsId, filepath):
    """
    Downloads a .TIF image from the ee server and stores it in a temp file
        
    Arguments:
    -----------
    image: ee.Image
        Image object to be downloaded
    polygon: list
        polygon containing the lon/lat coordinates to be extracted
        longitudes in the first column and latitudes in the second column
    bandsId: list of dict
        list of bands to be downloaded
    filepath: location where the temporary file should be saved
    
    Returns:
    -----------
    Downloads an image in a file named data.tif     
      
    """
    
    url = ee.data.makeDownloadUrl(ee.data.getDownloadId({
        'image': image.serialize(),
        'region': polygon,
        'bands': bandsId,
        'filePerBand': 'false',
        'name': 'data',
        }))
    local_zip, headers = urlretrieve(url)
    with zipfile.ZipFile(local_zip) as local_zipfile:
        return local_zipfile.extract('data.tif', filepath)


def retrieve_images(inputs):
    """
    Downloads all images from Landsat 5, Landsat 7, Landsat 8 and Sentinel-2 
    covering the area of interest and acquired between the specified dates. 
    The downloaded images are in .TIF format and organised in subfolders, divided 
    by satellite mission. The bands are also subdivided by pixel resolution.
    
    KV WRL 2018
        
    Arguments:
    -----------
    inputs: dict with the following keys
        'sitename': str
            name of the site
        'polygon': list
            polygon containing the lon/lat coordinates to be extracted,
            longitudes in the first column and latitudes in the second column,
            there are 5 pairs of lat/lon with the fifth point equal to the first point:
            ```
            polygon = [[[151.3, -33.7],[151.4, -33.7],[151.4, -33.8],[151.3, -33.8],
            [151.3, -33.7]]]
            ```
        'dates': list of str
            list that contains 2 strings with the initial and final dates in 
            format 'yyyy-mm-dd':
            ```
            dates = ['1987-01-01', '2018-01-01']
            ```
        'sat_list': list of str
            list that contains the names of the satellite missions to include: 
            ```
            sat_list = ['L5', 'L7', 'L8', 'S2']
            ```
        'filepath_data': str
            filepath to the directory where the images are downloaded
    
    Returns:
    -----------
    metadata: dict
        contains the information about the satellite images that were downloaded:
        date, filename, georeferencing accuracy and image coordinate reference system 
           
    """
    
    # initialise connection with GEE server
    ee.Initialize()
    
    # read inputs dictionnary
    sitename = inputs['sitename']
    polygon = inputs['polygon']
    dates = inputs['dates']
    sat_list= inputs['sat_list']
    filepath_data = inputs['filepath']
    
    # format in which the images are downloaded
    suffix = '.tif'
 
    # initialize metadata dictionnary (stores information about each image)       
    metadata = dict([])
    
    # create a new directory for this site
    if not os.path.exists(os.path.join(filepath_data,sitename)):
        os.makedirs(os.path.join(filepath_data,sitename))
        
    print('Downloading images:')
    
    #=============================================================================================#
    # download L5 images
    #=============================================================================================#
    
    if 'L5' in sat_list or 'Landsat5' in sat_list:
        
        satname = 'L5'
        # create a subfolder to store L5 images
        filepath = os.path.join(filepath_data, sitename, satname, '30m')
        filepath_meta = os.path.join(filepath_data, sitename, satname, 'meta')
        if not os.path.exists(filepath):
            os.makedirs(filepath)
        if not os.path.exists(filepath_meta):
            os.makedirs(filepath_meta)
            
        # Landsat 5 collection
        count_loop = 0
        while count_loop < 1:
            try:
                input_col = ee.ImageCollection('LANDSAT/LT05/C01/T1_TOA')
                # filter by location and dates
                flt_col = input_col.filterBounds(ee.Geometry.Polygon(polygon)).filterDate(dates[0],dates[1])
                # get all images in the filtered collection
                im_all = flt_col.getInfo().get('features')
                count_loop = 1
            except: 
                count_loop = 0   
                
        # remove very cloudy images (>95% cloud)
        cloud_cover = [_['properties']['CLOUD_COVER'] for _ in im_all]
        if np.any([_ > 95 for _ in cloud_cover]):
            idx_delete = np.where([_ > 95 for _ in cloud_cover])[0]
            im_col = [x for k,x in enumerate(im_all) if k not in idx_delete]
        else:
            im_col = im_all
        n_img = len(im_col)
        # print how many images there are
        print('%s: %d images'%(satname,n_img))
       
       # loop trough images
        timestamps = []
        acc_georef = []
        filenames = []
        all_names = []
        im_epsg = []
        for i in range(n_img):
            count_loop = 0
            while count_loop < 1:
                try:
                    # find each image in ee database
                    im = ee.Image(im_col[i]['id'])
                    count_loop = 1
                except: 
                    count_loop = 0 
            # read metadata
            im_dic = im_col[i]
            # get bands
            im_bands = im_dic['bands']
            # get time of acquisition (UNIX time)
            t = im_dic['properties']['system:time_start']
            # convert to datetime
            im_timestamp = datetime.fromtimestamp(t/1000, tz=pytz.utc)
            timestamps.append(im_timestamp)
            im_date = im_timestamp.strftime('%Y-%m-%d-%H-%M-%S')
            # get EPSG code of reference system
            im_epsg.append(int(im_dic['bands'][0]['crs'][5:]))
            # get geometric accuracy
            if 'GEOMETRIC_RMSE_MODEL' in im_dic['properties'].keys():
                acc_georef.append(im_dic['properties']['GEOMETRIC_RMSE_MODEL'])
            else:
                acc_georef.append(12) # default value of accuracy (RMSE = 12m)
            # delete dimensions key from dictionnary, otherwise the entire image is extracted
            for j in range(len(im_bands)): del im_bands[j]['dimensions']
            # bands for L5
            ms_bands = [im_bands[0], im_bands[1], im_bands[2], im_bands[3], im_bands[4], im_bands[7]]
            # filenames for the images
            filename = im_date + '_' + satname + '_' + sitename + suffix
            # if two images taken at the same date add 'dup' in the name (duplicate)
            if any(filename in _ for _ in all_names):
                filename = im_date + '_' + satname + '_' + sitename + '_dup' + suffix 
            all_names.append(filename)
            filenames.append(filename)
            # download .TIF image
            count_loop = 0
            while count_loop < 1:
                try:
                    local_data = download_tif(im, polygon, ms_bands, filepath)
                    count_loop = 1
                except: 
                    count_loop = 0   
            # update filename
            try:
                os.rename(local_data, os.path.join(filepath, filename))
            except:
                os.remove(os.path.join(filepath, filename))
                os.rename(local_data, os.path.join(filepath, filename))
            # write metadata in .txt file
            filename_txt = filename.replace('.tif','')
            metadict = {'filename':filename,'acc_georef':acc_georef[i],
                        'epsg':im_epsg[i]}
            with open(os.path.join(filepath_meta,filename_txt + '.txt'), 'w') as f:
                for key in metadict.keys():
                    f.write('%s\t%s\n'%(key,metadict[key]))
                    
            print('\r%d%%' % (int(((i+1)/n_img)*100)), end='')
        print('')
        
        # sort metadata (downloaded images are sorted by date in directory)
        timestamps_sorted = sorted(timestamps)
        idx_sorted = sorted(range(len(timestamps)), key=timestamps.__getitem__)
        acc_georef_sorted = [acc_georef[j] for j in idx_sorted]
        filenames_sorted = [filenames[j] for j in idx_sorted]
        im_epsg_sorted = [im_epsg[j] for j in idx_sorted]
        # save into dict
        metadata[satname] = {'dates':timestamps_sorted, 'acc_georef':acc_georef_sorted,
                'epsg':im_epsg_sorted, 'filenames':filenames_sorted}       
    
    #=============================================================================================#
    # download L7 images
    #=============================================================================================#
    
    if 'L7' in sat_list or 'Landsat7' in sat_list:
        
        satname = 'L7'
        # create subfolders (one for 30m multispectral bands and one for 15m pan bands)
        filepath = os.path.join(filepath_data, sitename, 'L7')
        filepath_pan = os.path.join(filepath, 'pan')
        filepath_ms = os.path.join(filepath, 'ms')
        filepath_meta = os.path.join(filepath, 'meta')
        if not os.path.exists(filepath_pan):
            os.makedirs(filepath_pan)
        if not os.path.exists(filepath_ms):
            os.makedirs(filepath_ms)
        if not os.path.exists(filepath_meta):
            os.makedirs(filepath_meta)
            
        # landsat 7 collection
        count_loop = 0
        while count_loop < 1:
            try:
                input_col = ee.ImageCollection('LANDSAT/LE07/C01/T1_RT_TOA')
                # filter by location and dates
                flt_col = input_col.filterBounds(ee.Geometry.Polygon(polygon)).filterDate(dates[0],dates[1])
                # get all images in the filtered collection
                im_all = flt_col.getInfo().get('features')
                count_loop = 1
            except: 
                count_loop = 0  

        # remove very cloudy images (>95% cloud)
        cloud_cover = [_['properties']['CLOUD_COVER'] for _ in im_all]
        if np.any([_ > 95 for _ in cloud_cover]):
            idx_delete = np.where([_ > 95 for _ in cloud_cover])[0]
            im_col = [x for k,x in enumerate(im_all) if k not in idx_delete]
        else:
            im_col = im_all
        n_img = len(im_col)
        # print how many images there are
        print('%s: %d images'%(satname,n_img)) 
        
        # loop trough images
        timestamps = []
        acc_georef = []
        filenames = []
        all_names = []
        im_epsg = []
        for i in range(n_img):
            
            count_loop = 0
            while count_loop < 1:
                try:
                    # find each image in ee database
                    im = ee.Image(im_col[i]['id'])
                    count_loop = 1
                except: 
                    count_loop = 0            
            # read metadata
            im_dic = im_col[i]
            # get bands
            im_bands = im_dic['bands']
            # get time of acquisition (UNIX time)
            t = im_dic['properties']['system:time_start']
            # convert to datetime
            im_timestamp = datetime.fromtimestamp(t/1000, tz=pytz.utc)
            timestamps.append(im_timestamp)
            im_date = im_timestamp.strftime('%Y-%m-%d-%H-%M-%S')
            # get EPSG code of reference system
            im_epsg.append(int(im_dic['bands'][0]['crs'][5:]))
            # get geometric accuracy
            if 'GEOMETRIC_RMSE_MODEL' in im_dic['properties'].keys():
                acc_georef.append(im_dic['properties']['GEOMETRIC_RMSE_MODEL'])
            else:
                acc_georef.append(12) # default value of accuracy (RMSE = 12m)  
            # delete dimensions key from dictionnary, otherwise the entire image is extracted
            for j in range(len(im_bands)): del im_bands[j]['dimensions']   
            # bands for L7
            pan_band = [im_bands[8]]
            ms_bands = [im_bands[0], im_bands[1], im_bands[2], im_bands[3], im_bands[4], im_bands[9]] 
            # filenames for the images
            filename_pan = im_date + '_' + satname + '_' + sitename + '_pan' + suffix
            filename_ms = im_date + '_' + satname + '_' + sitename + '_ms' + suffix  
            # if two images taken at the same date add 'dup' in the name
            if any(filename_pan in _ for _ in all_names):
                filename_pan = im_date + '_' + satname + '_' + sitename + '_pan' + '_dup' + suffix
                filename_ms = im_date + '_' + satname + '_' + sitename + '_ms' + '_dup' + suffix 
            all_names.append(filename_pan)
            filenames.append(filename_pan)
            # download .TIF image
            count_loop = 0
            while count_loop < 1:
                try:
                    local_data_pan = download_tif(im, polygon, pan_band, filepath_pan)
                    local_data_ms = download_tif(im, polygon, ms_bands, filepath_ms)
                    count_loop = 1
                except: 
                    count_loop = 0  
            # update filename
            try:
                os.rename(local_data_pan, os.path.join(filepath_pan, filename_pan))
            except:
                os.remove(os.path.join(filepath_pan, filename_pan))
                os.rename(local_data_pan, os.path.join(filepath_pan, filename_pan))
            try:
                os.rename(local_data_ms, os.path.join(filepath_ms, filename_ms))
            except:
                os.remove(os.path.join(filepath_ms, filename_ms))
                os.rename(local_data_ms, os.path.join(filepath_ms, filename_ms))
            # write metadata in .txt file
            filename_txt = filename_pan.replace('_pan','').replace('.tif','')
            metadict = {'filename':filename_pan,'acc_georef':acc_georef[i],
                        'epsg':im_epsg[i]}
            with open(os.path.join(filepath_meta,filename_txt + '.txt'), 'w') as f:
                for key in metadict.keys():
                    f.write('%s\t%s\n'%(key,metadict[key]))
                    
            print('\r%d%%' % (int(((i+1)/n_img)*100)), end='')
        print('')
            
        # sort metadata (dowloaded images are sorted by date in directory)
        timestamps_sorted = sorted(timestamps)
        idx_sorted = sorted(range(len(timestamps)), key=timestamps.__getitem__)
        acc_georef_sorted = [acc_georef[j] for j in idx_sorted]
        filenames_sorted = [filenames[j] for j in idx_sorted]
        im_epsg_sorted = [im_epsg[j] for j in idx_sorted]
        # save into dict
        metadata[satname] = {'dates':timestamps_sorted, 'acc_georef':acc_georef_sorted,
                'epsg':im_epsg_sorted, 'filenames':filenames_sorted}        
        
    #=============================================================================================#
    # download L8 images
    #=============================================================================================#
    
    if 'L8' in sat_list or 'Landsat8' in sat_list:

        satname = 'L8'  
        # create subfolders (one for 30m multispectral bands and one for 15m pan bands)
        filepath = os.path.join(filepath_data, sitename, 'L8')
        filepath_pan = os.path.join(filepath, 'pan')
        filepath_ms = os.path.join(filepath, 'ms')
        filepath_meta = os.path.join(filepath, 'meta')
        if not os.path.exists(filepath_pan):
            os.makedirs(filepath_pan)
        if not os.path.exists(filepath_ms):
            os.makedirs(filepath_ms)
        if not os.path.exists(filepath_meta):
            os.makedirs(filepath_meta)
            
        # landsat 8 collection
        count_loop = 0
        while count_loop < 1:
            try:
                input_col = ee.ImageCollection('LANDSAT/LC08/C01/T1_RT_TOA')
                # filter by location and dates
                flt_col = input_col.filterBounds(ee.Geometry.Polygon(polygon)).filterDate(dates[0],dates[1])
                # get all images in the filtered collection
                im_all = flt_col.getInfo().get('features')
                count_loop = 1
            except: 
                count_loop = 0 
                
        # remove very cloudy images (>95% cloud)
        cloud_cover = [_['properties']['CLOUD_COVER'] for _ in im_all]
        if np.any([_ > 95 for _ in cloud_cover]):
            idx_delete = np.where([_ > 95 for _ in cloud_cover])[0]
            im_col = [x for k,x in enumerate(im_all) if k not in idx_delete]
        else:
            im_col = im_all
        n_img = len(im_col)
        # print how many images there are
        print('%s: %d images'%(satname,n_img)) 
        
       # loop trough images
        timestamps = []
        acc_georef = []
        filenames = []
        all_names = []
        im_epsg = []
        for i in range(n_img):
            
            count_loop = 0
            while count_loop < 1:
                try:
                    # find each image in ee database
                    im = ee.Image(im_col[i]['id'])
                    count_loop = 1
                except: 
                    count_loop = 0             
            # read metadata
            im_dic = im_col[i]
            # get bands
            im_bands = im_dic['bands']
            # get time of acquisition (UNIX time)
            t = im_dic['properties']['system:time_start']
            # convert to datetime
            im_timestamp = datetime.fromtimestamp(t/1000, tz=pytz.utc)
            timestamps.append(im_timestamp)
            im_date = im_timestamp.strftime('%Y-%m-%d-%H-%M-%S')
            # get EPSG code of reference system
            im_epsg.append(int(im_dic['bands'][0]['crs'][5:]))
            # get geometric accuracy
            if 'GEOMETRIC_RMSE_MODEL' in im_dic['properties'].keys():
                acc_georef.append(im_dic['properties']['GEOMETRIC_RMSE_MODEL'])
            else:
                acc_georef.append(12) # default value of accuracy (RMSE = 12m)
            # delete dimensions key from dictionnary, otherwise the entire image is extracted
            for j in range(len(im_bands)): del im_bands[j]['dimensions']   
            # bands for L8    
            pan_band = [im_bands[7]]
            ms_bands = [im_bands[1], im_bands[2], im_bands[3], im_bands[4], im_bands[5], im_bands[11]]
            # filenames for the images
            filename_pan = im_date + '_' + satname + '_' + sitename + '_pan' + suffix
            filename_ms = im_date + '_' + satname + '_' + sitename + '_ms' + suffix  
            # if two images taken at the same date add 'dup' in the name
            if any(filename_pan in _ for _ in all_names):
                filename_pan = im_date + '_' + satname + '_' + sitename + '_pan' + '_dup' + suffix
                filename_ms = im_date + '_' + satname + '_' + sitename + '_ms' + '_dup' + suffix 
            all_names.append(filename_pan)  
            filenames.append(filename_pan)
            # download .TIF image
            count_loop = 0
            while count_loop < 1:
                try:
                    local_data_pan = download_tif(im, polygon, pan_band, filepath_pan)
                    local_data_ms = download_tif(im, polygon, ms_bands, filepath_ms)
                    count_loop = 1
                except: 
                    count_loop = 0 

            # update filename
            try:
                os.rename(local_data_pan, os.path.join(filepath_pan, filename_pan))
            except:
                os.remove(os.path.join(filepath_pan, filename_pan))
                os.rename(local_data_pan, os.path.join(filepath_pan, filename_pan))
            try:
                os.rename(local_data_ms, os.path.join(filepath_ms, filename_ms))
            except:
                os.remove(os.path.join(filepath_ms, filename_ms))
                os.rename(local_data_ms, os.path.join(filepath_ms, filename_ms))
            # write metadata in .txt file
            filename_txt = filename_pan.replace('_pan','').replace('.tif','')
            metadict = {'filename':filename_pan,'acc_georef':acc_georef[i],
                        'epsg':im_epsg[i]}
            with open(os.path.join(filepath_meta,filename_txt + '.txt'), 'w') as f:
                for key in metadict.keys():
                    f.write('%s\t%s\n'%(key,metadict[key]))
                
            print('\r%d%%' % (int(((i+1)/n_img)*100)), end='')
        print('')
    
        # sort metadata (dowloaded images are sorted by date in directory)
        timestamps_sorted = sorted(timestamps)
        idx_sorted = sorted(range(len(timestamps)), key=timestamps.__getitem__)
        acc_georef_sorted = [acc_georef[j] for j in idx_sorted]
        filenames_sorted = [filenames[j] for j in idx_sorted]
        im_epsg_sorted = [im_epsg[j] for j in idx_sorted]
        
        metadata[satname] = {'dates':timestamps_sorted, 'acc_georef':acc_georef_sorted,
                'epsg':im_epsg_sorted, 'filenames':filenames_sorted}

    #=============================================================================================#
    # download S2 images
    #=============================================================================================#
    
    if 'S2' in sat_list or 'Sentinel2' in sat_list:

        satname = 'S2' 
        # create subfolders for the 10m, 20m and 60m multipectral bands
        filepath = os.path.join(filepath_data, sitename, 'S2')
        if not os.path.exists(os.path.join(filepath, '10m')):
            os.makedirs(os.path.join(filepath, '10m'))
        if not os.path.exists(os.path.join(filepath, '20m')):
            os.makedirs(os.path.join(filepath, '20m'))
        if not os.path.exists(os.path.join(filepath, '60m')):
            os.makedirs(os.path.join(filepath, '60m'))
        filepath_meta = os.path.join(filepath, 'meta')
        if not os.path.exists(filepath_meta):
            os.makedirs(filepath_meta)
            
        # Sentinel2 collection
        count_loop = 0
        while count_loop < 1:
            try:
                input_col = ee.ImageCollection('COPERNICUS/S2')
                # filter by location and dates
                flt_col = input_col.filterBounds(ee.Geometry.Polygon(polygon)).filterDate(dates[0],dates[1])
                # get all images in the filtered collection
                im_all = flt_col.getInfo().get('features') 
                count_loop = 1
            except: 
                count_loop = 0 
 
        # remove duplicates in the collection (there are many in S2 collection)
        timestamps = [datetime.fromtimestamp(_['properties']['system:time_start']/1000,
                                             tz=pytz.utc) for _ in im_all]
        # utm zone projection
        utm_zones = np.array([int(_['bands'][0]['crs'][5:]) for _ in im_all])
        utm_zone_selected =  np.max(np.unique(utm_zones))
        # find the images that were acquired at the same time but have different utm zones
        idx_all = np.arange(0,len(im_all),1)
        idx_covered = np.ones(len(im_all)).astype(bool)
        idx_delete = []
        i = 0
        while 1:
            same_time = np.abs([(timestamps[i]-_).total_seconds() for _ in timestamps]) < 60*60*24
            idx_same_time = np.where(same_time)[0]
            same_utm = utm_zones == utm_zone_selected
            idx_temp = np.where([same_time[j] == True and same_utm[j] == False for j in idx_all])[0]
            idx_keep = idx_same_time[[_ not in idx_temp for _ in idx_same_time ]]
            # if more than 2 images with same date and same utm, drop the last ones
            if len(idx_keep) > 2: 
               idx_temp = np.append(idx_temp,idx_keep[-(len(idx_keep)-2):])
            for j in idx_temp:
                idx_delete.append(j)
            idx_covered[idx_same_time] = False
            if np.any(idx_covered):
                i = np.where(idx_covered)[0][0]
            else:
                break
        # update the collection by deleting all those images that have same timestamp and different
        # utm projection
        im_all_updated = [x for k,x in enumerate(im_all) if k not in idx_delete]
        
        # remove very cloudy images (>95% cloud)
        cloud_cover = [_['properties']['CLOUDY_PIXEL_PERCENTAGE'] for _ in im_all_updated]
        if np.any([_ > 95 for _ in cloud_cover]):
            idx_delete = np.where([_ > 95 for _ in cloud_cover])[0]
            im_col = [x for k,x in enumerate(im_all_updated) if k not in idx_delete]
        else:
            im_col = im_all_updated
        
        n_img = len(im_col)
        # print how many images there are
        print('%s: %d images'%(satname,n_img)) 
    
       # loop trough images
        timestamps = []
        acc_georef = []
        filenames = []
        all_names = []
        im_epsg = []
        for i in range(n_img):
            count_loop = 0
            while count_loop < 1:
                try:
                    # find each image in ee database
                    im = ee.Image(im_col[i]['id'])
                    count_loop = 1
                except: 
                    count_loop = 0             
            # read metadata
            im_dic = im_col[i]
            # get bands
            im_bands = im_dic['bands']
            # get time of acquisition (UNIX time)
            t = im_dic['properties']['system:time_start']
            # convert to datetime
            im_timestamp = datetime.fromtimestamp(t/1000, tz=pytz.utc)
            im_date = im_timestamp.strftime('%Y-%m-%d-%H-%M-%S')   
            # delete dimensions key from dictionnary, otherwise the entire image is extracted
            for j in range(len(im_bands)): del im_bands[j]['dimensions'] 
            # bands for S2
            bands10 = [im_bands[1], im_bands[2], im_bands[3], im_bands[7]]
            bands20 = [im_bands[11]]
            bands60 = [im_bands[15]]    
            # filenames for images
            filename10 = im_date + '_' + satname + '_' + sitename + '_' + '10m' + suffix
            filename20 = im_date + '_' + satname + '_' + sitename + '_' + '20m' + suffix
            filename60 = im_date + '_' + satname + '_' + sitename + '_' + '60m' + suffix
            # if two images taken at the same date skip the second image (they are the same)
            if any(filename10 in _ for _ in all_names):
                filename10 = filename10[:filename10.find('.')] + '_dup' + suffix
                filename20 = filename20[:filename20.find('.')] + '_dup' + suffix
                filename60 = filename60[:filename60.find('.')] + '_dup' + suffix
            all_names.append(filename10)  
            filenames.append(filename10)
            
            # download .TIF image and update filename
            count_loop = 0
            while count_loop < 1:
                try:
                    local_data = download_tif(im, polygon, bands10, os.path.join(filepath, '10m'))
                    count_loop = 1
                except: 
                    count_loop = 0  
                    
            try:
                os.rename(local_data, os.path.join(filepath, '10m', filename10))
            except:
                os.remove(os.path.join(filepath, '10m', filename10))
                os.rename(local_data, os.path.join(filepath, '10m', filename10))
                
            count_loop = 0
            while count_loop < 1:
                try:
                    local_data = download_tif(im, polygon, bands20, os.path.join(filepath, '20m'))
                    count_loop = 1
                except: 
                    count_loop = 0  
                    
            try:
                os.rename(local_data, os.path.join(filepath, '20m', filename20))
            except:
                os.remove(os.path.join(filepath, '20m', filename20))
                os.rename(local_data, os.path.join(filepath, '20m', filename20))
                
            count_loop = 0
            while count_loop < 1:
                try:
                    local_data = download_tif(im, polygon, bands60, os.path.join(filepath, '60m'))
                    count_loop = 1
                except: 
                    count_loop = 0 
                    
            try:
                os.rename(local_data, os.path.join(filepath, '60m', filename60))
            except:
                os.remove(os.path.join(filepath, '60m', filename60))
                os.rename(local_data, os.path.join(filepath, '60m', filename60))
                    
            # save timestamp, epsg code and georeferencing accuracy (1 if passed 0 if not passed)
            timestamps.append(im_timestamp)
            im_epsg.append(int(im_dic['bands'][0]['crs'][5:]))
            # Sentinel-2 products don't provide a georeferencing accuracy (RMSE as in Landsat)
            # but they have a flag indicating if the geometric quality control was passed or failed
            # if passed a value of 1 is stored if failed a value of -1 is stored in the metadata
            if 'GEOMETRIC_QUALITY_FLAG' in im_dic['properties'].keys():
                if im_dic['properties']['GEOMETRIC_QUALITY_FLAG'] == 'PASSED':
                    acc_georef.append(1)
                else:
                    acc_georef.append(-1)
            elif 'quality_check' in im_dic['properties'].keys():
                if im_dic['properties']['quality_check'] == 'PASSED':
                    acc_georef.append(1) 
                else:
                    acc_georef.append(-1)
            else:
                acc_georef.append(-1)
            # write metadata in .txt file
            filename_txt = filename10.replace('_10m','').replace('.tif','')
            metadict = {'filename':filename10,'acc_georef':acc_georef[i],
                        'epsg':im_epsg[i]}
            with open(os.path.join(filepath_meta,filename_txt + '.txt'), 'w') as f:
                for key in metadict.keys():
                    f.write('%s\t%s\n'%(key,metadict[key]))
                
            print('\r%d%%' % (int(((i+1)/n_img)*100)), end='')
        print('')
        # sort metadata (dowloaded images are sorted by date in directory)
        timestamps_sorted = sorted(timestamps)
        idx_sorted = sorted(range(len(timestamps)), key=timestamps.__getitem__)
        acc_georef_sorted = [acc_georef[j] for j in idx_sorted]
        filenames_sorted = [filenames[j] for j in idx_sorted]
        im_epsg_sorted = [im_epsg[j] for j in idx_sorted]
        
        metadata[satname] = {'dates':timestamps_sorted, 'acc_georef':acc_georef_sorted,
                'epsg':im_epsg_sorted, 'filenames':filenames_sorted} 
    
    # merge overlapping images (necessary only if the polygon is at the boundary of an image)
    if 'S2' in metadata.keys():
        metadata = merge_overlapping_images(metadata,inputs)

    # save metadata dict
    filepath = os.path.join(filepath_data, sitename)
    with open(os.path.join(filepath, sitename + '_metadata' + '.pkl'), 'wb') as f:
        pickle.dump(metadata, f)
    
    return metadata
        
            
def merge_overlapping_images(metadata,inputs):
    """
    Merge simultaneous overlapping images that cover the area of interest.
    When the area of interest is located at the boundary between 2 images, there 
    will be overlap between the 2 images and both will be downloaded from Google
    Earth Engine. This function merges the 2 images, so that the area of interest 
    is covered by only 1 image.
    
    KV WRL 2018
        
    Arguments:
    -----------
    metadata: dict
        contains all the information about the satellite images that were downloaded
    inputs: dict with the following keys
        'sitename': str
            name of the site
        'polygon': list
            polygon containing the lon/lat coordinates to be extracted,
            longitudes in the first column and latitudes in the second column,
            there are 5 pairs of lat/lon with the fifth point equal to the first point:
            ```
            polygon = [[[151.3, -33.7],[151.4, -33.7],[151.4, -33.8],[151.3, -33.8],
            [151.3, -33.7]]]
            ```
        'dates': list of str
            list that contains 2 strings with the initial and final dates in 
            format 'yyyy-mm-dd':
            ```
            dates = ['1987-01-01', '2018-01-01']
            ```
        'sat_list': list of str
            list that contains the names of the satellite missions to include: 
            ```
            sat_list = ['L5', 'L7', 'L8', 'S2']
            ```
        'filepath_data': str
            filepath to the directory where the images are downloaded
        
    Returns:
    -----------
    metadata_updated: dict
        updated metadata
            
    """
    
    # only for Sentinel-2 at this stage (not sure if this is needed for Landsat images)    
    sat = 'S2'
    filepath = os.path.join(inputs['filepath'], inputs['sitename'])
    filenames = metadata[sat]['filenames']
    # find the pairs of images that are within 5 minutes of each other
    time_delta = 5*60 # 5 minutes in seconds
    dates = metadata[sat]['dates'].copy()
    pairs = []
    for i,date in enumerate(metadata[sat]['dates']):
        # dummy value so it does not match it again
        dates[i] = pytz.utc.localize(datetime(1,1,1) + timedelta(days=i+1))
        # calculate time difference
        time_diff = np.array([np.abs((date - _).total_seconds()) for _ in dates])
        # find the matching times and add to pairs list
        boolvec = time_diff <= time_delta
        if np.sum(boolvec) == 0:
            continue
        else:
            idx_dup = np.where(boolvec)[0][0]
            pairs.append([i,idx_dup])
                
    # for each pair of image, create a mask and add no_data into the .tif file (this is needed before merging .tif files)
    for i,pair in enumerate(pairs):
        fn_im = []
        for index in range(len(pair)): 
            # get filenames of all the files corresponding to the each image in the pair
            fn_im.append([os.path.join(filepath, 'S2', '10m', filenames[pair[index]]),
                  os.path.join(filepath, 'S2', '20m',  filenames[pair[index]].replace('10m','20m')),
                  os.path.join(filepath, 'S2', '60m',  filenames[pair[index]].replace('10m','60m')),
                  os.path.join(filepath, 'S2', 'meta', filenames[pair[index]].replace('_10m','').replace('.tif','.txt'))])
            # read that image
            im_ms, georef, cloud_mask, im_extra, im_QA, im_nodata = SDS_preprocess.preprocess_single(fn_im[index], sat, False) 
            # im_RGB = SDS_preprocess.rescale_image_intensity(im_ms[:,:,[2,1,0]], cloud_mask, 99.9) 
            
            # in Sentinel2 images close to the edge of the image there are some artefacts, 
            # that are squares with constant pixel intensities. They need to be masked in the 
            # raster (GEOTIFF). It can be done using the image standard deviation, which 
            # indicates values close to 0 for the artefacts.      
            if len(im_ms) > 0:
                # calculate image std for the first 10m band
                im_std = SDS_tools.image_std(im_ms[:,:,0],1)
                # convert to binary
                im_binary = np.logical_or(im_std < 1e-6, np.isnan(im_std))
                # dilate to fill the edges (which have high std)
                mask10 = morphology.dilation(im_binary, morphology.square(3))
                # mask all 10m bands
                for k in range(im_ms.shape[2]):
                    im_ms[mask10,k] = np.nan
                # mask the 10m .tif file (add no_data where mask is True)
                SDS_tools.mask_raster(fn_im[index][0], mask10)
                
                # create another mask for the 20m band (SWIR1)
                im_std = SDS_tools.image_std(im_extra,1)
                im_binary = np.logical_or(im_std < 1e-6, np.isnan(im_std))
                mask20 = morphology.dilation(im_binary, morphology.square(3))     
                im_extra[mask20] = np.nan
                # mask the 20m .tif file (im_extra)
                SDS_tools.mask_raster(fn_im[index][1], mask20) 
                
                # use the 20m mask to create a mask for the 60m QA band (by resampling)
                mask60 = ndimage.zoom(mask20,zoom=1/3,order=0)
                mask60 = transform.resize(mask60, im_QA.shape, mode='constant', order=0,
                                          preserve_range=True)
                mask60 = mask60.astype(bool)
                # mask the 60m .tif file (im_QA)
                SDS_tools.mask_raster(fn_im[index][2], mask60)    
                            
            else:
                continue
            
            # make a figure for quality control
            # fig,ax= plt.subplots(2,2,tight_layout=True)
            # ax[0,0].imshow(im_RGB)
            # ax[0,0].set_title('RGB original')
            # ax[1,0].imshow(mask10)
            # ax[1,0].set_title('Mask 10m')
            # ax[0,1].imshow(mask20)  
            # ax[0,1].set_title('Mask 20m')
            # ax[1,1].imshow(mask60)
            # ax[1,1].set_title('Mask 60 m')
        
        # once all the pairs of .tif files have been masked with no_data, merge the using gdal_merge
        fn_merged = os.path.join(filepath, 'merged.tif')
        
        # merge masked 10m bands and remove duplicate file
        gdal_merge.main(['', '-o', fn_merged, '-n', '0', fn_im[0][0], fn_im[1][0]])
        os.chmod(fn_im[0][0], 0o777)
        os.remove(fn_im[0][0])
        os.chmod(fn_im[1][0], 0o777)
        os.remove(fn_im[1][0])
        os.chmod(fn_merged, 0o777)
        os.rename(fn_merged, fn_im[0][0])
        
        # merge masked 20m band (SWIR band)
        gdal_merge.main(['', '-o', fn_merged, '-n', '0', fn_im[0][1], fn_im[1][1]])
        os.chmod(fn_im[0][1], 0o777)
        os.remove(fn_im[0][1])
        os.chmod(fn_im[1][1], 0o777)
        os.remove(fn_im[1][1])
        os.chmod(fn_merged, 0o777)
        os.rename(fn_merged, fn_im[0][1])
    
        # merge QA band (60m band)
        gdal_merge.main(['', '-o', fn_merged, '-n', '0', fn_im[0][2], fn_im[1][2]])
        os.chmod(fn_im[0][2], 0o777)
        os.remove(fn_im[0][2])
        os.chmod(fn_im[1][2], 0o777)
        os.remove(fn_im[1][2])
        os.chmod(fn_merged, 0o777)
        os.rename(fn_merged, fn_im[0][2])
        
        # remove the metadata .txt file of the duplicate image
        os.chmod(fn_im[1][3], 0o777)
        os.remove(fn_im[1][3])
        
    print('%d pairs of overlapping Sentinel-2 images were merged' % len(pairs))
    
    # update the metadata dict
    metadata_updated = copy.deepcopy(metadata)
    idx_removed = []
    idx_kept = []
    for pair in pairs: idx_removed.append(pair[1])
    for idx in np.arange(0,len(metadata[sat]['dates'])):
        if not idx in idx_removed: idx_kept.append(idx)
    for key in metadata_updated[sat].keys():
        metadata_updated[sat][key] = [metadata_updated[sat][key][_] for _ in idx_kept]
        
    return metadata_updated  

def get_metadata(inputs):
    """
    Gets the metadata from the downloaded images by parsing .txt files located 
    in the \meta subfolder. 
    
    KV WRL 2018
        
    Arguments:
    -----------
    inputs: dict with the following fields
        'sitename': str
            name of the site
        'filepath_data': str
            filepath to the directory where the images are downloaded
    
    Returns:
    -----------
    metadata: dict
        contains the information about the satellite images that were downloaded:
        date, filename, georeferencing accuracy and image coordinate reference system 
           
    """
    # directory containing the images
    filepath = os.path.join(inputs['filepath'],inputs['sitename'])
    # initialize metadata dict
    metadata = dict([])
    # loop through the satellite missions
    for satname in ['L5','L7','L8','S2']:
        # if a folder has been created for the given satellite mission
        if satname in os.listdir(filepath):
            # update the metadata dict
            metadata[satname] = {'filenames':[], 'acc_georef':[], 'epsg':[], 'dates':[]}
            # directory where the metadata .txt files are stored
            filepath_meta = os.path.join(filepath, satname, 'meta')
            # get the list of filenames and sort it chronologically
            filenames_meta = os.listdir(filepath_meta)
            filenames_meta.sort()
            # loop through the .txt files
            for im_meta in filenames_meta:
                # read them and extract the metadata info: filename, georeferencing accuracy
                # epsg code and date
                with open(os.path.join(filepath_meta, im_meta), 'r') as f:
                    filename = f.readline().split('\t')[1].replace('\n','')
                    acc_georef = float(f.readline().split('\t')[1].replace('\n',''))
                    epsg = int(f.readline().split('\t')[1].replace('\n',''))
                date_str = filename[0:19]
                date = pytz.utc.localize(datetime(int(date_str[:4]),int(date_str[5:7]),
                                                  int(date_str[8:10]),int(date_str[11:13]),
                                                  int(date_str[14:16]),int(date_str[17:19])))
                # store the information in the metadata dict
                metadata[satname]['filenames'].append(filename)
                metadata[satname]['acc_georef'].append(acc_georef)
                metadata[satname]['epsg'].append(epsg)
                metadata[satname]['dates'].append(date)
                
    # save a .pkl file containing the metadata dict
    with open(os.path.join(filepath, inputs['sitename'] + '_metadata' + '.pkl'), 'wb') as f:
        pickle.dump(metadata, f)
    
    return metadata
                 