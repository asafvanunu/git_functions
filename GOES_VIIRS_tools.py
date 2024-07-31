# %% [markdown]
# # Tools for processing GOES and VIIRS data
# 
# **Name:** Asaf Vanunu  
# **Date:** July 30, 2024
# 
# ## Description
# In this notebook you can find all the functions for the project

# %%
import numpy as np
from PIL import Image
from numpy.lib.stride_tricks import as_strided
import pandas as pd
import geopandas as gpd
import shapely.geometry
import os
import itertools
import rasterio
import rioxarray
import xarray as xr
from datetime import datetime, timedelta, time
from rasterio import features
from rasterio.enums import MergeAlg
from rasterio.plot import show
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# %%
def open_VIIRS(path:"goes_data_fire path") -> "open_VIIRS shapefile":
    """Gets a path and returns a geodataframe of VIIRS
    Keyword arguments:
   path --  for example 'F:\\Project_data\\goes_data_creek
    """
    sub_folders = [os.path.join(path, folder) for folder in os.listdir(path)]
    
    VIIRS_index = sub_folders.index("{}\\{}".format(path, "VIIRS_pnt"))
    VIIRS_folder = sub_folders[VIIRS_index] ## Extract VIIRS folder

    VIIRS_files_list = [os.path.join(VIIRS_folder, file) for file in os.listdir(VIIRS_folder)] ## full path for VIIRS files
    VIIRS_file_index = os.listdir(VIIRS_folder).index("VIIRS.shp") ## get the index of VIIRS.shp
    VIIRS = gpd.read_file(VIIRS_files_list[VIIRS_file_index]) ## Read VIIRS
    
    return(VIIRS)

# %%
def open_bbox(path:"goes_data_fire path") -> "open_bbox shapefile":
    """Gets a path and return a geodataframe of bbox
    Keyword arguments:
   path --  for example 'F:\\Project_data\\goes_data_creek
    """
    sub_folders = [os.path.join(path, folder) for folder in os.listdir(path)]
    
    shp_index = sub_folders.index("{}\\{}".format(path, "shp"))
    shp_folder = sub_folders[shp_index] ## extract shapefile folder
    
    for index, file in enumerate(os.listdir(shp_folder)):
        if file.endswith(".shp"):
            shp_index = index ## Get BBOX index
        
    bbox_path = os.path.join(shp_folder, os.listdir(shp_folder)[shp_index])
    bbox = gpd.read_file(bbox_path) ## open BBOX
    
    
    return(bbox)

# %%
def get_time_from_VIIRS(VIIRS_file: "cliped VIIRS shapefile") -> "VIIRS file with a column of time/string":
    """Gets a VIIRS file returns a VIIRS file with a column of time string
    Keyword arguments:
    VIIRS_file -- VIIRS shapefile
    """
    VIIRS_dates = [] ## set empty list
    
    for index in range(len(VIIRS_file)): ## for index in VIIRS file
        date = VIIRS_file["ACQ_DATE"][index] ## take the date
        time = VIIRS_file["ACQ_TIME"][index] ## take the time
        
        ## set Hour and Minute
        H = time[:2]
        minute = time[2:]
        ## string it 
        date_time = "{} {}:{}".format(date, H, minute)
        
        ## append in the list
        VIIRS_dates.append(date_time)
        
    VIIRS_file["DATE_TIME"] = VIIRS_dates ## set as a new column
    
    return(VIIRS_file) ## return

# %%
def clip_VIIRS_to_bbox(path:"goes_data_fire path") -> "return a clipped VIIRS file":
    """Gets a path and returns a geodataframe of clipped VIIRS with a corrected time column into the fire bbox
    Keyword arguments:
   path --  for example 'F:\\Project_data\\goes_data_creek
    """
    ## open VIIRS and the bbox
    bbox = open_bbox(path)
    VIIRS = open_VIIRS(path)
    
    bbox = bbox.to_crs(VIIRS.crs) ## convert to the same CRS
    VIIRS_clip = VIIRS.overlay(bbox[['geometry']], how="intersection") ## clip to bbox
    
    VIIRS_time_column = get_time_from_VIIRS(VIIRS_clip)
    
    return(VIIRS_time_column)

# %%
def get_time_from_goes(file_name: "string of file path", form: "string\time", VIIRS_time: "Y/N", NetCDF:"Y/N") -> "time\string":
    
    """Gets a GOES file full name and returns the time in a string or time format
    Keyword arguments:
    file_name -- string file name for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_08_31_00_02.tif'
    form -- string to enter: string/time
    VIIRS_time - If it is a time coming from VIIRS then mark Y else mark N
    NetCDF - If it is a time coming from a NetCDF file then mark Y else mark N
    """
    
    if (form == "time") & (VIIRS_time == "Y"): ## if it is VIIRS time
        
        time_format = '%Y-%m-%d %H:%M' ## format
        
        converted_time = datetime.strptime(file_name, time_format) ## convert to time
        
        return(converted_time) ## return VIIRS time

    ## this is the GOES cropped tif part ##
    if (VIIRS_time == "N") & (NetCDF == "N") : ## if it is a cropped tif
        l = file_name.split("\\") ## get list of path parts
    
        string_date = l[-1].split(".")[0] ## get the string '2020_08_31_00_02'
    
        Y, M, D, H, minute = string_date.split("_")
    
        if (form == "string") & (VIIRS_time == "N") :
        
            string = "{}-{}-{} {}:{}".format(Y,M,D,H,minute)
        
            return(string)
    
        if (form == "time") & (VIIRS_time == "N"):
        
            string = "{}-{}-{} {}:{}".format(Y,M,D,H,minute)
        
            time_format = '%Y-%m-%d %H:%M' ## format
    
            converted_time = datetime.strptime(string, time_format) ## convert to time
        
            return(converted_time)
        
    ## This is the GOES NetCDF part ##
    if (VIIRS_time == "N") & (NetCDF == "Y"):
        
        ## Extract year, month, day, hour, minute from the file name
        x = xr.open_dataset(file_name)
        t = x.t.data
        ymd = str(t)[:10] ## subset the Year, month and day
        ymd = ymd.replace("-","_") ## fix syntax
        h = str(t)[11:16] ## subset h:m
        h = h.replace(":","_") ##fix syntax
        f_name = ymd +"_"+ h ## create file name
        
        Y, M, D, H, minute = f_name.split("_")
        ## return string
        if (form == "string"):
            
            string = "{}-{}-{} {}:{}".format(Y,M,D,H,minute)
        
            return(string)
        ## return time ##
        if (form == "time"):
            
            string = "{}-{}-{} {}:{}".format(Y,M,D,H,minute)
        
            time_format = '%Y-%m-%d %H:%M' ## format
    
            converted_time = datetime.strptime(string, time_format) ## convert to time
        
            return(converted_time)

# %%
def time_in_range(start, end, x): ## function that gets a start and end datetime and say if x time is inside True/False
    """Return true if x is in the range [start, end]"""
    if start <= end:
        return start <= x <= end
    else:
        return start <= x or x <= end

# %%
def filter_VIIRS(VIIRS_file: "cliped VIIRS shapefile", image_path: "image path", time_delta:"minutes") -> "VIIRS file filterd to image time":
    """Gets a VIIRS file returns a VIIRS file with a column of time string
    Keyword arguments:
    VIIRS_file -- VIIRS shapefile
    image_time_string -- image time string for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif' 
    time_delta -- number that describe the time delta for example 5 will be 5 minutes before and after the time of the image
    """
    im = image_path.split("\\")[-1]
    
    ## if it is a cropped GOES image (tif file)
    if (im.endswith(".tif") == True):
    
    ## We set the start and end time
        time_diff = time_delta
        start_time = get_time_from_goes(image_path, form="time", VIIRS_time="N", NetCDF="N") - timedelta(minutes = time_diff)
        end_time = get_time_from_goes(image_path, form="time", VIIRS_time="N", NetCDF="N") + timedelta(minutes = time_diff)
        
      ## set empty list
        VIIRS_filter_list = []
        for index in range(len(VIIRS_file)):
            d_t = VIIRS_file["DATE_TIME"].iloc[index]
            date_time = get_time_from_goes(d_t, form="time", VIIRS_time="Y", NetCDF="N")
            con = time_in_range(start_time, end_time, date_time)
            VIIRS_filter_list.append(con)
        
        VIIRS_filter = VIIRS_file[VIIRS_filter_list].reset_index(drop = True)
        return(VIIRS_filter)
    
    ## if it is a NetCDF GOES image (nc file)
    if (im.endswith(".nc") == True):
    
    ## We set the start and end time
        time_diff = time_delta
        start_time = get_time_from_goes(image_path, form="time", VIIRS_time="N", NetCDF="Y") - timedelta(minutes = time_diff)
        end_time = get_time_from_goes(image_path, form="time", VIIRS_time="N", NetCDF="Y") + timedelta(minutes = time_diff)
        
    ## set empty list
        VIIRS_filter_list = []
        for index in range(len(VIIRS_file)):
            d_t = VIIRS_file["DATE_TIME"].iloc[index]
            date_time = get_time_from_goes(d_t, form="time", VIIRS_time="Y", NetCDF="N")
            con = time_in_range(start_time, end_time, date_time)
            VIIRS_filter_list.append(con)
        
        VIIRS_filter = VIIRS_file[VIIRS_filter_list].reset_index(drop = True)
        return(VIIRS_filter)


# %%
def create_fp_list(path:"goes_data_fire path") -> "create GOES fire product list":
    """Gets a path and returns a full path list of GOES fire products
    Keyword arguments:
   path --  for example 'F:\\Project_data\\goes_data_creek
    """
    sub_folders = [os.path.join(path, folder) for folder in os.listdir(path)]
    
    fp_index = sub_folders.index("{}\\{}".format(path, "fire_product_out"))
    fp_folder = sub_folders[fp_index] ## extract cropped fire product folder
    
    fp_list = [os.path.join(fp_folder, fp) for fp in os.listdir(fp_folder)]

    return(fp_list)

# %%
def create_fp_NetCDF_list(path:"goes_data_fire path") -> "create GOES fire product list":
    """Gets a path and returns a full path list of GOES NetCDF fire products
    Keyword arguments:
   path --  for example 'F:\\Project_data\\goes_data_creek
    """
    sub_folders = [os.path.join(path, folder) for folder in os.listdir(path)]
    
    nc_index = sub_folders.index("{}\\{}".format(path, "raw_fire_product"))
    nc_folder = sub_folders[nc_index] ## extract raw fire product folder
    
    nc_list = [os.path.join(nc_folder, nc) for nc in os.listdir(nc_folder)]

    return(nc_list)

# %%
def filter_fp_list(path:"goes_data_fire path", 
                   VIIRS_file:"cliped VIIRS file",
                  time_delta:"minutes") -> "create filterd GOES fire product list with VIITS times":
    """Gets a path and returns a full path list of cropped tif GOES fire products after filtered with VIIRS times
    Keyword arguments:
   path --  for example 'F:\\Project_data\\goes_data_creek
   VIIRS_file -- "cliped VIIRS file"
   time_delta -- number that describes the time delta for example 5 will be 5 minutes before and after the time of the image
    """
    VIIRS_clip = VIIRS_file
    fp_list = create_fp_list(path)
    
    f_list = [] ## empty list
    for file in fp_list: ## for each file in the list
        VIIRS_filter = filter_VIIRS(VIIRS_clip, file, time_delta=time_delta) ## filter VIIRS points
        if (len(VIIRS_filter)>0) == True: ## If the VIIRS points length is more then zero then
            f_list.append(file) ## append the file name
    
    
    return(f_list)

# %%
def filter_NetCDF_list(path:"goes_data_fire path", 
                   VIIRS_file:"cliped VIIRS file",
                  time_delta:"minutes") -> "create a filtered NetCDF GOES fire product list with VIITS times":
    """Gets a path and returns a full path list of NetCDF GOES fire products after filtered with VIIRS times
    Keyword arguments:
   path --  for example 'F:\\Project_data\\goes_data_creek
   VIIRS_file -- "cliped VIIRS file"
   time_delta -- number that describes the time delta for example 5 will be 5 minutes before and after the time of the image
    """
    VIIRS_clip = VIIRS_file
    nc_list = create_fp_NetCDF_list(path)
    
    f_list = [] ## empty list
    for file in nc_list: ## for each file in the list
        VIIRS_filter = filter_VIIRS(VIIRS_clip, file, time_delta=time_delta) ## filter VIIRS points
        if (len(VIIRS_filter)>0) == True: ## If the VIIRS points length is more then zero then
            f_list.append(file) ## append the file name
    
    
    return(f_list)

# %%
def GOES_pixels(fp_file:"GOES cropped fire product path", CRS:"CRS") -> "return GOES pixels polygons":
    
    """Gets a GOES file path returns a geodataframe of the polygons and the values
    Keyword arguments:
   fp_file --  for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif'
   CRS -- a CRS of an image for exapmle "mcmi_src.rio.crs"
    """
    
    image_src = rioxarray.open_rasterio(fp_file) ## open image
    
    dataset = image_src ## set dataset

    data = dataset.values[0] ## take the values

    t = dataset.rio.transform() ## take the transformation matrix

    move_x = t[0]
    # t[4] is negative, as raster start upper left 0,0 and goes down
    # later for steps calculation (ymin=...) we use plus instead of minus
    move_y = t[4]

    height = dataset.shape[1] ## GOES dims for rows
    width = dataset.shape[2] ## GOES dims for cols

    polygons = [] ## empty list
    indices = list(itertools.product(range(width), range(height))) ## make each polygon in a list

    for x,y in indices:
        x_min, y_max = t * (x,y)
        x_max = x_min + move_x
        y_min = y_max + move_y
        polygons.append(shapely.geometry.box(x_min, y_min, x_max, y_max)) ## append geometry to each poly
        
    data_list = []
    for x,y in indices:
        data_list.append(data[y,x]) ## append data values for each poly
    
    ## create gdf
    gdf = gpd.GeoDataFrame(data=data_list, crs=CRS, geometry=polygons, columns=['value']) 
    
    gdf = gdf[gdf["value"]!= -1].reset_index(drop=False) ## remove nodata
    
    
    return(gdf)

# %%
def spatial_join_VIIRS_GOES(fp_file:"GOES cropped fire product path", pixels:"goes_pixels shp", filterd_VIIRS:"filted_VIIRS shp",
                            CRS:"crs", VIIRS_buffer:"buffer",GOES_buffer:"buffer",
                            pixels_group:"list of GOES fire product pixels group")-> "return false alarm table":
    """Gets a GOES pixels file and a filterd VIIRS file and returns a dataframe of false alarms for a specific image
    Keyword arguments:
   fp_file --  for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif' 
   pixels --  shapefile of GOES pixels
   filterd_VIIRS -- shapefile of VIIRS after it is filterd to a given image time
   CRS -- a CRS of an image for exapmle "mcmi_src.rio.crs"
   VIIRS_buffer -- a number in meters that takes into acount the buffer radius for VIIRS points - for example 2000 is 2km
   GOES_buffer -- a number in meters that takes into acount the buffer radius for GOES fire pixels - for example 2000 is 2km
   pixels_group --  list of GOES fire pixels for example group_1, group_2 and so on...
    """
    VIIRS_T = filterd_VIIRS.copy() ## set a copy of VIIRS filter
    VIIRS_T = filterd_VIIRS.to_crs(CRS) ## convert CRS
    VIIRS_T['counter'] = 1 ## set a counter column  
    
    ## So first we will do a spatial join between GOES pixels and VIIRS
    pixels_join = pixels.sjoin(VIIRS_T[["counter", "geometry"]], how = "left")
    
    ## we will remove some non-important columns
    pixels_join.pop("index_right")
    
    ## Now we will sum up the number of VIIRS pixels inside a GOES pixel
    pixels_dis = pixels_join.dissolve(by = "index",
                         aggfunc={
                             "index":"first",
                             "value":"first",
                             "counter":"sum"
                             })
    
    ## replace 0 in nan where the value is nan
    for index in range(len(pixels_dis)):
        if np.isnan(pixels_dis["value"][index]) == True:
            pixels_dis["counter"][index] = np.nan
            
    ## Creating buffers
    buff = gpd.GeoDataFrame(geometry=VIIRS_T.buffer(VIIRS_buffer), crs=CRS)
    ## Set counter column
    buff["buffer_counter"] = 1
    
    # Spatial Join
    pixels_buff_sj = pixels_dis.sjoin(buff[["buffer_counter", "geometry"]], how = "left").reset_index(drop=True)
    
    ## pop non important column
    pixels_buff_sj.pop("index_right")
    
    ## count buffers in each polygon
    pixels_buff_dis = pixels_buff_sj.dissolve(by = "index",
                         aggfunc={
                             "index":"first",
                             "value":"first",
                             "counter":"first",
                             "buffer_counter":"sum"
                             })
    
    ## replace 0 in nan where the value is nan
    for index in range(len(pixels_buff_dis)):
        if np.isnan(pixels_buff_dis["value"][index]) == True:
            pixels_buff_dis["buffer_counter"][index] = np.nan
            
    fire_pixels = pixels_group ## select fire pixels group
    
    con_filter = np.isin(pixels_buff_dis["value"], fire_pixels) ## take GOES fire pixels
    
    ## Take out the fire polygons
    GOES_fire_poly = pixels_buff_dis[con_filter][["value", "geometry"]].reset_index(drop=True)
    
    if (len(GOES_fire_poly)>0) == True: 
    ## take centroids and give them 3km buffer
        cent = GOES_fire_poly.centroid
        cent_buff= cent.buffer(GOES_buffer)

        ## create the GeoDataFrame of the GOES fire pixels
        cent_buff_poly = gpd.GeoDataFrame(data=GOES_fire_poly["value"], geometry=cent_buff, crs=CRS)

         ## Rename the column and just set a binary value for the buffer
        cent_buff_poly = cent_buff_poly.rename(columns={"value":"GOES_FP"})
        cent_buff_poly["GOES_FP"] = 1   

        ## spatial join between the pixels and GOES pixels buffer
        SJ = pixels_buff_dis.sjoin(cent_buff_poly[["GOES_FP", "geometry"]], how = "left").reset_index(drop=True)
        SJ.pop("index_right")

        ## Sum the number of intersections between GOES fire pixels and neighboring pixels
        SJ_dis = SJ.dissolve(by = "index",
                             aggfunc={
                                 "index":"first",
                                 "value":"first",
                                 "counter":"first",
                                 "buffer_counter":"first",
                                 "GOES_FP":"sum"
                             })

        ## replace 0 in nan where the value is nan
        for index in range(len(SJ_dis)):
            if np.isnan(SJ_dis["value"][index]) == True:
                SJ_dis["GOES_FP"][index] = np.nan

        # set empty column
        SJ_dis["false_alarm"] = 0
        ## set a copy
        tr = SJ_dis.copy()

        for index in range(len(tr)):
        ## ON
        ## If there is no fire pixels
        ##and there is no VIIRS points
        ## and there is no points buffer
        ##OR if there is no fire pixels
        ## and there is no VIIRS points
        ## and There are points buffer
        ## then it is ON
            if (np.isin(tr["value"][index], fire_pixels) == False) & \
                ((tr["counter"][index] > 0) == False) & \
                ((tr["buffer_counter"][index] > 0) == False) | \
                (np.isin(tr["value"][index], fire_pixels) == False) & \
                ((tr["counter"][index] > 0) == False) & \
                ((tr["buffer_counter"][index] > 0) == True):
                tr["false_alarm"][index] = 1

        ## OF
        ## if there is no fire pixel
        ## and there is VIIRS points
        ## and there is VIIRS buffer
        ## then OF
            if (np.isin(tr["value"][index], fire_pixels) == False) & \
                ((tr["counter"][index] > 0) == True) & \
                ((tr["buffer_counter"][index] > 0) == True):
                tr["false_alarm"][index] = 2

        ## CN
        ## If there is fire pixel
        ## and there is no VIIRS points
        ## and there is no VIIRS buffer
        ## then CN
            if (np.isin(tr["value"][index], fire_pixels) == True) & \
                ((tr["counter"][index] > 0) == False) & \
                ((tr["buffer_counter"][index] > 0) == False):
                tr["false_alarm"][index] = 3 

        ## CF
        ## if there is fire pixel
        ## and there is VIIRS buffer
        ## or if there is no fire pixel
        ## but there is GOES fire pixel buffer
        ## and there is VIIRS points
        ## then CF
            if (np.isin(tr["value"][index], fire_pixels) == True) & \
                ((tr["buffer_counter"][index] > 0) == True) | \
                (np.isin(tr["value"][index], fire_pixels) == False) & \
                ((tr["GOES_FP"][index] > 0) == True) & \
                ((tr["counter"][index] > 0) == True):
                tr["false_alarm"][index] = 4

        ## replace nan    
            if np.isnan(pixels_buff_dis["value"][index]) == True:
                tr["false_alarm"][index] = np.nan
        ## calculate false alarm
        df = calculate_false_alarm(spatial_join_gdf=tr, fp_file=fp_file)

        return(df)

# %%
def calculate_false_alarm(spatial_join_gdf:"GOES_pixels and VIIRS gdf spatial join", fp_file:"full image path")-> "false alarm df":
    """Gets a GOES pixels and VIIRS Spatial Join gdf and returns a false alarm gdf
    Keyword arguments:
   spatial_join_gdf --  shapefile of GOES pixels and VIIRS spatial join
   fp_file --  for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif'
   fire_name_path -- for example 'F:\\Project_data\\goes_data_creek
    """
    
    tr = spatial_join_gdf
    
    ## Calculate the number of each category
    n = 0
    ON = 0
    OF = 0
    CN = 0
    CF = 0

    for index in range(len(tr)):
        if np.isnan(tr["value"][index]) == False:
            n += 1
        if (tr["false_alarm"][index] == 1):
            ON +=1
        if (tr["false_alarm"][index] == 2):
            OF +=1
        if (tr["false_alarm"][index] == 3):
            CN +=1
        if (tr["false_alarm"][index] == 4):
            CF +=1
            
    prop_CE = CN/n
    prop_OE = OF/n
    CER = CN/(CN + CF)
    OER = OF/(OF + ON)
    
    fire_name = fp_file.split("\\")[-3]
    f_name = fp_file.split("\\")[-1]
    d_name = get_time_from_goes(file_name=fp_file, form="string", VIIRS_time="N")
    
    d = {"fire_name":[fire_name],"file_name":[f_name], "date":[d_name], "number_of_pixels_(n)":[n], "ON":[ON],
        "OF":[OF], "CN":[CN], "CF":[CF], "prop_comission_error":[prop_CE],
        "prop_omission_error":[prop_OE], "comission_error_rate":[CER], "omission_error_rate":[OER]}
    
    df = pd.DataFrame(data=d)
    
    return(df)

# %%
def rasterize_fp_nc(image_file:"GOES fire product path",filterd_VIIRS:"filted_VIIRS shp",
                   CRS:"CRS", VIIRS_points:"int number") -> "VIIRS rasterized image according to GOES dims":

    """Gets a GOES file and a filterd VIIRS file and return a rasterized VIIRS image
    Keyword arguments:
   image_file --  for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif' or NC 
   filterd_VIIRS -- shapefile of VIIRS after it is filterd to a given image time
   CRS -- a CRS of an image for exapmle "mcmi_src.rio.crs"
   VIIRS_points -- int number, the number of VIIRS points inside a GOES pixel that would count as fire, for example "5" would mean that a pixel below 5 VIIRS points would not be conciderd as fire
    """
    
    f_name = image_file.split("\\")[-1] ## take image path and split the last part - the file name
    
    VIIRS_T = filterd_VIIRS.copy() ## set a VIIRS copy
    
    VIIRS_T = VIIRS_T.to_crs(CRS) ## reproject to GOES
    
    VIIRS_T['counter'] = 1 ## set a counter column
    
    if f_name.endswith(".tif"): ## if the file is a cropped GOES
        
        src_fp = rioxarray.open_rasterio(image_file) ## open the file
        
        fp_rast = src_fp.values[0]
        
        geom_value = [(geom,value) for geom, value in zip(VIIRS_T.geometry, VIIRS_T['counter'])] ## take the geometry and values
        
        rasterized_fp = features.rasterize(geom_value, ## tkae the pointes
                                out_shape = src_fp.data[0].shape, ## take GOES shape
                                transform = src_fp.rio.transform(), ## take transfor matrix
                                all_touched = True, 
                                fill = 0,   # background value
                                merge_alg = MergeAlg.add, ## sum up the counter value
                                dtype = src_fp.data[0].dtype) ## set GOES dtype
        
        VIIRS_con = rasterized_fp<VIIRS_points ## set a condition take all the pixels that are smaller then VIITS point parm
        
        rasterized_fp[VIIRS_con] = 0 ## set the condition to zero. so if we enter "10" all the pixels below 10 will be set to zero
        
        rasterized_fp[np.isnan(fp_rast)] = np.nan ## set nan values as the original cropped image
        
        return(rasterized_fp) ## return the raserized VIIRS
    
    if f_name.endswith(".nc"): ## if it is a netcdf file 
        
        src_nc = rioxarray.open_rasterio(image_file) ## open it
        
        geom_value = [(geom,value) for geom, value in zip(VIIRS_T.geometry, VIIRS_T['counter'])] ## take the geometry and values
        
        rasterized_nc = features.rasterize(geom_value, ## tkae the pointes
                                out_shape = src_nc["Mask"][0].shape, ## take GOES shape
                                transform = src_nc.rio.transform(), ## take transfor matrix
                                all_touched = True, 
                                fill = 0,   # background value
                                merge_alg = MergeAlg.add, ## sum up the counter value
                                dtype = src_nc["Mask"][0].dtype) ## set GOES dtype
        
        VIIRS_con = rasterized_nc<VIIRS_points ## set a condition take all the pixels that are smaller then VIITS point parm
        
        rasterized_nc[VIIRS_con] = 0 ## set the condition to zero. so if we enter "10" all the pixels below 10 will be set to zero
        
        return(rasterized_nc) ## return NC file

# %%
def calculate_false_alarm_rasterized(rasterized_file:"raster file", distance:"int steps from fire pixel",
                                    image_path:"string file path",
                                    pixels_group:"GOES pixels group - list",
                                    return_image:"Y/N") -> "false alarm data frame":
    """Gets a rasterized file and return a false alarm gdf
    Keyword arguments:
   rasterize_file --  rasterized VIIRS file
   image_path --  for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif' or NC
   distance -- the steps from a fire pixel for calculating neighbors
   pixels_group --  list of GOES fire pixels for example group_1, group_2 and so on...
   return_image -- Y if we want to return image and N if not 
    """
    all_fire_pixels = [
    10, ## Processed fp (fire pixel)
    11, ## Saturated fp
    12, ## Cloud contaminated fp
    30, ## Temporally filtered Processed fp
    31, ## Temporally filtered Saturated fp
    32, ## Temporally filtered Cloud contaminated fp
    13, ## High probability fp
    33, ## Temporally filtered high probability fp
    14, ## medium probability fp
    34, ## Temporally filtered medium probability fp
    15, ## low probability fp
    35 ## Temporally filtered low probability fp 
    ]
    
    if distance == 1: ## set the size of the window 3x3
        kernel = 9
    
    if distance == 2: ## set the size of the window 5x5
        kernel = 25
        
    if distance == 3: ## set the size of the window 7x7
        kernel = 49
    
    f_name = image_path.split("\\")[-1] ## take file name
    
    if f_name.endswith(".tif"):
    
        src_fp = rioxarray.open_rasterio(image_path) ## open it with rasterio
    
        fp_rast = src_fp.values[0].copy() ## take the image 
        
        #process_con = np.any(np.isin(fp_rast,pixels_group)) ## this is only if there is this pixel group in the raster
        
        #if process_con == True: ## if it is true then:
        
        dummy_rast= fp_rast.copy() ## set a copy of the fire product
    
        row_i, col_j = np.where(rasterized_file > 0) ## get the index of each pixel in VIIRS raster that is larger then 0
    
        dummy_rast[np.isnan(dummy_rast) == False] = 1 ## set all the values in the dummy rast that are not nan as 1 (non-fire)
        
        for i, j in zip(row_i, col_j):## for row index i and col index j
            
            neighboors = cell_neighbors(arr=fp_rast, i= i, j= j, d=distance) ## take the GOES neighbors in the given distance
        
            con = np.any(np.isin(neighboors, all_fire_pixels)) ## if any one of the VIIRS pixel neighbors is in the all GOES pixels group
        
            if (con) == True: # if the con is true
                dummy_rast[i,j] = 1 ## set the dummy rast in that index as 1 (no-fire) (it used to be 4 i changed it 27/7/23)
            else: ## if there isnt any GOES fire pixels around VIIRS
                dummy_rast[i,j] = 2 ## set the dummy rast in that index as 2 (omission error by GOES)
                
        ### GOES PART ########
        
        GOES_row_i, GOES_col_j = np.where(np.isin(fp_rast,pixels_group)) ## take the row and col index for all GOES pixels that are fire
        ## in a given pixels group, for example g1, g2 etc etc...
        
        for i, j in zip(GOES_row_i, GOES_col_j): ## for row index i and col index j
            
            neighboors = cell_neighbors(arr=rasterized_file, i = i, j= j, d=distance) ## take the VIIRS neighbors in the given distance
                
            calc = np.sum(np.isin(neighboors,0)) + np.sum(np.isnan(neighboors)) ## calculate the sum of all the zero and nan of the neighbors
                
            if (np.any(neighboors>0)) == True: ## if any of the neighbors value is heigher than 0
                dummy_rast[i,j] = 4 ## set this pixel as fire
            if (calc == kernel) == True: ## if all of the VIIRS srounding neighbors dont have fire
                dummy_rast[i,j] = 3 ## set it to 3 - false alarm!
        
        ## calculate false alarm
        df = false_alarm_dummy_rast(dummy_rast=dummy_rast, image_path=image_path) 
        
        dummy_rast[np.isnan(fp_rast)] = np.nan ## set nan values according to fp raster
                
        if (return_image == "Y"): ## in case we want to also see the image to analyze it!
            return(df, dummy_rast)
        if (return_image == "N"):
            return(df)
    
    ######NETCDF#########
    
    if f_name.endswith(".nc"): ## if it is a NC file
    
        src_nc = rioxarray.open_rasterio(image_path) ## open it with rasterio
    
        nc_rast = src_nc["Mask"][0].values.copy() ## take the image 
        
        process_con = np.any(np.isin(nc_rast,pixels_group)) ## this is only if there is this pixel group in the raster
        
        if process_con == True:
        
            dummy_rast= nc_rast.copy()
    
            row_i, col_j = np.where(rasterized_file > 0) ## get the index of each pixel in VIIRS raster that is larger then 0
    
            dummy_rast[np.isnan(dummy_rast) == False] = 1 ## set all the values in the dummy rast that are not nan as 1 (non-fire)
        
            for i, j in zip(row_i, col_j):## for row index i and col index j
            
                neighboors = cell_neighbors(arr=nc_rast, i= i, j= j, d=distance) ## take the GOES neighbors in the given distance
        
                con = np.any(np.isin(neighboors, all_fire_pixels)) ## if any one of the VIIRS pixel neighbors is in the ALL GOES pixels group
        
                if (con) == True: # if the con is true
                    dummy_rast[i,j] = 1 ## set the dummy rast in that index as 4 (fire)
                else: ## if there isnt any GOES fire pixels around VIIRS
                    dummy_rast[i,j] = 2 ## set the dummy rast in that index as 2 (omission error by GOES)
                
        ### GOES PART ########
        
            GOES_row_i, GOES_col_j = np.where(np.isin(nc_rast,pixels_group)) ## take the row and col index for all GOES pixels that are fire
        ## in a given pixel group, for example g1, g2, g3 and so on... 
        
            for i, j in zip(GOES_row_i, GOES_col_j): ## for row index i and col index j
            
                neighboors = cell_neighbors(arr=rasterized_file, i = i, j= j, d=distance) ## take the VIIRS neighbors in the given distance
                
                calc = np.sum(np.isin(neighboors,0)) + np.sum(np.isnan(neighboors)) ## calculate the sum of all the zero and nan of the neighbors
                
                if (np.any(neighboors>0)) == True: ## if any of the neighbors value is heigher than 0
                    dummy_rast[i,j] = 4 ## set this pixel as fire
                if (calc == kernel) == True: ## if all of the VIIRS srounding neighbors dont have fire or nan
                    dummy_rast[i,j] = 3 ## set it to 3 - false alarm!
        
        ## calculate false alarm
            df = false_alarm_dummy_rast(dummy_rast=dummy_rast, image_path=image_path) 
        
            if (return_image == "Y"): ## in case we want to also see the image to analyze it!
                return(df, dummy_rast)
            if (return_image == "N"):
                return(df)

# %%
def calculate_false_alarm_rasterized_no_neighbors(rasterized_file:"raster file",
                                    image_path:"string file path",
                                    pixels_group:"GOES pixels group - list",
                                    return_image:"Y/N") -> "false alarm data frame without neighbors":
    """Gets a rasterized file and return a false alarm gdf that dosent take into account neighbors
    Keyword arguments:
   rasterize_file --  rasterized VIIRS file
   image_path --  for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif' or NC
   pixels_group --  list of GOES fire pixels for example group_1, group_2 and so on...
   return_image -- Y if we want to return image and N if not 
    """
    all_fire_pixels = [
    10, ## Processed fp (fire pixel)
    11, ## Saturated fp
    12, ## Cloud contaminated fp
    30, ## Temporally filtered Processed fp
    31, ## Temporally filtered Saturated fp
    32, ## Temporally filtered Cloud contaminated fp
    13, ## High probability fp
    33, ## Temporally filtered high probability fp
    14, ## medium probability fp
    34, ## Temporally filtered medium probability fp
    15, ## low probability fp
    35 ## Temporally filtered low probability fp 
    ]
     
        
    f_name = image_path.split("\\")[-1] ## take file name
    
    if f_name.endswith(".tif"):
    
        src_fp = rioxarray.open_rasterio(image_path) ## open it with rasterio
    
        fp_rast = src_fp.values[0].copy() ## take the image
        
        process_con = np.any(np.isin(fp_rast,pixels_group)) ## this is only if there is this pixel group in the raster
        
        if process_con == True:
        
            dummy_rast = fp_rast.copy() ## set a copy of the fire product
        
            dummy_rast[np.isnan(dummy_rast) == False] = 1 ## set all of the not nan values to one
        
            VIIRS_GOES_con_fire = (rasterized_file>0) & (np.isin(fp_rast, pixels_group)) ## if there is a GOES fire pixel and VIIRs fire pixel - before 03/09 it was all_fire_pixels instead of pixels group
        
            VIIRS_GOES_con_false = (rasterized_file==0) & (np.isin(fp_rast, pixels_group)) ## if there is GOES fire but no VIIRS fire
        
            VIIRS_GOES_con_omit = (rasterized_file>0) & (~np.isin(fp_rast, all_fire_pixels)) ## if there is no GOES fire and only VIIRS fire
        
            if np.any(VIIRS_GOES_con_fire): ## if there is a GOES fire pixel and VIIRS fire pixel
                dummy_rast[VIIRS_GOES_con_fire] = 4 #set it as fire
            if np.any(VIIRS_GOES_con_false): ## if there is GOES fire but no VIIRS fire
                dummy_rast[VIIRS_GOES_con_false] = 3 ## set it as false alarm - commission error
            if np.any(VIIRS_GOES_con_omit): ## if there is no GOES fire and only VIIRS fire
                dummy_rast[VIIRS_GOES_con_omit] = 2 ## set it as omission error

            ## calculate false alarm
            df = false_alarm_dummy_rast(dummy_rast=dummy_rast, image_path=image_path) 

            if (return_image == "Y"): ## in case we want to also see the image to analyze it!
                return(df, dummy_rast)
            if (return_image == "N"):
                return(df)
    
    ######NETCDF#########
    
    if f_name.endswith(".nc"): ## if it is a NC file
    
        src_nc = rioxarray.open_rasterio(image_path) ## open it with rasterio
    
        nc_rast = src_nc["Mask"][0].values.copy() ## take the image 
        
        process_con = np.any(np.isin(nc_rast,pixels_group)) ## this is only if there is this pixel group in the raster
        
        if process_con == True:
        
            dummy_rast= nc_rast.copy()
            
            dummy_rast[np.isnan(dummy_rast) == False] = 1 ## set all of the nan values to one
                
            VIIRS_GOES_con_fire = (rasterized_file>0) & (np.isin(nc_rast, pixels_group)) ## if there is a GOES fire pixel and VIIRs fire pixel
        
            VIIRS_GOES_con_false = (rasterized_file==0) & (np.isin(nc_rast, pixels_group)) ## if there is GOES fire but no VIIRS fire
        
            VIIRS_GOES_con_omit = (rasterized_file>0) & (~np.isin(nc_rast, all_fire_pixels)) ## if there is no GOES fire and only VIIRS fire
        
            if np.any(VIIRS_GOES_con_fire): ## if there is a GOES fire pixel and VIIRs fire pixel
                dummy_rast[VIIRS_GOES_con_fire] = 4 #set it as fire
            if np.any(VIIRS_GOES_con_false): ## if there is GOES fire but no VIIRS fire
                dummy_rast[VIIRS_GOES_con_false] = 3 ## set it as false alarm - comission error
            if np.any(VIIRS_GOES_con_omit): ## if there is no GOES fire and only VIIRS fire
                dummy_rast[VIIRS_GOES_con_omit] = 2 ## set it as omission error
    
            
        ## calculate false alarm
            df = false_alarm_dummy_rast(dummy_rast=dummy_rast, image_path=image_path) 
        
            if (return_image == "Y"): ## in case we want to also see the image to analyze it!
                return(df, dummy_rast)
            if (return_image == "N"):
                return(df)

# %%
def create_label_raster(rasterized_file:"raster file", distance:"int steps from fire pixel",
                                    image_path:"string file path",
                                    pixels_group:"GOES pixels group - list",
                                    save_image:"Y/N",
                                   dst:"out_dst of the image" = False) -> "labeled imgae":
    """Gets a rasterized file, image path and distance neighborehood and return a labeld image
    Keyword arguments:
   rasterize_file --  rasterized VIIRS file
   image_path --  for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif' or NC
   distance -- the steps from a fire pixel for calculating neighbors, distance = 1 is 3X3, distance = 2 is 5X5
   pixels_group --  list of GOES fire pixels for example group_1, group_2 and so on...
   save_image -- Y if we want to save the image and N if not. It will save it in the fire dir in a new folder cald "labeld_images"
   dst -- Y if we want to save the image and N if not. Just write the dir without file name for example 'F:\\Project_data\\goes_data_creek\\label'
    """
    
    if distance == 1: ## set the size of the window 3x3
        kernel = 9
    
    if distance == 2: ## set the size of the window 5x5
        kernel = 25
        
    if distance == 3: ## set the size of the window 7x7
        kernel = 49
    
    f_name = image_path.split("\\")[-1] ## take file name
    
    if f_name.endswith(".tif"):
    
        src_fp = rioxarray.open_rasterio(image_path) ## open it with rasterio
    
        fp_rast = src_fp.values[0].copy() ## take the image 
        
        process_con = np.any(np.isin(fp_rast,pixels_group)) ## this is only if there is this pixel group in the raster
        
        VIIRS_rast_con = np.any(rasterized_file>0) ## if there are VIIRS fire pixels
        
        if (process_con == True) & (VIIRS_rast_con == True): ## if it is true then:
        
            dummy_rast= fp_rast.copy() ## set a copy of the fire product
    
    
            dummy_rast[np.isnan(dummy_rast) == False] = 0 ## set all the values in the dummy rast as 0 (non-fire)
            dummy_rast[np.isnan(dummy_rast) == True] = 0 ## set all the values in the dummy rast as 0 (non-fire)
                
        ### GOES PART ########
        
            GOES_row_i, GOES_col_j = np.where(np.isin(fp_rast,pixels_group)) ## take the row and col index for all GOES pixels that are fire
        ## in a given pixels group, for example g1, g2 etc etc...
        
            for i, j in zip(GOES_row_i, GOES_col_j): ## for row index i and col index j
            
                neighboors = cell_neighbors(arr=rasterized_file, i = i, j= j, d=distance) ## take the VIIRS neighbors in the given distance
                
                if (np.any(neighboors>0)) == True: ## if any of the neighbors value is heigher than 0
                    dummy_rast[i,j] = 1 ## set this pixel as fire

            
            if (save_image == "Y"): ## if we want to save the image
                if os.path.exists(dst): ## if the dst exists
                    full_dst = os.path.join(dst, f_name) ## take the full path dst
                    im = Image.fromarray(dummy_rast) ## convert to imgae
                    im.save(full_dst) ## save
                    return(dummy_rast) ## and return the image
                else: ## if the dst not exists:
                    os.mkdir(dst)  ## create it
                    full_dst = os.path.join(dst, f_name) ## take the full file dst path
                    im = Image.fromarray(dummy_rast) ## convert to image
                    im.save(full_dst) ## save it
                    return(dummy_rast) ## return the image
    
            if (save_image == "N"): ## if save imgae is No
                return(dummy_rast) # Only return the image
    

# %%
def calculate_false_alarm_random_forest(VIIRS_rasterized_file:"raster file", distance:"int steps from fire pixel",
                                    RF_image_path:"string file path",
                                    return_image:"Y/N") -> "false alarm data frame":
    """Gets a rasterized file and return a false alarm gdf
    Keyword arguments:
   VIIRS_rasterize_file --  rasterized VIIRS file matching the same time of the random forest
   RF_image_path --  for example 'F:\\change_detection\\random_forest\\predictions_rf_temp\\fire_sawtooth\\2020_09_05_09_02.tif'
   distance -- the steps from a fire pixel for calculating neighbors. 1 is 3X3 and 2 is 5X5 and 3 is 7X7
   return_image -- Y if we want to return image and N if not 
    """
    
    if distance == 1: ## set the size of the window 3x3
        kernel = 9
    
    if distance == 2: ## set the size of the window 5x5
        kernel = 25
        
    if distance == 3: ## set the size of the window 7x7
        kernel = 49
    
    f_name = RF_image_path.split("\\")[-1] ## take file name
    
    if f_name.endswith(".tif"):
    
        src_rf = rioxarray.open_rasterio(RF_image_path) ## open the random forest with rasterio
    
        rf_rast = src_rf.values[0].copy() ## take the image 
        
        dummy_rast= rf_rast.copy() ## set a copy of the random forest image
        #### VIIRS PART#######
        row_i, col_j = np.where(VIIRS_rasterized_file > 0) ## get the index of each pixel in VIIRS raster that is larger then 0
    
        dummy_rast[np.isnan(dummy_rast) == False] = 1 ## set all the values in the dummy rast that are not nan as 1 (non-fire)
        
        for i, j in zip(row_i, col_j):## for row index i and col index j
            
            neighboors = cell_neighbors(arr=rf_rast, i= i, j= j, d=distance) ## take the random forest neighbors in the given distance
        
            con = np.any(neighboors>0) ## is any one of the VIIRS pixel neighbors is labeld fire by the random forest?
        
            if (con) == True: # if the con is true
                dummy_rast[i,j] = 1 ## set the dummy rast in that index as 4 (fire) ## changed at 3/8/23 to 1 - no fire
            else: ## if there isnt any random forst fire pixels around VIIRS
                dummy_rast[i,j] = 2 ## set the dummy rast in that index as 2 (omission error by random forest)
                
        ### RF PART ########
        
            RF_row_i, RF_col_j = np.where(rf_rast>0) ## take the row and col index for all Random Forest pixels that are fire
        
            for i, j in zip(RF_row_i, RF_col_j): ## for row index i and col index j
            
                neighboors = cell_neighbors(arr=VIIRS_rasterized_file, i = i, j= j, d=distance) ## take the VIIRS neighbors in the given distance
                
                calc = np.sum(np.isin(neighboors,0)) + np.sum(np.isnan(neighboors)) ## calculate the sum of all the zero and nan of the neighbors
                
                if (np.any(neighboors>0)) == True: ## if any of the neighbors value is higher than 0
                    dummy_rast[i,j] = 4 ## set this pixel as fire
                if (calc == kernel) == True: ## if all of the VIIRS surrounding neighbors don't have fire
                    dummy_rast[i,j] = 3 ## set it to 3 - false alarm!
        
        ## calculate false alarm
            df = false_alarm_dummy_rast(dummy_rast=dummy_rast, image_path=RF_image_path) 
        
            dummy_rast[np.isnan(rf_rast)] = np.nan ## set nan values according to fp raster
                
            if (return_image == "Y"): ## in case we want to also see the image to analyze it!
                return(df, dummy_rast)
            if (return_image == "N"):
                return(df)

# %%
def false_alarm_dummy_rast(dummy_rast:"raster named 'dummy_rast'",
                           image_path:"image path string") -> "gdf with false alarm":
    
    """Gets a dummy rasterized file and return a false alarm gdf
    Keyword arguments:
   dummy_rast --  raster dummy rast
   image_path --  for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif' or NC
    """
    
    n = np.sum(np.isnan(dummy_rast) == False) ## total number of pixels that are not nan
    
    ## Calcualte statistics
    
    ON = np.sum(dummy_rast == 1) 
    
    OF = np.sum(dummy_rast == 2)
    
    CN = np.sum(dummy_rast == 3)
    
    CF = np.sum(dummy_rast == 4)
    
    prop_CE = CN/n
    
    prop_OE = OF/n
    
    CER = CN/(CN+CF)
    
    OER = OF/(OF+ON)
    
    ## Get fire name and file name ##
    
    fire_name = image_path.split("\\")[-3]
    
    f_name = image_path.split("\\")[-1]
    
    if (f_name.endswith(".tif")): ## get the date if it is nc or tif
        
        d_name = get_time_from_goes(file_name=image_path, form="string", VIIRS_time="N",  NetCDF="N")
    
    if (f_name.endswith(".nc")):
        
        d_name = get_time_from_goes(file_name=image_path, form="string", VIIRS_time="N",  NetCDF="Y")
    
    # create dict
    d = {"fire_name":[fire_name],"file_name":[f_name], "date":[d_name], "number_of_pixels_(n)":[n], "ON":[ON],
        "OF":[OF], "CN":[CN], "CF":[CF], "prop_comission_error":[prop_CE],
        "prop_omission_error":[prop_OE], "comission_error_rate":[CER], "omission_error_rate":[OER]}
    
    # create df
    df = pd.DataFrame(data=d)
    
    return(df)

# %%
def sliding_window(arr:"array", window_size:"size of the window") -> "array sliced to the window":
    """ Construct a sliding window view of the array
    Keyword arguments:
   arr --  array, could be a raster image once it is open
   window_size --  for example "3" if we want a 3X3 window
    """
    arr = np.asarray(arr)
    window_size = int(window_size)
    if arr.ndim != 2:
        raise ValueError("need 2-D input")
    if not (window_size > 0):
        raise ValueError("need a positive window size")
    shape = (arr.shape[0] - window_size + 1,
             arr.shape[1] - window_size + 1,
             window_size, window_size)
    if shape[0] <= 0:
        shape = (1, shape[1], arr.shape[0], shape[3])
    if shape[1] <= 0:
        shape = (shape[0], 1, shape[2], arr.shape[1])
    strides = (arr.shape[1]*arr.itemsize, arr.itemsize,
               arr.shape[1]*arr.itemsize, arr.itemsize)
    return as_strided(arr, shape=shape, strides=strides)

# %%
def cell_neighbors(arr:"array", i:"int number index", j:"int number index", d:"int number distance") -> "cell neighbors":
    """Return d-th neighbors of cell (i, j)
     Keyword arguments:
   arr --  array, could be a raster image once it is open
   i --  row index
   j --  column index
   d --  the distance from the pixel. for example "1" will act like a 3X3 window and "2" will be a distance of 2 from the index
    """
    w = sliding_window(arr, 2*d+1)

    ix = np.clip(i - d, 0, w.shape[0]-1)
    jx = np.clip(j - d, 0, w.shape[1]-1)

    i0 = max(0, i - d - ix)
    j0 = max(0, j - d - jx)
    i1 = w.shape[2] - max(0, d - i + ix)
    j1 = w.shape[3] - max(0, d - j + jx)

    return w[ix, jx][i0:i1,j0:j1].ravel()

# %%
def img_is_color(img):

    if len(img.shape) == 3:
        # Check the color channels to see if they're all the same.
        c1, c2, c3 = img[:, : , 0], img[:, :, 1], img[:, :, 2]
        if (c1 == c2).all() and (c2 == c3).all():
            return True

    return False

def show_image_list(list_images, list_titles=None, list_cmaps=None, grid=True, num_cols=2, figsize=(20, 10), title_fontsize=30):
    '''
    Shows a grid of images, where each image is a Numpy array. The images can be either
    RGB or grayscale.

    Parameters:
    ----------
    images: list
        List of the images to be displayed.
    list_titles: list or None
        Optional list of titles to be shown for each image.
    list_cmaps: list or None
        Optional list of cmap values for each image. If None, then cmap will be
        automatically inferred.
    grid: boolean
        If True, show a grid over each image
    num_cols: int
        Number of columns to show.
    figsize: tuple of width, height
        Value to be passed to pyplot.figure()
    title_fontsize: int
        Value to be passed to set_title().
    '''

    assert isinstance(list_images, list)
    assert len(list_images) > 0
    assert isinstance(list_images[0], np.ndarray)

    if list_titles is not None:
        assert isinstance(list_titles, list)
        assert len(list_images) == len(list_titles), '%d imgs != %d titles' % (len(list_images), len(list_titles))

    if list_cmaps is not None:
        assert isinstance(list_cmaps, list)
        assert len(list_images) == len(list_cmaps), '%d imgs != %d cmaps' % (len(list_images), len(list_cmaps))

    num_images  = len(list_images)
    num_cols    = min(num_images, num_cols)
    num_rows    = int(num_images / num_cols) + (1 if num_images % num_cols != 0 else 0)

    # Create a grid of subplots.
    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize)
    
    # Create list of axes for easy iteration.
    if isinstance(axes, np.ndarray):
        list_axes = list(axes.flat)
    else:
        list_axes = [axes]

    for i in range(num_images):

        img    = list_images[i]
        title  = list_titles[i] if list_titles is not None else 'Image %d' % (i)
        cmap   = list_cmaps[i] if list_cmaps is not None else (None if img_is_color(img) else 'gray')
        
        list_axes[i].imshow(img, cmap=cmap)
        list_axes[i].set_title(title, fontsize=title_fontsize) 
        list_axes[i].grid(grid)
        
    for i in range(num_images, len(list_axes)):
        list_axes[i].set_visible(False)

    fig.tight_layout()
    _ = plt.show()

# %%
def create_i_j(i_list, j_list):
    """this function get a list of i and a list of j and return them as one list of i,j

    Args:
        i_list (list): a list of fire pixels locations in row (i)
        j_list (list): a list of fire pixels locations in col (j)
    """
    out_list = []
    for i,j in zip(i_list, j_list):
        loc = (i, j)
        out_list.append(loc)
    return(out_list)

# %%
def out_of_bounds_test(location, matrix_shape, test_type):
    """get a location and a matrix shape and determine if it is out of bounds, and return nan or the location of the fire

    Args:
        location (list): location list for example [2,4] ([row,col])
        matrix_shape (tupple): matrix shape for example (10,10)
        test_type (string): write if it's for row or col. e.g. "row", "col"
    """
    min_row = 0 ## min row
    min_col = 0 ## min col
    max_row = matrix_shape[0] - 1 ## take the max row and col
    max_col = matrix_shape[1] - 1 
    loc_row = location[0] ## take the row and col location
    loc_col = location[1]
    
    if test_type == "row": ## if the test is for row
        if (loc_row < min_row) | (loc_row > max_row): ## if it is out of bounds 
            return(np.nan) ## return nan
        else:
            return(location) ## else return the location
    if test_type == "col":
        if (loc_col < min_col) | (loc_col > max_col):
            return(np.nan)
        else:
            return(location)

# %%
def set_R(location_row, location_col, R_type, matrix_shape, matrix):
    """this function get location of row and location of col and the R type we want

    Args:
        location_row (int): location of row
        location_col (int): location of col
        R_type (string): for example "R1" or "R2" "R3" "R4"
        matrix_shape(tupple): for example (10,10)
        matrix(2d array): The VIIRS rasterize array
    """
    if R_type == "R0":
        return(matrix[location_row, location_col])
    
    if R_type == "R1":
     ## If V(i,j+1)≤V(i,j) then R(1) = V(i,j+1)+ V(i,j)   
        R1_con = out_of_bounds_test(location=(location_row, location_col+1), matrix_shape=matrix_shape, test_type="col")
        if np.all(np.isnan(R1_con)):
            return(np.nan)
        else:
            R1_matrix = matrix[location_row, location_col+1]
            regular_value = matrix[location_row, location_col]
            if  R1_matrix <= regular_value:
                return(R1_matrix + regular_value)
            else:
                return(regular_value)
            
    if R_type == "R2":
     ## If V(i,j-1)≤V(i,j) then R(2) = V(i,j-1)+ V(i,j)   
        R2_con = out_of_bounds_test(location=(location_row, location_col-1), matrix_shape=matrix_shape, test_type="col")
        if np.all(np.isnan(R2_con)):
            return(np.nan)
        else:
            R2_matrix = matrix[location_row, location_col-1]
            regular_value = matrix[location_row, location_col]
            if  R2_matrix <= regular_value:
                return(R2_matrix + regular_value)
            else:
                return(regular_value)
            
    if R_type == "R3":
     ## If V(i+1,j)≤V(i,j) then R(3) = V(i+1,j)+ V(i,j)   
        R3_con = out_of_bounds_test(location=(location_row+1, location_col), matrix_shape=matrix_shape, test_type="row")
        if np.all(np.isnan(R3_con)):
            return(np.nan)
        else:
            R3_matrix = matrix[location_row+1, location_col]
            regular_value = matrix[location_row, location_col]
            if  R3_matrix <= regular_value:
                return(R3_matrix + regular_value)
            else:
                return(regular_value)
            
    if R_type == "R4":
     ## If V(i-1,j)≤V(i,j) then R(4) = V(i-1,j)+ V(i,j)   
        R4_con = out_of_bounds_test(location=(location_row-1, location_col), matrix_shape=matrix_shape, test_type="row")
        if np.all(np.isnan(R4_con)):
            return(np.nan)
        else:
            R4_matrix = matrix[location_row-1, location_col]
            regular_value = matrix[location_row, location_col]
            if  R4_matrix <= regular_value:
                return(R4_matrix + regular_value)
            else:
                return(regular_value)

# %%
def Max_R(location_row, location_col, matrix_shape, matrix):
    """this function get the set R function and return max R

    Args:
        location_row (int): location of row
        location_col (int): location of col
        R_type (string): for example "R1" or "R2" "R3" "R4"
        matrix_shape(tupple): for example (10,10)
        matrix(2d array): The VIIRS rasterize array
    """ 
    R0 = set_R(location_row=location_row, location_col=location_col, R_type="R0", matrix_shape=matrix_shape, matrix=matrix)
    R1 = set_R(location_row=location_row, location_col=location_col, R_type="R1", matrix_shape=matrix_shape, matrix=matrix)
    R2 = set_R(location_row=location_row, location_col=location_col, R_type="R2", matrix_shape=matrix_shape, matrix=matrix)
    R3 = set_R(location_row=location_row, location_col=location_col, R_type="R3", matrix_shape=matrix_shape, matrix=matrix)
    R4 = set_R(location_row=location_row, location_col=location_col, R_type="R4", matrix_shape=matrix_shape, matrix=matrix)
    
    R = [R0, R1, R2, R3, R4]
    return(np.nanmax(R))

# %%
def update_W_matrix(original_matrix, threshold):
    """This function combines the functions above and produce an updated matrix

    Args:
        original_matrix (2d array): array
        threshold (int): an int for threshold
    """
    W = original_matrix.copy() # Set a copy of the original matrix
    matrix_shape = original_matrix.shape ## get the shape
    
    i_list, j_list = np.where(original_matrix>0) ## get where there are fire pixels
    
    i_j_list = create_i_j(i_list, j_list) ## create a list of these locations
    
    for loc in i_j_list: ## for each location
        loc_row = loc[0] ## take the row location
        loc_col = loc[1] ## take col location
        if W[loc_row, loc_col] < threshold:
        ## calculate max R
            R = Max_R(location_row = loc_row, location_col = loc_col, matrix_shape = matrix_shape, matrix = original_matrix)
            W[loc_row, loc_col] = R ## apply the new value
        
    W[W<threshold] = 0 ## everything below the threshold set to zero
        
    return(W)

# %%
def find_FN(corrected_rasterized_VIIRS, GOES_matrix, window_size, fire_class):
    """This function calculate the number of FN (false negative or omission) in a single GOES and VIIRS image. It returns a
        FN matrix where 1 means false negative pixel and zero means not false negative

    Args:
        corrected_rasterized_VIIRS (array): corrected rasterized VIIRS array
        GOES_matrix (array): GOES matrix
        window_size (int): a number that is the size of the window for example 1 for 3X3, 2 for 5X5 and 3 for 7X7 and 0 for 1x1
        fire_class (int): The GOES fire confidence. e.g. 1 is for saturated class, 2 is for high, 3 is for medium
                            4 is for low, and 5 is for all RF for random forest
    """
    ## Fire catagories
    if fire_class == 1:
        GOES_fire_label = [10, 11, 30, 31]
    elif fire_class == 2:
        GOES_fire_label = [13, 33]
    elif fire_class == 3:
        GOES_fire_label = [14,34]
    elif fire_class == 4:
        GOES_fire_label = [15,35]
    elif fire_class == 5:
        GOES_fire_label = [10,11,12,13,14,15,30,31,32,33,34,35]
    elif fire_class == "RF":
        GOES_fire_label = 1
    
    FN = np.zeros((GOES_matrix.shape)) ## set an empty matrix with GOES matrix shape
    
    ## VIIRS fire locations
    VIIRS_fire_i, VIIRS_fire_j = np.where(corrected_rasterized_VIIRS>0)
    ## For each VIIRS fire location:
    
    ## If the window size is 1X1 (no neighbors)
    if window_size == 0:
        for i_loc,j_loc in zip(VIIRS_fire_i, VIIRS_fire_j): ## for each fire pixel
            GOES_value_in_VIIRS_loc = GOES_matrix[i_loc, j_loc] ## take GOES value in VIIRS fire location
            con_GOES_value = np.isin(GOES_value_in_VIIRS_loc, GOES_fire_label) ## if this value is not inside the fire class
            if con_GOES_value == False:
                FN[i_loc, j_loc] = 1 ## Mark it as false negative
        return(FN) ## Return the matrix
                
    else: ## If the window size is bigger than 1X1    
        for i_loc,j_loc in zip(VIIRS_fire_i, VIIRS_fire_j):
        
       ## take the GOES neighbors in a window from a VIIRS fire location
            VIIRS_GOES_fire_neighbor = cell_neighbors(arr=GOES_matrix, i=i_loc, j=j_loc, d=window_size)
        ## if there isn't a  GOES fire pixel near VIIRS fire location
            con_fire = np.any(np.isin(VIIRS_GOES_fire_neighbor, GOES_fire_label))
        
            if con_fire == False:
                FN[i_loc, j_loc] = 1 ## Mark it as false negative
            
        return(FN) ## Return the matrix

# %%
def out_of_bounds_test_value_return(location, matrix_shape, test_type, matrix):
    """get a location and a matrix shape and determine if it is out of bounds, and return nan or the value of the location

    Args:
        location (list): location list for example [2,4] ([row,col])
        matrix_shape (tupple): matrix shape for example (10,10)
        test_type (string): write if it's for row or col. e.g. "row", "col"
        matrix (array): The array we want to get the value
    """
    min_row = 0 ## min row
    min_col = 0 ## min col
    max_row = matrix_shape[0] - 1 ## take the max row and col
    max_col = matrix_shape[1] - 1 
    loc_row = location[0] ## take the row and col location
    loc_col = location[1]
    
    if test_type == "row": ## if the test is for row
        if (loc_row < min_row) | (loc_row > max_row): ## if it is out of bounds 
            return(np.nan) ## return nan
        else:
            return(matrix[loc_row, loc_col]) ## else return the location
    if test_type == "col":
        if (loc_col < min_col) | (loc_col > max_col):
            return(np.nan)
        else:
            return(matrix[loc_row, loc_col])

# %%
def post_FN(corrected_rasterized_VIIRS, FN_matrix, threshold):
    """This function calculate the number of FN (false negative or omission) in a single GOES and VIIRS image after pre-process. 
        It returns FN matrix where 1 means false negative pixel and zero means not false negative and also a df with FN summary

    Args:
        corrected_rasterized_VIIRS (array): corrected rasterized VIIRS array
        FN_matrix (array): pre-processed FN matrix
        threshold (int): threshold of VIIRS fire pixels inside a GOES pixel
        
    """
    ## A new FN matrix in the size of the FN matrix
    ## FNnew = zeros(size(FN))
    FN_new = np.zeros(FN_matrix.shape)
    
    ## VIIRS fire locations
    ## For each location (i,j) not at the boundary
    VIIRS_fire_i, VIIRS_fire_j = np.where(corrected_rasterized_VIIRS>0)
    
    ## For each fire location
    for i_loc, j_loc in zip(VIIRS_fire_i, VIIRS_fire_j):
        
        ## If the value of the corrected VIIRS is higher than 2*threshold
        ## If W(i,j) >2t 
        if corrected_rasterized_VIIRS[i_loc, j_loc] > (2 * threshold):
            ## the new FN gets the value of the pre-process FN
            ## FNnew(i,j) = FN(i,j)
            FN_new[i_loc, j_loc] = FN_matrix[i_loc, j_loc]
         ## If the VIIRS fire value is equal or bigger from the t and also  smaller form 2*t 
         ## If W(i,,j) >= t and W(i,j) <= 2t  
        elif (corrected_rasterized_VIIRS[i_loc, j_loc] >= threshold) & (corrected_rasterized_VIIRS[i_loc, j_loc] <= (threshold * 2)):
            ## Cn = [W(i,j)==W(i+1,j),  W(i,j)==W(i-1,j), W(i,j)==W(i,j+1), W(i,j)==W(i,j-1)] 
            Cn = [] ## empty list for summing up neighbors
            W_i_j = corrected_rasterized_VIIRS[i_loc, j_loc] ## take the location of the VIIRS fire
            ## take the W(i+1,j)
            W_plus_row = out_of_bounds_test_value_return(location=[i_loc+1,j_loc], test_type="row",
                                            matrix_shape=FN_new.shape, matrix=corrected_rasterized_VIIRS)
            ## take the W(i-1,j)
            W_minus_row = out_of_bounds_test_value_return(location=[i_loc-1,j_loc], test_type="row",
                                            matrix_shape=FN_new.shape, matrix=corrected_rasterized_VIIRS)
            ## take the W(i,j+1)
            W_plus_col = out_of_bounds_test_value_return(location=[i_loc,j_loc+1], test_type="col",
                                            matrix_shape=FN_new.shape, matrix=corrected_rasterized_VIIRS)
            ## take the W(i,j-1)
            W_minus_col = out_of_bounds_test_value_return(location=[i_loc,j_loc-1], test_type="col",
                                            matrix_shape=FN_new.shape, matrix=corrected_rasterized_VIIRS)
            ## Append all
            if W_i_j == W_plus_row:
                Cn.append(1)
            else:
                Cn.append(0)
            if W_i_j == W_minus_row:
                Cn.append(1)
            else:
                Cn.append(0)
            if W_i_j == W_plus_col:
                Cn.append(1)
            else:
                Cn.append(0)
            if W_i_j == W_minus_col:
                Cn.append(1)
            else:
                Cn.append(0)
            
            #Fneighbor = [FN(i+1,j), FN(i-1,j), FN(i,j+1), FN(i,j-1)]
            Fneighbor = []
            ## take the F(i+1,j)
            F_plus_row = out_of_bounds_test_value_return(location=[i_loc+1,j_loc], test_type="row",
                                            matrix_shape=FN_new.shape, matrix=FN_matrix)
            ## take the F(i-1,j)
            F_minus_row = out_of_bounds_test_value_return(location=[i_loc-1,j_loc], test_type="row",
                                            matrix_shape=FN_new.shape, matrix=FN_matrix)
            ## take the F(i,j+1)
            F_plus_col = out_of_bounds_test_value_return(location=[i_loc,j_loc+1], test_type="col",
                                            matrix_shape=FN_new.shape, matrix=FN_matrix)
            ## take the F(i,j-1)
            F_minus_col = out_of_bounds_test_value_return(location=[i_loc,j_loc-1], test_type="col",
                                            matrix_shape=FN_new.shape, matrix=FN_matrix)
            
            Fneighbor.append(F_plus_row)
            Fneighbor.append(F_minus_row)
            Fneighbor.append(F_plus_col)
            Fneighbor.append(F_minus_col)
            ## If sum(Cn) = 0  
            if sum(Cn) == 0:
                FN_new[i_loc, j_loc] = FN_matrix[i_loc, j_loc]
            
            else:
                if FN_matrix[i_loc, j_loc] == 1:
                    ## FNnew(i,j) = ½ * sum( Cn * Fneighbor)
                    Cn_array = np.array(Cn)
                    Fneighbor_array = np.array(Fneighbor) 
                    FN_new[i_loc, j_loc] = 0.5 * np.nansum(Cn_array * Fneighbor_array)
                    
    return(FN_new)

# %%
def false_alarm_dummy_rast_non_max_sup(dummy_rast, image_path):
    """Gets a dummy rasterized file and return a false alarm gdf
    Keyword arguments

    Args:
        dummy_rast (array): raster dummy rast
        image_path (string): for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif'
    """
    n = np.sum(np.isnan(dummy_rast) == False) ## total number of pixels that are not nan
    
    ## Calcualte statistics
    
    ON = np.sum(dummy_rast == 0)
    
    OF =  np.sum(dummy_rast[(dummy_rast<88) & (dummy_rast>0)])
    
    CN = np.sum(dummy_rast == 88)
    
    CF = np.sum(dummy_rast == 99)
    
    prop_CE = CN/n
    
    prop_OE = OF/n
    
    CER = CN/(CN+CF)
    
    OER = OF/(OF+ON)
    
    ## Get fire name and file name ##
    
    fire_name = image_path.split("\\")[-3]
    
    f_name = image_path.split("\\")[-1]
    
    if (f_name.endswith(".tif")): ## get the date if it is nc or tif
        
        d_name = get_time_from_goes(file_name=image_path, form="string", VIIRS_time="N",  NetCDF="N")
    
    if (f_name.endswith(".nc")):
        
        d_name = get_time_from_goes(file_name=image_path, form="string", VIIRS_time="N",  NetCDF="Y")
    
    # create dict
    d = {"fire_name":[fire_name],"file_name":[f_name], "date":[d_name], "number_of_pixels_(n)":[n], "ON":[ON],
        "OF":[OF], "CN":[CN], "CF":[CF], "prop_comission_error":[prop_CE],
        "prop_omission_error":[prop_OE], "comission_error_rate":[CER], "omission_error_rate":[OER]}
    
    # create df
    df = pd.DataFrame(data=d)
    
    return(df)

# %%
def calculate_false_alarm_non_max_sup(corrected_rasterized_VIIRS, distance, image_path,
                                      fire_class_GOES, fire_class_FN, threshold, return_image):
    """Gets a rasterized file and return a false alarm df 

    Args:
        corrected_rasterized_VIIRS (array): rasterized VIIRS file
        distance (int): the steps from a fire pixel for calculating neighbors
        image_path (string): for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif' or NC
        fire_class_GOES (string): The GOES fire confidence. e.g. 1 is for saturated class, 2 is for high, 3 is for medium
                            4 is for low, and 5 is for all
        fire_class_FN (string): The GOES fire confidence. e.g. 1 is for saturated class, 2 is for high, 3 is for medium
                            4 is for low, and 5 is for all for the FN calculation!
        threshold (int): Number of VIIRS pixels threshold
        return_image (string): "Y" to return and "N" to return only df
    """
    ## Fire catagories
    if fire_class_GOES == 1:
        GOES_fire_label = [10, 11, 30, 31]
    elif fire_class_GOES == 2:
        GOES_fire_label = [13, 33]
    elif fire_class_GOES == 3:
        GOES_fire_label = [14,34]
    elif fire_class_GOES == 4:
        GOES_fire_label = [15,35]
    elif fire_class_GOES == 5:
        GOES_fire_label = [10,11,12,13,14,15,30,31,32,33,34,35]
    elif fire_class_GOES == "RF":
        GOES_fire_label = 1
        
        
    if distance == 1: ## set the size of the window 3x3
        kernel = 9
    
    if distance == 2: ## set the size of the window 5x5
        kernel = 25
        
    if distance == 3: ## set the size of the window 7x7
        kernel = 49
        
    f_name = image_path.split("\\")[-1] ## take file name
    
    if f_name.endswith(".tif"):
    
        src_fp = rioxarray.open_rasterio(image_path) ## open it with rasterio
    
        fp_rast = src_fp.values[0].copy() ## take the image 
        
        #process_con = np.any(np.isin(fp_rast,GOES_fire_label)) ## this is only if there is this pixel group in the raster
        
        #if process_con == True: ## if it is true then:
        
        dummy_rast = fp_rast.copy() ## set a copy of the fire product
            
    
        dummy_rast[np.isnan(dummy_rast) == False] = 0 ## set all the values in the dummy rast that are not nan as 0 (non-fire)
            
            ## Get a False negative matrix
        FN_matrix = find_FN(corrected_rasterized_VIIRS=corrected_rasterized_VIIRS,
                                GOES_matrix=fp_rast, window_size=distance, fire_class=fire_class_FN)
            ## correct it
        FN_corrected_matrix = post_FN(corrected_rasterized_VIIRS=corrected_rasterized_VIIRS, FN_matrix=FN_matrix,
                                          threshold=threshold)
            ## Take all the values of the FN corrected matrix and insert them to the dummy rast
        dummy_rast[FN_corrected_matrix != 0] = FN_corrected_matrix[FN_corrected_matrix!=0]
            
            ### GOES PART ########
        
        GOES_row_i, GOES_col_j = np.where(np.isin(fp_rast,GOES_fire_label)) ## take the row and col index for all GOES pixels that are fire
        ## in a given pixels group, for example g1, g2 etc etc...
        
        for i, j in zip(GOES_row_i, GOES_col_j): ## for row index i and col index j
                
                ## take the VIIRS neighbors in the given distance
            neighboors = cell_neighbors(arr=corrected_rasterized_VIIRS, i = i, j= j, d=distance) 
                
                ## calculate the sum of all the zero and nan of the neighbors
            calc = np.all(neighboors == 0) 
                
            if (np.any(neighboors>0)) == True: ## if any of the neighbors value is heigher than 0
                dummy_rast[i,j] = 99 ## set this pixel as fire (99 means fire)
            if calc == True: ## if all of the VIIRS srounding neighbors dont have fire
                dummy_rast[i,j] = 88 ## set it to 88 - false alarm!
        
        ## calculate false alarm
        df = false_alarm_dummy_rast_non_max_sup(dummy_rast=dummy_rast, image_path=image_path)
        
        dummy_rast[np.isnan(fp_rast)] = np.nan ## set nan values according to fp raster
                
        if (return_image == "Y"): ## in case we want to also see the image to analyze it!
            return(df, dummy_rast)
        if (return_image == "N"):
            return(df)
            

# %%
def calculate_false_alarm_nn_non_max_sup(corrected_rasterized_VIIRS, image_path,
                                      fire_class_GOES, fire_class_FN, threshold, return_image):
    """Gets a rasterized file and return a false alarm df with no neighbors

    Args:
        corrected_rasterized_VIIRS (array): rasterized VIIRS file
        image_path (string): for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif' or NC
        fire_class_GOES (string): The GOES fire confidence. e.g. 1 is for saturated class, 2 is for high, 3 is for medium
                            4 is for low, and 5 is for all
        fire_class_FN (string): The GOES fire confidence. e.g. 1 is for saturated class, 2 is for high, 3 is for medium
                            4 is for low, and 5 is for all for the FN calculation!
        threshold (int): Number of VIIRS pixels threshold
        return_image (string): "Y" to return and "N" to return only df
    """
    ## Fire catagories
    if fire_class_GOES == 1:
        GOES_fire_label = [10, 11, 30, 31]
    elif fire_class_GOES == 2:
        GOES_fire_label = [13, 33]
    elif fire_class_GOES == 3:
        GOES_fire_label = [14,34]
    elif fire_class_GOES == 4:
        GOES_fire_label = [15,35]
    elif fire_class_GOES == 5:
        GOES_fire_label = [10,11,12,13,14,15,30,31,32,33,34,35]
    elif fire_class_GOES == "RF":
        GOES_fire_label = 1
        
     
    f_name = image_path.split("\\")[-1] ## take file name
    
    if f_name.endswith(".tif"):
    
        src_fp = rioxarray.open_rasterio(image_path) ## open it with rasterio
    
        fp_rast = src_fp.values[0].copy() ## take the image 
        
        #process_con = np.any(np.isin(fp_rast,GOES_fire_label)) ## this is only if there is this pixel group in the raster
        
        #if process_con == True: ## if it is true then:
        
        dummy_rast = fp_rast.copy() ## set a copy of the fire product
            
    
        dummy_rast[np.isnan(dummy_rast) == False] = 0 ## set all the values in the dummy rast that are not nan as 0 (non-fire)
            
            ## Get a False negative matrix
        FN_matrix = find_FN(corrected_rasterized_VIIRS=corrected_rasterized_VIIRS,
                                GOES_matrix=fp_rast, window_size=0, fire_class=fire_class_FN)
            ## correct it
        FN_corrected_matrix = post_FN(corrected_rasterized_VIIRS=corrected_rasterized_VIIRS, FN_matrix=FN_matrix,
                                          threshold=threshold)
            ## Take all the values of the FN corrected matrix and insert them to the dummy rast
        dummy_rast[FN_corrected_matrix != 0] = FN_corrected_matrix[FN_corrected_matrix!=0]
            
            ### GOES PART ########
        
        GOES_row_i, GOES_col_j = np.where(np.isin(fp_rast,GOES_fire_label)) ## take the row and col index for all GOES pixels that are fire
        ## in a given pixels group, for example g1, g2 etc etc...
        
        for i, j in zip(GOES_row_i, GOES_col_j): ## for row index i and col index j
                
                ## take the VIIRS value in GOES fire locations
            VIIRS_value = corrected_rasterized_VIIRS[i,j]
                
            if (VIIRS_value>0) == True: ## if any of the neighbors value is heigher than 0
                dummy_rast[i,j] = 99 ## set this pixel as fire (99 means fire)
            else: ## if all of the VIIRS srounding neighbors dont have fire
                dummy_rast[i,j] = 88 ## set it to 88 - false alarm!
        
        ## calculate false alarm
        df = false_alarm_dummy_rast_non_max_sup(dummy_rast=dummy_rast, image_path=image_path)
        
        dummy_rast[np.isnan(fp_rast)] = np.nan ## set nan values according to fp raster
                
        if (return_image == "Y"): ## in case we want to also see the image to analyze it!
            return(df, dummy_rast)
        if (return_image == "N"):
            return(df)
            

# %%
def compare_GOES_VIIRS_image_time(GOES_fp_path, VIIRS_points, time_delta):
    """
    Get a GOES fire product path and a VIIRS unique time and return the GOES file name and VIIRS file name

    :GOES_fp_path: image time string for example 'F:\\Project_data\\goes_data_creek\\fire_product_out\\2020_09_05_09_02.tif'
    :VIIRS_points: VIIRS points file from "generate_VIIRS_points" function
    
    :return: return the GOES file name and VIIRS file name list
    """ 
    ## We set the start and end time
    time_diff = time_delta
    start_time = get_time_from_goes(GOES_fp_path, form="time", VIIRS_time="N", NetCDF="N") - timedelta(minutes = time_diff)
    end_time = get_time_from_goes(GOES_fp_path, form="time", VIIRS_time="N", NetCDF="N") + timedelta(minutes = time_diff)
    
    ## Get VIIRS time
    VIIRS_file_list = [] ## set empty VIIRS list
    VIIRS_GEO_list = [] ## empty list for GEO files
    GOES_list = []
    df_grouped = VIIRS_points.groupby(["DATE_TIME"]).first().reset_index() ## group by the unique time of the VIIRS image
    for i in range(len(df_grouped)): ## for every row
        d_t = df_grouped["DATE_TIME"].iloc[i] ## take the date_time of VIIRS
        date_time = get_time_from_goes(d_t, form="time", VIIRS_time="Y", NetCDF="N") ## convert into date format
        con = time_in_range(start_time, end_time, date_time) ## see if it is inside GOES time
        if con == True: ## If it is true
            VIIRS_file_name = df_grouped["Fire_file"].iloc[i] ## take the VIIRS file name
            VIIRS_file_list.append(VIIRS_file_name) ## append it to a list
            VIIRS_GEO_name = df_grouped["GEO_file"].iloc[i] ## take the VIIRS file name
            VIIRS_GEO_list.append(VIIRS_GEO_name) ## append it to a list
            GOES_list.append(GOES_fp_path)

    d = {"GOES_file":GOES_list, "VIIRS_fire_file":VIIRS_file_list, "VIIRS_GEO_file":VIIRS_GEO_list}
    df = pd.DataFrame(d)
    return(df)

# %%
def get_VIIRS_bounds_polygon(VIIRS_fire_path):
    """
    Get a VIIRS fire product path return a polygon of VIIRS image bounds

    :VIIRS_fire_path: VIIRS fire product for example "C:\\Asaf\\VIIRS\\VNP14IMG.A2020251.2042.001.2020258064138"
    
    :return: geodataframe of bounds in lat/lon WGS84 projection
    """ 
    VIIRS = xr.open_dataset(VIIRS_fire_path) ## Open file
    ## Now for the bbox bounds
    lat_bounds = VIIRS.attrs['GRingPointLatitude']
    lon_bounds = VIIRS.attrs['GRingPointLongitude']

    pol_bounds = [] ## Empty bounds list
    for i in range(len(lat_bounds)): ## for all bounds
        coord = (lon_bounds[i], lat_bounds[i]) ## Take lon/lat
        pol_bounds.append(coord) ## append
    pol_bounds.append((lon_bounds[0], lat_bounds[0])) ## append firt coords to "close" the polygon

    pol = shapely.Polygon(pol_bounds) ## Create polygon
    gdf_pol = gpd.GeoDataFrame({'geometry': [pol]}, crs=4326) ## Create as gdf

    return(gdf_pol)

# %%
def clip_GOES_to_VIIRS_bounds(GOES_path, VIIRS_image_path, GOES_CRS):
    """
    Get a VIIRS fire product path and a GOES fire product path and CRS and return a cliped GOES image
    
    :GOES_path: GOES fire product for example "F:\\Project_data2\\goes_data_false_alarm_cal\\fire_product_out\\2022_01_09_00_02.tif"
    :VIIRS_image_path: VIIRS fire product for example "C:\\Asaf\\VIIRS\\VNP14IMG.A2020251.2042.001.2020258064138"
    :GOES_CRS: GOES CRS code that can be extracted from an MCMI file
    
    :return: a clipped GOES fire product image
    """ 
    ## Open data
    CRS = GOES_CRS
    VIIRS_bounds_polygon = get_VIIRS_bounds_polygon(VIIRS_fire_path=VIIRS_image_path)
    GOES = rioxarray.open_rasterio(GOES_path)
    ## set CRS for GOES image and the polygon
    GOES = GOES.rio.set_crs(CRS)
    VIIRS_bounds_polygon = VIIRS_bounds_polygon.to_crs(CRS)
    ## Clip
    GOES_clip = GOES.rio.clip(VIIRS_bounds_polygon["geometry"])
    ## Convert -99 to nan
    con = GOES_clip.values[0] == -99
    GOES_clip.values[0][con] = np.nan
    return(GOES_clip)

# %%
def generate_VIIRS_time_from_image(VIIRS_fire_path, type_of_date:str) -> str:
    """
    Get a VIIRS fire product and return the date and time of the file

    :VIIRS_fire_path: VIIRS fire product for example 'C:\\Asaf\\VIIRS\\VNP14IMG.A2020251.2042.001.2020258064138.nc'
    :type_of_date: A string "date_time" for full date time, "date" just for date, "time" just for time
    :return: a string of Year-Month-Day Hour-Minute
    """ 
    VIIRS_fire = xr.open_dataset(VIIRS_fire_path)
    date_time = VIIRS_fire.attrs['StartTime']
    date = date_time.split(" ")[0]
    time_raw = date_time.split(" ")[1]
    time_raw_sub = time_raw.split(".")[0]
    hour = time_raw_sub.split(":")[0]
    minute = time_raw_sub.split(":")[1]
    correct_date_time = "{} {}:{}".format(date, hour, minute)
    only_time = "{}:{}".format(hour,minute)

    if type_of_date == "date_time":
        return(correct_date_time)
    elif type_of_date == "date":
        return(date)
    elif type_of_date == "time":
        return(only_time)


