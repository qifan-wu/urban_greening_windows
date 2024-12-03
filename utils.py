import rasterio
import xarray as xr
from rasterio.io import MemoryFile

# manually calculate the transform for nee file from flux.org
minx, miny, maxx, maxy = -180.0, -90.0, 180.0, 90.0
resolution_x = 360 / 4320
resolution_y = 180 / 2160
nee_transform = rasterio.transform.from_origin(minx, maxy, resolution_x, resolution_y)


# Write the ndarray to a GeoTIFF
def save_tiff(data, output_path, crs, transform):
    with rasterio.open(
        output_path,
        "w",
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=1,  # Number of bands
        dtype=data.dtype,
        crs=crs,
        transform=transform,
        nodata=np.nan, 
    ) as dst:
        # Write data to the first band
        dst.write(data, 1)

    print(f"GeoTIFF saved to {output_path}")

def create_in_memory_ds(data, crs, transform, return_file=False):
    '''
    data: 2d array
    '''
    memfile = MemoryFile()
    with memfile.open(
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype=data.dtype.name,
        crs=crs,
        transform=transform,
        nodata=np.nan
    ) as dst:
        # Write the data to the in-memory file
        dst.write(data, 1)
    if return_file:
        return memfile
    return dst

def read_nee(nee_file, nee_transform, nee_memory):

    nee_ds = xr.open_dataset(nee_file)
    time_length = len(nee_ds['NEE']) # should be 12
    for i in range(time_length - 10): #TODO: remove -10 after testing
        time_str = nee_ds.time[i].dt.strftime("%Y-%m-%d").values
        print("Writing NEE for time: ", time_str), "to memory"
        nee_time = nee_ds['NEE'][i]
        nee_crs = nee_time.rio.crs
        # Save to an in-memory GeoTIFF
        # memfile = MemoryFile()
        # with memfile.open(
        #     driver="GTiff",
        #     height=nee_time.shape[0],
        #     width=nee_time.shape[1],
        #     count=1,
        #     dtype=nee_time.dtype.name,
        #     crs=nee_crs,
        #     transform=nee_transform,
        # ) as dst:
        #     # Write the data to the in-memory file
        #     dst.write(nee_time.values, 1)
        memfile = create_in_memory_ds(nee_time, nee_crs, nee_transform, return_file=True)
        nee_memory.append(memfile)


from shapely.geometry import mapping
from rasterio.mask import mask
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject
import numpy as np

def pipe_read_gen_params(msa, gpp_file, nlcd_file, ua_file, nee_memory):
    with rasterio.open(nlcd_file) as nlcd_dstrd:
        nlcd_crs = nlcd_dstrd.crs
        geometries_aea = [mapping(geom) for geom in msa.geometry]

         # Clip the raster using the reprojected geometries
        try:
            nlcd_clip_image, nlcd_clip_transform = mask(nlcd_dstrd, geometries_aea, crop=True)
        except ValueError as e:
            print(f"No NLCD data for {msa['NAMELSAD'].values[0]}")
            nlcd_clip_image
            
        nlcd_msa = nlcd_clip_image[0]

        with rasterio.open(gpp_file) as gpp_dstrd:
            msa_crsgpp=msa.to_crs(gpp_dstrd.crs)
            msa_crs = msa_crsgpp.crs
            gpp_crs = gpp_dstrd.crs
            
            gpp_clip_image, gpp_clip_transform = mask(gpp_dstrd,msa_crsgpp.geometry,crop=True)
            

            # replace fillvalue 65535 with np.nan to avoid future calculation error
            gpp_clip_image = np.where(gpp_clip_image == 65535, np.nan, gpp_clip_image) 

            gpp_msa = gpp_clip_image[0]

            gpp_msa_rr = np.empty((nlcd_clip_image.shape[1], nlcd_clip_image.shape[2]), dtype=np.float32) #make it support nan

            reproject(
                source=gpp_msa,
                destination=gpp_msa_rr,
                src_transform=gpp_clip_transform,
                src_crs=gpp_crs,
                dst_transform=nlcd_clip_transform,
                dst_crs=nlcd_crs,
                resampling=Resampling.nearest, # use nearest to prevent unreasonable upscaling
                src_nodata=np.nan
            )
            

            with rasterio.open(ua_file) as ua_dstrd:
                ua_crs = ua_dstrd.crs

                ua_clip_image, ua_clip_transform = mask(ua_dstrd, geometries_aea, crop=True)
                ua_msa = ua_clip_image[0]

                ua_msa_rr = np.empty((nlcd_clip_image.shape[1], nlcd_clip_image.shape[2]), dtype=ua_dstrd.meta['dtype'])

                reproject(
                    source=ua_msa,
                    destination=ua_msa_rr,
                    src_transform=ua_clip_transform,
                    src_crs=ua_crs,
                    dst_transform=nlcd_clip_transform,
                    dst_crs=nlcd_crs,
                    resampling=Resampling.nearest
                )

            # TODO: replace nee_memory[0] with specific time
            # nee_dstrd = nee_memory[0]

            # nee_crs = gpp_crs
            
            # # clip nee raster with shape
            # nee_clip_image, nee_clip_transform = mask(nee_dstrd, [mapping(geom) for geom in msa_crsgpp.geometry], crop=True, all_touched=True)  # Include all touched pixels)
            # nee_clip_image = np.where(nee_clip_image == -9999, np.nan, nee_clip_image) 
            
            # nee_msa = nee_clip_image[0]

            memfile_nee = nee_memory[0]
            with memfile_nee.open() as nee_dstrd:
                nee_crs = gpp_crs
                
                # clip nee raster with shape
                nee_clip_image, nee_clip_transform = mask(nee_dstrd, [mapping(geom) for geom in msa_crsgpp.geometry], crop=True, all_touched=True)  # Include all touched pixels)
                nee_clip_image = np.where(nee_clip_image == -9999, np.nan, nee_clip_image) 
                
                nee_msa = nee_clip_image[0]

            # memfile.close()

    return {'gpp_msa_rr': gpp_msa_rr, 
            'ua_msa_rr': ua_msa_rr, 
            'nlcd_msa': nlcd_msa, 
            'nlcd_crs': nlcd_crs, 
            'nlcd_clip_transform': nlcd_clip_transform,
            'nee_msa': nee_msa, 
            'nee_clip_transform': nee_clip_transform, 
            'nee_crs': nee_crs}


def create_mask(gpp_msa_rr, ua_msa_rr, nlcd_msa):
    
    urban_mask = (ua_msa_rr != 0)
    suburban_mask = (ua_msa_rr == 0)
    forest_mask = (nlcd_msa == 41) | (nlcd_msa == 42) | (nlcd_msa == 43)
    shrub_mask = (nlcd_msa == 51) | (nlcd_msa == 52)
    grass_mask = (nlcd_msa > 70) & (nlcd_msa < 90)
    wetland_mask = (nlcd_msa >= 90) & (nlcd_msa < 100)

    nlcd_mask_dict = {
        'urban_forest': urban_mask & forest_mask,
        'suburban_forest': suburban_mask & forest_mask,
        'urban_shrub': urban_mask & shrub_mask,
        'suburban_shrub': suburban_mask & shrub_mask,
        'urban_grass': urban_mask & grass_mask,
        'suburban_grass': suburban_mask & grass_mask,
        'urban_wetland': urban_mask & wetland_mask,
        'suburban_wetland': suburban_mask & wetland_mask,
    }

    return nlcd_mask_dict


def gap_fill_gpp(gpp_msa_rr, ua_msa_rr, nlcd_msa, msa_name):
    valid_gpp_mask = ~np.isnan(gpp_msa_rr)
    nlcd_mask_dict = create_mask(gpp_msa_rr, ua_msa_rr, nlcd_msa)

    gpp_msa_rr_filled = np.copy(gpp_msa_rr)

    record = {'name': msa_name}

    water_mask = (nlcd_msa > 10) & (nlcd_msa < 20)
    developed_mask = (nlcd_msa > 20) & (nlcd_msa < 40)

    gpp_msa_rr_filled[~valid_gpp_mask & water_mask] = 0
    gpp_msa_rr_filled[~valid_gpp_mask & developed_mask] = 0
                    
    for category, landcover_mask in nlcd_mask_dict.items():

        category_gpp_mean = np.nanmean(gpp_msa_rr[landcover_mask])
        record[category] = category_gpp_mean

        # # debug
        # print(category)
        # import rasterio.plot
        # rasterio.plot.show(landcover_mask) 
        # print('----')

        gpp_msa_rr_filled[~valid_gpp_mask & landcover_mask] = category_gpp_mean

    gpp_mean_cat_data.append(record)

    print("     ", record)
    # print("Gap filled GPP values for each landcover category have been calculated")

    return gpp_msa_rr_filled


def get_gpp_coarse(gpp_msa_rr_filled, nee_msa, nlcd_clip_transform, nlcd_crs, nee_clip_transform, nee_crs):
    gpp_filled_coarse = np.empty_like(nee_msa, dtype=np.float32)

    # Reproject and resample GPP to match NEE's CRS, resolution, and extent
    reproject(
        source=gpp_msa_rr_filled, 
        destination=gpp_filled_coarse,  
        src_transform=nlcd_clip_transform, 
        src_crs=nlcd_crs,  
        dst_transform=nee_clip_transform, 
        dst_crs=nee_crs,  
        resampling=Resampling.average,  # Use nearest neighbor resampling
        src_nodata=np.nan,  
        dst_nodata=np.nan  
    )
    return gpp_filled_coarse

def get_nee_gpp_ratio_fine(gpp_msa_rr_filled, nee_msa, nlcd_clip_transform, nlcd_crs, nee_clip_transform, nee_crs):
    gpp_filled_coarse = get_gpp_coarse(gpp_msa_rr_filled, nee_msa, nlcd_clip_transform, nlcd_crs, nee_clip_transform, nee_crs)

    
    safe_gpp_filled_coarse = np.where(gpp_filled_coarse == 0, 0.1, gpp_filled_coarse) # Avoid division by zero
    nee_gpp_ratio_coarse = nee_msa / safe_gpp_filled_coarse

    nee_gpp_ratio_fine = np.empty_like(gpp_msa_rr_filled, dtype=np.float32)

    # reproject and resample NEE/GPP ratio to match gpp filled resolution
    reproject(
        source=nee_gpp_ratio_coarse,  
        destination=nee_gpp_ratio_fine,  
        src_transform=nee_clip_transform, 
        src_crs=nee_crs, 
        dst_transform=nlcd_clip_transform, 
        dst_crs=nlcd_crs,  
        resampling=Resampling.nearest,  # Use nearest neighbor resampling
        src_nodata=np.nan,  
        dst_nodata=np.nan  
    )
    return nee_gpp_ratio_fine


def get_downscaled_nee_msa(msa_ds, msa, gpp_file, nlcd_file, ua_file, nee_memory):
     # Subsetting to my AOI
    msa_name = msa['NAMELSAD'].values[0]
    pipe_output = pipe_read_gen_params(msa, gpp_file, nlcd_file, ua_file, nee_memory)
    gpp_msa_rr = pipe_output['gpp_msa_rr']
    ua_msa_rr = pipe_output['ua_msa_rr']
    nlcd_msa = pipe_output['nlcd_msa']
    nlcd_crs = pipe_output['nlcd_crs']
    nlcd_clip_transform = pipe_output['nlcd_clip_transform']
    nee_msa = pipe_output['nee_msa']
    nee_clip_transform = pipe_output['nee_clip_transform']
    nee_crs = pipe_output['nee_crs']

    gpp_msa_rr_filled = gap_fill_gpp(gpp_msa_rr, ua_msa_rr, nlcd_msa, msa_name)

    nee_gpp_ratio_fine = get_nee_gpp_ratio_fine(gpp_msa_rr_filled, nee_msa, nlcd_clip_transform, nlcd_crs, nee_clip_transform, nee_crs)
    downscaled_nee = nee_gpp_ratio_fine * gpp_msa_rr_filled
    # downscaled_nee[np.isnan(nee_gpp_ratio_fine) | np.isnan(gpp_msa_rr_filled)] = 0 #if nee or gpp does not have data, set downscaled nee to 0

    downscaled_nee_info = {
        'data': downscaled_nee,
        'crs': nlcd_crs,
        'transform': nlcd_clip_transform
    }
    return downscaled_nee_info

from rasterio.merge import merge
import rasterio
from rasterio.io import MemoryFile

def merge_in_batches(memfiles, batch_size):
    merged_files = []
    for i in range(0, len(memfiles), batch_size):
        batch = memfiles[i:i + batch_size]
        datasets = [mem.open() for mem in batch]
        merged_data, merged_transform = merge(datasets)
        for ds in datasets:
            ds.close()
        # Write intermediate result to a temporary MemoryFile
        temp_mem = MemoryFile()
        with temp_mem.open(
            driver="GTiff",
            height=merged_data.shape[1],
            width=merged_data.shape[2],
            count=1,
            dtype=merged_data.dtype,
            crs=datasets[0].crs,  # Assume same CRS
            transform=merged_transform,
        ) as dst:
            dst.write(merged_data[0], 1)
        merged_files.append(temp_mem)
    return merged_files



def main():
    import geopandas as gpd
    nee_memory = []
    nee_file = "../urban_greening/NEE.RS.FP-NONE.MLM-ALL.METEO-NONE.4320_2160.monthly.2015.nc"
    read_nee(nee_file, nee_transform, nee_memory)
    # print(nee_memory)

    # msa_file = '../urban_greening/msa/michiganMSA_reprojected.shp' # for michigan, crs: Albers Equal Area
    # msa_file = '../urban_greening/msa/msaUS/msaUS_aea.shp' # for US, crs: Albers Equal Area
    msa_file = '../urban_greening/msa/msaUS_mainland_aea.shp' # for US without MSA from Hawaii, Puerto Rico and Alaska (no NLCD or no carbon data), crs: Albers Equal Area

    gpp_file = "../urban_greening/nov.15/michigan_test/modis-250-gpp-2015001.tif" # EPSG:4326
    nee_file = "../urban_greening/NEE.RS.FP-NONE.MLM-ALL.METEO-NONE.4320_2160.monthly.2015.nc" # EPSG:4326, resolution 1/12 degree
    nlcd_file = "../urban_greening/nov.15/nlcd_2016_land_cover_l48_20210604.img" # crs: Albers Equal Area, resolution 30m
    ua_file = "../urban_greening/ua/ua_us_30_clip.tif" # crs: Albers Equal Area, resolution 30m

    msa_ds=gpd.read_file(msa_file)

    # from rasterio import plot
    # import matplotlib.pyplot as plt
    

    # msa_names = ['Grand Rapids-Wyoming-Kentwood, MI Metro Area', 'Lansing-East Lansing, MI Metro Area']
    # msas = [msa_ds[msa_ds['NAMELSAD']==msa_name] for msa_name in msa_names]

    import pandas as pd
    global gpp_mean_cat_data 
    gpp_mean_cat_data = []

    mem_downscaled_nee_list = []

    # ======== test with all msas ========
    for index, record in msa_ds.iterrows():
        msa_name = record['NAMELSAD']
        print(f'Generating downsclaed data for {msa_name}...')
        msa = msa_ds.loc[[index]]
        downscaled_nee_msa = get_downscaled_nee_msa(msa_ds, msa, gpp_file, nlcd_file, ua_file, nee_memory)
        mem = create_in_memory_ds(downscaled_nee_msa['data'], downscaled_nee_msa['crs'], downscaled_nee_msa['transform'], return_file=True)
        mem_downscaled_nee_list.append(mem)
        

    # gpp_mean_data_df = pd.DataFrame(gpp_mean_cat_data)
    # gpp_mean_data_df.to_csv('../output/gpp_mean_data2.csv', index=False)

    datasets = [mem.open() for mem in mem_downscaled_nee_list]
    
    # Step 2: Merge datasets
    print("Merging datasets. This might take a while...")
    merged_data, merged_transform = rasterio.merge.merge(datasets, nodata=np.nan)
    
    print("Merging completed")

    merged_raster = merged_data[0] #get the first band
    output_file = '../output/mergedNEE_US.tif'
    save_tiff(merged_raster, output_file, datasets[0].crs, merged_transform)

    # Step 4: Close datasets
    for ds in datasets:
        ds.close()
    # ======== test with all msas ========
    

    # ====== test with 1 msa ========
    # msa_name = 'Grand Rapids-Wyoming-Kentwood, MI Metro Area'
    # msa=msa_ds[msa_ds['NAMELSAD']==msa_name] # Subsetting to my AOI

    # pipe_output = pipe_read_gen_params(msa, gpp_file, nlcd_file, ua_file, nee_memory)
    # gpp_msa_rr = pipe_output['gpp_msa_rr']
    # ua_msa_rr = pipe_output['ua_msa_rr']
    # nlcd_msa = pipe_output['nlcd_msa']
    # nlcd_crs = pipe_output['nlcd_crs']
    # nlcd_clip_transform = pipe_output['nlcd_clip_transform']
    # nee_msa = pipe_output['nee_msa']
    # nee_clip_transform = pipe_output['nee_clip_transform']
    # nee_crs = pipe_output['nee_crs']

    # gpp_msa_rr_filled = gap_fill_gpp(gpp_msa_rr, ua_msa_rr, nlcd_msa, msa_name)
    # # save_tiff(gpp_msa_rr_filled, '../output/msagpptest.tif', nlcd_crs, nlcd_clip_transform)
    
    # nee_gpp_ratio_fine = get_nee_gpp_ratio_fine(gpp_msa_rr_filled, nee_msa, nlcd_clip_transform, nlcd_crs, nee_clip_transform, nee_crs)
    # downscaled_nee = nee_gpp_ratio_fine * gpp_msa_rr_filled
    # save_tiff(downscaled_nee, '../output/msaneetest.tif', nlcd_crs, nlcd_clip_transform)
    # ====== test with 1 msa ========
    


if __name__ == "__main__":
    main()