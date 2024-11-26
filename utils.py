import rasterio
import xarray as xr
from rasterio.io import MemoryFile

# manually calculate the transform for nee file from flux.org
minx, miny, maxx, maxy = -180.0, -90.0, 180.0, 90.0
resolution_x = 360 / 4320
resolution_y = 180 / 2160
nee_transform = rasterio.transform.from_origin(minx, maxy, resolution_x, resolution_y)


# Write the ndarray to a GeoTIFF
def save_tiff(output_file, output_path, crs, transform):
    with rasterio.open(
        output_path,
        "w",
        driver="GTiff",
        height=output_file.shape[0],
        width=output_file.shape[1],
        count=1,  # Number of bands
        dtype=output_file.dtype,
        crs=crs,
        transform=transform,
        # nodata=np.nan, # do not set this or there will be reprojet edge with value 0
    ) as dst:
        # Write data to the first band
        dst.write(output_file, 1)

    print(f"GeoTIFF saved to {output_path}")


def read_nee(nee_file, nee_transform, nee_memory):


    nee_ds = xr.open_dataset(nee_file)
    time_length = len(nee_ds['NEE']) # should be 12
    for i in range(time_length - 10): #TODO: remove -10 after testing
        time_str = nee_ds.time[i].dt.strftime("%Y-%m-%d").values
        print("Writing NEE for time: ", time_str), "to memory"
    # nee_ds['NEE'][0].rio.write_crs('EPSG:4326', inplace=True)
    # nee_ds['NEE'][0].rio.write_transform(nee_transform, inplace=True)
        nee_time = nee_ds['NEE'][i]
        # Save to an in-memory GeoTIFF
        memfile = MemoryFile()
        with memfile.open(
            driver="GTiff",
            height=nee_time.shape[0],
            width=nee_time.shape[1],
            count=1,
            dtype=nee_time.dtype.name,
            crs=nee_time.rio.crs,
            transform=nee_transform,
        ) as dst:
            # Write the data to the in-memory file
            dst.write(nee_time.values, 1)
        nee_memory.append(memfile)


from shapely.geometry import mapping
from rasterio.mask import mask
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject
import numpy as np

def pipe(msa, gpp_file, nlcd_file, ua_file, nee_memory):
    with rasterio.open(nlcd_file) as nlcd_dstrd:
        nlcd_crs = nlcd_dstrd.crs
        geometries_aea = [mapping(geom) for geom in msa.geometry]

         # Clip the raster using the reprojected geometries
        nlcd_clip_image, nlcd_clip_transform = mask(nlcd_dstrd, geometries_aea, crop=True)
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

                ua_msa_rr = np.empty((ua_clip_image.shape[1], ua_clip_image.shape[2]), dtype=ua_dstrd.meta['dtype'])

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
            memfile = nee_memory[0]
            with memfile.open() as nee_dstrd:
                nee_crs = gpp_crs
                
                # clip nee raster with shape
                nee_clip_image, nee_clip_transform = mask(nee_dstrd, [mapping(geom) for geom in msa_crsgpp.geometry], crop=True, all_touched=True)  # Include all touched pixels)
                nee_clip_image = np.where(nee_clip_image == -9999, np.nan, nee_clip_image) 
                
                nee_msa = nee_clip_image[0]

            memfile.close()
    return {'gpp_msa_rr': gpp_msa_rr, 
            'ua_msa_rr': gpp_msa_rr, 
            'nlcd_msa': nlcd_msa, 
            'nlcd_clip_image': nlcd_clip_image, 
            'nlcd_crs': nlcd_crs, 
            'nlcd_clip_transform': nlcd_clip_transform,
            'nee_msa': nee_msa, 
            'nee_clip_transform': nee_clip_transform, 
            'nee_crs': nee_crs}


def create_mask(gpp_msa_rr, ua_msa_rr, nlcd_msa):
    valid_gpp_mask = ~np.isnan(gpp_msa_rr)
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

    return valid_gpp_mask, nlcd_mask_dict


def gap_fill_gpp(gpp_msa_rr, ua_msa_rr, nlcd_msa):
    
    valid_gpp_mask, nlcd_mask_dict = create_mask(gpp_msa_rr, ua_msa_rr, nlcd_msa)
    gpp_msa_rr_filled = np.copy(gpp_msa_rr)

    category_gpp_mean_list = []
    record = {'name': 'test_msa'}

    water_mask = (nlcd_msa > 10) & (nlcd_msa < 20)
    developed_mask = (nlcd_msa > 20) & (nlcd_msa < 40)

    gpp_msa_rr_filled[~valid_gpp_mask & water_mask] = 0
    gpp_msa_rr_filled[~valid_gpp_mask & developed_mask] = 0
                    
    for category, landcover_mask in nlcd_mask_dict.items():

        category_gpp_mean = np.nanmean(gpp_msa_rr[landcover_mask])
        record[category] = category_gpp_mean
        
        gpp_msa_rr_filled[~valid_gpp_mask & landcover_mask] = category_gpp_mean

    print(record)
    print("Gap filled GPP values for each landcover category have been calculated")
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
    nee_gpp_ratio_fine = np.empty_like(gpp_msa_rr_filled, dtype=np.float32)
    gpp_filled_coarse = get_gpp_coarse(gpp_msa_rr_filled, nee_msa, nlcd_clip_transform, nlcd_crs, nee_clip_transform, nee_crs)

    # reproject and resample NEE/GPP ratio to match gpp filled resolution
    reproject(
        source=nee_msa / gpp_filled_coarse,  
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


def main():
    import geopandas as gpd
    nee_memory = []
    nee_file = "../urban_greening/NEE.RS.FP-NONE.MLM-ALL.METEO-NONE.4320_2160.monthly.2015.nc"
    read_nee(nee_file, nee_transform, nee_memory)
    # print(nee_memory)

    msa_file = '../urban_greening/msa/michiganMSA_reprojected.shp' # crs: Albers Equal Area
    gpp_file = "../urban_greening/nov.15/michigan_test/modis-250-gpp-2015001.tif" # EPSG:4326
    nee_file = "../urban_greening/NEE.RS.FP-NONE.MLM-ALL.METEO-NONE.4320_2160.monthly.2015.nc" # EPSG:4326, resolution 1/12 degree
    nlcd_file = "../urban_greening/nov.15/nlcd_2016_land_cover_l48_20210604.img" # crs: Albers Equal Area, resolution 30m
    ua_file = "../urban_greening/gis_processed/ua/ua_30.tif" # crs: Albers Equal Area, resolution 30m

    msa=gpd.read_file(msa_file)
    msa=msa[msa['NAMELSAD']=='Grand Rapids-Wyoming-Kentwood, MI Metro Area'] # Subsetting to my AOI



    # gpp_msa_rr, ua_msa_rr, nlcd_msa, nlcd_clip_image, nlcd_crs, nlcd_clip_transform, nee_msa, nee_clip_transform, nee_crs = pipe(msa, gpp_file, nlcd_file, ua_file, nee_memory).values()
    pipe_output = pipe(msa, gpp_file, nlcd_file, ua_file, nee_memory)
    gpp_msa_rr = pipe_output['gpp_msa_rr']
    ua_msa_rr = pipe_output['ua_msa_rr']
    nlcd_msa = pipe_output['nlcd_msa']
    nlcd_clip_image = pipe_output['nlcd_clip_image']
    nlcd_crs = pipe_output['nlcd_crs']
    nlcd_clip_transform = pipe_output['nlcd_clip_transform']
    nee_msa = pipe_output['nee_msa']
    nee_clip_transform = pipe_output['nee_clip_transform']
    nee_crs = pipe_output['nee_crs']
    

    # from rasterio import plot
    import matplotlib.pyplot as plt
    
    # # Show plog overlays
    # fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    # rasterio.plot.show(gpp_msa_rr, ax=ax, alpha=0.8, title="Overlaying Rasters")
    # cax = ax.imshow(gpp_msa_rr, cmap='viridis')  # Use a colormap of your choice
    # cbar = fig.colorbar(cax, ax=ax, orientation='vertical')
    # cbar.set_label("GPP Values")
    # # rasterio.plot.show(nlcd_msa, ax=ax, alpha=1)  # Adjust alpha for transparency
    # # rasterio.plot.show(ua_msa_rr, ax=ax, alpha=0.5)  # Adjust alpha for transparency
    # plt.show()
    gpp_msa_rr_filled = gap_fill_gpp(gpp_msa_rr, ua_msa_rr, nlcd_msa)

    # fig, ax = plt.subplots(figsize=(10, 6)) 
    # cax = ax.imshow(gpp_msa_rr_filled, cmap='viridis')  # Use a colormap of your choice
    # cbar = fig.colorbar(cax, ax=ax, orientation='vertical')
    # cbar.set_label("GPP Values")
    # plt.show()
    nee_gpp_ratio_fine = get_nee_gpp_ratio_fine(gpp_msa_rr_filled, nee_msa, nlcd_clip_transform, nlcd_crs, nee_clip_transform, nee_crs)
    downscaled_nee = nee_gpp_ratio_fine * gpp_msa_rr_filled
    save_tiff(downscaled_nee, '../output/downscaled_nee_check.tif', nlcd_crs, nlcd_clip_transform)
    



if __name__ == "__main__":
    main()