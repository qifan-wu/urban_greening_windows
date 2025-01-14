from utils_250_pipe import *

# MSA_FILE = '../gis/msa/msaUS_mland_aea1_M1.shp' # for US without MSA from Hawaii, Puerto Rico and Alaska (no NLCD or no carbon data), and without some msa in midwest (see removed_msa.txt); crs: Albers Equal Area
MSA_FILE = '../gis/msa/msaUS_mland_aea1_M1_all.shp' # for US without MSA from Hawaii, Puerto Rico and Alaska (no NLCD or no carbon data), crs: Albers Equal Area
UA_FILE = "../gis/ua/ua_us_30_clip1.tif" # crs: Albers Equal Area, resolution 30m


# def pipe_downscaled_nee_msa(msa_ds, msa, gpp_file, nlcd_file, ua_file, memfile_nee):
#      # Subsetting to my AOI
#     msa_name = msa['NAMELSAD'].values[0]

#     pipe_output = pipe_read_gen_params(msa, gpp_file, nlcd_file, ua_file, memfile_nee)
#     gpp_msa_rr = pipe_output['gpp_msa_rr']
#     ua_msa_rr = pipe_output['ua_msa_rr']
#     nlcd_msa = pipe_output['nlcd_msa']
#     nlcd_crs = pipe_output['nlcd_crs']
#     nlcd_clip_transform = pipe_output['nlcd_clip_transform']
#     nee_msa = pipe_output['nee_msa']
#     nee_clip_transform = pipe_output['nee_clip_transform']
#     nee_crs = pipe_output['nee_crs']


#     from rasterio.transform import Affine
#     target_transform = nlcd_clip_transform * Affine.scale(250 / 30)
#     gpp_msa_rr_filled_30m = gap_fill_gpp(gpp_msa_rr, ua_msa_rr, nlcd_msa, msa_name, save_mean_csv=True) # change to false if don't want mean gpp as csv
#     gpp_msa_rr_filled_250m = reproject_gpp_filled(gpp_msa_rr_filled_30m, nlcd_clip_transform, nlcd_crs, target_resolution=250, target_transform=target_transform)

#     nee_gpp_ratio_fine = get_nee_gpp_ratio_fine(gpp_msa_rr_filled_250m, nee_msa, target_transform, nlcd_crs, nee_clip_transform, nee_crs)
    
#     # testmem = create_in_memory_ds(nee_gpp_ratio_fine, nlcd_crs, target_transform, return_file=True) # test only, delete later
#     # test_ratio_list.append(testmem) # test only, delete later

#     downscaled_nee = nee_gpp_ratio_fine * gpp_msa_rr_filled_250m

#     downscaled_nee_info = {
#         'data': downscaled_nee,
#         'crs': nlcd_crs,
#         'transform': target_transform,  # Update with 250m transform
#     }

#     return downscaled_nee_info

    
def downscale_pipe(year_month, nee_memory):
    '''
    year_month: str, e.g. '201501'
    '''
    import geopandas as gpd

    
    year = int(year_month[:4])
    month = int(year_month[4:])
    

    # Prepare files
    gpp_file = f"../gis/GPP_monthly_mean/gpp_{year_month}.tif" # EPSG:4326
    
    if year < 2003:
        nlcd_year = '2001'
    elif year < 2006:
        nlcd_year = '2004'
    elif year < 2008:
        nlcd_year = '2006'
    elif year < 2010:
        nlcd_year = '2008'
    elif year < 2013:
        nlcd_year = '2011'
    elif year < 2015:
        nlcd_year = '2013'
    elif year >= 2015:
        nlcd_year = '2016'

    nlcd_file = f"../gis/NLCD/nlcd_{nlcd_year}_land_cover_l48_20210604.img" # crs: Albers Equal Area, resolution 30m

    msa_ds=gpd.read_file(MSA_FILE)

    

    memfile_nee = nee_memory[month - 1]

    # ======== test with all msas ========
    mem_downscaled_nee_list = []
    for index, record in msa_ds.iterrows(): #debug: change to msa_ds[:3}.iterrows()
        msa_name = record['NAMELSAD']
        print(f'Generating downsclaed data for {msa_name}...')
        msa = msa_ds.loc[[index]]
        
        downscaled_nee_msa = pipe_downscaled_nee_msa(msa_ds, msa, gpp_file, nlcd_file, UA_FILE, memfile_nee)
        mem = create_in_memory_ds(downscaled_nee_msa['data'], downscaled_nee_msa['crs'], downscaled_nee_msa['transform'], return_file=True)
        mem_downscaled_nee_list.append(mem)
        
        


    datasets_ratio = [mem.open() for mem in test_ratio_list]
    merged_data_ratio, merged_transform_ratio = rasterio.merge.merge(datasets_ratio, nodata=np.nan)
    save_tiff(merged_data_ratio[0], f'../gis/output/downscaleRatio/ratio_us_{year_month}.tif', datasets_ratio[0].crs, merged_transform_ratio)

    # Merge datasets
    datasets = [mem.open() for mem in mem_downscaled_nee_list]
    print("Merging datasets. This might take a while...")
    merged_data, merged_transform = rasterio.merge.merge(datasets, nodata=np.nan)
    
    print("Merging completed")

    merged_raster = merged_data[0] # get the first band
    output_file = f'../gis/output/downscaledNEE/downscaledNEE_US_{year_month}.tif'
    save_tiff(merged_raster, output_file, datasets[0].crs, merged_transform)

    # Close datasets
    for ds in datasets:
        ds.close()
    # ======== test with all msas ========
    


def main():
    print(f"==== Downscaling Start ====")
    for year in range(2003, 2016): #TODO change to (2001, 2016)
        print(f"Prepare Raw NEE for {year}...")
        nee_memory = []
        nee_file = f"../gis/NEE/NEE.RS.FP-NONE.MLM-ALL.METEO-NONE.4320_2160.monthly.{year}.nc" # EPSG:4326, resolution 1/12 degree
        read_nee(nee_file, nee_transform, nee_memory)

       
        for month in range(1,13): 
            print("Downscaling NEE for", f"{year:04d}{month:02d}...")

            # clear the global variable every month
            gpp_mean_cat_data.clear()
            test_ratio_list.clear # test only, delete later

            downscale_pipe(f"{year:04d}{month:02d}", nee_memory)

            # Save gpp_mean_values to csv
            gpp_mean_data_df = pd.DataFrame(gpp_mean_cat_data)
            gpp_mean_data_df.to_csv(f'../gis/output/statistics/gpp_mean_data_250_{year:04d}{month:02d}.csv', index=False)

        print(f"--- Finish downscaling NEE for {year} ---")
    print(f"==== Downscaling Finished ====")

if __name__ == "__main__":
    main()