from utils_250_pipe import *

# MSA_FILE = '../gis/msa/msaUS_mland_aea1_M1.shp' # for US without MSA from Hawaii, Puerto Rico and Alaska (no NLCD or no carbon data), and without some msa in midwest (see removed_msa.txt); crs: Albers Equal Area
MSA_FILE = '../gis/msa/msaUS_mland_aea1_M1_all.shp' # for US without MSA from Hawaii, Puerto Rico and Alaska (no NLCD or no carbon data), crs: Albers Equal Area
# MSA_FILE = '../gis/msa/msaUS_mland_aea1_M1_michigan.shp' # Michigan MSA for test
UA_FILE = "../gis/ua/ua_us_30_clip1.tif" # crs: Albers Equal Area, resolution 30m

    
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
        # print(f'Generating downscaled data for {msa_name}...')
        msa = msa_ds.loc[[index]]
        
        downscaled_nee_msa = pipe_downscaled_nee_msa(msa_ds, msa, gpp_file, nlcd_file, UA_FILE, memfile_nee)
        mem = create_in_memory_ds(downscaled_nee_msa['data'], downscaled_nee_msa['crs'], downscaled_nee_msa['transform'], return_file=True)
        mem_downscaled_nee_list.append(mem)
        

    # datasets_ratio = [mem.open() for mem in test_ratio_list]
    # merged_data_ratio, merged_transform_ratio = rasterio.merge.merge(datasets_ratio, nodata=np.nan)
    # save_tiff(merged_data_ratio[0], f'../gis/output/downscaleRatio/ratio_us_{year_month}.tif', datasets_ratio[0].crs, merged_transform_ratio)

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
    for year in range(2001, 2016): #TODO change to (2001, 2016)
        print(f"Prepare Raw NEE for {year}...")
        nee_memory = []
        nee_file = f"../gis/NEE/NEE.RS.FP-NONE.MLM-ALL.METEO-NONE.4320_2160.monthly.{year}.nc" # EPSG:4326, resolution 1/12 degree
        read_nee(nee_file, nee_transform, nee_memory)


        for month in range(1,13):
            print("Downscaling NEE for", f"{year:04d}{month:02d}...")

            # clear the global variable every month
            gpp_mean_cat_data.clear()

            downscale_pipe(f"{year:04d}{month:02d}", nee_memory)

            # Save gpp_mean_values to csv
            gpp_mean_data_df = pd.DataFrame(gpp_mean_cat_data)
            # gpp_mean_data_df.to_csv(f'../gis/output/statistics/gpp_mean_data_250_{year:04d}{month:02d}.csv', index=False)

            # Save ratio to tif
            datasets_ratio = [mem.open() for mem in test_ratio_list]
            merged_data_ratio, merged_transform_ratio = rasterio.merge.merge(datasets_ratio, nodata=np.nan)
            save_tiff(merged_data_ratio[0], f'../gis/output/downscaleRatio/ratio_us_{year:04d}_{month:02d}.tif', datasets_ratio[0].crs, merged_transform_ratio)

            # Close memory for ratio
            for memfile in test_ratio_list:
                memfile.close()
            test_ratio_list.clear() 

        print(f"--- Finish downscaling NEE for {year} ---")
        for memfile in nee_memory:
            memfile.close() 
        nee_memory.clear()

        
    print(f"==== Downscaling Finished ====")

if __name__ == "__main__":
    main()