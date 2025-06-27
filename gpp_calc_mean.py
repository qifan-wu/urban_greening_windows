import os
import re
import numpy as np
import rasterio
from rasterio.plot import show
from collections import defaultdict

# Directory containing .tif files
input_directory = "../gis/GPP/" 
output_directory = "../gis/GPP_monthly_mean/" 

def get_month_from_day(year, day):
    # Check if the year is a leap year
    is_leap_year = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    
    # Cumulative days for months in a regular year
    days_in_month = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
    
    # Adjust for leap year
    if is_leap_year:
        days_in_month = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
    
    # Determine the month
    for month in range(1, 13):
        if day <= days_in_month[month]:
            return month


def main():
    # Organize files by year and month
    files_by_month = defaultdict(list)

    for file_name in os.listdir(input_directory):
        if file_name.endswith(".tif"):
            year = int(file_name[-11:-7])
            day = int(file_name[-7:-4])
            month = get_month_from_day(year, day)
            month = f"{month:02d}"
            files_by_month[f"{year}{month}"].append(os.path.join(input_directory, file_name))

    # for k,v in files_by_month.items():
    #     print(k, len(v))
    print("Calculating monthly mean rasters begins...")
    # Process each year-month group
    for year_month, file_list in files_by_month.items():
        data_arrays = []
        for year_month, file_list in files_by_month.items():
            # Read data for the month
            for file_path in file_list:
                with rasterio.open(file_path) as src:
                    data = src.read(1)
                    data_arrays.append(data)
                    # Capture metadata from the first file
                    if len(data_arrays) == 1:
                        meta = src.meta.copy()

                # Calculate the mean for the month
                mean_array = np.mean(data_arrays, axis=0)

                # Cast to integer type (e.g., uint16)
                mean_array = mean_array.astype(np.uint16)

                # Update metadata for output
                meta.update(dtype=rasterio.uint16, count=1)

            # Output file name
            output_file = os.path.join(output_directory, f"gpp_{year_month}.tif")

            # Save the mean array as a new TIFF file
            with rasterio.open(output_file, 'w', **meta) as dst:
                dst.write(mean_array.astype(rasterio.float32), 1)

            print(f"Monthly mean raster saved as {output_file}")
    print("===Calculating monthly mean rasters finished===")

if __name__ == "__main__":
    main()