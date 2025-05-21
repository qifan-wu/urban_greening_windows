
MSA_FILE = r'c:\Users\qifanw\Documents\gis\msa\msaUS_mland_aea1_M1_all.shp'
MSA_REGION_FILE = r'c:\Users\qifanw\Documents\gis\msa\region\msa_all_region.shp'
MSA_NAME_REGION_TABLE = r'c:\Users\qifanw\Documents\data\intermedia\msa_name_region.csv'

# get NEE crs which is equal to GPP crs (EPSG:4326)
# gpp_file = f'../gis/GPP_monthly_mean/gpp_200101.tif'
# with rasterio.open(gpp_file) as gpp_dstrd:
#     gpp_crs = gpp_dstrd.crs
# nee_crs = gpp_crs
NEE_CRS = "EPSG:4326"
GPP_CRS = "EPSG:4326"

from rasterio.transform import Affine
# # get original NEE transform
# minx, miny, maxx, maxy = -180.0, -90.0, 180.0, 90.0
# nee_orig_resolution_x = 1/12
# nee_orig_resolution_y = 1/12
# nee_orig_transform = rasterio.transform.from_origin(minx, maxy, nee_orig_resolution_x, nee_orig_resolution_y)
# nee_orig_transform
# NEE_ORIG_TRANSFORM = Affine(*nee_orig_transform)
NEE_ORIG_TRANSFORM = Affine(*(0.08333333333333333, 0.0, -180.0,
       0.0, -0.08333333333333333, 90.0))

NEE_MSA_MONTHLY_MEAN = r'c:\Users\qifanw\Documents\data\intermedia\msa_nee_mean.csv'
NEE_MSA_MONTHLY_MEAN_XBASE = r'c:\Users\qifanw\Documents\data\intermedia\msa_nee_mean_fluxx.csv'

CLIMATE_ZONE = {
    'northeast': [
        'CT', 'DE', 'ME', 'MD', 'MA', 'NH', 'NJ', 'NY', 'PA', 'VT',
        'MA-NH', 'MD-WV', 'NY-NJ', 'PA-NJ', 'PA-NJ-DE-MD', 'RI-MA',
        'DC-VA-MD-WV', 'VA-WV'
    ],
    'upper_midwest': [
        'IA', 'MI', 'MN', 'WI',
        'MN-WI', 'WI-MN', 'SD-MN'
    ],
    'ohio_valley': [
        'IL', 'IN', 'KY', 'MO', 'OH', 'TN', 'WV',
        'KY-IL', 'KY-IN', 'MO-IL', 'OH-KY-IN', 'WV-OH', 'WV-KY-OH',
        'IA-IL', 'IL-IN', 'IN-MI', 'MO-KS', 'TN-KY'
    ],
    'southeast': [
        'AL', 'FL', 'GA', 'NC', 'SC', 'VA',
        'GA-AL', 'GA-SC', 'NC-SC', 'TN-GA', 'TN-VA', 'VA-NC'
    ],
    'northern_rockies_plains': [
        'MT', 'NE', 'ND', 'SD', 'WY',
        'IA-NE-SD', 'ND-MN', 'NE-IA'
    ],
    'south': [
        'AR', 'KS', 'LA', 'MS', 'OK', 'TX',
        'TX-AR', 'AR-OK', 'TN-MS-AR'
    ],
    'southwest': [
        'AZ', 'CO', 'NM', 'UT',
        'UT-ID'
    ],
    'northwest': [
        'ID', 'OR', 'WA',
        'ID-WA', 'OR-WA'
    ],
    'west': [
        'CA', 'NV'
    ]
}

# Meterological data MSA mean
PPT_MEAN_FILE = r'c:\Users\qifanw\Documents\data\intermedia\prism_msa_mean\ppt_msa_mean_2001-2015.csv'
TMEAN_MEAN_FILE = r'c:\Users\qifanw\Documents\data\intermedia\prism_msa_mean\tmean_msa_mean_2001-2015.csv'
TDMEAN_MEAN_FILE = r'c:\Users\qifanw\Documents\data\intermedia\prism_msa_mean\tdmean_msa_mean_2001-2015.csv'
SIF_MEAN_FILE = r'c:\Users\qifanw\Documents\data\intermedia\green_indices_msa_mean\sif_msa_mean_2001-2021.csv'

FF_MONTHLY_MEAN = r'c:\Users\qifanw\Documents\data\intermedia\ODIAC2000To2021mean.csv'

CS_TREND = r'c:\Users\qifanw\Documents\data\intermedia\trend\cs_trend.csv'