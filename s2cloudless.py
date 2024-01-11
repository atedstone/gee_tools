"""
Cloud masking Sentinel-2 imagery on GEE

Most functionality copied from 
https://developers.google.com/earth-engine/tutorials/community/sentinel-2-s2cloudless
"""

import ee

CLOUD_FILTER = 60 # %
CLD_PRB_THRESH = 40 # %
NIR_DRK_THRESH = 0.15 # refl
CLD_PRJ_DIST = 1 # km
BUFFER = 100 # metres

def get_s2_sr_cld_col(
    bounds,
    start_date, 
    end_date
    ):
    
    def _get_collection(
        coll, 
        start_date=start_date, 
        end_date=end_date, 
        bounds=bounds
        ):
        return ee.ImageCollection(coll) \
            .filterDate(start_date, end_date) \
            .filterBounds(bounds)
    
    # Previously I was using S2_SR, but this is not widely available over Greenland earlier in the time series.
    # 2019 NDWI scenes onwards were generated with S2_SR data. 
    # 2017, 2018 generated with S2 collection (which contains S2 L1C)
    refl = _get_collection('COPERNICUS/S2').filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', CLOUD_FILTER))

    clouds = _get_collection('COPERNICUS/S2_CLOUD_PROBABILITY')
    
    return ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply(**{
        'primary': refl,
        'secondary': clouds,
        'condition': ee.Filter.equals(**{
            'leftField': 'system:index',
            'rightField': 'system:index'
        })
    }))


def add_cloud_bands(img):
    # Get s2cloudless image, subset the probability band.
    cld_prb = ee.Image(img.get('s2cloudless')).select('probability')

    # Condition s2cloudless by the probability threshold value.
    is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds')

    # Add the cloud probability layer and cloud mask as image bands.
    return img.addBands(ee.Image([cld_prb, is_cloud]))


def add_shadow_bands(img, use_SCL=True):
    SR_BAND_SCALE = 1e4
    if use_SCL:
        # For L2A collection
        # Identify not-water pixels from the SCL band.
        not_water = img.select('SCL').neq(6)

    else:
        # For L1C collection
        # Identify not-water pixels based on simple detection threshold ...
        not_water = img.normalizedDifference(['B4', 'B2']).lt(0.2)

    # Identify dark NIR pixels that are not water (potential cloud shadow pixels).
    dark_pixels = img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE).multiply(not_water).rename('dark_pixels')   

    # Determine the direction to project cloud shadow from clouds (assumes UTM projection).
    shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')));

    # Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
    cld_proj = (img.select('clouds').directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
        .reproject(**{'crs': img.select(0).projection(), 'scale': 100})
        .select('distance')
        .mask()
        .rename('cloud_transform'))

    # Identify the intersection of dark pixels with cloud shadow projection.
    shadows = cld_proj.multiply(dark_pixels).rename('shadows')

    # Add dark pixels, cloud projection, and identified shadows as image bands.
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))


def add_cld_shdw_mask(img, use_SCL):
    # Add cloud component bands.
    img_cloud = add_cloud_bands(img)

    # Add cloud shadow component bands.
    img_cloud_shadow = add_shadow_bands(img_cloud, use_SCL=use_SCL)

    # Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
    is_cld_shdw = img_cloud_shadow.select('clouds').add(img_cloud_shadow.select('shadows')).gt(0)

    # Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
    # 20 m scale is for speed, and assumes clouds don't require 10 m precision.
    is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(BUFFER*2/20)
        .reproject(**{'crs': img.select([0]).projection(), 'scale': 20})
        .rename('cloudmask'))

    # Add the final cloud-shadow mask to the image.
    return img_cloud_shadow.addBands(is_cld_shdw)


def apply_cld_shdw_mask(img):
    # Subset the cloudmask band and invert it so clouds/shadow are 0, else 1.
    not_cld_shdw = img.select('cloudmask').Not()

    # Subset reflectance bands and update their masks, return the result.
    return img.select('B.*').updateMask(not_cld_shdw)