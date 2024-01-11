"""
Helper functions for GEE

AT, March 2022
"""
import datetime as dt
import ee
import pandas as pd

CRS='epsg:3413'

def get_date(im, attr='system:time_start'):
    """ Return Datetime of image timestamp 
    
    im : ee.Image
    attr : attribute in which date is stored

    returns : dt.datetime

    """
    def cb(im):
        return 0
    # Get timestamp in milliseconds POSIX format
    time_ms = im.get('system:time_start').getInfo(cb)
    # Convert to seconds POSIX time and convert to datetime
    return dt.datetime.fromtimestamp(time_ms * 1e-3)


def shapely_to_ee_poly(geom, crs):
    """ Given a shapely Polygon, return a ee Polygon. """
    xx, yy = geom.exterior.xy
    eep = ee.Geometry.Polygon([
        [[x, y] for x, y in zip(xx, yy)]
        ],
        geodesic=False, 
        proj=crs
        )
    return eep


def create_reduce_regions(collection, reducer, add_date=True, add_props=None, crs=CRS, scale=20):
    """
    Creates a function which can be mapped to images to extract region statistics from them.

    :param collection: FeatureCollection of regions over which to reducer
    :param reducer: an Earth Engine reducer
    :param add_date: if True (default), attach date to each returned region
    :param add_props: List of image properties to retrieve (don't list bands to be reduced here)
    :param crs: Coordinate reference system in which to apply the reducer
    :param scale: metres
    """
    def reduce_regions(img):
        """
        collection : FeatureCollection of image regions
        Call this function on an ImageCollection using its map() call.
        """
        d = img.reduceRegions(
            collection,
            reducer=reducer,
            crs=crs,
            scale=scale) #meters I believe

        def set_props(f):
            for kw in add_props:
                f = f.set({kw: img.get(kw)})
            return f
        if add_props is not None:
            d = d.map(set_props)
            
        # All the patches from this one image are returned within a FeatureCollection...
        # ...so add the (same) image date to each Feature in the FeatureCollection.
        if add_date:
            def set_date(f):
                return f.set({'millis': img.date().millis()}) # 
            dd = d.map(set_date)
            return dd
        
        return d
        
    return reduce_regions


def make_roi(xmin, ymin, xmax, ymax, crs=CRS):
    roi = ee.Geometry.Rectangle([ xmin, ymin, xmax, ymax ], 
            crs, 
            geodesic=False, evenOdd=True)
    return roi


def fc_to_dict(fc):
    """ Define a function to transfer FeatureCollection properties to an ee.Dict. 

    As of 09/2022, not in use for the Sentinel-1 retention project.
    """
    raise NotImplementedError
    prop_names = fc.first().propertyNames()
    prop_lists = fc.reduceColumns(
      reducer=ee.Reducer.toList().repeat(prop_names.size()),
      selectors=prop_names).get('list')
    return ee.Dictionary.fromLists(prop_names, prop_lists)


def fc_to_df(fc, selectors=None):
    """ Convert a FeatureCollection to a pd.DataFrame. """
    if selectors is None:
        selectors = fc.first().propertyNames().getInfo() # convert to local list, needed for Pandas later
        selectors_n = len(selectors)
    else:
        assert isinstance(selectors, list)
        selectors_n = len(selectors)

    as_list = fc.reduceColumns(reducer=ee.Reducer.toList(tupleSize=selectors_n), selectors=selectors).get('list').getInfo()
    df = pd.DataFrame(as_list, columns=selectors)
    if 'millis' in df.columns:
        df['timestamp'] = pd.to_datetime(df['millis'], unit='ms')
        df.index = df.timestamp
    return df


def create_df(dictobj):
    """ Unclear purpose, so deprecated as of 09/2022. """
    raise NotImplementedError
    tmp = fc_to_dict(dictobj).getInfo()
    df = pd.DataFrame(tmp)
    df['timestamp'] = pd.to_datetime(df['millis'], unit='ms')
    df.index = df.timestamp
    return df


def get_timeseries(collection, reduce_function, selectors=None):
    """ Convenience function to get a time series of region(s) in pandas format.

    :param collection: ImageCollection over which to map.
    :param reduce_function: reducer such as that produced by create_reduce_regions.
    :returns: pd.DataFrame
    """

    ts = collection.map(reduce_function)
    if selectors is None:
        selectors = list(ts.first().getInfo()['columns'].keys())
    ts_f = ts.flatten()
    ts_pd = fc_to_df(ts_f, selectors=selectors)
    return ts_pd