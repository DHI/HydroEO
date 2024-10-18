
import typing
import numpy as np
import shapely
import geopandas as gpd

def grid_bounds(bounds: tuple, res: float) -> list :
    
    # unpack the bounds
    xmin, ymin, xmax, ymax = bounds

    # create the cells in a loop
    grid_cells = []
    for x0 in np.arange(xmin, xmax+res, res ):
        for y0 in np.arange(ymin, ymax+res, res):
            # bounds
            x1 = x0-res
            y1 = y0+res
            grid_cells.append( shapely.geometry.box(x0, y0, x1, y1) )

    return grid_cells

def find_closest_geom(geom: shapely.Geometry , geoms: list) -> tuple :
    """finds the index of the geometry within the provided list that is closest to the main geometry provided

    Parameters
    ----------
    point : shapely.Geometry 
        point in which to compare all other points to.
    points : list
        list of points to compare to the main point

    Returns
    -------
    tuple
        tuple containg the index of the point within the list that is closest to the main point, the point object, and the distance
    """
    # get distances betewen the point and the individual lines
    distance_list = [geom.distance(sub_geom) for sub_geom in geoms]

    # find the index of the shortest distance
    index = np.argmin(distance_list)

    return (index, geoms[index], distance_list[index])