
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


def line_to_points(line: shapely.LineString, delta: float) -> tuple :
    """returns a list of points and a list of distances at a given interval along the provided linestring.
    Note that the last point will not be returned unless it lies along a multiple of the interval. This option coudl eb included if needed.

    Parameters
    ----------
    line : shapely.LineString
        the line to get points along
    delta : float
        the distance that should be between each point

    Returns
    -------
    tuple
        - list of points along linstring
        - list of distances associated with the points
    """
    points = list()
    distances = list()

    # calculate the length of the linestring
    length = line.length

    # loop through possible distances at that interval and record the point and distance
    for dist in np.arange(0, length+delta, delta):
        points.append(line.interpolate(dist))
        distances.append(dist)

    return points, distances


def format_coord_list(coords : list):

    if not hasattr(coords[0], '__iter__'):

        # make sure we have an even number of elements in the list
        if len(coords)%2 > 0:
            raise ValueError("The inputed list of elements is not even. Please retry the query with an even list of coordinates or a list of coordinate pairs")

        # unpack the 1d list into coordinate pairs
        new_list = list()
        while (len(coords) > 0):
            new_list.append((coords.pop[0], coords.pop(0)))
        coords = new_list

    return coords