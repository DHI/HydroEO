
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