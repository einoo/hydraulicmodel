"""
Discretization of the study area based on Cell Size
"""

import logging
import math
from itertools import product
import pandas as pd
import numpy as np
import geopandas as gpd
import shapely.geometry
from shapely.geometry import Point


def discretization_CellSize(src, n, outshp, outcsv):
    """
    Discretization of the study area based on the cell size

    Save the DEM of the study area as csv file including the x and y
    coordinates and the elevation in each point.  The spatial operation of
    the discretization is operated by geopandas.  The elevation in each grid
    is the mean value of the elevation of the points within the grid.

    Args:
        src (file): The source file of the elevation data including the (x, y)
                    coordinates and elevations
        n (int): Cell size in both x and y directions
        outshp (file): The output shape file of the grids with elevations
        outcsv (file): The output csv file of the grids with elevations

    Returns:
        NX (int):          The number of cells in the x-axis
        NY (int):          The number of cells in the y-axis
    """

    df = pd.read_csv(src, index_col=None)
    geometry = [Point(xy) for xy in zip(df.iloc[:, 0], df.iloc[:, 1])]

    # specify the WGS84 / UTM N54 coordinate system
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs='EPSG:32654')
    gdf = gdf.drop(columns=['POINT_X', 'POINT_Y'])
    logging.info('The coordinate system of the source file is %s.' % (gdf.crs))

    # total area for the grid
    xmin, ymin, xmax, ymax = gdf.total_bounds
    logging.info('The boundary of the study area is %s, %s, %s, %s.'
                 % (xmin, xmax, ymin, ymax))

    # calculate the number of cells in x and y direction
    cell_nox = math.ceil((xmax-xmin)/n)
    cell_noy = math.ceil((ymax-ymin)/n)
    logging.info('Numbers of cells in x and y are: %d, %d' %
                 (cell_nox, cell_noy))

    # create the cells in a loop
    grid_cells = []
    for x0, y0 in product(np.arange(xmin, xmin + cell_nox * n, n),
                          np.arange(ymin, ymin + cell_noy * n, n)):
        # polygon boundary
        x1 = x0 + n
        y1 = y0 + n
        # amend the data to the polygon
        grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))

    # convert the shape format to the geopandas format
    cell = gpd.GeoDataFrame(grid_cells, columns=['geometry'],
                            crs='EPSG:32654')
    logging.info('The cells of the study area are created.')

    # join the cell with DEM
    merged = gpd.sjoin(gdf, cell, how='right', op='within')
    merged = merged.drop(columns=['index_left'])

    # compute the average/mean DEM in each cell
    dissolve = merged.groupby(level=0).mean()

    # define the boundary NaN value to the maximum elevation
    MAX_ELE = np.max(dissolve)
    # along x-axias bottom
    for i in range(cell_noy):
        if np.isnan(dissolve['GRID_CODE'][i*cell_nox]):
            dissolve['GRID_CODE'][i*cell_nox] = MAX_ELE
    # along x-axis top
    for i in range(1, cell_noy + 1):
        if np.isnan(dissolve['GRID_CODE'][i*cell_nox-1]):
            dissolve['GRID_CODE'][i*cell_nox-1] = MAX_ELE
    # along y-axis left
    for i in range(cell_nox + 1):
        if np.isnan(dissolve['GRID_CODE'][i]):
            dissolve['GRID_CODE'][i] = MAX_ELE
    # along y-axis right
    for i in range(cell_nox * cell_noy - cell_noy, cell_nox * cell_noy):
        if np.isnan(dissolve['GRID_CODE'][i]):
            dissolve['GRID_CODE'][i] = MAX_ELE
    # the other values in the body
    for i in range(cell_nox * cell_noy):
        if np.isnan(dissolve['GRID_CODE'][i]):
            dissolve['GRID_CODE'][i] = dissolve['GRID_CODE'][i-1]

    # ensemble the average DEM to the grid cells
    cell.loc[dissolve.index, 'ave_dem'] = dissolve.GRID_CODE.values

    # save the discretized file to the shape file and csv file
    cell.to_file(outshp)
    cell.to_csv(outcsv)

    return cell_nox, cell_noy
