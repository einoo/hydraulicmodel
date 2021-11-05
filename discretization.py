"""
Discretization of the study area
"""

import pandas as pd
import numpy as np
import geopandas
import logging
import shapely.geometry
from shapely.geometry import Point


def discretization(src, n, outshp, outcsv):
    """
    Discretization of the study area based on the number of cells in x-axis

    Save the DEM of the study area as csv file including the x and y
    coordinates and the elevation in each point.  The spatial operation of
    the discretization is operated by geopandas.  The elevation in each grid
    is the mean value of the elevation of the points within the grid.

    Args:
        src (file): The source file of the elevation data including the (x, y)
                    coordinates and elevations
        n (int): Discreted number of cells in the x-axis
        outshp (file): The output shape file of the grids with elevations
        outcsv (file): The output csv file of the grids with elevations

    Returns:
        cell_size (float): The cell size of the discretization
        NX (int):          The number of cells in the x-axis
        NY (int):          The number of cells in the y-axis
    """

    df = pd.read_csv(src, index_col=None)
    geometry = [Point(xy) for xy in zip(df.iloc[:, 0], df.iloc[:, 1])]

    # specify the WGS84 / UTM N54 coordinate system
    gdf = geopandas.GeoDataFrame(df, geometry=geometry, crs='EPSG:32654')
    gdf = gdf.drop(columns=['POINT_X', 'POINT_Y'])
    logging.info('The coordinate system of the source file is %s.' % (gdf.crs))

    # total area for the grid
    xmin, ymin, xmax, ymax = gdf.total_bounds
    logging.info('The boundary of the study area is %s, %s, %s, %s.'
                 % (xmin, xmax, ymin, ymax))

    # calculate the cell size
    cell_size = (xmax-xmin)/n
    logging.info('The cell size is %0.1f.' % (cell_size))

    # create the cells in a loop
    grid_cells = []
    cell_no = 0
    cell_nox = 0
    for x0 in np.arange(xmin, xmax+cell_size, cell_size):
        cell_nox += 1
        for y0 in np.arange(ymin, ymax+cell_size, cell_size):
            # polygon boundary
            x1 = x0-cell_size
            y1 = y0+cell_size
            # amend the data to the point and polygon
            grid_cells.append(shapely.geometry.box(x0, y0, x1, y1))
            cell_no += 1
    logging.info('The total number of cells are %d.' % (cell_no))
    logging.info('The number of cells in horizontal axis are %d.' % (cell_nox))

    # convert the shape format to the geopandas format
    cell = geopandas.GeoDataFrame(grid_cells, columns=['geometry'],
                                  crs='EPSG:32654')
    logging.info('The cells of the study area are created.')

    # join the cell with DEM
    merged = geopandas.sjoin(gdf, cell, how='right', op='within')
    merged = merged.drop(columns=['index_left'])

    # compute the average/mean DEM in each cell
    dissolve = merged.groupby(level=0).mean()

    # define the boundary NaN value to the maximum elevation
    MAX_ELE = np.max(dissolve)
    cell_noy = int(cell_no/cell_nox)
    # along x-axias bottom
    for i in range(cell_noy):
        if np.isnan(dissolve['GRID_CODE'][i*89]):
            dissolve['GRID_CODE'][i*89] = MAX_ELE
    # along x-axis top
    for i in range(1, cell_noy + 1):
        if np.isnan(dissolve['GRID_CODE'][i*89-1]):
            dissolve['GRID_CODE'][i*89-1] = MAX_ELE
    # along y-axis left
    for i in range(cell_nox + 1):
        if np.isnan(dissolve['GRID_CODE'][i]):
            dissolve['GRID_CODE'][i] = MAX_ELE
    # along y-axis right
    for i in range(cell_no - cell_noy, cell_no):
        if np.isnan(dissolve['GRID_CODE'][i]):
            dissolve['GRID_CODE'][i] = MAX_ELE
    # the other values in the body
    for i in range(cell_no):
        if np.isnan(dissolve['GRID_CODE'][i]):
            dissolve['GRID_CODE'][i] = dissolve['GRID_CODE'][i-1]

    # ensemble the average DEM to the grid cells
    cell.loc[dissolve.index, 'ave_dem'] = dissolve.GRID_CODE.values

    # save the discretized file to the shape file and csv file
    cell.to_file(outshp)
    cell.to_csv(outcsv)

    return cell_size, cell_nox, cell_noy
