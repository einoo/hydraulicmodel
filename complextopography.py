"""
Create the topography map for the 1D channel flow and 2D overland flow
"""

import logging
import numpy as np
import pandas as pd
import geopandas as gpd
from itertools import product


def complextopography(line_src, polygon_src, NX, NY,
                      cell_size, out_river, out_overland):
    """
    Mark the location of river in the DEM map

    The coordinate of the rivers will be marked by R0, R1, ..., RN as the input
    of the flood wave. The other part are the DEM for the flood propagation.

    Args:
        line_src (file):     Shape file of the river network
        polygon_src (file):  Discritized grid file of the study area
        NX (int):            Number of cells in X direction
        NY (int):            Number of cells in Y direction
        cell_size (float):   The size of cell during the simulation
        out_river (file):    The coordinates of rivers
        out_overland (file): The coordinates of the surface elevation

    Returns:
        overland_topo (array): 2D array for the inundation simulation
    """

    # polygon layer
    poly = gpd.read_file(polygon_src)
    # line layer
    line = gpd.read_file(line_src)

    # spatially join with intersection options
    inter = gpd.sjoin(line, poly, how='right', op='intersects')

    # delete the NaN value of the intersections
    # index value was kept the same as the grid cell in the geometry column
    inter = inter[~np.isnan(inter['index_left'])]

    # create the topography for the overland flow simulation
    topo = np.zeros((NY, NX))
    for i in range(NX):
        for j in range(NY):
            topo[j, i] = poly['ave_dem'][i * NY + j]
            if int(i * NY + j) in inter.index:
                topo[j, i] = 0             # arbitary value for the river

    inter.to_file(out_river)

    # save the complex topography to the csv data file
    topo = topo.tolist()
    df = pd.DataFrame(topo)
    df.to_csv('../Test/topo_CellSize%d.csv' % (cell_size))

    # revise the river bank elevation
    # based on the official website, the critical river water level is 7.66 m
    RIVER_BANK = 7.66
    topo_rev = pd.read_csv('../Test/topo_CellSize%d.csv' % (cell_size),
                           index_col=0)
    x, y = topo_rev.values.shape
    for i, j in product(range(1, x-1), range(y)):
        if topo_rev.iloc[i, j] == 0. and topo_rev.iloc[i-1, j] < RIVER_BANK:
            topo_rev.iloc[i-1, j] = RIVER_BANK
        if topo_rev.iloc[i, j] == 0. and topo_rev.iloc[i+1, j] < RIVER_BANK:
            topo_rev.iloc[i+1, j] = RIVER_BANK
    '''
    # boundary lines
    for i in range(x):
        if topo_rev.iloc[i, 0] == 0.:
            topo_rev.iloc[i, 0] = RIVER_BANK
        if topo_rev.iloc[i, -1] == 0.:
            topo_rev.iloc[i, -1] = RIVER_BANK
    for i in range(y):
        if topo_rev.iloc[0, i] == 0.:
            topo_rev.iloc[0, i] = RIVER_BANK
        if topo_rev.iloc[-1, i] == 0.:
            topo_rev.iloc[-1, i] = RIVER_BANK
    '''

    topo_rev.to_csv('../Test/topo_rev_CellSize%d.csv' % (cell_size),
                    float_format='%0.2f')

    logging.info('River length in shp: %0.2f.' % (line['geometry'].length))
    logging.info('Model river len:%0.2f.' % (inter.values.shape[0]*cell_size))

    return topo_rev
