"""
Discretize the land use and cover data
"""

import geopandas as gpd
import logging

def landuse(grid, lulc, bldg, outlulc, outbldg):
    """
    Read the land use land cover, and building distributions

    Discretize the land use land cover data based on the grid cell.
    Add the building information to the geodataframe, for the risk assessment. 

    Args:
        grid (file): The grid cell shape file for this simulation (polygon)
        lulc (file): The land use land cover shape file (polygon)
        bldg (file): The building distribution shape file (polygon)
        outlulc (file): Discretized land use land cover shape file
        outbldg (file): Discretized building distribution shape file

    Returns:
        GeoDataFrame: With the land use land cover information as CODE3
        GeoDataFrame: With the building distributions as ORGGILVL
    """

    logging.info('The grid cell of this simulations is %s.' % (grid))
    logging.info('The land use land cover original data is %s.' % (lulc))
    logging.info('The building distribution original data is %s.' % (bldg))
    # grid cell polygon layer
    cell = gpd.read_file(grid)
    # land use land cover polygon layer
    landuse = gpd.read_file(lulc)
    # building distribution polygon layer
    building = gpd.read_file(bldg)

    luce = gpd.sjoin(landuse, cell, how='right', op='intersects')
    blce = gpd.sjoin(building, cell, how='right', op='intersects')

    luce.to_file(outlulc)
    blce.to_file(outbldg)

    return luce, blce
