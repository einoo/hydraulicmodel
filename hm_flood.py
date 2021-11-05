'''
The work flow of the Hydraulic model.
Including the
(1) discretization on DEM
(2) precipitation
(3) land use and cover map
(4) evapotranspiration
(5) river and overland flow distinguish
(6) ...

Author: Mao Ouyang
Date:   2021-08-02
Email:  einooumo@hotmail.com
'''

import logging  # Log File


def main():

    # Log information
    logging.basicConfig(
            filename='HydraulicModel_Flood.log',
            format='[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%m-%d %H:%M',
            level=logging.INFO)

    # (1) discretize the DME of study area
    logging.info('Start the study area discretization.')
    import discretization_CellSize as dc
    domain = '../Test/Yachiyo_DEM.csv'
    cell_size = 10
    out_shp = '../Test/Yachiyo_grid_CellSize%d.shp' % (cell_size)
    out_csv = '../Test/Yachiyo_grid_CellSize%d.csv' % (cell_size)
    NX, NY = dc.discretization_CellSize(domain, cell_size, out_shp, out_csv)
    logging.info('End the study area discretization.')

    # (5) river and overland flow distinguish
    logging.info('2D topography; river/floodplain separate.')
    import complextopography as ct
    river_network = '../Test/Yachiyo_river_epsg32654.shp'
    grid_file = '../Test/Yachiyo_grid_CellSize%d.shp' % (cell_size)
    out_river = '../Test/Yachiyo_river_network.shp'
    out_overland = '../Test/Yachiyo_overland_CellSize%d.csv' % (cell_size)
    topo = ct.complextopography(river_network, grid_file, NX, NY, cell_size,
                                out_river, out_overland)
    logging.info('End topography creation and river/floodplain separation.')


if __name__ == '__main__':
    main()
