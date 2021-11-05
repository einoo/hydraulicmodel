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
Date:   2021-01-31
Email:  einooumo@hotmail.com
'''

import logging                          # Log File
from datetime import datetime, timedelta    # Date Time package

def main():

    # Log information
    logging.basicConfig(
            filename='HydraulicModel.log',
            format='[%(asctime)s] %(levelname)s: %(message)s',
            datefmt='%m-%d %H:%M',
            level=logging.INFO)

    # (1) discretize the DME of study area
    logging.info('Start the study area discretization.')
    import discretization as dc
    domain = '../Test/Yachiyo_DEM.csv'
    n_cells = 100
    out_shp = '../Test/Yachiyo_grid.shp'
    out_csv = '../Test/Yachiyo_grid.csv'
    cell_size = dc.discretization(domain, n_cells, out_shp, out_csv)
    logging.info('End the study area discretization.')

    # (2) precipitation data
    logging.info('Read the precipitation data.')
    start_date = datetime(2019, 10, 25, 0, 0)
    end_date = datetime(2019, 10, 28, 0, 0)
    delta = timedelta(minutes=10)
    time_step = 0
    while start_date <= end_date:
        time_step += 1
        start_date += delta
    logging.info('The total time steps in the simulation is %s.' % (time_step))
    logging.info('Assume precp is uniformly distributed in study area.')
    import precipitation as prec
    prec_src = '../Test/Precipitation_Yachiyo.csv'
    pr = prec.precipitation(prec_src, time_step)
    logging.info('End the precipitation data input.')

    # (3) land cover and use map
    logging.info('Read the land use land cover, building information.')
    lulc_data = '../Test/Yachiyo_landform.shp'
    bldg_data = '../Test/Yachiyo_buildings.shp'
    out_lulc = '../Test/Yachiyo_landform_discrete.shp'
    out_bldg = '../Test/Yachiyo_buildings_discrete.shp'
    import landusecover as luc
    land = luc.landuse(out_shp, lulc_data, bldg_data, out_lulc, out_bldg)
    logging.info('End land use, building information discretization.')

    # (4) evapotranspiration
    logging.info('Calculate crop evaportranspiration (ETC) based on ET0')
    day_of_year = [298, 299, 300]  # corresponded to Oct 25, 26, 27, 2019
    latitude = 35.423348
    sunshine_hours = 11.0
    altitude = 8.80
    tmin = 14.0
    tmax = 22.0
    ws = 2.0
    z = 10.0
    import evapotranspiration as et
    eto = []
    for day in day_of_year:
        eto.append(et.daily_eto(day, latitude, sunshine_hours, altitude, tmax,
                                tmin, ws, z))
    logging.info('The reference et of the three days are %s.' % (eto))
    logging.info('Crop evaptr (ETC) is obtained based on land use patterns.')
    etc = []
    for i, j in enumerate(eto):
        etc.append(et.crop_et(out_lulc, j, i))
    logging.info('ETC in each step is recorded as shape file.')

    # from test.test_ab import hm_test as test
    # test()
    # logging.info('Finish test')


if __name__ == '__main__':
    main()
