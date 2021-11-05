"""
Calculate the crop evapotranspiration from the reference evaportranspiration
"""

import pyeto
import logging
import geopandas as gpd

def daily_eto(day_of_year, latitude, sunshine_hours, altitude, tmax, tmin, ws, z):
    """
    Estimate reference evapotranspiration (ETo) from a hypothetical short grass reference surface using the FAO-56 Penman-Monteith equation.
    Based on equation 6 in Allen et al (1998) and the code written by Mark Richards.

    Args:
        day_of_year (int): The day of year (1-365, 366)
        latitude (float): The latitude in degree sunshine_hours (float): The sunshine hours
        altitude (float): Elevation above the sea level
        tmax (float): Maximum temperature tmin (float): Minimum temperature ws (float): Wind speed at elevation (z) z (float): The elevation which measuring the wind speed
    Returns:
        ETo (float): Reference evapotranspiration [mm day-1] on the day of year
    """

    sol_dec = pyeto.sol_dec(day_of_year)
    latitude = pyeto.deg2rad(latitude)
    sha = pyeto.sunset_hour_angle(latitude, sol_dec)
    daylight_hours = pyeto.daylight_hours(sha)
    ird = pyeto.inv_rel_dist_earth_sun(day_of_year)
    et_rad = pyeto.et_rad(latitude, sol_dec, sha, ird)
    sol_rad = pyeto.sol_rad_from_sun_hours(daylight_hours, sunshine_hours, et_rad)
    net_in_sol_rad = pyeto.net_in_sol_rad(sol_rad, 0.23)
    cs_rad = pyeto.cs_rad(altitude, et_rad)
    avp = pyeto.avp_from_tmin(tmin)
    net_out_lw_rad = pyeto.net_out_lw_rad(tmin, tmax, sol_rad, cs_rad, avp)
    net_rad = pyeto.net_rad(net_in_sol_rad, net_out_lw_rad)
    t = (tmin+tmax)/2.0
    ws = pyeto.wind_speed_2m(ws, z)
    svp = pyeto.svp_from_t(t)
    delta_svp = pyeto.delta_svp(t)
    atm_pressure = pyeto.atm_pressure(altitude)
    psy = pyeto.psy_const(atm_pressure)

    eto = pyeto.fao56_penman_monteith(net_rad, t, ws, svp, avp, delta_svp, psy, 0.0)

    logging.info('Net radiation Rn is %s.' % (net_rad))
    logging.info('Air temperature at 2 m is %s.' % (t))
    logging.info('wind speed at 2m is %s.' % (ws))
    logging.info('Saturation vapor pressure (es) is %s.' % (svp))
    logging.info('Actual vapor pressure (ea) is %s.' % (avp))
    logging.info('Slope of saturation vapor pressure (Delta) is %s.' % (delta_svp))
    logging.info('Psychrometric constant is %s.' % (psy))
    logging.info('ET = %s.' % (eto))

    return eto

def crop_et(luc, eto, day):
    """
    Calculate the crop evapotranspiration/actual evapotranspiration based on the reference evapotranspiration

    ET0 is calculated from the pyeto package. 
    The land use land cover data is the discretized land form data with name of column 'CODE3'

    Args:
        luc (file): Discretized land form data
        eto (float): Daily average reference evapotranspiration

    Returns:
        etc (file): Discretized crop evapotranspiration in shape file format with column name 'etc'
    """
    gdf = gpd.read_file(luc)
    logging.info('The default kc coefficient for the crop evapotranspiration is 0.15.')
    gdf['kc'] = 0.15  # The default value of kc coefficient

    for i, luc in enumerate(gdf.CODE3):
        # print(luc, type(luc))
        if luc != None and int(luc) == 1:
            gdf.loc['kc', i] = 0.4
        elif luc != None and int(luc) == 2:
            gdf.loc['kc', i] = 1.0
        elif luc != None and int(luc) == 3:
            gdf.loc['kc', i] = 0.8
        elif luc != None and int(luc) == 4:
            gdf.loc['kc', i] = 0.4
        elif luc != None and int(luc) == 5:
            gdf.loc['kc', i] = 0.15
        elif luc != None and int(luc) == 6:
            gdf.loc['kc', i] = 0.2

    gdf.kc = gdf.kc * eto/24./6
    logging.info('The crop evapotranspiration of each time step (10 min/600 s) is saved as ../Test/Yachiyo_etc_%s.shp.' % (day))
    gdf.to_file('../Test/Yachiyo_etc_%s.shp' % (day))

    return gdf
