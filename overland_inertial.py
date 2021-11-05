"""
Floodplain simulation in Yachiyo area, Mobara City, Chiba, Japan
Hunter, et al. 2005; Bates, et al. 2010; de Almeida et al. 2012

Author: Mao Ouyang
Date:   2021-10-13
Email:  einooumo@hotmail.com
"""

import copy
import scipy as sp
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

cell_size = 10       # Cell size of the domain
N = 0.03             # Manning's coefficient
G = 9.81             # Gravity acceleration
ALPHA = 0.7          # Time step reduction coefficient
NT = 259200          # Total simulation period
#  NT = 159200          # Total simulation period
#  NT = 1900          # Total simulation period
DX = DY = cell_size  # Dimensions in the simulation
SDT = 60             # Saving interval

# Two-dimensional complex topography
COMP_TOPO = pd.read_csv('../Test/topo_rev_CellSize%d.csv' % (cell_size),
                        index_col=0)

# Initial condition
h0 = COMP_TOPO.to_numpy()

# Flood boundary condition / DT = 600 s / Measured river water level
MNT = 433
t = []
for i in range(MNT):
    t.append(i * 600)
ub_overland = pd.read_csv('../Test/ub_river.csv')
h_flood = ub_overland['Hayano'][0:MNT].astype(float)
h_flood = [i/100. for i in h_flood]  # Change unit from cm to m
h_flood = np.roll(h_flood, -65)  # Directly start from overtopping flow


def discharge(water_level, bed_level, water_disc,
              center_point, outer_point, manning_coef):
    '''
    Calculate the discharge from outer to the center

    Args:
        water_level (array):  water level of the previous step
        bed_level (array):    surface level of the study domain
        water_disc (array):   discharge of the previous step
        center_point (array): coordinate of the center point
        outer_point (array):  coordinate of the outer point
        manning_coef (float): manning's roughness coefficient

    Returns:
        Q (float): discharge to the center_point
    '''
    # flow direction
    if water_level[outer_point] >= water_level[center_point]:
        FLOW_DIRECTION = 1.
    else:
        FLOW_DIRECTION = -1.
    # effective flow depth
    hflow = max(water_level[center_point], water_level[outer_point]) - \
        max(bed_level[center_point], bed_level[outer_point])
    assert hflow >= 0., 'Effective flow depth cannot be negative'
    # Discharge based on inertial equation (semi-implicit)
    QINERT0 = water_disc[center_point]
    QINERT1 = G * hflow * DT * abs(water_level[center_point] -
                                   water_level[outer_point]) / DX
    if hflow ** (7./3.) > 0.:
        QINERT2 = G * DT * N ** 2. * abs(QINERT0) / (hflow ** (7./3.))
    else:
        QINERT2 = 0.
    Q = (QINERT0 + FLOW_DIRECTION * QINERT1) / (1. + QINERT2)
    return Q


def floodvolume(water_depth, dx, dy, time_interval):
    '''
    Calculate the water volume during the flood event

    Args:
        water_depth (float):   depth of water in each grid cell
        dx (float):            size in WE direction
        dy (float):            size in NS direction
        time_interval (float): time interval for this water depth

    Return:
        V (float):             flood volume in the time step
    '''
    V = water_depth * dx * dy * time_interval
    return V


# Numerical simulation (Bates et al. 2010)
# Fixed grid size and adaptive time step with inertial approximation
NX, NY = COMP_TOPO.values.shape  # Cell no. in X, Y directions
DT0 = 10.   # Initial time step
DT = DT0     # For the first step discharge calculation

qn_old = np.zeros((NX, NY))  # Previous discharge
qn_new = np.zeros((NX, NY))  # Updated dischange
hn_new = np.zeros((NX, NY))  # Updated water level
hn_old = copy.deepcopy(h0)   # Previous water level

# Main loop
TN = 0 + DT0  # Initial time step

# Record the adaptive time step evolution
adt_tn = []
adt_dt = []
adt_v = []
adt_tn.append(TN)
adt_dt.append(DT0)
adt_v.append(0.)

while TN <= NT:
    ##################################################################
    # First calculate the discharge (qn) from the previous time step #
    ##################################################################
    qn_we = np.zeros((NX, NY))  # calculated discharge in the WE direction
    qn_ns = np.zeros((NX, NY))  # calculated discharge in the NS direction
    # Discharge in WE direction
    for i, j in product(range(NX), range(1, NY)):
        cent = (i, j)
        outr = (i, j-1)
        qn_we[i, j] = discharge(hn_old, h0, qn_old, cent, outr, N)
    # Discharge in NS direction
    for i, j in product(range(NY), range(1, NX)):
        cent = (j, i)
        outr = (j-1, i)
        qn_ns[j, i] = discharge(hn_old, h0, qn_old, cent, outr, N)

    ###########################################################################
    # Calculate the water level according to the continue equation
    # Core body [0:-1, 0:-1]
    for i, j in product(range(NX-1), range(NY-1)):
        hn_new[i, j] = hn_old[i, j] + (qn_we[i, j] - qn_we[i, j+1] +
                                       qn_ns[i, j] - qn_ns[i+1, j]) *\
            DT / DX / DY
    # Corners
    hn_new[-1, -1] = hn_old[-1, -1] + (qn_ns[-1, -1] +
                                       qn_we[-1, -1]) * DT / DX / DY
    # Lines
    for i in range(NY-1):
        hn_new[-1, i] = hn_old[-1, i] + (qn_we[-1, i] - qn_we[-1, i+1] +
                                         qn_ns[-1, i]) * DT / DX / DY
    for i in range(NX-1):
        hn_new[i, -1] = hn_old[i, -1] + (qn_ns[i, -1] - qn_ns[i+1, -1] +
                                         qn_we[i, -1]) * DT / DX / DY

    ##################################################################
    #          Then specify the flood boundary condition             #
    ##################################################################
    FL = h_flood[int(TN/600)] + (TN % 600) / 600. *\
        (h_flood[int(TN/600)+1] - h_flood[int(TN/600)])
    depth = np.zeros((NX, NY))
    flood_volume = 0.
    for i, j in product(range(NX), range(NY)):
        if h0[i, j] == 0.:
            H_FLOOD = '{:0.5f}'.format(FL)
            hn_new[i, j] = H_FLOOD
            depth[i, j] = 0.0
        else:
            # Calculate the water level according to the continue equation
            hn_new[i, j] = '{:0.5f}'.format(hn_new[i, j])
            if hn_new[i, j] - h0[i, j] < 0.:
                hn_new[i, j] = h0[i, j]
            depth[i, j] = hn_new[i, j] - h0[i, j]
            flood_volume += floodvolume(depth[i, j], DX, DY, DT)

    # Update the new value to the old value for the next loop
    hn_old = copy.deepcopy(hn_new)
    qn_old = copy.deepcopy(qn_new)

    ##################################################################
    #              Save the data at time interval SDT                #
    ##################################################################
    if TN % SDT < DT0:
        T_SAVE = int(TN/SDT) * SDT
        H_SAVE = pd.DataFrame(hn_new-h0)
        Q_SAVE = pd.DataFrame(qn_new)
        H_SAVE.to_csv('../Test/Overland_inertial_H%d.csv' % (T_SAVE))
        #  Q_SAVE.to_csv('../Test/Overland_inertial_Q%d.csv' % (T_SAVE))

    ##################################################################
    #                 Calculate the next time step                   #
    ##################################################################
    if np.max(depth) == 0.:
        DT = DT0
    else:
        DT = ALPHA * DX / np.sqrt(G * np.max(depth))

    TN += DT
    adt_tn.append(TN)
    adt_dt.append(DT)
    adt_v.append(flood_volume)
    print(TN, FL, np.max(depth), DT)
    assert np.min(depth) >= 0., 'Inundation depth cannot be negative'
    assert np.max(depth) < 10., 'Inundation depth too high'

# Save some data
np.savetxt('../Test/Overland_TN.txt', adt_tn, fmt='%0.4f')
np.savetxt('../Test/Overland_V.txt', adt_v, fmt='%0.4f')
