"""
Floodplain flow simulation in Yachiyo area, Mobara City, Chiba, Japan
Hunter, et al. 2005; Bates, et al. 2010; de Almeida et al. 2012

Seperate the discharge calculation and the time step calculation

Author: Mao Ouyang
Date:   2021-07-21
Email:  einooumo@hotmail.com
"""

import copy
from itertools import product
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#  os.environ['PYTHONBREAKPOINT'] = '0'

plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

cell_size = 10
N = 0.02
NT = 259200
#  NT = 200
DX = DY = cell_size
HLIN = 0.010
SDT = 60  # The interval for saving the data

# Two-dimensional complex topography
COMP_TOPO = pd.read_csv('../Test/topo_rev_CellSize%d.csv' % (cell_size), index_col=0)

# Initial condition
h0 = COMP_TOPO.to_numpy()

# Flood boundary condition / DT = 600 s
MNT = 433
t = []
for i in range(MNT):
    t.append(i * 600)
ub_overland = pd.read_csv('../Test/ub_river.csv')
h_flood = ub_overland['Hayano'][0:MNT].astype(float)
h_flood = [i/100. for i in h_flood]  # Change unit from cm to m
h_flood = np.roll(h_flood, -65)


# To calculate the flood volume with sum of all positive values
def Positive_Sum(WaterLevel):
    s = 0.
    for i in WaterLevel:
        for j in i:
            if j > 0.:
                s += j
    return s


# Numerical simulation (Hunter et al. 2005)
DT = 10.   # First loop in the main
DT0 = 10.  # Initial time step
NX, NY = COMP_TOPO.values.shape  # Total number of cells in X, Y direction
q0 = h0 * 0.  # Initial discharge

# Main loop
TN = 0 + DT0  # Time step
adt_dt = []  # Record the time step
tn = []      # Record the time elapse
volume = []  # Record the flood volume

qn_new = np.zeros((NX, NY))
hn_new = np.zeros((NX, NY))
hn_old = copy.deepcopy(h0)
qn_old = copy.deepcopy(q0)
while TN <= NT:
    ###########################################################################
    # Four corners: 1. North-West
    if hn_old[0, 1] - hn_old[0, 0] >= 0.:
        FLOW_DIRECTION_X2 = 1.
    else:
        FLOW_DIRECTION_X2 = -1.
    if hn_old[0, 0] - hn_old[1, 0] >= 0.:
        FLOW_DIRECTION_Y2 = 1.
    else:
        FLOW_DIRECTION_Y2 = -1.
    # Manning's discharge
    hflow_x2 = abs(max(hn_old[0, 1], hn_old[0, 0]) -
                   max(h0[0, 1], h0[0, 0]))
    qn_x2 = hflow_x2 ** (5./3.) /\
        N * (abs(hn_old[0, 1] - hn_old[0, 0])/DX) ** (1./2.) * DY
    hflow_y2 = abs(max(hn_old[0, 0], hn_old[1, 0]) -
                   max(h0[0, 0], h0[1, 0]))
    qn_y2 = hflow_y2 ** (5./3.) /\
        N * (abs(hn_old[0, 0] - hn_old[1, 0])/DY) ** (1./2.) * DX
    # Linearized discharge
    lin_x2 = hflow_x2 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
        * (abs(hn_old[0, 1] - hn_old[0, 0]) / DX)
    lin_y2 = hflow_y2 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
        * (abs(hn_old[0, 0] - hn_old[1, 0]) / DY)
    # In X Direction
    if abs(hn_old[0, 1] - hn_old[0, 0]) <= HLIN:
        qn_x2 = abs(lin_x2) * FLOW_DIRECTION_X2
    else:
        qn_x2 = abs(qn_x2) * FLOW_DIRECTION_X2
    # In Y direction
    if abs(hn_old[0, 0] - hn_old[1, 0]) <= HLIN:
        qn_y2 = abs(lin_y2) * FLOW_DIRECTION_Y2
    else:
        qn_y2 = abs(qn_y2) * FLOW_DIRECTION_Y2
    # Summarize in the North-West
    qn_old[0, 0] = - qn_y2 + qn_x2

    ###########################################################################
    # Four corners: 2. South-West
    if hn_old[-1, 1] - hn_old[-1, 0] >= 0.:
        FLOW_DIRECTION_X2 = 1.
    else:
        FLOW_DIRECTION_X2 = -1.
    if hn_old[-2, 0] - hn_old[-1, 0] >= 0.:
        FLOW_DIRECTION_Y1 = 1.
    else:
        FLOW_DIRECTION_Y1 = -1.
    # Manning's discharge
    hflow_x2 = abs(max(hn_old[-1, 1], hn_old[-1, 0]) -
                   max(h0[-1, 1], h0[-1, 0]))
    qn_x2 = hflow_x2 ** (5./3.) /\
        N * (abs(hn_old[-1, 1] - hn_old[-1, 0])/DX) ** (1./2.) * DY
    hflow_y1 = abs(max(hn_old[-2, 0], hn_old[-1, 0]) -
                   max(h0[-2, 0], h0[-1, 0]))
    qn_y1 = hflow_y1 ** (5./3.) /\
        N * (abs(hn_old[-2, 0] - hn_old[-1, 0])/DY) ** (1./2.) * DX
    # Linearized discharge
    lin_x2 = hflow_x2 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
        * (abs(hn_old[-1, 1] - hn_old[-1, 0]) / DX)
    lin_y1 = hflow_y1 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
        * (abs(hn_old[-2, 0] - hn_old[-1, 0]) / DY)
    # In X Direction
    if abs(hn_old[-1, 1] - hn_old[-1, 0]) <= HLIN:
        qn_x2 = abs(lin_x2) * FLOW_DIRECTION_X2
    else:
        qn_x2 = abs(qn_x2) * FLOW_DIRECTION_X2
    # In Y direction
    if abs(hn_old[-2, 0] - hn_old[-1, 0]) <= HLIN:
        qn_y1 = abs(lin_y1) * FLOW_DIRECTION_Y1
    else:
        qn_y1 = abs(qn_y1) * FLOW_DIRECTION_Y1
    # Summarize in the South-West
    qn_old[-1, 0] = qn_y1 + qn_x2

    ###########################################################################
    # Four corners: 3. South-East
    if hn_old[-1, -1] - hn_old[-1, -2] >= 0.:
        FLOW_DIRECTION_X1 = 1.
    else:
        FLOW_DIRECTION_X1 = -1.
    if hn_old[-2, -1] - hn_old[-1, -1] >= 0.:
        FLOW_DIRECTION_Y1 = 1.
    else:
        FLOW_DIRECTION_Y1 = -1.
    # Manning's discharge
    hflow_x1 = abs(max(hn_old[-1, -1], hn_old[-1, -2]) -
                   max(h0[-1, -1], h0[-1, -2]))
    qn_x1 = hflow_x1 ** (5./3.) /\
        N * (abs(hn_old[-1, -1] - hn_old[-1, -2])/DX) ** (1./2.) * DY
    hflow_y1 = abs(max(hn_old[-2, -1], hn_old[-1, -1]) -
                   max(h0[-2, -1], h0[-1, -1]))
    qn_y1 = hflow_y1 ** (5./3.) /\
        N * (abs(hn_old[-2, -1] - hn_old[-1, -1])/DY) ** (1./2.) * DX
    # Linearized discharge
    lin_x1 = hflow_x1 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
        * (abs(hn_old[-1, -1] - hn_old[-1, -2]) / DX)
    lin_y1 = hflow_y1 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
        * (abs(hn_old[-2, -1] - hn_old[-1, -1]) / DY)
    # In X Direction
    if abs(hn_old[-1, -1] - hn_old[-1, -2]) <= HLIN:
        qn_x1 = abs(lin_x1) * FLOW_DIRECTION_X1
    else:
        qn_x1 = abs(qn_x1) * FLOW_DIRECTION_X1
    # In Y direction
    if abs(hn_old[-2, -1] - hn_old[-1, -1]) <= HLIN:
        qn_y1 = abs(lin_y1) * FLOW_DIRECTION_Y1
    else:
        qn_y1 = abs(qn_y1) * FLOW_DIRECTION_Y1
    # Summarize in the South-East
    qn_old[-1, -1] = qn_y1 - qn_x1

    ###########################################################################
    # Four corners: 4. North-East
    if hn_old[0, -1] - hn_old[0, -2] >= 0.:
        FLOW_DIRECTION_X1 = 1.
    else:
        FLOW_DIRECTION_X1 = -1.
    if hn_old[0, -1] - hn_old[1, -1] >= 0.:
        FLOW_DIRECTION_Y2 = 1.
    else:
        FLOW_DIRECTION_Y2 = -1.

    # Manning's discharge
    hflow_x1 = abs(max(hn_old[0, -1], hn_old[0, -2]) -
                   max(h0[0, -1], h0[0, -2]))
    qn_x1 = hflow_x1 ** (5./3.) /\
        N * (abs(hn_old[0, -1] - hn_old[0, -2])/DX) ** (1./2.) * DY

    hflow_y2 = abs(max(hn_old[0, -1], hn_old[1, -1]) -
                   max(h0[0, -1], h0[1, -1]))
    qn_y2 = hflow_y2 ** (5./3.) /\
        N * (abs(hn_old[0, -1] - hn_old[1, -1])/DY) ** (1./2.) * DX

    # Linearized discharge
    lin_x1 = hflow_x1 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
        * (abs(hn_old[0, -1] - hn_old[0, -2]) / DX)
    lin_y2 = hflow_y2 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
        * (abs(hn_old[0, -1] - hn_old[1, -1]) / DY)
    # In X Direction
    if abs(hn_old[0, -1] - hn_old[0, -2]) <= HLIN:
        qn_x1 = abs(lin_x1) * FLOW_DIRECTION_X1
    else:
        qn_x1 = abs(qn_x1) * FLOW_DIRECTION_X1
    # In Y direction
    if abs(hn_old[0, -1] - hn_old[1, -1]) <= HLIN:
        qn_y2 = abs(lin_y2) * FLOW_DIRECTION_Y2
    else:
        qn_y2 = abs(qn_y2) * FLOW_DIRECTION_Y2
    # Summarize in the North-East
    qn_old[0, -1] = - qn_y2 - qn_x1

    ###########################################################################
    # Fourt lines: 1. North [1:-1]
    for i in range(1, NY-1):
        # Calculate the flow directions
        if hn_old[0, i] - hn_old[1, i] >= 0.:
            FLOW_DIRECTION_Y2 = 1.
        else:
            FLOW_DIRECTION_Y2 = -1.
        if hn_old[0, i] - hn_old[0, i-1] >= 0.:
            FLOW_DIRECTION_X1 = 1.
        else:
            FLOW_DIRECTION_X1 = -1.
        if hn_old[0, i+1] - hn_old[0, i] >= 0.:
            FLOW_DIRECTION_X2 = 1.
        else:
            FLOW_DIRECTION_X2 = -1.
        # Manning's discharge
        hflow_y2 = abs(max(hn_old[0, i], hn_old[1, i]) -
                       max(h0[0, i], h0[1, i]))
        qn_y2 = hflow_y2 ** (5./3.) /\
            N * (abs(hn_old[0, i] - hn_old[1, i])/DY) ** (1./2.) * DX

        hflow_x1 = abs(max(hn_old[0, i], hn_old[0, i-1]) -
                       max(h0[0, i], h0[0, i-1]))
        qn_x1 = hflow_x1 ** (5./3.) /\
            N * (abs(hn_old[0, i] - hn_old[0, i-1])/DX) ** (1./2.) * DY

        hflow_x2 = abs(max(hn_old[0, i+1], hn_old[0, i]) -
                       max(h0[0, i+1], h0[0, i]))
        qn_x2 = hflow_x2 ** (5./3.) /\
            N * (abs(hn_old[0, i+1] - hn_old[0, i])/DX) ** (1./2.) * DY
        # Linearized discharge
        lin_y2 = hflow_y2 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
            * (abs(hn_old[0, i] - hn_old[1, i]) / DY)
        lin_x1 = hflow_x1 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
            * (abs(hn_old[0, i] - hn_old[0, i-1]) / DX)
        lin_x2 = hflow_x2 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
            * (abs(hn_old[0, i+1] - hn_old[0, i]) / DX)
        # In Y direction
        if abs(hn_old[0, i] - hn_old[1, i]) <= HLIN:
            qn_y2 = abs(lin_y2) * FLOW_DIRECTION_Y2
        else:
            qn_y2 = abs(qn_y2) * FLOW_DIRECTION_Y2
        # In X direction
        if abs(hn_old[0, i] - hn_old[0, i-1]) <= HLIN:
            qn_x1 = abs(lin_x1) * FLOW_DIRECTION_X1
        else:
            qn_x1 = abs(qn_x1) * FLOW_DIRECTION_X1
        if abs(hn_old[0, i+1] - hn_old[0, i]) <= HLIN:
            qn_x2 = abs(lin_x2) * FLOW_DIRECTION_X2
        else:
            qn_x2 = abs(qn_x2) * FLOW_DIRECTION_X2
        # Summarized in the North [1:-1]
        qn_old[0, i] = -qn_y2 + qn_x2 - qn_x1

        #######################################################################
        # Four lines: 2. South [1:-1]
        # Calculate the flow directions
        if hn_old[-2, i] - hn_old[-1, i] >= 0.:
            FLOW_DIRECTION_Y1 = 1.
        else:
            FLOW_DIRECTION_Y1 = -1.
        if hn_old[-1, i] - hn_old[-1, i-1] >= 0.:
            FLOW_DIRECTION_X1 = 1.
        else:
            FLOW_DIRECTION_X1 = -1.
        if hn_old[-1, i+1] - hn_old[-1, i] >= 0.:
            FLOW_DIRECTION_X2 = 1.
        else:
            FLOW_DIRECTION_X2 = -1.
        # Manning's discharge
        hflow_y1 = abs(max(hn_old[-2, i], hn_old[-1, i]) -
                       max(h0[-2, i], h0[-1, i]))
        qn_y1 = hflow_y1 ** (5./3.) /\
            N * (abs(hn_old[-2, i] - hn_old[-1, i])/DY) ** (1./2.) * DX

        hflow_x1 = abs(max(hn_old[-1, i], hn_old[-1, i-1]) -
                       max(h0[-1, i], h0[-1, i-1]))
        qn_x1 = hflow_x1 ** (5./3.) /\
            N * (abs(hn_old[-1, i] - hn_old[-1, i-1])/DX) ** (1./2.) * DY

        hflow_x2 = abs(max(hn_old[-1, i+1], hn_old[-1, i]) -
                       max(h0[-1, i+1], h0[-1, i]))
        qn_x2 = hflow_x2 ** (5./3.) /\
            N * (abs(hn_old[-1, i+1] - hn_old[-1, i])/DX) ** (1./2.) * DY
        # Linearized discharge
        lin_y1 = hflow_y1 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
            * (abs(hn_old[-2, i] - hn_old[-1, i]) / DY)
        lin_x1 = hflow_x1 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
            * (abs(hn_old[-1, i] - hn_old[-1, i-1]) / DX)
        lin_x2 = hflow_x2 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
            * (abs(hn_old[-1, i+1] - hn_old[-1, i]) / DX)
        # In Y direction
        if abs(hn_old[-2, i] - hn_old[-1, i]) <= HLIN:
            qn_y1 = abs(lin_y1) * FLOW_DIRECTION_Y1
        else:
            qn_y1 = abs(qn_y1) * FLOW_DIRECTION_Y1
        # In X direction
        if abs(hn_old[-1, i] - hn_old[-1, i-1]) <= HLIN:
            qn_x1 = abs(lin_x1) * FLOW_DIRECTION_X1
        else:
            qn_x1 = abs(qn_x1) * FLOW_DIRECTION_X1
        if abs(hn_old[-1, i+1] - hn_old[-1, i]) <= HLIN:
            qn_x2 = abs(lin_x2) * FLOW_DIRECTION_X2
        else:
            qn_x2 = abs(qn_x2) * FLOW_DIRECTION_X2
        # Summarized in the South [1:-1]
        qn_old[-1, i] = qn_y1 + qn_x2 - qn_x1

    ###########################################################################
    # Four lines: 3. West [1:-1]
    for i in range(1, NX-1):
        # Calculate the flow directions
        if hn_old[i, 1] - hn_old[i, 0] >= 0.:
            FLOW_DIRECTION_X2 = 1.
        else:
            FLOW_DIRECTION_X2 = -1.
        if hn_old[i-1, 0] - hn_old[i, 0] >= 0.:
            FLOW_DIRECTION_Y1 = 1.
        else:
            FLOW_DIRECTION_Y1 = -1.
        if hn_old[i, 0] - hn_old[i+1, 0] >= 0.:
            FLOW_DIRECTION_Y2 = 1.
        else:
            FLOW_DIRECTION_Y2 = -1.
        # Manning's discharge
        hflow_x2 = abs(max(hn_old[i, 1], hn_old[i, 0]) -
                       max(h0[i, 1], h0[i, 0]))
        qn_x2 = hflow_x2 ** (5./3.) /\
            N * (abs(hn_old[i, 1] - hn_old[i, 0])/DX) ** (1./2.) * DY

        hflow_y1 = abs(max(hn_old[i-1, 0], hn_old[i, 0]) -
                       max(h0[i-1, 0], h0[i, 0]))
        qn_y1 = hflow_y1 ** (5./3.) /\
            N * (abs(hn_old[i-1, 0] - hn_old[i, 0])/DY) ** (1./2.) * DX

        hflow_y2 = abs(max(hn_old[i, 0], hn_old[i+1, 0]) -
                       max(h0[i, 0], h0[i+1, 0]))
        qn_y2 = hflow_y2 ** (5./3.) /\
            N * (abs(hn_old[i, 0] - hn_old[i+1, 0])/DY) ** (1./2.) * DX
        # Linearized discharge
        lin_x2 = hflow_x2 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
            * (abs(hn_old[i, 1] - hn_old[i, 0]) / DX)
        lin_y1 = hflow_y1 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
            * (abs(hn_old[i-1, 0] - hn_old[i, 0]) / DY)
        lin_y2 = hflow_y2 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
            * (abs(hn_old[i, 0] - hn_old[i+1, 0]) / DY)
        # In X direction
        if abs(hn_old[i, 1] - hn_old[i, 0]) <= HLIN:
            qn_x2 = abs(lin_x2) * FLOW_DIRECTION_X2
        else:
            qn_x2 = abs(qn_x2) * FLOW_DIRECTION_X2
        # In Y direction
        if abs(hn_old[i-1, 0] - hn_old[i, 0]) <= HLIN:
            qn_y1 = abs(lin_y1) * FLOW_DIRECTION_Y1
        else:
            qn_y1 = abs(qn_y1) * FLOW_DIRECTION_Y1
        if abs(hn_old[i, 0] - hn_old[i+1, 0]) <= HLIN:
            qn_y2 = abs(lin_y2) * FLOW_DIRECTION_Y2
        else:
            qn_y2 = abs(qn_y2) * FLOW_DIRECTION_Y2
        # Summarized in the West [1:-1]
        qn_old[i, 0] = qn_y1 - qn_y2 + qn_x2

        #######################################################################
        # Four lines: 4. East [1:-1]
        # Calculate the flow directions
        if hn_old[i, -1] - hn_old[i, -2] >= 0.:
            FLOW_DIRECTION_X1 = 1.
        else:
            FLOW_DIRECTION_X1 = -1.
        if hn_old[i-1, -1] - hn_old[i, -1] >= 0.:
            FLOW_DIRECTION_Y1 = 1.
        else:
            FLOW_DIRECTION_Y1 = -1.
        if hn_old[i, -1] - hn_old[i+1, -1] >= 0.:
            FLOW_DIRECTION_Y2 = 1.
        else:
            FLOW_DIRECTION_Y2 = -1.
        # Manning's discharge
        hflow_x1 = abs(max(hn_old[i, -1], hn_old[i, -2]) -
                       max(h0[i, -1], h0[i, -2]))
        qn_x1 = hflow_x1 ** (5./3.) /\
            N * (abs(hn_old[i, -1] - hn_old[i, -2])/DX) ** (1./2.) * DY

        hflow_y1 = abs(max(hn_old[i-1, -1], hn_old[i, -1]) -
                       max(h0[i-1, -1], h0[i, -1]))
        qn_y1 = hflow_y1 ** (5./3.) /\
            N * (abs(hn_old[i-1, -1] - hn_old[i, -1])/DY) ** (1./2.) * DX

        hflow_y2 = abs(max(hn_old[i, -1], hn_old[i+1, -1]) -
                       max(h0[i, -1], h0[i+1, -1]))
        qn_y2 = hflow_y2 ** (5./3.) /\
            N * (abs(hn_old[i, -1] - hn_old[i+1, -1])/DY) ** (1./2.) * DX
        # Linearized discharge
        lin_x1 = hflow_x1 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
            * (abs(hn_old[i, -1] - hn_old[i, -2]) / DX)
        lin_y1 = hflow_y1 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
            * (abs(hn_old[i-1, -1] - hn_old[i, -1]) / DY)
        lin_y2 = hflow_y2 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
            * (abs(hn_old[i, -1] - hn_old[i+1, -1]) / DY)
        # In X direction
        if abs(hn_old[i, -1] - hn_old[i, -2]) <= HLIN:
            qn_x1 = abs(lin_x1) * FLOW_DIRECTION_X1
        else:
            qn_x1 = abs(qn_x1) * FLOW_DIRECTION_X1
        # In Y direction
        if abs(hn_old[i-1, -1] - hn_old[i, -1]) <= HLIN:
            qn_y1 = abs(lin_y1) * FLOW_DIRECTION_Y1
        else:
            qn_y1 = abs(qn_y1) * FLOW_DIRECTION_Y1
        if abs(hn_old[i, -1] - hn_old[i+1, -1]) <= HLIN:
            qn_y2 = abs(lin_y2) * FLOW_DIRECTION_Y2
        else:
            qn_y2 = abs(qn_y2) * FLOW_DIRECTION_Y2
        # Summarized in the East [1:-1]
        qn_old[i, -1] = qn_y1 - qn_y2 - qn_x1

    ###########################################################################
    # Core body [1:-1, 1:-1]
    for i, j in product(range(1, NX-1), range(1, NY-1)):
        # First calculate the discharge (qn) from the previous time step
        # 1.1 Comparison between the water levels in the neighour cells to
        #     determine flow directions
        if hn_old[i, j] - hn_old[i-1, j] >= 0.:
            FLOW_DIRECTION_X1 = 1.
        else:
            FLOW_DIRECTION_X1 = -1.
        if hn_old[i+1, j] - hn_old[i, j] >= 0.:
            FLOW_DIRECTION_X2 = 1.
        else:
            FLOW_DIRECTION_X2 = -1.
        if hn_old[i, j-1] - hn_old[i, j] >= 0.:
            FLOW_DIRECTION_Y1 = 1.
        else:
            FLOW_DIRECTION_Y1 = -1.
        if hn_old[i, j] - hn_old[i, j+1] >= 0.:
            FLOW_DIRECTION_Y2 = 1.
        else:
            FLOW_DIRECTION_Y2 = -1.
        # 1.2 Calculate the discharges in West (X1), East (X2), North (Y1),
        # and South (Y2) directions
        hflow_x1 = abs(max(hn_old[i, j], hn_old[i-1, j]) -
                       max(h0[i, j], h0[i-1, j]))
        hflow_x2 = abs(max(hn_old[i, j], hn_old[i+1, j]) -
                       max(h0[i, j], h0[i+1, j]))
        hflow_y1 = abs(max(hn_old[i, j], hn_old[i, j-1]) -
                       max(h0[i, j], h0[i, j-1]))
        hflow_y2 = abs(max(hn_old[i, j], hn_old[i, j+1]) -
                       max(h0[i, j], h0[i, j+1]))
        qn_x1 = FLOW_DIRECTION_X1 * hflow_x1 ** (5./3.) /\
            N * (abs(hn_old[i, j] - hn_old[i-1, j])/DX) ** (1./2.) * DY
        qn_x2 = FLOW_DIRECTION_X2 * hflow_x2 ** (5./3.) /\
            N * (abs(hn_old[i, j] - hn_old[i+1, j])/DX) ** (1./2.) * DY
        qn_y1 = FLOW_DIRECTION_Y1 * hflow_y1 ** (5./3.) /\
            N * (abs(hn_old[i, j] - hn_old[i, j-1])/DX) ** (1./2.) * DY
        qn_y2 = FLOW_DIRECTION_Y2 * hflow_y2 ** (5./3.) /\
            N * (abs(hn_old[i, j] - hn_old[i, j+1])/DX) ** (1./2.) * DY
        # 1.3 Employ linearized discharge in West (X1), East (X2), North (Y1),
        # and South (Y2) directions
        lin_x1 = hflow_x1 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
            * (abs(hn_old[i, j] - hn_old[i-1, j]) / DX)
        lin_x2 = hflow_x2 ** (5./3.) / N * (DX/HLIN) ** (1./2.)\
            * (abs(hn_old[i, j] - hn_old[i+1, j]) / DX)
        lin_y1 = hflow_y1 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
            * (abs(hn_old[i, j] - hn_old[i, j-1]) / DY)
        lin_y2 = hflow_y2 ** (5./3.) / N * (DY/HLIN) ** (1./2.)\
            * (abs(hn_old[i, j] - hn_old[i, j+1]) / DY)
        # 1.4 Comparison with the threshold to obtain the dischage and
        # adaptive time step
        if abs(hn_old[i, j] - hn_old[i-1, j]) <= HLIN:
            qn_x1 = abs(lin_x1) * FLOW_DIRECTION_X1
        else:
            qn_x1 = abs(qn_x1) * FLOW_DIRECTION_X1
        if abs(hn_old[i+1, j] - hn_old[i, j]) <= HLIN:
            qn_x2 = abs(lin_x2) * FLOW_DIRECTION_X2
        else:
            qn_x2 = abs(qn_x2) * FLOW_DIRECTION_X2
        if abs(hn_old[i, j-1] - hn_old[i, j]) <= HLIN:
            qn_y1 = abs(lin_y1) * FLOW_DIRECTION_Y1
        else:
            qn_y1 = abs(qn_y1) * FLOW_DIRECTION_Y1
        if abs(hn_old[i, j] - hn_old[i, j+1]) <= HLIN:
            qn_y2 = abs(lin_y2) * FLOW_DIRECTION_Y2
        else:
            qn_y2 = abs(qn_y2) * FLOW_DIRECTION_Y2
        # 1.5 Mark the cells of river
        if h0[i, j] == 0.:
            qn_x1 = qn_x2 = qn_y1 = qn_y2 = 0.
        # 1.6 Summarize all the discharge and time steps
        # The difference of discharge is directly stored in the list
        qn_old[i, j] = qn_y1 - qn_y2 + qn_x2 - qn_x1

    # Then specify the flood boundary condition
    FL = h_flood[int(TN/600)] + (TN % 600) / 600. *\
        (h_flood[int(TN/600)+1] - h_flood[int(TN/600)])
    for i, j in product(range(NX), range(NY)):
        if h0[i, j] == 0.:
            H_FLOOD = '{:0.5f}'.format(FL)
            hn_new[i, j] = H_FLOOD
        else:
            # Calculate the water level according to the continue equation
            HN_NEW = hn_old[i, j] + qn_old[i, j] * DT / DX / DY
            hn_new[i, j] = '{:0.5f}'.format(HN_NEW)
    # Update the water level in both corners, lines, and core body
    #  for i, j in product(range(NX), range(NY)):
    #      # Calculate the water level according to the continue equation
    #      hn_new[i, j] = hn_old[i, j] + qn_old[i, j] * DT / DX / DY
    print(TN, FL)

    ###########################################################################
    ###########################################################################
    # Save the data at time interval SDT
    if TN % SDT < DT0:
        T_SAVE = int(TN/SDT) * SDT
        H_SAVE = pd.DataFrame(hn_new-h0)
        H_SAVE.to_csv('../Test/Overland_%d.csv' % (T_SAVE))

    ###########################################################################
    ###########################################################################
    # Calculate the adaptive time step
    adt = []
    adt.append(DT0)  # Ensure the min() can obtain a value
    ###########################################################################
    # Four corners: 1. North-West
    hflow_x2 = abs(max(hn_old[0, 1], hn_old[0, 0]) -
                   max(h0[0, 1], h0[0, 0]))
    hflow_y2 = abs(max(hn_old[0, 0], hn_old[1, 0]) -
                   max(h0[0, 0], h0[1, 0]))
    # In X direction
    if hflow_x2 == 0.:
        dt_x2 = DT0
    elif abs(hn_old[0, 1] - hn_old[0, 0]) <= HLIN and hflow_x2 != 0.:
        dt_x2 = DX ** 2 / 4. * N / hflow_x2 ** (5./3.) * (HLIN/DX) ** (1./2.)
    elif abs(hn_old[0, 1] - hn_old[0, 0]) > HLIN and hflow_x2 != 0.:
        dt_x2 = DX ** 2 / 4. * 2. * N / hflow_x2 ** (5./3.)\
            * (abs(hn_old[0, 1] - hn_old[0, 0]) / DX) ** (1./2.)
    # In Y direction
    if hflow_y2 == 0.:
        dt_y2 = DT0
    elif abs(hn_old[0, 0] - hn_old[1, 0]) <= HLIN and hflow_y2 != 0.:
        dt_y2 = DX ** 2 / 4. * N / hflow_y2 ** (5./3.) * (HLIN/DY) ** (1./2.)
    elif abs(hn_old[0, 0] - hn_old[1, 0]) > HLIN and hflow_y2 != 0.:
        dt_y2 = DX ** 2 / 4. * 2. * N / hflow_y2 ** (5./3.)\
            * (abs(hn_old[0, 0] - hn_old[1, 0]) / DY) ** (1./2.)
    adt.append(min(dt_x2, dt_y2))
    #  if min(adt) < DT0:
    #      breakpoint()

    ###########################################################################
    # Four concers: 2. South-West
    hflow_x2 = abs(max(hn_old[-1, 1], hn_old[-1, 0]) -
                   max(h0[-1, 1], h0[-1, 0]))
    hflow_y1 = abs(max(hn_old[-2, 0], hn_old[-1, 0]) -
                   max(h0[-2, 0], h0[-1, 0]))
    # In X direction
    if hflow_x2 == 0.:
        dt_x2 = DT0
    elif abs(hn_old[-1, 1] - hn_old[-1, 0]) <= HLIN and hflow_x2 != 0.:
        dt_x2 = DX ** 2 / 4. * N / hflow_x2 ** (5./3.) * (HLIN/DX) ** (1./2.)
    elif abs(hn_old[-1, 1] - hn_old[-1, 0]) > HLIN and hflow_x2 != 0.:
        dt_x2 = DX ** 2 / 4. * 2. * N / hflow_x2 ** (5./3.)\
            * (abs(hn_old[-1, 1] - hn_old[-1, 0]) / DX) ** (1./2.)
    # In Y direction
    if hflow_y1 == 0.:
        dt_y1 = DT0
    elif abs(hn_old[-1, 0] - hn_old[-2, 0]) <= HLIN and hflow_y1 != 0.:
        dt_y1 = DX ** 2 / 4. * N / hflow_y1 ** (5./3.) * (HLIN/DY) ** (1./2.)
    elif abs(hn_old[-1, 0] - hn_old[-2, 0]) > HLIN and hflow_y1 != 0.:
        dt_y1 = DX ** 2 / 4. * 2. * N / hflow_y1 ** (5./3.)\
            * (abs(hn_old[-1, 0] - hn_old[-2, 0]) / DY) ** (1./2.)
    adt.append(min(dt_x2, dt_y1))
    #  if min(adt) < DT0:
    #      breakpoint()

    ###########################################################################
    # Four concers: 3. South-East
    hflow_x1 = abs(max(hn_old[-1, -1], hn_old[-1, -2]) -
                   max(h0[-1, -1], h0[-1, -2]))
    hflow_y1 = abs(max(hn_old[-2, -1], hn_old[-1, -1]) -
                   max(h0[-2, -1], h0[-1, -1]))
    # In X direction
    if hflow_x1 == 0.:
        dt_x1 = DT0
    elif abs(hn_old[-1, -1] - hn_old[-1, -2]) <= HLIN and hflow_x1 != 0.:
        dt_x1 = DX ** 2 / 4. * N / hflow_x1 ** (5./3.) * (HLIN/DX) ** (1./2.)
    elif abs(hn_old[-1, -1] - hn_old[-1, -2]) > HLIN and hflow_x1 != 0.:
        dt_x1 = DX ** 2 / 4. * 2. * N / hflow_x1 ** (5./3.)\
            * (abs(hn_old[-1, -1] - hn_old[-1, -2]) / DX) ** (1./2.)
    # In Y direction
    if hflow_y1 == 0.:
        dt_y1 = DT0
    elif abs(hn_old[-1, -1] - hn_old[-2, -1]) <= HLIN and hflow_y1 != 0.:
        dt_y1 = DX ** 2 / 4. * N / hflow_y1 ** (5./3.) * (HLIN/DY) ** (1./2.)
    elif abs(hn_old[-1, -1] - hn_old[-2, -1]) > HLIN and hflow_y1 != 0.:
        dt_y1 = DX ** 2 / 4. * 2. * N / hflow_y1 ** (5./3.)\
            * (abs(hn_old[-1, -1] - hn_old[-2, -1]) / DY) ** (1./2.)
    adt.append(min(dt_x1, dt_y1))
    #  if min(adt) < DT0:
    #      breakpoint()

    ###########################################################################
    # Four concers: 4. North-East
    hflow_x1 = abs(max(hn_old[0, -1], hn_old[0, -2]) -
                   max(h0[0, -1], h0[0, -2]))
    hflow_y2 = abs(max(hn_old[0, -1], hn_old[1, -1]) -
                   max(h0[0, -1], h0[1, -1]))
    # In X direction
    if hflow_x1 == 0.:
        dt_x1 = DT0
    elif abs(hn_old[0, -1] - hn_old[0, -2]) <= HLIN and hflow_x1 != 0.:
        dt_x1 = DX ** 2 / 4. * N / hflow_x1 ** (5./3.) * (HLIN/DX) ** (1./2.)
    elif abs(hn_old[0, -1] - hn_old[0, -2]) > HLIN and hflow_x1 != 0.:
        dt_x1 = DX ** 2 / 4. * 2. * N / hflow_x1 ** (5./3.)\
            * (abs(hn_old[0, -1] - hn_old[0, -2]) / DX) ** (1./2.)
    # In Y direction
    if hflow_y2 == 0.:
        dt_y2 = DT0
    elif abs(hn_old[0, -1] - hn_old[1, -1]) <= HLIN and hflow_y2 != 0.:
        dt_y2 = DX ** 2 / 4. * N / hflow_y2 ** (5./3.) * (HLIN/DY) ** (1./2.)
    elif abs(hn_old[0, -1] - hn_old[1, -1]) > HLIN and hflow_y2 != 0.:
        dt_y2 = DX ** 2 / 4. * 2. * N / hflow_y2 ** (5./3.)\
            * (abs(hn_old[0, -1] - hn_old[1, -1]) / DY) ** (1./2.)
    adt.append(min(dt_x1, dt_y2))
    #  if min(adt) < DT0:
    #      breakpoint()

    ###########################################################################
    # Four lines: 1. North [1:-1]
    for i in range(1, NY-1):
        hflow_y2 = abs(max(hn_old[0, i], hn_old[1, i]) -
                       max(h0[0, i], h0[1, i]))

        hflow_x1 = abs(max(hn_old[0, i], hn_old[0, i-1]) -
                       max(h0[0, i], h0[0, i-1]))

        hflow_x2 = abs(max(hn_old[0, i+1], hn_old[0, i]) -
                       max(h0[0, i+1], h0[0, i]))
        # In Y direction
        if hflow_y2 == 0.:
            dt_y2 = DT0
        elif abs(hn_old[0, i] - hn_old[1, i]) <= HLIN and hflow_y2 != 0.:
            dt_y2 = DX ** 2 / 4. * N / hflow_y2 ** (5./3.) *\
                (HLIN/DY) ** (1./2.)
        elif abs(hn_old[0, i] - hn_old[1, i]) > HLIN and hflow_y2 != 0.:
            dt_y2 = DX ** 2 / 4. * 2. * N / hflow_y2 ** (5./3.)\
                * (abs(hn_old[0, i] - hn_old[1, i]) / DY) ** (1./2.)
        # In X direction
        if hflow_x1 == 0.:
            dt_x1 = DT0
        elif abs(hn_old[0, i] - hn_old[0, i-1]) <= HLIN and hflow_x1 != 0.:
            dt_x1 = DX ** 2 / 4. * N / hflow_x1 ** (5./3.) *\
                (HLIN/DX) ** (1./2.)
        elif abs(hn_old[0, i] - hn_old[0, i-1]) > HLIN and hflow_x1 != 0.:
            dt_x1 = DX ** 2 / 4. * 2. * N / hflow_x1 ** (5./3.)\
                * (abs(hn_old[0, i] - hn_old[0, i-1]) / DX) ** (1./2.)

        if hflow_x2 == 0.:
            dt_x2 = DT0
        elif abs(hn_old[0, i+1] - hn_old[0, i]) <= HLIN and hflow_x2 != 0.:
            dt_x2 = DX ** 2 / 4. * N / hflow_x2 ** (5./3.) *\
                (HLIN/DX) ** (1./2.)
        elif abs(hn_old[0, i+1] - hn_old[0, i]) > HLIN and hflow_x2 != 0.:
            dt_x2 = DX ** 2 / 4. * 2. * N / hflow_x2 ** (5./3.)\
                * (abs(hn_old[0, i+1] - hn_old[0, i]) / DX) ** (1./2.)

        #  # Boundary cells connected with Rivers
        #  if h0[1, i] == 0.:
        #      dt_y2 = DT0
        adt.append(min(dt_y2, dt_x2, dt_x1))
        #  if min(adt) < DT0:
        #      breakpoint()

        #######################################################################
        # Four lines: 2. South [1:-1]
        hflow_y1 = abs(max(hn_old[-2, i], hn_old[-1, i]) -
                       max(h0[-2, i], h0[-1, i]))

        hflow_x1 = abs(max(hn_old[-1, i], hn_old[-1, i-1]) -
                       max(h0[-1, i], h0[-1, i-1]))

        hflow_x2 = abs(max(hn_old[-1, i+1], hn_old[-1, i]) -
                       max(h0[-1, i+1], h0[-1, i]))

        # In Y direction
        if hflow_y1 == 0.:
            dt_y1 = DT0
        elif abs(hn_old[-2, i] - hn_old[-1, i]) <= HLIN and hflow_y1 != 0.:
            dt_y1 = DX ** 2 / 4. * N / hflow_y1 ** (5./3.) *\
                (HLIN/DY) ** (1./2.)
        elif abs(hn_old[-2, i] - hn_old[-1, i]) > HLIN and hflow_y1 != 0.:
            dt_y1 = DX ** 2 / 4. * 2. * N / hflow_y1 ** (5./3.)\
                * (abs(hn_old[-2, i] - hn_old[-1, i]) / DY) ** (1./2.)
        # In X direction
        if hflow_x1 == 0.:
            dt_x1 = DT0
        elif abs(hn_old[-1, i] - hn_old[-1, i-1]) <= HLIN and hflow_x1 != 0.:
            dt_x1 = DX ** 2 / 4. * N / hflow_x1 ** (5./3.) *\
                (HLIN/DX) ** (1./2.)
        elif abs(hn_old[-1, i] - hn_old[-1, i-1]) > HLIN and hflow_x1 != 0.:
            dt_x1 = DX ** 2 / 4. * 2. * N / hflow_x1 ** (5./3.)\
                * (abs(hn_old[-1, i] - hn_old[-1, i-1]) / DX) ** (1./2.)

        if hflow_x2 == 0.:
            dt_x2 = DT0
        elif abs(hn_old[-1, i+1] - hn_old[-1, i]) <= HLIN and hflow_x2 != 0.:
            dt_x2 = DX ** 2 / 4. * N / hflow_x2 ** (5./3.) *\
                (HLIN/DX) ** (1./2.)
        elif abs(hn_old[-1, i+1] - hn_old[-1, i]) > HLIN and hflow_x2 != 0.:
            dt_x2 = DX ** 2 / 4. * 2. * N / hflow_x2 ** (5./3.)\
                * (abs(hn_old[-1, i+1] - hn_old[-1, i]) / DX) ** (1./2.)

        #  # Boundary cells connected with Rivers
        #  if h0[-2, i] == 0.:
        #      dt_y1 = DT0
        adt.append(min(dt_y1, dt_x2, dt_x1))
        #  if min(adt) < DT0:
        #      breakpoint()

    ###########################################################################
    # Four lines: 3. West [1:-1]
    for i in range(1, NX-1):
        hflow_x2 = abs(max(hn_old[i, 1], hn_old[i, 0]) -
                       max(h0[i, 1], h0[i, 0]))

        hflow_y1 = abs(max(hn_old[i-1, 0], hn_old[i, 0]) -
                       max(h0[i-1, 0], h0[i, 0]))

        hflow_y2 = abs(max(hn_old[i, 0], hn_old[i+1, 0]) -
                       max(h0[i, 0], h0[i+1, 0]))

        # In X direction
        if hflow_x2 == 0.:
            dt_x2 = DT0
        elif abs(hn_old[i, 1] - hn_old[i, 0]) <= HLIN and hflow_x2 != 0.:
            dt_x2 = DX ** 2 / 4. * N / hflow_x2 ** (5./3.) *\
                (HLIN/DX) ** (1./2.)
        elif abs(hn_old[i, 1] - hn_old[i, 0]) > HLIN and hflow_x2 != 0.:
            dt_x2 = DX ** 2 / 4. * 2. * N / hflow_x2 ** (5./3.)\
                * (abs(hn_old[i, 1] - hn_old[i, 0]) / DX) ** (1./2.)
        # In Y direction
        if hflow_y1 == 0.:
            dt_y1 = DT0
        elif abs(hn_old[i-1, 0] - hn_old[i, 0]) <= HLIN and hflow_y1 != 0.:
            dt_y1 = DX ** 2 / 4. * N / hflow_y1 ** (5./3.) *\
                (HLIN/DY) ** (1./2.)
        elif abs(hn_old[i-1, 0] - hn_old[i, 0]) > HLIN and hflow_y1 != 0.:
            dt_y1 = DX ** 2 / 4. * 2. * N / hflow_y1 ** (5./3.)\
                * (abs(hn_old[i-1, 0] - hn_old[i, 0]) / DY) ** (1./2.)

        if hflow_y2 == 0.:
            dt_y2 = DT0
        elif abs(hn_old[i, 0] - hn_old[i+1, 0]) <= HLIN and hflow_y2 != 0.:
            dt_y2 = DX ** 2 / 4. * N / hflow_y2 ** (5./3.) *\
                (HLIN/DY) ** (1./2.)
        elif abs(hn_old[i, 0] - hn_old[i+1, 0]) > HLIN and hflow_y2 != 0.:
            dt_y2 = DX ** 2 / 4. * 2. * N / hflow_y2 ** (5./3.)\
                * (abs(hn_old[i, 0] - hn_old[i+1, 0]) / DY) ** (1./2.)

        #  # Boundary cells connected with Rivers
        #  if h0[i, 1] == 0.:
        #      dt_x2 = DT0
        adt.append(min(dt_y1, dt_y2, dt_x2))
        #  if min(adt) < DT0:
        #      breakpoint()

        #######################################################################
        # Four lines: 4. East [1:-1]
        hflow_x1 = abs(max(hn_old[i, -1], hn_old[i, -2]) -
                       max(h0[i, -1], h0[i, -2]))

        hflow_y1 = abs(max(hn_old[i-1, -1], hn_old[i, -1]) -
                       max(h0[i-1, -1], h0[i, -1]))

        hflow_y2 = abs(max(hn_old[i, -1], hn_old[i+1, -1]) -
                       max(h0[i, -1], h0[i+1, -1]))

        # In X direction
        if hflow_x1 == 0.:
            dt_x1 = DT0
        elif abs(hn_old[i, -1] - hn_old[i, -2]) <= HLIN and hflow_x1 != 0.:
            dt_x1 = DX ** 2 / 4. * N / hflow_x1 ** (5./3.) *\
                (HLIN/DX) ** (1./2.)
        elif abs(hn_old[i, -1] - hn_old[i, -2]) > HLIN and hflow_x1 != 0.:
            dt_x1 = DX ** 2 / 4. * 2. * N / hflow_x1 ** (5./3.)\
                * (abs(hn_old[i, -1] - hn_old[i, -2]) / DX) ** (1./2.)
        # In Y direction
        if hflow_y1 == 0.:
            dt_y1 = DT0
        elif abs(hn_old[i-1, -1] - hn_old[i, -1]) <= HLIN and hflow_y1 != 0.:
            dt_y1 = DX ** 2 / 4. * N / hflow_y1 ** (5./3.) *\
                (HLIN/DY) ** (1./2.)
        elif abs(hn_old[i-1, -1] - hn_old[i, -1]) > HLIN and hflow_y1 != 0.:
            dt_y1 = DX ** 2 / 4. * 2. * N / hflow_y1 ** (5./3.)\
                * (abs(hn_old[i-1, -1] - hn_old[i, -1]) / DY) ** (1./2.)

        if hflow_y2 == 0.:
            dt_y2 = DT0
        elif abs(hn_old[i, -1] - hn_old[i+1, -1]) <= HLIN and hflow_y2 != 0.:
            dt_y2 = DX ** 2 / 4. * N / hflow_y2 ** (5./3.) *\
                (HLIN/DY) ** (1./2.)
        elif abs(hn_old[i, -1] - hn_old[i+1, -1]) > HLIN and hflow_y2 != 0.:
            dt_y2 = DX ** 2 / 4. * 2. * N / hflow_y2 ** (5./3.)\
                * (abs(hn_old[i, -1] - hn_old[i+1, -1]) / DY) ** (1./2.)

        #  # Boundary cells connected with Rivers
        #  if h0[i, -2] == 0.:
        #      dt_x1 = DT0
        adt.append(min(dt_y1, dt_y2, dt_x1))
        #  if min(adt) < DT0:
        #      breakpoint()

    ###########################################################################
    # Core body [1:-1, 1:-1]
    for i, j in product(range(1, NX-1), range(1, NY-1)):
        hflow_x1 = abs(max(hn_old[i, j], hn_old[i-1, j]) -
                       max(h0[i, j], h0[i-1, j]))
        hflow_x2 = abs(max(hn_old[i, j], hn_old[i+1, j]) -
                       max(h0[i, j], h0[i+1, j]))
        hflow_y1 = abs(max(hn_old[i, j], hn_old[i, j-1]) -
                       max(h0[i, j], h0[i, j-1]))
        hflow_y2 = abs(max(hn_old[i, j], hn_old[i, j+1]) -
                       max(h0[i, j], h0[i, j+1]))

        # In X direction
        if hflow_x1 == 0.:
            dt_x1 = DT0
        elif abs(hn_old[i, j] - hn_old[i-1, j]) <= HLIN and hflow_x1 != 0.:
            dt_x1 = DX ** 2 / 4. * N / hflow_x1 ** (5./3.) *\
                (HLIN/DX) ** (1./2.)
        elif abs(hn_old[i, j] - hn_old[i-1, j]) > HLIN and hflow_x1 != 0.:
            dt_x1 = DX ** 2 / 4. * 2. * N / hflow_x1 ** (5./3.)\
                * (abs(hn_old[i, j] - hn_old[i-1, j]) / DX) ** (1./2.)

        if hflow_x2 == 0.:
            dt_x2 = DT0
        elif abs(hn_old[i+1, j] - hn_old[i, j]) <= HLIN and hflow_x2 != 0.:
            dt_x2 = DX ** 2 / 4. * N / hflow_x2 ** (5./3.) *\
                (HLIN/DX) ** (1./2.)
        elif abs(hn_old[i+1, j] - hn_old[i, j]) > HLIN and hflow_x2 != 0.:
            dt_x2 = DX ** 2 / 4. * 2. * N / hflow_x2 ** (5./3.)\
                * (abs(hn_old[i, j] - hn_old[i+1, j]) / DX) ** (1./2.)

        # In Y direction
        if hflow_y1 == 0.:
            dt_y1 = DT0
        elif abs(hn_old[i, j-1] - hn_old[i, j]) <= HLIN and hflow_y1 != 0.:
            dt_y1 = DX ** 2 / 4. * N / hflow_y1 ** (5./3.) *\
                (HLIN/DY) ** (1./2.)
        elif abs(hn_old[i, j-1] - hn_old[i, j]) > HLIN and hflow_y1 != 0.:
            dt_y1 = DX ** 2 / 4. * 2. * N / hflow_y1 ** (5./3.)\
                * (abs(hn_old[i, j] - hn_old[i, j-1]) / DY) ** (1./2.)

        if hflow_y2 == 0.:
            dt_y2 = DT0
        elif abs(hn_old[i, j] - hn_old[i, j+1]) <= HLIN and hflow_y2 != 0.:
            dt_y2 = DX ** 2 / 4. * N / hflow_y2 ** (5./3.) *\
                (HLIN/DY) ** (1./2.)
        elif abs(hn_old[i, j] - hn_old[i, j+1]) > HLIN and hflow_y2 != 0.:
            dt_y2 = DX ** 2 / 4. * 2. * N / hflow_y2 ** (5./3.)\
                * (abs(hn_old[i, j] - hn_old[i, j+1]) / DY) ** (1./2.)

        # If the cells are river
        if h0[i, j] == 0.:
            dt_x1 = dt_x2 = dt_y1 = dt_y2 = DT0

        adt.append(min(dt_x1, dt_x2, dt_y1, dt_y2))
        #  if min(adt) < DT0:
        #      #  print(i, j, h0[i, j])
        #      breakpoint()

    # Fix the adaptive time step in space at each time step
    DT = min(adt)

    # Update the new value to the old value for the next loop
    hn_old = copy.deepcopy(hn_new)
    qn_old = copy.deepcopy(qn_new)
    TN += DT
    adt_dt.append(DT)

    # Record the evolution of flood volume
    tn.append(TN)
    #  volume.append(Positive_Sum(hn_new-h0))
    volume.append(np.sum(hn_new-h0))

# Save the final results of the water depth
H_SAVE = pd.DataFrame(hn_new-h0)
H_SAVE.to_csv('../Test/Overland_%d.csv' % (NT))

# One-column figure
# fig = plt.figure(constrained_layout=True, figsize=(4.5,3.0), dpi=100)
# Two-column figure
#  fig = plt.figure(constrained_layout=True, figsize=(7.2, 4.8), dpi=100)
#  spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)
plt.figure(figsize=(7.2, 4.8), dpi=100)

# Initial condition: h(x, 0)
ax1 = plt.subplot(2, 2, 1)
ax1 = sns.heatmap(h0, cmap='RdYlGn_r', linewidth=0.5, annot=False,
                  vmax=12, vmin=0)
ax1.tick_params(axis='both', which='both', length=0)
ax1.set(xlabel='DX = %s m' % (DX), ylabel='DY = %s m' % (DY))
ax1.title.set_text('Initial condition')

# Boundary condition: h(0, t)
ax2 = plt.subplot(2, 2, 2)
ax2.plot(t, h_flood, 'k-', label='Upstream boundary condition')
ax2.set_ylim(0, 15.0)
ax2.set_xlim(0, 259200)
ax2.set(xlabel='t (s)', ylabel='h (m)')
ax2.title.set_text('Upstream boundary condition')
ax2.legend(loc=0, numpoints=1, framealpha=1.0, edgecolor='black')

# Comparison between numerical and analytic solution
ax3 = plt.subplot(2, 2, 3)
ax3 = sns.heatmap(hn_new-h0, cmap='RdYlGn_r', linewidth=0.5, annot=False,
                  vmax=12, vmin=0)
ax3.tick_params(axis='both', which='both', length=0)
ax3.set(xlabel='DX = %s m' % (DX), ylabel='DY = %s m' % (DY))
ax3.title.set_text('Simulated water depth')

# Flood volume evolution
ax4 = plt.subplot(2, 2, 4)
ax4.plot(tn, volume, 'k-', label='Evolution of flood volume')
#  ax4.set_ylim(0, 20.0)
ax4.set_xlim(0, 7200)
ax4.set(xlabel='t (s)', ylabel='V (m$^3$)')
ax4.title.set_text('Evolution of flood volume')
ax4.legend(loc=0, numpoints=1, framealpha=1.0, edgecolor='black')

'''
# Evolution of time step during the solution
ax4 = plt.subplot(2, 2, 4)
ax4.plot(tn, adt_dt, 'k-', label='Adaptive time step')
ax4.set_yscale('log')
ax4.set_xlim(0, 7200)
ax4.set_ylim(0.001, 10)
ax4.set(xlabel='t (s)', ylabel='dt (s)')
ax4.title.set_text('Evolution of adaptive time step')
ax4.legend(loc=0, numpoints=1, framealpha=1.0, edgecolor='black')
'''

plt.savefig('../Test/Overland_adt.png')
plt.show()
