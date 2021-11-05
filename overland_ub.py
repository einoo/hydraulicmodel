'''
Upstream boundary condition of the overland flow
Using the mreasured river water level at Hayano station directly
'''

import pandas as pd

# Total simulation period: 3 days = 259200 seconds
# Total number of measured river water levels with 10 min interval
NT = 433
T = []
for i in range(NT):
    T.append(i * 600)

ub_overland = pd.read_csv('../Test/ub_river.csv')
FL = ub_overland['Hayano'][0:NT].astype(float)
FL = [i/100. for i in FL]
