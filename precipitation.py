"""
Read the precipitation data
"""

import pandas as pd
import numpy as np
import logging

def precipitation(src, n):
    """
    Read the precipitation data

    The precipitation can be uniformly or non uniformly distributed in the study area.
    Here is an example for the uniformly distribution rainfall

    Args:
        src (file): The source file of precipitation obtained from the meteorological stations
        n (int): The total time step for the simulation of floods

    Returns:
        list (float): The amount of precipitation (mm) in each time step (10 minutes here)
    """

    logging.info('The precipitation file is from %s.' % (src))
    data = pd.read_csv(src)
    # import pdb; pdb.set_trace()

    prec = np.zeros(n)
    for i in range(1, len(data)):
        s = (i-1)*6 + 1
        e = i*6 + 1
        prec[s:e] = data.iloc[i, 1]/6.0
    logging.info('The intensive precipitation is uniformly distributed in the time interval (10 min in this case)')
    # import pdb; pdb.set_trace()

    return prec
