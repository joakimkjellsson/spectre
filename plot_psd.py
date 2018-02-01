##
##
##

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

from settings import setup_analysis
from plot_functions_2 import *
from read_saved_data import *

setup = setup_analysis()

starttime = setup['starttime']
endtime   = setup['endtime']
step      = setup['outputStep']
ddir      = setup['outdir']
data_list = read_stored_data(setup['names'],setup['regions'],setup['prefix'],\
                             starttime,endtime,step,\
                             ddir)

plot_stored_data(data_list)

plt.show()
sys.exit()