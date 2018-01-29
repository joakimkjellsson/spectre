##
##
##

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

from plot_functions_2 import *

starttime = datetime(2009,1,5)
endtime   = datetime(2009,12,31)
step      = timedelta(days=5)
ddir      = '/Users/jkjellsson/Downloads/'
data_list = read_stored_data(['INALT10.L46-KJH0017-NEST1'],['agulhas-retro'],'test',\
                             starttime,endtime,step,\
                             ddir)

plot_stored_data(data_list)

plt.show()
sys.exit()