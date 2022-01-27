import numpy as np
from ypp_matdyn_parser import ypp_matdyn_parser
from matdyn_freq_parser import matdyn_freq_parser
#from interpolate import *
from yambopy import *
import sys
import matplotlib.pyplot as plt
np.set_printoptions(threshold=sys.maxsize)

# define qtilde_grid (or get it from matdyn ??)

# parse matdyn output
q_matdyn, freq_matdyn = matdyn_freq_parser()

# parse ypp output for the excitons energies interpolated on the matdyn fine grid
list_of_points_to_expand, list_of_exc_nrgs = ypp_matdyn_parser('/home/lechifflart/hBN_2D/test_code/test_data/o.excitons-1-9_grid9x9',[1,9])

# visualize which points are selected
car_uniques = red_car(uniques,lat.lat)
car_ibz = red_car(list_of_points_to_expand,lat.lat)
plt.scatter(list_of_points_to_expand[:,0],list_of_points_to_expand[:,1])
plt.scatter(uniques[:,0],uniques[:,1],marker='+',s = 80)
plt.show()

#
# find all qtilde points around each Q point (from BSE)
    # read data
with open('test_data/o-BSE.excitons_interpolated_IBZ','r') as f:
    data_BSE = np.genfromtxt(fname=f,comments="#")
exc_a = data_BSE[:,1:3]
qBSE = data_BSE[:,-3:]
    # get the FBZ original Q grid
for QQ in qBSE :
    print(QQ)

# analytic fit around Q = 0

# structure would be

"""
if Q==0:
    if beta<2:
        analytic fit and matdyn
    else ypp_interp and matdyn

else :
    ypp_interp and matdyn
"""
