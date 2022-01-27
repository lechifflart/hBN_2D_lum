import numpy as np
from ypp_output_parser import ypp_output_parser
#from matdyn_freq_parser import *
#from interpolate import *
from yambopy import *
import sys
import matplotlib.pyplot as plt
np.set_printoptions(threshold=sys.maxsize)

# define qtilde_grid (or get it from matdyn ??)

# parse matdyn output ( do we launch matdyn jobs from here ?)

# parse ypp output with interpolated grid 108x108x1, for all 20 excitons
        # check if each grid is the same

# expand IBZ exciton energies to FBZ (expand_kpoints from yambopy)
        # need the list of symmetries
list_of_points_to_expand, list_of_exc_nrgs = ypp_output_parser('/home/lechifflart/hBN_2D/test_code/o.exc-1-4_interpolated_IBZ_testgrid',[1,9])
lat = YamboLatticeDB.from_db_file(filename='SAVE/ns.db1')
#rec_sym = lat.sym_rec
rec_red_sym = lat.sym_rec_red

expanded_qtilde = expand_kpts(list_of_points_to_expand,rec_red_sym)
expanded_qtilde_idx = np.array(expanded_qtilde)[:,0]
expanded_qtilde_coord = np.array(expanded_qtilde)[:,1]

# discard the duplicates in the expanded points
new_array = [tuple(row.round(decimals=4)) for row in expanded_qtilde_coord]
uniques = np.unique(new_array,axis=0)
print(len(uniques))
print("uniques",uniques)

# visualize which points are selected
car_uniques = red_car(uniques,lat.lat)
car_ibz = red_car(list_of_points_to_expand,lat.lat)
plt.scatter(list_of_points_to_expand[:,0],list_of_points_to_expand[:,1])
plt.scatter(uniques[:,0],uniques[:,1],marker='+',s = 80)
plt.show()


sys.exit(1)

# find all qtilde points around each Q point (from BSE)
    # get the FBZ original Q grid
#expanded_Q = expand_kpts(,lat.sym_car)

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
