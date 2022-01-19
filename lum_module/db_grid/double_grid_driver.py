import numpy as np
from ypp_output_parser import ypp_output_parser
#from matdyn_freq_parser import *
#from interpolate import *
from yambopy import *
import sys

# define qtilde_grid (or get it from matdyn ??)

# parse matdyn output ( do we launch matdyn jobs from here ?)

# parse ypp output with interpolated grid 108x108x1, for all 20 excitons
        # check if each grid is the same

# expand IBZ exciton energies to FBZ (expand_kpoints from yambopy)
        # need the list of symmetries
list_of_points_to_expand, list_of_exc_nrgs = ypp_output_parser('/home/lechifflart/hBN_2D/test_code/o-BSE.exc-1-4_interpolated_IBZ',[1,4])
lat = YamboLatticeDB.from_db_file(filename='SAVE/ns.db1')

list_of_symmetries = lat.sym_car
print(expand_kpts(list_of_points_to_expand,lat.sym_car))
#print(lat.expand_kpoints)
sys.exit(1)
# find all qtilde points around each Q point (from BSE)

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
