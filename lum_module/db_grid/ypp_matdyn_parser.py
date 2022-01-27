import numpy as np

def ypp_output_parser(file,exc_range):
    # parse ypp output file
    # returns exciton interpolated energies on matdyn output grid
    #
    with open(file) as f:
        data = np.genfromtxt(fname=file,comments='#')
    qtilde = data[:3,:]
    exc_nrg = data[3:,:]

    return(qtilde,exc_nrg)

#qtilde, exc_nrgs = ypp_output_parser('/home/lechifflart/hBN_2D/test_code/o-BSE.excitons-1-4_interpolated_IBZ',[1,4])
#print(exc_nrgs[-1])
