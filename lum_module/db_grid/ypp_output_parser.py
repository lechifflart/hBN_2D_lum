import numpy as np

def ypp_output_parser(file,exc_range):
    # parse ypp output file
    # returns exciton interpolated energies on finer Grid
    #
    with open(file) as f:
        data = np.genfromtxt(fname=file,comments='#')
    qtilde = data[:,-3:]
    print("[WARNING] : this function will parse the q-grid several times\n \t you should check they are all the same")
    exc_nrg = data[:,exc_range[0]:exc_range[1]+1]

    return(qtilde,exc_nrg)

#qtilde, exc_nrgs = ypp_output_parser('/home/lechifflart/hBN_2D/test_code/o-BSE.excitons-1-4_interpolated_IBZ',[1,4])
#print(exc_nrgs[-1])
