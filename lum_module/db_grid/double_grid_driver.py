import numpy as np
from lum_module.db_grid.ypp_matdyn_parser import ypp_matdyn_parser
from lum_module.db_grid.matdyn_freq_parser import matdyn_freq_parser
from lum_module.db_grid.interpolate import interpolate
from yambopy import lattice
import sys
import matplotlib.pyplot as plt
np.set_printoptions(threshold=sys.maxsize)

def double_grid_driver(QQ,matdyn_freq_file,ypp_matdyn_file,exc_QQ_alpha,exc_beta,ph_freq):
    """
    This first version here is to integrate the energy denominators around Q=Gamma with a finer grid (dbg);
    finite Q integration will use db grid too, already started implementing it

    returns : one array for each QQ, containing energy denominators of dim nb_beta, nb_phonons
    """
    Gamma = np.array([0.,0.,0.])
    # parse matdyn output
    q_matdyn, freq_matdyn = matdyn_freq_parser(matdyn_freq_file)

    # parse ypp output for the excitons energies interpolated on the matdyn fine grid
    qtilde_dbg, exc_beta_dbg = ypp_matdyn_parser(ypp_matdyn_file)

    # expand BSE energies from IBZ to FBZ

    #
    # find all qtilde points around each Q point (from BSE)
    #
    #############################################################
    #               FOR NOW IT IS ONLY GAMMA                    #
    #############################################################
    #
    # read data
    full_path = "/home/lechifflart/hBN_2D/test_code/BSE_QDISP/"
    with open(full_path+"exciton_quad_lin.dat") as f:
        data_fit = np.genfromtxt(fname=f,comments='#')
    eV2ha = 1./27.211396132
    QBSE_fit = data_fit[:,0]
    exc1_fitQ = data_fit[:,1] * eV2ha
    exc2_fitQ = data_fit[:,2] * eV2ha

    # analytic fit around Q = 0
    # call interpolate.py, get fit parameters, define functions quad and lin

    if (QQ == Gamma).all():
        exc1_0 = exc1_fitQ[0]
        exc2_0 = exc2_fitQ[0]
        def quad(q,b):
            return(exc1_0 + b * q**2)

        def lin(q,b):
            return(exc2_0 + b*abs(q))

        analytic_fit_params = interpolate(QBSE_fit,exc1_fitQ,exc2_fitQ,quad,lin)

    # check if all qtilde are the same
    if not np.array(lattice.point_matching(qtilde_dbg,q_matdyn,double_check=False)==lattice.point_matching(q_matdyn,qtilde_dbg,double_check=False)).any() :
        print("some points are not matching between matdyn and ypp outputs")
    # here all the points are the same, in the same order

    # double grid integration
    energy_denominators = np.empty([len(ph_freq),len(exc_beta)],dtype=float)
    for mu in range(len(ph_freq)):
        for beta in range(len(exc_beta)):
            energy_denominators[mu,beta] =  np.sum((exc_beta_dbg[:,beta] - exc_QQ_alpha.real + freq_matdyn[:,mu])**2) / len(qtilde_dbg)

        if (QQ==Gamma).all():
            energy_denominators[mu,0] = np.sum([(quad(qt,analytic_fit_params[0]) - exc_QQ_alpha.real + freq_matdyn[iq,mu] )**2 for iq,qt in enumerate(qtilde_dbg)])
            energy_denominators[mu,1] = np.sum([(lin(qt,analytic_fit_params[1]) - exc_QQ_alpha.real + freq_matdyn[iq,mu] )**2 for iq,qt in enumerate(qtilde_dbg)])

    return(energy_denominators)

#double_grid_driver([0.,0.,0.],'test_datamatdyn_dbg.freq/','test_data/o.excitons-1-9_grid9x9',0.2,exc_beta,ph_freq)
