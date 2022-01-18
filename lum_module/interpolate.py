
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np
#from lum_module.parameters import eV2ha
import sys

# get the BSE data
full_path = "/home/lechifflart/hBN_2D/test_code/BSE_QDISP/"
with open(full_path+"exciton_quad_lin.dat") as f:
    data = np.genfromtxt(fname=f,comments='#')
eV2ha = 1./27.211396132
qpoints = data[:,0]
exc1 = data[:,1] #* eV2ha
exc2 = data[:,2] #* eV2ha

# functions to interpolate the first two exciton's dispersions
def quad(q,b):
    return(exc1[0] + b * q**2)

def lin(q,b):
    return(exc2[0] + b*abs(q))

def interpolate(qqq, eig1, eig2) :
    params1, par_cov1 = optimize.curve_fit(quad, qqq, eig1)
    print(params1)
    print(par_cov1)
    plt.scatter(qqq,eig1)
    plt.plot(qqq,quad(qqq,params1[0]))
    #plt.show()
    params2, par_cov2 =  optimize.curve_fit(lin, qqq, eig2)
    print(params2)
    print(par_cov2)
    plt.scatter(qqq,eig2)
    plt.plot(qqq,lin(qqq,params2[0]))
    plt.show()

interpolate(qpoints,exc1,exc2)

#sys.exit('fin du code pour l\'instant')

# params taken from the fits
quad_params = [5.39608951, 6.70135474]
lin_slope = [4.62603503]

# 3x3 regular grid centered on 0
K_min = -0.027778
q_grid = [[i,j] for i in range(-1,2) for j in range(-1,2)]
q_grid = np.array(q_grid) * abs(K_min) / 2.

print(q_grid)
ph_freqs=np.array([6.174255e-08, 3.330298e-06, 3.330298e-06, 1.433492e-05,
    4.394219e-05, 4.394219e-05])
# # double grid integration
# for alpha in range(2):
#     for beta in range(2):
#         dbg_sum = 0.
#         for q in q_grid:
#             if alpha == 1 :
#                 dbg_sum += 1. / (ebeta - quad(q,5.39608951, 6.70135474) + ph_freq)
#
# def db_grid_nrg(alpha,beta,mu):
#     dbg_sum = 0.
#     for q in q_grid:
#         dbg_sum += 1. / (ebeta - quad(q,5.39608951, 6.70135474) + ph_freq)
#     return(dbg_sum / np.shape(q_grid)[0])
