
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np
#from lum_module.parameters import eV2ha
import sys

def interpolate(qqq, eig1, eig2,quad,lin) :
    params1, par_cov1 = optimize.curve_fit(quad, qqq[:8], eig1[:8])
    # print(params1)
    # print(par_cov1)
    # plt.scatter(qqq,eig1)
    # plt.plot(qqq,quad(qqq,params1[0]))
    #plt.show()
    params2, par_cov2 = optimize.curve_fit(lin, qqq[:5], eig2[:5])
    # print(params2)
    # print(par_cov2)
    # plt.scatter(qqq,eig2)
    # plt.plot(qqq,lin(qqq,params2[0]))
    # plt.show()
    return(params1,params2)

#interpolate(qpoints,exc1,exc2)
