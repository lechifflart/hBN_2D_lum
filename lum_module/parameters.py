broadening = 0.01 /27.211396132 # ha
Nq = 1 #1296         # nb of q points in the full BZ

#################################################################################################
index_fortran_q = 1     # index of the q point in Fortran indexing, ie starts at 1
idq = str(index_fortran_q)

# constants and conversion factors
c_light  = 137.03599911   # speed of light, a.u.
#c_light     = 299792458e10   # angstrom / s
hbar        = 1. #6.5821e-16     # eV * s

k_boltzmann = 3.166811563e6 # ha / K
#k_boltzmann = 6.241506363142715e18# 8.617333e-5    # eV / K

ha2ev = 27.211396132
ha2cmm1 = ha2ev * 8065.5
ha2kel = 3.1577464e5  #11604.5077*ha2ev

eV2ha = 1./27.211396132
