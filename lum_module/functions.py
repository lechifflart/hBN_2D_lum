from lum_module.parameters import *
import numpy as np

def chi_res_static(omega,strength,eig):
    chi=0.+0.*1j
    for alpha in range(len(eig)):
        chi+= strength[alpha] / (omega - eig[alpha] + 1j*broadening)
    return(chi)

def renorm_factor(ph_freq,phfreq_indx,eig_a,eig_b,Nq,Gkkp_sq):
    # arguments : phonon frequencies, energies of incoming excitons, " of scattered excitons, nb of q points, exc-ph matrix elmts squared
    # to be at Gamma, iq1 = 0
    nb_a = len(eig_a)
    nb_b = len(eig_b)
    R_qa = np.zeros((Nq,nb_a),dtype=np.cfloat)
    #
    #print('eig_a',eig_a.real)
    #print('eig_b',eig_b.real)
    for iq1 in range(Nq):
        for alpha in range(nb_a-1):
            dummy_sum=0.
            for mu in range(phfreq_indx[0],phfreq_indx[1]):
                dummy_sum += np.sum(Gkkp_sq[mu,:,alpha] / (eig_b[:]-eig_a[alpha]+ph_freq[mu])**2)
                #print('ph_freq',ph_freq[mu]**2)
                #print('Gkkp_sq',Gkkp_sq[mu,:2,alpha])
                #print('dummy_sum',dummy_sum)
            R_qa[iq1,alpha] = dummy_sum
    return(1./Nq * R_qa)

def chi_one_dyn(omega,strength,eig_a,eig_b,ph_freq,phfreq_indx,Gkkp_sq,R_qa,Nq): # first-order correction in the dynamical part of the interaction
    # arguments : light freq, coupling strength, energies of incoming excitons, " of scattered excitons,exc-ph matrix elmts squared, phonon freqs, renorm_factor, nb of q points
    # coupling strength = |T_{0\alpha}|^2
    chi_stat=0.+0.*1j
    #chi_dyn=0.+0.*1j
    chi=0.+0.*1j
    nb_a = len(eig_a)
    nb_b = len(eig_b)

    for alpha in range(nb_a):
        for iq1 in range(Nq):   # at Gamma, Nq = 1
            chi_stat += abs(strength[alpha])*(1. - R_qa[iq1,alpha]) / (omega - eig_a[alpha] + 1j*broadening)
            chi_dyn = 0. + 0.*1j
            for mu in range(phfreq_indx[0],phfreq_indx[1]):
                chi_dyn += 1./Nq *abs(strength[alpha]) * np.sum(Gkkp_sq[mu,:,alpha] / (eig_b - eig_a[alpha] + ph_freq[mu])**2 / (omega - eig_b - ph_freq[mu] + 1j*broadening) )
    chi += chi_stat + chi_dyn
    return(chi_stat,chi_dyn,chi)

############################################################################
# Luminescence
def refrac_index(eps):
    return(np.sqrt(1./2.*np.sqrt(eps.real**2+eps.imag**2)+eps.real))

def spontaneous_rate(omega,eps,refrac,boltzmann):
    return(boltzmann * refrac * omega**3 * eps.imag / np.pi**2 / c_light**3 / hbar**4)

def satellite_rate(omega,sat,refrac,boltzmann, ph_freq, Nq):
    return(boltzmann * refrac * omega * (omega - 2*ph_freq)**2 * sat/ np.pi**2 / c_light**3 / hbar**4)

def boltzmann(omega, T, ref_E):
    return(np.exp(-(omega - ref_E)/(k_boltzmann * T)))

def excitonic_temperature(T):
    # Fit of experimental photoluminescence in bulk hBN
    # argument is lattice temperature ie temp of the sample in the experiment in Kelvin
    # output is excitonic temperature in Kelvin
    return((6.67858 + 1.78924 * T))

def dielectric_function(omega, strengths, eig_a, eig_b, ph_freq,phfreq_indx, Gkkp_sq, Z_qa, Nq, prefactor):
    # returns eps_static, eps_dyn, eps_full
    #    == eps_0, eps_sat, eps_full
    [eps_static, eps_dyn, eps_full] = np.asarray(chi_one_dyn(omega,strengths,
                             eig_a,eig_b,ph_freq,phfreq_indx,Gkkp_sq,Z_qa,1))
    return(1.-np.asarray([eps_static, eps_dyn, eps_full]) * prefactor)
