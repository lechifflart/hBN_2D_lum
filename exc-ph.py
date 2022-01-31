#!/bin/python3.8
#
# for now everything is calculated at Gamma ; just need one extra index in every loop to include finite q's
#
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
#
from lum_module import functions as funcs
from lum_module.parameters import *
from lum_module.db_grid import double_grid_driver as dbg

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--plot_abs', help='Plots the phonon-assisted absorption spectrum', action='store_true')
    parser.add_argument('--plot_lum', help='Plots the photoluminescence using the RS relation', action='store_true')
    parser.add_argument('--write_static', help='Writes the static dielectric function in static.dat file', action='store_true')
    parser.add_argument('--write_dynamic', help='Writes the dielectric function with a renormalization (NO SATELLITE)', action='store_true')

    args = parser.parse_args()

    ########################################################################
    #    Read netCDF databases written by yambo
    #
    with nc.Dataset('./BSE/ndb.cutoff') as data_cutoff:
        re_bare_qpg,im_bare_qpg = data_cutoff.variables['CUT_BARE_QPG'][:].T
        bare_qpg=re_bare_qpg + im_bare_qpg * 1j

    with nc.Dataset('ndb.elph_gkkp_expanded_fragment_1') as data_ph:
        ph_freq         =    data_ph.variables['PH_FREQS1'][:].T
        ph_freq = np.sqrt(ph_freq)

    # skip the first 3 modes (accoustic modes should have 0 frequency at Gamma)
    phfreq_indx = [3,6]
    for mu in range(phfreq_indx[0],phfreq_indx[1]):
        print(mu,ph_freq[mu])

    with nc.Dataset('./BSE/ndb.BS_diago_Q1_Lbar') as data_exc:
         #residuals
         rel,iml=data_exc.variables['BS_left_Residuals'][:].T
         rer,imr=data_exc.variables['BS_right_Residuals'][:].T
         l_residual = rel+iml*1j
         r_residual = rer+imr*1j

         #energies
         eig =  data_exc.variables['BS_Energies'][:]
         exc_eigenvalues = eig[:,0]+eig[:,1]*1j

    with nc.Dataset('./BSE_EXCPH/ndb.excph_Gkkp_fragment_1') as database:
            ph_modes         = database.dimensions['PH_modes'].size
            N_exc_alpha      = database.dimensions['N_exc_states'].size
            N_exc_beta       = database.dimensions['N_exc_sum'].size
            Gkkp_sq          = database.variables['EXCITON_PH_GKKP_SQUARED_Q'+idq][:].T
            #Gkkp             = database.variables['EXCITON_PH_GKKP_Q1'][:].T

    #Gkkp_sq *= ha2ev**2
    if (Gkkp_sq < 0.).any() :
        sys.exit('Gkkp negative')

    eig_a = exc_eigenvalues[:N_exc_alpha]
    eig_b = exc_eigenvalues[:N_exc_beta]

    energy_range = [5.00, 6.00]
    energy_n_steps = 2000
    d_en = (energy_range[1]-energy_range[0])/energy_n_steps
    omega_eV = [energy_range[0]+i*d_en for i in range(energy_n_steps)]
    omega = np.asarray(omega_eV)*eV2ha
    oscillator_strengths = l_residual*r_residual

    prefactor = 2.5512034E-05/ bare_qpg[0,0]**2 # Co_factor / bare_qpg(1,1) **2

    #################################################################################
    #  write static response function to be compared with the yambo output o-*.eps*
    #
    if args.write_static :
        with open('static.dat','w') as file:
            file.write("#E/ev \t EPS-Re \t EPS-Im \n")
            for i in range(energy_n_steps):
                eps_out=(1.- funcs.chi_res_static(omega[i],oscillator_strengths,exc_eigenvalues) )* prefactor
                file.write(f'{omega_eV[i]} \t {eps_out.real} \t  {eps_out.imag}\n')

    ########################################################################
    #  write response function with first-order dynamical correction
    #
    Z_qa = funcs.renorm_factor(ph_freq,phfreq_indx,eig_a,eig_b,1,Gkkp_sq)
    if (Z_qa.real > 1.).any() :
        #sys.exit('\n ERROR : Renorm factor > 1 \n')
        print(' \n')

    if args.plot_abs :
        eps_to_plot = np.zeros((energy_n_steps),dtype=np.cfloat)
        eps_stat = np.zeros((energy_n_steps),dtype=np.cfloat)
        for i in range(energy_n_steps):
            eps_to_plot[i]=1.- (funcs.chi_one_dyn(omega[i],
                                        oscillator_strengths, eig_a,eig_b,
                                        ph_freq,phfreq_indx,
                                        Gkkp_sq,
                                        Z_qa,1)[2]
                         ) * prefactor
            eps_stat[i]=1.-(funcs.chi_res_static(omega[i],oscillator_strengths,exc_eigenvalues) )*prefactor

        plt.xlim([5.25,5.7])
        plt.ylim([0.0,1.2e-5])
        ZO = exc_eigenvalues[0].real+ph_freq[3].real
        LO = exc_eigenvalues[1].real+ph_freq[5].real
        plt.axvline(x=ZO,ymin=0,ymax=1,color='black',lw=0.5,linestyle='--')
        plt.axvline(x=LO,ymin=0,ymax=1,color='black',lw=0.5,linestyle='--')

        plt.plot(omega_eV,eps_to_plot.imag,label='first order dynamic correction',color='orange')
        plt.plot(omega_eV,eps_stat.imag,color='black',linestyle='dotted', label='static')
        plt.xlabel('Frequency [eV]', fontsize = 16)
        plt.ylabel(r"$\epsilon_{2}(\omega)$", fontsize=16)

        plt.annotate("ZO", xy=(ZO,0.4e-5), xytext=(ZO-0.05, 0.5e-5),arrowprops=dict(arrowstyle="-"))
        plt.annotate("LO", xy=(LO,0.4e-5), xytext=(LO+0.05, 0.5e-5),arrowprops=dict(arrowstyle="-"))
        plt.legend()
        plt.show()

    if args.write_dynamic :
        with open('dynamic.dat','w') as file:
           file.write("#E/ev \t EPS-Re \t EPS-Im \n")
           for i in range(energy_n_steps):
               [eps_out_stat,eps_out_dyn,eps_out_full] = funcs.dielectric_function(omega[i],oscillator_strengths,eig_a,eig_b, ph_freq,phfreq_indx, Gkkp_sq, Z_qa,1,prefactor)
               file.write(f'{omega_eV[i]} \t {eps_out_stat.real} \t  {eps_out_stat.imag}\t {eps_out_dyn.real} \t  {eps_out_dyn.imag}\t {eps_out_full.real} \t  {eps_out_full.imag} \n')


    T_exc = funcs.excitonic_temperature(5)   # input is in Kelvin
    if args.plot_lum :
        R0_sp = np.empty(energy_n_steps,dtype=float)
        R_sat_sp = np.empty(energy_n_steps,dtype=float)
        R_full_sp = np.empty(energy_n_steps,dtype=float)

        boltzmann_occ = 1. #funcs.boltzmann(eig_a.real,T_exc, exc_eigenvalues[0].real)
        for w in range(energy_n_steps):
            eps_0, eps_sat, eps_full = funcs.dielectric_function(omega[w], oscillator_strengths,eig_a,eig_b, ph_freq,phfreq_indx, Gkkp_sq, Z_qa,1,prefactor)
            refrac = funcs.refrac_index(eps_0)

            R0_sp[w] =  funcs.spontaneous_rate(omega[w],eps_0,refrac,boltzmann_occ)

            for iq in range(Nq):
                for mu in range(phfreq_indx[0],phfreq_indx[1]):
                    sat = funcs.dielectric_function(omega[w] + 2*ph_freq[mu], oscillator_strengths,eig_a,eig_b, ph_freq,phfreq_indx, Gkkp_sq, Z_qa,1,prefactor)[1].imag
                    R_sat_sp[w] += funcs.satellite_rate(omega[w],sat,refrac,boltzmann_occ,ph_freq[mu],1)

            R_full_sp[w] = R0_sp[w] + 1./Nq * R_sat_sp[w]

        plt.xlabel('Frequency [eV]', fontsize = 16)
        plt.ylabel(r"R$^{sp}$", fontsize=16)
        plt.plot(omega_eV,R_full_sp, label='full',color='red')
        plt.plot(omega_eV,R_sat_sp,label='satellite only',color='blue')
        plt.plot(omega_eV,R0_sp,label='static',color='orange')
        plt.legend()
        plt.show()

    # print(Z_qa)
    energy_denominators = dbg.double_grid_driver([0.,0.,0.],'test_data/matdyn_dbg.freq','test_data/o.excitons-1-9_grid9x9',eig_a[0],eig_b[:9],ph_freq)
    print(energy_denominators[:2])

if __name__ == '__main__' :
    main()
