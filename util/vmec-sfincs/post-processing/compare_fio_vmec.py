import os, sys, glob
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import math

from netCDF4 import Dataset


def compare_FSA_quantities(vmec_wout1, vmec_wout2, ntheta, nphi, vmec1_descriptor=None, vmec2_descriptor=None):

    '''
    This function compares 2 VMEC files by constructing and plotting the flux surface average of |B|, g, B_theta, B_phi, B^theta, B^phi as a function of toroidal flux
    '''
    
    if (vmec1_descriptor == None):
        vmec1_descriptor = os.path.abspath(vmec_wout1)
    if (vmec2_descriptor == None):
        vmec2_descriptor = os.path.abspath(vmec_wout2)

    # Read first vmec file
    data_1 = Dataset(vmec_wout1, mode='r')    
    nfp_1 = data_1.variables['nfp'][:]
    mnmax_nyq_1 = data_1.variables['mnmax_nyq'][:]
    xm_nyq_1 = data_1.variables['xm_nyq'][:]
    xn_nyq_1 = data_1.variables['xn_nyq'][:]
    nsurf_1 = data_1.variables['ns'][:]
    tflux_1 = np.abs(data_1.variables['phi'][:])
    bmnc_1 = data_1.variables['bmnc'][:,:]
    bsubumnc_1 = data_1.variables['bsubumnc'][:,:]
    bsubvmnc_1 = data_1.variables['bsubvmnc'][:,:]
    bsubsmns_1 = data_1.variables['bsubsmns'][:,:]
    bsupumnc_1 = data_1.variables['bsupumnc'][:,:]
    bsupvmnc_1 = data_1.variables['bsupvmnc'][:,:]
    gmnc_1 = data_1.variables['gmnc'][:,:]
    rmnc_1 = data_1.variables['rmnc'][:,:]
    zmns_1 = data_1.variables['zmns'][:,:]
    mnmax_1 = data_1.variables['mnmax'][:]
    xm_1 = data_1.variables['xm'][:]
    xn_1 = data_1.variables['xn'][:]


    # Read second vmec file
    data_2 = Dataset(vmec_wout2, mode='r')    
    nfp_2 = data_2.variables['nfp'][:]
    mnmax_nyq_2 = data_2.variables['mnmax_nyq'][:]
    xm_nyq_2 = data_2.variables['xm_nyq'][:]
    xn_nyq_2 = data_2.variables['xn_nyq'][:]
    nsurf_2 = data_2.variables['ns'][:]
    tflux_2 = np.abs(data_2.variables['phi'][:])
    bmnc_2 = data_2.variables['bmnc'][:,:]
    bsubumnc_2 = data_2.variables['bsubumnc'][:,:]
    bsubvmnc_2 = data_2.variables['bsubvmnc'][:,:]
    bsubsmns_2 = data_2.variables['bsubsmns'][:,:]
    bsupumnc_2 = data_2.variables['bsupumnc'][:,:]
    bsupvmnc_2 = data_2.variables['bsupvmnc'][:,:]
    gmnc_2 = data_2.variables['gmnc'][:,:]
    rmnc_2 = data_2.variables['rmnc'][:,:]
    zmns_2 = data_2.variables['zmns'][:,:]
    mnmax_2 = data_2.variables['mnmax'][:]
    xm_2 = data_2.variables['xm'][:]
    xn_2 = data_2.variables['xn'][:]

 
    theta_1 = np.linspace(0, 2*np.pi, ntheta)
    phi_1 = np.linspace(0, 2*np.pi, nphi)
    mod_B_1 = np.zeros((nsurf_1,ntheta,nphi))
    B_sub_theta_1 = np.zeros((nsurf_1,ntheta,nphi))
    B_sub_phi_1 = np.zeros((nsurf_1,ntheta,nphi))
    B_sub_psi_1 = np.zeros((nsurf_1,ntheta,nphi))
    B_sup_theta_1 = np.zeros((nsurf_1,ntheta,nphi))
    B_sup_phi_1 = np.zeros((nsurf_1,ntheta,nphi))
    Jac_1 = np.zeros((nsurf_1,ntheta,nphi))
    R_1 = np.zeros((nsurf_1,ntheta,nphi))
    Z_1 = np.zeros((nsurf_1,ntheta,nphi))
    fsa_mod_B_1 = np.zeros(nsurf_1)
    fsa_B_sub_theta_1 = np.zeros(nsurf_1)
    fsa_B_sub_phi_1 = np.zeros(nsurf_1)
    fsa_B_sub_psi_1 = np.zeros(nsurf_1)
    fsa_B_sup_theta_1 = np.zeros(nsurf_1)
    fsa_B_sup_phi_1 = np.zeros(nsurf_1)
    fsa_Jac_1 = np.zeros(nsurf_1)
    dV_ds_1 = np.zeros(nsurf_1)

    theta_2 = np.linspace(0,2*np.pi,ntheta)
    phi_2 = np.linspace(0, 2*np.pi, nphi)
    mod_B_2 = np.zeros((nsurf_2,ntheta,nphi))
    B_sub_theta_2 = np.zeros((nsurf_2,ntheta,nphi))
    B_sub_phi_2 = np.zeros((nsurf_2,ntheta,nphi))
    B_sub_psi_2 = np.zeros((nsurf_2,ntheta,nphi))
    B_sup_theta_2 = np.zeros((nsurf_2,ntheta,nphi))
    B_sup_phi_2 = np.zeros((nsurf_2,ntheta,nphi))
    Jac_2 = np.zeros((nsurf_2,ntheta,nphi))
    lambda_2 = np.zeros((nsurf_2,ntheta,nphi))
    R_2 = np.zeros((nsurf_2,ntheta,nphi))
    Z_2 = np.zeros((nsurf_2,ntheta,nphi))
    fsa_mod_B_2 = np.zeros(nsurf_2)
    fsa_B_sub_theta_2 = np.zeros(nsurf_2)
    fsa_B_sub_phi_2 = np.zeros(nsurf_2)
    fsa_B_sub_psi_2 = np.zeros(nsurf_2)
    fsa_B_sup_theta_2 = np.zeros(nsurf_2)
    fsa_B_sup_phi_2 = np.zeros(nsurf_2)
    fsa_Jac_2 = np.zeros(nsurf_2)
    dV_ds_2 = np.zeros(nsurf_2)

    theta_1_grid, phi_1_grid = np.meshgrid(theta_1, phi_1, indexing='ij')
    theta_2_grid, phi_2_grid = np.meshgrid(theta_2, phi_2, indexing='ij')
    
    for isurf in range(nsurf_1):
        for imn in range(mnmax_nyq_1):
            mod_B_1[isurf] += bmnc_1[isurf,imn]*np.cos(xm_nyq_1[imn]*theta_1_grid - xn_nyq_1[imn]*phi_1_grid)
            B_sub_theta_1[isurf] += bsubumnc_1[isurf,imn]*np.cos(xm_nyq_1[imn]*theta_1_grid - xn_nyq_1[imn]*phi_1_grid)
            B_sub_phi_1[isurf] += bsubvmnc_1[isurf,imn]*np.cos(xm_nyq_1[imn]*theta_1_grid - xn_nyq_1[imn]*phi_1_grid)
            B_sub_psi_1[isurf] += bsubsmns_1[isurf,imn]*np.sin(xm_nyq_1[imn]*theta_1_grid - xn_nyq_1[imn]*phi_1_grid)
            B_sup_theta_1[isurf] += bsupumnc_1[isurf,imn]*np.cos(xm_nyq_1[imn]*theta_1_grid - xn_nyq_1[imn]*phi_1_grid)
            B_sup_phi_1[isurf] += bsupvmnc_1[isurf,imn]*np.cos(xm_nyq_1[imn]*theta_1_grid - xn_nyq_1[imn]*phi_1_grid)
            Jac_1[isurf] += gmnc_1[isurf,imn]*np.cos(xm_nyq_1[imn]*theta_1_grid - xn_nyq_1[imn]*phi_1_grid)

        fsa_mod_B_1[isurf] = flux_surface_average(mod_B_1[isurf,:,:], Jac_1[isurf,:,:], theta_1, phi_1)
        fsa_B_sub_theta_1[isurf] = flux_surface_average(B_sub_theta_1[isurf,:,:], Jac_1[isurf,:,:], theta_1, phi_1)
        fsa_B_sub_phi_1[isurf] = flux_surface_average(B_sub_phi_1[isurf,:,:], Jac_1[isurf,:,:], theta_1, phi_1)
        fsa_B_sub_psi_1[isurf] = flux_surface_average(B_sub_psi_1[isurf,:,:], Jac_1[isurf,:,:], theta_1, phi_1)
        fsa_B_sup_theta_1[isurf] = flux_surface_average(B_sup_theta_1[isurf,:,:], Jac_1[isurf,:,:], theta_1, phi_1)
        fsa_B_sup_phi_1[isurf] = flux_surface_average(B_sup_phi_1[isurf,:,:], Jac_1[isurf,:,:], theta_1, phi_1)
        fsa_Jac_1[isurf] = flux_surface_average(Jac_1[isurf,:,:], Jac_1[isurf,:,:], theta_1, phi_1)
        dV_ds_1[isurf] = flux_surface_average(Jac_1[isurf,:,:], Jac_1[isurf,:,:], theta_1, phi_1, volume=True)
        
    for isurf in range(nsurf_2):
        for imn in range(mnmax_nyq_2):
            mod_B_2[isurf] += bmnc_2[isurf,imn]*np.cos(xm_nyq_2[imn]*theta_2_grid - xn_nyq_2[imn]*phi_2_grid)
            B_sub_theta_2[isurf] += bsubumnc_2[isurf,imn]*np.cos(xm_nyq_2[imn]*theta_2_grid - xn_nyq_2[imn]*phi_2_grid)
            B_sub_phi_2[isurf] += bsubvmnc_2[isurf,imn]*np.cos(xm_nyq_2[imn]*theta_2_grid - xn_nyq_2[imn]*phi_2_grid)
            B_sub_psi_2[isurf] += bsubsmns_2[isurf,imn]*np.sin(xm_nyq_2[imn]*theta_2_grid - xn_nyq_2[imn]*phi_2_grid)
            B_sup_theta_2[isurf] += bsupumnc_2[isurf,imn]*np.cos(xm_nyq_2[imn]*theta_2_grid - xn_nyq_2[imn]*phi_2_grid)
            B_sup_phi_2[isurf] += bsupvmnc_2[isurf,imn]*np.cos(xm_nyq_2[imn]*theta_2_grid - xn_nyq_2[imn]*phi_2_grid)
            Jac_2[isurf] += gmnc_2[isurf,imn]*np.cos(xm_nyq_2[imn]*theta_2_grid - xn_nyq_2[imn]*phi_2_grid)
        fsa_mod_B_2[isurf] = flux_surface_average(mod_B_2[isurf,:,:], Jac_2[isurf,:,:], theta_2, phi_2)
        fsa_B_sub_theta_2[isurf] = flux_surface_average(B_sub_theta_2[isurf,:,:], Jac_2[isurf,:,:], theta_2, phi_2)
        fsa_B_sub_phi_2[isurf] = flux_surface_average(B_sub_phi_2[isurf,:,:], Jac_2[isurf,:,:], theta_2, phi_2)
        fsa_B_sub_psi_2[isurf] = flux_surface_average(B_sub_psi_2[isurf,:,:], Jac_2[isurf,:,:], theta_2, phi_2)
        fsa_B_sup_theta_2[isurf] = flux_surface_average(B_sup_theta_2[isurf,:,:], Jac_2[isurf,:,:], theta_2, phi_2)
        fsa_B_sup_phi_2[isurf] = flux_surface_average(B_sup_phi_2[isurf,:,:], Jac_2[isurf,:,:], theta_2, phi_2)
        fsa_Jac_2[isurf] = flux_surface_average(Jac_2[isurf,:,:], Jac_2[isurf,:,:], theta_2, phi_2)
        dV_ds_2[isurf] = flux_surface_average(Jac_2[isurf,:,:], Jac_2[isurf,:,:], theta_2, phi_2, volume=True)


            
    fig, ax = plt.subplots(2,3)
    ax[0,0].plot(tflux_1[1:], fsa_mod_B_1[1:], label=r"{}".format(vmec1_descriptor))
    ax[0,0].plot(tflux_2[1:], fsa_mod_B_2[1:], label=r"{}".format(vmec2_descriptor))
    ax[0,0].set_title(r"$<|B|>$")
    ax[0,0].set_xlabel(r"$\psi_t$")
    
    ax[0,1].plot(tflux_1[1:], fsa_B_sub_theta_1[1:], label=r"{}".format(vmec1_descriptor))
    ax[0,1].plot(tflux_2[1:], fsa_B_sub_theta_2[1:], label=r"{}".format(vmec2_descriptor))
    ax[0,1].set_title(r"$<B_{\theta}>$")
    ax[0,1].set_xlabel(r"$\psi_t$")
    
    ax[0,2].plot(tflux_1[1:], fsa_B_sub_phi_1[1:], label=r"{}".format(vmec1_descriptor))
    ax[0,2].plot(tflux_2[1:], fsa_B_sub_phi_2[1:], label=r"{}".format(vmec2_descriptor))
    ax[0,2].set_title(r"$<B_{\phi}>$")
    ax[0,2].set_xlabel(r"$\psi_t$")
    
    ax[1,0].plot(tflux_1[1:-1], dV_ds_1[1:-1], label=r"{}".format(vmec1_descriptor))
    ax[1,0].plot(tflux_2[1:-1], dV_ds_2[1:-1], label=r"{}".format(vmec2_descriptor))
    ax[1,0].set_title(r"$dV/ds$")
    ax[1,0].set_xlabel(r"$\psi_t$")
    
    ax[1,1].plot(tflux_1[1:], fsa_B_sup_theta_1[1:], label=r"{}".format(vmec1_descriptor))
    ax[1,1].plot(tflux_2[1:], fsa_B_sup_theta_2[1:], label=r"{}".format(vmec2_descriptor))
    ax[1,1].set_title(r"$<B^{\theta}>$")
    ax[1,1].set_xlabel(r"$\psi_t$")
    
    ax[1,2].plot(tflux_1[1:], fsa_B_sup_phi_1[1:], label=r"{}".format(vmec1_descriptor))
    ax[1,2].plot(tflux_2[1:], fsa_B_sup_phi_2[1:], label=r"{}".format(vmec2_descriptor))
    ax[1,2].set_title(r"$<B^{\phi}>$")
    ax[1,2].set_xlabel(r"$\psi_t$")
    
    plt.legend()
    plt.show()
    exit()


def flux_surface_average(quantity, jacobian, theta, phi, volume=False):

    denom = np.trapz(np.trapz(jacobian, x=theta, axis=0), x=phi, axis=0)
    numerator = np.trapz(np.trapz(quantity*jacobian, x=theta, axis=0), x=phi, axis=0)
    fsa_quantity = numerator / denom
    if (volume):
        return denom
    else:
        return fsa_quantity


vmec1 = sys.argv[1]
vmec2 = sys.argv[2]
ntheta = 200
nphi = 32

compare_FSA_quantities(vmec1, vmec2, ntheta, nphi, vmec1_descriptor="fusion-io VMEC", vmec2_descriptor="original VMEC")
    
