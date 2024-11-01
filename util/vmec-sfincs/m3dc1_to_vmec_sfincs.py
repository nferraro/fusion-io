#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%                                                                              %
#%                     Copyright (C) 2024 Thea Energy, Inc                      %
#%                           All rights reserved                                %
#%                                                                              %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import os, sys, glob
import numpy as np
from netCDF4 import Dataset
import pickle
from scipy.interpolate import CubicSpline
import subprocess

def m3dc1_to_vmec(
        m3dc1_outfile,
        timeslice,
        nfp,
        ntheta,
        nphi,
        nsurf_nominal,
        vmec_mpol,
        vmec_ntor,
        te_start=-1,
        te_end=-1,
        psi_start = -1,
        psi_end = -1,
        bootstrap = 0,
        usingT = 0,
        tol=1,
        R0=0.0,
        write_sfincs_input=True,
        vmec_file_descriptor=None
):

    '''
    MAIN DRIVING FUNCTION

    This function takes a C1.h5 output file from M3D-C1 and outputs an imitation VMEC file that can be used with SFINCS

    Input Variables:
    * m3dc1_outfile: C1.h5 output from M3D-C1
    * timeslice (int): Timeslice from M3D-C1 output (-1 is first timeslice)
    * nfp (int): Number of field periods in equilibrium
    * ntheta (int): Number of poloidal gridpoints for write_neo_input
    * nphi (int): Number of toroidal gridpoints for write_neo_input
    * vmec_mpol (int): Max poloidal Fourier resolution for VMEC reconstruction
    * vmec_ntor (int): Max toroidal Fourier resolution for VMEC reconstruction
    * te_start (float): Starting electron temperature for write_neo_input
    * te_end (float): End electron temperature for write_neo_input
    * psi_start (float): Starting psi_normal for write_neo_input
    * psi_end (float): End psi_normal for write_neo_input
    * usingT (int): flag to use psi or te for write_neo_input
    * bootstrap (int): flag to use bootstrap current output in write_neo_input if the model is on in the m3dc1 simulations 
    * tol (float): Tolerance for write_neo_input solve
    * dR0 (float): Guess for magnetic axis for write_neo_input
    * write_sfincs_input (bool): If True, create a SFINCS input file and corresponding runspec.dat file
    * vmec_file_descriptor (string): String to include in the imitation VMEC output filename

    Returns:
    * neo_input.nc: Output from running write_neo_input
    * imVMEC*.nc: Imitation VMEC wout file
    * imVMEC_data.pkl: Pickled dataset of imitation VMEC data
    * input.namelist: SFINCS input file (if write_sfincs_input == True)
    * runspec.dat: Data with profile information for a radial sfincsScan (if write_sfincs_input == True)

    '''

    # fusion-io DIR
    fio_bin = os.path.join(os.environ['FIO_INSTALL_DIR'], "bin")
    
    # Double VMEC resolution to get half and full mesh
    nsurf = 2*nsurf_nominal
    phi = np.linspace(0, 2*np.pi, nphi, endpoint=False)
    theta = np.linspace(0, 2*np.pi, ntheta, endpoint=False)
    if (vmec_file_descriptor) == None:
        imVMEC_outname = "wout_imVMEC.nc"
    else:
        imVMEC_outname = "wout_imVMEC-{}.nc".format(vmec_file_descriptor)

    
    if (os.path.exists("neo_input.nc")):
        print("neo_input.nc file already exists!")
    else:
        # fusion-io/install/bin should be in PATH
        
        if (os.path.exists(os.path.join(fio_bin, "write_neo_input"))):
            if (usingT ==1):
                cmd = "~/src/fusion-io/util/_stellar/write_neo_input -bootstrap {} -m3dc1 {} {} -te_start {} -te_end {} -nr {} -nphi {} -ntheta {} -R0 {} -tol {}".format(bootstrap,m3dc1_outfile, timeslice, te_start, te_end, nsurf, nphi, ntheta, R0, tol)
                subprocess.run(cmd, cwd=os.getcwd(), shell=True)
            else:
                cmd = "~/src/fusion-io/util/_stellar/write_neo_input -bootstrap {} -m3dc1 {} {} -psi_start {} -psi_end {} -nr {} -nphi {} -ntheta {} -R0 {} -tol {}".format(bootstrap,m3dc1_outfile, timeslice, psi_start, psi_end, nsurf, nphi, ntheta, R0, tol)
                subprocess.run(cmd, cwd=os.getcwd(), shell=True)
        else:
            print("Unable to find write_neo_input executable. Check that FIO_ROOT is set in your environment. Exiting...")
            exit()

    # if imitation VMEC file already exists, skip
    if (os.path.exists(imVMEC_outname) and os.path.exists("imVMEC_data.pkl")):
        print("Imitation VMEC file already exists!")
        pickled_imVMEC_data = open("imVMEC_data.pkl", "rb")
        imVMEC_data = pickle.load(pickled_imVMEC_data)

    else:
        
        modB = np.zeros((nsurf,nphi,ntheta))
        B_sub_psi = np.zeros((nsurf,nphi,ntheta))
        B_sub_theta = np.zeros((nsurf,nphi,ntheta))
        B_sub_phi = np.zeros((nsurf,nphi,ntheta))
        B_sup_theta = np.zeros((nsurf,nphi,ntheta))
        B_sup_phi = np.zeros((nsurf,nphi,ntheta))
        Jac_vmec = np.zeros((nsurf,nphi,ntheta))
        avg_minor_radius = [-1.0]
    
        calculate_field_in_vmec_coordinates("neo_input.nc", modB, B_sub_psi, B_sub_theta, B_sub_phi, B_sup_theta, B_sup_phi, Jac_vmec, nsurf, theta, phi, avg_minor_radius)

        imVMEC_data = put_quantities_in_vmec_format("neo_input.nc", modB, B_sub_psi, B_sub_theta, B_sub_phi, B_sup_theta, B_sup_phi, Jac_vmec, nsurf, theta, phi, vmec_mpol, vmec_ntor, nfp, avg_minor_radius)
    
        with open("imVMEC_data.pkl", "wb") as fp:
            pickle.dump(imVMEC_data, fp)

        write_imitation_VMEC(imVMEC_outname, imVMEC_data, s_start=0.0, s_end=1.0)
       

    # path stuff
    if (write_sfincs_input):
        toroidal_flux = 0.45
        write_SFINCS_input(imVMEC_data=imVMEC_data, imitation_vmec_file=imVMEC_outname, input_type="radial", )


        

def put_quantities_in_vmec_format(neo_input_file, modB, B_sub_psi, B_sub_theta, B_sub_phi, B_sup_theta, B_sup_phi, Jac_vmec, nsurf, theta, phi, vmec_mpol, vmec_ntor, nfp, avg_minor_radius):

    '''
    This function converts quantities from fusion-io and write_neo_input into the variables and formats used in VMEC. This is mainly performing FFTs on the geometry of the flux surfaces and on the various components of the magnetic field

    * neo_input_file: The netcdf file generated from running the write_neo_input script
    * modB, B_sub_psi, B_sub_theta, B_sub_phi, B_sup_theta, B_sup_phi: empty arrays of size 
    '''
    
    neo_file = Dataset(neo_input_file, mode="r")
    R = np.swapaxes(neo_file.variables['R'][:,:,:], 0,2)
    Z = np.swapaxes(neo_file.variables['Z'][:,:,:], 0,2)
    #Jac = np.swapaxes(neo_file.variables['Jac'][:,:,:], 0,2)
    B_R = np.swapaxes(neo_file.variables['B_R'][:,:,:], 0,2)
    B_Phi = np.swapaxes(neo_file.variables['B_Phi'][:,:,:], 0,2)
    B_Z = np.swapaxes(neo_file.variables['B_Z'][:,:,:], 0,2)
    psi_t = neo_file.variables['Psi_t'][:]
    q = neo_file.variables['q'][:]
    ne = neo_file.variables['ne0'][:]
    te = neo_file.variables['Te0'][:]
    ni = neo_file.variables['ni0'][:]
    ti = neo_file.variables['Ti0'][:]
    
    ntheta = len(theta)
    nphi = len(phi)
    theta_periodic = np.append(theta, 2*np.pi)
    phi_periodic = np.append(phi, 2*np.pi)
    
    xm, xn, xm_nyq, xn_nyq = create_xm_xn_arrays(vmec_mpol, vmec_ntor, nfp, nyq=True)
    mnmax = len(xm)
    mnmax_nyq = len(xm_nyq)
    iota = 1./q
    
    elementary_charge = 1.602176634e-19
    pressure = elementary_charge * (ne*te + ni*ti) # Pressure in Pa
    
    phi_vmec = psi_t
    phi_vmec_lcfs = phi_vmec[-1]
    s = phi_vmec / phi_vmec_lcfs
    dphids = np.zeros(nsurf)
    for i in range(nsurf-1):
        dphids[i+1] = (-1./(2*np.pi))*(phi_vmec[i+1]-phi_vmec[i])/(s[i+1]-s[i])

    nsurf_nominal = nsurf // 2
        
    imVMEC_data = {'nfp': nfp}
    imVMEC_data['mpol'] = vmec_mpol
    imVMEC_data['ntor'] = vmec_ntor
    imVMEC_data['xm'] = xm
    imVMEC_data['xn'] = xn
    imVMEC_data['mnmax'] = mnmax
    imVMEC_data['xm_nyq'] = xm_nyq
    imVMEC_data['xn_nyq'] = xn_nyq
    imVMEC_data['mnmax_nyq'] = mnmax_nyq
    imVMEC_data['nsurf'] = nsurf_nominal
    imVMEC_data['stell_sym'] = 0
    imVMEC_data['aminor_avg'] = avg_minor_radius[0]
    imVMEC_data['ne'] = np.zeros(nsurf_nominal)
    imVMEC_data['Te'] = np.zeros(nsurf_nominal)
    imVMEC_data['ni'] = np.zeros(nsurf_nominal)
    imVMEC_data['Ti'] = np.zeros(nsurf_nominal)
    imVMEC_data['presf'] = np.zeros((nsurf_nominal))
    imVMEC_data['phi'] = np.zeros((nsurf_nominal))
    imVMEC_data['phips'] = np.zeros((nsurf_nominal))
    imVMEC_data['iotas'] = np.zeros((nsurf_nominal))
    imVMEC_data['rmnc'] = np.zeros((nsurf_nominal, mnmax))
    imVMEC_data['rmns'] = np.zeros((nsurf_nominal, mnmax))
    imVMEC_data['zmnc'] = np.zeros((nsurf_nominal, mnmax))
    imVMEC_data['zmns'] = np.zeros((nsurf_nominal, mnmax))
    imVMEC_data['gmnc'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['gmns'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['bmnc'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['bmns'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['bsupumnc'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['bsupumns'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['bsupvmnc'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['bsupvmns'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['bsubsmnc'] = np.zeros((nsurf_nominal, mnmax_nyq)) # full mesh
    imVMEC_data['bsubsmns'] = np.zeros((nsurf_nominal, mnmax_nyq)) # full mesh
    imVMEC_data['bsubumnc'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['bsubumns'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['bsubvmnc'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['bsubvmns'] = np.zeros((nsurf_nominal, mnmax_nyq))
    imVMEC_data['fio_ntheta'] = ntheta
    imVMEC_data['fio_nphi'] = nphi   

    # It is assumed that of the nsurf surfaces computed, the even surfaces will be considered the 'full mesh' and the odd surfaces will be considered the 'half mesh'    
    for isurf in range(nsurf_nominal):
        half_mesh_idx = 2*isurf+1
        full_mesh_idx = 2*isurf
        imVMEC_data['presf'][isurf] = pressure[full_mesh_idx]
        imVMEC_data['phi'][isurf] = phi_vmec[full_mesh_idx]
        imVMEC_data['iotas'][isurf] = iota[half_mesh_idx]
        imVMEC_data['phips'][isurf] = dphids[half_mesh_idx]
        imVMEC_data['ne'][isurf] = ne[full_mesh_idx]
        imVMEC_data['Te'][isurf] = te[full_mesh_idx]
        imVMEC_data['ni'][isurf] = ni[full_mesh_idx]
        imVMEC_data['Ti'][isurf] = ti[full_mesh_idx]
        imVMEC_data['rmnc'][isurf,:], imVMEC_data['rmns'][isurf,:] = fourier_transform_2d(R[half_mesh_idx,:,:], theta_periodic, phi_periodic, xm, xn)
        imVMEC_data['zmnc'][isurf,:], imVMEC_data['zmns'][isurf,:] = fourier_transform_2d(Z[half_mesh_idx,:,:], theta_periodic, phi_periodic, xm, xn)
        imVMEC_data['gmnc'][isurf,:], imVMEC_data['gmns'][isurf,:] = fourier_transform_2d(Jac_vmec[half_mesh_idx,:,:], theta_periodic, phi_periodic, xm_nyq, xn_nyq)
        imVMEC_data['bmnc'][isurf,:], imVMEC_data['bmns'][isurf,:] = fourier_transform_2d(modB[half_mesh_idx,:,:], theta_periodic, phi_periodic, xm_nyq, xn_nyq)
        imVMEC_data['bsupumnc'][isurf,:], imVMEC_data['bsupumns'][isurf,:] = fourier_transform_2d(B_sup_theta[half_mesh_idx,:,:], theta_periodic, phi_periodic, xm_nyq, xn_nyq)
        imVMEC_data['bsupvmnc'][isurf,:], imVMEC_data['bsupvmns'][isurf,:] = fourier_transform_2d(B_sup_phi[half_mesh_idx,:,:], theta_periodic, phi_periodic, xm_nyq, xn_nyq)
        imVMEC_data['bsubumnc'][isurf,:], imVMEC_data['bsubumns'][isurf,:] = fourier_transform_2d(B_sub_theta[half_mesh_idx,:,:], theta_periodic, phi_periodic, xm_nyq, xn_nyq)
        imVMEC_data['bsubvmnc'][isurf,:], imVMEC_data['bsubvmns'][isurf,:] = fourier_transform_2d(B_sub_phi[half_mesh_idx,:,:], theta_periodic, phi_periodic, xm_nyq, xn_nyq)
        imVMEC_data['bsubsmnc'][isurf,:], imVMEC_data['bsubsmns'][isurf,:] = fourier_transform_2d(B_sub_psi[full_mesh_idx,:,:], theta_periodic, phi_periodic, xm_nyq, xn_nyq) # Full mesh
        '''
        imVMEC_data['aux'][isurf,:,:] = Jac[half_mesh_idx,:,:]
        imVMEC_data['B_sub_theta'][isurf,:,:] = B_sub_theta[half_mesh_idx,:,:]
        imVMEC_data['B_sub_psi'][isurf,:,:] = B_sub_psi[half_mesh_idx,:,:]
        imVMEC_data['B_sup_theta'][isurf,:,:] = B_sup_theta[half_mesh_idx,:,:]
        imVMEC_data['B_R'][isurf,:,:] = B_R[half_mesh_idx,:,:]
        imVMEC_data['B_Phi'][isurf,:,:] = B_Phi[half_mesh_idx,:,:]
        imVMEC_data['B_Z'][isurf,:,:] = B_Z[half_mesh_idx,:,:]
        '''
    return imVMEC_data

def calculate_field_in_vmec_coordinates(neo_input_file, modB, B_sub_psi, B_sub_theta, B_sub_phi, B_sup_theta, B_sup_phi, Jac_vmec, nsurf, theta, phi, avg_minor_radius):

    neo_file = Dataset(neo_input_file, mode="r")
    R = np.swapaxes(neo_file.variables['R'][:,:,:], 0,2)
    Z = np.swapaxes(neo_file.variables['Z'][:,:,:], 0,2)
    Jac = np.swapaxes(neo_file.variables['Jac'][:,:,:], 0,2)
    B_R = np.swapaxes(neo_file.variables['B_R'][:,:,:], 0,2)
    B_Phi = np.swapaxes(neo_file.variables['B_Phi'][:,:,:], 0,2)
    B_Z = np.swapaxes(neo_file.variables['B_Z'][:,:,:], 0,2)
    psi_t = neo_file.variables['Psi_t'][:]
    q = neo_file.variables['q'][:]
    ne = neo_file.variables['ne0'][:]
    te = neo_file.variables['Te0'][:]
    ni = neo_file.variables['ni0'][:]
    ti = neo_file.variables['Ti0'][:]

    ntheta = len(theta)
    nphi = len(phi)
    dtheta = theta[1]-theta[0]
    dphi = phi[1]-phi[0]
    
    dR_dpsi = np.zeros((nsurf,nphi,ntheta))
    dZ_dpsi = np.zeros((nsurf,nphi,ntheta))
    dR_dtheta = np.zeros((nsurf,nphi,ntheta))
    dZ_dtheta = np.zeros((nsurf,nphi,ntheta))
    dR_dphi = np.zeros((nsurf,nphi,ntheta))
    dZ_dphi = np.zeros((nsurf,nphi,ntheta))
    
    # Converting Jacobian derivates w.r.t. index from write_neo_input into VMEC coordinates

    np.copyto(Jac_vmec,Jac)
    Jac_vmec /= -(dtheta*dphi)
    psi_t_end = psi_t[-1]
    psi_t_norm = psi_t / psi_t_end

    #psi_t_end = psi_t[-1]   # If comparing to a VMEC file, hard code phi_t_LCFS to phi[-1] from the VMEC file
    #psi_t_norm = psi_t / 4.4245216949037 # Need to normalize wrt toroidal flux on boundary of original vmec file for comparison purposes

    psi_t_end=1#0.5144
    psi_t_norm = psi_t /psi_t_end # 1 for benchmarking case, 0.5144 for NCSX case, 0.4 for tcv comparison case If comparing to a VMEC file, hard code phi_t_LCFS to phi[-1] from the VMEC file
    print('psi_t_end')
    print(psi_t_end)
    for isurf in range(nsurf):
        for iphi in range(nphi):
            for itheta in range(ntheta):

                surf_idx_pos = isurf if isurf==nsurf-1 else isurf+1
                surf_idx_neg = 0 if isurf==0 else isurf-1
                numer_fac = 1.0 if isurf==0 or isurf==nsurf-1 else 2.0
                denom_fac = 1.0 if isurf==0 or isurf==nsurf-1 else 1.0
                dR_dpsi[isurf,iphi,itheta] = (R[surf_idx_pos,iphi,itheta] - R[surf_idx_neg,iphi,itheta])/(denom_fac*(psi_t_norm[surf_idx_pos]-psi_t_norm[surf_idx_neg]))
                dZ_dpsi[isurf,iphi,itheta] = (Z[surf_idx_neg,iphi,itheta] - Z[surf_idx_neg,iphi,itheta])/(denom_fac*(psi_t_norm[surf_idx_pos]-psi_t_norm[surf_idx_neg]))
                Jac_vmec[isurf,iphi,itheta] *= numer_fac/(denom_fac*(psi_t_norm[surf_idx_pos]-psi_t_norm[surf_idx_neg]))

                phi_idx_pos = 0 if iphi==nphi-1 else iphi+1
                phi_idx_neg = -1 if iphi==0 else iphi-1
                dR_dphi[isurf,iphi,itheta] = (R[isurf,phi_idx_pos,itheta] - R[isurf,phi_idx_neg,itheta])/(2.*dphi)
                dZ_dphi[isurf,iphi,itheta] = (Z[isurf,phi_idx_pos,itheta] - Z[isurf,phi_idx_neg,itheta])/(2.*dphi)

                theta_idx_pos = 0 if itheta==ntheta-1 else itheta+1
                theta_idx_neg = -1 if itheta==0 else itheta-1
                dR_dtheta[isurf,iphi,itheta] = (R[isurf,iphi,theta_idx_pos] - R[isurf,iphi,theta_idx_neg])/(2.*dtheta)
                dZ_dtheta[isurf,iphi,itheta] = (Z[isurf,iphi,theta_idx_pos] - Z[isurf,iphi,theta_idx_neg])/(2.*dtheta)
                

                # Covariant components of magnetic field in flux coordinates
                B_sub_psi[isurf,iphi,itheta] = B_R[isurf,iphi,itheta]*dR_dpsi[isurf,iphi,itheta] + B_Z[isurf,iphi,itheta]*dZ_dpsi[isurf,iphi,itheta]
                B_sub_theta[isurf,iphi,itheta] = B_R[isurf,iphi,itheta]*dR_dtheta[isurf,iphi,itheta] + B_Z[isurf,iphi,itheta]*dZ_dtheta[isurf,iphi,itheta]# + B_Phi[isurf,iphi,itheta]*dphi*dR_dtheta[isurf,iphi,itheta]
                B_sub_phi[isurf,iphi,itheta] = B_R[isurf,iphi,itheta]*dR_dphi[isurf,iphi,itheta] + B_Phi[isurf,iphi,itheta]*R[isurf,iphi,itheta] + B_Z[isurf,iphi,itheta]*dZ_dphi[isurf,iphi,itheta]

                # Metric tensor elements to calculate contravariant components
                g_sub_phi_theta = dR_dtheta[isurf,iphi,itheta]*dR_dphi[isurf,iphi,itheta] + dZ_dtheta[isurf,iphi,itheta]*dZ_dphi[isurf,iphi,itheta]
                g_sub_phi_phi = np.power(dR_dphi[isurf,iphi,itheta],2) + np.power(dZ_dphi[isurf,iphi,itheta],2) + np.power(R[isurf,iphi,itheta],2)
                g_sub_theta_theta = np.power(dR_dtheta[isurf,iphi,itheta],2) + np.power(dZ_dtheta[isurf,iphi,itheta],2)

                # Contravariant components of magnetic field in flux coordinates
                B_sup_theta[isurf,iphi,itheta] = (B_sub_theta[isurf,iphi,itheta] - (g_sub_phi_theta/g_sub_phi_phi)*B_sub_phi[isurf,iphi,itheta]) / (g_sub_theta_theta - np.power(g_sub_phi_theta,2)/g_sub_phi_phi)
                B_sup_phi[isurf,iphi,itheta] = (B_sub_phi[isurf,iphi,itheta] - g_sub_phi_theta*B_sup_theta[isurf,iphi,itheta])/g_sub_phi_phi

                # |B|
                modB[isurf,iphi,itheta] = np.sqrt(B_sub_theta[isurf,iphi,itheta]*B_sup_theta[isurf,iphi,itheta] + B_sub_phi[isurf,iphi,itheta]*B_sup_phi[isurf,iphi,itheta])


    # Calculate effective minor radius using geometry at outermost surface
    dZ_dtheta_edge = dZ_dtheta[-1,:,:]
    R_edge = R[-1,:,:]
    Abar = np.trapz(np.trapz(R_edge*dZ_dtheta_edge, x=phi, axis=0), x=theta, axis=0) / (2*np.pi)
    avg_minor_radius[0] = np.sqrt(Abar / np.pi)

    if (avg_minor_radius[0] < 0):
        print("Average minor radius was either calculated incorrectly or there is an issue with the geometry of the outermost surface. Exiting...")
        exit()

                
def write_SFINCS_runspec(psiN, ne_fio, Te_fio, ni_fio, Ti_fio, Tbar=1000, nbar=1.e20, vmec_compare=False, s_vmec=None):

    '''
    This function generates the SFINCS runspec file that is needed for a radial scan of SFINCS runs

    * ne_fio: Electron density from fusion-io
    * ni_fio: Ion density from fusion-io
    * Te_fio: Electron temperature from fusion-io
    * Ti_fio: Ion temperature from fusion-io
    * Tbar: Temperature normalization for SFINCS (default: 1000 eV)
    * nbar: Density normalization for SFINCS (default: 1.e20 m^{-3})

    The vmec_compare / s_vmec option are used for generating a 1-to-1 comparison with the FIO surfaces, which do not necessarily have the same range of toroidal flux values
    '''    

    if (vmec_compare):
        nradial = len(s_vmec)
    else:
        nradial = len(psiN)
        
    # Create splines from the data
    Ti = CubicSpline(psiN, Ti_fio/Tbar)
    dTi_dpsiN = Ti.derivative()
    Te = CubicSpline(psiN, Te_fio/Tbar)
    dTe_dpsiN = Te.derivative()
    ne = CubicSpline(psiN, ne_fio/nbar)
    dne_dpsiN = ne.derivative()    

    run_spec = open("runspec.dat", 'w')
    run_spec.write("! psiN_wish nHats_1 nHats_2 THats_1 THats_2 dnHatdpsiNs_1 dnHatdpsiNs_2 dTHatdpsiNs_1 dTHatdpsiNs_2\n")
    for irad in range(1,nradial):

        if (vmec_compare):
            line = "{} {} {} {} {} {} {} {} {}\n".format(s_vmec[irad], ne(s_vmec[irad]), ne(s_vmec[irad]), Te(s_vmec[irad]), Ti(s_vmec[irad]), dne_dpsiN(s_vmec[irad]), dne_dpsiN(s_vmec[irad]), dTe_dpsiN(s_vmec[irad]), dTi_dpsiN(s_vmec[irad]))
        else:
            line = "{} {} {} {} {} {} {} {} {}\n".format(psiN[irad], ne(psiN[irad]), ne(psiN[irad]), Te(psiN[irad]), Ti(psiN[irad]), dne_dpsiN(psiN[irad]), dne_dpsiN(psiN[irad]), dTe_dpsiN(psiN[irad]), dTi_dpsiN(psiN[irad]))
            
        run_spec.write(line)
        

    run_spec.close()
    
    
def write_SFINCS_input(imVMEC_data, imitation_vmec_file, input_type="surface", psiN=0.5, Tbar=1000, nbar=1.e20, vmec_compare=False, vmec_file=None):

    '''
    Based on profile information and VMEC output file, this function will output a SFINCS input file
    '''

    if (input_type == "surface"):
        print("Creating SFINCS input file for a single flux surface")
    elif (input_type == "radial"):
        print("Creating SFINCS input file and runspec.dat for a radial scan")
    else:
        print("SFINCS input_type must be either 'surface' or 'radial'. Exiting...")
        exit()
        
    desired_psi_t = psiN * imVMEC_data['phi'][-1]
    s = imVMEC_data['phi'] / imVMEC_data['phi'][-1]

    if (input_type == "radial"):
        if (vmec_compare):
            vmec_data = Dataset(sys.argv[1], mode='r')
            phi_vmec = vmec_data.variables['phi'][:]
            s_vmec = phi_vmec / phi_vmec[-1]
            idx = find_nearest_idx(phi_vmec, imVMEC_data['phi'][-1])
            s_vmec = s_vmec[:idx]
            write_SFINCS_runspec(s, imVMEC_data['ne'], imVMEC_data['Te'], imVMEC_data['ni'], imVMEC_data['Ti'], vmec_compare=True, s_vmec=s_vmec)
            exit()
        else:
            write_SFINCS_runspec(s, imVMEC_data['ne'], imVMEC_data['Te'], imVMEC_data['ni'], imVMEC_data['Ti'])        
    
    # Select flux surface closest to desired normalized toroidal flux input value
    id1, id2 = find_closest_indices(imVMEC_data['phi'], desired_psi_t)
    isurf = id1

    ne = imVMEC_data['ne'][id1] + (desired_psi_t - imVMEC_data['phi'][id1])*( (imVMEC_data['ne'][id2] - imVMEC_data['ne'][id1]) / ( imVMEC_data['phi'][id2] - imVMEC_data['phi'][id1] ))
    Te = imVMEC_data['Te'][id1] + (desired_psi_t - imVMEC_data['phi'][id1])*( (imVMEC_data['Te'][id2] - imVMEC_data['Te'][id1]) / ( imVMEC_data['phi'][id2] - imVMEC_data['phi'][id1] ))
    ni = imVMEC_data['ni'][id1] + (desired_psi_t - imVMEC_data['phi'][id1])*( (imVMEC_data['ni'][id2] - imVMEC_data['ni'][id1]) / ( imVMEC_data['phi'][id2] - imVMEC_data['phi'][id1] ))
    Ti = imVMEC_data['Ti'][id1] + (desired_psi_t - imVMEC_data['phi'][id1])*( (imVMEC_data['Ti'][id2] - imVMEC_data['Ti'][id1]) / ( imVMEC_data['phi'][id2] - imVMEC_data['phi'][id1] ))
    
    # From the equilibrium data, need to calculate derivatives of temp and density
    # Use SFINCS dnHatdpsiNs so derivative is w.r.t. the normalized toroidal flux
    dne_dpsi = (imVMEC_data['ne'][isurf+1] - imVMEC_data['ne'][isurf-1]) / (s[isurf+1] - s[isurf-1])
    dTe_dpsi = (imVMEC_data['Te'][isurf+1] - imVMEC_data['Te'][isurf-1]) / (s[isurf+1] - s[isurf-1])
    dni_dpsi = (imVMEC_data['ni'][isurf+1] - imVMEC_data['ni'][isurf-1]) / (s[isurf+1] - s[isurf-1])
    dTi_dpsi = (imVMEC_data['Ti'][isurf+1] - imVMEC_data['Ti'][isurf-1]) / (s[isurf+1] - s[isurf-1])

    # Write out SFINCS input file
    infile = open(os.path.join(os.environ['FIO_ROOT'],'util/vmec-sfincs/sfincs_input_template'),'r')
    outfile = open('./input.namelist','w')
    line = infile.readline()
    while line:
        if "!ss scanType" in line:
            if (input_type == "surface"):
                outfile.write("!ss scanType = 1\n")
            elif (input_type == "radial"):
                outfile.write("!ss scanType = 21\n")

        elif "psiN_wish =" in line:
            # This value is only used if doing a "surface" run
            outfile.write("  psiN_wish = {} ! Set by sfincsScan.\n".format(psiN))

        elif "equilibriumFile" in line:
            outfile.write("  equilibriumFile = '{}' ! this location is correct after this file is copied to subdirs\n".format(imitation_vmec_file))

        elif "nHats" in line:
            outfile.write("  nHats = {}  {}\n".format(ne/nbar, ni/nbar))

        elif "dnHatdpsiNs" in line:
            outfile.write("  dnHatdpsiNs = {}  {}\n".format(dne_dpsi/nbar, dni_dpsi/nbar))
        elif "THats" in line:
            outfile.write("  THats = {}  {}\n".format(Te/Tbar, Ti/Tbar))
            
        elif "dTHatdpsiNs" in line:
            outfile.write("  dTHatdpsiNs = {}  {}\n".format(dTe_dpsi/Tbar, dTi_dpsi/Tbar))
        elif "Zs =" in line:
            outfile.write("  Zs = {}  {}\n".format(-1, 1))

        elif "mHats =" in line:
            outfile.write("  mHats = 0.000548561  1.00723\n")

        else:
            outfile.write(line)

        line = infile.readline()

    outfile.close()
    infile.close()
    

def write_imitation_VMEC(outfile, imVMEC_data, s_start, s_end):

    '''
    This function will create a VMEC-like *.nc output file with only those fields required to initialize a SFINCS calculation
    (all other fields should probably be filled in with nulls or zeros...)
    '''

    head, tail = os.path.split(outfile)
    path = os.path.join(os.getcwd(), tail)
    
    file = Dataset(path, mode="w", format="NETCDF3_64BIT_OFFSET")

    s_start = imVMEC_data['phi'][0]
    # VMEC radial coordinate: s = rho^2 = Psi_t / Psi_t(LCFS)
    s_full = np.linspace(s_start, s_end, imVMEC_data['nsurf'])
    s_half = s_full[0:-1] + 0.5 / (imVMEC_data['nsurf'] - 1)
    
    # dimensions
    file.createDimension("radius", imVMEC_data['nsurf'])  # number of flux surfaces
    file.createDimension("fio_nphi", imVMEC_data['fio_nphi']) # number of toroidal grid points in fusion-io construction
    file.createDimension("fio_ntheta", imVMEC_data['fio_ntheta']) # number of poloidal grid points in fusion-io construction
    file.createDimension(
        "mn_mode", imVMEC_data['mnmax']
        #"mn_mode", (2 * N + 1)* M - N
    )  # number of Fourier modes
    file.createDimension(
        "mn_mode_nyq", imVMEC_data['mnmax_nyq']
        #"mn_mode_nyq", (2 * N_nyq + 1)* M_nyq - N_nyq
    )   # number of Nyquist Fourier modes
    file.createDimension("ntor_p1", imVMEC_data['ntor']+1)
    file.createDimension("dim_00001", 1)
    file.createDimension("dim_00020", 20)
    file.createDimension("dim_00100", 100)
    file.createDimension("dim_00200", 200)

    imitation_vmec = file.createVariable("IMITATION_VMEC", np.int32)
    imitation_vmec.long_name = "This dataset is meant to imitate a VMEC output file for use with the SFINCS neoclassical code. These data fields do NOT represent a VMEC equilibrium"
    imitation_vmec[:] = 1

    version_vmec = file.createVariable("version_", np.float64)
    version_vmec.long_name = "This is a fake version number that acts as a flag to tell SFINCS that this is a not a real VMEC file, and to handle it slightly differently"
    version_vmec[:] = -2.0
    
    lasym = file.createVariable("lasym__logical__", np.int32)
    lasym.long_name = "asymmetry logical (0 = stellarator symmetry)"
    lasym[:] = int(imVMEC_data['stell_sym'])
    
    nfp = file.createVariable("nfp", np.int32)
    nfp.long_name = "number of field periods"
    nfp[:] = imVMEC_data['nfp']

    ns = file.createVariable("ns", np.int32)
    ns.long_name = "number of flux surfaces"
    ns[:] = imVMEC_data['nsurf']

    mpol = file.createVariable("mpol", np.int32)
    mpol.long_name = "number of poloidal Fourier modes"
    mpol[:] = imVMEC_data['mpol']

    ntor = file.createVariable("ntor", np.int32)
    ntor.long_name = "number of positive toroidal Fourier modes"
    ntor[:] = imVMEC_data['ntor']

    mnmax = file.createVariable("mnmax", np.int32)
    mnmax.long_name = "total number of Fourier modes"
    mnmax[:] = file.dimensions["mn_mode"].size

    xm = file.createVariable("xm", np.float64, ("mn_mode",))
    xm.long_name = "poloidal mode numbers"
    xm[:] = imVMEC_data['xm']

    xn = file.createVariable("xn", np.float64, ("mn_mode",))
    xn.long_name = "toroidal mode numbers"
    xn[:] = imVMEC_data['xn']

    mnmax_nyq = file.createVariable("mnmax_nyq", np.int32)
    mnmax_nyq.long_name = "total number of Nyquist Fourier modes"
    mnmax_nyq[:] = file.dimensions["mn_mode_nyq"].size

    xm_nyq = file.createVariable("xm_nyq", np.float64, ("mn_mode_nyq",))
    xm_nyq.long_name = "poloidal Nyquist mode numbers"
    xm_nyq[:] = imVMEC_data['xm_nyq']

    xn_nyq = file.createVariable("xn_nyq", np.float64, ("mn_mode_nyq",))
    xn_nyq.long_name = "toroidal Nyquist mode numbers"
    xn_nyq[:] = imVMEC_data['xn_nyq']

    phi = file.createVariable("phi", np.float64, ("radius",))
    phi.long_name = "toroidal flux on full mesh"
    phi.units = "Wb"
    phi[:] = imVMEC_data['phi']

    iotas = file.createVariable("iotas", np.float64, ("radius",))
    iotas.long_name = "rotational transform on half mesh"
    iotas.units = "None"
    iotas[:] = imVMEC_data['iotas']
    #iotas[0] = 0
    #iotas[1:] = -i_half  # negative sign for negative Jacobian

    presf = file.createVariable("presf", np.float64, ("radius",))
    presf.long_name = "pressure on full mesh"
    presf.units = "Pa"
    presf[:] = imVMEC_data['presf']

    phips = file.createVariable("phips", np.float64, ("radius",))
    phips.long_name = "d(phi)/ds * -1/2pi: toroidal flux derivative on half mesh"
    phips[:] = imVMEC_data['phips']
    #phips[0] = 0
    #phips[1:] = -phipf[1:] / (2 * np.pi)

    signgs = file.createVariable("signgs", np.int32)
    signgs.long_name = "sign of coordinate system Jacobian"
    signgs[:] = -1  # VMEC always uses a negative Jacobian

    Aminor_p = file.createVariable("Aminor_p", np.float64)
    Aminor_p.long_name = "minor radius"
    Aminor_p.units = "m"
    Aminor_p[:] = imVMEC_data['aminor_avg']
    
    # R
    rmnc = file.createVariable("rmnc", np.float64, ("radius", "mn_mode"))
    rmnc.long_name = "cos(m*t-n*p) component of cylindrical R, on full mesh"
    rmnc.units = "m"
    rmnc[:,:] = imVMEC_data['rmnc'][:,:]
    if not imVMEC_data['stell_sym'] == 0:
        rmns = file.createVariable("rmns", np.float64, ("radius", "mn_mode"))
        rmns.long_name = "sin(m*t-n*p) component of cylindrical R, on full mesh"
        rmns.units = "m"
        rmns[:,:] = imVMEC_data['rmns'][:,:]

    # Z
    zmns = file.createVariable("zmns", np.float64, ("radius", "mn_mode"))
    zmns.long_name = "sin(m*t-n*p) component of cylindrical Z, on full mesh"
    zmns.units = "m"
    zmns[:,:] = imVMEC_data['zmns'][:,:]
    if not imVMEC_data['stell_sym'] == 0:
        zmnc = file.createVariable("zmnc", np.float64, ("radius", "mn_mode"))
        zmnc.long_name = "cos(m*t-n*p) component of cylindrical Z, on full mesh"
        zmnc.units = "m"
        zmnc[:,:] = imVMEC_data['zmnc'][:,:]
        
    # Jacobian
    gmnc = file.createVariable("gmnc", np.float64, ("radius", "mn_mode_nyq"))
    gmnc.long_name = "cos(m*t-n*p) component of Jacobian, on half mesh"
    gmnc.units = "m"
    gmnc[0, :] = 0
    gmnc[1:, :] = imVMEC_data['gmnc'][1:,:]
    if not imVMEC_data['stell_sym'] == 0:
        gmns = file.createVariable("gmns", np.float64, ("radius", "mn_mode_nyq"))
        gmns.long_name = "sin(m*t-n*p) component of Jacobian, on half mesh"
        gmns.units = "m"
        gmns[0, :] = 0
        gmns[1:, :] = imVMEC_data['gmns'][1:,:]

    # |B|
    bmnc = file.createVariable("bmnc", np.float64, ("radius", "mn_mode_nyq"))
    bmnc.long_name = "cos(m*t-n*p) component of |B|, on half mesh"
    bmnc.units = "T"
    bmnc[0, :] = 0
    bmnc[1:, :] = imVMEC_data['bmnc'][1:,:]
    if not imVMEC_data['stell_sym'] == 0:
        bmns = file.createVariable("bmns", np.float64, ("radius", "mn_mode_nyq"))
        bmns.long_name = "sin(m*t-n*p) component of |B|, on half mesh"
        bmns.units = "T"
        bmns[0, :] = 0
        bmns[1:, :] = imVMEC_data['bmns'][1:,:]

    # B^theta
    bsupumnc = file.createVariable("bsupumnc", np.float64, ("radius", "mn_mode_nyq"))
    bsupumnc.long_name = "cos(m*t-n*p) component of B^theta, on half mesh"
    bsupumnc.units = "T/m"
    bsupumnc[0, :] = 0
    bsupumnc[1:, :] = imVMEC_data['bsupumnc'][1:,:]
    if not imVMEC_data['stell_sym'] == 0:
        bsupumns = file.createVariable("bsupumns", np.float64, ("radius", "mn_mode_nyq"))
        bsupumns.long_name = "sin(m*t-n*p) component of B^theta, on half mesh"
        bsupumns.units = "T/m"
        bsupumns[0, :] = 0
        bsupumns[1:, :] = -imVMEC_data['bsupumns'][1:,:]

    # B^zeta
    bsupvmnc = file.createVariable("bsupvmnc", np.float64, ("radius", "mn_mode_nyq"))
    bsupvmnc.long_name = "cos(m*t-n*p) component of B^zeta, on half mesh"
    bsupvmnc.units = "T/m"
    bsupvmnc[0, :] = 0
    bsupvmnc[1:, :] = imVMEC_data['bsupvmnc'][1:,:]
    if not imVMEC_data['stell_sym'] == 0:
        bsupvmns = file.createVariable("bsupvmns", np.float64, ("radius", "mn_mode_nyq"))
        bsupvmns.long_name = "sin(m*t-n*p) component of B^zeta, on half mesh"
        bsupvmns.units = "T/m"
        bsupvmns[0, :] = 0
        bsupvmns[1:, :] = imVMEC_data['bsupvmns'][1:,:]

    # B_psi
    bsubsmns = file.createVariable("bsubsmns", np.float64, ("radius", "mn_mode_nyq"))
    bsubsmns.long_name = "sin(m*t-n*p) component of B_psi, on full mesh"
    bsubsmns.units = "T*m"
    bsubsmns[:] = imVMEC_data['bsubsmns']
    if not imVMEC_data['stell_sym'] == 0:
        bsubsmnc = file.createVariable("bsubsmnc", np.float64, ("radius", "mn_mode_nyq"))
        bsubsmnc.long_name = "cos(m*t-n*p) component of B_psi, on full mesh"
        bsubsmnc.units = "T*m"
        bsubsmnc[:] = imVMEC_data['bsubsmnc']

    # B_theta
    bsubumnc = file.createVariable("bsubumnc", np.float64, ("radius", "mn_mode_nyq"))
    bsubumnc.long_name = "cos(m*t-n*p) component of B_theta, on half mesh"
    bsubumnc.units = "T*m"
    bsubumnc[0, :] = 0
    bsubumnc[1:, :] = imVMEC_data['bsubumnc'][1:,:]
    if not imVMEC_data['stell_sym'] == 0:
        bsubumns = file.createVariable("bsubumns", np.float64, ("radius", "mn_mode_nyq"))
        bsubumns.long_name = "sin(m*t-n*p) component of B_theta, on half mesh"
        bsubumns.units = "T*m"
        bsubumns[0, :] = 0
        bsubumns[1:, :] = -imVMEC_data['bsubumns'][1:,:]

    # B_zeta
    bsubvmnc = file.createVariable("bsubvmnc", np.float64, ("radius", "mn_mode_nyq"))
    bsubvmnc.long_name = "cos(m*t-n*p) component of B_zeta, on half mesh"
    bsubvmnc.units = "T*m"
    bsubvmnc[0, :] = 0
    bsubvmnc[1:, :] = imVMEC_data['bsubvmnc'][1:,:]
    if not imVMEC_data['stell_sym'] == 0:
        bsubvmns = file.createVariable("bsubvmns", np.float64, ("radius", "mn_mode_nyq"))
        bsubvmns.long_name = "sin(m*t-n*p) component of B_zeta, on half mesh"
        bsubvmns.units = "T*m"
        bsubvmns[0, :] = 0
        bsubvmns[1:, :] = imVMEC_data['bsubvmns'][1:,:]

    # Start of creating empty arrays and auxiliary arrays
    raxis_cc = file.createVariable("raxis_cc", np.float64, ("ntor_p1",))
    raxis_cc[:] = 0.0

    zaxis_cs = file.createVariable("zaxis_cs", np.float64, ("ntor_p1",))
    zaxis_cs[:] = 0.0

    # Some potential debugging variables
    '''
    # B_sub_theta
    B_sub_theta = file.createVariable("B_sub_theta", np.float64, ("radius", "fio_nphi", "fio_ntheta"))
    B_sub_theta.long_name = "B_sub_theta"
    B_sub_theta.units = "m"
    B_sub_theta[:, :, :] = imVMEC_data['B_sub_theta'][:,:,:]  # negative sign for negative Jacobian

    # B_sub_psi
    B_sub_psi = file.createVariable("B_sub_psi", np.float64, ("radius", "fio_nphi", "fio_ntheta"))
    B_sub_psi.long_name = "B_sub_psi"
    B_sub_psi.units = "m"
    B_sub_psi[:, :, :] = imVMEC_data['B_sub_psi'][:,:,:]  # negative sign for negative Jacobian

    # B_sup_theta
    B_sup_theta = file.createVariable("B_sup_theta", np.float64, ("radius", "fio_nphi", "fio_ntheta"))
    B_sup_theta.long_name = "B_sup_theta"
    B_sup_theta.units = "m"
    B_sup_theta[:, :, :] = imVMEC_data['B_sup_theta'][:,:,:]  # negative sign for negative Jacobian
    
    # B_R
    B_R = file.createVariable("B_R", np.float64, ("radius", "fio_nphi", "fio_ntheta"))
    B_R.long_name = "B_R"
    B_R.units = "m"
    B_R[:, :, :] = imVMEC_data['B_R'][:,:,:]  # negative sign for negative Jacobian

    # B_Phi
    B_Phi = file.createVariable("B_Phi", np.float64, ("radius", "fio_nphi", "fio_ntheta"))
    B_Phi.long_name = "B_Phi"
    B_Phi.units = "m"
    B_Phi[:, :, :] = imVMEC_data['B_Phi'][:,:,:]  # negative sign for negative Jacobian

    # B_Z
    B_Z = file.createVariable("B_Z", np.float64, ("radius", "fio_nphi", "fio_ntheta"))
    B_Z.long_name = "B_Z"
    B_Z.units = "m"
    B_Z[:, :, :] = imVMEC_data['B_Z'][:,:,:]  # negative sign for negative Jacobian
    '''
    

def create_xm_xn_arrays(mpol, ntor, nfp, nyq=False):

    '''
    Calculates the arrays of the ordered Fourier modes used in VMEC
    * mpol (int): Number of poloidal Fourier modes
    * ntor (int): Number of toroidal Fourier modes (NOT multiplied by nfp)
    * nfp (int): Number of field periods
    * nyq (bool): Write arrays for Nyquist modes
    '''

    mnmax = (2*ntor+1)*mpol - ntor
    all_toroidal_modes = nfp*np.linspace(-ntor, ntor, 2*ntor+1)
    positive_toroidal_modes = nfp*np.linspace(0, ntor, ntor+1)
    for ipol in range(mpol):
        if (ipol == 0):
            xm = np.zeros(ntor+1)
            xn = positive_toroidal_modes
        else:
            xm = np.concatenate((xm,ipol*np.ones(2*ntor+1)))
            xn = np.concatenate((xn,all_toroidal_modes))

    if nyq: # Write out Nyquist modes
        mpol_nyq = mpol + 4
        ntor_nyq = ntor + 2 if ntor > 0 else 0
        mnmax_nyq = (2*ntor_nyq+1)*mpol_nyq - ntor_nyq   
        all_toroidal_modes_nyq = nfp*np.linspace(-ntor_nyq, ntor_nyq, 2*ntor_nyq+1)
        positive_toroidal_modes_nyq = nfp*np.linspace(0, ntor_nyq, ntor_nyq+1)
        for ipol in range(mpol_nyq):
            if (ipol == 0):
                xm_nyq = np.zeros(ntor_nyq+1)
                xn_nyq = positive_toroidal_modes_nyq
            else:
                xm_nyq = np.concatenate((xm_nyq,ipol*np.ones(2*ntor_nyq+1)))
                xn_nyq = np.concatenate((xn_nyq,all_toroidal_modes_nyq))

        return xm, xn, xm_nyq, xn_nyq
    
    else:
        return xm, xn

    
def fourier_transform_2d(X_array, theta_vals, phi_vals, xm, xn):

    '''
    Performs a 2D Fourier transform over theta,phi for the doubly periodic array, X_array
    Assumes equally spaced theta and phi coordinates. Should be replaced with FFT once normalizations in scipy are figured out
    '''

    xmnc = np.zeros(len(xm))
    xmns = np.zeros(len(xm))

    # Ensuring periodicity of the 2D array
    X_array = np.concatenate((X_array, np.reshape(X_array[0,:], (1, -1))), axis=0)
    X_array = np.concatenate((X_array, np.reshape(X_array[:,0], (-1, 1))), axis=1)
    X_array = np.transpose(X_array)
    
    phi, theta = np.meshgrid(phi_vals, theta_vals)
    
    norm_factor = 4*np.pi*np.pi
    for imn in range(len(xm)):
        m = float(xm[imn])
        n = float(xn[imn])
        if ( (int(m) == 0) and (int(n) == 0) ):
            xmnc_temp = (1./norm_factor) * np.trapz(np.trapz(X_array, x=theta_vals, axis=0), x=phi_vals, axis=0)
            xmns_temp = 0.0
        else:
            xmnc_temp = (2./norm_factor) * np.trapz(np.trapz(X_array * np.cos(m*theta - n*phi), x=theta_vals, axis=0), x=phi_vals, axis=0)
            xmns_temp = (2./norm_factor) * np.trapz(np.trapz(X_array * np.sin(m*theta - n*phi), x=theta_vals, axis=0), x=phi_vals, axis=0)

        xmnc[imn] = xmnc_temp
        xmns[imn] = xmns_temp

    return xmnc, xmns



def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def find_closest_indices(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, idx+1

        
