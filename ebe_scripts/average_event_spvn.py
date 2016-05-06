#! /usr/bin/env python
"""
     This script performs event averaging for particle 
     spectra and anisotropic flow coefficients calculated 
     from event-by-event simulations

     v_n is analyzed up to n = 6

     Format for particle_XXX_vndata.dat file:
     n_order  real_part  real_part_err  imag_part  imag_part_err
     
     Format for particle_XXX_vndata_diff.dat file:
     pT(GeV)  pT_err(GeV)  dN/(2pi dy pT dpT)(GeV^-2)  dN/(2pi dy pT dpT)_err(GeV^-2)
     vn_real  vn_real_err  vn_imag  vn_imag_err

     All the errors are only statistic errors
"""

from sys import argv, exit
from os import path, mkdir
from glob import glob
from numpy import *
import shutil

# define colors
purple = "\033[95m"
green = "\033[92m"
blue = "\033[94m"
yellow = "\033[93m"
red = "\033[91m"
normal = "\033[0m"

try:
    working_folder = path.abspath(argv[1])
    avg_folder = path.join(path.abspath(argv[2]), working_folder.split('/')[-1])
    if(path.isdir(avg_folder)):
        shutil.rmtree(avg_folder)
    mkdir(avg_folder)
except IndexError:
    print("Usage: average_event_spvn.py working_folder results_folder")
    exit(1)

particle_list = ['211', '-211', '321', '-321', '2212', '-2212', '9999']
particle_name_list = ['pion_p', 'pion_m', 'kaon_p', 'kaon_m', 'proton', 
                      'anti_proton', 'charged_hadron']

n_order = 6

def calcualte_inte_vn(pT_low, pT_high, data):
    npT = 50
    pT_inte_array = linspace(pT_low, pT_high, npT)
    dN_event = data[:, 2]
    pT_event = data[:, 0]
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event+1e-30)))
    temp_vn_array = []
    for iorder in range(1, n_order):
        vn_real_event = temp_data[:, 4*iorder]
        vn_imag_event = temp_data[:, 4*iorder+2]
        vn_real_interp = interp(pT_inte_array, pT_event, vn_real_event)
        vn_imag_interp = interp(pT_inte_array, pT_event, vn_imag_event)
        vn_real_inte = (
            sum(vn_real_interp*dN_interp*pT_inte_array)
            /sum(dN_interp*pT_inte_array))
        vn_imag_inte = (
            sum(vn_imag_interp*dN_interp*pT_inte_array)
            /sum(dN_interp*pT_inte_array))
        vn_inte = vn_real_inte + 1j*vn_imag_inte
        temp_vn_array.append(vn_inte)
    return(temp_vn_array)


def calcualte_diff_vn_SP(pT_ref_low, pT_ref_high, data):
    npT = 50
    pT_inte_array = linspace(pT_ref_low, pT_ref_high, npT)
    dN_event = data[:, 2]
    pT_event = data[:, 0]
    dN_interp = exp(interp(pT_inte_array, pT_event, log(dN_event+1e-30)))
    dN_ref = sum(dN_interp*pT_inte_array)
    temp_vn_real_array = []
    temp_vn_imag_array = []
    temp_vn_denorm_array = []
    for iorder in range(1, n_order):
        vn_real_event = temp_data[:, 4*iorder]
        vn_imag_event = temp_data[:, 4*iorder+2]
        vn_real_interp = interp(pT_inte_array, pT_event, vn_real_event)
        vn_imag_interp = interp(pT_inte_array, pT_event, vn_imag_event)
        vn_real_inte = (
            sum(vn_real_interp*dN_interp*pT_inte_array)
            /sum(dN_interp*pT_inte_array))
        vn_imag_inte = (
            sum(vn_imag_interp*dN_interp*pT_inte_array)
            /sum(dN_interp*pT_inte_array))
        vn_ref = vn_real_inte + 1j*vn_imag_inte
        vn_pt = vn_real_event + 1j*vn_imag_event
        numerator_real = real(dN_event*vn_pt*dN_ref*conj(vn_ref))
        numerator_imag = imag(dN_event*vn_pt*dN_ref*conj(vn_ref))
        denorm = dN_event*dN_ref
        temp_vn_real_array.append(numerator_real)
        temp_vn_imag_array.append(numerator_imag)
    temp_vn_denorm_array.append(denorm)
    return(temp_vn_real_array, temp_vn_imag_array, temp_vn_denorm_array)


file_folder_list = glob(path.join(working_folder, '*'))
nev = len(file_folder_list)
for ipart, particle_id in enumerate(particle_list):
    print("processing %s ..." % particle_name_list[ipart])
    
    # first particle yield dN/dy
    file_name = 'particle_%s_vndata.dat' % particle_id

    dN_dy = []
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        temp_data = loadtxt(path.join(results_folder, file_name))

        dN_dy.append(temp_data[0, 1])

    dN_dy = array(dN_dy)
    dN_dy_avg = mean(dN_dy)
    dN_dy_avg_err = std(dN_dy)/sqrt(nev)

    # then <pT>, vn, dN/(2pi dy pT dpT), vn{SP}(pT)
    file_name = 'particle_%s_vndata_diff.dat' % particle_id
   
    pT_array = []
    dN_array = []
    vn_phenix_array = []
    vn_star_array = []
    vn_alice_array = []
    vn_cms_array = []
    vn_atlas_array = []
    vn_diff_phenix_real = []; vn_diff_phenix_imag = []; vn_diff_phenix_denorm = []
    vn_diff_star_real = []; vn_diff_star_imag = []; vn_diff_star_denorm = []
    vn_diff_alice_real = []; vn_diff_alice_imag = []; vn_diff_alice_denorm = []
    vn_diff_cms_real = []; vn_diff_cms_imag = []; vn_diff_cms_denorm = []
    vn_diff_atlas_real = []; vn_diff_atlas_imag = []; vn_diff_atlas_denorm = []
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        temp_data = loadtxt(path.join(results_folder, file_name))
        
        dN_event = temp_data[:, 2]  # dN/(2pi dy pT dpT)
        pT_event = temp_data[:, 0]

        # record particle spectra
        pT_array.append(pT_event)
        dN_array.append(dN_event)

        # pT-integrated vn
        # vn with PHENIX pT cut
        temp_vn_array = calcualte_inte_vn(0.2, 2.0, temp_data)
        vn_phenix_array.append(temp_vn_array)

        # vn with STAR pT cut
        temp_vn_array = calcualte_inte_vn(0.15, 2.0, temp_data)
        vn_star_array.append(temp_vn_array)

        # vn with ALICE pT cut
        temp_vn_array = calcualte_inte_vn(0.2, 3.0, temp_data)
        vn_alice_array.append(temp_vn_array)
        
        # vn with CMS pT cut
        temp_vn_array = calcualte_inte_vn(0.3, 3.0, temp_data)
        vn_cms_array.append(temp_vn_array)
        
        # vn with ATLAS pT cut
        temp_vn_array = calcualte_inte_vn(0.5, 3.0, temp_data)
        vn_atlas_array.append(temp_vn_array)

        # pT-differential vn using scalar-product method
        # vn{SP}(pT) with PHENIX pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff = (
                                     calcualte_diff_vn_SP(0.15, 2.0, temp_data))
        vn_diff_phenix_real.append(temp_vn_diff_real);
        vn_diff_phenix_imag.append(temp_vn_diff_imag);
        vn_diff_phenix_denorm.append(temp_dn_diff);

        # vn{SP}(pT) with STAR pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff = (
                                     calcualte_diff_vn_SP(0.15, 2.0, temp_data))
        vn_diff_star_real.append(temp_vn_diff_real);
        vn_diff_star_imag.append(temp_vn_diff_imag);
        vn_diff_star_denorm.append(temp_dn_diff);
        
        # vn{SP}(pT) with ALICE pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff = (
                                     calcualte_diff_vn_SP(0.2, 3.0, temp_data))
        vn_diff_alice_real.append(temp_vn_diff_real);
        vn_diff_alice_imag.append(temp_vn_diff_imag);
        vn_diff_alice_denorm.append(temp_dn_diff);
        
        # vn{SP}(pT) with CMS pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff = (
                                     calcualte_diff_vn_SP(0.3, 3.0, temp_data))
        vn_diff_cms_real.append(temp_vn_diff_real);
        vn_diff_cms_imag.append(temp_vn_diff_imag);
        vn_diff_cms_denorm.append(temp_dn_diff);
        
        # vn{SP}(pT) with ATLAS pT cut
        temp_vn_diff_real, temp_vn_diff_imag, temp_dn_diff = (
                                     calcualte_diff_vn_SP(0.5, 3.0, temp_data))
        vn_diff_atlas_real.append(temp_vn_diff_real);
        vn_diff_atlas_imag.append(temp_vn_diff_imag);
        vn_diff_atlas_denorm.append(temp_dn_diff);

    # now we perform event average
    dN_array = array(dN_array)
    pT_array = array(pT_array)
    n_pT = len(pT_array[0, :])
    pT_spectra = zeros([n_pT])
    for ipT in range(len(pT_array[0, :])):
        dN_temp = sum(dN_array[:, ipT]*pT_array[:, ipT])
        if(dN_temp > 0):
            pT_spectra[ipT] = sum(pT_array[:, ipT]**2.*dN_array[:, ipT])/dN_temp
        else:
            pT_spectra[ipT] = mean(pT_array[:, ipT])
    dN_spectra = mean(pT_array*dN_array, 0)/pT_spectra   # dN/(2pi dy pT dpT)
    dN_spectra_err = std(pT_array*dN_array, 0)/pT_spectra/sqrt(nev)

    # calculate mean pT
    pT_interp = linspace(0.05, 2.95, 30)
    dN_interp = exp(interp(pT_interp, pT_spectra, log(dN_spectra+1e-30)))
    dN_interp_err = interp(pT_interp, pT_spectra, dN_spectra_err)
    mean_pT = sum(pT_interp**2.*dN_interp)/sum(pT_interp*dN_interp)
    mean_pT_upper = (sum(pT_interp**2.*(dN_interp+dN_interp_err))
                     /sum(pT_interp*(dN_interp+dN_interp_err)))
    mean_pT_lower = (sum(pT_interp**2.*(dN_interp-dN_interp_err))
                     /sum(pT_interp*(dN_interp-dN_interp_err)))
    mean_pT_err = max(abs(mean_pT_upper - mean_pT), 
                      abs(mean_pT - mean_pT_lower))

    # calcualte vn{2}
    vn_phenix_array = array(vn_phenix_array)
    vn_phenix_2 = sqrt(mean(abs(vn_phenix_array)**2., 0)) + 1e-30
    vn_phenix_2_err = std(abs(vn_phenix_array)**2., 0)/sqrt(nev)/2./vn_phenix_2

    vn_star_array = array(vn_star_array)
    vn_star_2 = sqrt(mean(abs(vn_star_array)**2., 0)) + 1e-30
    vn_star_2_err = std(abs(vn_star_array)**2., 0)/sqrt(nev)/2./vn_star_2

    vn_alice_array = array(vn_alice_array)
    vn_alice_2 = sqrt(mean(abs(vn_alice_array)**2., 0)) + 1e-30
    vn_alice_2_err = std(abs(vn_alice_array)**2., 0)/sqrt(nev)/2./vn_alice_2
    
    vn_cms_array = array(vn_cms_array)
    vn_cms_2 = sqrt(mean(abs(vn_cms_array)**2., 0)) + 1e-30
    vn_cms_2_err = std(abs(vn_cms_array)**2., 0)/sqrt(nev)/2./vn_cms_2
    
    vn_atlas_array = array(vn_atlas_array)
    vn_atlas_2 = sqrt(mean(abs(vn_atlas_array)**2., 0)) + 1e-30
    vn_atlas_2_err = std(abs(vn_atlas_array)**2., 0)/sqrt(nev)/2./vn_atlas_2

    # calcualte vn{SP}(pT)
    vn_diff_phenix_real = array(vn_diff_phenix_real)
    vn_diff_phenix_imag = array(vn_diff_phenix_imag)
    vn_diff_phenix_denorm = array(vn_diff_phenix_denorm) + 1e-30
    vn_denorm = vn_phenix_2.reshape(len(vn_phenix_2), 1)
    vn_denorm_err = vn_phenix_2_err.reshape(len(vn_phenix_2_err), 1)
    vn_diff_SP_phenix = (
        mean(vn_diff_phenix_real, 0)/mean(vn_diff_phenix_denorm, 0)/vn_denorm)
    vn_diff_SP_phenix_err = sqrt(
        ( std(vn_diff_phenix_real, 0)/sqrt(nev)/mean(vn_diff_phenix_denorm, 0)
          /vn_denorm)**2.
        + (vn_diff_SP_phenix*vn_denorm_err/vn_denorm)**2.)
        
    vn_diff_star_real = array(vn_diff_star_real)
    vn_diff_star_imag = array(vn_diff_star_imag)
    vn_diff_star_denorm = array(vn_diff_star_denorm) + 1e-3
    vn_denorm = vn_star_2.reshape(len(vn_star_2), 1)
    vn_denorm_err = vn_star_2_err.reshape(len(vn_star_2_err), 1)
    vn_diff_SP_star = (
        mean(vn_diff_star_real, 0)/mean(vn_diff_star_denorm, 0)/vn_denorm)
    vn_diff_SP_star_err = sqrt(
        ( std(vn_diff_star_real, 0)/sqrt(nev)/mean(vn_diff_star_denorm, 0)
          /vn_denorm)**2.
        + (vn_diff_SP_star*vn_denorm_err/vn_denorm)**2.)
    
    vn_diff_alice_real = array(vn_diff_alice_real)
    vn_diff_alice_imag = array(vn_diff_alice_imag)
    vn_diff_alice_denorm = array(vn_diff_alice_denorm) + 1e-30
    vn_denorm = vn_alice_2.reshape(len(vn_alice_2), 1)
    vn_denorm_err = vn_alice_2_err.reshape(len(vn_alice_2_err), 1)
    vn_diff_SP_alice = (
        mean(vn_diff_alice_real, 0)/mean(vn_diff_alice_denorm, 0)/vn_denorm)
    vn_diff_SP_alice_err = sqrt(
        ( std(vn_diff_alice_real, 0)/sqrt(nev)/mean(vn_diff_alice_denorm, 0)
          /vn_denorm)**2.
        + (vn_diff_SP_alice*vn_denorm_err/vn_denorm)**2.)
    
    vn_diff_cms_real = array(vn_diff_cms_real)
    vn_diff_cms_imag = array(vn_diff_cms_imag)
    vn_diff_cms_denorm = array(vn_diff_cms_denorm) + 1e-30
    vn_denorm = vn_cms_2.reshape(len(vn_cms_2), 1)
    vn_denorm_err = vn_cms_2_err.reshape(len(vn_cms_2_err), 1)
    vn_diff_SP_cms = (
        mean(vn_diff_cms_real, 0)/mean(vn_diff_cms_denorm, 0)/vn_denorm)
    vn_diff_SP_cms_err = sqrt(
        ( std(vn_diff_cms_real, 0)/sqrt(nev)/mean(vn_diff_cms_denorm, 0)
          /vn_denorm)**2.
        + (vn_diff_SP_cms*vn_denorm_err/vn_denorm)**2.)
    
    vn_diff_atlas_real = array(vn_diff_atlas_real)
    vn_diff_atlas_imag = array(vn_diff_atlas_imag)
    vn_diff_atlas_denorm = array(vn_diff_atlas_denorm) + 1e-30
    vn_denorm = vn_atlas_2.reshape(len(vn_atlas_2), 1)
    vn_denorm_err = vn_atlas_2_err.reshape(len(vn_atlas_2_err), 1)
    vn_diff_SP_atlas = (
        mean(vn_diff_atlas_real, 0)/mean(vn_diff_atlas_denorm, 0)/vn_denorm)
    vn_diff_SP_atlas_err = sqrt(
        ( std(vn_diff_atlas_real, 0)/sqrt(nev)/mean(vn_diff_atlas_denorm, 0)
          /vn_denorm)**2.
        + (vn_diff_SP_atlas*vn_denorm_err/vn_denorm)**2.)
    
    # then particle rapidity distribution
    if particle_id == '9999':
        file_name = 'particle_%s_dNdeta.dat' % particle_id
    else:
        file_name = 'particle_%s_dNdy.dat' % particle_id

    eta_array = []
    dN_array = []
    vn_array = []
    for ifolder in range(nev):
        results_folder = path.abspath(file_folder_list[ifolder])
        temp_data = loadtxt(path.join(results_folder, file_name))

        eta_array.append(temp_data[:, 0])
        dN_array.append(temp_data[:, 1])
        temp_vn_array = []
        for iorder in range(1, n_order):
            vn_real = temp_data[:, 6*iorder-3]
            vn_imag = temp_data[:, 6*iorder-1]
            vn = vn_real + 1j*vn_imag
            temp_vn_array.append(vn)
        vn_array.append(temp_vn_array)

    eta_array = array(eta_array)
    dN_array = array(dN_array)
    vn_array = array(vn_array)

    eta_point = mean(eta_array, 0)
    dNdeta = mean(dN_array, 0)
    dNdeta_err = std(dN_array, 0)/sqrt(nev)
    vn_eta = sqrt(mean(abs(vn_array)**2., 0))
    vn_eta_err = std(abs(vn_array)**2., 0)/sqrt(nev)/2./(vn_eta + 1e-15)
    vn_eta_real = mean(real(vn_array), 0)
    vn_eta_real_err = std(real(vn_array), 0)/sqrt(nev)
    
    # finally, output all the results
    output_filename = "%s_integrated_observables.dat" % particle_name_list[ipart]
    f = open(output_filename, 'w')
    f.write("dN/dy= %.10e +/- %.10e\n" % (dN_dy_avg, dN_dy_avg_err))
    f.write("<pT>= %.10e +/- %.10e\n" % (mean_pT, mean_pT_err))
    for iorder in range(1, n_order):
        f.write("v_%d{2}(phenix)= %.10e +/- %.10e\n" 
                % (iorder, vn_phenix_2[iorder-1], vn_phenix_2_err[iorder-1]))
        f.write("v_%d{2}(STAR)= %.10e +/- %.10e\n" 
                % (iorder, vn_star_2[iorder-1], vn_star_2_err[iorder-1]))
        f.write("v_%d{2}(ALICE)= %.10e +/- %.10e\n" 
                % (iorder, vn_alice_2[iorder-1], vn_alice_2_err[iorder-1]))
        f.write("v_%d{2}(CMS)= %.10e +/- %.10e\n" 
                % (iorder, vn_cms_2[iorder-1], vn_cms_2_err[iorder-1]))
        f.write("v_%d{2}(ATLAS)= %.10e +/- %.10e\n" 
                % (iorder, vn_atlas_2[iorder-1], vn_atlas_2_err[iorder-1]))
    f.close()
    shutil.move(output_filename, avg_folder)
    
    output_filename = ("%s_differential_observables_PHENIX.dat" 
                       % particle_name_list[ipart])
    f = open(output_filename, 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_phenix[iorder-1, ipT], 
                                        vn_diff_SP_phenix_err[iorder-1, ipT]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)

    output_filename = ("%s_differential_observables_STAR.dat" 
                       % particle_name_list[ipart])
    f = open(output_filename, 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_star[iorder-1, ipT], 
                                        vn_diff_SP_star_err[iorder-1, ipT]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)
    
    output_filename = ("%s_differential_observables_ALICE.dat" 
                       % particle_name_list[ipart])
    f = open(output_filename, 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_alice[iorder-1, ipT], 
                                        vn_diff_SP_alice_err[iorder-1, ipT]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)
    
    output_filename = ("%s_differential_observables_CMS.dat" 
                       % particle_name_list[ipart])
    f = open(output_filename, 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_cms[iorder-1, ipT], 
                                        vn_diff_SP_cms_err[iorder-1, ipT]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)
    
    output_filename = ("%s_differential_observables_ATLAS.dat" 
                       % particle_name_list[ipart])
    f = open(output_filename, 'w')
    f.write("#pT  dN/(2pi dy pT dpT)  dN/(2pi dy pT dpT)_err  vn{SP}  vn{SP}_err\n")
    for ipT in range(len(pT_spectra)):
        f.write("%.10e  %.10e  %.10e  "
                % (pT_spectra[ipT], dN_spectra[ipT], dN_spectra_err[ipT]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  " % (vn_diff_SP_atlas[iorder-1, ipT], 
                                        vn_diff_SP_atlas_err[iorder-1, ipT]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)

    output_filename = ("%s_rapidity_distribution.dat" 
                       % particle_name_list[ipart])
    f = open(output_filename, 'w')
    if(particle_id == '9999'):
        f.write("#eta  dN/deta  dN/deta_err  vn{2}(eta)  vn{2}(eta)_err\n")
    else:
        f.write("#y  dN/dy  dN/dy_err  vn{2}(y)  vn{2}(y)_err\n")
    for ieta in range(len(eta_point)):
        f.write("%.10e  %.10e  %.10e  "
                % (eta_point[ieta], dNdeta[ieta], dNdeta_err[ieta]))
        for iorder in range(1, n_order):
            f.write("%.10e  %.10e  %.10e  %.10e  "
                    % (vn_eta[iorder-1, ieta], vn_eta_err[iorder-1, ieta],
                       vn_eta_real[iorder-1, ieta],
                       vn_eta_real_err[iorder-1, ieta]))
        f.write("\n")
    f.close()
    shutil.move(output_filename, avg_folder)

print "Analysis is done."

