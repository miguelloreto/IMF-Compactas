#!/usr/bin/env python
##############################################################################
#
# This PPXF_POPULATION_EXAMPLE routine shows how to study stellar population with
# the procedure PPXF, which implements the Penalized Pixel-Fitting (pPXF) method by
# Cappellari M., & Emsellem E., 2004, PASP, 116, 138.
#
# MODIFICATION HISTORY:
#   V1.0.0: Adapted from PPXF_KINEMATICS_EXAMPLE.
#       Michele Cappellari, Oxford, 12 October 2011
#   V1.1.0: Made a separate routine for the construction of the templates
#       spectral library. MC, Vicenza, 11 October 2012
#   V1.1.1: Includes regul_error definition. MC, Oxford, 15 November 2012
#   V2.0.0: Translated from IDL into Python. MC, Oxford, 6 December 2013
#   V2.0.1: Fit SDSS rather than SAURON spectrum. MC, Oxford, 11 December 2013
#   V2.0.2: Support both Python 2.6/2.7 and Python 3.x. MC, Oxford, 25 May 2014
#   V2.0.3: Explicitly sort template files as glob() output may not be sorted.
#       Thanks to Marina Trevisan for reporting problems under Linux.
#       MC, Sydney, 4 February 2015
#   ********
#   V3.0.0: Adapted from PPXF_POPULATION_EXAMPLE_SDSS.
#
##############################################################################

from __future__ import print_function

import pyfits
import math
import csv
from scipy import ndimage
import numpy as np
import glob
import matplotlib.pyplot as plt
from time import clock

from ppxf import ppxf
import ppxf_util as util
from numpy import *

def setup_spectral_library(file_base, nAges, nMetal, FWHM_tem, lnorm, velscale, FWHM_gal, lamRange_gal, base_dir):

    # Read the list of filenames from the Single Stellar Population library
    # 

    file = file_base
    basefiles = loadtxt(file, unpack=True, usecols=[0], dtype=str)
    ages, mets, mstar = loadtxt(file, unpack=True, usecols=[1,2,3])

    # Extract the wavelength range and logarithmically rebin one spectrum
    # to the same velocity scale of the SAURON galaxy spectrum, to determine
    # the size needed for the array which will contain the template spectra.
    #
    hdu = pyfits.open(base_dir + '/' + basefiles[0])
    ssp = hdu[0].data
    h2 = hdu[0].header
    lamRange_stars = h2['CRVAL1'] + np.array([0.,h2['CDELT1']*(h2['NAXIS1']-1)])
    sspNew, logLam_stars, velscale_t = util.log_rebin(lamRange_stars, ssp, velscale=velscale)

    wave2 = np.exp(logLam_stars)

    #sspNew = sspNew[(wave2 >= lamRange_gal[0]) & (wave2 <= lamRange_gal[1])]
    #wave2 = wave2[(wave2 > lamRange_gal[0]) & (wave2 < lamRange_gal[1])]

    # Create a three dimensional array to store the
    # two dimensional grid of model spectra
    #
    templates = np.empty((sspNew.size, nAges, nMetal))
    L_M = np.empty((nAges, nMetal))

    # Convolve the whole Vazdekis library of spectral templates
    # with the quadratic difference between the SAURON and the
    # Vazdekis instrumental resolution. Logarithmically rebin
    # and store each template as a column in the array TEMPLATES.

    # Quadratic sigma difference in pixels Vazdekis --> SAURON
    # The formula below is rigorously valid if the shapes of the
    # instrumental spectral profiles are well approximated by Gaussians.
    #
    FWHM_dif = np.sqrt(FWHM_gal**2 - FWHM_tem**2)
    sigma = FWHM_dif/2.355/h2['CDELT1'] # Sigma difference in pixels

    # Here we make sure the spectra are sorted in both [M/H]
    # and Age along the two axes of the rectangular grid of templates.
    # A simple alphabetical ordering of Vazdekis's naming convention
    # does not sort the files by [M/H], so we do it explicitly below
    #
    #print(min(wave2, key=lambda x:abs(x-4020)))
    normf = min(wave2, key=lambda x:abs(x-lnorm)) 

    for k in range(0, nMetal):
      for j in range(0, nAges):
         hdu = pyfits.open(base_dir + basefiles[k * nAges + j])
         ssp = hdu[0].data
         ssp = ndimage.gaussian_filter1d(ssp,sigma)
         sspNew_temp, logLam_stars, velscale_t = util.log_rebin(lamRange_stars, ssp, velscale=velscale)
         sspNew = sspNew_temp#[(wave2 >= lamRange_gal[0]) & (wave2 <= lamRange_gal[1])]
         templates[:,j,k] = sspNew/sspNew_temp[wave2 == normf]
         L_M[j, k] = sspNew_temp[wave2 == normf]

    #lamRange_stars[0] = np.exp(min(logLam_stars))
    #lamRange_stars[1] = np.exp(max(logLam_stars))
    #logLam_stars = logLam_stars[(wave2 >= lamRange_gal[0]) & (wave2 <= lamRange_gal[1])]
    return templates, lamRange_stars, logLam_stars, ages, mets, mstar, L_M

#------------------------------------------------------------------------------

def format_vector(v, fmt="%9.5f"):
    return " ".join([fmt % x for x in v])

def ppxf_population_gas_example(file, file_in, lmin, ngas = 0, m_stars = 2, v_stars = 0.0, s_stars = 150., m_gas1 = 2, v_gas1 = 0., s_gas1 = 50., v_gas2 = 0., s_gas2 = 100.):

    z = 0.0           # OBSERVED SPECTRA ARE ASSUMED TO BE IN REST FRAME
    
    print('Running pPXF for ', file_in)
    # Read input parameters
    #
    f = open(file_in, "r")
    obs_dir = f.readline().split(' ')[0]
    base_dir = f.readline().split(' ')[0]
    out_dir = f.readline().split(' ')[0]
    file_base = f.readline().split(' ')[0]
    nAges     = int(f.readline().split(' ')[0])
    nMetal    = int(f.readline().split(' ')[0])
    FWHM_tem  = float(f.readline().split(' ')[0])
    lmin_0      = float(f.readline().split(' ')[0])
    lmax      = float(f.readline().split(' ')[0])
    lnorm     = float(f.readline().split(' ')[0])
    FWHM_gal  = float(f.readline().split(' ')[0])
    regul_err = float(f.readline().split(' ')[0])
    regul     = float(f.readline().split(' ')[0])
    ngas_0      = int(f.readline().split(' ')[0])

    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('Running pPXF for ', file)
    print('Number of gas components ', ngas)
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++')

    # Read galaxy spectrum 
    #
    with open(obs_dir+file) as f:
      reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
      first_row = next(reader)
      num_cols = len(first_row)
      if num_cols == 1:
        reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)

    if num_cols == 3:
      wave, galaxy, noise = loadtxt(obs_dir+file, unpack=True)
    else:
      wave, galaxy, noise, flag = loadtxt(obs_dir+file, unpack=True)

    # Only use the wavelength range in common between galaxy and stellar library.
    #
    lamRange_in = [lmin, lmax]   # MILES 3540 - 7409 A  
    mask = (wave >= lamRange_in[0]) & (wave <= lamRange_in[1])
    galaxy = galaxy[mask]  #/np.median(galaxy[mask])  # Normalize spectrum to avoid numerical issues
    noise = noise[mask]
    wave = wave[mask]
    
    lamRange_gal = lamRange_in
    #print(lamRange_in)
    lamRange_gal[0] = min(wave)
    lamRange_gal[1] = max(wave)
    #print(lamRange_gal)

    temp = noise[noise < 0]
    if temp.size > 0:
      print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      print('Negative noise values! (', temp.size, ' pixels)')
      noise[noise < 0] = -1 * noise[noise < 0]
      print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')

    # print(min(wave), max(wave))
    # noise = galaxy*0 + 0.01528           # Assume constant noise per pixel here

    # The velocity step was already chosen by the SDSS pipeline
    # and we convert it below to km/s
    #
    c = 299792.458 # speed of light in km/s
    velscale = np.log(wave[1]/wave[0])*c

    galaxyNew, logLam_gal, velscale = util.log_rebin(lamRange_gal, galaxy, velscale=velscale)
    noiseNew, logLam_gal, velscale = util.log_rebin(lamRange_gal, noise, velscale=velscale)

    stars_templates, lamRange_stars, logLam_stars, ages, mets, mstar, L_M = setup_spectral_library(file_base, 
             nAges, nMetal, FWHM_tem, lnorm, velscale, FWHM_gal, lamRange_gal, base_dir)

    waveNew = np.exp(logLam_gal)
    
    #print(min(waveNew), max(waveNew))  
    #print(np.exp(min(logLam_stars)), np.exp(max(logLam_stars)))   
    #print(lamRange_stars) 

    npix = galaxyNew.shape
    npmax = npix[0]-10
    galaxyNew = galaxyNew[0:npmax]
    noiseNew = noiseNew[0:npmax]
    waveNew = waveNew[0:npmax]

    norm_gal = float(galaxyNew[waveNew == min(waveNew, key=lambda x:abs(x-lnorm))])

    # The galaxy and the template spectra do not have the same starting wavelength.
    # For this reason an extra velocity shift DV has to be applied to the template
    # to fit the galaxy spectrum. We remove this artificial shift by using the
    # keyword VSYST in the call to PPXF below, so that all velocities are
    # measured with respect to DV. This assume the redshift is negligible.
    # In the case of a high-redshift galaxy one should de-redshift its
    # wavelength to the rest frame before using the line below as described
    # in PPXF_KINEMATICS_EXAMPLE_SAURON.
    #
    c = 299792.458
    #dv = (np.log(lamRange_stars[0])-np.log(waveNew[0]))*c # km/s
    #print('*****************************************************************')
    #print(dv)
    dv = (logLam_stars[0]-logLam_gal[0])*c # km/s
    #print(dv)
    #dv = (np.log(lamRange_stars[0])-np.log(lamRange_gal[0]))*c # km/s
    #print(dv)
    vel = c*z            # Initial estimate of the galaxy velocity in km/s
    goodpixels = util.determine_goodpixels(np.log(waveNew),lamRange_gal,vel)

    # Here the actual fit starts. The best fit is plotted on the screen.
    #
    # IMPORTANT: Ideally one would like not to use any polynomial in the fit
    # as the continuum shape contains important information on the population.
    # Unfortunately this is often not feasible, due to small calibration
    # uncertainties in the spectral shape. To avoid affecting the line strength of
    # the spectral features, we exclude additive polynomials (DEGREE=-1) and only use
    # multiplicative ones (MDEGREE=10). This is only recommended for population, not
    # for kinematic extraction, where additive polynomials are always recommended.
    #
    start_stars = [v_stars, s_stars]
    if ngas > 0:
      start_gas1 = [v_gas1, s_gas1] # (km/s), starting guess for [V,sigma]
    if ngas == 2:
      start_gas2 = [v_gas2, s_gas2] # (km/s), starting guess for [V,sigma]

    #------------------- Setup templates -----------------------

    #stars_templates, lamRange_temp, logLam_temp = \
    #    setup_spectral_library(velscale, FWHM_gal)

    # The stellar templates are reshaped into a 2-dim array with each spectrum
    # as a column, however we save the original array dimensions, which are
    # needed to specify the regularization dimensions
    #
    reg_dim = stars_templates.shape[1:]
    stars_templates = stars_templates.reshape(stars_templates.shape[0],-1)

    # See the pPXF documentation for the keyword REGUL,
    # for an explanation of the following two lines
    #
    #norm_stars = np.median(stars_templates)
    #stars_templates /= np.median(stars_templates) # Normalizes stellar templates by a scalar
    #regul_err = 0.004 # Desired regularization error

    # Construct a set of Gaussian emission line templates.
    # Estimate the wavelength fitted range in the rest frame.
    #
    if ngas > 0:
      #lamRange_gal = np.array([np.min(waveNew), np.max(waveNew)])/(1 + z)
      gas_templates, line_names, line_wave = \
          util.emission_lines(logLam_stars, lamRange_gal, FWHM_gal)
      nLines = gas_templates.shape[1]



    # Combines the stellar and gaseous templates into a single array
    # during the PPXF fit they will be assigned a different kinematic
    # COMPONENT value
    #
    if ngas == 0:
      templates = stars_templates
    if ngas == 1:
      templates = np.hstack([stars_templates, gas_templates])
    if ngas == 2:
      templates = np.hstack([stars_templates, gas_templates, gas_templates])

    #print(line_names)
    #print(lamRange_stars)
    #print(lamRange_gal)
    #print(logLam_stars.shape)
    #print(stars_templates.shape)
    #print(gas_templates.shape)

    # Assign component=0 to the stellar templates and
    # component=1 to the gas emission lines templates.
    # One can easily assign different kinematic components to different gas species
    # e.g. component=1 for the Balmer series, component=2 for the [OIII] doublet, ...)
    # Input a negative MOMENTS value to keep fixed the LOSVD of a component.
    #
    nTemps = stars_templates.shape[1]

    if ngas == 0:
      component = [0]*nTemps
      moments = m_stars # 4 --> fit (V,sig,h3,h4) for the stars and 2 (V,sig)
      start = start_stars
    if ngas == 1:
      component = [0]*nTemps + [1]*nLines 
      moments = [m_stars, 2] # [4, 2] --> fit (V,sig,h3,h4) for the stars and (V,sig) for the gas
      start = [start_stars, start_gas1] # adopt the same starting value for both gas and stars
    if ngas == 2:
      component = [0]*nTemps + [1]*nLines + [2]*nLines
      moments = [m_stars, m_gas1, 2] # [4, 2, 2] --> fit (V,sig,h3,h4) for the stars and (V,sig) for the gas
      start = [start_stars, start_gas1, start_gas2] # adopt the same starting value for both gas and stars

    t = clock()

    plt.clf()
    plt.subplot(211)

    if ngas == 0:
      pp = ppxf(templates, galaxyNew, noiseNew, velscale, start, 
              plot=True, moments=moments, degree=-1, goodpixels = goodpixels,
              vsyst=dv, clean=True, mdegree=-1, regul=regul, reg_dim=reg_dim,
              lam=waveNew, component=component, reddening=True)
    else:
      pp = ppxf(templates, galaxyNew, noiseNew, velscale, start, 
              plot=True, moments=moments, degree=-1,
              vsyst=dv, clean=False, mdegree=-1, regul=regul, reg_dim=reg_dim,
              lam=waveNew, component=component, reddening=True)

    #print(pp.weights)
    # When the two numbers below are the same, the solution is the smoothest
    # consistent with the observed spectrum.
    #
    #print('Desired Delta Chi^2: %.4g' % np.sqrt(2*goodpixels.size))
    #print('Current Delta Chi^2: %.4g' % ((pp.chi2 - 1)*goodpixels.size))
    print('Elapsed time in PPXF: %.2f s' % (clock() - t))

    if ngas > 0:
      w = np.where(np.array(component) == 1)[0] # Extract weights of gas emissions
      print('++++++++++++++++++++++++++++++')
      print('Gas V= %.4g and sigma= %.2g km/s' % (pp.sol[1][0], pp.sol[1][1]))
      print('Emission lines peak intensity:')
      for name, weight, line in zip(line_names, pp.weights[w], pp.matrix[:,w].T):
        print('%12s: %.3g' % (name, weight*np.max(line)))
      print('------------------------------')
     
    if ngas == 2:
      w = np.where(np.array(component) == 2)[0] # Extract weights of gas emissions
      print('++++++++++++++++++++++++++++++')
      print('Gas V= %.4g and sigma= %.2g km/s' % (pp.sol[2][0], pp.sol[2][1]))
      print('Emission lines peak intensity:')
      for name, weight, line in zip(line_names, pp.weights[w], pp.matrix[:,w].T):
        print('%12s: %.3g' % (name, weight*np.max(line)))
      print('------------------------------')

    plt.subplot(212)
    #s = reg_dim
    #print(reg_dim)
    #weights = pp.weights.reshape(s[1],s[2])/pp.weights.sum()
    weights = pp.weights[:np.prod(reg_dim)].reshape(reg_dim) #/pp.weights.sum()
    weights /= weights.sum()
    print(weights.shape)
    plt.imshow(np.rot90(weights * 100), interpolation='nearest', cmap='Reds', aspect='auto',
               extent=(np.log10(ages[0]), np.log10(ages[nAges*nMetal-1]), 
                       np.log10(mets[0]/0.019), np.log10(mets[nAges*nMetal-1]/0.019)))
    plt.colorbar()
    plt.title("Light Fraction (%)")
    plt.xlabel("log$_{10}$ Age (yr)")
    plt.ylabel("[M/H]")
    plt.tight_layout()

    sfh_name = out_dir + file + '_' + str(int(ngas)) + 'comp.png'
    plt.savefig(sfh_name)

    # plt.show()   

    file = open(out_dir + file + '_' + str(int(ngas)) + 'comp.out', "w")

    file.write("Best Fit:       V     sigma        h3        h4        h5        h6\n")
    file.write("comp. 0")
    if ngas > 0:
      for j in pp.sol[0]:
        file.write("{:10.3g} ".format(j))
      file.write("\ncomp. 1")
      for j in pp.sol[1]:
        file.write("{:10.3g} ".format(j))
      if ngas == 2:
        file.write("\ncomp. 2")
        for j in pp.sol[2]:
          file.write("{:10.3g} ".format(j))
      if ngas == 1:
        file.write("\ncomp. 2      0.0       0.0")
    else:
      for j in pp.sol:
        file.write("{:10.3g} ".format(j))
      file.write("\ncomp. 1      0.0       0.0")
      file.write("\ncomp. 2      0.0       0.0")
 
    file.write("\nchi2/DOF: %.4g" % pp.chi2)
    #file.write('Function evaluations:', ncalls)

    file.write('\nReddening E(B-V): {:.3f}'.format(pp.reddening))
    nw = nAges * nMetal
    file.write('\nNonzero Templates: {:.0f} / {:.0f}'.format(sum(weights > 0), nw))

    #file.write('\nDesired Delta Chi^2: %.4g' % np.sqrt(2*goodpixels.size))
    #file.write('\nCurrent Delta Chi^2: %.4g' % ((pp.chi2 - 1)*goodpixels.size))
    file.write('\nElapsed time in PPXF: %.2f s' % (clock() - t))

    # Mcor, Mini
    # ------------------------------------- 
    Mini = 0.0
    Mcor = 0.0

    for i in range(0,nMetal):
      for j in range(0, nAges):
        jj = i * nAges + j
        x_j = weights[j, i]
        L_M_j = L_M[j, i]
        Mini = Mini + x_j/L_M_j 
        Mcor = Mcor + x_j/L_M_j * mstar[jj]

    print('Norm factor: ', norm_gal)
    print('Mini: ', (Mini * norm_gal))
    print('Mcor: ', (Mcor * norm_gal))

    file.write('\nSpec Norm factor: {:.4f}'.format(norm_gal))
    file.write('\nSpec Norm lambda: {:.2f}'.format(lnorm))
    file.write('\nMini: {:.4f}'.format(Mini * norm_gal))
    file.write('\nMcor: {:.4f}'.format(Mcor * norm_gal))

    # Emission lines
    # ------------------------------------- 
    if ngas > 0:
      w = np.where(np.array(component) == 1)[0] # Extract weights of gas emissions
      file.write('\n++++++++++++++++++++++++++++++')
      file.write('\nGas V= {:.2f} and sigma= {:.2f} km/s'.format(pp.sol[1][0], pp.sol[1][1]))
      file.write('\nNumber of emission lines: {:.0f}'.format(line_names.size))
      file.write('\nEmission lines peak intensity:')
      for name, weight, line in zip(line_names, pp.weights[w], pp.matrix[:,w].T):
        file.write('\n{:s} {:.2e}'.format(name, weight*np.max(line)))
      file.write('\n------------------------------')

    if ngas == 2:
      w = np.where(np.array(component) == 2)[0] # Extract weights of gas emissions
      file.write('\n++++++++++++++++++++++++++++++')
      file.write('\nGas V= {:.2f} and sigma= {:.2f} km/s'.format(pp.sol[2][0], pp.sol[2][1]))
      file.write('\nNumber of emission lines: {:.0f}'.format(line_names.size))
      file.write('\nEmission lines peak intensity:')
      for name, weight, line in zip(line_names, pp.weights[w], pp.matrix[:,w].T):
        file.write('\n{:s} {:.2e}'.format(name, weight*np.max(line)))
      file.write('\n------------------------------')


    # Weights
    # ------------------------------------- 

    file.write("\n\n\n {:>5s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s} {:>15s}".format('# j', 
                'x_j(%)', 'Mini(%)', 'Mcor(%)', 'age_j', 'met_j', '(L/M)_j', 'Mstar'))

    for i in range(0,nMetal):
      for j in range(0, nAges):
        jj = i * nAges + j
        x_j = weights[j, i] 
        L_M_j = L_M[j, i]
        Mini_j = (x_j/L_M_j) / Mini
        Mcor_j = (x_j/L_M_j) * mstar[jj] / Mcor

        file.write("\n {:5.0f} {:15.5f} {:15.5f} {:15.5f} {:15.5e} {:15.5f} {:15.5e} {:15.5f}".format(jj+1, 
                   x_j*100, Mini_j*100, Mcor_j*100, ages[jj], mets[jj], L_M_j, mstar[jj]))


    #pp.weights.reshape(s[1],s[2])  #/pp.weights.sum()
    weights = pp.weights[:np.prod(reg_dim)].reshape(reg_dim) #/pp.weights.sum()
    weights /= weights.sum()


#    FWHM_gal = pp.sol[1]/c * 6000 * 2.355
#    templates, lamRange_temp, age, met, mstar, L_M = setup_spectral_library(file_base, 
#           nAges, nMetal, FWHM_tem, lnorm, velscale, FWHM_gal, lamRange, base_dir)
#    templates /= median(templates) 

    
#    spec = templates[:, 1, 1]
#    spec[:] = 0.0
#    for a in range(0,(s[1])):
#      for m in range(0,(s[2])):
#        spec = spec + templates[:, a, m] * weights[a, m]


    file.write("\n\n\n  {:>12s} {:>10s} {:>10s} {:>10s}".format('# lambda', 'obs', 'bestfit', 'noise'))
    s = galaxyNew.shape
    for i in range(1,(s[0])):
       file.write("\n  {:12.5f} {:10.5f} {:10.5f} {:10.5f}".format(waveNew[i], galaxyNew[i], pp.bestfit[i], noiseNew[i]))

    file.close()

#------------------------------------------------------------------------------

#if __name__ == '__main__':
    #ppxf_population_gas_example('0710-52203-0351.cxt', 'grid_in.dat', 3704, 1)

