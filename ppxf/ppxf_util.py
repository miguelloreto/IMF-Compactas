#######################################################################
#
# Copyright (C) 2001-2014, Michele Cappellari
# E-mail: cappellari_at_astro.ox.ac.uk
#
# This software is provided as is without any warranty whatsoever.
# Permission to use, for non-commercial purposes is granted.
# Permission to modify for personal or internal use is granted,
# provided this copyright and disclaimer are included unchanged
# at the beginning of the file. All other rights are reserved.
#
#######################################################################
#
# NAME:
#   LOG_REBIN
#
# PURPOSE:
#   Logarithmically rebin a spectrum, while rigorously conserving the flux.
#   Basically the photons in the spectrum are simply ridistributed according
#   to a new grid of pixels, with non-uniform size in the spectral direction.
#
#   This routine makes the `standard' zero-order assumption that the spectrum
#   is *constant* within each pixels. It is possible to perform log-rebinning
#   by assuming the spectrum is represented by a piece-wise polynomial of
#   higer degree, while still obtaining a uniquely defined linear problem,
#   but this reduces to a deconvolution and amplifies noise.
#
#   This same routine can be used to compute approximate errors
#   of the log-rebinned spectrum. To do this type the command
#
#       LOG_REBIN, lamRange, err^2, err2New
#
#   and the desired errors will be given by SQRT(err2New).
#   NB: This rebinning of the error-spectrum is very *approximate* as
#   it does not consider the correlation introduced by the rebinning!
#
# CALLING SEQUENCE:
#   LOG_REBIN, lamRange, spec, specNew, logLam, $
#       OVERSAMPLE=oversample, VELSCALE=velScale, /FLUX
#
# INPUTS:
#   LAMRANGE: two elements vector containing the central wavelength
#       of the first and last pixels in the spectrum, which is assumed
#       to have constant wavelength scale! E.g. from the values in the
#       standard FITS keywords: LAMRANGE = CRVAL1 + [0,CDELT1*(NAXIS1-1)].
#       It must be LAMRANGE[0] < LAMRANGE[1].
#   SPEC: input spectrum.
#
# OUTPUTS:
#   SPECNEW: logarithmically rebinned spectrum.
#   LOGLAM: log(lambda) (*natural* logarithm: ALOG) of the central
#       wavelength of each pixel. This is the log of the geometric
#       mean of the borders of each pixel.
#
# KEYWORDS:
#   FLUX: Set this keyword to preserve total flux. In this case the
#       log rebinning changes the pixels flux in proportion to their
#       dLam so the following command will show large differences
#       beween the spectral shape before and after LOG_REBIN:
#
#           plot, exp(logLam), specNew  # Plot log-rebinned spectrum
#           oplot, range(lamRange[0],lamRange[1],n_elements(spec)), spec
#
#       By defaul, when this keyword is *not* set, the above two lines
#       produce two spectra that almost perfectly overlap each other.
#   OVERSAMPLE: Oversampling can be done, not to loose spectral resolution,
#       especally for extended wavelength ranges and to avoid aliasing.
#       Default: OVERSAMPLE=1 ==> Same number of output pixels as input.
#   VELSCALE: velocity scale in km/s per pixels. If this variable is
#       not defined, then it will contain in output the velocity scale.
#       If this variable is defined by the user it will be used
#       to set the output number of pixels and wavelength scale.
#
# MODIFICATION HISTORY:
#   V1.0.0: Using interpolation. Michele Cappellari, Leiden, 22 October 2001
#   V2.0.0: Analytic flux conservation. MC, Potsdam, 15 June 2003
#   V2.1.0: Allow a velocity scale to be specified by the user.
#       MC, Leiden, 2 August 2003
#   V2.2.0: Output the optional logarithmically spaced wavelength at the
#       geometric mean of the wavelength at the border of each pixel.
#       Thanks to Jesus Falcon-Barroso. MC, Leiden, 5 November 2003
#   V2.2.1: Verify that lamRange[0] < lamRange[1].
#       MC, Vicenza, 29 December 2004
#   V2.2.2: Modified the documentation after feedback from James Price.
#       MC, Oxford, 21 October 2010
#   V2.3.0: By default now preserve the shape of the spectrum, not the
#       total flux. This seems what most users expect from the procedure.
#       Set the keyword /FLUX to preserve flux like in previous version.
#       MC, Oxford, 30 November 2011
#   V3.0.0: Translated from IDL into Python. MC, Santiago, 23 November 2013
#   V3.1.0: Fully vectorized log_rebin. Typical speed up by two orders of magnitude.
#       MC, Oxford, 4 March 2014
#
#----------------------------------------------------------------------

from __future__ import print_function

import numpy as np

def log_rebin(lamRange, spec, oversample=False, velscale=None, flux=False):
    """
    Logarithmically rebin a spectrum, while rigorously conserving the flux.
    Basically the photons in the spectrum are simply ridistributed according
    to a new grid of pixels, with non-uniform size in the spectral direction.
    
    When the flux keyword is set, this program performs an exact integration 
    of the original spectrum, assumed to be a step function within the 
    linearly-spaced pixels, onto the new logarithmically-spaced pixels. 
    The output was tested to agree with the analystic solution.

    """
    lamRange = np.asarray(lamRange)
    if len(lamRange) != 2:
        raise ValueError('lamRange must contain two elements')
    if lamRange[0] >= lamRange[1]:
        raise ValueError('It must be lamRange[0] < lamRange[1]')
    s = spec.shape
    if len(s) != 1:
        raise ValueError('input spectrum must be a vector')
    n = s[0]
    if oversample:
        m = int(n*oversample)
    else:
        m = int(n)

    dLam = np.diff(lamRange)/(n - 1.)        # Assume constant dLam
    lim = lamRange/dLam + [-0.5, 0.5]        # All in units of dLam
    borders = np.linspace(*lim, num=n+1)     # Linearly
    logLim = np.log(lim)

    c = 299792.458                           # Speed of light in km/s
    if velscale is None:                     # Velocity scale is set by user
        velscale = np.diff(logLim)/m*c       # Only for output
    else:
        logScale = velscale/c
        m = int(np.diff(logLim)/logScale)    # Number of output pixels
        logLim[1] = logLim[0] + m*logScale

    newBorders = np.exp(np.linspace(*logLim, num=m+1)) # Logarithmically
    k = (newBorders - lim[0]).clip(0, n-1).astype(int)

    specNew = np.add.reduceat(spec, k)[:-1]  # Do analytic integral
    specNew *= np.diff(k) > 0    # fix for design flaw of reduceat()
    specNew += np.diff((newBorders - borders[k])*spec[k])

    if not flux:
        specNew /= np.diff(newBorders)

    # Output log(wavelength): log of geometric mean
    logLam = np.log(np.sqrt(newBorders[1:]*newBorders[:-1])*dLam)

    return specNew, logLam, velscale

#----------------------------------------------------------------------
#
# PPXF_DETERMINE_GOODPIXELS: Example routine to generate the vector of goodPixels
#     to be used as input keyword for the routine PPXF. This is useful to mask
#     gas emission lines or atmospheric absorptions.
#     It can be trivially adapted to mask different lines.
#
# INPUT PARAMETERS:
# - LOGLAM: Natural logarithm ALOG(wave) of the wavelength in Angstrom
#     of each pixel of the log rebinned *galaxy* spectrum.
# - LAMRANGETEMP: Two elements vectors [lamMin2,lamMax2] with the minimum and
#     maximum wavelength in Angstrom in the stellar *template* used in PPXF.
# - VEL: Estimate of the galaxy velocity in km/s.
#
# V1.0.0: Michele Cappellari, Leiden, 9 September 2005
# V1.0.1: Made a separate routine and included additional common emission lines.
#   MC, Oxford 12 January 2012
# V2.0.0: Translated from IDL into Python. MC, Oxford, 10 December 2013
# V2.0.1: Updated line list. MC, Oxford, 8 January 2014

def determine_goodpixels(logLam, lamRangeTemp, vel):
    """
    Generates a list of goodpixels to mask a given set of gas emission
    lines. This is meant to be used as input for PPXF.

    """
#                     -----[OII]-----    Hdelta   Hgamma   Hbeta   -----[OIII]-----   [OI]    -----[NII]-----   Halpha   -----[SII]-----
    lines = np.array([3726.03, 3728.82, 4101.76, 4340.47, 4861.33, 4958.92, 5006.84, 6300.30, 6548.03, 6583.41, 6562.80, 6716.47, 6730.85])
    dv = lines*0 + 800 # width/2 of masked gas emission region in km/s
    c = 299792.458 # speed of light in km/s

    flag = logLam < 0  # empy mask
    for j in range(lines.size):
        flag |= (logLam > np.log(lines[j]) + (vel - dv[j])/c) \
              & (logLam < np.log(lines[j]) + (vel + dv[j])/c)

    flag |= logLam < np.log(lamRangeTemp[0]) + (vel + 900.)/c  # Mask edges of
    flag |= logLam > np.log(lamRangeTemp[1]) + (vel - 900.)/c  # stellar library

    #flag |= (logLam > np.log(4500)) & (logLam < np.log(4700)) # Mask region between 4500 and 4700 A

    return np.where(flag == 0)[0]

#------------------------------------------------------------------------------
# V1.0.0: Michele Cappellari, Oxford, 7 January 2014
# V1.1.0: Fixes [OIII] and [NII] doublets to the theoretical flux ratio.
#       Returns line names together with emission lines templates.
#       MC, Oxford, 3 August 2014
# V1.1.1: Only returns lines included within the estimated fitted wavelength range.
#       This avoids identically zero gas templates being included in the PPXF fit
#       which can cause numearical instabilities in the solution of the system.
#       MC, Oxford, 3 September 2014

def emission_lines(logLam_temp, lamRange_gal, FWHM_gal):
    """
    Generates an array of Gaussian emission lines to be used as templates in PPXF.
    Additional lines can be easily added by editing this procedure.
    - logLam_temp is the natural log of the wavelength of the templates in Angstrom.
      logLam_temp should be the same as that of the stellar templates.
    - lamRange_gal is the estimated rest-frame fitted wavelength range
      Typically lamRange_gal = np.array([np.min(wave), np.max(wave)])/(1 + z),
      where wave is the observed wavelength of the fitted galaxy pixels and
      z is an initial very rough estimate of the galaxy redshift.
    - FWHM_gal is the instrumantal FWHM of the galaxy spectrum under study in
      Angstrom. Here it is assumed constant. It could be a function of wavelength.
    - The [OI], [OIII] and [NII] doublets are fixed at theoretical flux ratio~3.

    """
    lam = np.exp(logLam_temp)
    sigma = FWHM_gal/2.355 # Assumes instrumental sigma is constant in Angstrom

    # Balmer Series:      H12       H10       H10      H9       H8       Hepsilon Hdelta   Hgamma   Hbeta    Halpha
    line_wave = np.array([3750.158, 3770.637, 3797.90, 3835.39, 3889.06, 3970.08, 4101.76, 4340.47, 4861.33, 6562.82])
    line_names = np.array(['H12', 'H11', 'H10', 'H9', 'H8', 'Hepsilon', 'Hdelta', 'Hgamma', 'Hbeta', 'Halpha'])
    emission_lines = np.exp(-0.5*((lam[:,np.newaxis] - line_wave)/sigma)**2)

    #                 -----[OII]-----    -----[SII]-----
    lines = np.array([3726.03, 3728.82, 6716.44, 6730.81])
    names = np.array(['[OII]3726', '[OII]3729', '[SII]6716', '[SII]6731'])
    gauss = np.exp(-0.5*((lam[:,np.newaxis] - lines)/sigma)**2)
    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, lines)

    #                 -----[SII]----- 
    lines = np.array([4068.600])
    names = np.array(['[SII]4068'])
    gauss = np.exp(-0.5*((lam[:,np.newaxis] - lines)/sigma)**2)
    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, lines)

    #                 -----[SIII]----- 
    lines = np.array([6312.060])
    names = np.array(['[SIII]6312'])
    gauss = np.exp(-0.5*((lam[:,np.newaxis] - lines)/sigma)**2)
    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, lines)


    #                 -----[NIII]----- 
    lines = np.array([3868.76, 3967.47])
    names = np.array(['[NIII]3868', '[NIII]3967'])
    gauss = np.exp(-0.5*((lam[:,np.newaxis] - lines)/sigma)**2)
    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, lines)

    #                 -----[HeI]----- 
    lines = np.array([3888.647, 4026.190, 4143.761, 4471.479, 5875.62, 6678.15])
    names = np.array(['[HeI]3888', '[HeI]4026', '[HeI]4143', '[HeI]4471', '[HeI]5875', '[HeI]6678'])
    gauss = np.exp(-0.5*((lam[:,np.newaxis] - lines)/sigma)**2)
    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, lines)

    #                 -----[OIII]----- 
    lines = np.array([4363.21])
    names = np.array(['[OIII]4363'])
    gauss = np.exp(-0.5*((lam[:,np.newaxis] - lines)/sigma)**2)
    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, lines)

    #                 -----[ClIII]----- 
    lines = np.array([5517.709])
    names = np.array(['[ClIII]5517'])
    gauss = np.exp(-0.5*((lam[:,np.newaxis] - lines)/sigma)**2)
    emission_lines = np.column_stack([emission_lines, gauss])
    line_names = np.append(line_names, names)
    line_wave = np.append(line_wave, lines)

    #                 -----[OIII]-----
    lines = np.array([4958.92, 5006.84])
    doublet = np.exp(-0.5*((lam - lines[1])/sigma)**2) + 0.35*np.exp(-0.5*((lam - lines[0])/sigma)**2)
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OIII]5007d') # single template for this doublet
    line_wave = np.append(line_wave, lines[1])

    #                  -----[OI]-----
    lines = np.array([6363.776, 6300.304])
    doublet = np.exp(-0.5*((lam - lines[1])/sigma)**2) + 0.33*np.exp(-0.5*((lam - lines[0])/sigma)**2)
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[OI]6300d') # single template for this doublet
    line_wave = np.append(line_wave, lines[1])

    #                 -----[NII]-----
    lines = np.array([6548.05, 6583.46])
    doublet = np.exp(-0.5*((lam - lines[1])/sigma)**2) + 0.34*np.exp(-0.5*((lam - lines[0])/sigma)**2)
    emission_lines = np.column_stack([emission_lines, doublet])
    line_names = np.append(line_names, '[NII]6583d') # single template for this doublet
    line_wave = np.append(line_wave, lines[1])

    # Only include lines falling within the estimated fitted wavelength range.
    # This is important to avoid instabilities in the PPXF system solution
    #
    w = (line_wave > lamRange_gal[0]) & (line_wave < lamRange_gal[1])
    emission_lines = emission_lines[:, w]
    line_names = line_names[w]
    line_wave = line_wave[w]

    print('Emission lines included in gas templates:')
    print(line_names)

    return emission_lines, line_names, line_wave

#------------------------------------------------------------------------------
