remove(list = ls())

#print(Sys.time())

#-----------------------------------
# Functions
#-----------------------------------

source('interp.R')

for(model_type in c('ref')){
  #-----------------------------------
  # SOME DEFINITIONS
  #-----------------------------------
  fwhm_in = 0.5    # resolution of the input spectrum, in A
  sigma_out_kms = 160  # Desired resolution, in km/s
  delta_sampl = 1  # Sampling of the output spectrum, in Angstroms
  folder = 'tiO1' 
  model_type = 'ref'  # ref, imf or afe
  
  #-----------------------------------
  # Input/output files
  #-----------------------------------
  #if(model_type == 'ref'){
  input_spec = read.table(sprintf('%s/SSP_Fe+0.13_a+0.00_C+0.00_N+0.00_O+0.00_Mg+0.00_Si+0.00_Ca+0.00_Ti+0.00_Na+0.00_Al+0.00_Ba+0.00_Eu+0.00_age09.0_slope1.30', folder))
  #}
  #if(model_type == 'imf'){
  #  input_spec = read.table(sprintf('%s/SSP_Fe+0.13_a+0.00_C+0.00_N+0.00_O+0.00_Mg+0.00_Si+0.00_Ca+0.00_Ti+0.00_Na+0.00_Al+0.00_Ba+0.00_Eu+0.00_age09.0_slope3.30', folder))
  #}
  #if(model_type == 'afe'){
  #  input_spec = read.table(sprintf('%s/SSP_Fe+0.13_a+0.40_C+0.00_N+0.00_O+0.40_Mg+0.40_Si+0.40_Ca+0.40_Ti+0.40_Na+0.00_Al+0.00_Ba+0.00_Eu+0.00_age09.0_slope1.30', folder))
  #}
  output_file = sprintf('%s/spec_%s_conv%.0f.txt', folder, model_type, sigma_out_kms)
  
  
  # input_spec = read.table('CaT/SSP_Fe+0.40_a+0.00_C+0.00_N+0.00_O+0.00_Mg+0.00_Si+0.00_Ca+0.00_Ti+0.00_Na+0.00_Al+0.00_Ba+0.00_Eu+0.00_age09.0_slope1.30')
  # input_spec = read.table('CaT/SSP_Fe+0.40_a+0.00_C+0.00_N+0.00_O+0.00_Mg+0.00_Si+0.00_Ca+0.00_Ti+0.00_Na+0.00_Al+0.00_Ba+0.00_Eu+0.00_age09.0_slope3.30')
  # input_spec = read.table('CaT/SSP_Fe+0.40_a+0.40_C+0.00_N+0.00_O+0.40_Mg+0.40_Si+0.40_Ca+0.40_Ti+0.40_Na+0.00_Al+0.00_Ba+0.00_Eu+0.00_age09.0_slope1.30')
  # output_file = sprintf('%s/conv_feh_0.40/spec_%s_conv%.0f.txt', folder, model_type, sigma_out_kms)
  
  #-----------------------------------
  # SPECTRUM
  #-----------------------------------
  lambda = input_spec$V1
  Npix = length(lambda)
  flux = input_spec$V2
  
  #-----------------------------------
  # SIGMAS
  #-----------------------------------
  sigma_in = fwhm_in / (2 * sqrt(2 * log(2)))
  sigma_out = sigma_out_kms / 3e5 * lambda
  sigma_conv = sqrt(sigma_out**2 - sigma_in**2)
  sigma_conv[(sigma_out**2 - sigma_in**2) <= 0] = 0.001
  print(sigma_conv[1])
  #-----------------------------------
  # CONVOLUTION
  #-----------------------------------
  flux_conv = flux; flux_conv[] = 0
  for(n in 1:Npix){
    dl2conv = 20 * sigma_conv[n]
    pixels2sum = which(lambda >= (lambda[n] - dl2conv) & lambda <= (lambda[n] + dl2conv))
    psf = exp(-((lambda[n] - lambda[pixels2sum]) / sigma_conv[n])**2 / 2)
    norm_psf = sum(psf)
    print(norm_psf)
    for(m in pixels2sum){
      if((m+1) <= Npix){
        psf = exp(-((lambda[n] - lambda[m]) / sigma_conv[n])**2 / 2)
        flux_conv[n] = flux_conv[n] + flux[m] * psf / norm_psf
      }
    #print(sigma_conv)
    }
    # ------------------------------------------------------------------
    # METHOD 2
    # ------------------------------------------------------------------
    # THE METHOD 2 BELOW WORKS WELL FOR SIGMA_OUT > 1.01 * SIGMA_IN
    #   WHEN SIGMA_OUT <~ SIGMA_IN, IT DOES NOT WORK
    #   DIFFERENCES BETWEEN METHODS
    #   ARE WITHIN < 0.5% FOR SIGMA_OUT > 1.01 * SIGMA_IN
    #   SEE DIRECTORIES 
    #        conv_test_fortran_code AND 
    #        Python
    # ------------------------------------------------------------------
    # dl2conv = 10 * sigma_conv[n] 
    # pixels2sum = which(lambda >= (lambda[n] - dl2conv) & lambda <= (lambda[n] + dl2conv))
    # for(m in pixels2sum){
    #   psf = exp(-((lambda[n] - lambda[m]) / sigma_conv[n])**2 / 2) / (sigma_conv[n] * sqrt(2 * pi)) * (lambda[m+1] - lambda[m])
    #   flux_conv[n] = flux_conv[n] + flux[m] * psf
    # }
  #print(m)
  }
  print(head(flux_conv,5))
  print(tail(flux_conv,5))
  
  c#-----------------------------------
  # REBIN
  #-----------------------------------
  lambda_out = seq(min(lambda), max(lambda), delta_sampl)
  flux_conv_interp = interp(lambda, lambda_out, flux_conv)
  
  #-----------------------------------
  # SAVE CONVOLVED SPECTRUM
  #-----------------------------------
  output_spec = data.frame(l = lambda_out, f = flux_conv_interp)
  write.table(output_spec, file = output_file, row.names = F, col.names = F)
  
  #-----------------------------------
  # PLOT SPECTRA
  #-----------------------------------
  plot(lambda, flux, type = 'l')
  lines(lambda, flux_conv, col = 'red')
  lines(lambda_out, flux_conv_interp, col = 'green')
  
  #-----------------------------------
  # COMPARE WITH OBSERVATIONS
  #-----------------------------------
  # t = read.table('Stack_100-160.dat')
  # t = read.table('Stack_240-400.dat')
  # xx.xx = t$V1 > 8400 & t$V1 < 8500
  # plot(t$V1, t$V2/median(t$V2[xx.xx]), type = 'l', xlab = expression(paste(lambda, ' (', ring(A), ')')),  
  #      ylab = 'Flux (arbitrary units)')
  # xx.xx = lambda > 8400 & lambda < 8500
  # lines(lambda, flux_conv/median(flux_conv[xx.xx]), col = 'red')
  # lines(lambda, flux/median(flux[xx.xx]))
  # legend('bottomright', c('Input', 'Output convolved'), col = c('black', 'red'), bty = 'n', cex = 1.6, lty = 1)
  
  #print(Sys.time())
}

