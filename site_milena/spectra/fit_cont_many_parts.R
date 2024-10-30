remove(list = ls())

source('interp.R')

star_name = 'caT'

# -----------------------------------------------
# LOG FILE
# -----------------------------------------------
log_file = sprintf('%s.log', star_name)
write(sprintf('Star name: %s', star_name), file = log_file)

#################################################
# -----------------------------------------------
# PARAMETERS FOR NORMALIZATION
# -----------------------------------------------
#################################################
l_breaks = c(8500, 8600, 8700, 8800) #n de cortes no espectro
Nparts = length(l_breaks) - 1 
flag_norm = vector(length = Nparts); flag_norm[] = 1
nsig_up_b = vector(length = Nparts); nsig_up_b[] = 6 
nsig_low_b = vector(length = Nparts); nsig_low_b[] = 0.3
degree_b = vector(length = Nparts); degree_b[] = 1
span_b = vector(length = Nparts); span_b[] = 0.5
Nit_b = vector(length = Nparts); Nit_b[] = 3

nsig_up_b[c(1, 4, 5, 6, 8, 12)] = 1.2
span_b[c(3, 4, 5)] = 0.9
parts_to_correct = c(1 : Nparts)
factor_number = 0.037

# -----------------------------------------------
#################################################
spec_in_file = 'SSP_Fe+0.13_a+0.00_C+0.00_N+0.00_O+0.00_Mg+0.00_Si+0.00_Ca+0.00_Ti+0.00_Na+0.00_Al+0.00_Ba+0.00_Eu+0.00_age09.0_slope1.30'
spec_out_file = 'caTbase.out'
fig_out_file = 'caTbase.pdf' 

spec_in = read.table(spec_in_file)
spec_in$V3 <- rep(1e-6, times=length(spec_in$V1)) 
lmin = l_breaks[1]
lmax = l_breaks[Nparts + 1]
xx.xx = spec_in$V1 > lmin & spec_in$V1 < lmax
spec_in = spec_in[xx.xx, ]

flux_cont = vector()
flux_norm = vector()
err_norm = vector()
flux_cont_lmin = vector()
flux_cont_lmax = vector()

# -----------------------------------------------
# INITIALIZE FIGURE
# -----------------------------------------------
cairo_pdf(fig_out_file, width = 12, height = 8)
par(mfrow = c(3, 1), mar = c(5, 5, 1, 1), cex.lab = 1.5, cex.axis = 1.5)
par(font = list(family = 'Times'))
plot(spec_in$V1, spec_in$V2, type = 'l', log = 'y', #xlim = c(9500, 11000), #ylim = c(0.1, 0.2),
     xlab = expression(paste(lambda, ' (', ring(A), ')')), ylab = 'Flux',
     xaxt = 'n')
abline(v = l_breaks, lty = 5, col = 'red', lwd = 0.5)

# -----------------------------------------------
# FIT CONTINUUM IN PARTS
# -----------------------------------------------
for(b in 1:Nparts){
  lmin = l_breaks[b]
  lmax = l_breaks[b + 1]
  xx.xx = spec_in$V1 > lmin & spec_in$V1 < lmax
  spec_in_b = spec_in[xx.xx, ]
  
  print(c(min(spec_in_b$V1), max(spec_in_b$V1)))
  
  legend(lmin-300, factor_number * max(spec_in$V2) / min(spec_in$V2) * min(spec_in$V2), sprintf('%.0f', b), bty = 'n', cex = 1.2)
  if(flag_norm[b] == 1){
    ##################
    # 
    ##################
    nsig_up = nsig_up_b[b]
    nsig_low = nsig_low_b[b]
    degree = degree_b[b]
    span = span_b[b]
    Nit = Nit_b[b]
    
    write(sprintf('--------- Wavelenght range %.0f ---------', b), file = log_file, append = T)
    write(sprintf('lmin = %.1f, lmax = %.1f', lmin, lmax), file = log_file, append = T)
    write(sprintf('nsig_up = %.1f, nsig_low = %.1f, degree = %.0f, span = %.2f, Nit = %.0f', 
                  nsig_up, nsig_low, degree, span, Nit), file = log_file, append = T)
    
    x2_temp = spec_in_b$V2
    x1_temp = spec_in_b$V1
    
    for(nn in 1:Nit){
      temp <- loess(x2_temp ~ x1_temp, degree = degree, span = span)
      xx.xx_temp1 = x2_temp - predict(temp) <= nsig_up * sd(x2_temp - predict(temp)) & x2_temp - predict(temp) >= 0
      xx.xx_temp2 = abs(x2_temp - predict(temp)) <= nsig_low * sd(x2_temp - predict(temp)) & x2_temp - predict(temp) < 0
      
      xx.xx_temp = xx.xx_temp1 | xx.xx_temp2
      
      x1_temp = x1_temp[xx.xx_temp]
      x2_temp = x2_temp[xx.xx_temp]
    }
    
    temp <- loess(x2_temp ~ x1_temp, degree = degree, span = span)
    x1_temp = seq(lmin, lmax, 1)
    flux_cont_temp = predict(temp, newdata = data.frame(x1_temp = x1_temp))
    
    xx.xx = !is.na(flux_cont_temp)
    temp = interp(x1_temp[xx.xx], spec_in_b$V1, flux_cont_temp[xx.xx])
    flux_cont = c(flux_cont, temp)
    flux_cont_lmin = c(flux_cont_lmin, temp[1])
    flux_cont_lmax = c(flux_cont_lmax, temp[length(temp)])
  }else{
    
    flux_norm = c(flux_norm, spec_in_b$V2 * 0 + 1)
    if(obs_data){err_norm = c(err_norm, spec_in_b$V3 * 0)}
    
  }
}

flux_norm = spec_in$V2 / flux_cont
err_norm = spec_in$V3 / flux_cont

##################
# PLOT
##################
write('-----------------------------------', file = log_file, append = T)
write(paste('Flux corrected for interval: ', parts_to_correct), file = log_file, append = T)

# -----------------------------------------------
# CORRECT "STEPS" IN FLUX
# -----------------------------------------------
flux_corr = spec_in
flux_cont_corr = flux_cont
if(!(0 %in% parts_to_correct)){
  for(b in parts_to_correct){
    ##################
    # 
    ##################
    lmin = l_breaks[b]
    lmax = l_breaks[b + 1]
    xx.xx = spec_in$V1 > lmin & spec_in$V1 < lmax
    #xx.xx = spec_in$V1 < lmax
    spec_in_b = spec_in[xx.xx, ]
    m_factor = flux_cont_lmin[b + 1] / flux_cont_lmax[b]

    
    flux_corr$V2[xx.xx] = flux_corr$V2[xx.xx] * m_factor
    flux_corr$V3[xx.xx] = flux_corr$V3[xx.xx] * m_factor
    flux_cont_corr[xx.xx] = flux_cont_corr[xx.xx] * m_factor
  }
}
xx.xx = spec_in$V1 < lmax
lines(flux_corr$V1[xx.xx], flux_corr$V2[xx.xx], col = 'gray', lty = 1, lwd = 0.8)
lines(flux_corr$V1[xx.xx], flux_cont_corr[xx.xx], col = rgb(1, 0.5, 0.5), lty = 1, lwd = 0.1)
lines(spec_in$V1, spec_in$V2)
lines(spec_in$V1, flux_cont, col = 'red', lwd = 0.5)

leg1 = 'Original spectra'
leg2 = 'Flux-corrected spectra'
leg3 = 'Continuum fit'
leg4 = 'Intervals for continuum fit'
legend('topright', c(leg1, leg2, leg3, leg4), col = c('black', 'gray', 'red', 'red'), lty = c(1, 1, 1, 5), bty = 'n', cex = 1.2)

# -----------------
# AXIS
# -----------------
axis(1, seq(6000, 30000, 1000), tcl = -0.6, labels = T)
axis(1, seq(6000, 30000, 500), tcl = -0.3, labels = F)


# -----------------------------------------------
# PLOT NORMALIZED SPECTRUM
# -----------------------------------------------
plot(spec_in$V1, flux_norm, type = 'l', #xlim = c(11000, 14000), ####################################################
     xlab = expression(paste(lambda, ' (', ring(A), ')')), ylab = 'Normalized Flux',
     xaxt = 'n')
abline(h = 1, lty = 5, col = 'red')
abline(h = c(0.98, 1.02), lty = 3, col = 'red')
abline(v = l_breaks, lty = 5, col = 'red', lwd = 0.5)

# -----------------
# AXIS
# -----------------
axis(1, seq(6000, 30000, 1000), tcl = -0.6, labels = T)
axis(1, seq(6000, 30000, 500), tcl = -0.3, labels = F)

#################################################
# -----------------------------------------------
# WRITE OUTPUT FILE WITH NORMALIZED SPECTRUM
# -----------------------------------------------
#################################################
spec_out = spec_in

spec_out$V4 = flux_corr$V2
spec_out$V5 = flux_corr$V3
spec_out$V6 = flux_norm
spec_out$V7 = err_norm
spec_out$V8 = flux_cont
spec_out$V9 = flux_cont_corr
names(spec_out) = c('lambda', 'flux_orig', 'err_orig', 'flux_corr', 'err_corr', 'flux_norm', 'err_norm', 'flux_cont', 'flux_cont_corr')

write.table(na.omit(spec_out), file = spec_out_file, row.names = F, col.names = T)
# -----------------------------------------------
#################################################

# -----------------------------------------------
# PLOT CORRECTED SPECTRUM (NO "STEPS")
# -----------------------------------------------
t = read.table(spec_out_file, header = T)
plot(t$lambda, t$flux_corr, type = 'l', log = 'y', # xlim = c(9500, 11000), ylim = c(0.1, 0.2),
     xlab = expression(paste(lambda, ' (', ring(A), ')')), ylab = 'Flux')
lines(t$lambda, t$flux_cont_corr, col = 'red', lwd = 0.5)

leg1 = 'Flux-corrected spectra'
leg2 = 'Continuum fit'
legend('topright', c(leg1, leg2), col = c('black', 'red'), lty = 1, bty = 'n', cex = 1.2)

par(font = 2)
legend('bottomleft', star_name, bty = 'n', cex = 2)
par(font = 1)
dev.off()