spec_ref = read.table('na8190/SSP_Fe+0.13_a+0.00_C+0.00_N+0.00_O+0.00_Mg+0.00_Si+0.00_Ca+0.00_Ti+0.00_Na+0.00_Al+0.00_Ba+0.00_Eu+0.00_age09.0_slope1.30')
spec_imf = read.table('na8190/SSP_Fe+0.13_a+0.00_C+0.00_N+0.00_O+0.00_Mg+0.00_Si+0.00_Ca+0.00_Ti+0.00_Na+0.00_Al+0.00_Ba+0.00_Eu+0.00_age09.0_slope3.30')
spec_afe = read.table('na8190/SSP_Fe+0.13_a+0.40_C+0.00_N+0.00_O+0.00_Mg+0.00_Si+0.00_Ca+0.00_Ti+0.00_Na+0.00_Al+0.00_Ba+0.00_Eu+0.00_age09.0_slope1.30')

xx.xx = spec_ref$V1 < 1000000
temp = spec_ref$V2[xx.xx]
for(i in 1:10){
  xx.xx = temp > median(temp) - 0.5 * sd(temp)
  temp = temp[xx.xx]
}
cont_ref = median(temp)

xx.xx = spec_imf$V1 < 1000000
temp = spec_imf$V2[xx.xx]
for(i in 1:10){
  xx.xx = temp > median(temp) - 0.5 * sd(temp)
  temp = temp[xx.xx]
}
cont_imf = median(temp)

xx.xx = spec_afe$V1 < 1000000
temp = spec_afe$V2[xx.xx]
for(i in 1:10){
  xx.xx = temp > median(temp) - 0.5 * sd(temp)
  temp = temp[xx.xx]
}
cont_afe = median(temp)

spec_ref$V3 = spec_ref$V2 / cont_ref
spec_imf$V3 = spec_imf$V2 / cont_imf
spec_afe$V3 = spec_afe$V2 / cont_afe
par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
plot(spec_ref$V1, spec_ref$V2 / cont_ref, type = 'l')
lines(spec_imf$V1, spec_imf$V2 / cont_imf, col = 'red')
lines(spec_afe$V1, spec_afe$V2 / cont_afe, col = 'blue')

plot(spec_imf$V1, (spec_imf$V3 - spec_ref$V3) / spec_ref$V3 * 100, type = 'l', col = 'red')
lines(spec_afe$V1, (spec_afe$V3 - spec_ref$V3) / spec_ref$V3 * 100,  col = 'blue')
abline(h = 0)
