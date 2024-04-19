t = read.table('BASE.IAA2012.All', stringsAsFactors = F)
Nbase = nrow(t)

age = t$V2
LM = vector()
for(i in 1:Nbase){
  temp = read.table(t$V1[i])
  LM[i] = temp$V2[temp$V1 == 4020]
}

plot(age, LM, log = 'xy')