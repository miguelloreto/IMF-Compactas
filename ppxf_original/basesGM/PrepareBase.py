import numpy as np


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

spec,code=np.loadtxt('BASE.IAA2012.All.edt',usecols=(0,3),dtype='str',unpack=True)
age,z,ms,Yav,aFe=np.loadtxt('BASE.IAA2012.All.edt',usecols=(1,2,4,5,6),unpack=True)


agesForBase=[1E6,5E6,10E6,15E6,20E6,30E6,50E6,100E6,200E6,300E6,400E6,500E6,600E6,700E6,800E6,900E6,1E9,2E9,5E9,9E9,13E9]
baseagesIn=[]
for a in agesForBase:
       baseagesIn.append(find_nearest(age,a))
       
tobase_spec=[]
tobase_ages=[]
tobase_met=[]
tobase_code=[]
tobase_Mstar=[]
tobase_Yav=[]
tobase_alf=[]
for t in baseagesIn:
#      cont=0
#      for aa in age:
 #         if aa == t:
            #print spec[cont], '{:.2E}'.format(age[cont]), '{:.6F}'.format(z[cont]),code[cont],'{:.4F}'.format(ms[cont]),'{:.2F}'.format(Yav[cont]),'{:.2F}'.format(aFe[cont])
            tt=spec[t==age]
            [tobase_spec.append(tt[i]) for i in np.arange(0,(len(tt)+0))]
            tt=age[t==age]
            [tobase_ages.append(tt[i]) for i in np.arange(0,(len(tt)+0))]
            tt=z[t==age]
            [tobase_met.append(tt[i]) for i in np.arange(0,(len(tt)+0))]
            tt=code[t==age]
            [tobase_code.append(tt[i]) for i in np.arange(0,(len(tt)+0))]
            tt=ms[t==age]
            [tobase_Mstar.append(tt[i]) for i in np.arange(0,(len(tt)+0))]
            tt=Yav[t==age]
            [tobase_Yav.append(tt[i]) for i in np.arange(0,(len(tt)+0))]
            tt=aFe[t==age]
            [tobase_alf.append(tt[i]) for i in np.arange(0,(len(tt)+0))]


f=open('BaseGMManga2','w+')

f.write(str(len(tobase_ages))+'   [N_base]\n')

indexes=np.array(tobase_met).argsort()

for i in indexes:
   s=(30-len(tobase_spec[i]))*' '
   ss=(20-len(tobase_code[i]))*' '
   f.write(tobase_spec[i]+s+'{:.2E}'.format(tobase_ages[i])+'   {:.8F}'.format(tobase_met[i])+ss
          +tobase_code[i]+'   {:.5F}'.format(tobase_Mstar[i])+'   {:.2F}'.format(tobase_Yav[i])+'   {:.2F}'.format(tobase_alf[i])+'\n')


f.close()

