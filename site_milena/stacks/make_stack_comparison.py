import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
import pandas as pd
import os
import math

c=3e5 #Velocidade da luz em Km/s
pi=np.pi
sigmas=[160, 300]

stack_path=os.getcwd()
stack_dir=os.listdir(stack_path)
index=[i for i in stack_dir if not i.__contains__('.')]
stacks=[i for i in stack_dir if i.__contains__('Stack')]
stacks.sort()

def gaussiana(lam, lam_0, dellambda):
  return np.exp(-np.power((lam - lam_0)/dellambda, 2)/2)

def extract_data(data):
  low_s, high_s=data[0], data[1]
  low=[np.array(low_s['lam']), np.array(low_s['high'])]
  high=[np.array(high_s['lam']), np.array(high_s['high'])]
  return(low, high)
  
def convolui_espectro(modelo, sigma_out):
  lam, fluxo=modelo['lam'], modelo['fluxo']
  lam=np.array(lam)
  fluxo=np.array(fluxo)
  sigma_in = 0.5 #Resolução do espectro de input em Angstrons
  sigma_in= sigma_in/(2 * math.sqrt(2 * np.log(2)))
  sigma_out=sigma_out/c *lam
  sigma_conv=np.sqrt(np.array(sigma_out)**2 - np.array(sigma_in)**2)
  sigma_conv[(sigma_out**2 - sigma_in**2) <= 0] = 0.001
  convoluido=np.zeros(len(lam))
  for i in range(0, len(lam)):
    intervalo=10*sigma_conv[i]
    j=np.where(np.logical_and(lam<=(lam[i]+intervalo), lam>=(lam[i]-intervalo)))
    j=np.array(j[0])
    psf=gaussiana(lam[i], lam[j], sigma_conv[i])
    norm_psf=sum(psf) 
    for k in j:
      if((k+1) <= len(lam)):
        psf=gaussiana(lam[i], lam[k], sigma_conv[i])
        convoluido[i]= convoluido[i] +fluxo[k]*psf/norm_psf
  modelo['fluxo']=convoluido
  return(modelo)


def normaliza_continuo_simples(blue_min, blue_max, red_min, red_max, espectro):
  bluecont=espectro['lam'].between(blue_min, blue_max, inclusive='both')
  redcont=espectro['lam'].between(red_min, red_max, inclusive='both')
  x=[(espectro['lam'][bluecont.values]).median(),(espectro['lam'][redcont.values]).median()]
  y=[(espectro['fluxo'][bluecont.values]).median(),(espectro['fluxo'][redcont.values]).median()]
  a,b=np.polyfit(x,y,1)
  lam=np.array(espectro['lam'])
  cont= lambda a, b, lam: a*lam + b
  cont=cont(a,b,lam)
  espectro['fluxo']=np.array(espectro['fluxo']/cont)
  return(espectro)

def plota_spec(observado, base, alfa, imf, info, nome_linha):
  centralmin, centralmax, bluemin, bluemax, redmin, redmax= info['info'].iloc[0],info['info'].iloc[1],info['info'].iloc[2],info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5],
  observado_l, observado_h=extract_data(observado)
  base_l, base_h=extract_data(base)
  imf_l, imf_h=extract_data(imf)
  alfa_l, alfa_h=extract_data(alfa)


  fig, (ax1, ax2)=plt.subplots(nrows=2,ncols=1)
  fig.suptitle(nome_linha)
  ax1.plot(observado_l[0], observado_l[0], color='xkcd:green', ls='dashed')
  ax1.plot(base_l[0], base_l[1], color='xkcd:black')
  ax1.plot(alfa_l[0], alfa_l[1], color='xkcd:primary blue')
  ax1.plot(imf_l[0], imf_l[1], color='xkcd:cherry red')
  ax1.set_xlim(bluemin-20, redmax+20)
  ax1.legend(labels=[r'Base', r'$[α/Fe]$', r'IMF'] ,loc="upper right")
  ax1.set_xlabel(r' $\lambda$($\AA$)')
  ax1.set_ylabel(r'Fluxo Normalizado')
  ax1.set_ylim(0.6, 1)
  ax1.fill_betweenx([0,1.2], bluemin, bluemax, facecolor='xkcd:gunmetal', edgecolor='black', alpha=0.4)
  ax1.fill_betweenx([0,1.2], centralmin, centralmax, facecolor='xkcd:greyish', edgecolor='black', alpha=0.4)
  ax1.fill_betweenx([0,1.2], redmin, redmax, facecolor='xkcd:gunmetal',edgecolor='black', alpha=0.4)

  ax2.plot(observadox, observadoy, color='xkcd:green', ls='dashed')
  ax2.plot(basex, basey, color='xkcd:black')
  ax2.plot(chemx, chemy, color='xkcd:primary blue')
  ax2.plot(imfx, imfy, color='xkcd:cherry red')
  ax2.set_xlim(bluemin-20, redmax+20)
  ax2.legend(labels=[r'Base', r'$[α/Fe]$', r'IMF'] ,loc="upper right")
  ax2.set_xlabel(r' $\lambda$($\AA$)')
  ax2.set_ylabel(r'Fluxo Normalizado')
  ax2.set_ylim(0.6, 1)
  ax2.fill_betweenx([0,1.2], bluemin, bluemax, facecolor='xkcd:gunmetal', edgecolor='black', alpha=0.4)
  ax2.fill_betweenx([0,1.2], centralmin, centralmax, facecolor='xkcd:greyish', edgecolor='black', alpha=0.4)
  ax2.fill_betweenx([0,1.2], redmin, redmax, facecolor='xkcd:gunmetal',edgecolor='black', alpha=0.4)
  plt.subplots_adjust(hspace=0.4)
  fig.set_figheight(5)
  fig.set_figwidth(20)
  fig.savefig(nome_linha+'.png', bbox_inches='tight')


for i in range(0, len(stacks)):
  stackado=pd.read_csv(stacks[i], sep='\s+', header=None, usecols=[0,1], names=['lam', 'fluxo'])
  stackado=stackado[stackado.fluxo > 0]
  stacks[i]=stackado

for i in index:
    index_path=stack_path+'/'+i
    index_dir=os.listdir(index_path)
    ssps=[i for i in index_dir if i.__contains__('SSP')]
    if len(ssps) == 0:
        continue
    else:
      os.chdir(index_path)
      info=[i for i in index_dir if i.__contains__('.txt')]
      info=pd.read_csv(info[0], header=None, names=['name_info', 'info'])
      alfa=pd.read_csv(ssps[0], sep='\s+', skiprows=1, header=None, names=['lam', 'fluxo'])
      imf=pd.read_csv(ssps[1], sep='\s+', skiprows=1, header=None, names=['lam', 'fluxo'])
      base=pd.read_csv(ssps[2], sep='\s+', skiprows=1, header=None, names=['lam', 'fluxo'])
    
    lambda_interpol=np.array(base['lam'])

    base_sigma=[]
    alfa_sigma=[]
    imf_sigma=[]

    print(i)
    for j in sigmas:
      base_sigma.append(convolui_espectro(base, j))
      alfa_sigma.append(convolui_espectro(alfa, j))
      imf_sigma.append(convolui_espectro(imf, j))

    for j in range(0, len(base_sigma)):
      base_sigma[j]=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5], base_sigma[j])
      alfa_sigma[j]=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5], alfa_sigma[j])
      imf_sigma[j]=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5], imf_sigma[j])

    stacks_interpolados=[]
    for stack in stacks:
      stack=stack[(stack.lam > lambda_interpol[0]) & (stack.lam < lambda_interpol[-1])]
      interpolador=CubicSpline(np.array(stack.lam), np.array(stack.fluxo))
      fluxo_interpol=interpolador(lambda_interpol)
      stack_interpolado=pd.DataFrame(data={'lam': lambda_interpol, 'fluxo':fluxo_interpol})
      stack_interpolado=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5], stack_interpolado)
      stacks_interpolados.append(stack_interpolado)
    
    plota_spec(stacks_interpolados, base_sigma, alfa_sigma, imf_sigma, info, i)