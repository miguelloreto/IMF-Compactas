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

def extract_data(data):
  low_s, high_s=data[0], data[1]
  low=[np.array(low_s['lam']), np.array(low_s['fluxo'])]
  high=[np.array(high_s['lam']), np.array(high_s['fluxo'])]
  return(low, high)
  
def convolui_espectro(modelo, sigma_out):
  '''Convolui um espectro com sigma de saída sigma_out E dá um novo vetor lambda conforme deltalambda'''
  lam, fluxo=modelo['lam'], modelo['fluxo']
  lam=np.array(lam)
  fluxo=np.array(fluxo)
  #Rebin no lambda:
  deltalambda=1
  lambda_interp=np.arange(lam.min(), lam.max(), deltalambda)
  #Convolução:
  sigma_in = 0.5 #Resolução do espectro de input em Angstrons
  sigma_in= sigma_in/(2 * math.sqrt(2 * np.log(2)))
  sigma_out=sigma_out/c *lam
  sigma_conv=np.sqrt(np.array(sigma_out)**2 - np.array(sigma_in)**2)
  sigma_conv[(sigma_out**2 - sigma_in**2) <= 0] = 0.001
  convoluido=np.zeros(len(lam))
  for i in range(0, len(lam)):
    intervalo=20*sigma_conv[i]
    pixels2sum=np.where(np.logical_and(lam>=(lam[i]-intervalo), lam<=(lam[i]+intervalo)))[0]
    psf=np.exp(-((lam[i] - lam[pixels2sum]) / sigma_conv[i])**2 / 2)
    norm_psf=sum(psf) 
    for k in pixels2sum:
      if((k+2) <= len(lam)):
        psf=np.exp(-((lam[i] - lam[k]) / sigma_conv[i])**2 / 2)
        convoluido[i]= convoluido[i] +fluxo[k]*psf/norm_psf
  fluxo_interp=np.interp(lambda_interp,lam,convoluido)
  modelo=pd.DataFrame(data={'lam': lambda_interp, 'fluxo':fluxo_interp})
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

def plota_spec(ref, observado, base, alfa, imf, info, nome_linha):
  centralmin, centralmax, bluemin, bluemax, redmin, redmax= info['info'].iloc[0],info['info'].iloc[1],info['info'].iloc[2],info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5],
  observado_l, observado_h=extract_data(observado)
  base_l, base_h=extract_data(base)
  imf_l, imf_h=extract_data(imf)
  alfa_l, alfa_h=extract_data(alfa)
  ref_x, ref_y=ref['lam'], ref['fluxo'] 
  fig, (ax1, ax2)=plt.subplots(nrows=2,ncols=1)
  fig.suptitle(nome_linha)
  ax1.plot(ref_x, ref_y, color='xkcd:dusty blue', lw=0.8)
  ax1.plot(base_l[0], base_l[1], color='xkcd:black')
  ax1.plot(alfa_l[0], alfa_l[1], color='xkcd:primary blue')
  ax1.plot(imf_l[0], imf_l[1], color='xkcd:cherry red')
  ax1.plot(observado_l[0], observado_l[1], color='xkcd:green', ls='dotted')
  ax1.set_xlim(bluemin-5, redmax+5)
  ax1.legend(labels=[r'Observado', r'Base', r'$[α/Fe]$', r'IMF', r'Sintético']  ,loc="lower right")
  ax1.set_xlabel(r' $\lambda$($\AA$)')
  ax1.set_ylabel(r'Fluxo Normalizado')
  ax1.set_ylim(0.55, 1.1)
  ax1.fill_betweenx([0,1.2], bluemin, bluemax, facecolor='xkcd:gunmetal', edgecolor='black', alpha=0.4)
  ax1.fill_betweenx([0,1.2], centralmin, centralmax, facecolor='xkcd:greyish', edgecolor='black', alpha=0.4)
  ax1.fill_betweenx([0,1.2], redmin, redmax, facecolor='xkcd:gunmetal',edgecolor='black', alpha=0.4)
  ax1.text(bluemax+5, 0.6, r'σ=160km/s', fontsize=12)

  ax2.plot(ref_x, ref_y, color='xkcd:dusty blue', lw=0.8)
  ax2.plot(base_h[0], base_h[1], color='xkcd:black')
  ax2.plot(alfa_h[0], alfa_h[1], color='xkcd:primary blue')
  ax2.plot(imf_h[0], imf_h[1], color='xkcd:cherry red')
  ax2.plot(observado_h[0], observado_h[1], color='xkcd:green', ls='dotted')
  ax2.set_xlim(bluemin-5, redmax+5)
  ax2.legend(labels=[r'Observado', r'Base', r'$[α/Fe]$', r'IMF', r'Sintético']  ,loc="lower right")
  ax2.set_xlabel(r' $\lambda$($\AA$)')
  ax2.set_ylabel(r'Fluxo Normalizado')
  ax2.set_ylim(0.55, 1.1)
  ax2.fill_betweenx([0,1.2], bluemin, bluemax, facecolor='xkcd:gunmetal', edgecolor='black', alpha=0.4)
  ax2.fill_betweenx([0,1.2], centralmin, centralmax, facecolor='xkcd:greyish', edgecolor='black', alpha=0.4)
  ax2.fill_betweenx([0,1.2], redmin, redmax, facecolor='xkcd:gunmetal',edgecolor='black', alpha=0.4)
  ax2.text(bluemax+5,0.6, r'σ=300km/s', fontsize=12)
  plt.subplots_adjust(hspace=0.4)
  fig.set_figheight(5)
  fig.set_figwidth(20)
  fig.savefig(nome_linha+'.png', bbox_inches='tight')


for i in range(0, len(stacks)):
  stackado=pd.read_csv(stacks[i], sep='\s+', header=None, usecols=[0,1], names=['lam', 'fluxo'])
  stackado=stackado[stackado.fluxo > 0]
  stacks[i]=stackado

#print(index)
#index=[index[5], 'mg4780']
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
      base=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5], base)
      base_sigma[j]=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5], base_sigma[j])
      alfa_sigma[j]=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5], alfa_sigma[j])
      imf_sigma[j]=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5], imf_sigma[j])

    stacks_normalizados=[]
    for stack in stacks:
      stack_norm=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5], stack)
      stacks_normalizados.append(stack_norm)

    plota_spec(base, stacks_normalizados, base_sigma, alfa_sigma, imf_sigma, info, i)