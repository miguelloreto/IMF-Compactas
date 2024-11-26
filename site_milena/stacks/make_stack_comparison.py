import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
import pandas as pd
import os
import math

c=3e5 #Velocidade da luz em Km/s
pi=np.pi
sigmas=[160, 180, 200, 220, 240, 400]

stack_path=os.getcwd()
stack_dir=os.listdir(stack_path)
index=[i for i in stack_dir if not i.__contains__('.')]
stacks=[i for i in stack_dir if i.__contains__('Stack')]
stacks.sort()

def gaussiana(lam, lam_0, dellambda):
  return 1/(np.sqrt(2*pi)*dellambda)*np.exp(-np.power((lam - lam_0)/dellambda, 2)/2)

def convolui_espectro(modelo, sigma):
  lam, fluxo=modelo['lam'], modelo['fluxo']
  lam=np.array(lam)
  fluxo=np.array(fluxo)
  deltalambda=lam*(sigma/c) / (2 * math.sqrt(2 * np.log(2)))
  convoluido=np.zeros(len(lam))
  vetor_corte=[]
  for i in range(0, len(lam)):
    intervalo=2*deltalambda[i]
    j=np.where(np.logical_and(lam<=(lam[i]+intervalo), lam>=(lam[i]-intervalo)))
    psf=gaussiana(lam[i], lam[j], deltalambda[i])
    convoluido[i]=np.sum(psf*fluxo[j])
  modelo['fluxo']=convoluido
  return(modelo)


def normaliza_continuo_simples(lnorm_min, lnorm_max, espectro):
  bluecont=espectro['lam'].between(lnorm_min, lnorm_max, inclusive='both')
  fluxo_bluecont=espectro['fluxo'][bluecont.values]
  for i in range(1, 6):
    bluecont= fluxo_bluecont > (fluxo_bluecont.median() - 0.5*fluxo_bluecont.std()) #pd.Series de boleanos onde fluxobluecont é maior do que a mediana de fluxo blue cont menos metade de seu dp
    fluxo_bluecont=fluxo_bluecont[bluecont.values]
  cont=fluxo_bluecont.median()
  espectro['fluxo']=espectro['fluxo']/cont
  return(espectro)

def plota_spec(base, chem, imf, info_linha, nome_linha):
  centralmin, centralmax, bluemin, bluemax, redmin, redmax= info['info'].iloc[0],info['info'].iloc[1],info['info'].iloc[2],info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5],
  basex, basey = base['lam'], base['fluxo']
  chemx, chemy = chem['lam'], chem['fluxo']
  imfx, imfy = imf['lam'], imf['fluxo']
  minimo=np.min(chemy)
  maximo=np.max(chemy)
  comp_chem=(chemy-basey)/chemy
  comp_imf=(imfy-basey)/imfy
  fig, (ax1, ax2)=plt.subplots(nrows=2,ncols=1)
  fig.suptitle(nome_linha)
  ax1.plot(basex, basey, color='xkcd:black', lw=1.8)
  ax1.plot(chemx, chemy, color='xkcd:primary blue')
  ax1.plot(imfx, imfy, color='xkcd:cherry red')
  ax1.set_xlim(bluemin-20, redmax+20)
  ax1.legend(labels=[r'Base', r'$[α/Fe]$', r'IMF'] ,loc="upper right")
  ax1.set_xlabel(r' $\lambda$($\AA$)')
  ax1.set_ylabel(r'Fluxo Normalizado')
  ax1.set_ylim(minimo -0.15*minimo, maximo +0.05*maximo)
  ax1.fill_betweenx([0,1.2], bluemin, bluemax, facecolor='xkcd:gunmetal', edgecolor='black', alpha=0.4)
  ax1.fill_betweenx([0,1.2], centralmin, centralmax, facecolor='xkcd:greyish', edgecolor='black', alpha=0.4)
  ax1.fill_betweenx([0,1.2], redmin, redmax, facecolor='xkcd:gunmetal',edgecolor='black', alpha=0.4)

  ax2.plot(chemx, comp_chem, color='xkcd:primary blue')
  ax2.plot(imfx, comp_imf, color='xkcd:cherry red')
  ax2.set_xlim(bluemin-20, redmax+20)
  ax2.legend(labels=[r'$[α/Fe]$', r'IMF'], loc="upper right")
  ax2.set_xlabel(r' $\lambda$($\AA$)')
  ax2.set_ylabel(r' $\Delta$Fluxo / Fluxo')
  ax2.set_ylim(-0.18, 0.18)
  ax2.set_yticks([-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15])
  ax2.axhline(-0.05, ls='--', lw=1.3, color='xkcd:dark')
  ax2.axhline(0.05, ls='--', lw=1.3, color='xkcd:dark')
  ax2.axhline(0, ls='dotted', lw=1, color='xkcd:dark')
  ax2.fill_betweenx([-1,1], bluemin, bluemax, facecolor='xkcd:gunmetal',edgecolor='black', alpha=0.4)
  ax2.fill_betweenx([-1,1], centralmin, centralmax, facecolor='xkcd:greyish', edgecolor='black',alpha=0.4)
  ax2.fill_betweenx([-1,1], redmin, redmax, facecolor='xkcd:gunmetal',edgecolor='black', alpha=0.4)
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

    for j in sigmas:
      base_sigma.append(convolui_espectro(base, j))
      alfa_sigma.append(convolui_espectro(alfa, j))
      imf_sigma.append(convolui_espectro(imf, j))
    
    for j in range(0, len(base_sigma)):
      base_sigma[j]=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3], base_sigma[j])
      alfa_sigma[j]=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3], alfa_sigma[j])
      imf_sigma[j]=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3], imf_sigma[j])

    stacks_interpolados=[]
    for stack in stacks:
      stack=stack[(stack.lam > lambda_interpol[0]) & (stack.lam < lambda_interpol[-1])]
      interpolador=CubicSpline(np.array(stack.lam), np.array(stack.fluxo))
      fluxo_interpol=interpolador(lambda_interpol)
      stack_interpolado=pd.DataFrame(data={'lam': lambda_interpol, 'fluxo':fluxo_interpol})
      stack_interpolado=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3], stack_interpolado)
      stacks_interpolados.append(stack_interpolado)
    
    for j in range(0, len(stacks_interpolados)): plota_spec(stacks_interpolados[j], alfa_sigma[j], imf_sigma[j], info, i+'_'+str(sigmas[j]))