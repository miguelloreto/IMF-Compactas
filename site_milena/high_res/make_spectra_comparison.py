import matplotlib.pyplot as plt
import numpy as np
import sys
import shutil
import pandas as pd
import os

spectra_path=os.getcwd()
spectra_dir=os.listdir(spectra_path)
index=[i for i in spectra_dir if not i.__contains__('.')]

def plota_spec(base, chem, imf, info_linha, nome_linha):
        
  centralmin, centralmax, bluemin, bluemax, redmin, redmax= info['info'].iloc[0],info['info'].iloc[1],info['info'].iloc[2],info['info'].iloc[3],info['info'].iloc[4],info['info'].iloc[5],
  basex, basey = base['lambda'], base['fluxo']
  chemx, chemy = chem['lambda'], chem['fluxo']
  imfx, imfy = imf['lambda'], imf['fluxo']
  minimo=np.min(chemy)
  maximo=np.max(chemy)
  comp_chem=(chemy-basey)/basey
  comp_imf=(imfy-basey)/basey
  minimo_comp=np.min(comp_chem)
  maximo_comp=np.max(comp_chem)
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

def normaliza_continuo_simples(lnorm_min, lnorm_max, espectro):
  espectro.columns=['lambda','fluxo']
  bluecont=espectro['lambda'].between(lnorm_min, lnorm_max, inclusive='both')
  fluxo_bluecont=espectro['fluxo'][bluecont.values]
  for i in range(1, 10):
    bluecont= fluxo_bluecont > (fluxo_bluecont.median() - 0.5*fluxo_bluecont.std()) #pd.Series de boleanos onde fluxobluecont é maior do que a mediana de fluxo blue cont menos metade de seu dp
    fluxo_bluecont=fluxo_bluecont[bluecont.values]
  cont=fluxo_bluecont.median()
  espectro['fluxo']=espectro['fluxo']/cont
  return(espectro)

for i in index:
    index_path=spectra_path+'/'+i
    index_dir=os.listdir(index_path)
    ssps=[i for i in index_dir if i.__contains__('SSP')]
    if len(ssps) == 0:
        continue
    else:
      os.chdir(index_path)
      info=[i for i in index_dir if i.__contains__('.txt')]
      info=pd.read_csv(info[0], header=None, names=['name_info', 'info'])
      alfa=pd.read_csv(ssps[0], sep='\s+', skiprows=1, header=None)
      imf=pd.read_csv(ssps[1], sep='\s+', skiprows=1, header=None)
      base=pd.read_csv(ssps[2], sep='\s+', skiprows=1, header=None)
    base=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3], base)
    alfa=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3], alfa)
    imf=normaliza_continuo_simples(info['info'].iloc[2], info['info'].iloc[3], imf)
    plota_spec(base, alfa, imf, info, i)