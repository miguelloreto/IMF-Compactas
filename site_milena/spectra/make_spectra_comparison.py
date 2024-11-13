import matplotlib.pyplot as plt
import numpy as np
import sys
import shutil
import pandas as pd
import os

spectra_path=os.getcwd()
spectra_dir=os.listdir(spectra_path)
index=[i for i in spectra_dir if not i.__contains__('.')]

def plota_spec(base, chem, imf, il, fl, nome_linha):
  basex, basey = base['lambda'], base['flux_norm']
  chemx, chemy = chem['lambda'], chem['flux_norm']
  imfx, imfy = imf['lambda'], imf['flux_norm']
  comp_chem=(chemy-basey)/basey
  comp_imf=(imfy-basey)/basey
  comp=chemy-imfy
  g1=plt.figure(figsize=(8,8))
  plt.title(nome_linha)
  plt.xlabel('$\lambda$($\AA$)')
  plt.ylabel('Diferença absoluta') #($erg/cm^2/s/A$)
  plt.xlim([il,fl])
  plt.ylim(-0.1, 0.1)
  plt.plot(chemx,comp_chem, color='green', label='Abundância',linewidth=0.5)
  plt.plot(imfx, comp_imf ,color='blue', label='IMF',linewidth=0.5)
  plt.plot(imfx,comp, color='orange', label='IMF-Abundância', linewidth=0.5)
  plt.legend()
  plt.savefig(nome_linha+'.png')

for i in index:
    index_path=spectra_path+'/'+i
    index_dir=os.listdir(index_path)
    norms=[i for i in index_dir if i.__contains__('.out')]
    if len(norms) == 0:
        continue
    else:
      os.chdir(index_path)
      info=[i for i in index_dir if i.__contains__('.txt')]
      info=pd.read_csv(info[0], header=None, names=['name_info', 'info'])
      alfa=pd.read_csv(norms[0], sep=' ')
      base=pd.read_csv(norms[1], sep=' ')
      imf=pd.read_csv(norms[2], sep=' ')
    plota_spec(base, alfa, imf, info['info'].iloc[0], info['info'].iloc[1], i)