import matplotlib.pyplot as plt
import numpy as np

def plota_spec(base, chem, imf, il, fl, linha):
  basex, basey = np.loadtxt(base, unpack=True)
  chemx, chemy = np.loadtxt(chem, unpack=True)
  imfx, imfy = np.loadtxt(imf, unpack=True)
  g1=plt.figure(figsize=(8,8))
  plt.title(linha)
  plt.xlabel('$\lambda$($\AA$)')
  plt.ylabel('Fluxo ($erg/cm^2/s/A$)')
  plt.xlim([il,fl])
  plt.plot(basex, basey, color='red', label='Base')
  plt.plot(chemx,chemy, color='green', label='Abundância')
  plt.plot(imfx, imfy ,color='blue', label='IMF')
  plt.legend()
  return g1

print(plota_spec('tio1b.txt', 'tio1c.txt', 'tio1i.txt', 5930, 6000, "TiO1"))
print(plota_spec('tio2b.txt', 'tio2c.txt', 'tio2i.txt', 6180, 6275, "TiO2" ))
print(plota_spec('nab.txt', 'nac.txt', 'nai.txt', 8175, 8205, "Na8190"))
print(plota_spec('cab.txt', 'cac.txt', 'cai.txt', 8480, 8520, "Ca1"))