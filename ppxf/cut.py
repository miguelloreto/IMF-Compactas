import os
import numpy as np
import pandas as pd

list_spectra=os.listdir(os.getcwd())
list_spectra.sort()
outra_lista=set(list_spectra)
conjunto_rep={'spec-0530-52026-0278.fits.cxt','spec-0544-52201-0278.fits.cxt','spec-0656-52148-0523.fits.cxt','spec-0848-52669-0279.fits.cxt','spec-1417-53141-0522.fits.cxt','spec-1671-53446-0522.fits.cxt','spec-1694-53472-0278.fits.cxt', 'spec-2590-54175-0250.fits.cxt', 'spec-1237-52762-0249.fits.cxt','spec-1306-52996-0250.fits.cxt', 'spec-1368-53084-0249.fits.cxt', 'spec-1399-53172-0249.fits.cxt', 'spec-1646-53498-0175.fits.cxt'}

list_spectra=list(outra_lista-conjunto_rep)
list_spectra.sort()

spectra=[]
for i in list_spectra:
    if(i[0] == 's'): spectra.append(pd.read_csv(i, header=0, sep='\s+'))

j=1
for i in spectra:
    i=i.set_axis(['lam', 'flux', 'a', 'b'], axis='columns', copy=False)
    if(i.loc[6400]['flux'] == 0):print(list_spectra[j])
    j=j+1 
