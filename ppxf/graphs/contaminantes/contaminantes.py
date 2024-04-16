import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
plt.rcParams.update({'figure.max_open_warning': 0})


tabela_contaminantes=pd.read_csv('contaminates.csv')
tabela_mags=pd.read_csv("myTable1_out_miguel_loreto.csv")

#Procedimento padrão de ler tabelas
tabela_mags.drop_duplicates(inplace=True)
tabela_contaminantes.drop_duplicates(inplace=True)
tabela_mags=tabela_mags.sort_values(by=['plate', 'mjd', 'fiberid'])
tabela_contaminantes=tabela_contaminantes.sort_values(by=['aid'])
tabela_mags = pd.DataFrame(data=tabela_mags.to_numpy(), columns=tabela_mags.columns)
tabela_contaminantes = pd.DataFrame(data=tabela_contaminantes.to_numpy(), columns=tabela_contaminantes.columns)

#Problema dos nomes, de novo:
nomes_uniforme=[]
#Conservando o nome antigo!
tabela_nomes_spec=tabela_contaminantes[['aid', 'ra', 'dec']]
tabela_nomes_spec=tabela_nomes_spec.sort_values(by=['ra', 'dec'])
for i in tabela_contaminantes['aid']:
    splitado_hifen=i.split('-')
    splitado_hifen.pop(0)
    nome=[j.lstrip('0') for j in splitado_hifen]
    nome=[j.lstrip('0') for j in nome]
    nomes_uniforme.append('-'.join(nome))
tabela_contaminantes.aid=nomes_uniforme

tabela_mags.plate=tabela_mags.plate.astype(int)
tabela_mags.mjd=tabela_mags.mjd.astype(int)
tabela_mags.fiberid=tabela_mags.fiberid.astype(int)
tabela_mags.plate=tabela_mags.plate=tabela_mags.plate.astype(str)
tabela_mags.mjd=tabela_mags.mjd.astype(str)
tabela_mags.fiberid=tabela_mags.fiberid.astype(str)
nome=tabela_mags.plate.str.cat(tabela_mags.mjd, sep='-', join='right')
nome=nome.str.cat(tabela_mags.fiberid, sep='-', join='right')
tabela_mags.drop(labels=['plate', 'mjd', 'fiberid'], axis=1, inplace=True)
tabela_mags.insert(0, "aid", nome)

#Enxugando a tabela de contaminantes:
tabela_contaminantes.drop(['e', 'ng', 'Kq', 'Kn', 'lim_sup', 'z'], axis=1, inplace=True)
tabela_contaminantes.rename(columns={'Mdyn':'log_M_dyn'}, inplace=True)
tabela_contaminantes=tabela_contaminantes[['aid','ra','dec','Rhlr','log_M', 'log_M_dyn', 'sigma_e']]

#Coincidindo as tabelas:
tabela_mags=tabela_mags[tabela_mags.aid.isin(tabela_contaminantes.aid)]
#CONCAT EH SEMPRE PRA ADICIONAR MAIS DADOS EM UMA MESMA COLUNA:
#tabela_contaminantes=pd.concat([tabela_contaminantes,tabela_mags], axis=0, ignore_index=True)
# JOIN, adiciona um dataframe no outro
tabela_contaminantes=tabela_contaminantes.join(tabela_mags.set_index('aid'), on='aid')
tabela_contaminantes['fator_mag']=10**((-0.4)*(tabela_contaminantes['fibermag_r']-tabela_contaminantes['petromag_r']))


tabela_contaminantes.to_csv("coord_contaminantes.csv", columns=['ra', 'dec'], header=None, index=False,sep=',')
#Via Vizier, de um jeito bem idiota, serio, vtnc, consigo essa escala:
tabela_scale=pd.read_csv("scales.tsv", header=None, sep=';',skiprows=38)
tabela_scale.columns=['z', 'scale', 'ra', 'dec']
tabela_scale.drop_duplicates(inplace=True)
tabela_contaminantes=tabela_contaminantes.sort_values(by=['ra', 'dec'])
tabela_scale=tabela_scale.sort_values(by=['ra', 'dec'])
tabela_contaminantes.insert(len(tabela_contaminantes.columns), "scale", np.array(tabela_scale.scale))
tabela_contaminantes = pd.DataFrame(data=tabela_contaminantes.to_numpy(), columns=tabela_contaminantes.columns)

r_petro=[1.76, 1.99, 2.56, 2.18, 30.69,2.15,2.50,2.28, 2.68, 2.32,2.42,2.27,80,2.31, 7.85,2.35,1.58,1.66,2.78,2,1.81,1.79,2.83,5.12,32.52,2.14,5.92,1.73,1.64,2.31,2.28,1.9,20.69,1.99,1.91,1.89]
tabela_contaminantes.insert(len(tabela_contaminantes.columns), "r_petro", np.array(r_petro))
tabela_contaminantes.r_petro=tabela_contaminantes.r_petro*tabela_contaminantes.scale
tabela_scale=tabela_scale.sort_values(by=['ra', 'dec'])
print(tabela_contaminantes)
print(tabela_nomes_spec)
tabela_contaminantes.aid=np.array(tabela_nomes_spec.aid)
print(tabela_contaminantes)
tabela_contaminantes=tabela_contaminantes.sort_values(by=['aid'])
tabela_contaminantes.drop(['z'],inplace=True, axis=1)
tabela_contaminantes = pd.DataFrame(data=tabela_contaminantes.to_numpy(), columns=tabela_contaminantes.columns)
print(tabela_contaminantes)

#Imprimindo o .pdf:
tabela_latex=tabela_contaminantes.to_latex(float_format="{:.2f}".format,index=False)
tabela_latex=tabela_latex.replace("\\\n", "\\ \hline\n")
tabela_latex=tabela_latex.replace("\\toprule", "\hline")
tabela_latex=tabela_latex.replace("\\midrule", "")
tabela_latex=tabela_latex.replace("\\bottomrule", "")
print(tabela_latex)
dir_contaminantes=os.getcwd()
dir_ppxf=dir_contaminantes +'/output_ppxf_ch'
dir_images=dir_contaminantes +'/images'
dir_fits=dir_contaminantes +'/fits'
list_ppxf=os.listdir(dir_ppxf)
list_images=os.listdir(dir_images)
list_fits=os.listdir(dir_fits)
list_ppxf.sort()
list_images.sort()
list_fits.sort()

for index, row in tabela_contaminantes.iterrows():
    print(r'\newpage')
    print(r'\begin{figure}')
    print('\centering')
    print(r'\begin{tabular}{cc}')
    os.chdir(dir_ppxf)
    print(r'\multicolumn{2}{c}{\includegraphics[width=0.7\hsize]{output_ppxf_ch/'+list_ppxf[index] +r'}} \\')
    os.chdir(dir_fits)
    print(r'\includegraphics[width=0.45\hsize]{fits/'+list_fits[index] +r'}  & ')
    os.chdir(dir_images)
    print(r'\includegraphics[width=0.45\hsize]{images/'+list_images[index] +r'} \\')
    print('\end{tabular}')
    print('\end{figure}')
    linha_latex=row.to_frame().T.to_latex(header=['Aid', 'Ra', 'Dec', 'R_hlr', 'mstar', 'mdyn', 'sig', 'fiber', 'petro','fator','scale','rpetro'], index=False, float_format="{:.2f}".format, column_format='|c|c|c|c|c|c|c|c|c|c|c|c|c|')
    linha_latex=linha_latex.replace("\\\n", "\\ \hline\n")
    linha_latex=linha_latex.replace("\\toprule", "\hline")
    linha_latex=linha_latex.replace("\\midrule", "")
    linha_latex=linha_latex.replace("\\bottomrule", "")
    print(linha_latex)

#Exemplo de printar espectros, caso venha a ser útil de novo:
'''
dir_contaminantes=os.getcwd()
dir_spectra=dir_contaminantes +'/spectra'
os.chdir(dir_spectra)
espectros=os.listdir(dir_spectra)

for i in espectros:
    espectro=np.loadtxt(i)
    espectro=pd.DataFrame(espectro)
    espectro.columns=['lam', 'fluxo', 'erro', 'erro2']
    espectro= espectro[espectro.fluxo != 0]
    print(espectro)
    plt.figure(figsize=(6,6))
    plt.plot(np.array(espectro.lam),np.array(espectro.fluxo), color='xkcd:red', label='Espectro')
    plt.xlabel('$\lambda$($\AA$)')
    plt.ylabel('Fluxo ($erg/cm^2/s/A$)')
    plt.legend()
    plt.savefig('specter_'+ str(i) +'.jpg', dpi=900)

'''