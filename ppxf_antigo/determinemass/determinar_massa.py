import pandas as pd
import os
import numpy as np

z_solar=0.019

def calcula_media_ponderada(lista, peso):
    '''Calcula a media ponderda de uma lista[N] sendo que cada elemento tem peso peso[N]'''
    if(len(lista)!=len(peso)):
        raise Exception()
    media= (peso*lista).sum()/peso.sum()
    return media

def calcula_media_ponderada_log(lista, peso):
    '''Calcula a media ponderda de uma log(lista[N]) sendo que cada elemento tem peso: peso[N]'''
    import numpy as np
    if(len(lista)!=len(peso)):
        raise Exception()
    media= (peso*np.log10(lista)).sum()/peso.sum()
    return media

def calcula_cosmodist(z, h_0, omega_l, omega_m, dz=1e-6):
    '''Calcula a distancia DL para um dado z'''
    c = 299792
    DH_Mpc = c / h_0     # Mpc
    tH_Gyr = 978 / h_0   # Gyr
    
    
    E_z=np.sqrt(omega_m * (1 +np.arange(0,z,dz))**3 + omega_l)
    DC_Mpc = DH_Mpc * sum(dz/E_z)
    DM_Mpc = DC_Mpc       # Omega_z = 0 (flat Universe)
    DL_Mpc = (1 + z) * DM_Mpc

    return DL_Mpc

cwd_path=os.getcwd()
ppxf_path= cwd_path + '/ppxf'
determinemass_path= ppxf_path + '/determinemass'
output_path=ppxf_path + '/output'
os.chdir(determinemass_path) 
tabela_propriedades=pd.read_csv("table_compact_stellar_population_properties.csv")
tabela_mags=pd.read_csv("myTable_out_miguel_loreto.csv")
os.chdir(ppxf_path) 
tabela_bases=open('BaseGM_LCGs', 'r')
#Assumindo que no BasesGM todos tem o mesmo slope e todos estao no formato do MILES (assumir isso é ok porque eu que criei o arquivo):
slope_esperado=0.8
slope=0.8
for i in tabela_bases:
    i=i.split()
    if(slope==slope_esperado):
        slope=float(i[0][3:7])
    else:
        print('Os slopes da base e do output não coincidem!')
os.chdir(cwd_path)


os.chdir(output_path)
lista_df_espectros=[]
nome_espectros=[]
mcor=[]
mini=[]
fator_norm=[]
for i in os.listdir(output_path):
    i=i.split(".")
    if (i[-1] == 'out'):
        i='.'.join(i)
        spec=open(output_path+'/'+ i, 'r')
        linhas=spec.readlines()
        mini.append(float(linhas[10].split(" ")[-1]))
        mcor.append(float(linhas[11].split(" ")[-1]))
        fator_norm.append(float(linhas[8].split(" ")[-1]))
        df=pd.read_csv(i, index_col=0, skiprows=40,names=['j','x_j', 'M_ini', 'Mcor', 'age', 'met', 'L/M', 'Mstar'], nrows=112, delim_whitespace=True, usecols=[0,1,2,3,4,5,6,7])
        lista_df_espectros.append(df)
        nome_espectros.append(i)

tabela_massas=pd.DataFrame(index=[i for i in range(0,len(lista_df_espectros))], columns=['Nome', 'Slope', 'Idade Média(L)', 'Idade Média Log(L)', 'Metalicidade Média(L)', 'Metalicidade Média Log(L)' 'Idade Média(M)', 'Idade Média Log(M)', 'Metalicidade Média(M)', 'Metalicidade Média Log(M)' , 'Mcor', 'Mini'])

os.chdir(cwd_path)


idade_media_luminosidade=[]
idade_media_log_luminosidade=[]
metalicidade_media_luminosidade=[]
metalicidade_media_log_luminosidade=[]
idade_media_massa=[]
idade_media_log_massa=[]
metalicidade_media_massa=[]
metalicidade_media_log_massa=[]

for i in lista_df_espectros:
    idade_media_luminosidade.append(calcula_media_ponderada(i['age'], i['x_j']))
    idade_media_log_luminosidade.append(calcula_media_ponderada_log(i['age'], i['x_j']))
    metalicidade_media_luminosidade.append(calcula_media_ponderada(i['met']/z_solar, i['x_j']))
    metalicidade_media_log_luminosidade.append(calcula_media_ponderada_log(i['met']/z_solar, i['x_j']))
    idade_media_massa.append(calcula_media_ponderada(i['age'], i['Mcor']))
    idade_media_log_massa.append(calcula_media_ponderada_log(i['age'], i['Mcor']))
    metalicidade_media_massa.append(calcula_media_ponderada(i['met']/z_solar, i['Mcor']))
    metalicidade_media_log_massa.append(calcula_media_ponderada_log(i['met']/z_solar, i['Mcor']))

#Fator norm:
mcor=np.array(mcor)
mini=np.array(mini)
fator_norm=np.array(fator_norm)
mini=mini/fator_norm
mcor=mcor/fator_norm

for i in range(0, len(tabela_massas)):
    tabela_massas.loc[i]['Nome']=nome_espectros[i]
    tabela_massas.loc[i]['Slope']=slope
    tabela_massas.loc[i]['Idade Média(L)']=idade_media_luminosidade[i]/1e9
    tabela_massas.loc[i]['Idade Média Log(L)']=(10**idade_media_log_luminosidade[i])/1e9
    tabela_massas.loc[i]['Metalicidade Média(L)']=np.log10(metalicidade_media_luminosidade[i])
    tabela_massas.loc[i]['Metalicidade Média Log(L)']=metalicidade_media_log_massa[i]
    tabela_massas.loc[i]['Idade Média(M)']=idade_media_massa[i]/1e9
    tabela_massas.loc[i]['Idade Média Log(M)']=(10**idade_media_log_massa[i])/1e9
    tabela_massas.loc[i]['Metalicidade Média(M)']=np.log10(metalicidade_media_massa[i])
    tabela_massas.loc[i]['Metalicidade Média Log(M)']=metalicidade_media_log_massa[i]
    tabela_massas.loc[i]['Mcor']=mcor[i]
    tabela_massas.loc[i]['Mini']=mini[i]

tabela_massas=tabela_massas.sort_values(by='Nome')
tabela_mags=tabela_mags.sort_values(['plate', 'mjd', 'fiberid'])

fator_z=[]
for i in tabela_mags.z:
    fator_z.append(1e-17*4*np.pi*((calcula_cosmodist(i, 70, 0.7, 0.3)*3.086e24)**2)/(3.826e33))

tabela_mags.insert(6, "fator_z", fator_z)
tabela_mags['fator_mag']=10**((-0.4)*(tabela_mags['fibermag_r']-tabela_mags['petromag_r']))

tabela_massas['Mcor']=tabela_massas['Mcor']*tabela_mags['fator_z']/tabela_mags['fator_mag']
tabela_massas['Mini']=tabela_massas['Mini']*tabela_mags['fator_z']/tabela_mags['fator_mag']

print(tabela_mags)
print(tabela_massas)