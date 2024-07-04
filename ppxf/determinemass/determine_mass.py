import pandas as pd
import os
import numpy as np
import sys
import shutil

z_solar=0.019
c = 299792

#Qual foi o slope utilizado ao construir as bases?
slope_esperado=float(sys.argv[1])

#Qual foi o regul utilizado no ppxf?
regul=100

#Qual foi a amostra utilizda?
amostra = 'csgs'

cortadas_mcgs=['spec-0303-51615-0278']
#Quais galáxias tem seu espectro cortados (output do ppxf vai ser diferente):
cortadas=['spec-0544-52201-0278','spec-0656-52148-0523', 'spec-0848-52669-0279', 'spec-1417-53141-0522', 'spec-1671-53446-0522', 'spec-1694-53472-0278']

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

def calcula_cosmodist(z, h_0, omega_l, omega_m, type_z, dz=1e-6):
    '''Calcula a distancia DL para um dado z, retorna um float ou uma lista de distâncias luminosidade'''
    DH_Mpc = c / h_0     # Mpc
    #tH_Gyr = 978 / h_0   # Gyr
    if type_z=='float':
        E_z=np.sqrt(omega_m * (1 +np.arange(0,z,dz))**3 + omega_l)
        DC_Mpc = DH_Mpc * sum(dz/E_z)
        DM_Mpc = DC_Mpc       # Omega_z = 0 (flat Universe)
        DL_Mpc = (1 + z) * DM_Mpc
        return DL_Mpc

    elif type_z=='list':
        DL_Mpc_list=[]
        for i in z:
            E_z=np.sqrt(omega_m * (1 +np.arange(0,i,dz))**3 + omega_l)
            DC_Mpc = DH_Mpc * sum(dz/E_z)
            DM_Mpc = DC_Mpc       # Omega_z = 0 (flat Universe)
            DL_Mpc = (1 + i) * DM_Mpc
            DL_Mpc_list.append(DL_Mpc)
        return np.array(DL_Mpc_list)
    else:
        raise TypeError("Z deve ser uma lista ou um float!")

#Lendo arquivos:
ppxf_path=os.getcwd()
determinemass_path=ppxf_path + '/determinemass'
output_path=ppxf_path + '/output'
amostra_path= output_path + '/' + amostra
spectra_path=ppxf_path + '/spectra'
os.chdir(determinemass_path) 
tabela_mags=pd.read_csv("table_"+amostra+"_cassjobs_mags.csv")
os.chdir(ppxf_path) 
tabela_bases=open('BaseGM_LCGs', 'r')

# Assumindo que no BasesGM todos tem o mesmo slope e todos estao no formato do MILES (assumir isso é ok porque eu que criei o arquivo):
for i in tabela_bases:
    i=i.split()
    slope=float(i[0][3:7])
    if(slope!=slope_esperado):
        raise('Os slopes da base e do output não coincidem!')

os.chdir(spectra_path)

# Renormalizar as massas, pegando
l_norm_min = 5590      #Re-normalizando Mini and Mcor
l_norm_max = 5680    
fator_norm_spectra=[]
espectros=os.listdir(spectra_path)
espectros.sort()
for i in espectros:
    if i[0]!='c':
        espectro=pd.read_csv(i, delimiter='\s+', index_col=False)
        espectro.columns=['lam', 'flux', 'sn', '?']
        espectro=espectro[(espectro.lam >l_norm_min) & (espectro.lam <l_norm_max)]
        fator_norm_spectra.append(espectro.flux.median())
    else:
        continue
        
os.chdir(output_path)
lista_df_espectros=[]
nome_espectros=[]
mcor=[]
mini=[]
fator_norm_ppxf=[]

outputs=os.listdir(output_path)
outputs.sort()
for i in outputs:
    i=i.split(".")
    if (i[-1] == 'out'):
        if(i[0] in cortadas):
            selecionado='.'.join(i)
            spec=open(output_path+'/'+ selecionado, 'r')
            linhas=spec.readlines()
            mini.append(float(linhas[10].split(" ")[-1]))
            mcor.append(float(linhas[11].split(" ")[-1]))
            fator_norm_ppxf.append(float(linhas[8].split(" ")[-1]))
            df=pd.read_csv(selecionado, skiprows=34,names=['x_j', 'M_ini', 'Mcor', 'age', 'met', 'L/M', 'Mstar'], nrows=112, sep='\s+', usecols=[1,2,3,4,5,6,7])
            lista_df_espectros.append(df)
            nome_espectros.append(selecionado)
        else:
            selecionado='.'.join(i)
            spec=open(output_path+'/'+ selecionado, 'r')
            linhas=spec.readlines()
            mini.append(float(linhas[10].split(" ")[-1]))
            mcor.append(float(linhas[11].split(" ")[-1]))
            fator_norm_ppxf.append(float(linhas[8].split(" ")[-1]))
            df=pd.read_csv(selecionado, skiprows=40,names=['x_j', 'M_ini', 'Mcor', 'age', 'met', 'L/M', 'Mstar'], nrows=112, sep='\s+', usecols=[1,2,3,4,5,6,7])
            lista_df_espectros.append(df)
            nome_espectros.append(selecionado)
    else:
        continue

if(len(fator_norm_ppxf) > 5574): raise Exception('sei la cara')

#Movendo os arquivos para bons lugares:
for i in outputs:
    i=i.split("-")
    if(i[0]=='spec'):shutil.move("-".join(i), amostra_path)

os.chdir(amostra_path)
outputs_amostra=os.listdir(amostra_path)
outputs_amostra.sort()

check_dir=False
for j in outputs_amostra: 
    if(j=='output_regul'+str(regul)+'_total-'+str(slope)): check_dir=True
if(check_dir==False): os.mkdir("output_regul"+str(regul)+"_total-" + str(slope))

for k in outputs_amostra:
    k=k.split("-")
    if(k[0]=='spec'):shutil.move("-".join(k), amostra_path+'/output_regul'+str(regul)+'_total-'+str(slope))

#Re-Normalizando:
fator_norm=np.array(fator_norm_spectra)/np.array(fator_norm_ppxf)
mini=np.array(mini)*fator_norm
mcor=np.array(mcor)*fator_norm

tabela_massas=pd.DataFrame(index=[i for i in range(0,len(lista_df_espectros))], columns=['Nome', 'Slope', 'Idade Média(L)', 'Idade Média Log(L)', 'Metalicidade Média(L)', 'Metalicidade Média Log(L)' 'Idade Média(M)', 'Idade Média Log(M)', 'Metalicidade Média(M)', 'Metalicidade Média Log(M)' , 'Mcor', 'Mini'])


os.chdir(determinemass_path)
#Certamente tem um jeito melhor de fazer esse processo, OH WELL...
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


tabela_massas['Mcor']=mcor
tabela_massas['Mini']=mini
tabela_mags.drop_duplicates(inplace=True)
tabela_massas.drop_duplicates(inplace=True)
tabela_massas=tabela_massas.sort_values(by='Nome')
tabela_mags=tabela_mags.sort_values(by=['plate', 'mjd', 'fiberid'])
# O Pandas gosta de multiplicar usando os indices, truque interessante pra reseta-los:
tabela_massas = pd.DataFrame(data=tabela_massas.to_numpy(), columns=tabela_massas.columns)
tabela_mags = pd.DataFrame(data=tabela_mags.to_numpy(), columns=tabela_mags.columns)

# Uniformizando os nomes nas tabelas
splitado_ponto=tabela_massas.Nome.str.split('.')
primeiro_nome=[]
for i in splitado_ponto:
    primeiro_nome.append(i[0])
nomes_uniforme=[]
for i in primeiro_nome:
    splitado_hifen=i.split('-')
    splitado_hifen.pop(0)
    nome=[j.lstrip('0') for j in splitado_hifen]
    nome=[j.lstrip('0') for j in nome]
    nomes_uniforme.append('-'.join(nome))
tabela_massas.Nome=nomes_uniforme

tabela_mags.plate=tabela_mags.plate.astype(int)
tabela_mags.mjd=tabela_mags.mjd.astype(int)
tabela_mags.fiberid=tabela_mags.fiberid.astype(int)
tabela_mags.plate=tabela_mags.plate=tabela_mags.plate.astype(str)
tabela_mags.mjd=tabela_mags.mjd.astype(str)
tabela_mags.fiberid=tabela_mags.fiberid.astype(str)
nome=tabela_mags.plate.str.cat(tabela_mags.mjd, sep='-', join='right')
nome=nome.str.cat(tabela_mags.fiberid, sep='-', join='right')
tabela_mags.drop(labels=['plate', 'mjd', 'fiberid'], axis=1, inplace=True)
tabela_mags.insert(0, "Nome", nome)

#Fazendo todos cálculos:
tabela_mags['fator_z']=1e-17*4*np.pi*((calcula_cosmodist(tabela_mags.z, 70, 0.7, 0.3, type_z='list')*3.086e24)**2)/(3.826e33)
tabela_mags['fator_mag']=10**((-0.4)*(tabela_mags['fibermag_r']-tabela_mags['petromag_r']))
tabela_massas['Mcor']=tabela_massas['Mcor']*tabela_mags['fator_z']/tabela_mags['fator_mag']
tabela_massas['Mini']=tabela_massas['Mini']*tabela_mags['fator_z']/tabela_mags['fator_mag']
tabela_massas['Mini']=tabela_massas['Mini'].astype(float)
tabela_massas['Mini']=np.log10(tabela_massas['Mini'])
tabela_massas['Mcor']=tabela_massas['Mcor'].astype(float)
tabela_massas['Mcor']=np.log10(tabela_massas['Mcor'])


#Escrevendo arquivos na pasta dos gráficos:
os.chdir(ppxf_path + '/graphs/todas_galaxias/'+ amostra + '/regul' + str(regul))
tabela_massas.to_csv('nomes_amostra_total.csv', index=False, columns=['Nome'], header='Nomes')
tabela_massas.to_csv('massas_{}_amostra_total.csv'.format(slope), index=False, columns=['Nome',"Mcor"], header=['Nome',str(slope)])
#if (slope==1.3): tabela_massas.to_csv('comparar_salim.csv'.format(slope),index=False, columns=['Nome', 'Mcor'], header=['Amostra', 'Mcor'])
#Opção para comparar com as massas do salim