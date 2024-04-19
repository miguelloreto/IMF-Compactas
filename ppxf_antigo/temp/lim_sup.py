import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
import os
import math as m 

"Programa que dado os dataframes de cada IMF, ajusta um polinômio para cada galáxia e acha o limite superior de IMF através da interseção do polinômio com a reta da massa dinâmica"
G=4.3009/1e6

def Kq (q):
    #Função da literatura para corrigir a massa dinâmica devido a geometria da galáxia
    return (0.87 +0.38*np.exp(-3.78*(1-q)))**2
def Kn (n):
    #Função da literatura para corrigir a massa dinâmica devido a geometria da galáxia
    return 8.87 - 0.831*n + 0.0241*n*n
def resolve_eq_seg_ordem(abc):
    #Resolve a equação do segundo grau e retorna um número positivo.
    xmais=(-abc[1]+np.sqrt(abc[1]**2-4*abc[0]*abc[2]))/(2*abc[0])
    xmenos=(-abc[1]-np.sqrt(abc[1]**2-4*abc[0]*abc[2]))/(2*abc[0])
    if (xmais < 0):
        return xmenos.real
    else:
        return xmais.real

#Lendo arquivos e montando o dataframe
os.chdir(os.getcwd()+ "/ppxf/temp")
direstelarmasses=list(os.listdir(os.getcwd()))
tabela_mestelar=pd.DataFrame()
for i in direstelarmasses:
    if (i[0:3]=='gal'):
        massa_imf=pd.read_csv(i)
        index_coluna=massa_imf.columns[0]
        tabela_mestelar.insert(0,column=index_coluna,value=massa_imf[index_coluna])
    if (i=='nomes_amostra_total.csv'):
        nomes=pd.read_csv(i)
    if (i=='info_MCGs_Mdyn.csv'):
        tabela_mdyn=pd.read_csv(i)
#Dando sort no dataframe para que slopes fiquem em ordem crescente e colocando a coluna nomes no final. 
#Reservando um vetor x que tem os slopes.
tabela_mestelar=tabela_mestelar.reindex(sorted(tabela_mestelar.columns), axis=1)
x=[float(i) for i in tabela_mestelar]
#index_coluna_nomes=nomes.columns[0]
#tabela_mestelar_N=tabela_mestelar.insert(0,column=index_coluna_nomes,value=nomes[index_coluna_nomes])
#Calculando a massa dinâmica
tabela_mdyn['Kq']=Kq(-1*(tabela_mdyn['e']-1))
tabela_mdyn['Kn']=Kn(tabela_mdyn['ng'])
tabela_mdyn['Mdyn']=np.log10((tabela_mdyn['Kn']*tabela_mdyn['Kq']*tabela_mdyn['Rhlr']*(tabela_mdyn['sigma_e'])**2)/G)
#Removendo as duas galáxias que não calculei IMF:
indexs_dropados=list()
indexs_dropados.append(tabela_mdyn[(tabela_mdyn.aid =="spec-0303-51615-0278")].index)
indexs_dropados.append(tabela_mdyn[(tabela_mdyn.aid =="spec-1949-53433-0280")].index)
for i in indexs_dropados:
    tabela_mdyn.drop(i[0], inplace=True)


massas_dinamicas=np.array(tabela_mdyn["Mdyn"])
limsup=[]
for i in range(0, len(tabela_mestelar)):
    y=np.array(tabela_mestelar.loc[i])
    abc=np.polyfit(x, y, deg=2)
    abc[2]=abc[2]-massas_dinamicas[i]
    limsup.append(resolve_eq_seg_ordem(abc))
tabela_mdyn.insert(6, "lim_sup", limsup)

plt.figure(figsize=(4,4))
plt.scatter(np.array(tabela_mdyn["log_M"]), np.array(tabela_mdyn["lim_sup"]), color='xkcd:azure', s=0.9)
plt.title("Limite Superior IMF x $M_{*}$")
plt.ylabel("Limite Superior IMF")
plt.xlabel("$\log(M_{*})$ ($M_{\odot}$)")
plt.savefig('plot_comparativo_LIMSUPXMESTELAR.jpg', dpi=900)

plt.figure(figsize=(4,4))
plt.scatter(np.array(tabela_mdyn["sigma_e"]), np.array(tabela_mdyn["lim_sup"]), color='xkcd:azure', s=0.9)
plt.title("Limite Superior IMF x $\sigma$")
plt.ylabel("Limite Superior IMF")
plt.xlabel("$\sigma$ (km/s)")
plt.savefig('plot_comparativo_LIMSUPXDISPVELOCIDADE.jpg', dpi=900)