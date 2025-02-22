import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

"Programa que dado os dataframes de cada IMF, ajusta um polinômio para cada galáxia e acha o limite superior de IMF através da interseção do polinômio com a reta da massa dinâmica"
G=4.3009/1e6
P=0.006

def Kq (q):
    #Função da literatura para corrigir a massa dinâmica devido a geometria da galáxia
    return (0.87 +0.38*np.exp(-3.78*(1-q)))**2
def Kn (n):
    #Função da literatura para corrigir a massa dinâmica devido a geometria da galáxia
    return 8.87 - 0.831*n + 0.0241*n*n
def escolhe_raiz(raizes):
    #Função que escolhe a raiz positiva de duas soluções da equação do segundo grau em forma de vetor
    if(type(raizes[0])==np.complex128 and type(raizes[1])==np.complex128): return 0
    elif(raizes[1]>raizes[0]):
        return raizes[1]
    elif(raizes[0]>raizes[1]): 
        return raizes[0]
    else:
        # Cair nessa condiçaõ quer dizer que o limite superior é imaginário. Isso que eu chamo de contaminante. 
        return 0
    
#Qual foi o regul utilizado no ppxf?
regul=10

#Qual foi a amostra utilizda?
amostra = 'mcgs'

initial_path=os.getcwd()
amostra_path=initial_path+'/'+amostra
os.chdir(amostra_path)
tabela_mdyn=pd.read_csv("table_"+amostra+"_mdyn.csv")
#Lendo arquivos e montando o dataframe
os.chdir(amostra_path+"/regul"+str(regul))
direstelarmasses=list(os.listdir(os.getcwd()))
direstelarmasses.sort()
tabela_mestelar=pd.DataFrame()
flag=0
for i in direstelarmasses:
    if (i[0:2]=='ma'):
        massa_imf=pd.read_csv(i)
        nomes=massa_imf.Nome
        index_slope=massa_imf.columns[1]
        tabela_mestelar.insert(flag,column=index_slope,value=massa_imf[index_slope])
        flag+=1

#Reservando um vetor x que tem os slopes.
x=[float(i) for i in tabela_mestelar]

#Calculando a massa dinâmica
tabela_mdyn['Kq']=Kq(-1*(tabela_mdyn['e']-1))
tabela_mdyn['Kn']=Kn(tabela_mdyn['ng'])
tabela_mdyn['Mdyn']=np.log10((tabela_mdyn['Kn']*tabela_mdyn['Kq']*tabela_mdyn['Rhlr']*(tabela_mdyn['sigma_e'])**2)/G)

#Dropando as linhas cuja massa dinâmica não pode ser calculada (não tem o Rhlr ou algum outro parâmetro), resetando índices
for i in range(0, len(tabela_mdyn)):
    if(np.isnan(tabela_mdyn['Mdyn'][i]) == True):tabela_mdyn.drop(axis=0,index=i, inplace=True), tabela_mestelar.drop(axis=0,index=i, inplace=True)
tabela_mdyn = pd.DataFrame(data=tabela_mdyn.to_numpy(), columns=tabela_mdyn.columns)
tabela_mestelar = pd.DataFrame(data=tabela_mestelar.to_numpy(), columns=tabela_mestelar.columns)

massas_dinamicas=np.array(tabela_mdyn["Mdyn"])
limsup=[]
flag=0

os.chdir(initial_path)
#Problema: todos são acusados como imaginarios, i.e., lim==0, problema deve ser na função escolhe_raizes HAHAHAHAAHAH
for i in range(0, len(tabela_mdyn)):
    y=np.array(tabela_mestelar.loc[i])
    abc=np.polyfit(x, y, deg=2)
    coef=np.polyfit(x, y, deg=2)
    abc[2]=abc[2]-massas_dinamicas[i]
    raizes=np.roots(abc)
    lim=float(escolhe_raiz(raizes))
    if(i == 4):
        xis=np.arange(0, x[-1], 0.001)
        plt.rcParams.update({'font.size': 15})
        plt.figure(figsize=(8,8))
        plt.scatter(x,np.array(y), color='xkcd:azure', s=35, label= 'Massa Calculada')
        plt.plot(xis,coef[0]*xis**2 + coef[1]*xis + coef[2], color='xkcd:red', label='Fit')
        plt.axhline(massas_dinamicas[i], color='xkcd:purple', label=r'$\log(M_{dyn}/M_\odot)$')
        plt.axvline(lim, color='xkcd:orange', label=r'$\Gamma_{max}$', linestyle='dashed')
        plt.xticks(np.arange(0, 3.6, 0.5), )
        plt.yticks([10.9, 11, 11.1, 11.2, 11.3, 11.4])
        plt.ylabel(r"$\log(M_{\star}/M_\odot)$", fontsize=25)
        plt.xlabel(r"$\Gamma$", fontsize=25)
        plt.legend(prop = { "size":15 })
        plt.grid()
        plt.savefig('fit_example.png', dpi=1200)
    limsup.append(lim)

'''
tabela_mdyn.insert(12, "lim_sup", limsup)
tabela_mestelar.insert(0,column="Nome",value=nomes)
tabela_contaminantes=tabela_mdyn[tabela_mdyn['lim_sup'] <= 0]
tabela_mdyn=tabela_mdyn[tabela_mdyn['lim_sup'] != 0]
tabela_contaminantes.to_csv("contaminates.csv", index=False)
tabela_mdyn.to_csv("tabela_lim_sup.csv", index=False)
'''