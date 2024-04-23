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
        #print(raizes)
        return 0
def cria_bins(dados):
    #Cria o número de bins apropriado para dos dados
    max=dados.max()
    min=dados.min()
    dp=dados.std()
    return int((max-min)/(dp*0.5))

#Lendo arquivos e montando o dataframe
os.chdir(os.getcwd()+'/regul0')
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
    if (i=='info_MCGs_Mdyn.csv'):
        tabela_mdyn=pd.read_csv(i)

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

'#Problema: todos são acusados como imaginarios, i.e., lim==0, problema deve ser na função escolhe_raizes'
for i in range(0, len(tabela_mdyn)):
    y=np.array(tabela_mestelar.loc[i])
    abc=np.polyfit(x, y, deg=2)
    coef=np.polyfit(x, y, deg=2)
    abc[2]=abc[2]-massas_dinamicas[i]
    raizes=np.roots(abc)
    lim=float(escolhe_raiz(raizes))
    '''
    if(lim <= 0):
        xis=np.arange(0, x[-1], 0.001)
        plt.figure(figsize=(6,6))
        plt.scatter(x,np.array(y), color='xkcd:azure', s=6, label= 'Massa Calculada')
        plt.plot(xis,coef[0]*xis**2 + coef[1]*xis + coef[2], color='xkcd:red', label='Fit')
        plt.axhline(massas_dinamicas[i], color='xkcd:purple', label=r'$\log(M_{Dyn})$')
        plt.ylabel(r"$\log(M_\ast)$")
        plt.xlabel("Slope IMF")
        plt.legend()
        plt.grid()
        plt.savefig('fit_example'+ str(tabela_mdyn['aid'][i]) +'.jpg', dpi=500)
    '''
    limsup.append(lim)

tabela_mdyn.insert(12, "lim_sup", limsup)
tabela_mestelar.insert(0,column="Nome",value=nomes)
tabela_contaminantes=tabela_mdyn[tabela_mdyn['lim_sup'] <= 0]
tabela_mdyn=tabela_mdyn[tabela_mdyn['lim_sup'] != 0]
print(tabela_contaminantes)
'''
plt.figure(figsize=(5,5))
plt.scatter(np.array(tabela_mdyn["log_M"]), np.array(tabela_mdyn["lim_sup"]), color='xkcd:azure', s=0.9)
plt.title("Limite Superior IMF x $M_{*}$")
plt.ylabel("Limite Superior IMF")
plt.xlabel("$\log(M_{*})$ ($M_{\odot}$)")
plt.savefig('plot_comparativo_LIMSUPXMESTELAR.jpg', dpi=900)

plt.figure(figsize=(5,5))
plt.scatter(np.array(tabela_mdyn["sigma_e"]), np.array(tabela_mdyn["lim_sup"]), color='xkcd:azure', s=0.9)
plt.title("Limite Superior IMF x $\sigma$")
plt.ylabel("Limite Superior IMF")
plt.xlabel("$\sigma$ (km/s)")
plt.savefig('plot_comparativo_LIMSUPXDISPVELOCIDADE.jpg', dpi=900)

plt.figure(figsize=(5,5))
plt.hist(np.array(tabela_mdyn['lim_sup']),bins=cria_bins(tabela_mdyn['lim_sup']), density=True, color = 'xkcd:lightblue',edgecolor='black', label='Distribuição')
plt.axvline(tabela_mdyn['lim_sup'].mean(), color='xkcd:indigo', ls='dashed', label="Média")
plt.legend()
#plt.plot(array_x_kde_plot, kde(array_x_kde_plot),color='xkcd:blue', label="KDE")
plt.ylabel("Frequência")
plt.xlabel("Limites Superiors")
plt.savefig('histograma_limsup.jpg', dpi=900)

tabela_contaminantes.to_csv("contaminates.csv", index=False)
'''