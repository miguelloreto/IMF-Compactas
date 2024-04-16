import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np
import os

G=4.3009/1e6
cores=['xkcd:cobalt', 'xkcd:pale orange', 'xkcd:kelly green', 'xkcd:deep red', 'xkcd:sea', 'xkcd:medium purple', 'xkcd:raw sienna', 'xkcd:bubblegum pink','xkcd:grey', 'xkcd:sand yellow']

def Kq (q):
    return (0.87 +0.38*np.exp(-3.78*(1-q)))**2
def Kn (n):
    return 8.87 - 0.831*n + 0.0241*n*n

slopes=[1.3, 1.8, 2.3, 3.3]
direstelarmasses=list(os.listdir(os.getcwd()))
tabelas=[]

for i in direstelarmasses:
    if (i[0:3]=='gal' and i[15]!='c'):
        tabelas.append(i)
    if (i=='nomes.csv'):
        nomes=pd.read_csv(i)
    if (i=='info_MCGs_Mdyn.csv'):
        tabela_mdyn=pd.read_csv(i)
tabelas.sort()

tabela_slope_massa=pd.DataFrame()
cont=0
for i in tabelas:
    i=pd.read_csv(i)
    tabela_slope_massa.insert(cont, str(slopes[cont]), i, allow_duplicates=True)
    cont+=1

nomes=list(nomes['Nome'])
for i in range(0, len(nomes)):
    nomes[i]=nomes[i].split('.')
    nomes[i]=nomes[i].pop(0)
    

tabela_slope_massa.index=nomes

tabela_mdyn['Kq']=Kq(-1*(tabela_mdyn['e']-1))
tabela_mdyn['Kn']=Kn(tabela_mdyn['ng'])
tabela_mdyn['Mdyn']=(tabela_mdyn['Kn']*tabela_mdyn['Kq']*tabela_mdyn['Rhlr']*(tabela_mdyn['sigma_e'])**2)/G
tabela_mdyn['Mdyn']=np.log10(tabela_mdyn['Mdyn'])
tabela_mdyn=tabela_mdyn.head(10)
tabela_mdyn.index=tabela_mdyn['aid']
tabela_mdyn.drop(columns='aid', axis=1, inplace=True)
tabela_mdyn.index.name=None

nomes_provisorio=nomes.copy()
'''
for i in nomes:
    if(tabela_slope_massa.loc[i]['1.3'] < tabela_mdyn.loc[i]['Mdyn']):
        tabela_slope_massa.drop(labels=i, axis=0, inplace=True)
        tabela_mdyn.drop(labels=i, axis=0, inplace=True)
        nomes_provisorio.remove(i)
        cores.pop()

nomes=nomes_provisorio
'''
tabela_slope_massa=tabela_slope_massa.transpose()
tabela_slope_massa.insert(0, 'Slopes', slopes)
tabela_slope_massa.reset_index(drop=True, inplace=True)


plt.figure()
ax1=tabela_slope_massa.plot(kind='scatter', x='Slopes', y=nomes[0], title='Gráfico corrigido', color=cores[0], label=nomes[0])
for i in range(1, len(nomes)):
    ax2=tabela_slope_massa.plot(kind='scatter', x='Slopes', y=nomes[i], title='Gráfico corrigido', color=cores[i],label=nomes[i], ax=ax1)

j=0
for i in np.array(tabela_mdyn['Mdyn']):
    plt.axhline(i, color=cores[j])
    j=j+1
plt.xticks(slopes)
plt.xlabel("Slopes")
plt.ylabel("Log de Massa Estelar (MCor)")
plt.legend(prop = { "size": 6}, loc ="lower right")
plt.savefig('g1_final.jpg', dpi=400)

