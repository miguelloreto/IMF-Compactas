import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
slopes=[0.8, 1.3, 1.8, 2.3, 3.3]
direstelarmasses=list(os.listdir(os.getcwd()))
tabelas=[]
for i in direstelarmasses:
    if (i[0:3]=='gal'):
        tabelas.append(i)
    if (i=='nomes.csv'):
        nomes=pd.read_csv(i)
tabelas.sort()

tabelona=pd.DataFrame()
cont=0
for i in tabelas:
    i=pd.read_csv(i)
    tabelona.insert(cont, str(slopes[cont]), i, allow_duplicates=True)
    cont+=1

nomes=list(nomes['Nome'])
for i in range(0, len(nomes)):
    nomes[i]=nomes[i].split('.')
    nomes[i]=nomes[i].pop(0)

tabelona.index=nomes
tabelona=tabelona.transpose()
print(tabelona)

plt.figure()
tabelona.plot(use_index=True, y=nomes, title='Gráfico para 10 espectros').get_figure()
plt.xlabel("Slopes")
plt.ylabel("Log de Massa Estelar (MCor)")
plt.legend(prop = { "size": 8 }, loc ="lower right")
plt.savefig('g1.jpg', dpi=400)

