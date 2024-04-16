import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats  

#Note que aqui o slope tem que ser o 1.3 para fazer a comparação com o Salim.
tabela_massas_calculadas=pd.read_csv("comparar_salim.csv")
tabela_info_MCGs=pd.read_csv("info_MCGs_Mdyn.csv")
tabela_info_MCGs.drop_duplicates(inplace=True)
tabela_massas_calculadas["M_Salim"]=tabela_info_MCGs["log_M"]
tabela_massas_calculadas["Diferença_Salim"]=tabela_massas_calculadas["Mcor"]-tabela_massas_calculadas["M_Salim"]

difmax=tabela_massas_calculadas["Diferença_Salim"].max()
difmin=tabela_massas_calculadas["Diferença_Salim"].min()
desvio_padrao_dif_salim=tabela_massas_calculadas["Diferença_Salim"].std()
print(difmax, difmin, desvio_padrao_dif_salim)
bins=int((difmax-difmin)/(desvio_padrao_dif_salim*0.5))
array_mcor=np.array(tabela_massas_calculadas["Mcor"])
array_salim=np.array(tabela_massas_calculadas["M_Salim"])
array_dif_salim=np.array(tabela_massas_calculadas['Diferença_Salim'])

kde=scipy.stats.gaussian_kde(array_dif_salim)
array_x_kde_plot=np.linspace(difmin, difmax, 1000)

plt.figure(figsize=(6,6))
plt.hist(array_dif_salim,bins=bins, density=True, color = 'xkcd:lightblue',edgecolor='black')
plt.plot(array_x_kde_plot, kde(array_x_kde_plot),color='xkcd:blue', label="KDE")
plt.legend()
plt.ylabel(r"$\log(M_\ast)$ (Chabrier)")
plt.xlabel(r"$\log(M_\ast) - \log(M_{Salim}$")
plt.savefig('histograma_final_ch.jpg', dpi=900)

plt.figure(figsize=(6,6))
plt.scatter(array_salim, array_mcor, color='xkcd:azure', s=0.5)
plt.plot(array_salim, array_salim, color='black')
plt.plot(array_salim, array_salim-0.5,color='black', linewidth=0.4, label="Dex $\pm 0.5$")
plt.plot(array_salim, array_salim+0.5,color='black', linewidth=0.4)
plt.legend()
plt.ylabel(r"$\log(M_\ast)$ (Chabrier)")
plt.xlabel(r"$\log(M_{Salim})$")
plt.savefig('plot_comparativo_final_ch.jpg', dpi=900)



