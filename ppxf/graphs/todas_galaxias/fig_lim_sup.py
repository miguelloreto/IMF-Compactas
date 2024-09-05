import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import scipy

#Em qual regul quer-se fazer o gráfico?
regul=10

def cria_bins(dados):
    #Cria o número de bins apropriado para dos dados
    max=dados.max()
    min=dados.min()
    dp=dados.std()
    return int((max-min)/(dp*0.5))

def kde_line(dataframe):
    difmax=dataframe["lim_sup"].max()
    difmin=dataframe["lim_sup"].min()
    kde=scipy.stats.gaussian_kde(np.array(dataframe["lim_sup"]))
    array_x_kde_plot=np.linspace(difmin, difmax, 1000)
    return array_x_kde_plot, kde(array_x_kde_plot)

def get_binned(dados, n_bins,medida):
  percentiles=np.arange(1/n_bins, 1, 1/n_bins)
  quantis=[]
  mediana_medida=[]
  mediana_gamma=[]
  erro=[]
  for i in percentiles: quantis.append(dados[medida].quantile(i))
  for i in quantis:
    dados_quantil=dados[dados[medida] < i]
    mediana_medida.append(np.median(dados_quantil[medida]))
    mediana_gamma.append(np.median(dados_quantil.lim_sup))
    erro.append(np.std(dados_quantil.lim_sup))
    dados_upper=dados[dados[medida]> i]
    dados=dados_upper
  return mediana_medida, mediana_gamma, erro

n_bins_hist=15
n_bins_rel=6

initial_path=os.getcwd()
mcgs_path=initial_path+'/mcgs/regul'+str(regul)
csgs_path=initial_path+'/csgs/regul'+str(regul)
os.chdir(mcgs_path)
tabela_mcgs_ls=pd.read_csv("tabela_lim_sup.csv")
tabela_mcgs_ls=tabela_mcgs_ls[tabela_mcgs_ls['lim_sup']  < 5]
tabela_mcgs_ls=tabela_mcgs_ls[tabela_mcgs_ls['lim_sup']  > 0]
os.chdir(csgs_path)
tabela_csgs_ls=pd.read_csv("tabela_lim_sup.csv")
tabela_csgs_ls=tabela_csgs_ls[tabela_csgs_ls['lim_sup']  < 5]
tabela_csgs_ls=tabela_csgs_ls[tabela_csgs_ls['lim_sup']  > 0]

ks=scipy.stats.ks_2samp(np.array(tabela_mcgs_ls['lim_sup']), np.array(tabela_csgs_ls['lim_sup']))

r_pearson_massa_mcgs=tabela_mcgs_ls["log_M"].astype(float).corr(tabela_mcgs_ls['lim_sup'].astype(float), method='pearson')
r_pearson_sigma_mcgs=tabela_mcgs_ls["sigma_e"].astype(float).corr(tabela_mcgs_ls['lim_sup'].astype(float), method='pearson')
tau_kendall_massa_mcgs=tabela_mcgs_ls["log_M"].astype(float).corr(tabela_mcgs_ls['lim_sup'].astype(float), method='kendall')
tau_kendall_sigma_mcgs=tabela_mcgs_ls["sigma_e"].astype(float).corr(tabela_mcgs_ls['lim_sup'].astype(float), method='kendall')

r_pearson_massa_csgs=tabela_csgs_ls["log_M"].astype(float).corr(tabela_csgs_ls['lim_sup'].astype(float), method='pearson')
r_pearson_sigma_csgs=tabela_csgs_ls["sigma_e"].astype(float).corr(tabela_csgs_ls['lim_sup'].astype(float), method='pearson')
tau_kendall_massa_csgs=tabela_csgs_ls["log_M"].astype(float).corr(tabela_csgs_ls['lim_sup'].astype(float), method='kendall')
tau_kendall_sigma_csgs=tabela_csgs_ls["sigma_e"].astype(float).corr(tabela_csgs_ls['lim_sup'].astype(float), method='kendall')

mcgs_m=get_binned(tabela_mcgs_ls, n_bins_rel, 'log_M')
csgs_m=get_binned(tabela_csgs_ls, n_bins_rel, 'log_M')
mcgs_s=get_binned(tabela_mcgs_ls, n_bins_rel, 'sigma_e')
csgs_s=get_binned(tabela_csgs_ls, n_bins_rel, 'sigma_e')

os.chdir(initial_path)

plt.figure(figsize=(5,5))
plt.scatter(np.array(tabela_mcgs_ls["log_M"]), np.array(tabela_mcgs_ls["lim_sup"]), color='xkcd:cherry', s=0.7,label=r'MCGs, r={0:.2f}, $\tau={0:.2f}$ '.format(r_pearson_massa_mcgs, tau_kendall_massa_mcgs))
plt.scatter(np.array(tabela_csgs_ls["log_M"]), np.array(tabela_csgs_ls["lim_sup"]), color='xkcd:green', s=0.7,label=r'CSGs, r={0:.2f}, $\tau={0:.2f}$ '.format(r_pearson_massa_csgs, tau_kendall_massa_csgs))
plt.plot(mcgs_m[0],mcgs_m[1], linestyle='-.', marker='s', color='black')
plt.plot(csgs_m[0],csgs_m[1], linestyle='-.', marker='o', color='black')
plt.ylim([0,4])
plt.ylabel(r"$\Gamma_{max}$")
plt.xlabel(r"$\log(M_{\star})$ ($M_{\odot}$)")
plt.legend()
plt.savefig('limsupxlogm_'+str(regul)+'.png', dpi=900)

plt.figure(figsize=(5,5))
plt.scatter(np.log10(np.array(tabela_mcgs_ls["sigma_e"])), np.array(tabela_mcgs_ls["lim_sup"]), color='xkcd:cherry', s=0.7, label=r'MCGs, r={0:.2f}, $\tau={0:.2f}$'.format(r_pearson_sigma_mcgs, tau_kendall_sigma_mcgs))
plt.scatter(np.log10(np.array(tabela_csgs_ls["sigma_e"])), np.array(tabela_csgs_ls["lim_sup"]), color='xkcd:green', s=0.7, label=r'CSGs, r={0:.2f}, $\tau={0:.2f}$'.format(r_pearson_sigma_csgs, tau_kendall_sigma_csgs))
plt.plot(np.log10(np.array(mcgs_s[0])), mcgs_s[1], linestyle='-.', marker='s', color='black')
plt.plot(np.log10(np.array(csgs_s[0])), csgs_s[1], linestyle='-.', marker='o', color='black')
plt.ylim([0,4])
plt.ylabel(r"$\Gamma_{max}$")
plt.xlabel(r"$\log(\sigma_e)$ (km/s)")
plt.legend()
plt.savefig('limsupxlogsigma_'+str(regul)+'.png', dpi=900)

plt.rcParams.update({'font.size': 15})
plt.figure(figsize=(8,8))
plt.hist(np.array(tabela_mcgs_ls['lim_sup']),bins=n_bins_hist, density=True, color = 'xkcd:cherry',edgecolor='black', label='Distribuição MCGs')
plt.axvline(tabela_mcgs_ls['lim_sup'].mean(), color='xkcd:salmon', ls='dashed', label="Média MCGs",linewidth = 3)
plt.hist(np.array(tabela_csgs_ls['lim_sup']),bins=n_bins_hist, density=True, color = 'xkcd:green',edgecolor='black', label='Distribuição CSGs', alpha=0.7)
plt.axvline(tabela_csgs_ls['lim_sup'].mean(), color='xkcd:lightgreen', ls='dashed', label="Média CSGs", linewidth = 3)
plt.plot([], [], ' ', label=r"KS p-value < {0:.0e}".format(1e-10))
plt.legend()
plt.xlim([0,4])
plt.plot(kde_line(tabela_mcgs_ls)[0],kde_line(tabela_mcgs_ls)[1],color='xkcd:salmon', label="KDE MCGs", linewidth = 3)
plt.plot(kde_line(tabela_csgs_ls)[0],kde_line(tabela_csgs_ls)[1],color='xkcd:lightgreen', label="KDE CSGs", linewidth = 3)
plt.ylabel("Frequência")
plt.xlabel(r"$\Gamma_{max}$")
plt.savefig('histograma_'+str(regul)+'.png', dpi=900)
