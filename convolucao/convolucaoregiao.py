# -*- coding: utf-8 -*-
#---------
# PACOTES
#---------
import matplotlib.pyplot as plt
import numpy as np
#------------------------------------
# DEFINIR CONSTANTES, CARREGAR DADOS
#------------------------------------
c=3e5 #Velocidade da luz em Km/s
pi=np.pi
lambdamodelo, fluxomodelo=np.loadtxt('imf30.txt', unpack='True')
lambdaobservado, fluxoobservado, nada, nada =np.loadtxt('galaxia_observada.txt', unpack='True')

'''----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'''
#------------------------------------
# DEFININDO UMA REGIÃO DO ESPECTRO
#------------------------------------

#Cortando os espectros simulados
def corte(vetornaocortadox, vetornaocortadoy, inicio, fim):
  vetorcortadox=[]
  vetorcortadoy=[]
  for i in np.arange(0, len(vetornaocortadox), 1):
    if vetornaocortadox[i] > inicio and vetornaocortadox[i] < fim: 
      vetorcortadox.append(vetornaocortadox[i])
      vetorcortadoy.append(vetornaocortadoy[i])
  return vetorcortadox, vetorcortadoy

lambdamodelocortado, fluxomodelocortado=corte(lambdamodelo, fluxomodelo, 5800, 6100)
lambdaobservadocortado, fluxoobservadocortado=corte(lambdaobservado, fluxoobservado, 5800, 6100)

#Definido os lambdas a serem usados na gaussiana
def compespectro(inicio, fim, sigma):
  lam=np.arange(inicio,fim, 0.001)
  lam_0=((lam[-1]-lam[0])/2) + lam[0]
  dellambda=sigma*lam/c
  return lam, lam_0, dellambda    

lambdacriado, lambdacriado_0, deltalambdacriado = compespectro(5800, 6100, 165)

'''----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'''
#------------------------------------
# INTERPOLAÇÃO
#------------------------------------
#Interpolarando dados criados e do modelo para os dois terem mesma dimensão
#fluxomodelointerpolado=np.interp(lambdacriado, lambdamodelo, fluxomodelo)
#Testando para ver se a interpolação e os cortes funcionaram:
#plt.plot(lambdacriado, fluxomodelointerpolado)
#plt.plot(lambdamodelocortado, fluxomodelocortado)

#Função de Interpolação
def interpolacao_linear(modelx, modely, x):

  i = 0
  #Enquanto o vetor original tiver um elemento menor que x, não se tem interesse porque está fora da reta que une dois pontos no vetor original 
  while(modelx[i]<x):
    i+=1
  x0 = modelx[i-1]
  y0 = modely[i-1]
  x1 = modelx[i]
  y1 = modely[i]

  #Uma vez dentro da reta, define-se a inclinação
  a = (y1 - y0)/(x1 - x0)
  #Fórmula da Interpolação
  return (y0 + a*(x-x0))

#Agora isso tem que ser repetido para todo x.
fluxomodelointerpolado2=[]
for i in lambdacriado:
  fluxomodelointerpolado2.append(interpolacao_linear(lambdamodelo, fluxomodelo, i))

plt.plot(lambdacriado, fluxomodelointerpolado2)
plt.plot(lambdamodelocortado, fluxomodelocortado)
plt.show
'''----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'''
#-----------
# GAUSSIANA
#-----------
# Montando a gaussiana que vamos usar na convolução.
def gaussiana(lam, lam_0, dellambda):
  return 1/(np.sqrt(2*pi)*dellambda)*np.exp(-np.power((lam - lam_0)/dellambda, 2)/2)
#Testando para ver se a Gaussiana Funciona: 
#plt.plot(lambdacriado, gaussiana(lambdacriado, lambdacriado_0, deltalambdacriado))

'''----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'''
#------------
# CONVOLUÇÃO
#------------

#Função de convolução do numpy, recebe dois sinais de argumento, no modo 'full' len(convolucao) = len(fluxo) + len(gaussiana) - 1.
#convolucao1=np.convolve(gaussiana(lambdacriado, lambdacriado_0, deltalambdacriado), fluxomodelointerpolado, mode='full')
#Por isso o ajuste nos lambdas: Dobro -1 de pontos
#plt.plot(np.arange(5800, 6100-0.0005, 0.0005), convolucao1)

#Outro modo 'same': len(convolucao) = len(fluxo) = len(gaussiana)
convolucao2=np.convolve(gaussiana(lambdacriado, lambdacriado_0, deltalambdacriado), fluxomodelointerpolado2, mode='same')
#plt.plot(lambdacriado, convolucao2)

'''----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'''
#-------------
#NORMALIZAÇÃO
#-------------

def fatordenormalizacao(fluxo1, fluxo2):
  fator=np.mean(fluxo1)/np.mean(fluxo2)
  return fator

'''----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'''
#---------
# OUTPUT
#---------

#dadosconvolucao1=np.column_stack([np.arange(5800, 6100 - 0.0005, 0.0005), convolucao1])
#np.savetxt('convolucao1.txt', dadosconvolucao1, fmt=['%.6e','%.6e'])
dadosconvolucao2=np.column_stack([lambdacriado, convolucao2])
np.savetxt('convolucao2.txt', dadosconvolucao2, fmt=['%.6e','%.6e'])

#Gráfico final:
g1=plt.figure(figsize=(8,8))
plt.title('Espectros na região de 5800 até 6100')
plt.xlabel('$\lambda$($\AA$)')
plt.ylabel('Fluxo ($erg/cm^2/s/A$)')
plt.xlim([5800,6100])
plt.plot(lambdacriado, fluxomodelointerpolado2, color='red', label='IMF30')
plt.plot(lambdacriado, convolucao2/fatordenormalizacao(convolucao2, fluxomodelointerpolado2), color='green', label=('Convolução'))
plt.plot(lambdaobservadocortado, fluxoobservadocortado/fatordenormalizacao(fluxoobservadocortado, fluxomodelointerpolado2),color='blue', label=('Observado'))
plt.legend()
plt.savefig('5800_6100.png', dpi=300)