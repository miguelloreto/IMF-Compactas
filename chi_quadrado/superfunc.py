from matplotlib.figure import Figure
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

def corte(vetornaocortadox, vetornaocortadoy, inicio, fim):
  vetorcortadox=[]
  vetorcortadoy=[]
  for i in np.arange(0, len(vetornaocortadox), 1):
    if vetornaocortadox[i] > inicio and vetornaocortadox[i] < fim: 
      vetorcortadox.append(vetornaocortadox[i])
      vetorcortadoy.append(vetornaocortadoy[i])
  return vetorcortadox, vetorcortadoy

def compespectro(inicio, fim, sigma):
  lam=np.arange(inicio,fim, 0.001)
  lam_0=((lam[-1]-lam[0])/2) + lam[0]
  dellambda=sigma*lam/c
  return lam, lam_0, dellambda   

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

def gaussiana(lam, lam_0, dellambda):
  return 1/(np.sqrt(2*pi)*dellambda)*np.exp(-np.power((lam - lam_0)/dellambda, 2)/2)

def fatordenormalizacao(fluxo1, fluxo2):
  fator=np.mean(fluxo1)/np.mean(fluxo2)
  return fator

#Vetor dados
dados=['01.txt','02.txt', '03.txt', '04.txt', '05.txt', '06.txt', '07.txt', '08.txt']
observado='galaxia_observada.txt'

def superfunc(modelos, observado,  inicio, fim, sigma):
  lamo, fluxo, nada, nada=np.loadtxt(observado, unpack='True')

  lam=[]
  flu=[]
  for i in modelos:
    l, f=np.loadtxt(i, unpack='True')
    lam.append(l)
    flu.append(f)

  lamc=[]
  fluc=[]
  for i in range(0, len(modelos), 1):
    lc, fc=(corte(lam[i], flu[i], inicio, fim))
    lamc.append(lc)
    fluc.append(fc)

  lgaus, lgaus_0, delgaus = compespectro(inicio, fim, sigma)

  flucI=[]
  for i in range(0, len(modelos), 1):
    interpolado=np.interp(lgaus, lamc[i], fluc[i])
    flucI.append(interpolado)

  convol=[]
  for i in range(0, len(modelos), 1):
    convoluido=np.convolve(gaussiana(lgaus, lgaus_0, delgaus), flucI[i], mode='same')
    convol.append(convoluido)
    
  def make_graph(numero, lam, fluxointerpolado, fluxoconvoluido, lamoc, fluxoc):
    graph=plt.figure(figsize=(8,8))
    plt.title('Espectros na região de 5800 até 6100')
    plt.xlabel('$\lambda$($\AA$)')
    plt.ylabel('Fluxo ($erg/cm^2/s/A$)')
    plt.xlim([5800,6100])
    plt.plot(lam, fluxointerpolado , color='red', label='Modelo{numero}'.format(numero=numero+1))
    plt.plot(lam, fluxoconvoluido/fatordenormalizacao(fluxoconvoluido, fluxointerpolado), color='green', label=('Convolução'))
    plt.plot(lamoc, fluxoc/fatordenormalizacao(fluxoc, fluxointerpolado),color='blue', label=('Observado'))
    plt.legend()
    plt.savefig('Modelo {numero}.png'.format(numero=numero), dpi=300)
    return graph
    
  supergraph=[]

  for i in range(0, len(modelos), 1):
    grafico= make_graph(i, lgaus, flucI[i], convol[i], lamo, fluxo)
    supergraph.append(grafico)

  return supergraph


listagraficos=superfunc(dados, observado, 5800, 6100, 300)