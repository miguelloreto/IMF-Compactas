#---------
# PACOTES
#---------
import matplotlib.pyplot as plt
import numpy as np
import math

sigma=300
c=3e5 #Velocidade da luz em Km/s
pi=np.pi
lambdamodelo, fluxomodelo=np.loadtxt('imf30.txt', unpack='True')
lambdaobservado, fluxoobservado, nada, nada =np.loadtxt('galaxia_observada.txt', unpack='True')


#Criando o vetor dos lambdas para todo o espectro
lambdacriado=np.arange(lambdamodelo[0], lambdamodelo[-1], 0.01)
#deltalambda/lambda=R=v/c então deltalambda=lambda*v/c. A divisão pelo (2 * math.sqrt(2 * np.log(2))) é a largura meia-altura da gaussiana
deltalambda=lambdacriado*(sigma/c) / (2 * math.sqrt(2 * np.log(2)))
#Interpolação do Python. Dar uma olhada na interpolação da Marina
fluxomodelointerpolado=np.interp(lambdacriado, lambdamodelo,fluxomodelo)

def fatordenormalizacao(fluxo1, fluxo2):
  fator=np.mean(fluxo1)/np.mean(fluxo2)
  return fator

#Gaussiana já normalizada
def gaussiana(lam, lam_0, dellambda):
  return 1/(np.sqrt(2*pi)*dellambda)*np.exp(-np.power((lam - lam_0)/dellambda, 2)/2)

#Definindo um espectro convoluido que tem len()=len(fluxomodelointerpolado)=len(lambdacriado)
convoluido=np.zeros(len(lambdacriado))
#Definindo um vetor "de corte" j
j=[]

for i in range(0, len(lambdacriado)):
  #Definindo um intervalo de corte
  intervalo=20*deltalambda[i]
  #Vetor j definido para cortar a gaussiana. Ela vale 0 para valores muito distantes de lambda_0 por isso vamos centrando ela dentro de um intervalo.
  #Tentei fazer isso com um if mas é muito demorado e não funciona.
  #Essa próxima linha de código?
  j=np.where(np.logical_and(lambdacriado<=(lambdacriado[i]+intervalo), lambdacriado>=(lambdacriado[i]-intervalo)))
  #Nossa psf é nossa gaussiana:
  psf=gaussiana(lambdacriado[i], lambdacriado[j], deltalambda[i])
  #Essa é a integral. len(psf)=len(fluxointerpolado[j]). Multiplica os vetores (definição de convolução) e integra (soma o valor da função multiplicada naquele ponto) com o np.sum().
  convoluido[i]=np.sum(psf*fluxomodelointerpolado[j])


#Gráfico convolução total
g1=plt.figure(figsize=(8,8))
plt.title('Espectros na região de 3500 até 7400')
plt.xlabel('$\lambda$($\AA$)')
plt.ylabel('Fluxo ($erg/cm^2/s/A$)')
plt.xlim([3500,7400])
plt.plot(lambdacriado, fluxomodelointerpolado, color='red', label='IMF30')
plt.plot(lambdacriado, convoluido/fatordenormalizacao(convoluido, fluxomodelointerpolado), color='green', label=('Convolução'))
plt.plot(lambdaobservado, fluxoobservado/fatordenormalizacao(fluxoobservado, fluxomodelointerpolado),color='blue', label=('Observado'))
plt.legend()
plt.savefig('convolucaototalimf30.png', dpi=300)