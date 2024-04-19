# MÓDULOS
import os 
import numpy as np
import shutil
import sys

# Programa que cria um arquivo de base chamado "BaseGM_LCGs" e deixa o ppxf pronto para rodar com uma nova base de arquivos.
# Usa os espectros baixados do MILES e uma tabela de massas remanescentes, os dois em uma pasta chamada newbase.
# A nova base pode ter qualquer tipo de IMF, qualquer slope, qualquer metalicidade...
# O formato do arquivo de bases se dá em 4 colunas: nome do espectro, idade, metalicidade e massa remascente.


#Lidando com caminhos de diretório:
ppxf_path=os.getcwd()
newbase_path=ppxf_path +'/newbase/'
basesGM_path=ppxf_path+'/basesGM/'


#Elemento recursivo: eliminando o que foi feito com a base anterior
espectros_anteriores=os.listdir(newbase_path)
os.chdir(newbase_path)
for i in espectros_anteriores:
    if i[0]=='M': os.remove(i)
os.chdir(ppxf_path)

#Primeiro, seleciona-se o slope e a IMF desejada para a nova base.
slope=float(sys.argv[1])
imf_type=sys.argv[2]
#Depois, olhamos para as bases possíveis no basesGM:
espectros_possiveis=os.listdir(basesGM_path)
for i in espectros_possiveis:
    if(i[0]!='M'):
        continue
    if(float(i[3:7])==slope and i[1:3]==imf_type):
        shutil.copy(basesGM_path + i, newbase_path)

if (len(espectros_possiveis) == 0): raise("A base de dados não contém esse slope")

# Cria o novo arquivo das bases no diretório newbase
os.chdir(newbase_path)
base=open('BaseGM_LCGs', 'w')
os.chdir(ppxf_path)
# Guarda o nome dos espectros no newbase. Esses serão os espectros utilizados para fazermos nossa base.
espectros=os.listdir(newbase_path)
for i in espectros:
    if i[0]!='M':
        espectros.remove(i)
    else:
        continue
# Abre a tabela de massas. 
tabela=open(newbase_path+'/tabelamassas'+imf_type+'baseFe.txt', 'r')

# O desafio agora é casar a informação do nome dos espectros com as idade, metalicidade e massas remascentes possíveis.
# Quais são as idades e metalicidades possíveis?
metalicidades=np.array([-0.96, -0.66, -0.35, -0.25, 0.06, 0.15, 0.26, 0.4]) # Em [M/H]
metalicidades=10**(metalicidades + np.log10(0.019)) # Conversão de [M/H] para Z
metalicidades=tuple(metalicidades)
idades=np.arange(np.log10(0.03*1e9), np.log10(15*1e9), 0.2) # Seleção de idades, conforme instruido pela Marina
idades=(10**np.array(idades))/1e9
idades=np.append(idades, 13.5)
idades=tuple(idades)
# Sei minhas idades e metalicadades possíveis. Essas vão ser 

# Ok, agora tenho que extrair as informações dos espectros. Vou organizar todas informações dele em uma classe:
class specter:
    def __init__(espectro, nome, imftipo, imfslope, metal, idade, mreman):
        espectro.nome=nome
        espectro.tipo=imftipo
        espectro.slope=imfslope
        espectro.metal=metal
        espectro.idade=idade
        espectro.mreman=mreman

# Fazer uma função que pega os dados do nome dos arquivos, que é uma string:
def get_info_specter(nome):
    '''Função que extrai os dados do nome do arquivo
    Note que todo espectro gerado pelo MILES tem o MESMO FORMATO DE NOME DE ARQUIVO, então isso funcionaria pra qualquer espectro.'''
    tipo=nome[1:3]
    slope=float(nome[3:7])
    if nome[8]=='m':
        metal=-1*float(nome[9:13])
        metal=10**(metal + np.log10(0.019)) # Conversão de [M/H] para Z
        idade=float(nome[14:21])
    else:
        metal=float(nome[9:13])
        metal=10**(metal + np.log10(0.019)) # Conversão de [M/H] para Z
        idade=float(nome[14:21])    

    return specter(nome, tipo, slope, metal, idade, False)

# Armazenando todos os espectros do diretório em um vetor espectros diretório
espectros_diretorio=[]
for i in espectros:
    espectros_diretorio.append(get_info_specter(i))

# Fazer uma função que pega os dados da tabela e seleciona com qual slope a nova base será criada. Note que a seleção do tipo de IMF está a cargo da tabela!
def get_info_tabela (tabela, slope):
    espectros_tabela=[]
    j=0
    for i in tabela:
        a=[]
        a = i.split()
        if len(a) == 0:
            continue
        else:
            j=j+1
            nome='Linha {}'.format(j)
            tipo=a[0].lower()
            slopetabela=float(a[1])
            metal=float(a[2])
            metal=10**(metal + np.log10(0.019)) # Conversão de [M/H] para Z
            idade=float(a[3])
            mreman=float(a[5])
            if (slopetabela== slope):
                espectro=specter(nome, tipo, slope, metal, idade, mreman)
                espectros_tabela.append(espectro)
    return espectros_tabela

espectros_tabela=get_info_tabela(tabela, slope)

# Agora tenho que organizar as idades e metalicidades de todos espectros de forma a escolher os espectros que melhor se aproximem das minhas tuplas idades e metalicidades.
# Tenho que fazer isso pro diretório e pra tabela 

def seletor(vetora, vetorb):
    ''' 
    Função que dado um vetor A de referencia, seleciona os valores de um vetor B que melhor se aproximam dos valores de A.
    Retorna um vetor no shape do vetor A, com os valores de B.
    '''
    vetor_convertido=[]
    for i in vetora:
        absold=1e1000
        for j in vetorb:
            absnew=abs(i-j)
            if (absnew<absold):
                absold=absnew
                selecionado=j    
        vetor_convertido.append(selecionado)
    return vetor_convertido

def aproximador(valor, vetor):
    '''
    Função para ser usada juntamente com a função seletor, aproxima um valor de uma lista de vetores possiveis
    '''
    absold=1e10
    for i in vetor:
        absnew=abs(valor-i)
        if (absnew<absold):
            absold=absnew
            traduzido=i
    return traduzido
############
# DIRETÓRIO
############
idades_selecionadas=[]
for i in espectros_diretorio:
    idades_selecionadas.append(i.idade)
idades_selecionadas=list(set(idades_selecionadas))
idades_selecionadas=seletor(idades, idades_selecionadas)

metalicidades_selecionadas=[]
for i in espectros_diretorio:
    metalicidades_selecionadas.append(i.metal)
metalicidades_selecionadas=list(set(metalicidades_selecionadas))
metalicidades_selecionadas=seletor(metalicidades, metalicidades_selecionadas)

espectros_selecionados_diretorio=[]
for i in espectros_diretorio:
    if(idades_selecionadas.count(i.idade) ==1 and metalicidades_selecionadas.count(i.metal) == 1):
        i.metal=aproximador(i.metal, metalicidades)
        espectros_selecionados_diretorio.append(i)

#########
# TABELA
#########
idades_selecionadas=[]
for i in espectros_tabela:
    idades_selecionadas.append(i.idade)
idades_selecionadas=list(set(idades_selecionadas))
idades_selecionadas=seletor(idades, idades_selecionadas)

metalicidades_selecionadas=[]
for i in espectros_tabela:
    metalicidades_selecionadas.append(i.metal)
metalicidades_selecionadas=list(set(metalicidades_selecionadas))
metalicidades_selecionadas=seletor(metalicidades, metalicidades_selecionadas)

espectros_selecionados_tabela=[]
for i in espectros_tabela:
    if(idades_selecionadas.count(i.idade) ==1 and metalicidades_selecionadas.count(i.metal) == 1):
        i.metal=aproximador(i.metal, metalicidades)
        espectros_selecionados_tabela.append(i)


# Parte onde eu escrevo tudo no novo arquivo BaseGM_LCGs. O objetivo aqui é achar um metalicade remanescente (que tá na tabela) e casar ela com os dados do espectro.
# Claro, a massa remascente da tabela, deve ter o mesmo slope, idade e metalicidade da do espectro. 
espectros_finais=[]
for i in espectros_selecionados_tabela:
    for j in espectros_selecionados_diretorio:
        if (i.tipo != j.tipo):
            continue
        if(i.slope == j.slope and i.idade == j.idade and i.metal == j.metal):
            final=specter(j.nome, j.tipo, j.slope, j.metal, j.idade, i.mreman)
            espectros_finais.append(final)
espectros_finais_sorted=sorted(espectros_finais, key=lambda specter:specter.idade)

for i in espectros_finais_sorted:
    base.write('{nome}\t{idade:.4e}\t{metalicidade:.7f}\t{mreman:.4f}\n'.format(nome=i.nome, idade=i.idade*1e9, metalicidade=i.metal, mreman=i.mreman))

# Fecha a tabela os espectros
tabela.close()
base.close()

# Se tiver um BaseGM_LCGs, descarte ele. Mova para a pasta fora da que estamos trabalhando, tomando o lugar do BaseGM_LCGs original.
dirppxf=os.listdir(ppxf_path)
for i in dirppxf:
    if i == 'BaseGM_LCGs':
        os.remove(ppxf_path +'/'+ i)
shutil.move(newbase_path+'BaseGM_LCGs', ppxf_path)

# Os especrtros que vão servir de base estão na pasta basesGM, então tenho que passar todos os espectros da minha nova base para lá.
# Mas o programa roda só com os espectros selecionados pelo arquivo BasesGM_LCGs