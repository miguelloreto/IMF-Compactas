import pandas as pd

tabela_propriedades=pd.read_csv("table_compact_stellar_population_properties.csv")

plate=[]
mjd=[]
fiberid=[]
for i in tabela_propriedades['aid']:
    nome_fibra_separado=i.split("-")
    if(nome_fibra_separado[1][0]=='0'):
        nome_fibra_separado[1]=nome_fibra_separado[1][1]+nome_fibra_separado[1][2]+nome_fibra_separado[1][3]
    print(nome_fibra_separado[1])
    plate.append(nome_fibra_separado[1])
    mjd.append(nome_fibra_separado[2])
    fiberid.append(nome_fibra_separado[3])

dicionario={'plate': plate, 'mjd': mjd, 'fiberid': fiberid}
tabela_upload_cassjobs=pd.DataFrame(data=dicionario)
tabela_upload_cassjobs.to_csv('myTable1_in_miguel_loreto.csv', index=False)
    
    