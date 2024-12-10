import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
import numpy as np
import os

arquivo = 'protein-peptides_2024_mAb_antiS1_5_enzimas_ISL_PEAKS61.csv'

df = pd.read_csv(arquivo)

ids = df['Protein ID'].unique()

id = input('Insira o ID da proteína: ')

mask = (df['Protein Accession'].values == id)
start = (df['Start'].values)[mask]
end = (df['End'].values)[mask]

contador = np.full(max(end)-min(start)+1, 0)

for j in range(len(start)):
    contador[start[j]-1:end[j]] += 1

linhas = list(range(min(start), max(end)+1))
data = {'Line': linhas, 'Count': contador}

result = pd.DataFrame(data=data)

os.makedirs(f'{id}', exist_ok=True)

result.to_csv(f'{id}/result_{id}.csv', index=False)

plt.figure(figsize=(30,8))
plt.style.use('seaborn-darkgrid')

plt.bar(linhas, contador, color ='maroon')
plt.xlabel('Linha')
plt.ylabel('Contagem')

plt.title(f'Contagem da proteína {id}')

plt.savefig(f'{id}/Cobertura_{id}.png')