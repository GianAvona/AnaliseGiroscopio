# -*- coding: utf-8 -*-
"""

NOME: GIANFRANCO AVONA COVELLO          TIA: 3185692-6
NOME: RAMON SOUZA ARAUJO                TIA: 3190715-6

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

### Importando dados da WEB, arquivo CSV
url = 'https://raw.githubusercontent.com/pjasimoes/PIPythonData/main/giroscopio1.csv'
data = pd.read_csv(url, sep= ",")

### Nomeando com variaveis as colunas dos dados obtidos
k = data.keys()
t = k[0]; x = k[1]; y = k[2]; z = k[3]; a = k[4]

"""
###Plotando gráficos das leituras de cada Eixo, sem filtro algum
fig, ax = plt.subplots(4, 1, figsize = (18,10))

ax[0].set_title('Velocidade Angular com Ruído', fontsize = '24', color = 'black')

ax[0].plot(data[t], data[x], color = 'blue')
ax[0].set_ylabel(x)

ax[1].plot(data[t], data[y], color = 'red')
ax[1].set_ylabel(y)

ax[2].plot(data[t], data[z], color = 'green')
ax[2].set_ylabel(z)

ax[3].plot(data[t], data[a], color = 'k')
ax[3].set_ylabel(a)
ax[3].set_xlabel(t)

for i in ax:
    i.grid()
    
"""    
### Removendo o ruído do sinal, aplicando filtros

wx = (data[x].rolling(100, min_periods=0,center=True).mean())
wy = (data[y].rolling(100, min_periods=0,center=True).mean())
wz = (data[z].rolling(100, min_periods=0,center=True).mean())
wa = (data[a].rolling(100, min_periods=0,center=True).mean())

### Removendo os ruídos e movimentos indesejados antes dos 10s,
### Onde nosso sensor esta parado e não deveria registrar movimentos
ind = np.where(data[t] < 10)
## Aplicando a média onde t < 10s
mx = wx.iloc[ind].mean()
my = wy.iloc[ind].mean()
mz = wz.iloc[ind].mean()
ma = wa.iloc[ind].mean()

wxf = wx - mx
wyf = wy - my
wzf = wz - mz
waf = wa - ma

### Gráficos dos eixos plotados na mesma figura, do sinal com e sem ruído e filtrado
### Antes dos 10s de medição
fig, b = plt.subplots(4, 1, figsize = (18, 10))
b[0].set_title('Velocidade Angular', fontsize = 24)

b[0].plot(data[t], data[x], label = 'Com ruído', color = 'blue')
b[0].legend(loc = 'upper left')
b[0].plot(data[t], wxf, label = 'Filtrado', color = 'orange')
b[0].legend(loc = 'upper left')
b[0].set_ylabel(x)

b[1].plot(data[t], data[y], label = 'Com ruído', color = 'red')
b[1].legend(loc = 'upper left')
b[1].plot(data[t], wyf, label = 'Filtrado', color = 'orange')
b[1].legend(loc = 'upper left')
b[1].set_ylabel(y)

b[2].plot(data[t], data[z], label = 'Com ruído', color = 'green')
b[2].legend(loc = 'upper left')
b[2].plot(data[t], wzf, label = 'Filtrado', color = 'orange')
b[2].legend(loc = 'upper left')
b[2].set_ylabel(z)

b[3].plot(data[t], data[a], label = 'Com ruído', color = 'k')
b[3].legend(loc = 'upper left')
b[3].plot(data[t], waf, label = 'Filtrado', color = 'orange')
b[3].legend(loc = 'upper left')
b[3].set_ylabel(a)
b[3].set_xlabel(t)

for i in b:
    i.grid()

# ------ Dataframe para Array ------- #
## Dados transformados foram os sem ruídos 
# ------ Acel. Ang em função do tempo de cada coordenada ------- #

grad_x = np.gradient(wxf, data[t])
grad_y = np.gradient(wyf, data[t])
grad_z = np.gradient(wzf, data[t])
grad_a = np.gradient(waf, data[t])

fig, acel_ang = plt.subplots(3, 1, figsize = (18,10))

## acel_ang = Aceleração ANGULAR
acel_ang[0].set_title('Aceleração Angular', fontsize = '24', color = 'black')

acel_ang[0].plot(data[k[0]], grad_x, color = 'blue')
acel_ang[0].set_ylabel([k[1]])

acel_ang[1].plot(data[k[0]], grad_y, color = 'orange')
acel_ang[1].set_ylabel([k[2]])

acel_ang[2].plot(data[k[0]], grad_z, color = 'green')
acel_ang[2].set_ylabel([k[3]])

acel_ang[2].set_xlabel([k[0]], fontsize = '18', color = 'black')

for i in acel_ang:
    i.grid()


### APLICANDO A INTEGRAL, PARA PODER ENCONTRAR O DESLOCAMENTO DE CADA EIXO NO TEMPO

from scipy.integrate import trapz
import math

dx = trapz(wxf, data[t])
dx = round(dx * (180 / math.pi), 5)
print('Movimento em Graus - Eixo X: ', dx,'º')

dy = trapz(wyf, data[t])
dy = round(dy * (180 / math.pi), 5)
print('Movimento em Graus - Eixo Y: ', dy,'º')

dz = trapz(wzf, data[t])
dz = round(dz * (180 / math.pi), 5)
print('Movimento em Graus - Eixo Z: ', dz,'º')

da = trapz(waf, data[t])
da = round(da * (180 / math.pi), 5)
print('Movimento em Graus - Eixo Absoluto: ', da,'º')

### APLICANDO A DERIVADA, PARA PODER ENCONTRAR A VARIAÇÃO DA POSIÇÃO ANGULAR NO TEMPO
from scipy.integrate import cumtrapz
dt = data['Time (s)']
dvx = cumtrapz(wxf, dt)
dvy = cumtrapz(wyf, dt)
dvz = cumtrapz(wzf, dt)
dva = cumtrapz(waf, dt)

fig, c = plt.subplots(4, 1, figsize = (18, 10))
c[0].set_title('Posição Angular no Tempo', fontsize = 24)

c[0].plot(dt[1::], dvx, label = 'x', color = 'blue')
c[0].legend(loc = 'upper left')
c[0].set_ylabel(x)

c[1].plot(dt[1::], dvy, label = 'y', color = 'red')
c[1].legend(loc = 'upper left')
c[1].set_ylabel(y)

c[2].plot(dt[1::], dvz, label = 'z', color = 'green')
c[2].legend(loc = 'upper left')
c[2].set_ylabel(z)

c[3].plot(dt[1::], dva, label = 'abs', color = 'k')
c[3].legend(loc = 'upper left')
c[3].set_ylabel(a)
c[3].set_xlabel(t)

for i in c:
    i.grid()

print("")
print("")
resp_a = '6a- Visualizando os gráficos, vemos que os movimentos foram no eixo Z, eixo X e eixo Y, respectivamente'
print(resp_a)
print("")
resp_b = '6b- A maior velocidade foi no Eixo Y e registrou', round(wyf.max(), 5), 'rad/s'
print(resp_b)
print("")
resp_c = '6c- O eixo Y tambem teve o maior deslocamento angular, deslocando em graus' , dy, 'graus'
print(resp_c)
print("")
resp_d = '6d- A maior aceleração também foi no Eixo Y e registrou', round(grad_y.max(), 5), 'rad/s^2'
print(resp_d)
print("")













