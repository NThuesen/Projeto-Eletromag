# Vamos importar os módulos que precisamos
from cmath import *
from numpy import linalg
import numpy as np
import matplotlib.pyplot as plt

# Parâmetros do circuito
R1 = 0.5
R2 = 0.5
L1 = 1e-6
L2 = 1e-6
C1 = 1e-6
C2 = 1e-6
f = 36000 # Frequência
V1 = 5*2/np.pi
w = 2*np.pi*f
ka = 0.2 # fator de acoplamento
M = ka*((L1*L2)**(1/2))

# Reatâncias
XL1 = 1j*w*L1
XL2 = 1j*w*L2
XM = 1j*w*M
XC1 = 1/(1j*w*C1)
XC2 = 1/(1j*w*C2)

# Valores da capacitância
capacitancias = [1.5e-7, 1e-7, 4.7e-7, 1e-6, 4.7e-6]

# Valores da frequência
frequencias = np.arange(0.01, 100e3, 10)

# Funções para calcular a indutância e frequência angular
def CalculaIndutância(C, w):
    L = 1/((w**2)*C)
    return L

def CalculaFrequenciaAngular(f):
    return 2*np.pi*f

# Funções para calcular transformador em série e paralelo
def CalcularTransformadorSerie(Uf, Rc, Rbobina1, Rbobina2, XL, XC, XM):
    Z = np.array([[Rbobina1, -XM], [-XM, Rbobina2 + Rc ]])
    V = np.array([Uf, 0])
    i = np.dot(linalg.inv(Z), V)
    return i[0], i[1]

def CalcularTransformadorParalelo(Uf, Rc, Rbobina1, Rbobina2, XL, XC, XM, ZP):
    Z = np.array([[Rbobina1, -XM], [-XM, Rbobina2 + ZP + XL]])
    V = np.array([Uf, 0])
    i = np.dot(linalg.inv(Z), V)
    return i[0], i[1]

# Inicializando variáveis para o gráfico
Cores = ['b', 'r', 'g', 'purple', 'black']

# Plot para a configuração em série
plt.figure(figsize=(12, 6))
fig, Eixo_X1_S = plt.subplots(figsize=(15, 8))
Eixo_X2_S = Eixo_X1_S.twinx()

for i, C in enumerate(capacitancias):
    Tensoes_S = []
    Eficiências_S = []
    for f in frequencias:
        w = CalculaFrequenciaAngular(f)
        L = CalculaIndutância(C, w)
        M = ka * L
        XL = 1j * w * L
        XC = 1 / (1j * w * C)
        XM = 1j * w * M
        Rbobina1 = R1 + (R1 / 100e3) * f
        Rbobina2 = R2 + (R2 / 100e3) * f

        # Calculando Tensão
        i1, i2 = CalcularTransformadorSerie(V1, 5, Rbobina1, Rbobina2, XL, XC, XM)
        TensaoS = abs(i2 * 5)
        Tensoes_S.append(TensaoS)

        # Calculando eficiência
        P_saida = (TensaoS * i2.conjugate()) / 2
        P_entrada = (V1 * i1.conjugate()) / 2
        N_S = np.real(P_saida) / np.real(P_entrada)
        Eficiências_S.append(N_S * 100)

    # Plotando tensão
    Eixo_X1_S.plot(frequencias, Tensoes_S, linestyle='-', label=f'C = {C * 1e6:.2f} uF', color=Cores[i])
    Eixo_X1_S.set_xlabel("Frequência de ressonância (Hz)")
    Eixo_X1_S.set_ylabel("Tensão (V)")

    # Plotando eficiência
    Eixo_X2_S.plot(frequencias, Eficiências_S, linestyle="--", color=Cores[i])
    Eixo_X2_S.set_ylabel('Eficiência (%)')

plt.title('Eficiência do WPT em Série para diferentes capacitâncias')
fig.legend(loc='upper right', bbox_to_anchor=(1, 0.85))
plt.show()

# Plot para a configuração em paralelo
plt.figure(figsize=(12, 6))
fig, Eixo_X1_P = plt.subplots(figsize=(15, 8))
Eixo_X2_P = Eixo_X1_P.twinx()

for i, C in enumerate(capacitancias):
    Tensoes_P = []
    Eficiências_P = []
    for f in frequencias:
        w = CalculaFrequenciaAngular(f)
        L = CalculaIndutância(C, w)
        M = ka * L
        XL = 1j * w * L
        XC = 1 / (1j * w * C)
        XM = 1j * w * M
        ZP = (XC * 5) / (XC + 5)
        Rbobina1 = R1 + (R1 / 100e3) * f
        Rbobina2 = R2 + (R2 / 100e3) * f

        # Calculando Tensão
        i1, i2 = CalcularTransformadorParalelo(V1, 5, Rbobina1, Rbobina2, XL, XC, XM, ZP)
        TensaoP = abs(i2 * ZP)
        Tensoes_P.append(TensaoP)

        # Calculando eficiência
        P_saidaP = (TensaoP * i2.conjugate()) / 2
        P_entradaP = (V1 * i1.conjugate()) / 2
        N_P = np.real(P_saidaP) / np.real(P_entradaP)
        Eficiências_P.append(N_P * 100)

    # Plotando tensão
    Eixo_X1_P.plot(frequencias, Tensoes_P, linestyle='-', label=f'C = {C * 1e6:.2f} uF', color=Cores[i])
    Eixo_X1_P.set_xlabel("Frequência de ressonância (Hz)")
    Eixo_X1_P.set_ylabel("Tensão (V)")

    # Plotando eficiência
    Eixo_X2_P.plot(frequencias, Eficiências_P, linestyle="--", color=Cores[i])
    Eixo_X2_P.set_ylabel('Eficiência (%)')

plt.title('Eficiência do WPT em Paralelo para diferentes capacitâncias')
fig.legend(loc='upper right', bbox_to_anchor=(1, 0.85))
plt.show()
