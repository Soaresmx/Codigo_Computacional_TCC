# -*- coding: utf-8 -*-

"""
Universidade Tecnológica Federal do Paraná - UTFPR-CT 2023
Trabalho de Conclusão de Curso - TCC02
Aluno : Daniel Souza Soares 
Orientador: Victor Frencl ; Coorientador: Alexandre Tuoto. 
Código para Sintonia de PID - Planta Didática de Vazão e pH.
""" 

# Importando as bibliotecas para as funções no Python.
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import control as ctl

plt.close("all")

# Método para escrita da FT da planta de uma maneira mais ágil 

s = ctl.TransferFunction.s

# Incializando os parâmetros do controlador PID zerados

kp = 0; ki = 0; kd = 0; R = 0; L = 0;

# FUNCOES PARA OS CONTROLES ---------------------------------------------------#

# Parâmetros para a sintonia de Ganhos de Controlador PID por ZN em MA, primeiro método

# Controlador PID, FT: G(s) = Kp + Ki/s + Kd*s
def Clt_PID_ZN(t, L):
    Kp1 = (1.2*(t/L));
    Ti1 = 2*L;
    Td1 = (0.5*L);
    ZN = Kp1*(1+(1/Ti1*s)+Td1*s);
    return ZN 

# Controlador PID ZN Resp. Dom. Freq. (Técnica para um sistema em MF), segundo método. 
def Clt_PID_ZN_freq(K_cr,PR_cr):
    ZNF = (0.6*K_cr) * (1+ (1/(0.5*PR_cr)*s) + (PR_cr*0.125*s))
    return ZNF

# Parâmetros para a sintonia de Ganhos de Controlador PID por CC em MA
def Clt_PID_CC(theta, tau, Kcc):
    Kp3 = (1.25+0.25*(theta/tau))*(tau)/(Kcc*theta);
    Ti3 = ((1.35+0.25*(theta/tau))/(0.54+0.33*(theta/tau)))*theta;
    Td3 = (0.5*theta) /(1.35+0.25*(theta/tau));
    CC = Kp3*(1+ (1/Ti3*s) + Td3*s);
    return CC

#----------------------------------------------------------------------------#

# Obtendo a FT do Matlab Resp ao Degrau 
num = []
den = []

# Selecionando o modelo de ref para gerar a minha FT do processo 
control = open("FT_planta.csv", "r")

# Realizando a leitura dos parâmetros da Tabela em CSV da FT do Matlab Resp ao Degrau
lines = control.readlines()

for line in lines :
    parameters = line.split(",")
    num.append(float(parameters[0]))
    den.append(float(parameters[1]))

control.close()

# Converter as strings em números float usando compreensões de lista
num = [float(string) for string in num]
den = [float(string) for string in den]

#--------------------------------------------------------------------------------- #

# Criando a FT do Processo
Ps_TF = signal.TransferFunction(num, den)

#------------------------------------------------------------------------
# Para o arquivo  -----> control = open("FT_planta.csv", "r")

num_aux = (num[0]*s**2 + num[1]*s**1 + num[2]*s**0)
den_aux = (den[0]*s**3 + den[1]*s**2 + den[2]*s**1 + den[3]*s**0)

#------------------------------------------------------------------------

#tempo de Amostragem do sinal 
Tsim = 1000
# Sistema de Realimentação 
H_s =  1

# Ajustes para a FT da Planta em uso para uso das bibliotecas <signal, ctl>
ps_tf = num_aux/den_aux
print('FT da Planta: ',ps_tf)

# Resp ao Degrau da Própria FT 
T, out_pid = ctl.step_response(ps_tf,Tsim);

# Tempo de amostragem para o sinal degrau 
t = np.linspace(0, Tsim, 1000)  # Cria um vetor de tempo de 0 a 5 segundos

# Amplitude do degrau usada no teste de identificação da FT da planta 
amplitude = 1
# Função degrau
degrau = np.where(t >= 1, amplitude, 0)  

# --------------------------------------------------------------------------------- #

# Resp ao Degrau do PID Tune usando o método de ZN - Resp em Freq. 

#calculo do Ganho Crítico: K_cr, e o tempo de oscilação: T_cr
out = ctl.stability_margins(ps_tf)

# Valores obtidos pela análise no MATLAB
K_cr= 1.6143 
PR_cr= (34.3196) 

W_cr = 2* np.pi/PR_cr
G = Clt_PID_ZN_freq(K_cr,PR_cr)
G_s = ctl.series(G,ps_tf);
G1_s=ctl.feedback(G_s, H_s, sign= -1);

# Sintoniai ZN Resposta em Frequência
print ("-------------Sintonia de ZN Resp em Frequência----------")
print ('Kp= ',(0.6*K_cr) ,'Ti= ',(PR_cr*0.5) , 'Td= ', (0.125*PR_cr))
print('Ganho crítico = ',K_cr )
print('Período crítico = ',PR_cr ) 
print('Frequência crítica = ',W_cr )
print ("--------------------------------------------------------")
T1, out_pid1 = ctl.step_response(G1_s,Tsim);

# --------------------------------------------------------------------------------- #

# Sintonia ZN Resposta ao Degrau, curva em "S"
M = 1.36 ; tx = (29-12); L = 12;
#M = 1.36 ; tx = (29-12); L = 12;
#M = 1.36 ; tx = 31; L = 10;
R = M/tx; 

# Sintonia PID  ZN Resposta ao deg, curva "S"
G2 = Clt_PID_ZN(tx, L);
G2_s = ctl.series(G2,ps_tf);
G22_s=ctl.feedback(G2_s, H_s, sign= -1);

print ("-------------Sintonia de PID - ZN Resp ao Degrau----------")
print('Kp = ', ((1.2*(tx/L))),'Ti = ',(2*L), 'Td= ',(0.5*L))
print ("--------------------------------------------------------")
T2, out_pid2 = ctl.step_response(G22_s,Tsim);

# --------------------------------------------------------------------------------- #

# Sintonia CC Resposta ao Degrau, curva de reação em formato de "S"
tau = (67-(29+12)); theta = 12; FI = theta/tau; Kcc= 1.36
#tau = 17; theta = 12; FI = theta/tau; Kcc= 1.36

# Sintoniai CC Resposta em MA, curva "S"
G4 = Clt_PID_CC(theta, tau, Kcc)
G4_s = ctl.series(G4,ps_tf);
G44_s=ctl.feedback(G4_s, H_s, sign= -1);

print ("-------------Sintonia de CC Resp ao Degrau----------")
print ('Kp= ',((1.25+0.25*(theta/tau))*(tau)/(Kcc*theta)) ,'Ti= ',(((1.35+0.25*(theta/tau))/(0.54+0.33*(theta/tau)))*theta) , 'Td= ', ((0.5*theta) /(1.35+0.25*(theta/tau)))
, 'FI = ', (theta/tau))
print ("--------------------------------------------------------")

T4, out_pid4 = ctl.step_response(G44_s,Tsim);

#--------------------------------------------------------------------------------- #

# Sintonia Via Matlab Resposta ao Degrau, curva em "S"

Kp1 = 0.876 ; ki1 = 0.02 ; kd1 = 0.0;

# FT que o Matlab sintoniza o PID: Kp + Ki/s + Kd*s

G3 = (Kp1 + ki1/s + kd1*s);
G3_s = ctl.series(G3,ps_tf);
G33_s=ctl.feedback(G3_s, H_s, sign= -1);

print ("-------------Sintonia de PID - Matlab Resp ao Degrau----------")
print('Kp = ', Kp1,'Ti = ',ki1, 'Td= ',kd1)
print ("--------------------------------------------------------")
T3, out_pid3 = ctl.step_response(G33_s,Tsim);

# --------------------------------------------------------------------------------- #

#Melhor Sintonia Encontada via refino da sugestão do MATLAB: 
    
print ("-------------Sintonia de PID - Final----------")
#print('Kp = ', 0.3912,'Ti = ',0.0632, 'Td= ',0.1837)
print('Kp = ', 0.8,'Ti = ',0.086, 'Td= ',0.020)
print ("--------------------------------------------------------")   

#--------------------------   Plot dos Gráficos: ------------------------------# 
    
#---- Primeiro Método de ZN

#Plot da Figura utiliza para ZN - Resp ao degrau
plt.figure(1)
 
plt.plot (T, out_pid,'b-', label=" Resp. ao Degrau da FT da planta.");
plt.plot(T2, out_pid2, color = 'green',label="Resp. ao degrau usando ZN - primeiro método.");
plt.plot(t, degrau, 'y-', label = 'Entrada Degrau.')
#plt.plot(t, degrau, 'y-', label = 'Entrada Degrau')
plt.legend()
plt.xlabel('Tempo')
plt.ylabel('Amplitude')
plt.title('Aplicação do método de ZN, primeiro método.')
plt.grid(True)
plt.show()
 
#---- Segundo Método de ZN
    
#Plot da Figura utiliza para ZN - Resp em Frequencia
plt.figure(2)
plt.plot (T, out_pid,'b-', label=" Resp. ao Degrau da FT da planta.");
plt.plot (T1,out_pid1,'r-', label=" Resp. usando ZN - segundo método.");
plt.plot(t, degrau, 'y-', label = 'Entrada Degrau')
#plt.plot(t, degrau, 'y-', label = 'Entrada Degrau')
plt.legend()
plt.xlabel('Tempo')
plt.ylabel('Amplitude')
plt.title('Aplicação do método de ZN, segundo método.')
plt.grid(True)
plt.show()

#---- Método de CC

#Plot da Figura utiliza para CC - Resp ao degrau
plt.figure(3)

plt.plot (T, out_pid,'b-', label=" Resp. ao Degrau da FT da planta.");
plt.plot (T4,out_pid4,color = 'orange', label=" Resp. usando CC.");
plt.plot(t, degrau, 'y-', label = 'Entrada Degrau.')
#plt.plot(t, degrau, 'y-', label = 'Entrada Degrau')
plt.legend()
plt.xlabel('Tempo')
plt.ylabel('Amplitude')
plt.title('Aplicação do método de CC.')
plt.grid(True)
plt.show()


#---- Método de sugerido pelo MATLAB

#Plot da Figura utiliza para sintoniai via MATLAB - Resp ao degrau
plt.figure(4)
plt.plot (T, out_pid,'b-', label=" Resposta ao Degrau da FT da planta.");
plt.plot (T3,out_pid3,color = 'orange', label=" Resp. usando MATLAB.");
plt.plot(t, degrau, 'y-', label = 'Entrada Degrau.')
plt.legend()
plt.xlabel('Tempo')
plt.ylabel('Amplitude')
plt.title('Sugestão de sintonia obtida pelo Toolbox Control utilizando o MATLAB.')
plt.grid(True)
plt.show() 