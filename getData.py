# -*- coding: utf-8 -*-

"""
Universidade Tecnológica Federal do Paraná - UTFPR-CT 2023
Trabalho de Conclusão de Curso - TCC02
Aluno : Daniel Souza Soares 
Orientador: Victor Frencl ; Coorientador: Alexandre Tuoto. 
Código para Sintonia de PID - Planta Didática de Vazão e pH.
""" 
# Importando as bibliotecas para as funções no Python.

import Adafruit_ADS1x15 as serial
import RPi.GPIO as gpio
from time import time as now

#*** Entradas e saídas do Raspberry ***
gpio.setwarnings(False)
gpio.setmode(gpio.BOARD)
#Saída 1: pino 33
#Entrada 1: A0 do ADC
pino_saida1=33;

#Inicializa a saída
gpio.setup(pino_saida1,gpio.OUT)

saida1=gpio.PWM(pino_saida1,300)

saida1.start(0)


# FUNCOES ------------------------------------------

# amostra os dados sampleData(ADC)
def sampleData(ADC):
    data = []
    for i in range(4):
        #data.append(i)
        data.append(ADC.read_adc(i, gain=1))
        #data.append(randomNumber(0, 10000)) #remover esta linha depois
    return data

# calcula tempo que passou
def getElapsed(startTime):
    return now() - startTime

# calcula o tempo de aquisicao
def getSampleTime(timestamp, lastTimestamp):
    return timestamp - lastTimestamp

# escreve dados no arquivo
def writeData(file, data):
    j = len(data)
    for i in range(j):
        file.write(str(data[i]))
        if i < (j-1) :
            file.write(",")
    file.write("\n")

# --------------------------------------------------

# INICIALIZACOES -----------------------------------

ADC = serial.ADS1115()

# tempo entre medicoes (em segundos)
sampleTime  = 0

# marca o momento inicial da execucao do script
startTime = now()

# marca o momento da ultima amostragem
lastSampleTime = startTime

# objeto do arquivo
file = open("output.csv", "w")


intensity = []
changeTime = []
control = open("sigControl.data", "r")

lines = control.readlines()

for line in lines :
    parameters = line.split(",")
    changeTime.append(float(parameters[0]))
    intensity.append(float(parameters[1]))

control.close()

# Caso o usuário queira Ver o que está sendo printado
#print(intensity)
#print(waitTime)

# tempo total de aquisicao de dados (em segundos)
totalTime = changeTime[-1] + 1

# -------------------------------------------------

# MAIN --------------------------------------------

A3_min = 5300 
A3_max = 26500
A3_value = 0
data_list = []

lastTimestamp = getElapsed(startTime)
current = 0

while current < len(changeTime) :
    
    # Timestamp da leitura
    timestamp = getElapsed(startTime)
        
    # Realiza a leitura do ADC
    leitura = ADC.read_adc(3, gain=1)
    # Converte o valor do ADC
    A3_value = ((leitura-A3_min)/(A3_max-A3_min))*100
    
    # Salva os dados para gravar em um arquivo
    data = [timestamp, getSampleTime(timestamp, lastTimestamp), leitura, A3_value, intensity[current]]
    data_list.append(data)
    
    lastTimestamp = timestamp
    
    
    if timestamp > changeTime[current] :
        print(timestamp)
        saida1.ChangeDutyCycle(intensity[current])
        current = current + 1
#
print("Teste finalizado")

for dados in data_list:
    writeData(file, dados)

file.close()
# Ordem dos vetores salvos no arquivo "output.csv"
#[timestamp, sampleTime, ADC_A3, vazao(%), inversor(%)]
print("Escrita Finalizada")

