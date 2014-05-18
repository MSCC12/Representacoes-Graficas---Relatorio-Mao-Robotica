Representacoes-Graficas---Relatorio-Mao-Robotica
================================================
# Relatório (Report) - Projecto 'Mão Robótica'

## Objectivo (Goal)

No contexto do desenvolvimento do projecto da tese de mestrado de um dos elementos do grupo, recorrendo à utilização do sistema Bitalino, procedeu-se à aquisição do sinal electromiográfico gerado pela contracção isométrica do músculo braquiorradial do antebraço.
Através dos dados extraídos do sinal adquirido pelo sistema em aquisição, pretende-se a extrapolação e modelação posterior dos mecanismos de contracção subsequentes ao movimento da mão.


## Equipa (Team)

Melissa Cruz
Joana Castro


## Processo (Process)

Para aquisição do sinal pretendido e anteriormente enunciado, procedeu-se à execução dos passos cuja descrição é apresentada seguidamente.


### Step 1 - Primeiro sinal adquirido (First sinal recorded)

O primeiro sinal adquirido, disponibilizado através do link que se segue abaixo.

https://cloud.githubusercontent.com/assets/7483292/2949655/214d5e30-da11-11e3-9d4d-82cc7e09fb94.jpg

Considerando, através de uma primeira análise, que os dados obtidos são compatíveis com o tratamento necessário para aplicação no desenvolvimento do projecto em discussão, não se procederam a mais aquisições de sinal, tendo-se iniciado então o processamento do mesmo.


### Step 2 - Sinal Melhorado (Improved signal)

Para melhoramento do sinal originalmente obtido, isto é, para supressão da interferência de fontes de ruído, procedeu-se à implementação dos seguintes passos:
- Implementação de um filtro Notch, para remoção do ruído associado à ligação do sistema de aquisição à rede eléctrica nacional (banda de remoção de frequências entre os 49,5 Hz e os 50,5 Hz);
- Implementação de um filtro Passa-Alto, para remoção de ruído associado a interferência do sinal electrocardiográfico;
- Implementação de um filtro Passa-Baixo, para remoção das componentes DC-offset e ruído associado a movimento;

As representações das três implementações consecutivas dos filtros anteriormente mencionados, são disponibilizadas através do link que seguidamente se apresenta.

https://cloud.githubusercontent.com/assets/7483292/2949636/f76fafbe-da10-11e3-960d-d0ba109a0a14.jpeg


### Step 3 - Definição dos passos a implementar para processamento do sinal (Defining the signal processing stages)

Para processamento do sinal adquirido, foram aplicadas ferramentas para:
- Análise estatística do sinal (valor máximo, médio, mínimo e desvio padrão do sinal);
- Análise de amplitude (valores de amplitude do sinal, valores dos picos do sinal, valores dos picos máximo e mínimo e momentos no tempo em que estas ocorrências se verificaram)
- Análise do início da contracção muscular (determinação do momento de ocorrência do primeiro pico registado no sinal)
- Análise em frequência (determinação da taxa de amostragem e histograma de frequências do sinal)

### Step 4 - Código de processamento

O código aplicado no processamento de sinal consta do exemplo abaixo apresentado.

from pylab import *
from scipy.signal import filtfilt
from scipy import signal
from numpy import array, clip, argsort, sort
from pylab import find
from scipy.signal import argrelmax

data = loadtxt('EMG_mao_robotica.txt')

EMGb=data[:,5]

figure()
subplot(2,1,1)
plot(EMGb)
xlim([5000,20000])
ylim([270,700])
title('Sinal Electromiografico Original')
xlabel('Sample #')
ylabel('ADC')

t = arange(len(EMGb))/1000.

nbits=10
Vcc=3.3

G=1000.

EMGv = (EMGb/(2**nbits-1)-.5)*Vcc/2./G
EMGmv = (EMGv-EMGv.mean())*G

subplot(2,1,2)
plot(t, EMGmv)
xlim([0,30])
ylim([-0.45,0.55])
title('Sinal Electromiografico Original')
xlabel('Tempo(ms)')
ylabel('Tensao(mV)')
show()
"""savefig('EMG_Mão_Robótica_Sinal_Original.jpg')"""


sinal=EMGmv
fs=1000.0
order=2
t=arange(len(EMGb))/1000.

figure()
subplot(4,1,1)
f1=49.5
f2=50.5
b,a=signal.butter(order, [f1 * 2 / fs, f2 * 2 / fs], btype='bandstop') 
s1=filtfilt(b, a, EMGmv)   
plot(t,s1)
xlim([0,30])
ylim([-0.45,0.55])
title('Implementacao de Filtro Notch')
xlabel('Tempo(s)')
ylabel('Tensao(mV)')
    
subplot(4,1,2)
f=1000
b,a=signal.butter(order, f/fs)
s2=filtfilt(b, a, s1)
plot(t,s2)
xlim([0,30])
ylim([-0.45,0.55])
title('Implementacao de Filtro Passa-Baixo')
xlabel('Tempo(s)')
ylabel('Tensao(mV)')
      
subplot(4,1,3)
f=12.5
b,a=signal.butter(order, f * 2 / fs, btype='highpass')
s3=filtfilt(b, a, s2)
plot(t,s3)
xlim([0,30])
ylim([-0.45,0.55])
title('Implementacao de Filtro Passa-Alto')
xlabel('Tempo(s)')
ylabel('Tensao(mV)')

"""savefig('EMG_Mão_Robótica_Sinal_Melhorado.jpg')"""

s4=abs(s3)
subplot(4,1,4)
plot(t,s4)
xlim([0,30])
ylim([-0.45,0.55])
title('Sinal Rectificado')
xlabel('Tempo(s)')
ylabel('Tensao(mV)')

print "Descrição estatística do sinal filtrado:\n"
print "- Valor máximo do sinal filtrado: " + str(max(s4))
print "- Valor mínimo do sinal filtrado: " + str(min(s4))
print "- Desvio padrão do sinal filtrado: " + str(std(s4))
print "- Valor médio do sinal filtrado: " + str(mean(s4)) + "\n"

amplitude_pp=(max(s4)-min(s4))

print "Análise de amplitude do sinal filtrado:\n"
print "- Valor máximo de amplitude pico a pico sinal filtrado: " + str(amplitude_pp)

def peaks(s4, tol=None):
    if (tol is None):
        tol = min(s4)
        pks = argrelmax(clip(s4, tol, s4.max()))
        return pks[0]
        
picos=peaks(s4)
t_amplitude_max=max(peaks(s4))
t_amplitude_med=mean(peaks(s4))
t_amplitude_min=min(peaks(s4))

print "- Momentos em que ocorrem picos de amplitude do sinal filtrado: "+ str(picos) + "\n"

print "- Momento, em ms, em que ocorre o valor máximo de amplitude pico a pico sinal filtrado: " + str(t_amplitude_max) + "\n"

amp_max=s4[31196]
print " Valor máximo de amplitude pico a pico sinal filtrado: " + str(amp_max) + "\n"

print "- Momento, em ms, em que ocorre o valor médio de amplitude pico a pico sinal filtrado: " + str(t_amplitude_med ) + "\n"

amp_med=s4[15543.3585046]
print " Valor médio de amplitude pico a pico sinal filtrado: " + str(amp_med) + "\n"

print "- Momento, em ms, em que ocorre o valor mínimo de amplitude pico a pico sinal filtrado: " + str(t_amplitude_min )+ "\n"

amp_min=s4[8]
print " Valor mínimo de amplitude pico a pico sinal filtrado: " + str(amp_min) + "\n"

contracao=t[peaks(s4, tol=None)[1]]
print "- Momento, em ms, em se iniciou a contracçao muscular " + str(contracao) + "\n"

taxaAmostragem= mean(1/diff(t))
print "- Taxa de amostragem do sinal (Hz):" + str(taxaAmostragem) + "\n"


figure()
hist(s4,100,orientation='vertical')
title('Histograma de frequências')
xlabel('Frequências')
show()


### Step 5 - Conclusões (Final considerations)

Através da realização do presente trabalho, foi possível a compreensão da necessidade da aplicação de ferramentas de processamento específicas ao processamento de sinal electromiográfico para que assim seja possível a extracção de informação pertinente que permita a avaliação dos processos de contracção, relaxamento e fatiga muscular, e respectivas diferenças, em indivíduos saudáveis ou afectados por patologias do foro muscular, ou entre outras características cuja relação se pretenda estabelecer entre as mesmas e os processos anteriormente mencionados. 
