#importando as bibliotecas
import control.matlab as ctl
import matplotlib.pyplot as plt
import numpy as np
import math

#fechando gráficos abertos
plt.close('all')

#1-
#criando a função de transferência
numGp = np.array([13500])
denGp = np.array([1,-120,3600])
Gp = ctl.tf(numGp,denGp)

#gerando a FTMF
FTMF = ctl.feedback(Gp,1)
print('FTMF sem controlador: ',FTMF)

#fazendo a resposta ao degrau
TQL, t = ctl.step(FTMF)

plt.figure()
plt.title('Resposta ao degrau sem compensação')
plt.plot(t,TQL,label='SMF')
plt.plot([0,0.0001,0.1],[0,1,1],label='Referência')
plt.grid()
plt.ylabel('Tensão proporcional à quantidade de levitação')
plt.xlabel('Tempo')
plt.legend(loc='best')
plt.show()

#2-
#Especificações do projeto
Mp = 0.2
print('Mp = ',Mp)
ts = 1
print('ts = ',ts,' s')
zeta = abs(math.log(Mp))/math.sqrt((math.pi**2)+(math.log(Mp)**2))
print('zeta = ',zeta)
wn = 4/(zeta*ts)
print('wn = ',wn,' rad/s')

sigma_i = wn*zeta
wi = math.sqrt((wn**2)-(sigma_i**2))
print('si = -',sigma_i,'+j',wi)

#Projetando o controlador PID

# Avaliação da condição angular
alfaID = math.atan(wi/sigma_i)*180/math.pi
theta_PPID = 180 - alfaID

theta_P0 = math.atan(wi/((60-sigma_i)))*180/math.pi

alfa1 = math.atan(wi/((60+sigma_i)))*180/math.pi
theta_P1 = 180 - alfa1

theta_z = (theta_PPID+theta_P1+theta_P0)-180
theta_z1 = theta_z/2
theta_z2 = theta_z - theta_z1
print('theta z1 = ',theta_z1)
print('theta z2 = ',theta_z2)

#Definindo os zeros
x1 = wi/(math.tan(theta_z1*math.pi/180))
z1 = -sigma_i - x1
print('z1 = ',z1)

x2 = wi/(math.tan(theta_z2*math.pi/180))
z2 = -sigma_i - x2
print('z2 = ',z2)

zero1 = -z1
zero2 = -z2

#Definindo os ganhos
ApPID = math.sqrt((sigma_i)**2+(wi)**2)
Ap0 = math.sqrt((60-sigma_i)**2+(wi)**2)
Ap1 = math.sqrt((60+sigma_i)**2+(wi)**2)
Az1 = math.sqrt((z1)**2+(wi)**2)
Az2 = math.sqrt((z2)**2+(wi)**2)

#Definindo os K's do controlador
Kd = (ApPID*Ap1*Ap0)/(Az1*Az2*13500)
print('Kd = ',Kd)
Kp = Kd*(zero1+zero2)
print('Kp = ',Kp)
Ki = Kd*(zero1*zero2)
print('Ki = ',Ki)


# Controlador PID projetado
C=ctl.tf(np.array([Kd,Kd*(zero1+zero2),Kd*zero1*zero2]),np.array([1,0])) 
print('Controlador PID',C)


#Gerando o SMF compensado
FTMAcomp = ctl.series(C,Gp)
FTMFcomp = ctl.feedback(FTMAcomp,1)

#3-
#Gerando o gráfico do Root Locus
rlist, klist = ctl.rlocus(FTMAcomp, print_gain=True, grid=True)
plt.grid()
plt.show()

#4-
#Gerando o gráficocom os três sinais
TQL_comp, t_comp = ctl.step(FTMFcomp)

plt.figure()
plt.title('Resposta ao degrau')
plt.plot(t_comp,TQL_comp,label='SMFcomp')
plt.plot([0,0.0001,2.5],[0,1,1],label='Referência')
plt.plot(t,TQL,label='SMFuncomp')
plt.grid()
plt.ylabel('Tensão proporcional à quantidade de levitação')
plt.xlabel('Tempo')
plt.legend(loc='best')
plt.show()

#Pegando os parâmetros do sistema compensado para o relatório
print('\nParâmetros do sistema compensado: ',ctl.stepinfo(FTMFcomp))
