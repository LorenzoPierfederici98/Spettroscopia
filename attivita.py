import numpy as np
import matplotlib.pyplot as plt

A0 = 74000
T = 30.17 #tempo di dimezzamento in anni
s = 60*60*24*365 #secondi in un anno
T = T * s #tempo di dimezzamento in secondi
l = (1/T) * np.log(2) #costante di decadimento [1/s]
t = 17.75 #tempo trascorso dal 2/2005 in anni
t = t * s #t in secondi
A_theo = A0*np.exp(-l*t)

a = 2.540 #raggio del rivelatore in cm
da = 0.005
r = 19.9 #distanza rivelatore-sorgente
dr = 0.1
Omega = np.pi*a**2/r**2 #angolo solido sotteso dal rivelatore
dOmega = np.sqrt((2*Omega/a)**2*da**2 + (Omega/r)**2*dr**2)
epsilon = 0.25 #efficienza intrinseca di picco
t_live = 1380 #tempo live di acquisizione
BR = 0.851
N_net = 48367
dN = 464

A = N_net*4*np.pi/(t_live*Omega*epsilon*BR)
dA = np.sqrt((A/N_net)**2*dN**2 + (A/Omega)**2*dOmega**2)
print(f'Attività misurata: {A/1000:.3f} +- {dA/1000:.3f} kBq\n')
print(f'Attività prevista: {A_theo/1000:.3f} kBq')

plt.errorbar(t/s, A/1000, dA/1000, marker = 'o', linestyle = '', label = 'Attività misurata')
plt.plot(t/s, A_theo/1000, marker = 'o', linestyle = '', label = 'Attività prevista')
plt.xlabel('t dal 2/2005 [y]')
plt.ylabel('A [$kBq$]')
plt.legend()
plt.show()