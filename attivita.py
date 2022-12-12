import numpy as np
import matplotlib.pyplot as plt

#misura di attività per lo spettro di Cs137_3
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
r = 19.9 #distanza rivelatore-sorgente in cm
dr = 0.1
Omega = 2*np.pi*(1 - r/np.sqrt(r**2 + a**2)) #angolo solido sotteso dal rivelatore
dOmega = 2*np.pi*np.sqrt((a**2/((np.sqrt(a**2 + r**2))**3))**2*dr**2 + (a*r/((np.sqrt(a**2 + r**2))**3))**2*da**2)

epsilon = 0.25 #efficienza intrinseca di picco
depsilon = 0.1*epsilon
t_live = 1380 #tempo live di acquisizione
BR = 0.851

N_net_continuum = 48367
dN_continuum = 464

N_net_fondo = 53847
dN_fondo = 250

print(f'Ang. solido = {Omega:.5f} +- {dOmega:.5f}\n')
A = N_net_continuum*4*np.pi/(t_live*Omega*epsilon*BR)
dA = np.sqrt((A/N_net_continuum)**2*dN_continuum**2 + (A/Omega)**2*dOmega**2 + (A/epsilon)**2*depsilon**2)
print(f'Attività misurata: {A/1000:.3f} +- {dA/1000:.3f} kBq\n')
print(f'Attività prevista: {A_theo/1000:.3f} kBq')

plt.errorbar(t/s, A/1000, dA/1000, marker = 'o', linestyle = '', label = 'Attività misurata\n')
plt.plot(t/s, A_theo/1000, marker = 'o', linestyle = '', label = 'Attività prevista\n')
plt.xlabel('t dal 2/2005 [y]')
plt.ylabel('A [$kBq$]')
plt.legend()
plt.minorticks_on()
plt.show()