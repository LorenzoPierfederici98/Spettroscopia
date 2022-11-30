import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#fit per trovare il coefficiente di attenuazione
#le aree sotto ai picchi sono state prese sottraendo il continuum ed il fondo naturale
#scalato con l'attenuazione del Pb.

A = np.array([3393, 1202, 1278, 852, 732])
dA = np.array([275, 164, 200, 205, 222])

x = np.array([0, 2*0.150, 3*0.150, 4*0.150, 5*0.150])
dx = np.array([0, 2*0.005, 3*0.005, 4*0.005, 5*0.005])

t_live = np.array([1207, 544, 594, 713, 808])
A = A/t_live #intensità conteggi/s
dA = dA/t_live
r_Al = 2.70 #densità alluminio g/cm^3
r_Pb = 11.34 #densità piombo g/cm^3

mu_Al = 0.073215 #coeff. att. alluminio cm^2/g
mu_Pb = 0.10675 #coeff. att. piombo cm^2/g

t = r_Pb*x #spessori massici g/cm^2
dt = r_Pb*dx
t_Pb = 5.17*r_Pb #spessore massico del piombo g/cm^2

def esponenziale(l, a, b):
    return a*np.exp(-b*l)

init_values = [A[0], mu_Pb]
pars, covm = curve_fit(esponenziale, t, A, init_values, dA)
a0, b0 = pars
da, db = np.sqrt(covm.diagonal())

dy = np.sqrt(dA**2 + (-b0*a0*np.exp(-b0*t))**2*dt**2)
pars1, covm1 = curve_fit(esponenziale, t, A, init_values, dy)
a1, b1 = pars1
da1, db1 = np.sqrt(covm1.diagonal())
chisq = (((A - esponenziale(t, *pars1))/dy)**2).sum()
ndof = len(A) - 2
print(f'chisq/ndof = {chisq:.3f}/{ndof}\n')
print(f'mu_Pb = {b1:.3f} +- {db1:.3f}\n')

plt.errorbar(t, A, dA, dt, marker = 'o', linestyle = '')
xx = np.linspace(t[0], t[-1], 1000)
plt.plot(xx, esponenziale(xx, *pars))
plt.xlabel('Spessori massici [$\dfrac{g}{cm^2}$]')
plt.ylabel('Intensità [$\dfrac{conteggi}{s}$]')
plt.minorticks_on()
plt.show()