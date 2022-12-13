import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

"""
================================================================================================
Fit lineare per la calibrazione canali-energia, con le sorgenti di Am241, Ba133, Cs137_2, Co61, Na22
================================================================================================
"""


Energy, Channel, dChannel = np.loadtxt('C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/linear.txt', unpack = True)
dEnergy = np.zeros(len(Energy))

def linear_fit(x, a, b):
    return a*x + b

init_values = [1., 10.]
pars, covm = curve_fit(linear_fit, Channel, Energy, init_values, dChannel)
a0, b0 = pars
da, db = np.sqrt(covm.diagonal())

chisq = (((Energy - linear_fit(Channel, *pars))/dChannel)**2).sum()
ndof = len(Energy) - 2

if __name__ == '__main__':
    print(f'a = {a0:.3f} +- {da:.3f} keV/canale')
    print(f'b = {b0:.3f} +- {db:.3f} keV')
    print(f'chisq/ndof = {chisq:.3f}/{ndof}')
    print(covm)

    plt.errorbar(Channel, Energy, dEnergy, dChannel, marker = 'o', linestyle = '')
    x = np.linspace(Channel[0], Channel[-1], 100)
    plt.plot(x, linear_fit(x, *pars), color = 'r')
    plt.ylabel('Energy [KeV]')
    plt.xlabel('Channels')
    plt.show()
