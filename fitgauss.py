import logging
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

"""
======================================================================================================
Fit gaussiano nell'intorno dei picchi, per ricavare media (posizione del canale), deviazione standard,
ampiezza ed offset.
======================================================================================================
"""

PATH = 'C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/spettri'  # percorso dei file .txt
NOME_SPETTRO = 'Am241_1.txt'  # modificare con il nome del file
PATH = os.path.join(PATH, NOME_SPETTRO)

# mette i risultati del fit nel file NOME_SPETTROlog.txt
# logging.basicConfig(filename=NOME_SPETTRO.replace('txt', '_bckg_log.txt'),
# level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')

# salta i commenti ed acquisice i conteggi dei canali 0-2047
counts = np.loadtxt(PATH, skiprows=12, max_rows=2048, unpack=True)
channels = np.array([i for i in range(0, 2048)],
                    dtype=float)  # numero di canali

#ch1 Am: 69-110, Ba: 34-65, 98-147, 364-441, 441-545
#Co: 1420-1584, 1600-1790, Na: 600-750, 1550-1711 
#Cs_1: 790-940 

# canali vicino al picco, da n a n_max-1
channels1 = np.array([channels[i] for i in range(69, 110)], dtype=float)
counts1 = np.array([counts[i] for i in range(69, 110)], dtype=float)

PATH_LOGS = 'C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/logs'
NOME_LOG = NOME_SPETTRO.replace('.txt', '.log.txt')
PATH_LOGS = os.path.join(PATH_LOGS, NOME_LOG)
init_values = [np.loadtxt(PATH_LOGS, skiprows=1, usecols=(0, 2, 4, 6),
                          max_rows=1, unpack=True)]

def gaussiana(x, mu, sigma, A, B):
    """Funzione per fit gaussiano channels-counts. A è l'ampiezza della gaussiana
    B è l'offset."""
    return A*(1/(sigma*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-mu)/sigma)**2) + B


def log_results(channel, mu, dmu, sigma, dsigma, A, dA, B, dB):
    """Inserisce i risultati del fit nel file NOME_SPETTRO_bckg.log.txt ."""
    logging.info(f'Range di canali: {channel[0]}-{channel[-1]}\n')
    logging.info(f'media = {mu:.3f} +- {dmu:.3f}\n')
    logging.info(f'dev. std = {sigma:.3f} +- {dsigma:.3f}\n')
    logging.info(f'A = {A:.3f} +- {dA:.3f}\n')
    logging.info(f'B = {B:.3f} +- {dB:.3f}\n')


def plot_results(channels, counts, mu, sigma, A, B, NOME_SPETTRO):
    """Plotta i risultati del fit nell'intervallo counts1, channels1."""
    plt.plot(channels, gaussiana(channels, mu, sigma, A, B), color='red')
    plt.plot(channels, counts, marker='o')
    NOME_SPETTRO = NOME_SPETTRO.replace('_1.txt', '')
    plt.title('Channels vs counts' + ' ' + NOME_SPETTRO)
    plt.xlabel('Channels')
    plt.ylabel('Counts')
    plt.minorticks_on()
    plt.show()


class Fit_iterator:
    """Iteratore di FitGauss."""

    def __init__(self, classe):
        self.classe = classe
        self.i = 0

    def __next__(self):
        if self.i < len(self.classe.Fit()[0]):
            result = self.classe.Fit()[0][self.i]
            self.i += 1
            return result
        raise StopIteration


class FitGauss:
    """Classe per il fit gaussiano. Prende in input gli array di canali,
    conteggi e valori iniziali e restituisce i parametri ottimali e la matrice
    di correlazione."""

    def __init__(self, x, y, init):
        self.init = init
        self.x = x
        self. y = y
        self.pars = np.array([])
        self.covm = np.array([], [])

    def Fit(self):
        self.pars, self.covm = curve_fit(gaussiana, self.x, self.y, self.init)
        return self.pars, self.covm

    def __iter__(self):
        return Fit_iterator(self)


def risultati(F):
    """"Restituisce i risultati del fit, in ordine: media, sigma, A, B.
    Prende in input la classe per il fit."""
    elem = np.array([])
    iterator = iter(F)
    while True:
        try:
            # Get next element from TeamIterator object using iterator object
            elem = np.append(elem, [next(iterator)])
        except StopIteration:
            break
    return elem


if __name__ == '__main__':

    F = FitGauss(channels1, counts1, init_values)
    risultati = risultati(F)
    mu0 = risultati[0]
    sigma0 = risultati[1]
    A0 = risultati[2]
    B0 = risultati[3]
    dm, dsigma, dA, dB = np.sqrt(F.covm.diagonal())
    print(f'{A0} +- {dA}')
    print(f'{B0} +- {dB}')
    #log_results(channels1, mu0, dm, sigma0, dsigma, A0, dA, B0, dB)
    plot_results(channels1, counts1, mu0, sigma0, A0, B0, NOME_SPETTRO)
