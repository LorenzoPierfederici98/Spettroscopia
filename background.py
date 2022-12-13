import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

PATH = 'C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/spettri/'  # percorso dei file .txt
NOME_SPETTRO = 'Cs137_3.txt'  # modificare con il nome del file
PATH_SPETTRO = os.path.join(PATH, NOME_SPETTRO)

# salta i commenti ed acquisice i conteggi dei canali 0-2047
counts = np.loadtxt(PATH_SPETTRO, skiprows=12, max_rows=2048, unpack=True)
channels = np.array([i for i in range(0, 2048)],
                    dtype=float)  # numero di canali

NOME_BCKG = 'fondo_54437.txt'
PATH_BCKG = os.path.join(PATH, NOME_BCKG)

background = np.loadtxt(PATH_BCKG, skiprows=12, max_rows=2048, unpack=True)

NOME_SPETTRO = NOME_SPETTRO.replace('_3.txt', '')
#NOME_SPETTRO = NOME_SPETTRO.replace('_1.txt', '')

#live time di acquisizione del background e dei radionuclidi
LIVE_TIME_BCKG = 54437
LIVE_TIME = {'Am241': 271, 'Ba133': 233,
             'Co60': 344, 'Na22': 1517, 'Cs137': 194}
LIVE_TIME_Cs3 = 1380
#il background è riscalato con il rapporto fra il live time della
#misura dello spettro ed il live time della misura del fondo senza
#sorgente

#background = background * LIVE_TIME[NOME_SPETTRO]/LIVE_TIME_BCKG
background = background * LIVE_TIME_Cs3/LIVE_TIME_BCKG

def find_delta(channels, counts):
    """Trova la delta ottimale per il filtro gaussiano di convoluzione, ovvero quella per cui 
    la doppia convoluzione ha un minimo."""
    delta_min = 1
    delta_max = 40
    min_convolv = np.array([])

    for delta in range(delta_min, delta_max + 1):
        C1 = Iconvolution(channels, counts, delta)
        C2 = IIconvolution(channels, counts, delta, C1)
        min_convolv = np.append(min_convolv, [np.amin(C2)])
    return np.argmin(min_convolv) + 1


def Iconvolution(channels, counts, delta):
    """Prima convoluzione con la derivata prima della gaussiana."""
    convolved = np.zeros(len(counts))
    for j in range(len(counts)):
        convolved[j] = (-counts*((channels - j)/delta**2) *
                        np.exp(-((channels - j)**2/(2*delta**2)))).sum()
    return convolved


def IIconvolution(channels, counts, delta, convolved):
    """Seconda convoluzione della convoluzione precedente con la derivata prima della gaussiana."""
    convolved2 = np.zeros(len(counts))
    for j in range(len(counts)):
        convolved2[j] = (-convolved*((channels - j) /
                         delta**2)*np.exp(-((channels - j)**2/(2*delta**2)))).sum()
    return convolved2


def find_neg_index(delta, channels, counts):
    """Trova tutti gli indici neg_ind per cui la doppia convoluzione è negativa e maggiore in modulo 
    di una certa prominence (distanza dalla base), i picchi corrispondono
    ai minimi della doppia convoluzione."""
    convolved = Iconvolution(channels, counts, delta)
    convolved2 = IIconvolution(channels, counts, delta, convolved)
    [neg_ind, _] = find_peaks(-convolved2, prominence=100)
    return [neg_ind, _]


def lin_interpol(delta, channels, counts):
    """Interpolazione lineare per trovare gli zeri di channels-convolved2."""
    convolved = Iconvolution(channels, counts, delta)
    convolved2 = IIconvolution(channels, counts, delta, convolved)
    neg_ind = find_neg_index(delta, channels, counts)
    zc_indx = np.where(np.diff(np.sign(convolved2))
                       )  # indice prima di passare lo zero
    conv_zero = np.array([])  # array che contiene gli indici degli zeri
    for zc_i in zc_indx:
        t1 = channels[zc_i]
        t2 = channels[zc_i + 1]
        a1 = convolved2[zc_i]
        a2 = convolved2[zc_i + 1]
        conv_zero = np.append(
            (conv_zero), [(t1 + (0 - a1) * ((t2 - t1) / (a2 - a1)))])

    return conv_zero


def fondo(counts, sigma):
    """Metodo per trovare il background. m corrisponde alla FWHM del picco più largo, da trovare con
    l'algoritmo per i picchi."""
    m = int(np.floor(2.35*sigma))
    k = np.array([i for i in counts])
    z = np.zeros(len(counts))
    for p in range(1, m+1):
        for i in range(p-1, len(counts)-p):
            z[i] = min(k[i], 0.5*(k[i-p] + k[i+p]))
        for i in range(p-1, len(counts)-p):
            k[i] = z[i]
    return k


if __name__ == '__main__':

    delta = find_delta(channels, counts)
    print(f'Delta ottimale = {delta}\n')

    convolved = Iconvolution(channels, counts, delta)
    convolved2 = IIconvolution(channels, counts, delta, convolved)

    [neg_ind, _] = find_neg_index(delta, channels, counts)
    conv_zero = lin_interpol(delta, channels, counts)


    #calcolo della sigma_left e right per ogni picco, calcolate come 
    # sqrt(|(zero-picco)**2 - 2*delta_ottimale**2|)
    #sigma_left = np.array(
    #    np.sqrt(np.abs((conv_zero[:2*len(neg_ind):2] - neg_ind)**2 - 2*delta**2)))
    #sigma_right = np.array(
    #    np.sqrt(np.abs((conv_zero[1:2*len(neg_ind):2] - neg_ind)**2 - 2*delta**2)))
    
    conv_zero = conv_zero[:2*len(neg_ind)]

    #convolved2_zeros corrisponde a zero nei canali all'interno di conv_zero
    #è usato per misurare l'area sotto ai picchi
    convolved2_zeros = convolved2
    for i in range(len(conv_zero[:2*len(neg_ind)])):
        convolved2_zeros[int(conv_zero[i])] = 0

    #misura dell'area fra gli estremi di convolved2_zeros, per cui la
    #doppia convoluzione è nulla
    a = 0
    area = np.array([])
    i = 0
    while i <(2*len(neg_ind) - 1):
        for j in range(int(conv_zero[i]), int(conv_zero[i+1] + 1)):
            a += convolved2_zeros[j]
        a = -a
        area = np.append(area, a)
        a = 0
        i += 2

    #sigma è l'array in cui ci sono le deviazioni standard dei picchi 
    #viene stimato con la formula data dall'articolo
    #sigma_i = Area_i*sqrt(6*pi)/(9*conteggio di picco), le quantità sono
    #convertite in positivo
    sigma = np.array([])
    for i in range(len(neg_ind)):
        s = area[i]*(np.sqrt(6*np.pi))/(-9*convolved2[neg_ind[i]])
        sigma = np.append(sigma, s)

    print(f'{sigma}\n')
    #sigma_max che serve per trovare il fondo deve essere adattato alla larghezza
    #di ogni picco per trovare il fondo+continuum sotto allo stesso picco
    sigma_max = np.amax(sigma)
    #print(sigma_max)
    fondo = fondo(counts, sigma_max)
    print(f'Canali dei picchi: {neg_ind}\n')
    print(
        f'Zeri della funzione di doppia convoluzione: {np.floor(conv_zero[0:2*len(neg_ind)])}\n')

    net_counts = np.array([i for i in counts]) - fondo

    #print(f'sigma_left: {sigma_left}\n')
    #print(f'sigma_right: {sigma_right}\n')

    plt.title('Counts and background' + ' ' +
              NOME_SPETTRO + ' ' + f'$\delta$ = {delta}')
    plt.xlabel('Channels [UA]')
    plt.ylabel('Counts [UA]')
    #plt.plot(conv_zero, np.zeros(len(conv_zero)), 'o')
    plt.plot(neg_ind, convolved2[neg_ind], marker='P', label='Picchi')
    plt.plot(channels, convolved2, marker='o',
             label='Doppia convoluzione con $-\dfrac{x_{channel}-y}{\delta^2}e^{-\dfrac{(x_{channel}-y)^2}{2\delta^2}}$')
    plt.plot(channels, counts, marker='o', color='b', label='Dati')
    plt.plot(channels, fondo, marker='o', color='r', label='Fondo stimato')
    #plt.plot(channels, net_counts, marker='o', color='skyblue', label='Net counts')
    plt.plot(channels, background, marker = 'o', label = 'Fondo true')
    plt.minorticks_on()
    plt.legend()
    plt.show()
