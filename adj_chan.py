import os
import matplotlib.pyplot as plt
import fitgauss as fit
import numpy as np
"""
------------------------------------------------------------
Metodo per trovare il fondo con i canali adiacenti al picco.
------------------------------------------------------------
"""

NOME_BCKG = 'fondo_54437.txt'
PATH_BCKG = os.path.join('C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/spettri/', NOME_BCKG)

background = np.loadtxt(PATH_BCKG, skiprows=12, max_rows=2048, unpack=True)

NOME_SPETTRO = fit.NOME_SPETTRO

PATH = 'C:/Users/Lorenzo/Desktop/Lab/Spettroscopia/spettri'
PATH = os.path.join(PATH, NOME_SPETTRO)
channels = fit.channels
counts = fit.counts

channels1 = fit.channels1
counts1 = fit.counts1

NOME_SPETTRO = NOME_SPETTRO.replace('_1.txt', '')

#live time di acquisizione del background e dei radionuclidi
LIVE_TIME_BCKG = 54437
LIVE_TIME = {'Am241': 271, 'Ba133': 233,
             'Co60': 344, 'Na22': 1517, 'Cs137': 194}
#LIVE_TIME_Cs3 = 1380

#il background è riscalato con il rapporto fra il live time della
#misura dello spettro ed il live time della misura del fondo senza
#sorgente
background = background * LIVE_TIME[NOME_SPETTRO]/LIVE_TIME_BCKG
#background = background * LIVE_TIME_Cs3/LIVE_TIME_BCKG
init_values = fit.init_values

F = fit.FitGauss(channels1, counts1, init_values)
risultati = fit.risultati(F)
mu0 = risultati[0]
sigma0 = risultati[1]
A0 = risultati[2]
B0 = risultati[3]
dm, dsigma, dA, dB = np.sqrt(F.covm.diagonal())


#canali a 3 sigma dal picco, trovato con il fit
A = int(np.floor(mu0 - 3*sigma0))
B = int(np.floor(mu0 + 3*sigma0))

#numero di canali da mediare prima di A e dopo B
m = 5

#media dei conteggi degli m canali prima di A e dopo B
mu_i = 0
mu_f = 0
for i in range(m+1):
    mu_i += counts[A-i]
    mu_f += counts[B+i]
mu_i = mu_i/m
mu_f = mu_f/m

mu_i = int(np.floor(mu_i))
mu_f = int(np.floor(mu_f))

#media degli m canali prima di A e dopo B
chan_i = A
chan_f = B
for i in range(m):
    chan_i += - 1/m
    chan_f +=  1/m

chan_i = int(np.floor(chan_i))
chan_f = int(np.floor(chan_f))

AREA_TOT = 0
for i in range(chan_i, chan_f):
    AREA_TOT += counts[i]

VAR_AREA_TOT = AREA_TOT
SIGMA_AREA_TOT = np.sqrt(VAR_AREA_TOT)

print(f'Area totale = {AREA_TOT:.3f} +- {SIGMA_AREA_TOT:.3f}\n')

N = abs(int(chan_i - chan_f))
AREA_CONTINUUM = 0.5*(mu_i + mu_f)*N
VAR_FONDO = AREA_CONTINUUM * N/(2*m)
SIGMA_FONDO = np.sqrt(VAR_FONDO)

BCKG = 0
for i in range(chan_i, chan_f):
    BCKG += background[i]

#area ottenuta sottraendo continuum e fondo ai dati
AREA_NETTA = AREA_TOT - BCKG - AREA_CONTINUUM
SIGMA_NETTA = np.sqrt(VAR_AREA_TOT + BCKG + VAR_FONDO)

def retta(x, y, A, B):
    """Interpolazione lineare fra i punti A e B del continuum."""
    return ((x - A)/(B - A))*(y[B] - y[A]) + y[A]

x = np.linspace(chan_i, chan_f, N)
retta = retta(x, counts, chan_i, chan_f)

channels_restrict = np.array([])#canali ristretti fra A e B, per plottare la retta
counts_restrict = np.array([])#conteggi ristretti fra A e B, per il fit allo spettro netto
spettro_netto = np.array([])#spettro a cui si sottraggono fondo e continuum


for i in range(chan_i, chan_f):
    spettro_netto = np.append(spettro_netto, counts[i] - retta[i - chan_i] - background[i])

for i in range(chan_i, chan_f):
    channels_restrict = np.append(channels_restrict, channels[i]) 
    counts_restrict = np.append(counts_restrict, counts[i])

init_values = [mu0, sigma0, A0, 0]
F1 = fit.FitGauss(channels_restrict, spettro_netto, init_values)
risultati1 = fit.risultati(F1)
mu1 = risultati1[0]
sigma1 = risultati1[1]
A1 = risultati1[2]
B1 = risultati1[3]

dm1, dsigma1, dA1, dB1 = np.sqrt(F1.covm.diagonal())
FWHM = 2.35*sigma1

print(f'Area continuum = {AREA_CONTINUUM:.3f} +- {SIGMA_FONDO:.3f}\n')
print(f'Area netta = {AREA_NETTA:.3f} +- {SIGMA_NETTA:.3f}\n')
#print(f'media = {mu0:.3f} +- {dm:.3f}\n')
#print(f'sigma = {sigma0:.3f} +- {dsigma:.3f}\n')
#print(f'A = {A0:.3f} +- {dA:.3f}\n')
#print(f'B = {B0:.3f} +- {dB:.3f}\n')
print(f'media netta = {mu1:.3f} +- {dm1:.3f}\n')
print(f'sigma netta = {sigma1:.3f} +- {dsigma1:.3f}\n')

plt.title(NOME_SPETTRO + ' ' + 'con continuum e background')
plt.xlabel('Channels [UA]')
plt.ylabel('Counts [UA]')
plt.plot(channels, counts, marker='o', color='b', label='Dati')
plt.plot(channels, background, marker = 'o', label = 'Fondo true')
plt.plot(x, retta, marker = 'o', label = 'Retta fra i punti a $\pm$ $3\sigma$ dal picco')
plt.plot(channels_restrict, spettro_netto, marker = 'o', label = 'Spettro al netto di fondo e continuum')
plt.minorticks_on()
plt.legend()
plt.show()
