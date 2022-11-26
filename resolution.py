import numpy as np
import logging
import linearfit as linear
import adj_chan as adj


"""
Calcolo della risoluzione energetica dalla media e deviazione standard date
dal fit gaussiano allo spettro al netto di fondo e continuum, in adj_chan.
La media trovata viene convertita in energia [keV] attraverso la calibrazione,
data dal fit in linearfit.
"""

NOME_SPETTRO = adj.NOME_SPETTRO + '_res.txt'
logging.basicConfig(filename=NOME_SPETTRO,
                    level=logging.INFO, format='%(asctime)s:%(levelname)s:%(message)s')

FWHM = adj.FWHM
sigma = adj.sigma1
dsigma = adj.dsigma1

pars = linear.pars
covm = linear.covm
da, db = np.sqrt(covm.diagonal())
correlation = covm[0][1]/(da*db)  # correlazione fra a e b

# canali che si vogliono convertire in energia, da usare per il calcolo della risoluzione
# è il mu1 di adj_chan, media dello spettro al netto di fondo e continuum
x = adj.mu1
dx = sigma

ener_calibration = pars[0]*x + pars[1]  # energia ottenuta dalla calibrazione
dener_calibration = np.sqrt(
    covm[0][0]*x**2 + covm[1][1] + 2*x*correlation*da*db)  # incertezza sull'energia
# considerando la correlazione fra i parametri
#RISOLUZIONE = FWHM/ener_calibration
#DRISOLUZIONE = 2.35*np.sqrt((dsigma/ener_calibration) **
#                            2 + (sigma/(ener_calibration**2))**2*dener_calibration**2)

RISOLUZIONE = FWHM/x
DRISOLUZIONE = 2.35*np.sqrt((dsigma/x) **
                            2 + (sigma/(x**2))**2*dx**2)

NOME_SPETTRO = NOME_SPETTRO.replace('_res.txt', '')

logging.info(f'calibrazione energia {ener_calibration:.3f} +- {dener_calibration:.3f}\n')
logging.info(f'Risoluzione energetica {NOME_SPETTRO} = {100*RISOLUZIONE:.3f} +- {100*DRISOLUZIONE:.3f} %\n')

print(f'calibrazione energia {ener_calibration:.3f} +- {dener_calibration:.3f}\n')
print(f'Risoluzione energetica {NOME_SPETTRO} = {100*RISOLUZIONE:.3f} +- {100*DRISOLUZIONE:.3f} %\n')
