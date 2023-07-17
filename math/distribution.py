import numpy as np

def gaussian1D(x, mu, sigma, amplitude):
    y = (amplitude / (sigma*np.sqrt(2*np.pi))) * np.exp(-0.55 * ((x-mu) / sigma)**2)
    return y