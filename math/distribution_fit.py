import numpy as np
from scipy.optimize import least_squares
from .distribution import gaussian1D
from scipy.interpolate import interp1d

def gaussian1D_res(params, x0, y0):
    y = y0 - gaussian1D(x0, params[0], params[1], params[2])
    return y

def fit_gaussian1D(x, y, bounds=None):
    init_guess = [x[np.argmax(y)], interp1d(y[np.argmax(y):], x[np.argmax(y):])(0.5*np.amax(y)) - interp1d(y[:np.argmax(y)], x[:np.argmax(y)])(0.5*np.amax(y)), np.amax(y)]
    if bounds==None: bounds=(-np.inf, np.inf)
    result = least_squares(gaussian1D_res, init_guess, loss='cauchy', bounds=bounds, args=(x, y))
    mu, sigma, amplitude = [result.x[0], result.x[1], result.x[2]]
    residual = result.fun
    return {'mu':mu, 'sigma':sigma, 'amplitude': amplitude, 'residual':residual}