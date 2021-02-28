import numpy as np
from scipy.optimize import curve_fit
from astropy.time import Time
from matplotlib import pyplot as plt

mjd = np.arange(41685, 59218)
t = Time(mjd, format='mjd', scale='tt')
tt, dt_tt_ut1 = t.tt.mjd, t.tt.mjd-t.ut1.mjd

def func(t, a, b, c, A, w, p):
    return a*t**2+b*t+c+A*np.sin(w*t+p)

lbounds = [-1e-12, 1e-9, -np.inf,
           -np.inf, 2*np.pi/8000, -np.inf]
ubounds = [-1e-16, 1e-6, +np.inf,
           +np.inf, 2*np.pi/6000, +np.inf]
popt, pcov = curve_fit(func, tt, dt_tt_ut1, bounds=(lbounds, ubounds))
print(popt)
print(pcov)

plt.plot(tt, dt_tt_ut1)
plt.plot(tt, func(tt, *popt))
plt.show()
