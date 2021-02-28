import numpy as np
from scipy.fftpack import fft
from scipy.optimize import curve_fit
from astropy.time import Time
from astropy.utils import iers
from matplotlib import pyplot as plt

iers_b = iers.IERS_B.open(iers.IERS_B_URL)

mjd = np.arange(37666, 59217)
t = Time(mjd, format='mjd', scale='tt')
xp, yp = iers_b.pm_xy(t)
tt, xp, yp = t.tt.mjd, xp.to_value(), yp.to_value()

plt.plot(xp, yp)
plt.axis('square')
plt.show()

f_xp, f_yp = abs(fft(xp)), abs(fft(yp))

plt.plot(np.arange(200), f_xp[:200], c='green', linestyle=':')
plt.plot(np.arange(200), f_yp[:200], c='black', linestyle=':')
plt.show()

def func(t, a, b, A1, w1, p1, A2, w2, p2):
    return a*t+b+A1*np.sin(w1*t+p1)+A2*np.sin(w2*t+p2)

lbounds = [-np.inf, -np.inf,
           -np.inf, 2*np.pi/((mjd.max()-mjd.min())/49), -np.inf,
           -np.inf, 2*np.pi/((mjd.max()-mjd.min())/58), -np.inf]
ubounds = [+np.inf, +np.inf,
           +np.inf, 2*np.pi/((mjd.max()-mjd.min())/51), +np.inf,
           +np.inf, 2*np.pi/((mjd.max()-mjd.min())/60), +np.inf]
popt_xp, pcov_xp = curve_fit(func, tt, xp, bounds=(lbounds, ubounds))
popt_yp, pcov_yp = curve_fit(func, tt, yp, bounds=(lbounds, ubounds))
print(popt_xp)
print(pcov_xp)
print(popt_yp)
print(pcov_yp)

plt.plot(tt, xp)
plt.plot(tt, func(tt, *popt_xp))
plt.show()

plt.plot(tt, yp)
plt.plot(tt, func(tt, *popt_yp))
plt.show()
