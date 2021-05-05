import numpy as np
from scipy.optimize import curve_fit
from astropy.time import Time
from astropy.utils import iers
from matplotlib import pyplot as plt

iers_b = iers.IERS_B.open(iers.IERS_B_URL)

mjd = np.arange(37666, 59217)
t = Time(mjd, format='mjd', scale='tt')
xp, yp = iers_b.pm_xy(t)
tt, xp, yp = t.tt.mjd, xp.to_value(), yp.to_value()

'''
plt.plot(xp, yp)
plt.axis('square')
plt.show()
'''

def func(t, a, b, c, A1, w1, p1, A2, w2, p2):
    return a*t**2+b*t+c+A1*np.sin(w1*t+p1)+A2*np.sin(w2*t+p2)

lbounds = [-np.inf, -np.inf, -np.inf,
           -np.inf, 2*np.pi/440, -np.inf,
           -np.inf, 2*np.pi/370, -np.inf]
ubounds = [+np.inf, +np.inf, +np.inf,
           +np.inf, 2*np.pi/430, +np.inf,
           +np.inf, 2*np.pi/360, +np.inf]
popt_xp, pcov_xp = curve_fit(func, tt, xp, bounds=(lbounds, ubounds))
popt_yp, pcov_yp = curve_fit(func, tt, yp, bounds=(lbounds, ubounds))
print(popt_xp)
# print(pcov_xp)
print(np.sqrt(np.diag(pcov_xp)))
print(popt_yp)
# print(pcov_yp)
print(np.sqrt(np.diag(pcov_yp)))

plt.plot(tt, xp)
plt.plot(tt, func(tt, *popt_xp))
plt.plot(tt, np.polyval(popt_xp[0:3], tt))
plt.show()

plt.plot(tt, yp)
plt.plot(tt, func(tt, *popt_yp))
plt.plot(tt, np.polyval(popt_yp[0:3], tt))
plt.show()
