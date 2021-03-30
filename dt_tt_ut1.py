import numpy as np
from astropy.time import Time
from matplotlib import pyplot as plt

def mjd2y(mjd):
    y = 1961+(mjd-(37300+365.25/24))/365.25
    return y

def dt_1961_1986(y):
    dy = y-1975
    dt = 45.45+1.067*dy-3.84615e-3*dy**2-1.13927577e-3*dy**3
    dt = dt/86400
    return dt

mjd = np.arange(41685, 46434)
t = Time(mjd, format='mjd', scale='tt')
tt, dt_tt_ut1, dt = t.tt.mjd, t.tt.mjd-t.ut1.mjd, dt_1961_1986(mjd2y(mjd))
print(np.std((dt-dt_tt_ut1)*86400, ddof=1))
plt.plot(tt, dt_tt_ut1, tt, dt)
plt.show()
