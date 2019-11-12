#E.M. Vavagiakis 11/12/19
#Adapted from E. Schaan 11/01/19
#From Kravtsov 2014: 1401.7329

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import pandas as pd

#Variables from Kravtsov Appendix A Table 3
log10M1 = 11.39
log10e = -1.685
alpha = -1.740   # minus sign missing in Table 3
delta = 4.335
gamma = 0.531

def f(x):
      """Fitting function, Kravtsov+14, eq A4
      """
      result = np.log10(1.+np.exp(x))**gamma
      result *= delta
      result /= 1. + np.exp(10.**(-x))
      result += -np.log10(10.**(alpha*x) + 1.)
      return result

def fmStar(mVir):
      """Computes stellar mass [M_sun]
      from halo mass Mvir [Msun].
      """
      result = log10e + log10M1
      result += f(np.log10(mVir) - log10M1)
      result -= f(0.)
      result = 10.**result
      return result

#To get function to interpolate 
mVir_sample = np.logspace(np.log10(1.e5), np.log10(1.e20), 501, 10.)
mStar_sample = np.array(map(fmStar, mVir_sample))  # [M_sun]

#Interpolating
fmStarTomVir = interpolate.interp1d(mStar_sample, mVir_sample, kind='cubic', bounds_error=True, fill_value=0.)
fmVirTomStar = interpolate.interp1d(mVir_sample, mStar_sample, kind='cubic', bounds_error=True, fill_value=0.)

#Loading DR15 data
data=pd.read_csv('DR15_actplanck_150and90div_PS5and8arcmin_catalog_wbestObjID_20190501_EMV_CUTS_20191009.csv', comment = '#')
lum=np.array(data["lum"])
lumselect=[3.0*lum[x] for x in range(len(lum))] #Factor of 3.0 converts luminosity to stellar mass as per Kravtsov

massselect=fmStarTomVir(lumselect) #masses for the sample
print 'Length of new mass array:'
print len(massselect)

plt.plot(np.log10(lumselect),np.log10(massselect),'o',label='DR15, M*/L=3.0')
plt.xlabel('log10(M*)')
plt.ylabel('log10(Mass)')
plt.legend()
plt.show()

