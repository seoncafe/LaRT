#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

#fname = 'rin0.1_Vrot100_NHI18'
fname = 'rin0.1_Vrot300_NHI18'

#fname = 'rin0.1_Vrot100_NHI18_unif'
#fname = 'rin0.1_Vrot300_NHI18_unif'

c  = fits.open(fname+'.fits.gz')
c1 = fits.open(fname+'_obs_001.fits.gz')
c2 = fits.open(fname+'_obs_002.fits.gz')
c3 = fits.open(fname+'_obs_003.fits.gz')
c4 = fits.open(fname+'_obs_004.fits.gz')

vel   = c[1].data['velocity']
spec1 = np.sum(c1[0].data, axis=(0,1))
spec2 = np.sum(c2[0].data, axis=(0,1))
spec3 = np.sum(c3[0].data, axis=(0,1))
spec4 = np.sum(c4[0].data, axis=(0,1))
#spec1 = np.sum(c1[0].data, axis=(0,1)) + np.sum(c1[1].data, axis=(0,1))
#spec2 = np.sum(c2[0].data, axis=(0,1)) + np.sum(c2[1].data, axis=(0,1))
#spec3 = np.sum(c3[0].data, axis=(0,1)) + np.sum(c3[1].data, axis=(0,1))
#spec4 = np.sum(c4[0].data, axis=(0,1)) + np.sum(c4[1].data, axis=(0,1))

beta1 = c1[0].header['beta']
beta2 = c2[0].header['beta']
beta3 = c3[0].header['beta']
beta4 = c4[0].header['beta']

plt.plot(vel, spec1, label=r'$\beta=%d^\circ$' % beta1)
plt.plot(vel, spec2, label=r'$\beta=%d^\circ$' % beta2)
plt.plot(vel, spec3, label=r'$\beta=%d^\circ$' % beta3)
plt.plot(vel, spec4, label=r'$\beta=%d^\circ$' % beta4)

plt.xlim(-500.0, 500.0)
plt.xlabel(r'Velocity [km s$^{-1}$]')
plt.ylabel(r'Spectrum')

plt.legend()
plt.show()
