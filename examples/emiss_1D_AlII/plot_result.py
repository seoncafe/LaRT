#!/usr/bin/env python
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../python')
from radial_profile import *

fname = 'AlII_ex'
a = fits.open(fname+'.fits.gz')
b = fits.open(fname+'_obs.fits.gz')

wavelength = a[1].data.wavelength

im1 = np.sum(b[0].data, axis=2)
im2 = np.sum(b[1].data, axis=2)
im = im1 + im2

r, prof = radial_profile(im)
rmax = 12.0
r    = r*rmax
prof = prof/np.amax(prof)

spec1 = np.sum(b[0].data, axis=(0,1))
spec2 = np.sum(b[1].data, axis=(0,1))
spec = spec1 + spec2
spec = spec/np.amax(spec)

plt.subplot(211)
plt.plot(r,prof)
plt.xlim(0.0, rmax)
plt.ylim(0.0, 1.1)
plt.xlabel('radius (pc)')
plt.ylabel('Normalized Intensity Profile')

plt.subplot(212)
plt.plot(wavelength,spec)
plt.xlim(np.amin(wavelength), np.amax(wavelength))
plt.ylim(0.0, 1.1)
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Normalized Spectrum')

plt.tight_layout()
plt.show()
