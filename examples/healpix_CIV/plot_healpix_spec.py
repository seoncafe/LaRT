#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import healpy as hp 

model_name = 'CIV_test'
a = fits.open(model_name + '_obs.fits.gz')
b = fits.open(model_name + '.fits.gz')

# wavelength
wavelength = b[1].data.wavelength

# scattered and direct light
I_scatt = a[0].data
I_direc = a[1].data

# number of pixels
npix = a[0].header['NAXIS2']

# integrate over frequency
spec_scatt = np.sum(I_scatt, 0)/npix
spec_direc = np.sum(I_direc, 0)/npix
spec_tot   = spec_scatt + spec_direc

plt.plot(wavelength, spec_tot)
plt.xlabel(r'Wavelength (\AA)')
plt.ylabel(r'Specrum')

pdf_name = model_name + '_spec.pdf'
plt.savefig(pdf_name, dpi=1200, bbox_inches='tight', pad_inches=0.02)

plt.show()
