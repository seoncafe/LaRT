#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import healpy as hp 
import utils

model_name = 'CIV_test'
a = fits.open(model_name + '_obs.fits.gz')

# delta frequency
dxfreq = a[0].header['DXFREQ']
# scattered and direct light
I_scatt = a[0].data
I_direc = a[1].data

# integrate over frequency
map_scatt = np.sum(I_scatt, 1) * dxfreq
map_direc = np.sum(I_direc, 1) * dxfreq
map_tot   = map_scatt + map_direc

vmin, vmax = utils.find_histogram_limits(map_tot, log=True)
hp.mollview(map_tot, norm='log', title=model_name, min=vmin, max=vmax)
hp.graticule()

pdf_name = model_name + '_map.pdf'
plt.savefig(pdf_name, dpi=1200, bbox_inches='tight', pad_inches=0.02)

plt.show()
