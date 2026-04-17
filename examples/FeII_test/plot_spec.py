#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import copy
import matplotlib as mpl
from   astropy.io import fits
from   matplotlib.ticker import AutoMinorLocator, MultipleLocator
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
from   square_subplots import square_subplots
from   scipy.ndimage import gaussian_filter

np.seterr(divide='ignore', invalid='ignore')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':14})
plt.rcParams.update({'lines.linewidth':0.9})
fontsize = 11

cmap  = 'seismic'

pdf_file = 'pol_spec.pdf'
fig, ax = plt.subplots(1,1,figsize=(6,3.5))

dir       = './'
fname_arr = ['FeII_UV3']

arr = np.arange(len(fname_arr))
for idx in arr:
   fname  = fname_arr[idx]
   hdu1   = fits.open(dir+fname+'.fits.gz')
   hdu3   = fits.open(dir+fname+'_stokes.fits.gz')

   wavl   = hdu1[1].data['WAVELENGTH']
   vel    = hdu1[1].data['VELOCITY']
   spec   = hdu1[1].data['JOUT']

   hdr1   = hdu1[1].header
   tau0   = hdr1['TAUMAX']
   Ngas   = hdr1['NGASMAX']
   Vexp   = hdr1['VEXP']
   temp   = hdr1['TEMP']
   vtherm = 0.12843374*np.sqrt(temp)

   ax.plot(wavl, spec)
   ax.set_xlim(2339.0, 2384.0)
   ax.set_ylim(0.0, 2.1)
   ax.set_xlabel(r'Wavelength [\rm\AA]')
   ax.set_ylabel(r'Spectrum')

pdf_file = 'FeII_UV3_spec.pdf'
plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
plt.show()
