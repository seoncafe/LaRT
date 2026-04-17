#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

def plot_spec(ax, fname, ext=None,ymax=None):
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

   beta1 = c1[0].header['beta']
   beta2 = c2[0].header['beta']
   beta3 = c3[0].header['beta']
   beta4 = c4[0].header['beta']

   ax.plot(vel, spec1, label=r'$\beta=%d^\circ$' % beta1)
   ax.plot(vel, spec2, label=r'$\beta=%d^\circ$' % beta2)
   ax.plot(vel, spec3, label=r'$\beta=%d^\circ$' % beta3)
   ax.plot(vel, spec4, label=r'$\beta=%d^\circ$' % beta4)

   ax.set_xlim(-500.0, 500.0)
   ax.set_xlabel(r'Velocity [km s$^{-1}$]')
   ax.set_ylabel(r'Spectrum')
   ii = fname.index('Vrot')+4
   Vrot = np.float64(fname[ii:ii+3])
   if ext == None:
      ax.set_title(r'V$_{\rm rot}$ = %d km s$^{-1}$' % Vrot)
   else:
      ax.set_title(r'V$_{\rm rot}$ = %d km s$^{-1}$ (%s)' % (Vrot, ext))

   ax.legend()

fig, axs = plt.subplots(ncols=2,nrows=2, figsize=(10,10))
axs = axs.flat

ext   = r'n$_x$=201, N$_{\rm HI}$ = 10$^{18}$'
fname = 'nx201/rin0.1_Vrot100_NHI18'
ax    = axs[0]
plot_spec(ax,fname,ext=ext)

fname = 'nx201/rin0.1_Vrot300_NHI18'
ax    = axs[1]
plot_spec(ax,fname,ext=ext)

ext   = r'n$_x$=801, N$_{\rm HI}$ = 10$^{18}$'
fname = 'rin0.1_Vrot100_NHI18'
ax    = axs[2]
plot_spec(ax,fname,ext=ext)

fname = 'rin0.1_Vrot300_NHI18'
ax    = axs[3]
plot_spec(ax,fname,ext=ext)

pdf_file = 'spec_compare.pdf'
plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
plt.show()
