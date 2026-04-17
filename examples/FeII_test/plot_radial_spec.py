#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import copy
import matplotlib as mpl
from   astropy.io import fits
from   matplotlib.ticker import AutoMinorLocator, MultipleLocator
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
from   scipy.ndimage  import gaussian_filter
from   radial_profile import *

def plot_data(fname, rbins=None, mode=0):
   ncol = 1
   nrow = 1
   fig, ax = plt.subplots(nrow,ncol,figsize=(8*ncol,6*nrow))

   dir = './'

   hdu   = fits.open(dir+fname+'.fits.gz')
   wavl  = hdu[1].data['WAVELENGTH'].squeeze()
   vel   = hdu[1].data['VELOCITY'].squeeze()
   dwavl = np.abs(wavl[1] - wavl[0])
   dvel  = np.abs(vel[1]  - vel[0])
   hdu.close()

   hdu_obs     = fits.open(dir+fname+'_obs.fits.gz')
   peel_scatt  = hdu_obs[0].data
   peel_direc  = hdu_obs[1].data
   peel_tot    = peel_scatt + peel_direc
   #try:
   #   peel_direc0 = hdu_obs[2].data
   #except:
   #   pass
   hdu_obs.close()

   ny, nx, nwavelength = peel_tot.shape

   #rbins = np.array([0.0, 0.05, 0.15, 0.25, 0.35, 0.55, 0.6, 0.7, 1.0])
   #rbins = np.array([0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 1.0])
   if rbins == None:
      #rbins = np.array([0.0, 0.1, 0.4, 0.7, 1.0])
      #rbins = np.array([0.0, 0.4, 1.0])
      rbins = np.array([0.0, 1.0])

   r_spec,  spec       = radial_spectrum(peel_tot,   rbins = rbins)
   r_spec1, spec_scatt = radial_spectrum(peel_scatt, rbins = rbins)

   vel1  = np.amin(vel)  - dvel/2.0
   vel2  = np.amax(vel)  + dvel/2.0
   wavl1 = np.amin(wavl) - dwavl/2.0
   wavl2 = np.amax(wavl) + dwavl/2.0

   scale      = np.median(spec[0,:])
   spec       = spec/scale
   spec_scatt = spec_scatt/scale

   for i in np.arange(r_spec.size):
      if mode == 0:
         ax.plot(wavl, spec[i,:])
         ax.plot(wavl, spec_scatt[i,:])
      else:
         ax.plot(vel,  spec[i,:])
         ax.plot(vel,  spec_scatt[i,:])

   if mode == 0:
      ax.set_xlim(wavl1, wavl2)
      ax.set_xlabel(r'Wavelength [$\rm\AA$]')
   else:
      ax.set_xlim(vel1, vel2)
      ax.set_xlabel(r'Velocity [km $s^{-1}$]')

   ax.set_ylim(0.0, np.nanmax(spec)*1.02)
   #ax.set_ylim(0.0, 1.1)
   ax.set_ylabel(r'Spectrum')

   plt.subplots_adjust(wspace=0.45)
   #plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
   plt.show()

def main():
   import sys
   nargs = len(sys.argv)
   fname = sys.argv[1]
   plot_data(fname)

if __name__ == "__main__":
   main()
