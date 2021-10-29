#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib as mpl
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from   astropy.io import fits
import glob
from   radial_profile import radial_profile
import os

np.seterr(divide='ignore', invalid='ignore')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':12})
plt.rcParams.update({'lines.linewidth':0.9})
plt.rcParams.update({'legend.fontsize':9})
fontsize = 12

cmap = mpl.cm.get_cmap("nipy_spectral").copy()
cmap.set_bad('white')

def find_minmax(array, nbins=500, frac=0.10, log_scale=True):
   arr             = array.reshape(array.size)
   if log_scale == True:
      w            = np.where(arr > 0.0)
      arr          = np.log10(arr[w])
   else:
      w            = np.where(arr > 0.0)
      arr          = arr[w]
   hist, bin_edges = np.histogram(arr, bins=nbins, density=True)
   wh              = (np.where(hist == np.amax(hist)))[0][0]
   xarr            = (bin_edges[0:-1] + bin_edges[1:])/2.0
   if wh == 0:
      xmin         = np.amax(hist)*frac
   else:
      xmin         = np.interp(np.amax(hist)*frac, hist[:wh], xarr[:wh])
   xmax            = np.interp(np.amax(hist)*frac, np.flip(hist[wh:]), np.flip(xarr[wh:]))
   return xmin, xmax

def plot_spectrum(fname, maxI=None, maxI0=None, display=True):
   hdu     = fits.open(fname+'.fits.gz')
   wavelength  = hdu[1].data['LAMBDA']
   dwavelength = hdu[1].header['DLAMBDA']
   zmax        = hdu[1].header['ZMAX']
   rmax        = zmax
   hdu.close()
   
   hdu_obs = fits.open(fname+'_obs.fits.gz')
   hdr     = hdu_obs[0].header

   dxim    = hdr['CD1_1'] * np.pi/180.0
   dyim    = hdr['CD2_2'] * np.pi/180.0
   try:
      flux_factor = hdr['FLUXFAC']
   except:
      flux_factor = 1.0
   try:
      intensity_unit = hdr['I_UNIT']
   except:
      intensity_unit = 0
   if intensity_unit == 0: bin_unit = hdr['DXFREQ']
   if intensity_unit == 1: bin_unit = hdr['DLAMBDA']

   scatt   = hdu_obs[0].data * flux_factor
   direc   = hdu_obs[1].data
   #-- integrate over frequency or wavelength.
   Spec_scatt  = np.sum(scatt,  axis=(0,1)) * (dxim * dyim)
   Spec_direc  = np.sum(direc,  axis=(0,1)) * (dxim * dyim)
   Spec_tot    = Spec_scatt + Spec_direc
   try:
      direc0  = hdu_obs[2].data
      Spec_direc0 = np.sum(direc0, axis=(0,1)) * (dxim * dyim)
      draw_direc0 = True
   except:
      draw_direc0 = False
   hdu_obs.close()

   #-- plot images
   if draw_direc0 == True:
      ncol = 2
   else:
      ncol = 1
   nrow = 1
   fig = plt.figure(figsize=(4.2*ncol,4.2*nrow))
   gap_w = 0.25/ncol
   gap_h = 0.25
   box_w = 1.0/ncol - gap_w
   box_h = 1.0/nrow - gap_h
   px    = gap_w/2.0
   py    = gap_h/2.0
   
   wav1 = np.amin(wavelength)
   wav2 = np.amax(wavelength)
   #### (1) scattered light
   ax   = fig.add_axes([px,py,box_w,box_h])

   ex       = np.floor(np.log10(np.nanmax(Spec_tot)))
   scal     = 10.0**(-ex)
   Spec_tot   = Spec_tot   * scal
   Spec_direc = Spec_direc * scal
   Spec_scatt = Spec_scatt * scal

   ax.plot(wavelength, Spec_tot,   color='k', label=r'Spec$_{\rm tot}$')
   ax.plot(wavelength, Spec_direc, color='r', label=r'Spec$_{\rm direc}$')
   ax.plot(wavelength, Spec_scatt, color='b', label=r'Spec$_{\rm scatt}$ (FLUXFAC = %0.3e)' % flux_factor)
   ax.set_title(r'%s' % fname.replace("_",'\\textunderscore{}'))
   ax.set_xlabel(r'Wavelength (\AA)')
   ax.set_ylabel(r'Intensity [10$^{%d}$ s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]' % ex)
   ax.set_xlim(wav1, wav2)
   ax.set_ylim(0.0, np.amax(Spec_tot) * 1.05)
   ax.legend(loc='upper left', frameon=False)

   #### (2) direct0 light
   if draw_direc0:
      ax    = fig.add_axes([ax.get_position().x1 + gap_w, py, box_w, box_h])

      ex        = np.floor(np.log10(np.nanmax(Spec_direc0)))
      scal      = 10.0**(-ex)
      Spec_direc0 = Spec_direc0 * scal

      ax.plot(wavelength, Spec_direc0, color='k', label=r'Spec$_{\rm direc0}$')
      ax.set_xlabel(r'Wavelength (\AA)')
      ax.set_ylabel(r'Intensity [10$^{%d}$ s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]' % ex)
      ax.set_xlim(wav1, wav2)
      ax.set_ylim(0.0, np.amax(Spec_direc0) * 1.05)
      ax.legend(loc='upper left', frameon=False)
   
   pdf_file = fname + '_spectrum.pdf'
   plt.subplots_adjust(wspace=0.45)
   plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
   print(' saved %s' % pdf_file)
   if display == True:
      plt.show()
   else:
      plt.close()

#-------------------------
#display = False
display = True
#flist = list(set(glob.glob('*_obs.fits.gz')) - set(glob.glob('*diffuse*_obs.fits.gz')))
flist = glob.glob('*_obs.fits.gz')
flist.sort()
nf    = len(flist)

for i in np.arange(nf):
    fname = flist[i][0:-12]
    print('plotting %s...' % fname)
    plot_spectrum(fname, display=display)

#----------------------------------------------
#--- combine output pdf files.
pdf_list = ''
for i in np.arange(nf):
   fname    = flist[i][0:-12]
   pdf_list = pdf_list+' '+fname+'_spectrum.pdf'

if nf > 1:
   all_pdf = 'all_spectrum.pdf'
   os.system("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s" % (all_pdf, pdf_list))
   #-- remove the individual pdf files.
   for i in np.arange(nf):
      pdf_file = flist[i][0:-12]+'_spectrum.pdf'
      os.remove(pdf_file)
