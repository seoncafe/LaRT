#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib as mpl
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import glob
from   radial_profile import radial_profile
from   lart_io       import load_lart, find_lart_file, glob_lart
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
   # Load the spectrum (Spectrum BinTable HDU in FITS, /Spectrum group in HDF5).
   main_path = find_lart_file(fname)
   if main_path is None:
      raise FileNotFoundError(f'No LaRT output found for {fname!r}')
   main = load_lart(main_path)
   spec = main.section('Spectrum')
   wavelength = spec.col('wavelength')

   obs_path = find_lart_file(fname, suffix='_obs')
   if obs_path is not None:
      obs = load_lart(obs_path)
      scattered = obs.section('Scattered')
      direct    = obs.section('Direct')
      dxim = float(scattered.attr('CD1_1', 1.0)) * np.pi/180.0
      dyim = float(scattered.attr('CD2_2', 1.0)) * np.pi/180.0
      flux_factor = float(scattered.attr('FLUXFAC', 1.0))

      scatt = scattered.data * flux_factor
      direc = direct.data
      Spec_scatt  = np.sum(scatt, axis=(0,1)) * (dxim * dyim)
      Spec_direc  = np.sum(direc, axis=(0,1)) * (dxim * dyim)
      Spec_tot    = Spec_scatt + Spec_direc
      direct0 = obs.section('Direct0')
      if direct0 is not None:
         Spec_direc0 = np.sum(direct0.data, axis=(0,1)) * (dxim * dyim)
         draw_direc0 = True
      else:
         draw_direc0 = False
      peeloff = True
   else:
      Spec_tot = spec.col('Jout')
      jin = spec.col('Jin')
      if jin is not None:
         Spec_direc0 = jin
         draw_direc0 = True
      else:
         draw_direc0 = False
      peeloff = False

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
   Spec_tot = Spec_tot   * scal
   ax.plot(wavelength, Spec_tot,   color='k', label=r'Spec$_{\rm tot}$')

   if peeloff == True:
      Spec_direc = Spec_direc * scal
      Spec_scatt = Spec_scatt * scal
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

def _strip_ext(path):
   lower = path.lower()
   for ext in ('.fits.gz', '.hdf5', '.fits', '.h5'):
      if lower.endswith(ext):
         return path[: -len(ext)]
   return path

flist = glob_lart('*', suffix='_obs')
if len(flist) == 0:
   flist   = glob_lart('*')
   peeloff = False
else:
   peeloff = True
flist.sort()
nf    = len(flist)

for i in np.arange(nf):
   stem = _strip_ext(flist[i])
   fname = stem[:-4] if peeloff else stem   # strip '_obs' suffix
   print('plotting %s...' % fname)
   plot_spectrum(fname, display=display)

#----------------------------------------------
#--- combine output pdf files.
pdf_list = ''
for i in np.arange(nf):
   stem = _strip_ext(flist[i])
   fname = stem[:-4] if peeloff else stem
   pdf_list = pdf_list+' '+fname+'_spectrum.pdf'

if nf > 1:
   all_pdf = 'all_spectrum.pdf'
   os.system("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s" % (all_pdf, pdf_list))
   #-- remove the individual pdf files.
   for i in np.arange(nf):
      stem = _strip_ext(flist[i])
      fname = stem[:-4] if peeloff else stem
      pdf_file = fname+'_spectrum.pdf'
      os.remove(pdf_file)
