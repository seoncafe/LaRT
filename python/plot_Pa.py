#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import copy
import matplotlib as mpl
from   matplotlib.ticker import AutoMinorLocator, MultipleLocator
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
from   scipy.ndimage  import gaussian_filter
from   radial_profile import radial_profile
from   lart_io        import load_lart, find_lart_file, glob_lart
import glob
import os

np.seterr(divide='ignore', invalid='ignore')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':12})
plt.rcParams.update({'lines.linewidth':0.9})
plt.rcParams.update({'legend.fontsize':9})
fontsize = 12

#cmap = mpl.cm.get_cmap("jet").copy()
#cmap = mpl.cm.get_cmap("nipy_spectral").reversed().copy()
cmap = mpl.cm.get_cmap("nipy_spectral").copy()
#cmap = mpl.cm.get_cmap("viridis").copy()
cmap.set_bad('white')

def find_minmax(array, nbins=500, frac=0.05, log_scale=True):
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

def plot_Pa(model_type, gauss_sigma=1.0, frac=0.05, display=True, Lum=5e39):
   print('--- Plotting %s ' % model_type)
   dir = './'
   fname = model_type

   #--- Read P-alpha
   main_path = find_lart_file(dir+fname)
   if main_path is None:
      raise FileNotFoundError(f'No LaRT output found for {dir+fname!r}')
   main = load_lart(main_path)
   spec = main.section('Spectrum')
   xfreq = np.asarray(spec.col('Xfreq')).squeeze()
   #--- Pa is the next section after Spectrum.  In FITS this was HDU 2 (an
   #    image extension); in HDF5 it shows up as the second group at root
   #    (typically /Pa_2D or /Pa_3D depending on geometry_JPa).
   pa_sec = None
   for s in main.sections:
      if s.name.lower().startswith('pa'):
         pa_sec = s
         break
   if pa_sec is None or pa_sec.data is None:
      raise ValueError(f'No P-alpha image found in {main_path!r}')
   Pa = pa_sec.data
   zmax = float(spec.attr('ZMAX', 1.0))
   rmax = zmax

   obs_path = find_lart_file(dir+fname, suffix='_obs')
   if obs_path is not None:
      obs = load_lart(obs_path)
      distance_unit = float(obs.section('Scattered').attr('DIST_CM', 1.0))
   else:
      distance_unit = 1.0

   nz, nr = Pa.shape
   if (nr//2)*2 == nr:
      Pa1 = np.flip(Pa,       axis=1)
   else:
      Pa1 = np.flip(Pa[:,1:], axis=1)

   Pa     = np.concatenate((Pa1, Pa), axis=1)
   nz, nr = Pa.shape

   #---
   #--- Update to use a correct luminosity.
   #---
   Pa  = Pa * Lum

   w = np.where(Pa == 0.0)
   if w[0].size > 0:
      Pa[w]  = np.nan

   #---- plot figures
   fig = plt.figure(figsize=(4.2,4.2))
   gap_w = 0.35
   gap_h = 0.35
   box_w = 1.0 - gap_w
   box_h = box_w
   px    = gap_w/2.0
   py    = gap_h/2.0

   #--- Pa image
   ax   = fig.add_axes([px,py,box_w,box_h])

   ext1 = [-rmax,rmax,-rmax,rmax]
   #vmin1, vmax1 = find_minmax(Pa, frac=frac, log_scale=True)
   vmin1, vmax1 = find_minmax(gaussian_filter(Pa,sigma=2.0), frac=frac, log_scale=True)

   image = ax.imshow(Pa, extent=ext1, cmap=cmap, origin='lower', vmin=0.0, vmax=10.0**vmax1)
   ex    = np.floor(np.log10(Lum))
   scal  = 10.0**(-ex)
   ax.set_title(r'($L = %3.2f\times10^{%d}$) - %s' % (Lum*scal, ex, model_type.replace("_","\\textunderscore{}")))

   ax.set_xlabel(r'$r_{\rm cyl}/r_{\rm planet}$')
   ax.set_ylabel(r'$z/r_{\rm planet}$')

   tick_interval = 5.0
   if (rmax <=  3.0): tick_interval = 1.0
   if (rmax <=  5.0): tick_interval = 2.0
   if (rmax >= 20.0): tick_interval = 10.0
   ax.xaxis.set_major_locator(MultipleLocator(tick_interval))
   ax.yaxis.set_major_locator(MultipleLocator(tick_interval))

   cax  = inset_axes(ax, width="5%", height="100%", loc="lower left",
                     bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
   cbar = fig.colorbar(image, cax=cax)
   cbar.set_label(r'P$\alpha$ [s$^{-1}$]')

   #--- Save
   pdf_file = model_type+'_Pa.pdf'
   # the amount of width reserved for space between subplots,
   # expressed as a fraction of the average axis width (default is 0.2).
   plt.subplots_adjust(wspace=0.45)
   plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
   print("    saved %s" % pdf_file)
   if display == True:
      plt.show()
   else:
      plt.close()

#------------------------
#display = False
display = True

def _strip_ext(path):
   lower = path.lower()
   for ext in ('.fits.gz', '.hdf5', '.fits', '.h5'):
      if lower.endswith(ext):
         return path[: -len(ext)]
   return path

flist = glob_lart('*', suffix='_obs')
flist.sort()
nf    = len(flist)

for i in np.arange(nf):
    fname = _strip_ext(flist[i])[:-4]   # strip '_obs'
    plot_Pa(fname, display=display)

#----------------------------------------------
#--- combine output pdf files.
pdf_list = ''
for i in np.arange(nf):
   fname    = _strip_ext(flist[i])[:-4]
   pdf_list = pdf_list+' '+fname+'_Pa.pdf'

if nf > 1:
   all_pdf = 'all_Pa.pdf'
   os.system("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s" % (all_pdf, pdf_list))
   #-- remove the individual pdf files.
   for i in np.arange(nf):
      fname    = _strip_ext(flist[i])[:-4]
      pdf_file = fname+'_Pa.pdf'
      os.remove(pdf_file)
