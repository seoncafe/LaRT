#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib as mpl
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import glob
from   lart_io import load_lart, find_lart_file, glob_lart
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

def plot_intensity_profile(fname, maxI=None, maxI0=None, display=True):
   main_path = find_lart_file(fname)
   if main_path is None:
      raise FileNotFoundError(f'No LaRT output found for {fname!r}')
   main = load_lart(main_path)
   spec = main.section('Spectrum')
   wavelength  = spec.col('wavelength')
   dwavelength = float(spec.attr('DWAVE', 1.0))
   xmax        = float(spec.attr('XMAX', 1.0))
   rmax        = xmax

   obs_path = find_lart_file(fname, suffix='_obs')
   if obs_path is None:
      raise FileNotFoundError(f'No peel-off (_obs) file found for {fname!r}')
   obs       = load_lart(obs_path)
   scattered = obs.section('Scattered')
   direct    = obs.section('Direct')
   direct0   = obs.section('Direct0')

   dxim = float(scattered.attr('CD1_1', 1.0))
   dyim = float(scattered.attr('CD2_2', 1.0))
   intensity_unit = int(scattered.attr('I_UNIT', 0))
   if intensity_unit == 0:
      bin_unit = float(scattered.attr('DXFREQ', 1.0))
   else:
      bin_unit = float(scattered.attr('DWAVE', 1.0))
   distance_unit = float(scattered.attr('DIST_CM', 1.0))
   distance_in   = float(scattered.attr('DISTANCE', 0.0))
   obsx = float(scattered.attr('OBSX', 0.0))
   obsy = float(scattered.attr('OBSY', 0.0))
   obsz = float(scattered.attr('OBSZ', 0.0))
   distance_star_to_planet = float(scattered.attr('DIST_S2P', 0.0))
   flux_factor    = float(scattered.attr('FLUXFAC', 1.0))
   limb_darkening = int(scattered.attr('LIMBDARK', 0))

   # Now (LaRT_v1.33a), distance is defined to be the distance from the rotation center to the observer.
   distance = np.sqrt(obsx**2 + obsy**2 + (obsz + distance_star_to_planet)**2)

   scatt  = scattered.data * flux_factor
   direc  = direct.data
   direc0 = direct0.data if direct0 is not None else np.zeros_like(scatt)
   #-- integrate over frequency or wavelength.
   Iscatt  = np.sum(scatt,  axis=2) * bin_unit
   Idirec  = np.sum(direc,  axis=2) * bin_unit
   Idirec0 = np.sum(direc0, axis=2) * bin_unit
   Itot    = Iscatt + Idirec

   #-- Note that rmax == rp_max if only par%nxim and par%nyim are given in the input file.
   #             (in this case, par%dxim and par%dyim are automatically determined.)
   #        but  rmax /= rp_max if par%dxim and par%dyim are also specified in the input file.
   #-- In general, rmax /= rp_max
   #-- In gneral, xp_max /= yp_max if nxim*dxim /= nyim*dyim
   #-- xp_max, yp_max = maximum sizes in the detector plane.
   nyim, nxim = Itot.shape
   xp_max     = distance * np.sin(dxim*(nxim/2.0) * (np.pi/180.0))
   yp_max     = distance * np.sin(dyim*(nyim/2.0) * (np.pi/180.0))
   #-- spherical model.
   rp_max     = np.amin([xp_max, yp_max])

   r0, Itot_r    = radial_profile(Itot)
   r0, Iscatt_r  = radial_profile(Iscatt)
   r0, Idirec_r  = radial_profile(Idirec)
   r0, Idirec0_r = radial_profile(Idirec0)
   r = r0 * rp_max

   #-- plot images
   ncol  = 2
   nrow  = 1
   fig   = plt.figure(figsize=(4.2*ncol,4.2*nrow))
   gap_w = 0.25/ncol
   gap_h = 0.25
   box_w = 1.0/ncol - gap_w
   box_h = 1.0/nrow - gap_h
   px    = gap_w/2.0
   py    = gap_h/2.0
   
   #### (1) scattered, direc, total light
   ax       = fig.add_axes([px,py,box_w,box_h])

   ex       = np.floor(np.log10(np.nanmax(Itot_r)))
   scal     = 10.0**(-ex)
   Itot_r   = Itot_r   * scal
   Idirec_r = Idirec_r * scal
   Iscatt_r = Iscatt_r * scal

   ax.plot(r, Itot_r,   color='k', label=r'$I_{\rm tot}$')
   ax.plot(r, Iscatt_r, color='r', label=r'$I_{\rm scatt}$ (FLUXFAC = %0.3e)' % flux_factor)
   ax.plot(r, Idirec_r, color='b', label=r'$I_{\rm direc}$')
   ax.set_title(r'%s' % fname.replace("_",'\\textunderscore{}'))
   ax.set_xlabel(r'$r_{\rm proj.}/r_{\rm planet}$')
   ax.set_ylabel(r'Intensity [10$^{%d}$ sr$^{-1}$ s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]' % ex)
   ax.set_xlim(0.0, rp_max)
   ax.set_ylim(0.0, np.amax(Itot_r))
   ax.legend(loc='upper left', frameon=False)

   #### (2) direct0 light
   ax        = fig.add_axes([ax.get_position().x1 + gap_w, py, box_w, box_h])

   # rotation about (0,0,0) gives zero intensity of Idirec0_r because the star is centered in the image plane.
   if np.nansum(Idirec0_r) > 0.0:
      ex     = np.floor(np.log10(np.nanmax(Idirec0_r)))
      scal   = 10.0**(-ex)
   Idirec0_r = Idirec0_r * scal

   ax.plot(r, Idirec0_r, color='k', label=r'$I_{\rm direc0}$')
   ax.set_xlabel(r'$r_{\rm proj.}/r_{\rm planet}$')
   ax.set_ylabel(r'Intensity [10$^{%d}$ sr$^{-1}$ s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]' % ex)
   ax.set_xlim(0.0, rp_max)
   ymax = np.nanmax(Idirec0_r)*1.1
   ax.set_ylim(0.0, ymax)
   ax.legend(loc='upper left', frameon=False)

   #---- Theoretical Profile.
   try:
      rstar = hdr['RSTAR']
   except:
      w     = np.where(Idirec0_r > 0.0)
      rstar = np.nanmax(r[w])
   sina  = r / distance
   cosa  = np.sqrt(1.0 - sina**2)
   cost  = np.sqrt((rstar/distance)**2 - sina**2) / (rstar/distance)
   if limb_darkening == 0:
      f_limb = cosa / cost
   elif limb_darkening == 1:
      f_limb = cosa
      w      = np.where(r > rstar)
      f_limb[w] = 0.0
   elif limb_darkening == 2:
      f_limb = cosa*(1.5*cost + 1.0)
   else:
      coeff  = [0.55, 0.12, 0.33]
      f_limb = cosa*(coeff[0] + coeff[1]*cost + coeff[2]*cost**2)
   w = np.where(np.isfinite(f_limb) == False)
   f_limb[w] = 0.0

   w      = np.where((r > 0.05) & (r < rstar*0.95))
   f_limb = f_limb/np.nansum(f_limb[w]) * np.nansum(Idirec0_r[w])
   ax.plot(r, f_limb, color='b', linestyle='', marker='o', ms=2)
   ax.set_title(r'\small blue: theoretical profile', color='b')
   
   #------------------------------------------------
   pdf_file = fname + '_intensity_profile.pdf'
   plt.subplots_adjust(wspace=0.45)
   plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
   print(' saved %s' % pdf_file)
   if display == True:
      plt.show()
   else:
      plt.close()

#-------------------------
display = False
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
   fname = _strip_ext(flist[i])[:-4]  # strip '_obs'
   print('plotting %s...' % fname)
   plot_intensity_profile(fname, display=display)

#----------------------------------------------
#--- combine output pdf files.
pdf_list = ''
for i in np.arange(nf):
   fname    = _strip_ext(flist[i])[:-4]
   pdf_list = pdf_list+' '+fname+'_intensity_profile.pdf'

if nf > 1:
   all_pdf = 'all_intensity_profile.pdf'
   os.system("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s" % (all_pdf, pdf_list))
   #-- remove the individual pdf files.
   for i in np.arange(nf):
      fname    = _strip_ext(flist[i])[:-4]
      pdf_file = fname+'_intensity_profile.pdf'
      os.remove(pdf_file)
