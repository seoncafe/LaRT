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

def plot_intensity_profile(fname, maxI=None, maxI0=None, display=True):
   hdu         = fits.open(fname+'.fits.gz')
   wavelength  = hdu[1].data['LAMBDA']
   dwavelength = hdu[1].header['DLAMBDA']
   xmax        = hdu[1].header['XMAX']
   rmax        = xmax
   hdu.close()
   
   hdu_obs = fits.open(fname+'_obs.fits.gz')
   hdr     = hdu_obs[0].header
   dxim    = hdr['CD1_1']
   dyim    = hdr['CD2_2']
   try:
      intensity_unit = hdr['I_UNIT']
   except:
      intensity_unit = 0
   if intensity_unit == 0: bin_unit = hdr['DXFREQ']
   if intensity_unit == 1: bin_unit = hdr['DLAMBDA']
   distance_unit = hdr['DIST_CM']
   distance_in   = hdr['DISTANCE']
   obsx          = hdr['OBSX']
   obsy          = hdr['OBSY']
   obsz          = hdr['OBSZ']
   try:
      distance_star_to_planet = hdr['DIST_S2P']
   except:
      distance_star_to_planet = 0.0
   try:
      flux_factor = hdr['FLUXFAC']
   except:
      flux_factor = 1.0
   try:
      limb_darkening = hdr['LIMBDARK']
   except:
      limb_darkening = 0

   # Now (LaRT_v1.33a), distance is defined to be the distance from the rotation center to the observer.
   distance = np.sqrt(obsx**2 + obsy**2 + (obsz + distance_star_to_planet)**2)

   scatt   = hdu_obs[0].data * flux_factor
   direc   = hdu_obs[1].data
   direc0  = hdu_obs[2].data
   #-- integrate over frequency or wavelength.
   Iscatt  = np.sum(scatt,  axis=2) * bin_unit
   Idirec  = np.sum(direc,  axis=2) * bin_unit
   Idirec0 = np.sum(direc0, axis=2) * bin_unit
   Itot    = Iscatt + Idirec
   hdu_obs.close()

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
#flist = list(set(glob.glob('*_obs.fits.gz')) - set(glob.glob('*diffuse*_obs.fits.gz')))
flist = glob.glob('*_obs.fits.gz')
flist.sort()
nf    = len(flist)

for i in np.arange(nf):
   fname = flist[i][0:-12]
   print('plotting %s...' % fname)
   plot_intensity_profile(fname, display=display)

#----------------------------------------------
#--- combine output pdf files.
pdf_list = ''
for i in np.arange(nf):
   fname    = flist[i][0:-12]
   pdf_list = pdf_list+' '+fname+'_intensity_profile.pdf'

if nf > 1:
   all_pdf = 'all_intensity_profile.pdf'
   os.system("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s" % (all_pdf, pdf_list))
   #-- remove the individual pdf files.
   for i in np.arange(nf):
      pdf_file = flist[i][0:-12]+'_intensity_profile.pdf'
      os.remove(pdf_file)
