#!/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib as mpl
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from   scipy.ndimage import gaussian_filter
from   astropy.io import fits
import glob
import os

np.seterr(divide='ignore', invalid='ignore')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':10})
plt.rcParams.update({'lines.linewidth':0.9})
plt.rcParams.update({'legend.fontsize':9})

# to make set_bad work, interpolation='none' must be used in imshow.
cmap = mpl.cm.get_cmap("nipy_spectral").copy()
cmap.set_bad(color='white')

def find_minmax(array, nbins=500, frac=0.10, log_scale=True):
   arr             = array.reshape(array.size)
   if log_scale == True:
      w            = np.where(arr > 0.0)
      arr          = np.log10(arr[w])
   else:
      w            = np.where(arr > 0.0)
      arr          = arr[w]
   hist, bin_edges = np.histogram(arr, bins=nbins, density=True)
   wh              = (np.where(hist == np.nanmax(hist)))[0][0]
   xarr            = (bin_edges[0:-1] + bin_edges[1:])/2.0
   if wh == 0:
      xmin         = np.nanmax(hist)*frac
   else:
      xmin         = np.interp(np.nanmax(hist)*frac, hist[:wh], xarr[:wh])
   xmax            = np.interp(np.nanmax(hist)*frac, np.flip(hist[wh:]), np.flip(xarr[wh:]))
   return xmin, xmax

#-- plot images
ncol  = 4
nrow  = 3
xsize = 3.0*ncol
ysize = 3.0*nrow
#fig, ax = plt.subplots(nrow,ncol, figsize=(xsize,ysize), sharex=True,sharey=True)
fig, ax = plt.subplots(nrow,ncol, figsize=(xsize,ysize))

fname0      = 'stars2'
hdu         = fits.open(fname0+'.fits.gz')
wavelength  = hdu[1].data['LAMBDA']
dwavelength = hdu[1].header['DLAMBDA']
zmax        = hdu[1].header['ZMAX']
rmax        = zmax
hdu.close()

flist = glob.glob(fname0+'_obs_2D??.fits.gz')
flist.sort()
nf    = len(flist)

vmin  =  np.inf
vmax  = -np.inf

#-- find min, max intensities.
for kk in np.arange(nf):
   fname   = flist[kk]
   hdu_obs = fits.open(fname)
   hdr     = hdu_obs[0].header

   Iscatt = hdu_obs[0].data
   Idirec = hdu_obs[1].data
   Itot   = Iscatt + Idirec
   hdu_obs.close()

   vmin1, vmax1 = find_minmax(gaussian_filter(Itot, sigma=2.0), frac=0.05, log_scale=True)
   if vmin1 < vmin: vmin = vmin1
   if vmax1 > vmax: vmax = vmax1

#for kk in np.arange(nf):
for kk in np.arange(ncol*nrow):
   j = kk // ncol 
   i = kk - j*ncol

   if kk < nf:
      fname   = flist[kk]
      hdu_obs = fits.open(fname)
      hdr     = hdu_obs[0].header
      dxim    = hdr['CD1_1']
      dyim    = hdr['CD2_2']
      distance_unit = hdr['DIST_CM']
      distance_in   = hdr['DISTANCE']
      obsx          = hdr['OBSX']
      obsy          = hdr['OBSY']
      obsz          = hdr['OBSZ']
      alpha         = hdr['ALPHA']
      beta          = hdr['BETA']
      gamma         = hdr['GAMMA']
      distance      = np.sqrt(obsx**2 + obsy**2 + obsz**2)
 
      Iscatt = hdu_obs[0].data
      Idirec = hdu_obs[1].data
      Itot   = Iscatt + Idirec
      hdu_obs.close()

      nyim, nxim = Itot.shape
      xp_max     = distance * np.sin(dxim*(nxim/2.0) * (np.pi/180.0))
      yp_max     = distance * np.sin(dyim*(nyim/2.0) * (np.pi/180.0))
      ext1       = [-xp_max,xp_max,-yp_max,yp_max]
      #vmin1, vmax1 = find_minmax(gaussian_filter(Itot, sigma=2.0), frac=0.05, log_scale=True)

      w1 = np.where(Idirec <= 0.0)
      w2 = np.where(Iscatt <= 0.0)
      w3 = np.where(Itot   <= 0.0)
      if len(w1[0]) > 0: Idirec[w1] = np.nan
      if len(w2[0]) > 0: Iscatt[w2] = np.nan
      if len(w3[0]) > 0: Itot[w3]   = np.nan
   
      #--
      print(kk,j,i)
      image = ax[j,i].imshow(Itot, interpolation='none', extent=ext1, cmap=cmap, origin='lower', vmin=0.0, vmax=10.0**vmax)
      ax[j,i].set_title(r'$\alpha = %3d^{\circ}$' % alpha)

      if (i == ncol-1) | (kk == nf-1):
         cax    = inset_axes(ax[j,i], width="5%", height="100%", loc="lower left",
                             bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax[j,i].transAxes, borderpad=0,)
         cbar   = fig.colorbar(image, cax=cax)
         cbar.set_label(r'Intensity [sr$^{-1}$ s$^{-1}$ cm$^{-2}$]')
   else:
      ax[j,i].axis('off')

   ####
   #pdf_file = fname + '_map.pdf'
   #plt.subplots_adjust(wspace=0.45)
   #plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
   #print('  saving %s' % pdf_file)

plt.subplots_adjust(wspace=0.45)
plt.show()
