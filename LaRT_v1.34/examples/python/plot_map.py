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

def plot_intensity_map(fname, maxI=None, maxI0=None, display=False):
   hdu         = fits.open(fname+'.fits.gz')
   wavelength  = hdu[1].data['LAMBDA']
   dwavelength = hdu[1].header['DLAMBDA']
   zmax        = hdu[1].header['ZMAX']
   rmax        = zmax
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

   distance = np.sqrt(obsx**2 + obsy**2 + obsz**2)
 
   scatt   = hdu_obs[0].data
   direc   = hdu_obs[1].data
   #-- integrate over frequency or wavelength.
   Iscatt  = np.sum(scatt,  axis=2) * bin_unit
   Idirec  = np.sum(direc,  axis=2) * bin_unit
   Itot    = Iscatt + Idirec
   try:
      direc0  = hdu_obs[2].data
      Idirec0 = np.sum(direc0, axis=2) * bin_unit
      draw_direc0 = True
   except:
      draw_direc0 = False
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

   ext1 = [-xp_max,xp_max,-yp_max,yp_max]
   #vmin1, vmax1 = find_minmax(Itot,    log_scale=True)
   vmin1, vmax1 = find_minmax(gaussian_filter(Itot,   sigma=2.0), log_scale=True)
   print('  max total  Intensity = %15.7e' % (10.0**vmax1))

   w1 = np.where(Idirec  <= 0.0)
   w2 = np.where(Iscatt  <= 0.0)
   w3 = np.where(Itot    <= 0.0)
   if len(w1[0]) > 0: Idirec[w1]  = np.nan
   if len(w2[0]) > 0: Iscatt[w2]  = np.nan
   if len(w3[0]) > 0: Itot[w3]    = np.nan

   #-- plot images
   if draw_direc0:
      ncol  = 4
      nrow  = 1
      extra = 0.04 ## extra space for colorbar
   else:
      ncol  = 3
      nrow  = 1
      extra = 0.04 ## extra space for colorbar

   xsize = 3.0*(nxim/nyim)
   ysize = 3.0
   fig   = plt.figure(figsize=(xsize*(1.0+extra)*ncol,ysize*nrow))
   gap_w = 0.35/ncol
   gap_h = 0.35/nrow
   box_w = (1.0-extra)/ncol - gap_w
   box_h = 1.0/nrow - gap_h
   px    = gap_w/2.0
   if nrow == 1:
      py = gap_h/2.0
   else:
      py = 1.0 - gap_h/2.0 - box_h
   
   #### (1) total light
   ax   = fig.add_axes([px,py,box_w,box_h])
   if (maxI != None): vmax1 = np.log10(maxI)
   image = ax.imshow(Itot, interpolation='none', extent=ext1, cmap=cmap, origin='lower', vmin=0.0, vmax=10.0**vmax1)
   ax.set_title(r'$I_{\rm tot}$ - {\small %s}' % fname.replace('atm_','').replace('_','\\textunderscore{}'))

   cax    = inset_axes(ax, width="5%", height="100%", loc="lower left",
                       bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
   cbar   = fig.colorbar(image, cax=cax)
   cbar.set_label(r'Intensity [sr$^{-1}$ s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]')

   #### (2) direct light
   if nrow == 1:
      ax   = fig.add_axes([ax.get_position().x1 + gap_w + extra/ncol, py, box_w, box_h])
   else:
      ax   = fig.add_axes([px, ax.get_position().y0 - gap_h - box_h, box_w, box_h])
   image = ax.imshow(Idirec, interpolation='none', extent=ext1, cmap=cmap, origin='lower', vmin=0.0, vmax=10.0**vmax1)
   ax.set_title(r'$I_{\rm direc}$')

   cax    = inset_axes(ax, width="5%", height="100%", loc="lower left",
                       bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
   cbar   = fig.colorbar(image, cax=cax)
   cbar.set_label(r'Intensity [sr$^{-1}$ s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]')

   #### (3) scattered light
   if nrow == 1:
      ax    = fig.add_axes([ax.get_position().x1 + gap_w + extra/ncol, py, box_w, box_h])
   else:
      ax   = fig.add_axes([px, ax.get_position().y0 - gap_h - box_h, box_w, box_h])
   image = ax.imshow(Iscatt, interpolation='none', extent=ext1, cmap=cmap, origin='lower', vmin=0.0, vmax=10.0**vmax1)
   ax.set_title(r'$I_{\rm scatt}$')

   cax    = inset_axes(ax, width="5%", height="100%", loc="lower left",
                       bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
   cbar   = fig.colorbar(image, cax=cax)
   cbar.set_label(r'Intensity [sr$^{-1}$ s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]')

   #### (4) direct0 light
   if draw_direc0:
      #vmin2, vmax2 = find_minmax(Idirec0, log_scale=True)
      vmin2, vmax2 = find_minmax(gaussian_filter(Idirec0,sigma=2.0), log_scale=True)
      print('  max direc0 Intensity = %15.7e' % (10.0**vmax2))
      w4 = np.where(Idirec0 <= 0.0)
      if len(w4[0]) > 0: Idirec0[w4] = np.nan
      if nrow == 1:
         ax    = fig.add_axes([ax.get_position().x1 + gap_w + extra/ncol, py, box_w, box_h])
      else:
         ax   = fig.add_axes([px, ax.get_position().y0 - gap_h - box_h, box_w, box_h])
      if (maxI != None): vmax2 = np.log10(maxI0)
      image = ax.imshow(Idirec0, interpolation='none', extent=ext1, cmap=cmap, origin='lower', vmin=0.0, vmax=10.0**vmax2)
      ax.set_title(r'$I_{\rm direc0}$')

      cax    = inset_axes(ax, width="5%", height="100%", loc="lower left",
                          bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
      cbar   = fig.colorbar(image, cax=cax)
      cbar.set_label(r'Intensity [sr$^{-1}$ s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]')
   
   ####
   pdf_file = fname + '_map.pdf'
   plt.subplots_adjust(wspace=0.45)
   plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
   print('  saving %s' % pdf_file)
   if display == True:
      plt.show()
   else:
      plt.close()

#-------------------------
#display = False
display = True

flist = glob.glob('*_obs.fits.gz')
flist.sort()
nf    = len(flist)

for i in np.arange(nf):
    fname = flist[i][0:-12]
    print('\n>>> plotting %s...' % fname)
    try:
       plot_intensity_map(fname, display=display)
    except:
       plot_intensity_map(fname, display=display)

#----------------------------------------------
#--- combine output pdf files.
pdf_list = ''
for i in np.arange(nf):
   fname    = flist[i][0:-12]
   pdf_list = pdf_list+' '+fname+'_map.pdf'

if nf > 1:
   all_pdf = 'all_map.pdf'
   os.system("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s" % (all_pdf, pdf_list))
   #-- remove the individual pdf files.
   for i in np.arange(nf):
      pdf_file = flist[i][0:-12]+'_map.pdf'
      os.remove(pdf_file)
