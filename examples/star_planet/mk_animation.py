#!/usr/bin/env python
import sys, getopt
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from   astropy.io import fits
from   matplotlib.ticker import AutoMinorLocator, MultipleLocator
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
from   scipy.ndimage import gaussian_filter
#from   astropy.convolution import convolve
#from   astropy.convolution import Gaussian2DKernel
#from   square_subplots import square_subplots

np.seterr(divide='ignore', invalid='ignore')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':10})
plt.rcParams.update({'lines.linewidth':0.9})
plt.rcParams.update({'legend.fontsize':9})

# to make set_bad work, interpolation='none' must be used in imshow.
#cmap = mpl.cm.get_cmap("jet").copy()
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

def find_minmax_from_flist(flist):
   #-- find min, max intensities.
   vmin, vmax = np.inf, -np.inf
   nf = len(flist)
   for kk in np.arange(nf):
      fname   = flist[kk]
      hdu_obs = fits.open(fname)
      hdr     = hdu_obs[0].header
      try:
         flux_factor = hdr['FLUXFAC']
      except:
         flux_factor = 1.0
      Iscatt = hdu_obs[0].data * flux_factor
      Idirec = hdu_obs[1].data
      Itot   = Iscatt + Idirec
      hdu_obs.close()

      vmin1, vmax1 = find_minmax(gaussian_filter(Itot, sigma=2.0), frac=0.05, log_scale=True)
      if vmin1 < vmin: vmin = vmin1
      if vmax1 > vmax: vmax = vmax1
   return vmin, vmax

def make_animation(argv):
   fname0 = argv
   np.seterr(divide='ignore', invalid='ignore')

   hdu   = fits.open(fname0+'.fits.gz')
   flist = glob.glob(fname0+'_obs_2D??.fits.gz')
   flist.sort()
   nf    = len(flist)

   vmin, vmax = find_minmax_from_flist(flist)
   alpha = np.zeros(nf)
   beta  = np.zeros(nf)
   gamma = np.zeros(nf)

   #======================
   for i in np.arange(nf):
      fname    = flist[i]
      hdu_obs  = fits.open(fname)
      hdr      = hdu_obs[0].header
      alpha[i] = hdr['ALPHA']
      beta[i]  = hdr['BETA']
      gamma[i] = hdr['GAMMA']
      try:
         flux_factor = hdr['FLUXFAC']
      except:
         flux_factor = 1.0
      Iscatt = hdu_obs[0].data * flux_factor
      Idirec = hdu_obs[1].data
      Itot   = Iscatt + Idirec
      ww     = np.where(Itot <= 0.0)
      if len(ww[0]) > 0: Itot[ww] = np.nan
      if i == 0:
         dxim          = hdr['CD1_1']
         dyim          = hdr['CD2_2']
         distance_unit = hdr['DIST_CM']
         distance_in   = hdr['DISTANCE']
         obsx          = hdr['OBSX']
         obsy          = hdr['OBSY']
         obsz          = hdr['OBSZ']
         distance      = np.sqrt(obsx**2 + obsy**2 + obsz**2)
         nyim, nxim    = Itot.shape
         xp_max        = distance * np.sin(dxim*(nxim/2.0) * (np.pi/180.0))
         yp_max        = distance * np.sin(dyim*(nyim/2.0) * (np.pi/180.0))
         ext1          = [-xp_max,xp_max,-yp_max,yp_max]
         I_cub = np.zeros((nf,nyim,nxim))
      I_cub[i,:,:] = Itot[:,:]
      hdu_obs.close()

   # Normalize
   I_cub = I_cub / 10.0**vmax
   #======================

   #--------
   fig, ax = plt.subplots(1,1,figsize=(3.5,3.5))

   #---- Image data
   I_img = ax.imshow(I_cub[0], animated=True, origin='lower',cmap=cmap, extent=ext1, vmin=0.0,vmax=1.0)

   ax.set_xticks([])
   ax.set_yticks([])
   cax  = inset_axes(ax, width="5%", height="100%", loc="lower left",
                     bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax.transAxes, borderpad=0,)
   cbar = fig.colorbar(I_img, cax=cax)
   cbar.set_label(r'Brightness')
   cbar.ax.minorticks_off()
   I_text = ax.set_title(r'\alpha $=%3.0f^\circ$' % alpha[i])

   #--------
   def init():
       I_img.set_array(np.zeros(I_cub[0].shape))
       return I_img,
   def updatefig(iframe,I_img):
       I_img.set_array(I_cub[iframe])
       label = r'\alpha $=%3.0f^\circ$'% alpha[iframe]
       I_text.set_text(label)
       return I_img,
   #--------

   #---- Video resolution
   w_in_inches, h_in_inches = fig.get_size_inches()
   w_resol = 1280
   dpi     = np.int32(w_resol/w_in_inches)
   print('  width in inches    = ', w_in_inches)
   print('  dot per inch (dpi) = ', dpi)
   print('  width resolution   = ', w_resol)

   plt.subplots_adjust(wspace=0.25, hspace=0.25, left=0.1,right=0.9,bottom=0.1,top=0.9)

   ani = animation.FuncAnimation(fig, updatefig, init_func=init, fargs=(I_img,),
                                 frames=nf, interval=200, blit=True, repeat_delay=2000)

   #-- make mp4 file.
   #-- (to do this, you need to install "ffmpeg". "sudo port install ffmpeg" in macos, "sudo apt install ffmpeg" in Linux.
   Writer   = animation.writers['ffmpeg']
   metadata = dict(title=r'%s' % fname0, artist='K.-I. Seon')
   duration = 10.0   # in sec
   #fps      = np.int32(nf/duration)
   fps      = nf/duration
   writer   = Writer(metadata=metadata, fps=fps)
   #ani.save(fname0+'_ani.mp4', dpi=dpi,writer=writer)
   ani.save(fname0+'_ani.mp4', writer=writer)

   #-- or show
   #plt.show()

if __name__ == "__main__":
   base_name = 'star_planet'
   make_animation(base_name)
