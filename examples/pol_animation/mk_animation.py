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
from   astropy.convolution import convolve
from   astropy.convolution import Gaussian2DKernel
from   square_subplots import square_subplots

#cmap  = 'viridis'
cmap  = 'jet'

def make_animation(argv):
   fname = argv
   np.seterr(divide='ignore', invalid='ignore')

   hdu1  = fits.open(fname+'.fits.gz')
   hdu3  = fits.open(fname+'_stokes.fits.gz')

   xfreq    = hdu1[1].data['Xfreq'].squeeze()
   Jout     = hdu1[1].data['Jout'].squeeze()
   #Vexp     = hdu1[1].header['Vexp']
   w        = (np.where(Jout == np.amax(Jout)))[0][0]
   freq_min = np.interp(np.amax(Jout)*1e-3, Jout[:w], xfreq[:w])
   freq_max = np.interp(np.amax(Jout)*1e-3, np.flip(Jout[w:]), np.flip(xfreq[w:]))
   k1       = np.amax(np.where(xfreq <= freq_min))
   k2       = np.amin(np.where(xfreq >= freq_max))

   Iim     = hdu3[0].data
   Qim     = hdu3[1].data
   Uim     = hdu3[2].data
   hdu1.close()
   hdu3.close()

   I_cub   = Iim
   pol_cub = np.sqrt(Qim**2 + Uim**2)/Iim
   w = np.where(Iim <= 0.0)
   pol_cub[w] = 0.0

   I_cub   = np.transpose(I_cub,(2,0,1))
   pol_cub = np.transpose(pol_cub,(2,0,1))
   I_cub   = I_cub[k1:k2+1,:,:]
   pol_cub = pol_cub[k1:k2+1,:,:]
   freq    = xfreq[k1:k2+1]
   Jout    = Jout[k1:k2+1]

   #--- convolution & normalization of Image data cube
   sigma  = 1.0
   kernel = Gaussian2DKernel(sigma)
   for jj in np.arange(I_cub.shape[0]):
       I_cub[jj] = convolve(I_cub[jj], kernel)
       I_cub[jj] = I_cub[jj]/np.amax(I_cub[jj])
   #scale = np.median(I_cub[np.where(I_cub > 0.0)])
   #I_cub = np.arcsinh(I_cub/scale)
   #I_cub = I_cub/np.amax(I_cub)

   #--- convolution of polarization data cube
   sigma  = 2.0
   kernel = Gaussian2DKernel(sigma)
   for jj in np.arange(pol_cub.shape[0]):
       pol_cub[jj] = convolve(pol_cub[jj], kernel)

   #--------
   fig, ax = plt.subplots(1,3,figsize=(10.5,3.5))

   #---- spectrum
   Jout  = Jout/np.amax(Jout)
   ax[0].plot(freq, Jout)
   ax[0].set_xlim(np.amin(freq),np.amax(freq))
   ax[0].set_ylim(0.0,np.amax(Jout)*1.15)
   ax[0].set_xlabel(r'$x$ (frequency)')
   ax[0].set_ylabel(r'$J_x$')
   x1, x2 = ax[0].get_xlim()
   y1, y2 = ax[0].get_ylim()
   ax[0].text(x1+(x2-x1)*0.05, y2-(y2-y1)*0.1, r'%s' % fname, va='center', ha='left', fontsize=12)

   Jpoint, = ax[0].plot([freq[0]], [Jout[0]], marker='o', markersize=4,color='red')

   #---- Image data
   vmin  = 0.0
   vmax  = np.amax(I_cub)
   I_img = ax[1].imshow(I_cub[0], animated=True, origin='lower',cmap=cmap,
                         extent=[-1,1,-1,1], vmin=vmin,vmax=vmax)

   ax[1].set_xticks([])
   ax[1].set_yticks([])
   ax[1].set_title(r'$J$ ($J^{\rm max} = 1$ at each $x$)', fontsize=12)
   cax  = inset_axes(ax[1], width="5%", height="100%", loc="lower left",
                     bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax[1].transAxes, borderpad=0,)
   cbar = fig.colorbar(I_img, cax=cax, ticks=np.arange(6)*0.2)
   cbar.ax.minorticks_off()
   I_text = ax[1].text(-0.95,-0.95, r'$\mathbf{x = %4.1f}$' % freq[0], fontsize=14, color='white', va='bottom',ha='left')

   #---- polarization data
   vmin  = 0.0
   vmax  = np.amax(pol_cub)
   pol_img = ax[2].imshow(pol_cub[0], animated=True, origin='lower',cmap=cmap,
                          extent=[-1,1,-1,1], vmin=vmin,vmax=vmax)

   ax[2].set_xticks([])
   ax[2].set_yticks([])
   ax[2].set_title(r'Polarization', fontsize=12)
   cax  = inset_axes(ax[2], width="5%", height="100%", loc="lower left",
                     bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax[2].transAxes, borderpad=0,)
   cbar = fig.colorbar(pol_img, cax=cax, ticks=np.arange(6)*0.2)
   cbar.ax.minorticks_off()
   pol_text = ax[2].text(-0.95,-0.95, r'$x = %4.1f$' % freq[0], fontsize=14, color='white', va='bottom',ha='left')

   #--------
   def init():
       Jpoint.set_data([],[])
       I_img.set_array(np.zeros(I_cub[0].shape))
       pol_img.set_array(np.zeros(pol_cub[0].shape))
       return (I_img,pol_img)

   def updatefig(i,I_img,pol_img):
       x = [freq[i]]
       y = [Jout[i]]
       Jpoint.set_data(x,y)
       I_img.set_array(I_cub[i])
       pol_img.set_array(pol_cub[i])
       label = r'$x = %4.1f$' % freq[i]
       I_text.set_text(label)
       pol_text.set_text(label)
       return I_img,pol_img
   #--------

   #---- Video resolution
   w_in_inches, h_in_inches = fig.get_size_inches()
   w_resol = 1280
   dpi     = np.int32(w_resol/w_in_inches)
   print('  width in inches    = ', w_in_inches)
   print('  dot per inch (dpi) = ', dpi)
   print('  width resolution   = ', w_resol)

   plt.subplots_adjust(wspace=0.35, hspace=0.25, left=0.1,right=0.9,bottom=0.1,top=0.9)
   s = square_subplots(fig)

   ani = animation.FuncAnimation(fig, updatefig, init_func=init, fargs=(I_img,pol_img),
                                 frames=I_cub.shape[0], interval=200, blit=True, repeat_delay=2000)

   #-- make mp4 file.
   #-- (to do this, you need to install "ffmpeg". "sudo port install ffmpeg" in macos, "sudo apt install ffmpeg" in Linux.
   Writer = animation.writers['ffmpeg']
   metadata = dict(title=r'%s' % fname, artist='K.-I. Seon')
   duration = 10.0   # in sec
   fps    = np.int32(I_cub.shape[0]/duration)
   writer = Writer(metadata=metadata, fps=fps)
   ani.save(fname+'_ani.mp4', dpi=dpi,writer=writer)

   #-- or show
   #plt.show()

if __name__ == "__main__":
   opts, args = getopt.getopt(sys.argv[1:], 'ai:', ['all'])
   if len(opts) == 0 and len(args) == 0:
      print('')
      print('Usage:')
      print('   mk_animation.py --all or')
      print('   mk_animation.py -i basic_file_name')
      print('')
      sys.exit(2)

   #-- why this does not work.
   #try:
   #   opts, args = getopt.getopt(sys.argv[1:], 'ai:', ['all'])
   #except getopt.GetoptError as err:
   #   print('"mk_anim.py --all" or "mk_anim.py -i basic_file_name"')
   #   sys.exit(2)

   for opt, arg in opts:
      if opt in ("-a","--all"):
         flist = glob.glob('*_stokes.fits.gz')
         flist.sort()
         for fname in flist:
            i0        = fname.find('_stokes.fits.gz')
            base_name = fname[0:i0]
            print('\n---> making animation for %s' % base_name)
            make_animation(base_name)
      elif opt in ("-i"):
         base_name = arg
         print('\n---> making animation for %s' % base_name)
         make_animation(base_name)
