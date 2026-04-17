#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import copy
import matplotlib as mpl
from   astropy.io import fits
from   matplotlib.ticker import AutoMinorLocator, MultipleLocator
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
from   square_subplots import square_subplots
from   scipy.ndimage import gaussian_filter

np.seterr(divide='ignore', invalid='ignore')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':14})
plt.rcParams.update({'lines.linewidth':0.9})
fontsize = 11

cmap  = 'seismic'

#cmap0 = mpl.cm.nipy_spectral.reversed()
#cmap0 = mpl.cm.jet
cmap0 = copy.copy(mpl.cm.get_cmap("jet"))
cmap0.set_bad('white')

pdf_file = 'pol_map.pdf'
fig, ax = plt.subplots(1,3,figsize=(10,3.5))

dir       = './'
#fname_arr = ['FeII_UV3']
fname_arr = ['FeII_UV3_V300']

arr = np.arange(len(fname_arr))
for idx in arr:
   fname  = fname_arr[idx]
   hdu1   = fits.open(dir+fname+'.fits.gz')
   hdu3   = fits.open(dir+fname+'_stokes.fits.gz')

   wavl   = hdu1[1].data['WAVELENGTH']
   vel    = hdu1[1].data['VELOCITY']
   #wvel,  = np.where((vel > 2800.0) & (vel < 3200.0))
   #wvel,  = np.where((vel > 1200.0) & (vel < 1600.0))
   # For FeII UV3
   wvel,  = np.where((vel > 3400.0) & (vel < 3800.0))
   # For SiII
   #wvel,  = np.where(wavl > 1263.0)

   hdr1   = hdu1[1].header
   tau0   = hdr1['TAUMAX']
   Ngas   = hdr1['NGASMAX']
   Vexp   = hdr1['VEXP']
   temp   = hdr1['TEMP']
   vtherm = 0.12843374*np.sqrt(temp)

   #--- Image for Q and U
   print(hdu3[0].data[:,:,wvel].shape)
   Iim = hdu3[0].data[:,:,wvel].sum(axis=2)
   Qim = hdu3[1].data[:,:,wvel].sum(axis=2)
   Uim = hdu3[2].data[:,:,wvel].sum(axis=2)

   #---- First, find zero locations.
   w   = np.where(Iim == 0.0)

   #---- Second, smoothing
   gauss_sigma = 1.0
   Iim = gaussian_filter(Iim, sigma=gauss_sigma, mode='nearest')
   Qim = gaussian_filter(Qim, sigma=gauss_sigma, mode='nearest')
   Uim = gaussian_filter(Uim, sigma=gauss_sigma, mode='nearest')

   #---- Third, normalize
   Qim = Qim/Iim
   Uim = Uim/Iim

   if w[0].size > 0:
      Iim[w] = np.nan
      Qim[w] = 0.0
      Uim[w] = 0.0

   theta  = 0.5 * np.arctan2(Uim, Qim)
   polim    = np.sqrt(Qim**2 + Uim**2) * 100.0
   polim[w] = np.nan

   #--- Radial Profile for polarization
   r      = hdu3[4].data['R']
   pol    = hdu3[4].data['POL']
   pmax   = np.int32(np.int32(np.amax(pol)*100)/5.0+1)*5.0/100.0
   print('pmax = %10.6e' % pmax)

   Imin  = 0.0
   #Imax  = np.amax(Iim[np.where(Iim > 0.0)])
   Icen  = np.median(Iim[np.where(Iim > 0.0)])
   Imax  = Icen*5.0
   levels= np.arange(17)*0.01 * Imax
   imgc  = ax[0].contour(Iim,  levels, extent=[-1,1,-1,1],colors='grey', origin='lower', vmin=Imin, vmax=Imax, linewidths=0.1)
   img0  = ax[0].contourf(Iim, levels, extent=[-1,1,-1,1],cmap=cmap0,    origin='lower', vmin=Imin, vmax=Imax)

   cax0  = inset_axes(ax[0], width="5%", height="100%", loc="lower left",
                      bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax[0].transAxes, borderpad=0,)
   cbar0 = fig.colorbar(img0, cax=cax0, ticks=levels)
   cbar0.ax.set_yticklabels(['0.00','','','','','0.05','','','','','0.10','','','','','0.15',''])
   cbar0.ax.minorticks_off()
   cbar0.ax.tick_params(labelsize=12, color='grey', length=7.5, width=0.1)

   #cbar0 = fig.colorbar(img0, cax=cax0, ticks=np.arange(5)*0.05)
   #cbar0.ax.minorticks_on()
   #cbar0.ax.tick_params(labelsize=12)

   #lab = r'$I$ $(T = 10^{%d} {\rm K}, \tau_0 = 10^{%d})$' % (np.log10(temp), np.log10(tau0))
   #ax[0].set_title(lab, fontsize=12)
   ax[0].xaxis.set_major_locator(MultipleLocator(0.2))
   ax[0].yaxis.set_major_locator(MultipleLocator(0.2))
   ax[0].xaxis.set_minor_locator(MultipleLocator(0.1))
   ax[0].yaxis.set_minor_locator(MultipleLocator(0.1))
   ax[0].xaxis.set_ticklabels([])
   ax[0].yaxis.set_ticklabels([])

   img1  = ax[1].imshow(Qim, extent=[-1,1,-1,1],cmap=cmap, origin='lower', vmin=-pmax, vmax=pmax)
   cax1  = inset_axes(ax[1], width="5%", height="100%", loc="lower left",
                      bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax[1].transAxes, borderpad=0,)
   cbar1 = fig.colorbar(img1, cax=cax1)
   cbar1.ax.tick_params(labelsize=12)
   ax[1].set_title(r'$Q/I$', fontsize=12)
   ax[1].xaxis.set_major_locator(MultipleLocator(0.2))
   ax[1].yaxis.set_major_locator(MultipleLocator(0.2))
   ax[1].xaxis.set_minor_locator(MultipleLocator(0.1))
   ax[1].yaxis.set_minor_locator(MultipleLocator(0.1))
   ax[1].xaxis.set_ticklabels([])
   ax[1].yaxis.set_ticklabels([])

   img2  = ax[2].imshow(Uim, extent=[-1,1,-1,1],cmap=cmap, origin='lower', vmin=-pmax, vmax=pmax)
   cax2  = inset_axes(ax[2], width="5%", height="100%", loc="lower left",
                      bbox_to_anchor=(1.05, 0.0, 1, 1), bbox_transform=ax[2].transAxes, borderpad=0,)
   cbar2 = fig.colorbar(img2, cax=cax2)
   cbar2.ax.tick_params(labelsize=12)
   ax[2].set_title(r'$U/I$', fontsize=12)
   ax[2].xaxis.set_major_locator(MultipleLocator(0.2))
   ax[2].yaxis.set_major_locator(MultipleLocator(0.2))
   ax[2].xaxis.set_minor_locator(MultipleLocator(0.1))
   ax[2].yaxis.set_minor_locator(MultipleLocator(0.1))
   ax[2].xaxis.set_ticklabels([])
   ax[2].yaxis.set_ticklabels([])

   xvec = -polim * np.sin(theta)
   yvec =  polim * np.cos(theta)

   xmax = 1.0
   ymax = 1.0
   xmin = -xmax
   ymin = -ymax
   nx   = np.size(Iim, 0)
   ny   = np.size(Iim, 1)
   dx   = 2.0*xmax/nx
   dy   = 2.0*ymax/ny
   x    = np.linspace(xmin+dx/2.0,xmax-dx/2.0,nx)
   y    = np.linspace(ymin+dy/2.0,ymax-dy/2.0,ny)

   scale = 150.0
   step  = 10
   j0    = np.int32(step/2)
   if (j0/2)*2 != 2:
      j0 = j0-1

   if (nx/2)*2 != nx:
      nxcen = np.int32((nx-1)/2)
      id1   = np.concatenate([np.flip(np.arange(nxcen,0,-step)), np.arange(nxcen+step,nx,step)])
   else:
      id1   = np.arange(j0,nx,step)

   if (ny/2)*2 != ny:
      nycen = np.int32((ny-1)/2)
      id2   = np.concatenate([np.flip(np.arange(nycen,0,-step)), np.arange(nycen+step,ny,step)])
   else:
      id2   = np.arange(j0,ny,step)

   width = 0.007
   x    = x[id1]
   y    = y[id2]
   xvec = xvec[id1,:]
   yvec = yvec[id1,:]
   xvec = xvec[:,id2]
   yvec = yvec[:,id2]
   ax[0].quiver(x,y,xvec,yvec,headaxislength=0,headlength=0,headwidth=0,pivot='mid',scale=scale,width=width,color='k')

   ax[0].xaxis.set_major_locator(MultipleLocator(0.2))
   ax[0].yaxis.set_major_locator(MultipleLocator(0.2))
   ax[0].xaxis.set_minor_locator(MultipleLocator(0.1))
   ax[0].yaxis.set_minor_locator(MultipleLocator(0.1))
   ax[0].xaxis.set_ticklabels([])
   ax[0].yaxis.set_ticklabels([])

   pol_ref = 20.0
   x0  = np.array([-0.92])
   y0  = np.array([-0.9])
   xv0 = np.array([pol_ref])
   yv0 = np.array([0.0])
   ax[0].quiver(x0,y0,xv0,yv0,headaxislength=0,headlength=0,headwidth=0,pivot='tail',scale=scale,width=width,color='k')
   lab = r'%d%s' % (pol_ref, '\%')
   ax[0].text(x0[0],y0[0]+0.03,lab,fontsize=9, va='bottom', ha='left', color='k')

plt.subplots_adjust(wspace=0.4, hspace=0.2, left=0.1,right=0.9,bottom=0.1,top=0.9)
s = square_subplots(fig)
plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
#plt.show()
