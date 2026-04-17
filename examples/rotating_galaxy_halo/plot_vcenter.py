#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Disable zero division warnings
np.seterr(divide='ignore')

cmap = cm.coolwarm

fname = 'rin0.1_Vrot100_NHI18'
#fname = 'rin0.1_Vrot300_NHI18'

#fname = 'rin0.1_Vrot100_NHI18_unif'
#fname = 'rin0.1_Vrot300_NHI18_unif'

c    = fits.open(fname+'.fits.gz')
vel  = c[1].data['velocity']
ii   = fname.index('Vrot')+4
vmin = -np.float64(fname[ii:ii+3])
vmax =  np.float64(fname[ii:ii+3])

c1 = fits.open(fname+'_obs_001.fits.gz')
c2 = fits.open(fname+'_obs_002.fits.gz')
c3 = fits.open(fname+'_obs_003.fits.gz')
c4 = fits.open(fname+'_obs_004.fits.gz')

#TEST
#ny,nx,nl = c1[0].data.shape
#im1 = np.zeros([ny,nx])
#for j in np.arange(ny):
#  for i in np.arange(nx):
#    im1[j,i] = np.sum(c1[0].data[j,i,:] * vel)/np.sum(c1[0].data[j,i,:])
#plt.imshow(im1)

im1 = np.sum(c1[0].data * vel[np.newaxis,np.newaxis,:], axis=2)/np.sum(c1[0].data, axis=2)
im2 = np.sum(c2[0].data * vel[np.newaxis,np.newaxis,:], axis=2)/np.sum(c2[0].data, axis=2)
im3 = np.sum(c3[0].data * vel[np.newaxis,np.newaxis,:], axis=2)/np.sum(c3[0].data, axis=2)
im4 = np.sum(c4[0].data * vel[np.newaxis,np.newaxis,:], axis=2)/np.sum(c4[0].data, axis=2)

im1 = np.swapaxes(im1, 0,1)
im2 = np.swapaxes(im2, 0,1)
im3 = np.swapaxes(im3, 0,1)
im4 = np.swapaxes(im4, 0,1)

#w1 = np.isinf(im1)
#w2 = np.isinf(im2)
#w3 = np.isinf(im3)
#w4 = np.isinf(im4)
#im1[w1] = 0.0
#im2[w2] = 0.0
#im3[w3] = 0.0
#im4[w4] = 0.0

#sigma = 0.5
#im1 = gaussian_filter(im1, sigma, mode='nearest')
#im2 = gaussian_filter(im2, sigma, mode='nearest')
#im3 = gaussian_filter(im3, sigma, mode='nearest')
#im4 = gaussian_filter(im4, sigma, mode='nearest')
#im1[w1] = np.inf
#im2[w2] = np.inf
#im3[w3] = np.inf
#im4[w4] = np.inf

beta1 = c1[0].header['beta']
beta2 = c2[0].header['beta']
beta3 = c3[0].header['beta']
beta4 = c4[0].header['beta']

fig, axs = plt.subplots(ncols=2,nrows=2,figsize=(10,10))
axs = axs.flat

ax  = axs[0]
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
pos = ax.imshow(im1, origin='lower', vmin=vmin,vmax=vmax, cmap=cmap)
div = make_axes_locatable(ax)
cax = div.append_axes("right", size="5%", pad=0.1)
fig.colorbar(pos,cax=cax)
ax.set_title(r'$\beta=%d^\circ$' % beta1)

ax  = axs[1]
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
pos = ax.imshow(im2, origin='lower', vmin=vmin,vmax=vmax, cmap=cmap)
div = make_axes_locatable(ax)
cax = div.append_axes("right", size="5%", pad=0.1)
fig.colorbar(pos,cax=cax)
ax.set_title(r'$\beta=%d^\circ$' % beta2)

ax  = axs[2]
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
pos = ax.imshow(im3, origin='lower', vmin=vmin,vmax=vmax, cmap=cmap)
div = make_axes_locatable(ax)
cax = div.append_axes("right", size="5%", pad=0.1)
fig.colorbar(pos,cax=cax)
ax.set_title(r'$\beta=%d^\circ$' % beta3)

ax  = axs[3]
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
pos = ax.imshow(im4, origin='lower', vmin=vmin,vmax=vmax, cmap=cmap)
div = make_axes_locatable(ax)
cax = div.append_axes("right", size="5%", pad=0.1)
fig.colorbar(pos,cax=cax)
ax.set_title(r'$\beta=%d^\circ$' % beta4)

#pdf_file = 'vcenter_map.pdf'
#plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
plt.show()
