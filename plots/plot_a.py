#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from   astropy.io import fits
from   matplotlib.ticker import AutoMinorLocator, MultipleLocator

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':14})
plt.rcParams.update({'lines.linewidth':0.9})
fontsize = 12

pdf_file = 'fig08b.pdf'
fig, ax = plt.subplots(1,2,figsize=(7,3.5))

dir       = '../sphere_allph/'
#fname_arr = ['t4tau3','t4tau4','t4tau5','t4tau6','t4tau7']
fname_arr = ['t4tau3']

cc  = ['r','b','g','c','m']
arr = np.arange(len(fname_arr))
for idx in arr:
   fname  = fname_arr[idx]
   hdu1   = fits.open(dir+fname+'.fits.gz')
   hdu3   = fits.open(dir+fname+'_stokes.fits.gz')
   hdu4   = fits.open(dir+fname+'_allph.fits.gz')

   rarr   = hdu3[4].data['r'].squeeze()
   dr     = rarr[2] - rarr[1]
   r1     = rarr - dr/2.0
   r1[0]  = 0.0
   r2     = rarr + dr/2.0
   for k in np.arange(rarr.size):
      print(r1[k],r2[k])
