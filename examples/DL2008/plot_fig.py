#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from   astropy.io import fits
from   matplotlib.ticker import AutoMinorLocator, MultipleLocator

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':14})
plt.rcParams.update({'lines.linewidth':0.9})
fontsize = 11

ps_file = 'fig_test.ps'
fig, ax = plt.subplots(1,3,figsize=(10,3.5))
#fig, ax = plt.subplots(1,3,figsize=(10,3))

dir = '../DL2008/'
fname_arr = ['DL20e','DL19e']

x1 =  60.0
x2 = -100.0

cc  = ['k','r','k','b','r','g','y','k','b','r','g','y']
arr = np.arange(len(fname_arr))
for idx in arr:
   if idx < 2 :
      dir = '../DL2008_FS/'
   else:
      dir = '../DL2008_FS2/'

   fname  = fname_arr[idx]
   hdu1   = fits.open(dir+fname+'.fits.gz')
   hdu2   = fits.open(dir+fname+'_obs.fits.gz')
   hdu3   = fits.open(dir+fname+'_stokes.fits.gz')

   im     = hdu2[0].data.sum(axis=2) + hdu2[0].data.sum(axis=2)
   Fscatt = hdu2[0].data.sum(axis=(0,1))
   Fdirec = hdu2[1].data.sum(axis=(0,1))
   Ftot   = Fscatt + Fdirec

   r      = hdu3[4].data['r'].squeeze()
   inten  = hdu3[4].data['I'].squeeze()
   Q      = hdu3[4].data['Q'].squeeze()
   U      = hdu3[4].data['U'].squeeze()
   pol    = np.sqrt(Q**2 + U**2)/inten * 100.0

   xfreq  = hdu1[1].data['xfreq'].squeeze()
   Jin    = hdu1[1].data['Jin'].squeeze()

   hdr1   = hdu1[1].header
   hdr2   = hdu2[0].header

   N_HI   = hdr1['N_HI']
   temp   = hdr1['TEMP']
   vtherm = 0.12843374*np.sqrt(temp)
   dxfreq = hdr2['DXFREQ']
   dxim   = hdr2['CD1_1']
   dyim   = hdr2['CD2_2']
   nxim   = hdr2['NAXIS2']
   nyim   = hdr2['NAXIS3']
   dist   = hdr2['DISTANCE']
   npix   = np.any(im > 0.0)
   area   = dist**2 * np.tan(dxim*np.pi/180.0)*np.tan(dyim*np.pi/180.0) * npix
   scal   = 4.0*np.pi*area
   areaJ  = 4.0*np.pi
   scalJ  = 2.0*np.pi * areaJ

   Ftot   = Ftot * scal * 100.0
   Jin    = Jin  * scalJ * 100.0
   #ymax   = np.max(Ftot)*1.1
   ymax   = 0.04 * 100.0

   if idx == 0:
      ax[0].plot(xfreq, Jin, 'b')
      ax[0].set_xlim(x1,x2)
      ax[0].set_ylim(0.0,ymax)
      ax[0].xaxis.set_minor_locator(MultipleLocator(10))
      ax[0].xaxis.set_major_locator(MultipleLocator(50))
      ax[0].set_xlabel(r'Frequency ($x$)')
      ax[0].set_ylabel(r'Ly$\alpha$ flux')
      ax[0].text(-10, 0.035*100, r'$N_{\rm HI} = 10^{19}$ cm$^{-2}$', va='bottom',ha='left',  color='r', fontsize=fontsize)
      ax[0].text(-10, 0.032*100, r'$N_{\rm HI} = 10^{20}$ cm$^{-2}$', va='bottom',ha='left',  color='k', fontsize=fontsize)
      ax[0].text( 20, 0.015*100, r'input',                            va='bottom',ha='right', color='b', fontsize=fontsize)
      ax2 = ax[0].twiny()
      ax2.set_xlim(-x1*vtherm, -x2*vtherm)
      ax2.xaxis.set_major_locator(MultipleLocator(500))
      ax2.set_xlabel(r'Velocity [km/s]', fontsize=fontsize)
      ax2.xaxis.set_tick_params(labelsize=fontsize)

      ax[0].plot(np.array([0,7])-40.0, np.array([1,1])*0.030*100, 'k')
      ax[0].plot(np.array([0,7])-40.0, np.array([1,1])*0.027*100, 'k', linestyle='', marker=".", markersize=3.0)
      ax[0].text(-50.0, 0.030*100, r'w/o dust', va='center',ha='left', color='k', fontsize=fontsize)
      ax[0].text(-50.0, 0.027*100, r'w\ \ \ \ dust', va='center',ha='left', color='k', fontsize=fontsize)

      ax[1].set_xlim(0,1)
      ax[1].set_yscale('log')
      ax[1].set_ylim(0.03,3)
      ax[1].xaxis.set_minor_locator(MultipleLocator(0.1))
      ax[1].xaxis.set_major_locator(MultipleLocator(0.2))
      ax[1].set_xlabel(r'$\alpha/\alpha_{\rm max}$')
      ax[1].set_ylabel(r'Surface Brightness')
      scal_SB = 1.0/max(inten[1:])

      ax[2].set_xlim(0,1)
      ax[2].set_ylim(0,50.0)
      ax[2].xaxis.set_minor_locator(MultipleLocator(0.1))
      ax[2].xaxis.set_major_locator(MultipleLocator(0.2))
      ax[2].set_xlabel(r'$\alpha/\alpha_{\rm max}$')
      ax[2].set_ylabel(r'Polarization (\%)')

   if idx <= 1:
      #ax[0].plot(xfreq, Ftot,          color=cc[idx], drawstyle='steps-mid')
      ax[0].plot(xfreq, Ftot,          color=cc[idx])
      ax[1].plot(r,     inten*scal_SB, color=cc[idx])
      ax[2].plot(r,     pol,           color=cc[idx])
   else:
      ax[0].plot(xfreq, Ftot,          linestyle='',color=cc[idx], marker=".", markersize=3.0, markevery=2)
      ax[1].plot(r,     inten*scal_SB/1.1, linestyle='',color=cc[idx], marker=".", markersize=3.0, markevery=2)
      ax[2].plot(r,     pol,           linestyle='',color=cc[idx], marker=".", markersize=3.0, markevery=2)

plt.tight_layout()
plt.savefig(ps_file, papertype = 'a3')
plt.show()
