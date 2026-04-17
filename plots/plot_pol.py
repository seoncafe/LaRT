#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import copy
from   astropy.io import fits
from   matplotlib.ticker import AutoMinorLocator, MultipleLocator
from   mpl_toolkits.axes_grid1.inset_locator import inset_axes
from   square_subplots import square_subplots
from   scipy.ndimage import gaussian_filter

# Do not print out warning for a zero by zero division.
np.seterr(divide='ignore', invalid='ignore')

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':14})
plt.rcParams.update({'lines.linewidth':0.9})
plt.rcParams.update({'legend.fontsize':9})
fontsize = 12

fig   = plt.figure(figsize=(4.2,4.2))
fig_w, fig_h = fig.get_size_inches()

gap_w = 0.22
gap_h = 0.12
box_w = 1.0 - gap_w
box_h = box_w

#---
#dir       = '../DL2008_FS_allph/'
#fname_arr = ['DL19e','DL20e','DL20e_dust']

#---
dir       = '../sphere_allph/'
fname_arr = ['t4tau4','t4tau5','t4tau6','t4tau7']
#fname_arr = ['t1tau4','t1tau5','t1tau6','t1tau7']

#---
#dir       = '../vel_effect_allph/'
#fname_arr = ['t4NHI2_20_V0020','t4NHI2_20_V0200','t4NHI2_20_V0500','t4NHI2_20_V2000']
#fname_arr = ['t4NHI2_20_V0020','t4NHI2_20_V0100','t4NHI2_20_V0200','t4NHI2_20_V0500',
#             't4NHI2_20_V1000','t4NHI2_20_V1500','t4NHI2_20_V2000','t4NHI2_20_V3000']

#---
#dir       = '../SSH_MUSE/'
#fname_arr = ['MUSE1185', 'MUSE0082', 'MUSE6905', 'MUSE1343', 'MUSE0053', 'MUSE0171', 'MUSE0547', 'MUSE0364']

nf    = len(fname_arr)
#rbin1 = np.array([0.0, 0.1, 0.2, 0.4, 0.6, 0.8])
#rbin2 = np.array([0.1, 0.2, 0.4, 0.6, 0.8, 1.0])
#nr    = rbin1.size

nr   = 20
rbin = np.arange(nr)/(nr-1.0)
nx   = 10
cc = [plt.cm.rainbow(i) for i in np.linspace(0,1,nx)]
#cc = [plt.cm.rainbow(i) for i in np.linspace(0,1,nr)]
#cc = [plt.cm.winter(i) for i in np.linspace(0,1,nr)]

if fname_arr[0][0:5] == 't4NHI': pdf_file = 'fig_pol_vel.pdf'
if fname_arr[0][0:5] == 't4tau': pdf_file = 'fig_pol_sph4.pdf'
if fname_arr[0][0:5] == 't1tau': pdf_file = 'fig_pol_sph1.pdf'
if fname_arr[0][0:2] == 'DL':    pdf_file = 'fig_pol_DL2008.pdf'
if fname_arr[0][0:2] == 'MU':    pdf_file = 'fig_pol_muse.pdf'

arr = np.arange(nf)
for idx in arr:
   fname  = fname_arr[idx]
   hdu1   = fits.open(dir+fname+'.fits.gz')
   hdu4   = fits.open(dir+fname+'_allph.fits.gz')

   hdr1   = hdu1[1].header
   tau0   = hdr1['TAUMAX']
   try:
      N_HI   = hdr1['N_HI']
   except:
      N_HI   = hdr1['N_HIMAX']
   Vexp   = hdr1['VEXP']
   temp   = hdr1['TEMP']
   vtherm = 0.12843374*np.sqrt(temp)
   xJ     = hdu1[1].data['xfreq'].squeeze()
   #Jin    = hdu1[1].data['Jin'].squeeze()
   Jout   = hdu1[1].data['Jout'].squeeze()
   Jout   = Jout/np.amax(Jout)

   xfreq  = hdu4[1].data['xfreq2'].squeeze()
   r      = hdu4[1].data['r'].squeeze()

   #nscatt = hdu4[1].data['nscatt_HI'].squeeze()
   #Iph    = hdu4[1].data['I'].squeeze()
   Q      = hdu4[1].data['Q'].squeeze()
   U      = hdu4[1].data['U'].squeeze()

   #lab1 = r'$N_{\rm HI} = 2\times10^{20}$ cm$^{-2}$'
   #lab2 = r'$V_{\rm exp} = %d$ km/s' % (Vexp)
   if fname == 't4NHI2_20_V3000':
      x1 = -150.0
      x2 =   5.0
   if fname == 't4NHI2_20_V2000':
      x1 = -320.0
      x2 =   10.0
   if fname == 't4NHI2_20_V1500':
      x1 = -180.0
      x2 =   10.0
   if fname == 't4NHI2_20_V1000':
      x1 = -180.0
      x2 =   10.0
   if fname == 't4NHI2_20_V0500':
      x1 = -170.0
      x2 =   10.0
   if fname == 't4NHI2_20_V0200':
      x1 = -110.0
      x2 =   10.0
   if fname == 't4NHI2_20_V0100':
      x1 =  -85.0
      x2 =   10.0
   if fname == 't4NHI2_20_V0020':
      x1 =  -50.0
      x2 =   35.0

   if fname == 't4tau4':
      x1 =  -8.0
      x2 =   8.0
   if fname == 't4tau5':
      x1 =  -12.0
      x2 =   12.0
   if fname == 't4tau6':
      x1 =  -22.0
      x2 =   22.0
   if fname == 't4tau7':
      x1 =  -45.0
      x2 =   45.0

   if fname == 't1tau4':
      x1 =  -20.0
      x2 =   20.0
   if fname == 't1tau5':
      x1 =  -40.0
      x2 =   40.0
   if fname == 't1tau6':
      x1 =  -75.0
      x2 =   75.0
   if fname == 't1tau7':
      x1 =  -125.0
      x2 =   125.0

   if fname[0:5] == 'DL19e':
      x1 =  -100.0
      x2 =   60.0
   if fname[0:5] == 'DL20e':
      x1 =  -150.0
      x2 =   60.0

   if fname[0:4] == 'MUSE':
      x1 =  -120.0
      x2 =   40.0

   if (idx == 0):
      px = gap_w/2.0
      py = gap_h/2.0
      pw = box_w
      ph = box_h
   else:
      pos = ax.get_position()
      px  = pos.x1 + gap_w
      py  = gap_h/2.0
      pw  = box_w
      ph  = box_h

   ax = fig.add_axes([px,py,pw,ph])

   x1    = 0.0
   dx    = (x2-x1)/nx
   xbin1 = np.arange(nx)*dx + x1
   xbin2 = xbin1 + dx
   ymax  = -999.9
   #for ix in np.arange(nx-2):
   for ix in np.arange(nx):
      #w  = np.where((xfreq >= xbin1[ix]) & (xfreq < xbin2[ix]))
      w  = np.where((np.abs(xfreq) >= xbin1[ix]) & (np.abs(xfreq) < xbin2[ix]))
      nw = w[0].size
      yQ,     x_edges = np.histogram(r[w], bins=rbin, weights=Q[w])
      yU,     x_edges = np.histogram(r[w], bins=rbin, weights=U[w])
      yp              = np.sqrt(yQ**2 + yU**2)
      yy_num, x_edges = np.histogram(r[w], bins=rbin)
      wz   = np.where(yy_num > 0.0)
      yy   = np.zeros(yp.size)
      if wz[0].size > 0:
         yy[wz]   = yp[wz]/yy_num[wz]
      xx   = (x_edges[0:-1]+x_edges[1:])/2.0
      ymax = np.amax([np.amax(yy),ymax])

      print(xbin1[ix], xbin2[ix])
      lab  = r'$%3.1f < x < %3.1f$' % (xbin1[ix], xbin2[ix])
      ax.plot(xx, yy, color=cc[ix], drawstyle='steps-mid', label=lab)

   ax.plot(xJ, Jout, color='k')
   ax.set_xlabel(r'$r_{\rm p}/r_{\rm max}$')
   ax.set_ylabel(r'Polarization')
   ax.set_xlim(0.0,1.0)
   ax.set_ylim(0.0,np.amin([1.0,ymax*1.1]))
   ax.xaxis.set_major_locator(MultipleLocator(0.2))
   ax.legend(loc="upper left", frameon=False, handlelength=0.7, prop={'size':8}, labelspacing=0.2, handletextpad=0.3,
             labelcolor='linecolor')

plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
#plt.show()
