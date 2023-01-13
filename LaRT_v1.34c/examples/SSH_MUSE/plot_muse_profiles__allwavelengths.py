#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from   astropy.io import fits
from   matplotlib.ticker import AutoMinorLocator, MultipleLocator
from   astropy.convolution import Moffat2DKernel, Gaussian1DKernel, convolve
from   radial_profile import radial_profile
from   square_subplots import square_subplots
import os

# --- 2020.10.14.
# This python code is to compare the best-fit models with the MUSE data.
# ---

# Do not print out warning for a zero by zero division.
np.seterr(divide='ignore', invalid='ignore')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size':10})
plt.rcParams.update({'lines.linewidth':0.9})
fontsize  = 9

cmap = 'seismic'

#-----------------
dir   = './'

#--- Read best-fit data
result_table = './Leclercq/muse_bestfit.txt'
idarr, rsUV,rsHI,rpeak,Vpeak,DelV,logtau,DGR,zarr, r_lim,lambda1,lambda2 = np.loadtxt(result_table, skiprows=3, unpack=True)

#--- Open a file to record the scaling factors for spectrum and surface brightness.
fileout = open('muse_scaling_factors_all.txt','w')

def plot_profiles(fname):
   #--- should be reset
   fig, ax = plt.subplots(1,3,figsize=(10.5,3.5))

   #--- set best-fit parameters and limits for wavelength and radial size.
   museid    = int(fname[4:])
   wid       = np.where(idarr == museid)[0]
   redshift  = zarr[wid]
   rscale_UV = rsUV[wid]
   lam1      = lambda1[wid]
   lam2      = lambda2[wid]
   rlim      = r_lim[wid]
   DGR1      = DGR[wid]

   #--- read simulation data
   hdu1   = fits.open(dir+fname+'.fits.gz')
   hdu2   = fits.open(dir+fname+'_stokes_2D.fits.gz')
   hdu3   = fits.open(dir+fname+'_stokes.fits.gz')

   #--- Parameters
   hdr1   = hdu1[1].header
   hdr2   = hdu2[0].header
   hdr3   = hdu3[0].header
   temp   = hdr1['TEMP']
   #---
   dist   = hdr3['DISTANCE']
   nxfreq = hdr3['NAXIS1']
   #nx     = hdr3['NAXIS2']
   #ny     = hdr3['NAXIS3']
   nx     = hdr2['NAXIS1']
   ny     = hdr2['NAXIS2']
   dxfreq = hdr1['DXFREQ']
   dx     = hdr3['CD1_1'] * np.pi/180.
   dy     = hdr3['CD2_2'] * np.pi/180.
   vtherm = 0.12843374*np.sqrt(temp)
   nrmax  = (nx+1)/2 - 0.5
   nrmax1 = np.int32(nrmax) + 1

   #--- Wavelength
   wavl0  = 1215.67
   speedc = 2.99792458e5
   resol  = 3000.0
   fwhm   = wavl0*(1.0 + redshift)/resol
   stddev = fwhm / (2.0 * np.sqrt(2.0*np.log10(2.0)))

   xfreq  = hdu1[1].data['Xfreq'].squeeze()
   wavl   = wavl0 * (1.0 + (-xfreq * vtherm)/speedc)
   wavlz  = wavl*(1.0 + redshift)
   dwavl  = abs(wavl[1]  - wavl[0])
   dwavlz = abs(wavlz[1] - wavlz[0])

   #---- image PSF (Moffat) and spectrum LSF (Gaussian) kernel
   moffat_pow   = 2.8
   rmax_arcsec  = 5.0
   moffat_fwhm  = 0.7
   moffat_width = ((moffat_fwhm/rmax_arcsec) /2.0)/np.sqrt(2.0**(1.0/moffat_pow)-1.0)
   moffat_width_pix = moffat_width * nrmax

   #kern_size = nx/2
   kern_size = np.int32(nx/2)
   if (kern_size/2)*2 == kern_size: kern_size = kern_size + 1
   kernel = Moffat2DKernel(moffat_width_pix, moffat_pow, x_size=kern_size, y_size=kern_size)
   gauss  = Gaussian1DKernel(stddev/dwavlz)

   #--- Image
   I_im   = hdu2[0].data
   Q_im   = hdu2[1].data
   U_im   = hdu2[2].data

   Idat   = hdu3[0].data
   Qdat   = hdu3[1].data
   Udat   = hdu3[2].data
   Itot   = np.sum(Idat) * dxfreq * (dx*dy)*dist**2 * 4.0*np.pi
   print('--- %04d: DGR = %3.1f' % (museid, DGR1))
   print('   Itot = %10.5e  (Note this should be ~ 1 if DGR = 0)' % Itot)
   fileout.write('--- %04d: DGR = %3.1f\n' % (museid, DGR1))
   fileout.write('   Itot = %10.5e  (Note this should be ~ 1 if DGR = 0)\n' % Itot)

   #wlam   = np.where((wavlz > lam1) & (wavlz < lam2))[0]
   ##I_im   = Idat[:,:,wlam].sum(axis=2)
   ##Q_im   = Qdat[:,:,wlam].sum(axis=2)
   ##U_im   = Udat[:,:,wlam].sum(axis=2)
   #I_im   = Idat.sum(axis=2)
   #Q_im   = Qdat.sum(axis=2)
   #U_im   = Udat.sum(axis=2)
   pol_im = np.sqrt(Q_im**2 + U_im**2)/I_im
   wz     = np.where(I_im <= 0.0)
   pol_im[wz] = 0.0

   #---- input radial profile
   center   = [nx/2.0, ny/2.0]
   y, x     = np.indices((ny,nx)) + 0.5
   r        = np.sqrt((x - center[0])**2 + (y - center[1])**2)
   r        = r/nrmax
   Iin      = np.exp(-r/rscale_UV)
   wrz      = np.where(r > 1.0)
   Iin[wrz] = 0.0

   #--- spectrum
   #Jout = np.zeros(nxfreq)
   #wr   = np.where(r <= rlim)
   #for idx in np.arange(nxfreq):
   #    Itmp      = Idat[:,:,idx]
   #    Jout[idx] = Itmp[wr].sum()
   Jout   = Idat.sum(axis=(0,1))

   spec   = convolve(Jout, gauss)

   wavl   = np.flip(wavl)
   wavlz  = np.flip(wavlz)
   Jout   = np.flip(Jout)
   spec   = np.flip(spec)
   gauss  = np.flip(gauss)

   #---- radial profiles
   I_conv   = convolve(I_im, kernel)
   Q_conv   = convolve(Q_im, kernel)
   U_conv   = convolve(U_im, kernel)
   pol_conv = np.sqrt(Q_conv**2 + U_conv**2)/I_conv
   Iin_conv = convolve(Iin, kernel)

   rr, prof_I     = radial_profile(I_im)
   rr, prof_Iconv = radial_profile(I_conv)
   rr, prof_pol0  = radial_profile(pol_im)
   rr, prof_pol   = radial_profile(pol_conv)
   rr, prof_Iin   = radial_profile(Iin_conv)
   raxis          = rr
   prof_Iin       = prof_Iin / np.sum(prof_Iin) * np.sum(prof_Iconv)

   # Read Spectrum and SB data.
   dir1  = './Leclercq/'
   fspec = dir1+'%04d_spec.txt' % museid
   fSB   = dir1+'%04d_SB.txt'   % museid
   lam_obs, spec_obs, spec_err = np.loadtxt(fspec, skiprows=0, unpack=True)
   r_obs,   lnSB_obs, lnSB_err = np.loadtxt(fSB,   skiprows=0, unpack=True)
   lam_res = lam_obs / (1.0 + redshift)
   r_obs   = r_obs / rmax_arcsec
   SB_err  = 10.0**(lnSB_obs + lnSB_err) - 10.0**lnSB_obs
   SB_up   = 10.0**(lnSB_obs + lnSB_err)
   SB_low  = 10.0**(lnSB_obs - lnSB_err)
   SB_obs  = 10.0**lnSB_obs

   #--- Scale the model data to the physical units of the observational data.
   #--- We need to place the target to the distance of the model for SB.
   Jout       = Jout * (dxfreq / dwavlz) * (dx*dy*dist**2) * (4.0*np.pi)
   spec       = spec * (dxfreq / dwavlz) * (dx*dy*dist**2) * (4.0*np.pi)
   prof_I     = prof_I     * dxfreq * (4.0*np.pi * dx*dy*dist**2) / (4.0*np.pi*dist**2 *(180.*60.*60./np.pi)**2 )
   prof_Iconv = prof_Iconv * dxfreq * (4.0*np.pi * dx*dy*dist**2) / (4.0*np.pi*dist**2 *(180.*60.*60./np.pi)**2 )
   prof_Iin   = prof_Iin   * dxfreq * (4.0*np.pi * dx*dy*dist**2) / (4.0*np.pi*dist**2 *(180.*60.*60./np.pi)**2 )

   #--- scaling factor
   spec1      = np.interp(lam_obs, wavlz, spec)
   spec_scale = np.sum(spec_obs * spec1 / spec_err**2)/np.sum(spec1**2 / spec_err**2)
   spec       = spec * spec_scale
   Jout       = Jout * spec_scale
 
   prof_Iconv1 = np.interp(r_obs, raxis, prof_Iconv)
   wSB         = np.where(SB_err > 0.0)
   # Note that SB should be multiplied by radius.
   SB_scale    = np.sum(SB_obs[wSB]*prof_Iconv1[wSB]*r_obs[wSB]**2 / SB_err[wSB]**2) \
                  / np.sum(prof_Iconv1[wSB]**2*r_obs[wSB]**2/ SB_err[wSB]**2)
   prof_Iconv  = prof_Iconv * SB_scale
   prof_I      = prof_I     * SB_scale
   prof_Iin    = prof_Iin   * SB_scale

   print('   SB_scale              = %6.3e' % (SB_scale))
   print('   spec_scale            = %6.3e' % (spec_scale))
   print('   SB_scale / spec_scale = %6.3e' % (SB_scale / spec_scale))
   fileout.write('   SB_scale              = %6.3e\n' % (SB_scale))
   fileout.write('   spec_scale            = %6.3e\n' % (spec_scale))
   fileout.write('   SB_scale / spec_scale = %6.3e\n' % (SB_scale / spec_scale))
   #----
   w         = (np.where(spec == np.amax(spec)))[0][0]
   wavl_min  = np.interp(np.amax(spec)*1e-4, spec[:w], wavl[:w])
   wavl_max  = np.interp(np.amax(spec)*1e-4, np.flip(spec[w:]), np.flip(wavl[w:]))
   wavlz_min = wavl_min * (1.0 + redshift)
   wavlz_max = wavl_max * (1.0 + redshift)
   # note that xp-sequence should be increasing in numpy.interp(x,xp,fp)

   hdu1.close
   hdu3.close

   ax[0].plot(wavlz, Jout, color='b')
   ax[0].plot(wavlz, spec, linewidth=2, color='r')
   ax[0].set_xlim(wavlz_min, wavlz_max)
   ax[0].set_ylim(0.0, np.amax(spec)*1.2)
   ax[0].xaxis.set_major_locator(MultipleLocator(10.0))
   ax[0].set_xlabel(r'$\lambda_{\rm obs}$ (\AA)')
   ax[0].set_ylabel(r'$10^{-20}$ erg cm$^{-2}$ s$^{-1}$ \AA$^{-1}$')

   ax[0].errorbar(lam_obs, spec_obs, yerr=spec_err, fmt='o', ms=2.0, capsize=1.0, linewidth=1.0, color='k')

   ax2 = ax[0].twiny()
   ax2.set_xlim(wavl_min, wavl_max)
   ax2.xaxis.set_tick_params(labelsize=10)
   ax2.set_xlabel(r'$\lambda_{\rm rest}$ (\AA)')

   #---
   ax[1].plot(raxis, prof_I, color='b')
   ax[1].plot(raxis, prof_Iconv, linewidth=2, color='r')
   ax[1].plot(raxis, prof_Iin, '--', color='g')
   ax[1].set_yscale('log')
   ax[1].set_xlim(0.0, 1.0)
   ymin0 = np.amin(prof_Iconv[:nrmax1])
   ymax0 = np.amax(prof_Iconv[:nrmax1])
   ymin  = ymin0*10.0**(-0.1*np.log10(ymax0/ymin0))
   ymax  = ymax0*10.0**( 0.2*np.log10(ymax0/ymin0))
   ax[1].set_ylim(ymin, ymax)
   ax[1].xaxis.set_major_locator(MultipleLocator(0.2))
   ax[1].set_title('MUSE %s' % museid, fontsize=fontsize)
   ax[1].set_xlabel(r'$r/r_{\rm max}$')
   ax[1].set_ylabel(r'erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2}$')

   ax[1].errorbar(r_obs, SB_obs, yerr=SB_err, fmt='o', ms=2.0, capsize=1.0, linewidth=1.0, color='k')

   #---
   ax[2].plot(raxis, prof_pol0, color='b')
   ax[2].plot(raxis, prof_pol, linewidth=2, color='r')
   ax[2].set_xlim(0.0, 1.0)
   ax[2].set_ylim(0.0, np.amax(prof_pol0)*1.05)
   ax[2].set_title('MUSE %s' % museid, fontsize=fontsize)
   ax[2].set_xlabel(r'$r/r_{\rm max}$')
   ax[2].set_ylabel(r'Polarization ($P_{\rm L}$)')

   #--- make pdf file
   pdf_file = fname+'_profile_all.pdf'
   plt.subplots_adjust(wspace=0.4, hspace=0.2, left=0.1,right=0.9,bottom=0.1,top=0.9)
   s = square_subplots(fig)
   plt.savefig(pdf_file, dpi=1200, bbox_inches='tight', pad_inches=0.02)
   #plt.show()
   
#flist = ['MUSE6905']
flist = ['MUSE1185',
         'MUSE0082',
         'MUSE6905',
         'MUSE1343',
         'MUSE0053',
         'MUSE0171',
         'MUSE0547',
         'MUSE0364']

nf   = len(flist)
npdf = 0

for idx in np.arange(nf):
   fname = flist[idx]
   if os.path.exists(fname+'.fits.gz'):
      plot_profiles(fname)
      npdf = npdf + 1
fileout.close

for idx in np.arange(nf):
   if os.path.exists(flist[idx]+'.fits.gz'):
      if idx == 0:
         in_pdf_str = flist[idx]+'_profile_all.pdf'
      else:
         in_pdf_str = in_pdf_str+' '+flist[idx]+'_profile_all.pdf'

if npdf > 1:
   all_pdf = 'muse_profiles_all.pdf'
   os.system("gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=%s %s" % (all_pdf, in_pdf_str))
   os.system("rm MUSE????_profile_all.pdf")
