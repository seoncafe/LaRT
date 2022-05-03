#!/usr/bin/env python
import numpy as np

def make_data(fname='clumpy',fraction_hot=0.4,nx=33,temperature=1e4,xyz_symmetry=False,Vexp=0.0):
   from astropy.io import fits
   #----------------------
   # generate density, temperature (and velocity) data file.
   # (2021.06.26)
   #----------------------
 
   ny = nx
   nz = nx
   z, y, x = np.indices((nz, ny, nx))

   if xyz_symmetry == False or (nx//2)*2 == nx: x = x + 0.5
   if xyz_symmetry == False or (ny//2)*2 == ny: y = y + 0.5
   if xyz_symmetry == False or (nz//2)*2 == nz: z = z + 0.5

   if xyz_symmetry == True:
      xcen = 0.0
      ycen = 0.0
      zcen = 0.0
   else:
      xcen = nx/2.0
      ycen = ny/2.0
      zcen = nz/2.0
      x    = x - xcen
      y    = y - xcen
      z    = z - xcen

   rmax = nx/2.0
   x    = x / rmax
   y    = y / rmax
   z    = z / rmax
   r    = np.sqrt(x**2 + y**2 + z**2)

   #--- assuming the pressure equilibrium.
   dens  = np.zeros((nz, ny, nx), dtype='float32') + 1.0
   temp  = np.zeros((nz, ny, nx), dtype='float32') + temperature
   w_hot = np.where(np.random.random((nz, ny, nz)) < fraction_hot)
   dens[w_hot] = 0.0
   w           = np.where(r > 1.0)
   dens[w]     = 0.0
   temp[w]     = 0.0

   hdu_dens = fits.PrimaryHDU(dens)
   hdu1     = fits.HDUList([hdu_dens])
   hdu1.writeto(fname+'_dens.fits.gz', overwrite=True)
   hdu_temp = fits.PrimaryHDU(temp)
   hdu2     = fits.HDUList([hdu_temp])
   hdu2.writeto(fname+'_temp.fits.gz', overwrite=True)

   if np.abs(Vexp) > 0.0:
      vfx  = np.zeros((nz, ny, nz), dtype='float32')
      vfy  = np.zeros((nz, ny, nz), dtype='float32')
      vfz  = np.zeros((nz, ny, nz), dtype='float32')
      velo = np.zeros((nz, ny, nz, 3), dtype='float32')
      w    = np.where((r > 0.0) & (r < 1.0))
      vfx[w] = x[w]/r[w] * Vexp
      vfy[w] = y[w]/r[w] * Vexp
      vfz[w] = z[w]/r[w] * Vexp
      velo[:,:,:,0] = vfx
      velo[:,:,:,1] = vfy
      velo[:,:,:,2] = vfz
      hdu_velo = fits.PrimaryHDU(velo)
      hdu2     = fits.HDUList([hdu_velo])
      hdu2.writeto(fname+'_velo.fits.gz', overwrite=True)

def make_input(fname='clumpy',NHI=1e20,source_rmax=0.0,nx=33,temperature=1e4,Vexp=0.0,xyz_symmetry=False,
               lambda_min=1214.0, lambda_max=1218.0,nlambda=200, no_photons=1e5,
               nxim=129,nyim=129,
               obsx=0.0, obsy=0.0, obsz=1.0):
    temp       = 1e4
    temp0      = 1e4
    file_in    = fname+'.in'
    dens_file  = fname+'_dens.fits.gz'
    temp_file  = fname+'_temp.fits.gz'
    velo_file  = fname+'_velo.fits.gz'
    rmax       = 1.0

    ny = nx
    nz = nx

    f = open(file_in, 'w')
    f.write("&parameters\n")
    a = int(np.log10(no_photons))
    b = no_photons/np.power(10,a)
    f.write(" par%%no_photons   = %0.1fe%d\n" % (b,a))
    a = int(np.log10(temp))
    b = temp/np.power(10,a)
    f.write(" par%%temperature  = %0.1fe%d\n" % (b,a))
    a = int(np.log10(temp0))
    b = temp0/np.power(10,a)
    f.write(" par%%temperature0 = %0.1fe%d\n" % (b,a))
    a = int(np.log10(NHI))
    b = NHI/np.power(10.0,a)
    f.write(" par%%N_HIhomo    = %0.1fe%d\n" % (b,a))
    f.write(" par%DGR         = 0.0\n")
    f.write(" par%%Vexp        = %0.1f\n" % (Vexp))
    if xyz_symmetry == True: f.write(" par%xyz_symmetry    = .true.\n")
    #f.write(" par%comoving_source = .true.\n")
    f.write(" par%comoving_source = .false.\n")
    f.write(" par%recoil          = .true.\n")
    f.write(" par%source_geometry = 'sphere'\n")
    f.write(" par%spectral_type   = 'voigt0'\n")
    f.write(" par%%source_rmax     = %3.1f\n" % source_rmax)
    f.write(" par%%dens_file       = %s\n" % dens_file)
    f.write(" par%%temp_file       = %s\n" % temp_file)
    if np.abs(Vexp) > 0.0:
        f.write(" par%%velo_file       = %s\n" % velo_file)
    f.write(" par%%nx = %d\n" %nx)
    f.write(" par%%ny = %d\n" %ny)
    f.write(" par%%nz = %d\n" %nz)
    f.write(" par%%rmax    = %5.1f\n" %rmax)
    f.write(" par%%lambda_min  = %8.2f\n"  %lambda_min)
    f.write(" par%%lambda_max  = %8.2f\n"  %lambda_max)
    f.write(" par%%nlambda  = %d\n"    %nlambda)
    f.write(" par%%nxim     = %d\n" %nxim)
    f.write(" par%%nyim     = %d\n" %nyim)
    f.write(" par%%obsx     = %5.2f\n" %obsx)
    f.write(" par%%obsy     = %5.2f\n" %obsy)
    f.write(" par%%obsz     = %5.2f\n" %obsz)
    f.write(" par%nprint   = 1000000\n")
    f.write("/\n")
    f.close()

Vexp = 20.0
make_data(Vexp=Vexp)
make_input(Vexp=Vexp, obsx=1.0, obsy=1.0, obsz=1.0)
