#!/usr/bin/env python
import numpy as np

def make_data(line='MgII',fname='fiducial',nH0=0.1,abund=10.0**(-5.47),rinner=1.0,router=20.0,nx=300,ny=300,nz=300,T=1e4,
              Vexp=1000.0,save_velo_file=False):
   from astropy.io import fits
   #----------------------
   # generate density, temperature (and velocity) data file.
   # (2021.07.03)
   #----------------------
   fname = line[0:4]+'_'+fname

   nMg0 = abund * nH0
   z, y, x = np.indices((nz, ny, nx))

   x = x + 0.5
   y = y + 0.5
   z = z + 0.5

   xcen = nx/2.0
   ycen = ny/2.0
   zcen = nz/2.0
   x    = (x - xcen)/(nx/2.0) * router
   y    = (y - xcen)/(ny/2.0) * router
   z    = (z - xcen)/(nz/2.0) * router

   r    = np.sqrt(x**2 + y**2 + z**2)

   # assume that the medium is fully ionized in r < rinner.
   dens        = np.zeros((nz, ny, nx), dtype='float32')
   w           = np.where((r >= rinner) & (r <= router))
   dens[w]     = nMg0 * (rinner/r[w])**2

   hdu_dens = fits.PrimaryHDU(dens)
   hdu1     = fits.HDUList([hdu_dens])
   hdu1.writeto(fname+'_dens.fits.gz', overwrite=True)

   if (save_velo_file == True) & (Vexp > 0.0):
      vfx  = np.zeros((nz, ny, nz), dtype='float32')
      vfy  = np.zeros((nz, ny, nz), dtype='float32')
      vfz  = np.zeros((nz, ny, nz), dtype='float32')
      velo = np.zeros((nz, ny, nz, 3), dtype='float32')
      w    = np.where((r >= rinner) & (r <= router))
      vfx[w] = x[w] * Vexp / router
      vfy[w] = y[w] * Vexp / router
      vfz[w] = z[w] * Vexp / router
      velo[:,:,:,0] = vfx
      velo[:,:,:,1] = vfy
      velo[:,:,:,2] = vfz
      hdu_velo = fits.PrimaryHDU(velo)
      hdu2     = fits.HDUList([hdu_velo])
      hdu2.writeto(fname+'_velo.fits.gz', overwrite=True)

def make_input(line='MgII',fname='fiducial',rinner=1.0,router=20.0,nx=300,ny=300,nz=300,
               bturb=15.0,Vexp=1000.0,save_velo_file=False):
    if line == 'MgII':
       lambda_min = 2784.0
       lambda_max = 2815.0
       nlambda    = 620
       line_id    = 'MgII_2796'
    if line == 'FeII_UV1':
       lambda_min = 2578.0
       lambda_max = 2634.0
       #nlambda    = 800
       nlambda    = 1120
       line_id    = line
    if line == 'FeII_UV2':
       lambda_min = 2373.0
       lambda_max = 2398.0
       nlambda    = 500
       line_id    = line
    if line == 'FeII_UV3':
       lambda_min = 2340.0
       lambda_max = 2385.0
       nlambda    = 900
       line_id    = line

    fname1 = line+'_'+fname
    no_photons = 1e8
    #bturb      = 15.0  # km/s
    temp       = (bturb/(0.12843374/np.sqrt(24.0)))**2
    file_in    = fname1+'.in'
    dens_file  = line[0:4]+'_'+fname+'_dens.fits.gz'
    temp_file  = line[0:4]+'_'+fname+'_temp.fits.gz'
    velo_file  = line[0:4]+'_'+fname+'_velo.fits.gz'

    source_rmax = 0.5
    xmax = router
    ymax = router
    zmax = router
    rmin = rinner
    rmax = router
    nxim = 129
    nyim = 129

    f = open(file_in, 'w')
    f.write("&parameters\n")
    f.write(" par%%line_id = '%s'\n" % line_id)
    a = int(np.log10(no_photons))
    b = no_photons/np.power(10,a)
    f.write(" par%%no_photons   = %0.3fe%d\n" % (b,a))
    a = int(np.log10(temp))
    b = temp/np.power(10,a)
    f.write(" par%%temperature  = %0.3fe%d\n" % (b,a))
    f.write(" par%DGR          = 0.0\n")
    if (Vexp > 0.0) & (save_velo_file == False):
       f.write(" par%velocity_type = 'hubble'\n")
       f.write(" par%%Vexp        = %0.1f\n" % (Vexp))
    f.write(" par%comoving_source = .false.\n")
    f.write(" par%recoil          = .true.\n")
    f.write(" par%source_geometry = 'sphere'\n")
    f.write(" par%%source_rmax     = %3.1f\n" % source_rmax)
    f.write(" par%spectral_type   = 'continuum'\n")
    f.write(" par%%dens_file       = '%s'\n" % dens_file)
    if (Vexp > 0.0) & (save_velo_file == True):
       f.write(" par%%velo_file       = '%s'\n" % velo_file)
    f.write(" par%%nx               = %d\n" % nx)
    f.write(" par%%ny               = %d\n" % ny)
    f.write(" par%%nz               = %d\n" % nz)
    f.write(" par%%xmax             = %3.1f\n" % xmax)
    f.write(" par%%ymax             = %3.1f\n" % ymax)
    f.write(" par%%zmax             = %3.1f\n" % zmax)
    f.write(" par%%rmin             = %3.1f\n" % rmin)
    f.write(" par%%rmax             = %3.1f\n" % rmax)
    f.write(" par%distance_unit    = 'kpc'\n")
    f.write(" par%%nlambda  = %d\n"    %nlambda)
    f.write(" par%%lambda_min = %5.1f\n" % lambda_min)
    f.write(" par%%lambda_max = %5.1f\n" % lambda_max)
    f.write(" par%%nxim = %d\n" % nxim)
    f.write(" par%%nyim = %d\n" % nyim)
    f.write(" par%out_bitpix = -64\n")
    f.write(" par%use_stokes = .false.\n")
    f.write(" par%nprint   = 1e7\n")
    f.write("/\n")
    f.close()

make_data(line='MgII', fname='a',abund=10.0**(-5.47),Vexp=1000.0,save_velo_file=False)
make_data(line='MgII', fname='b',abund=10.0**(-5.47)/2.0,Vexp=1000.0,save_velo_file=False)

make_input(line='MgII',fname='a',Vexp=1000.0,save_velo_file=False)
make_input(line='MgII',fname='b',Vexp=1000.0,save_velo_file=False)

#-------------------
make_data(line='FeII_UV1', fname='a',abund=10.0**(-5.47),Vexp=1000.0,save_velo_file=False)
make_data(line='FeII_UV1', fname='b',abund=10.0**(-5.47)/2.0,Vexp=1000.0,save_velo_file=False)
make_data(line='FeII_UV1', fname='c',abund=10.0**(-5.47)/4.0,Vexp=1000.0,save_velo_file=False)

make_input(line='FeII_UV1',fname='a',Vexp=1000.0,save_velo_file=False)
make_input(line='FeII_UV1',fname='b',Vexp=1000.0,save_velo_file=False)
make_input(line='FeII_UV1',fname='c',Vexp=1000.0,save_velo_file=False)

make_input(line='FeII_UV2',fname='a',Vexp=1000.0,save_velo_file=False)
make_input(line='FeII_UV2',fname='b',Vexp=1000.0,save_velo_file=False)
make_input(line='FeII_UV2',fname='c',Vexp=1000.0,save_velo_file=False)

make_input(line='FeII_UV3',fname='a',Vexp=1000.0,save_velo_file=False)
make_input(line='FeII_UV3',fname='b',Vexp=1000.0,save_velo_file=False)
make_input(line='FeII_UV3',fname='c',Vexp=1000.0,save_velo_file=False)

#-------------
#make_data(fname='fiducial2',abund=10.0**(-5.47),Vexp=1000.0,router=10.0,save_velo_file=False)
#make_input(fname='fiducial2',Vexp=1000.0,router=10.0,save_velo_file=False)
