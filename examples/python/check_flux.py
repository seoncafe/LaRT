#!/usr/bin/env python
from   astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import glob

flist = glob.glob('*_obs*.fits.gz')
flist.sort()
nf    = len(flist)

for idx in np.arange(nf):
   # peeled-off intensity image
   fname = flist[idx]
   hdu   = fits.open(fname)
   hdr   = hdu[0].header
   dxim  = hdr['CD1_1']
   dyim  = hdr['CD2_2']
   distance_unit = hdr['DIST_CM']
   distance_in   = hdr['DISTANCE']
   obsx          = hdr['OBSX']
   obsy          = hdr['OBSY']
   obsz          = hdr['OBSZ']
   try:
      distance_star_to_planet = hdr['DIST_S2P']
   except:
      distance_star_to_planet = 0.0

   print('\n>>> ', fname)

   try:
      intensity_unit = hdr['I_UNIT']
   except:
      intensity_unit = 0
   if intensity_unit == 0: bin_unit = hdr['DXFREQ']
   if intensity_unit == 1: bin_unit = hdr['DWAVE']

   try:
      # Now (LaRT_v1.33a), distance is defined to be the distance from the rotation center to the observer.
      flux_factor = hdr['FLUXFAC']
   except:
      flux_factor = 1.0
   distance = np.sqrt(obsx**2 + obsy**2 + (obsz + distance_star_to_planet)**2)

   im1 = hdu[0].data
   im2 = hdu[1].data
   im0 = hdu[2].data
   im  = im1 + im2
   hdu.close()


   #-- the output intensity is in unit of per dimensionless frequency [x = (nu - nu0)/(nu0*(vthermal/c))] if intensity_unit = 0.
   #-- the ouptut is in unit of per wavelength (angstrom) if intensity_unit = 1.
   #-- scale factor to convert the intensity (counts/cm^2/s/sr^2/angstrom) to total luminosity.
   #scale = 4.0*np.pi * (distance * np.tan(dxim * np.pi/180.0)) * (distance * np.tan(dyim * np.pi/180.0)) * (distance_unit**2) * bin_unit
   scale = 4.0*np.pi * (distance * dxim * np.pi/180.0) * (distance * dyim * np.pi/180.0) * (distance_unit**2) * bin_unit
   ftot  = im.sum()  * scale
   fsca  = im1.sum() * scale
   fdir  = im2.sum() * scale
   fdir0 = im0.sum() * scale
   ftot1 =  fsca * flux_factor + fdir

   print("  distance unit in cm                      : %15.7e" % distance_unit)
   print("  distance given in the model input        : %15.7e" % distance_in)
   print("  distance between the star and planet     : %15.7e" % distance_star_to_planet)
   if distance_star_to_planet == 0.0:
      print("  distance between the planet and observer : %15.7e" % distance)
   else:
      print("  distance between the star   and observer : %15.7e" % distance)
      print("  FLUX_FACTOR                              : %15.7e" % flux_factor)
   print("  observer coordinates (X, Y, Z)           : %15.7e, %15.7e, %15.7e" % (obsx, obsy, obsz))

   print('     TOTAL : %10.7f (%10.7f)' % (ftot1, ftot))
   print('     SCATT : %10.7f (%10.7f)' % (fsca, fsca*flux_factor))
   print('     DIREC : %10.7f' % fdir)
   print('     DIREC0: %10.7f' % fdir0)
