#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
from   lart_io import load_lart, glob_lart

flist = glob_lart('*', suffix='_obs')
flist.sort()
nf    = len(flist)

for idx in np.arange(nf):
   # peeled-off intensity image
   fname = flist[idx]
   obs   = load_lart(fname)
   scattered = obs.section('Scattered')
   direct    = obs.section('Direct')
   direct0   = obs.section('Direct0')

   dxim          = float(scattered.attr('CD1_1', 1.0))
   dyim          = float(scattered.attr('CD2_2', 1.0))
   distance_unit = float(scattered.attr('DIST_CM', 1.0))
   distance_in   = float(scattered.attr('DISTANCE', 0.0))
   obsx          = float(scattered.attr('OBSX', 0.0))
   obsy          = float(scattered.attr('OBSY', 0.0))
   obsz          = float(scattered.attr('OBSZ', 0.0))
   distance_star_to_planet = float(scattered.attr('DIST_S2P', 0.0))

   print('\n>>> ', fname)

   intensity_unit = int(scattered.attr('I_UNIT', 0))
   if intensity_unit == 0:
      bin_unit = float(scattered.attr('DXFREQ', 1.0))
   else:
      bin_unit = float(scattered.attr('DWAVE', 1.0))

   flux_factor = float(scattered.attr('FLUXFAC', 1.0))
   distance = np.sqrt(obsx**2 + obsy**2 + (obsz + distance_star_to_planet)**2)

   im1 = scattered.data
   im2 = direct.data
   im0 = direct0.data if direct0 is not None else np.zeros_like(im1)
   im  = im1 + im2

   #-- the output intensity is in unit of per dimensionless frequency [x = (nu - nu0)/(nu0*(vthermal/c))] if intensity_unit = 0.
   #-- the ouptut is in unit of per wavelength (angstrom) if intensity_unit = 1.
   #-- scale factor to convert the intensity (counts/cm^2/s/sr^2/angstrom) to total luminosity.
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
