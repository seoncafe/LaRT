#!/usr/bin/env python
import glob, os
import matplotlib.pyplot as plt
import numpy as np
from   astropy.io import fits
from   matplotlib.ticker import AutoMinorLocator

#dir = '../sphere/'
dir = './'
file_list  = glob.glob(dir+"t1tau??*.fits.gz")
file_list1 = glob.glob(dir+"t1tau?.fits.gz")
file_list2  = glob.glob(dir+"t4tau??*.fits.gz")
file_list3 = glob.glob(dir+"t4tau?.fits.gz")
file_list.sort()
file_list1.sort()
file_list2.sort()
file_list3.sort()
file_list.extend(file_list1)
file_list.extend(file_list2)
file_list.extend(file_list3)

f = open("Nscatt_sphere_uniform_voigt.txt","w")
for file in file_list:
   hdu     = fits.open(file)
   hdr     = hdu[1].header
   tau0    = hdr['TAUMAX']
   temp    = hdr['TEMP']
   try:
      nscatt  = hdr['NSC_GAS']
   except:
      nscatt  = 0
   print("sphere %5.0e %5.0e %9.6e" % (temp, tau0, nscatt))
   f.write("%5.0e %5.0e %9.6e\n" % (temp, tau0, nscatt))
f.close()

#dir = '../slab/'
#file_list = glob.glob(dir+"t?tau?.fits.gz")
#file_list.sort()
#
#f = open("Nscatt_slab.txt","w")
#for file in file_list:
#   hdu     = fits.open(file)
#   hdr     = hdu[1].header
#   tau0    = hdr['TAUMAX']
#   temp    = hdr['TEMP']
#   try:
#      nscatt  = hdr['NSC_HI']
#   except:
#      nscatt  = 0
#   print("slab   %5.0e %5.0e %9.6e" % (temp, tau0, nscatt))
#   f.write("%5.0e %5.0e %9.6e\n" % (temp, tau0, nscatt))
#f.close()
#
#dir = '../slab1/'
#file_list = glob.glob(dir+"t?tau?.fits.gz")
#file_list.sort()

#f = open("Nscatt_slab1.txt","w")
#for file in file_list:
#   hdu     = fits.open(file)
#   hdr     = hdu[1].header
#   tau0    = hdr['TAUMAX']
#   temp    = hdr['TEMP']
#   try:
#      nscatt  = hdr['NSC_HI']
#   except:
#      nscatt  = 0
#   print("slab1  %5.0e %5.0e %9.6e" % (temp, tau0, nscatt))
#   #f.write("%5.0e %5.0e %9.6e\n" % (temp, tau0, nscatt))
##f.close()
