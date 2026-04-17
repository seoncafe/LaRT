#!/usr/bin/env python
import numpy as np

#-- a toy model
#-- density, velocity, temperature vs radius
nr   = 91
rmin = 1.0
rmax = 10.0
dr   = (rmax-rmin)/(nr-1)
r    = np.arange(nr)*dr + rmin

rscale = 2.0
#dens   = 5e8   * np.exp(-(r-rmin)/rscale)
dens   = 1e8   * np.exp(-(r-rmin)/rscale)
velo   = 140.0 * (r-rmin) / (rmax-rmin)
temp   = 8e3   * (r-rmin)/(rmax-rmin) * np.exp(-(r-rmin)/(rmax-rmin)) + 8e3

f_dens = open('dens_profile.txt','w')
f_velo = open('velo_profile.txt','w')
f_temp = open('temp_profile.txt','w')
for i in np.arange(nr):
   #print("%4.1f  %0.5e  %0.5e  %0.5e" % (r[i], dens[i], velo[i], temp[i]))
   f_dens.write("%4.1f  %0.5e\n" % (r[i], dens[i]))
   f_velo.write("%4.1f  %0.5e\n" % (r[i], velo[i]))
   f_temp.write("%4.1f  %0.5e\n" % (r[i], temp[i]))
f_dens.close()
f_velo.close()
f_temp.close()

#-- line profile
nwav = 101
wav1 = 1214.0
wav2 = 1217.4
dwav = (wav2-wav1)/(nwav-1)
wav  = np.arange(nwav)*dwav + wav1
wav0 = 1215.668
sig  = 0.30
line = np.exp(-(wav-wav0)**2/(2.0*sig**2))

f_line = open('line_profile.txt','w')
for i in np.arange(nwav):
   f_line.write("%8.3f  %0.5e\n" % (wav[i], line[i]))
f_line.close()
