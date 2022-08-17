#!/usr/bin/env python
import os
import socket
import glob
import numpy as np

no_photons = 1e6
no_print   = 1e5
use_reduced_wgt     = '.true.'
use_stokes          = '.true.'
save_sightline_tau  = '.true.'
save_peeloff_2D     = '.true.'
comoving_source = '.false.'
out_merge       = '.true.'
save_Jin        = '.true.'
save_Jabs       = '.true.'
source_geometry = 'ssh'
velocity_type   = 'ssh'

temp           = 1e4
scatt_mat_file =  '../../data/mueller_Lyalpha.dat'
spectral_type  = 'voigt'

xmax = 1.0
ymax = 1.0
zmax = 1.0
rmax = 1.0
nx   = 201
ny   = 201
nz   = 201
nxim = 129
nyim = 129
nxfreq    = 401
xfreq_min = -150.0
xfreq_max =   50.0
distance  = rmax*1e3

DGR1 = 0.0

#-----------
fname = './Leclercq/muse_bestfit.txt'
id, rsUV, rsHI, rpeak, Vpeak, DeltaV, log10tau, DGR, z, r_lim,lambda1,lambda2 = np.loadtxt(fname, skiprows=3, unpack=True)
nfile = id.size
#----------- run file to use the MPI version.
file_sh  = 'run.sh'
fsh      = open(file_sh,'w')
fsh.write('#!/bin/bash\n')
fsh.write('exec < /dev/null 2>&1\n')
fsh.write('trap "" HUP\n')
fsh.write("\n")
fsh.write("EXEC=../../LaRT_FS.x\n")
fsh.write("HOST=all_hosts\n")
fsh.write("\n")
#----------- run file to use the MPI+Openmp hybrid version.
file_sh2 = 'run_hybrid.sh'
fsh2     = open(file_sh2,'w')
fsh2.write('#!/bin/bash\n')
fsh2.write('exec < /dev/null 2>&1\n')
fsh2.write('trap "" HUP\n')
fsh2.write("\n")
fsh2.write("EXEC=../../LaRT_FS.x\n")
fsh2.write("HOSTS=mocafe,lart1,lart2,lart3\n")
fsh2.write("\n")
#-----------
for idx in np.arange(nfile):
   file_in  = "MUSE%04d" % (id[idx]) + ".in"
   taumax   = 10.0**log10tau[idx]

   f = open(file_in,'w')
   f.write("&parameters\n")
   a = int(np.log10(no_photons))
   b = no_photons/10.0**a
   f.write(" par%%no_photons      = %0.1fe%d\n" % (b,a))
   f.write(" par%%scatt_mat_file  = '%s'\n" %scatt_mat_file)
   f.write(" par%%use_stokes      = %s\n" %use_stokes)
   f.write(" par%%use_reduced_wgt = %s\n" %use_reduced_wgt)
   f.write(" par%%comoving_source = %s\n" %comoving_source)
   f.write(" par%%save_sightline_tau   = %s\n" %save_sightline_tau)
   f.write(" par%%save_peeloff_2D      = %s\n" %save_peeloff_2D)
   f.write(" par%%save_Jin        = %s\n" %save_Jin)
   f.write(" par%%save_Jabs       = %s\n" %save_Jabs)
   f.write(" par%%temperature     = %6.1e\n" %temp)
   f.write(" par%%spectral_type   = '%s'\n" %spectral_type)
   f.write(" par%%source_geometry = '%s'\n" %source_geometry)
   f.write(" par%%velocity_type   = '%s'\n" %velocity_type)
   f.write(" par%%Vpeak           = %6.1f\n" %Vpeak[idx])
   f.write(" par%%rpeak           = %6.3f\n" %rpeak[idx])
   f.write(" par%%DeltaV          = %6.1f\n" %DeltaV[idx])
   f.write(" par%%density_rscale  = %6.4f\n" %rsHI[idx])
   f.write(" par%%source_rscale   = %6.4f\n" %rsUV[idx])
   f.write(" par%%DGR             = %5.3f\n" %DGR[idx])
   f.write(" par%%taumax          = %9.3e\n" % (taumax))
   f.write("\n")
   f.write(" par%distance_unit   = ''\n")
   f.write(" par%%xmax    = %5.1f\n" %xmax)
   f.write(" par%%ymax    = %5.1f\n" %ymax)
   f.write(" par%%zmax    = %5.1f\n" %zmax)
   f.write(" par%%rmax    = %5.1f\n" %rmax)
   f.write(" par%%nx = %d\n" %nx)
   f.write(" par%%ny = %d\n" %ny)
   f.write(" par%%nz = %d\n" %nz)
   f.write("\n")
   f.write(" par%%nxim = %d\n" %nxim)
   f.write(" par%%nyim = %d\n" %nyim)
   f.write(" par%%nxfreq    = %d\n" %nxfreq)
   f.write(" par%%xfreq_min = %d\n" %xfreq_min)
   f.write(" par%%xfreq_max = %d\n" %xfreq_max)

   #f.write(" par%obsx      = 0.0  0.0  0.0  0.0  1.0 -1.0\n")
   #f.write(" par%obsy      = 0.0  0.0  1.0 -1.0  0.0  0.0\n")
   #f.write(" par%obsz      = 1.0 -1.0  0.0  0.0  0.0  0.0\n")

   #f.write(" par%obsx      = 1.0   1.0   1.0  -1.0  -1.0  -1.0   1.0  -1.0\n")
   #f.write(" par%obsy      = 1.0   1.0  -1.0   1.0  -1.0   1.0  -1.0  -1.0\n")
   #f.write(" par%obsz      = 1.0  -1.0   1.0   1.0   1.0  -1.0  -1.0  -1.0\n")

   #f.write(" par%obsx      = 0.0  0.0  0.0  0.0  1.0 -1.0  1.0   1.0   1.0  -1.0  -1.0  -1.0   1.0  -1.0\n")
   #f.write(" par%obsy      = 0.0  0.0  1.0 -1.0  0.0  0.0  1.0   1.0  -1.0   1.0  -1.0   1.0  -1.0  -1.0\n")
   #f.write(" par%obsz      = 1.0 -1.0  0.0  0.0  0.0  0.0  1.0  -1.0   1.0   1.0   1.0  -1.0  -1.0  -1.0\n")

   f.write(" par%%distance  = %0.1e\n" % distance)
   f.write("\n")
   f.write(" par%use_master_slave  = .true.\n")
   a = int(np.log10(no_print))
   b = no_print/10.0**a
   f.write(" par%%no_print  = %6.1e\n" % (no_print))
   f.write(" par%%out_merge = %s\n" %out_merge)
   f.write("/\n")
   f.close()

   fsh.write("mpirun -machinefile $HOST $EXEC %s\n" % (file_in))
   fsh2.write("mpirun -hosts $HOSTS -ppn 1 $EXEC %s\n" % (file_in))
fsh.close()
fsh2.close()
os.system("chmod +x "+file_sh)
os.system("chmod +x "+file_sh2)
