#!/usr/bin/env python
import numpy as np
import os
import socket

run_file    = 'run.sh'

def make_runfile(run_file='run.sh', file_in='test.in',hybrid=False):
    import os
    hosts  = 'lart4,lart3,lart2,lart1,mocafe'

    global f2
    try:
       if hybrid == True:
          f2.write('mpirun -hosts $HOSTS -ppn 1 $EXEC %s\n' % file_in)
       else:
          f2.write('mpirun -machinefile $host_file $EXEC %s\n' % file_in)
    except:
       f2 = open(run_file,'w')
       f2.write('#!/bin/bash\n')
       f2.write('exec < /dev/null 2>&1\n')
       f2.write('trap "" HUP\n')
       f2.write('\n')
       f2.write('EXEC=%s\n' % EXEC)
       if hybrid == True:
          f2.write('HOSTS=%s\n' % hosts)
          f2.write('\n')
          f2.write('mpirun -hosts $HOSTS -ppn 1 $EXEC %s\n' % file_in)
       else:
          f2.write('HOSTS=%s\n' % hosts)
          f2.write('\n')
          f2.write('#====== Do not touch starting from here =========\n')
          f2.write('host_file=/tmp/host_file_$RANDOM\n')
          f2.write('\n')
          f2.write('array=$(echo $HOSTS | tr "," "\\n")\n')
          f2.write('for host in $array\n')
          f2.write('do\n')
          f2.write('   if [[ $host = "mocafe" ]]\n')
          f2.write('   then\n')
          f2.write('      num=88\n')
          f2.write('   else\n')
          f2.write('      num=72\n')
          f2.write('   fi\n')
          f2.write('   echo $host:$num >> $host_file\n')
          f2.write('done\n')
          f2.write('\n')
          f2.write('echo ""\n')
          f2.write('echo "Running $EXEC on $HOSTS"\n')
          f2.write('echo "   with the machinefile $host_file"\n')
          f2.write('#====== Do not touch up to here =========\n')
          f2.write('\n')
          f2.write('mpirun -machinefile $host_file $EXEC %s\n' % file_in)

       os.system('chmod +x %s' % run_file)
       #--------------------------


def make_input(line='ly_alpha',spectral_type='voigt',
               tau=1.0,temp=1e4, nx=201,ny=201,nz=201,
               no_photons=None, geometry='sphere', source_geometry='uniform_sphere', DGR_HI=0.0):

    if no_photons == None:
       #if tau <= 1e2:
       #   no_photons = 1e8
       #elif tau == 1e3:
       #   no_photons = 1e7
       if tau <= 1e4:
          no_photons = 1e8
       elif tau <= 1e5:
          no_photons = 1e7
       elif tau <= 1e6:
          no_photons = 1e6
       elif tau <= 1e7:
          no_photons = 1e5
       elif tau == 1e8:
          no_photons = 1e4
       elif tau == 1e9:
          no_photons = 1e3

    num_send_at_once = 100
    if num_send_at_once > no_photons/1000:
       num_send_at_once = np.int32(no_photons/1000)
    if num_send_at_once < 1: num_send_at_once = 1

    line_id    = line
    nxfreq     = 201
    nvelocity  = 201
    atom_num   = 1.0
    abundance  = 1.0
    DGR  = DGR_HI / abundance

    file_in = 't%dtau%d.in' % (np.log10(temp), np.log10(tau))

    xfreq_max = None
    #if tau <= 1e1: xfreq_max = 5.0
    #if tau == 1e2 and temp == 1e1: xfreq_max = 6.0
    #if tau == 1e2 and temp == 1e4: xfreq_max = 5.0
    #if tau == 1e3 and temp == 1e4: xfreq_max = 6.0
    #if xfreq_max != None: xfreq_min = -xfreq_max

    velocity_max = None
    velocity_min = -80.0
    velocity_max = 40.0

    xmax = 1.0
    ymax = 1.0
    zmax = 1.0
    rmax = 1.0
    source_rmax = 1.0

    f = open(file_in, 'w')
    f.write("&parameters\n")
    f.write(" par%%line_id = '%s'\n" % line_id)
    f.write(" par%%num_send_at_once = %d\n" % num_send_at_once)

    f.write("\n")
    f.write(" par%%no_photons   = %.5e\n" % no_photons)
    f.write(" par%%temperature  = %.5e\n" % temp)
    f.write(" par%%taumax       = %.5e\n" % tau)
    f.write(" par%%DGR             = %.4f\n" % DGR)
    if DGR > 0.0:
       f.write(" par%%albedo          = %.4f\n" % albedo)
       f.write(" par%%hgg             = %.4f\n" % hgg)
       f.write(" par%use_reduced_wgt = .true.\n")
       f.write(" par%use_stokes = .false.\n")

    #f.write(" par%xyz_symmetry    = .true.\n")
    f.write(" par%comoving_source = .false.\n")
    f.write(" par%save_all        = .false.\n")
    #f.write(" par%save_direc0     = .true.\n")
    f.write(" par%recoil          = .false.\n")
    f.write(" par%%geometry        = '%s'\n" % geometry)
    f.write(" par%%source_geometry = '%s'\n" % source_geometry)
    f.write(" par%%source_rmax     = %6.4f\n" % source_rmax)
    f.write(" par%%spectral_type   = '%s'\n" % spectral_type)
    f.write(" par%%nx               = %d\n" % nx)
    f.write(" par%%ny               = %d\n" % ny)
    f.write(" par%%nz               = %d\n" % nz)
    f.write(" par%%rmax             = %6.4f\n" % rmax)
    f.write(" par%%xmax             = %6.4f\n" % xmax)
    f.write(" par%%ymax             = %6.4f\n" % ymax)
    f.write(" par%%zmax             = %6.4f\n" % zmax)
    #f.write(" par%%nxfreq    = %d\n"    % nxfreq)
    if xfreq_max != None:
       f.write(" par%%xfreq_min = %6.1f\n" % xfreq_min)
       f.write(" par%%xfreq_max = %6.1f\n" % xfreq_max)
    f.write(" par%%nvelocity    = %d\n"    % nxfreq)
    if velocity_max != None:
       f.write(" par%%velocity_min = %6.1f\n" % velocity_min)
       f.write(" par%%velocity_max = %6.1f\n" % velocity_max)

    f.write(" par%%nprint     = %6.1e\n" % (no_photons/10.0))
    #f.write(" par%out_bitpix = -64\n")
    f.write("/\n")
    f.close()

    #--- run file.
    HOST = socket.gethostname()
    make_runfile(run_file='run.sh',file_in=file_in,hybrid=hybrid)

#--------------------------------
global EXEC
hybrid = False
#--------------------------------
line            = 'HeI_10833'
geometry        = 'sphere'
source_geometry = 'uniform_sphere'
DGR_HI          = 0.0
spectral_type   = 'voigt'

tau_arr  = [1e0, 1e1, 1e2, 1e3]
temp_arr = [1e4]

EXEC='../../LaRT.x'
no_photons = None

tau_arr  = np.array(tau_arr)
temp_arr = np.array(temp_arr)
ntau     = tau_arr.size
ntemp    = temp_arr.size

for j in np.arange(ntau):
   for i in np.arange(ntemp):
      tau  = tau_arr[j]
      temp = temp_arr[i]

      make_input(line=line, spectral_type=spectral_type,
                 tau=tau, temp=temp, geometry=geometry, source_geometry=source_geometry, no_photons=no_photons, DGR_HI=DGR_HI)
