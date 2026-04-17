#!/usr/bin/env python
import numpy as np
import os
import socket

run_file    = 'run.sh'
run_polaris = 'run_polaris.sh'

def float_format(real_num):
    if real_num >= 10.0:
       a = int(np.log10(real_num))
       b = real_num/np.power(10,a)
    else:
       a = 0
       b = real_num
    return "%1de%+1d" % (b,a)
    #return "%0.3fe%02d" % (b,a)
    #return "%6.1e" % real_num
    #return "%4.0e" % real_num

def make_runfile(run_file='run.sh', file_in='test.in',hybrid=True):
    import os
    if hybrid == True:  hosts  = 'lart4'
    #if hybrid == True:  hosts  = 'mocafe,lart1,lart2,lart3,lart4'
    #if hybrid == True:  hosts  = 'lart1,lart2,lart3,lart4'
    if hybrid == False: hosts = 'all_hosts'
    global f2
    try:
       if hybrid == True:
          f2.write('mpirun -hosts $HOSTS -ppn 1 $EXEC %s\n' % file_in)
       else:
          f2.write('mpirun -machinefile $HOST $EXEC %s\n' % file_in)
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
          f2.write('HOST=%s\n' % hosts)
          f2.write('\n')
          f2.write('mpirun -machinefile $HOST $EXEC %s\n' % file_in)
       os.system('chmod +x %s' % run_file)
       #--------------------------
       if hybrid == False:
          f1 = open(hosts, "w")
          f1.write('mocafe:88\n')
          f1.write('lart1:72\n')
          f1.write('lart2:72\n')
          f1.write('lart3:72\n')
          f1.write('lart4:72\n')
          f1.close()
       #--------------------------

def make_runfile_polaris(jobname='MgII',run_file='run_polaris.sh',file_in='test.in',outfile='out_polaris.txt',hybrid=True):
    import os
    number_of_nodes = 21
    ntasks_per_node = 16
    global f3
    try:
       if hybrid == True:
          f3.write('mpirun -ppn 1 $EXEC %s\n' % file_in)
       else:
          f3.write('mpirun $EXEC %s\n' % file_in)
    except:
       f3 = open(run_file,'w')
       f3.write('#!/bin/bash\n')
       f3.write('#SBATCH --job-name="%s"    # Job name\n' % jobname)
       f3.write('#SBATCH --output="%s"        # Name of stdout output file\n' % outfile)
       f3.write('#SBATCH --open-mode=append    # open the output and error files using append mode (default is truncate mode)\n')
       f3.write('#SBATCH --time=7-00:00:00     # time limit, 1 hour, for this job.\n')
       f3.write('#SBATCH --partition=ibm       # select sm partition (or queue), which is composed of supermicro servers\n')
       f3.write('#SBATCH --nodes=%d            # requested number of nodes\n' % number_of_nodes)
       f3.write('#SBATCH --ntasks-per-node=%d   # requested number of mpi tasks per node\n' % ntasks_per_node)
       f3.write('module purge                  # remove all the modules\n')
       f3.write('export I_MPI_FABRICS=shm:ofi  # use shared memory and infiniband as interconnects of the MPI job\n')
       f3.write('unset I_MPI_SHM_LMT\n')
       f3.write('cd $SLURM_SUBMIT_DIR          # change your working directory where you launch your job.\n')
       f3.write('\n')
       f3.write('EXEC=%s\n' % EXEC)
       f3.write('\n')
       if hybrid == True:
          f3.write('mpirun -ppn 1 $EXEC %s\n' % file_in)
       else:
          f3.write('mpirun $EXEC %s\n' % file_in)
       os.system('chmod +x %s' % run_file)

def make_input(line='CIV',use_stokes=True, spectral_type='continuum',
               tau=1.0,bturb=15.0, nx=200,ny=200,nz=1, Vexp=0.0,
               no_photons=1e7, geometry='sphere', source_geometry='uniform', DGR_HI=0.0):
    file_in = 'tau%s_V%03d.in' % (float_format(tau), Vexp)

    num_send_at_once = 1000

    # abundances are taken from Table 9.5 in Draine's book (for zeta Oph).
    if line == 'MgII':
       line_id    = 'MgII_2796'
       nlambda    = 800
       atom_num   = 12.0
       lambda_min = 2790.0
       lambda_max = 2810.0
       # abundance in WNM (F_star = 0.1)
       abundance  = 1.78e-5
       albedo     = 0.57146498
       hgg        = 0.54936271
       mat_file   = 'mueller_2796.dat'
    elif line == 'CIV' :
       line_id    = 'CIV_1548'
       nlambda    = 400
       atom_num   = 6.0
       lambda_min = 1546.0
       lambda_max = 1554.0
       # abundance in WIM (F_star = -0.1)
       abundance  = 1.14e-4
       albedo     = 0.38960096
       hgg        = 0.65915936
       mat_file   = 'mueller_1548.dat'
    elif line == 'SiII_1193':
       line_id    = 'SiII_1193'
       nlambda    = 240
       atom_num   = 14.0
       lambda_min = 1188.0
       lambda_max = 1200.0
       # abundance in WNM (F_star = 0.1)
       abundance  = 1.87e-5
       albedo     = 0.320
       hgg        = 0.674
       mat_file   = 'mueller_1193.dat'
   # abundance of Fe (WNM F_star = 0.1)
   # abundance = 2.9e-6
   # abundance = 5.92e-4 (Oxygen, WNM F_star = 0.1)

    temp = (bturb/(0.12843374/np.sqrt(2.0*atom_num)))**2
    DGR  = DGR_HI / abundance

    xmax = 1.0
    ymax = 1.0
    zmax = 1.0
    rmax = 1.0
    source_rmax = 1.0

    nxim = 100
    nyim = 100
    distance = 1e3
    dxim     = np.arctan(xmax/(distance-xmax))/(nxim/2.0) * (180.0/np.pi)
    dyim     = np.arctan(ymax/(distance-ymax))/(nyim/2.0) * (180.0/np.pi)

    f = open(file_in, 'w')
    f.write("&parameters\n")
    f.write(" par%%line_id = '%s'\n" % line_id)
    f.write(" par%%num_send_at_once = %d\n" % num_send_at_once)

    f.write("\n")
    f.write(" par%%no_photons   = %.5e\n" % no_photons)
    f.write(" par%%temperature  = %.5e\n" % temp)
    f.write(" par%%taumax       = %.5e\n" % tau)
    f.write(" par%velocity_type = 'hubble'\n")
    f.write(" par%%Vexp         = %.1f\n" % Vexp)

    f.write(" par%%DGR             = %.4f\n" % DGR)
    if DGR > 0.0:
       f.write(" par%%albedo          = %.4f\n" % albedo)
       f.write(" par%%hgg             = %.4f\n" % hgg)
       f.write(" par%use_reduced_wgt = .true.\n")
       if use_stokes == True:
          f.write(" par%%scatt_mat_file  = '../../LaRT_v1.36_hybrid/data/%s'\n" % mat_file)
    if use_stokes == True:
       f.write(" par%use_stokes = .true.\n")
    else:
       f.write(" par%use_stokes = .false.\n")

    #f.write(" par%comoving_source = .false.\n")
    f.write(" par%save_direc0     = .true.\n")
    f.write(" par%recoil          = .true.\n")
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
    f.write(" par%distance_unit    = ''\n")
    f.write(" par%%nlambda    = %d\n"    % nlambda)
    f.write(" par%%lambda_min = %6.1f\n" % lambda_min)
    f.write(" par%%lambda_max = %6.1f\n" % lambda_max)

    f.write("\n")
    f.write(" par%%distance = %6.1e\n" % distance)
    f.write(" par%%nxim     = %d\n" % nxim)
    f.write(" par%%nyim     = %d\n" % nyim)
    f.write(" par%%dxim     = %0.14e\n" % dxim)
    f.write(" par%%dyim     = %0.14e\n" % dyim)
    #f.write(" par%beta     = 0.0, 30.0, 45.0, 60.0, 70.0, 75.0, 80.0, 85.0, 88.0, 90.0\n")
    #f.write(" par%alpha    = 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0\n")
    #f.write(" par%gamma    = 0.0, 0.0, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0\n")

    #f.write(" par%out_bitpix = -64\n")
    f.write(" par%%nprint     = %6.1e\n" % (no_photons/10.0))
    f.write("/\n")
    f.close()

    #--- run file.
    HOST = socket.gethostname()
    if HOST == 'polaris':
       jobname = line + '_b%02d' % bturb
       if source_rscale <= 0.0: jobname = jobname + 'p'
       make_runfile_polaris(jobname=jobname,run_file='run_polaris.sh',file_in=file_in,hybrid=hybrid)
    else:
       make_runfile(run_file='run.sh',file_in=file_in,hybrid=hybrid)

#--------------------------------
global EXEC
#--------------------------------
line       = 'SiII_1193'
hybrid     = False
use_stokes = True
geometry        = 'sphere'
source_geometry = 'point'
DGR_HI          = 0.0
spectral_type   = 'continuum'

tau_arr  = [1.0, 2.0, 5e0, 1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3]
Vexp_arr = [0.0, 50.0, 100.0, 200.0, 400.0]

#if hybrid == True:
#   EXEC='../../LaRT_v1.40_hybrid/LaRT.x'
#else:
#   EXEC='../../LaRT_v1.40/LaRT.x'
EXEC='../LaRT.x'

tau_arr  = np.array(tau_arr)
Vexp_arr = np.array(Vexp_arr)
ntau     = tau_arr.size
nVexp    = Vexp_arr.size

for j in np.arange(ntau):
   for i in np.arange(nVexp):
      tau  = tau_arr[j]
      Vexp = Vexp_arr[i]
      no_photons = 1e8

      make_input(line=line, use_stokes=use_stokes, spectral_type=spectral_type,
                 tau=tau, Vexp=Vexp, geometry=geometry, source_geometry=source_geometry, no_photons=no_photons, DGR_HI=DGR_HI)
