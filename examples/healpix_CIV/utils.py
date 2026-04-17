#!/usr/bin/env python
import numpy as np

def find_histogram_limits(x, bins=200, log=True, perc=0.01, histogram=False):
   if log == True:
      w = np.where(x > 0.0)
      hist, bin_edges = np.histogram(np.log10(x[w]), bins=bins, density=True)
   else:
      hist, bin_edges = np.histogram(x, bins=bins, density=True)
   hist_cum = np.cumsum(hist)
   hist_cum = hist_cum / np.amax(hist_cum)
   x        = (bin_edges[1:] + bin_edges[0:-1])/2.0
   xmin     = np.interp(perc/2.0, hist_cum, x)
   xmax     = np.interp(1.0-perc/2.0, hist_cum, x)
   if log == True:
      xmin = 10.0**xmin
      xmax = 10.0**xmax
      bin_edges = 10.0**bin_edges
   #print(xmin, xmax, hist_cum[0], hist_cum[-1])
   if histogram == True:
      return hist, bin_edges
   else:
      return xmin, xmax

def my_histogram(x, bins=200, log=True, perc=0.01):
   hist, bin_edges = find_histogram_limits(x, bins=bins, log=log, perc=perc, histogram=True)
   bin_centers     = (bin_edges[1:]+bin_edges[0:-1])/2.0
   return hist, bin_centers

def print_time_stamp():
   from datetime import datetime
   global started
   try:
      from mpi4py import MPI
      comm = MPI.COMM_WORLD
      rank = comm.Get_rank()
   except:
      rank = 0

   if (rank == 0):
      try:
         started
      except:
         print('\nStart: %s' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'), flush=True)
         started = 1
      else:
         print('End  : %s\n' % datetime.now().strftime('%Y-%m-%d %H:%M:%S'), flush=True)
         del started

def get_time():
   global first_time
   try:
      from mpi4py import MPI
      wtime = MPI.Wtime()
   except:
      import time
      wtime = time.time()

   try:
      first_time
   except:
      first_time = wtime
      delta_time = 0.0
   else:
      delta_time = wtime - first_time
      del first_time

   return delta_time
