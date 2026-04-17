#!/usr/bin/env python
import numpy as np

nph   = 36
alpha = 360.0 * np.arange(nph)/nph
beta  = 60.0  * np.ones(nph)
gamma = 90.0  * np.ones(nph)

#----
ncol = 15
nrow = nph // ncol
n0   = nph - nrow*ncol

#fmtstr = '{:5.1f} ' * len(pha)
#fmtstr = '{:4.0f} ' * len(pha)

nloop = nrow + 1
for i in np.arange(nloop):
   k1, k2 = i*ncol, (i+1)*ncol
   if k2 > nph: k2 = nph
   fmtstr = '{:5.1f} ' * (k2 - k1)
   if i == 0:
      print(' par%alpha =', fmtstr.format(*alpha[k1:k2]))
   else:
      print('            ', fmtstr.format(*alpha[k1:k2]))

for i in np.arange(nloop):
   k1, k2 = i*ncol, (i+1)*ncol
   if k2 > nph: k2 = nph
   fmtstr = '{:5.1f} ' * (k2 - k1)
   if i == 0:
      print(' par%beta  =', fmtstr.format(*beta[k1:k2]))
   else:
      print('            ', fmtstr.format(*beta[k1:k2]))

for i in np.arange(nloop):
   k1, k2 = i*ncol, (i+1)*ncol
   if k2 > nph: k2 = nph
   fmtstr = '{:5.1f} ' * (k2 - k1)
   if i == 0:
      print(' par%gamma =', fmtstr.format(*gamma[k1:k2]))
   else:
      print('            ', fmtstr.format(*gamma[k1:k2]))
