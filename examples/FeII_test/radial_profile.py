#!/usr/bin/env python
import numpy as np

def radial_profile(data, center=None, normalize=True, whole_area=False, rbins=None):
    #---------------------------
    # Note that boundaries for an array (0,1,...,nx-1) is assumed to be (0,1,...,nx).
    # so that the centeral coordinates of a pixel (i,j) is (i+0.5,j+0.5).
    # 2020.10.23, updated for an option of normalize, whole_area
    # 2020.10.09, initial version KI Seon.
    #---------------------------

    ny, nx = data.shape
    if center == None:
       center = [nx/2.0, ny/2.0]
    if np.int32(nx/2)*2 != nx and np.int32(ny/2)*2 != ny:
       roff = -0.5
    else:
       roff =  0.0

    y, x = np.indices((data.shape)) + 0.5
    r0   = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r0   = r0 - roff

    r     = np.floor(r0).astype(np.int32)
    rarr  = np.unique(r)
    rmax  = np.amax(x - center[0])
    rmax0 = rmax

    if normalize == True:
       rarr = rarr / rmax
       rmax = 1.0

    if np.any(rbins) != True:
       tbin  = np.bincount(r.ravel(), weights=data.ravel())
       t2bin = np.bincount(r.ravel(), weights=data.ravel()**2)
       nr    = np.bincount(r.ravel())
       radialprofile = tbin / nr
       errorprofile  = np.sqrt(t2bin / nr - radialprofile**2)
       if whole_area == False:
          w    = np.where(rarr <= rmax)
          rarr = rarr[w]
          radialprofile = radialprofile[w]
          errorprofile  = errorprofile[w]
    else:
       r     = r0 / rmax0
       nbins = len(rbins)-1
       rarr  = (rbins[0:-1] + rbins[1:])/2.0
       if rbins[0] <= 0.0: rarr[0] = 0.0

       radialprofile = np.zeros(nbins)
       errorprofile  = np.zeros(nbins)
       for i in np.arange(nbins):
          w                = np.where((r >= rbins[i]) & (r < rbins[i+1]))
          ncount           = w[0].size
          radialprofile[i] = np.sum(data[w]) / ncount
          p2               = np.sum(data[w]**2)
          errorprofile[i]  = np.sqrt(p2 / ncount - radialprofile[i]**2)

    return rarr, radialprofile, errorprofile

def radial_spectrum(data, center=None, normalize=True, whole_area=False, rbins=None):
    #----
    ny, nx, nwavl = data.shape
    for i in np.arange(nwavl):
       rarr, prof, err = radial_profile(data[:,:,i],rbins=rbins)
       if i == 0:
          nr   = rarr.size
          spec = np.zeros([nr, nwavl])
       spec[:,i] = prof

    return rarr, spec
