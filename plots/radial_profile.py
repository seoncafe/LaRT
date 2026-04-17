#!/usr/bin/env python
import numpy as np
def radial_profile(data, center=None, normalize=True, whole_area=False):
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
    r    = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r    = np.floor(r - roff).astype(np.int)
    rarr = np.unique(r)
    rmax = np.amax(x - center[0])

    if normalize == True:
       rarr = rarr / rmax
       rmax = 1.0

    tbin = np.bincount(r.ravel(), weights=data.ravel())
    nr   = np.bincount(r.ravel())
    radialprofile = tbin / nr

    if whole_area == False:
       w    = np.where(rarr <= rmax)
       rarr = rarr[w]
       radialprofile = radialprofile[w]
    return rarr, radialprofile
