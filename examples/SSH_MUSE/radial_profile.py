#!/usr/bin/env python
import numpy as np
def radial_profile(data, center=None, normalize=True, whole_area=False):
    #---------------------------
    # Note that boundaries for an array (0,1,...,nx-1) is assumed to be (0,1,...,nx).
    # so that the centeral coordinates of a pixel (i,j) is (i+0.5,j+0.5).
    # 2021.09.12, updated to deal with the case where x-axis and y-axis are different.
    # 2020.10.23, updated for an option of normalize, whole_area
    # 2020.10.09, initial version KI Seon.
    #---------------------------

    ny, nx = data.shape
    if center == None:
       center = [nx/2.0, ny/2.0]
    if (nx//2)*2 != nx and (ny//2)*2 != ny:
       roff = -0.5
    else:
       roff =  0.0

    y, x = np.indices((data.shape)) + 0.5
    r    = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r    = np.floor(r - roff).astype(np.int)
    rarr = np.unique(r)
    xmax = np.amax(x - center[0])
    ymax = np.amax(y - center[1])
    rmax = np.amax([xmax, ymax])

    tbin = np.bincount(r.ravel(), weights=data.ravel())
    nr   = np.bincount(r.ravel())
    radialprofile = tbin / nr

    if whole_area == False:
       rmax = np.amin([xmax, ymax])
       w    = np.where(rarr <= rmax)
       rarr = rarr[w]
       radialprofile = radialprofile[w]

    if normalize == True:
       rarr = rarr / rmax
       rmax = 1.0

    return rarr, radialprofile
