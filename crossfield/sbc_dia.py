"""SBC's version of the Bramich 2008 (MNRAS, 386, L77) Difference Image Analysis Algorithm"""

import numpy as np
import pylab as py
import dia
import sys

def sbc_dia(r, i, k, w=None, verbose=False, noback=False):

    """
    Computational tool for Difference Image Analysis (DIA)
    
    :INPUTS:
       R -- reference image.  This should have the highest possible
            signal-to-noise and the sharpest PSF.

       I -- Current image to be analysed.

       k -- 2D kernel basis mask: 1 for pixels to be used, 0 for
            pixels to be ignored

    :OPTIONS:
       w       -- weights of the pixel values in I; typically (sigma)^-2

       noback  -- do not fit for a variable background; assume constant.

       verbose -- Print output statements and make a plot or two

    :OUTPUTS:       (M, K, B, C):
       M -- R, convolved to match I

       K -- kernel used in convolution

       B -- background offset

       C -- chisquared of fit. If no weights were specified, weights
            are set to unity for this calculation.

    Borrows heavily from Ian Crossfield code, but I think corrects some bugs.
    (http://www.mpia-hd.mpg.de/homes/ianc/python/index.html)
    """
    
    if noback:
        if verbose: print "Not fitting for a variable background..."

    tol = 1e-10  # tolerance for singularity

    r = np.array(r, copy=True)
    i = np.array(i, copy=True)
    k = np.array(k, copy=True, dtype=bool)

    # Couple of quick sanity checks
    if not (r.shape == i.shape):
        print "Bad user: input arrays must be of the same size!"
        sys.exit(1)
    if not (k.shape[0] % 2 == 1):
        print "Bad user: input kernel must be odd in size!"
        sys.exit(1)

    Nrx = r.shape[0]; Nry = r.shape[1]; Nr  = Nrx * Nry
    Nkx = k.shape[0]; Nky = k.shape[1]; Nk  = k.sum()
    # Nk = k.sum instead of Nkx * Nky in case you want to
    # weight some pixels as 0 in the kernel

    dx  = int(py.floor(Nkx/2)) + 1; dy  = int(py.floor(Nky/2)) + 1
    ix  = Nrx - Nkx + 1; iy  = Nry - Nky + 1
    # These indices will leave out the edges (where the convolution
    # is not well defined.
    
    # These define the l,m values that we will loop over
    pvec = py.find(k.ravel())
    pl = pvec/Nkx - (Nkx - 1) / 2 
    pm = (pvec % Nky) - (Nky - 1) / 2 

    if w==None:
        w = np.ones(i.shape, dtype=float)

    if verbose: 
        print "Nrx,Nry,Nr>>" + str((Nrx,Nry,Nr))
        print "r>>" + str(r)
        print "i>>" + str(i)
        print "Nkx,Nky,Nk>>" + str((Nkx,Nky,Nk))
        print "k>>" + str(k)
        print "pvec>>" + str(pvec)
        print "pl>>" + str(pl)
        print "pm>>" + str(pm)
        print "dx,dy>>" + str((dx,dy))
        print "ix,iy>>" + str((ix,iy))

    # Construct the "b" matrix (in Ua = b)
    b = np.zeros(Nk+1, dtype=float)
    # Remove edges (convolution not well defined there
    itemp = i[dx:ix,dy:iy]
    wtemp = w[dx:ix,dy:iy]
    if verbose: print "itemp>>" + str(itemp)
    # For all l,m values, want I_ij * w_ij * R_(i+l)(j+m)
    # Get all possible combinations of l,m from pl and pm
    for ii in np.arange(Nk):
        b[ii] = (itemp * wtemp * r[dx+pl[ii]:ix+pl[ii],dy+pm[ii]:iy+pm[ii]]).sum()
    # Final value is just I_ij * w_ij
    b[Nk] = (wtemp * itemp).sum()
    
    if verbose: print "b>>" + str(b)
    if verbose: print "Made it through 'b' Calculation!"

    # Now construct the "U" matrix
    U = np.zeros((Nk+1,Nk+1))

    # Will do this inefficiently with two for loops for now
    for ii in np.arange(Nk):
        l1 = pl[ii]; m1 = pm[ii]
        for jj in np.arange(Nk):
            l2 = pl[jj]; m2 = pm[jj]
            m = r[dx+l1:ix+l1,dy+m1:iy+m1] * r[dx+l2:ix+l2,dy+m2:iy+m2] * wtemp
            U[ii, jj] = m.sum()
        m2 = r[dx+l1:ix+l1,dy+m1:iy+m1] * wtemp
        U[ii, Nk] = m2.sum()
        U[Nk, ii] = m2.sum()
        
    # Slightly different formula for last row / column
    U[Nk, Nk] = wtemp.sum()

    if verbose: print "U>>" + str(U)

    if noback:
        U = U[0:Nk, 0:Nk]
        b = b[0:Nk]

    # Can we invert U?
    detU = np.linalg.det(U)
    if verbose: print "det(U) is:  " + str(detU)

    if detU<tol:
        if verbose: print "Singular matrix: det(U) < tol.  Using pseudoinverse..."
        a = np.dot(np.linalg.pinv(U), b)
    else:
        a = np.dot(np.linalg.inv(U), b)
    
    if noback:
        K = a
        B0 = 0.0
    else:
        K = a[0:len(a)-1]
        B0 = a[-1]
    if verbose: print "K>>" + str(K)
    if verbose: print "B0>>" + str(B0)
    if verbose: print "find(k.ravel())>>" + str(py.find(k.ravel()))

    kernel = np.zeros(Nkx*Nky, dtype=float)
    kernel[py.find(k.ravel())] = K

    kernel = kernel.reshape((Nkx, Nky))
                        
    # Bramich's convolution is not Numpy-standard, so we do it by hand:
    m = dia.rconvolve2d(r, kernel, mode='valid') + B0

    if verbose: print "itemp.shape>>" + str(itemp.shape)
    if verbose: print "wtemp.shape>>" + str(wtemp.shape)
    if verbose: print "r.shape>>" + str(r.shape)
    if verbose: print "m.shape>>" + str(m.shape)
    if verbose: print "m[dx:ix,dy:iy].shape>>" + str(m[dx:ix,dy:iy].shape)

    chisq  = ( wtemp * (itemp - m[dx:ix,dy:iy])**2 ).sum()
    chisq0 = ( wtemp * (itemp - r[dx:ix,dy:iy])**2 ).sum()

    #if verbose: print "Kernel is:  " + str(kernel)
    if verbose: print "Background: " + str(B0)
    if verbose: print "Phot. scaling: " + str(kernel.sum())
    if verbose: print "For the (" + str(Nr) + " - " + str(Nk+1) + ") = " + str(Nr-Nk-1) + " DOF:"
    if verbose: print "Red. Chisquared (I-R): " + str(chisq0/(Nr-Nk-1))
    if verbose: print "Red. Chisquared (I-M): " + str(chisq/(Nr-Nk-1))
        
    return (m, kernel, B0, chisq)

###########################################



