""" Python scripts to perform 1D and 2D deconvolution.

  :CONTENTS:

           :func:`dsa`      --  1D: Difference Spectral Analysis

           :func:`dsamulti` --  1D: DSA for sets of multiple spectra

           :func:`dia`      --  2D: Difference Image Analysis

           :func:`dialin`   --  2D DIA, with a linearly-varying kernel


  :REQUIREMENTS:
      :doc:`analysis` (for :func:`dsamulti`)

      :doc:`nsdata` (if verbose options are used in :func:`dia` or :func:`dialin`)

      :doc:`numpy`,   :doc:`matplotlib`, :doc:`sys`

  :HISTORY:
    2008-11-19 08:54 IJC: Created @ UCLA

    2012-04-02 14:26 IJMC: Updated documentation; code unchanged.
"""

from numpy import sum, array, arange, linalg, tile, concatenate, floor, ones, zeros, dot
from pylab import find, figure, plot, subplot, grid, legend, title, show, inv, colorbar, clim, imshow, meshgrid

def dsa(r, i, Nk, **kw): #, w=None, verbose=False, noback=False):
    """
    Computational tool for Difference Spectral Analysis (DSA)
    
    :INPUTS:
       R -- reference spectrum.  This should have the highest possible
            signal-to-noise and the highest spectral resolution.

       I -- Current spectrum to be analysed.

       Nk -- number of pixels in the desired convolution kernel

    :OPTIONS:
       w       -- weights of the pixel values in I; typically (sigma)^-2

       noback  -- do not fit for a variable background; assume constant.

       tol=1e-10 -- if matrix determinant is less than tol, use
                    pseudoinverse rather than straight matrix
                    inversion

       verbose -- Print output statements and make a plot or two

       retinv -- return a fourth output, the Least Squares inverse
                 matrix (False by default)

    :OUTPUTS:       (M, K, B, C):
       M -- R, convolved to match I

       K -- kernel used in convolution

       B -- background offset

       C -- chisquared of fit. If no weights were specified, weights
            are set to unity for this calculation.

    :OPTIONS:
       I -- inverse matrix


    :NOTES:
        Best results are obtained with proper registration of the spectra.

        Also, beware of edge effects.  As a general rule, anything within
        a kernel width of the edges is suspect.


    :SEE_ALSO:  
       :func:`dsamulti`    """
#    """Testing Bramich's algorithm for 1D spectra."""
#    Based on the 2D Bramich (2008) DIA algorithm
#    -----
#    2008-11-14 10:56 IJC: Created @ UCLA.
#    2008-11-18 11:12 IJC: Registration now works correctly
#    2008-12-09 16:10 IJC: Somewhat optimized
#    2009-02-26 22:06 IJC: Added retinv, changed optional input format

    defaults = dict(verbose=False, w=None, noback=False, tol=1e-10, \
                        retinv=False)
    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]
    verbose = bool(kw['verbose'])
    noback = bool(kw['noback'])
    retinv = bool(kw['retinv'])
    w = kw['w']
    if verbose:
        print "kw>>" + str(kw)

    if noback:
        if verbose: print "Not fitting for a variable background..."

    tol = 1e-10  # tolerance for singularity

    r = array(r, copy=True)
    i = array(i, copy=True)
    dx = int(floor(Nk/2))
    Nk = int(Nk)  # length of kernel

    if w==None:
        w = ones(len(r), dtype=float)

    Nr = len(r)  # length of Referene
    ind = arange(Nr-Nk+1, dtype=int)
    wind = w[ind]
    if verbose: 
        #print "r>>" + str(r)
        #print "i>>" + str(i)
        #print "ind>>" + str(ind)
        print ""
        
    if noback:    
        U = zeros((Nk,Nk), dtype=float)
        b = zeros(Nk, dtype=float)
    else:
        U = zeros((Nk+1,Nk+1), dtype=float)
        b = zeros(Nk+1, dtype=float)

    # Build the b vector and U matrix
    tempval0 = w[ind+dx] * i[ind+dx]
    for p in range(Nk):
        b[p] = (tempval0 * r[ind+p]).sum()
        tempval2 = wind*r[ind+p]
        for q in range(p, Nk):
            U[p,q] = (tempval2 * r[ind+q]).sum()
            U[q,p] = U[p,q]

    if not noback:
        b[Nk] = (w[ind+dx] * i[ind+dx]).sum()
        for q in range(Nk):
            U[Nk, q] = (wind * r[ind+q]).sum()
            U[q, Nk] = U[Nk, q]

        U[Nk,Nk] = wind.sum()
    
    detU = linalg.det(U)
    if verbose: print "det(U) is:  " + str(detU)

    if detU<tol:
        print "Singular matrix: det(U) < tol.  Using pseudoinverse..."
        if verbose: 
            print 'U>>',U
        invmat = linalg.pinv(U)
    else:
        invmat = inv(U)

    a = dot(invmat, b)

    if noback:
        K = a
        B0 = 0.0
    else:
        K = a[0:len(a)-1]
        B0 = a[-1]

    m = rconvolve1d(r, K, mode='valid') + B0

    chisq  = ( wind * (i[ind] - m[ind])**2 ).sum()

    if verbose:
        chisq0 = ( wind * (i[ind] - r[ind])**2 ).sum()
        #print "Kernel is:  " + str(K)
        print "Background: " + str(B0)
        print "For the (" + str(Nr) + " - " + str(Nk+1) + ") = " + str(Nr-Nk-1) + " DOF:"
        print "Red. Chisquared (I-R): " + str(chisq0/(Nr-Nk-1))
        print "Red. Chisquared (I-M): " + str(chisq/(Nr-Nk-1))
    
        figure(); subplot(311)
        plot(r, '--'); plot(i, '-x'); plot(m, '-..'); legend('RIM'); 
        subplot(312); plot(r - i, '--'); plot(m - i); legend(['R-I', 'M-I'])
        subplot(313); plot(K, '-o'); grid('on'); legend(['Kernel']); 

    if retinv:
        return (m, K, B0, chisq, invmat)
    else:
        return (m, K, B0, chisq)

def dsamulti(mod, obs, N, **kw):
    """
    Run LSD using a given set of spectra and a model delta-function spectrum.

    :INPUTS:
       mod -- length L array containing a model delta-functon spectrum,
               made with a model linelist and nsdata.linespec

       obs -- L X M array containing M observations of spectra of length L

       N   -- size of deconvolution kernel to solve for

    :OPTIONS:
       w=None         -- L x M array containing the statistical weights of 'obs'

       meansub=False  -- subtract a weighted mean from obs (along axis=1)

    :OUTPUT:
       lsdbox -- M x N array of computer deconvolution kernels

       varbox -- M x N array of statistical variance-uncertainties on 'lsdbox'

       chisq  -- length M array with the chi-squared values from the LSD

    :SEE_ALSO: 
      :func:`dsa`
    """
    #2009-03-01 09:51 IJC: Created to compartmentalize computations

    from numpy import diag, sqrt, ones
    from analysis import wmean, wstd
    import sys

    obs = array(obs).copy()
    mod = array(mod).copy()

    defaults = dict(w=None, meansub=False, verbose=False, noback=False)
    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]
    verbose = bool(kw['verbose'])
    meansub = bool(kw['meansub'])
    noback = bool(kw['noback'])

    osh = obs.shape
    if len(osh)==1:
        osh = osh.reshape(osh[0], 1)
    nlam = osh[0]
    nobs = osh[1]
    if kw['w']==None:
        w = ones((nlam, nobs), float)
    else:
        w = kw['w']

    if len(mod)<>nlam:
        print "Model and observed spectra must have the same length"
        return -1
    if obs.shape<>w.shape:
        print "Observed spectra and weights array must be the same size"
        return -1

    if meansub:
        meanobs = wmean(obs, w, 1) 
        meanerr = wstd(obs, w, 1)/sqrt(nobs)
        if verbose:
            print "meanobs.shape>> " + str(meanobs.shape)
            print "meanerr.shape>> " + str(meanerr.shape)
            print "obs.shape>> " + str(obs.shape)
            print "meanerr>>" + str(meanerr)
            print "w>>" + str(w)

        obs = obs - meanobs
        w = (1./w + meanerr**2)**(-1)   # Add variances in quadrature
        if verbose:
            print "w>>" + str(w)
        
    lsdbox = zeros((nobs, N), float)   # LSD profiles
    varbox = zeros((nobs, N), float)   # variances
    chisq  = zeros(nobs, float)
    for ii in range(nobs):
        ret = dsa(mod, obs[:,ii], N, noback=noback, \
                      verbose=verbose, w=w[:,ii], retinv=True)
        lsdbox[ii,:] = ret[1]
        varbox[ii,:] = diag(ret[4])
        chisq[ii] = ret[3]
        print str(ii+1) + '/' + str(nobs),
        sys.stdout.flush()
        

    return (lsdbox, varbox, chisq)

def dialin(r, i, k, w=None, verbose=True, noback=False):
    """
    Computational tool for Difference Image Analysis (DIA) with a
       linearly-varying kernel
    
    :INPUTS:
       R -- reference image.  This should have the highest possible
            signal-to-noise and the sharpest PSF.

       I -- Current image to be analysed.

       k -- pixel mask for kernel; 1 for pixels to be used, 0 for
            pixels to be ignored

    :OPTIONS:
       w       -- weights of the pixel values in I; typically (sigma)^-2

       noback  -- do not fit for a variable background; assume constant.

       verbose -- Print output statements and make a plot or two

    :OUTPUTS:       (M, K, X, Y, B, C):
       M -- R, convolved to match I

       K -- kernel used in convolution

       B -- background offset

       C -- chisquared of fit. If no weights were specified, weights
            are set to unity for this calculation.


    :NOTES:
        Best results are obtained with proper registration of the images.

        Also, beware of edge effects.  As a general rule, anything within
        a kernel width of the edges is suspect.

        Originally based on the 2D Bramich (2008) DIA algorithm

        2008-11-19 21:50 IJC: Trying it out
    """
#    """Testing Bramich's algorithm for 2D DIA."""

    from nsdata import imshow

    if noback:
        if verbose: print "Not fitting for a variable background..."

    tol = 1e-10  # tolerance for singularity

    r = array(r, copy=True)
    i = array(i, copy=True)
    k = array(k, copy=True, dtype=bool)
    if len(k.ravel())==1:
        k = k.reshape((1,1))
    elif len(k.shape)==1:
        k = k.reshape((k.shape[0], 1))

    Nrx = r.shape[0]
    Nry = r.shape[1]
    Nr  = Nrx * Nry
    Nkx = k.shape[0]
    Nky = k.shape[1]
    Nk  = k.sum()
    dx  = int(floor(Nkx/2))
    dy  = int(floor(Nky/2))
    ix  = Nrx - Nkx + 1
    iy  = Nry - Nky + 1

    pvec = find(k.ravel())
    pl = pvec/Nkx
    pm = pvec % Nky

    if w==None:
        w = ones(i.shape, dtype=float)

    #ind = arange(Nr-Nk+1, dtype=int)

    if verbose: print "Nrx,Nry,Nr>>" + str((Nrx,Nry,Nr))
    if verbose: print "r>>" + str(r)
    if verbose: print "i>>" + str(i)
    if verbose: print "Nkx,Nky,Nk>>" + str((Nkx,Nky,Nk))
    if verbose: print "k>>" + str(k)
    if verbose: print "pvec>>" + str(pvec)
    if verbose: print "pl>>" + str(pl)
    if verbose: print "pm>>" + str(pm)
    if verbose: print "dx,dy>>" + str((dx,dy))
    if verbose: print "ix,iy>>" + str((ix,iy))

    xcoords = arange(Nrx) - floor(Nrx/2)
    ycoords = arange(Nry) - floor(Nry/2)
    xx,yy = meshgrid(xcoords, ycoords)
    if verbose: print "xx>>" + str(xx)
    if verbose: print "yy>>" + str(yy)
    
    b = zeros(3*Nk+1, dtype=float)
    itemp = i[dx:dx+ix,dy:dy+iy]
    wtemp = w[dx:dx+ix,dy:dy+iy]
    iwtemp = itemp * wtemp
    #xtemp = xx[dx:dx+ix,dy:dy+iy]
    #ytemp = yy[dx:dx+ix,dy:dy+iy]
    
    if verbose: print "itemp>>" + str(itemp)
    #if verbose: print "xtemp>>" + str(ytemp)
    #if verbose: print "ytemp>>" + str(xtemp)
    for ii in arange(Nk):
        l = pl[ii]
        m = pm[ii]
        b[ii]      = (iwtemp * r[l:l+ix,m:m+iy]).sum()
        b[ii+Nk]   = (iwtemp * r[l:l+ix,m:m+iy] * xx[l:l+ix,m:m+iy]).sum()
        b[ii+2*Nk] = (iwtemp * r[l:l+ix,m:m+iy] * yy[l:l+ix,m:m+iy]).sum()

    b[3*Nk] = iwtemp.sum()

    if verbose: print "b>>" + str(b)
    if verbose: print "Made it through 'b' Calculation!"
    
    # Construct the U_pq matrix; here p is "ii" and q is "jj"
    U = zeros((3*Nk+1,3*Nk+1))
    for ii in arange(Nk):
        l = pl[ii]
        m = pm[ii]
        xtemp = xx[l:l+ix,m:m+iy]
        ytemp = yy[l:l+ix,m:m+iy]
        rtemp = r[l:l+ix, m:m+iy]
        #if verbose: print "xtemp.shape>>" + str(xtemp.shape)
        #if verbose: print "ytemp.shape>>" + str(ytemp.shape)

        # Compute the final row and column:
        U[ii,-1]      = ( wtemp * rtemp ).sum()
        U[-1,ii]      = U[ii,-1]
        U[ii+Nk,-1]   = ( wtemp * rtemp * xtemp ).sum()
        U[-1,ii+Nk]   = U[ii+Nk,-1]
        U[ii+2*Nk,-1] = ( wtemp * rtemp * ytemp ).sum()
        
        # Compute the guts of the matrix:
        for jj in arange(ii+1):
            l2 = pl[jj]
            m2 = pm[jj]
            if verbose: print "ii,jj, l, m, l2, m2>>" + str((ii,jj,l,m,l2,m2))
            xtemp2 = xx[l2:l2+ix,m2:m2+iy]
            ytemp2 = yy[l2:l2+ix,m2:m2+iy]
            rtemp2 = r[l2:l2+ix, m2:m2+iy]
            #if verbose: print "xtemp2.shape>>" + str(xtemp2.shape)
            #if verbose: print "ytemp2.shape>>" + str(ytemp2.shape)

            wrr2 = wtemp * rtemp * rtemp2

            U[ii,jj] =           ( wrr2 ).sum()
            U[ii+Nk,jj]        = ( wrr2 * xtemp ).sum()
            U[ii+Nk,jj+Nk]     = ( wrr2 * xtemp * xtemp2 ).sum()
            U[ii+2*Nk,jj]      = ( wrr2 * ytemp ).sum()
            U[ii+2*Nk,jj+Nk]   = ( wrr2 * ytemp * xtemp2 ).sum()
            U[ii+2*Nk,jj+2*Nk] = ( wrr2 * ytemp * ytemp2 ).sum()

            if ii>jj:   # Fill out the upper diagonal:
                if verbose: print "Filling upper diag. component " + str((ii,jj))
                U[jj     ,ii     ] = U[ii     ,jj     ]
                U[jj     ,ii+Nk  ] = U[ii+Nk  ,jj     ]
                U[jj+Nk  ,ii+Nk  ] = U[ii+Nk  ,jj+Nk  ]
                U[jj     ,ii+2*Nk] = U[ii+2*Nk,jj     ]
                U[jj+Nk  ,ii+2*Nk] = U[ii+2*Nk,jj+Nk  ]
                U[jj+2*Nk,ii+2*Nk] = U[ii+2*Nk,jj+2*Nk]

    U[-1,-1] = wtemp.sum()

    if verbose: print "U>>" + str(U)

    if noback:
        U = U[0:3*Nk, 0:3*Nk]
        b = b[0:3*Nk]

    detU = linalg.det(U)
    if verbose: print "det(U) is:  " + str(detU)

    if detU<tol:
        if verbose: print "Singular matrix: det(U) < tol.  Using pseudoinverse..."
        a = dot(linalg.pinv(U), b)
    else:
        a = dot(inv(U), b)
    
    if noback:
        K = a
        B0 = 0.0
    else:
        K = a[0:len(a)-1]
        B0 = a[-1]
    if verbose: print "K>>" + str(K)
    if verbose: print "B0>>" + str(B0)
    if verbose: print "find(k.ravel())>>" + str(find(k.ravel()))

    kernel  = zeros(Nkx*Nky, dtype=float)
    xkernel = zeros(Nkx*Nky, dtype=float)
    ykernel = zeros(Nkx*Nky, dtype=float)
    kernel[find(k.ravel())] = K[0:Nk]
    xkernel[find(k.ravel())] = K[Nk:2*Nk]
    ykernel[find(k.ravel())] = K[2*Nk:3*Nk]

    kernel  =  kernel.reshape((Nkx, Nky))
    xkernel = xkernel.reshape((Nkx, Nky))
    ykernel = ykernel.reshape((Nkx, Nky))
                        
    # Bramich's convolution is not Numpy-standard, so we do it by hand:
    #m = rconvolve2d(r, kernel, mode='valid') + B0

    #chisq  = ( wtemp * (itemp - m[dx:dx+ix,dy:dy+iy])**2 ).sum()
    #chisq0 = ( wtemp * (itemp - r[dx:dx+ix,dy:dy+iy])**2 ).sum()

    if verbose: print "Kernel is:  " + str(kernel)
    if verbose: print "X-Kernel is:  " + str(xkernel)
    if verbose: print "Y-Kernel is:  " + str(ykernel)
    if verbose: print "Background: " + str(B0)
    if verbose: print "Phot. scaling: " + str(kernel.sum())
    if verbose: print "For the (" + str(Nr) + " - " + str(3*Nk+1) + ") = " + str(Nr-3*Nk-1) + " DOF:"
    #if verbose: print "Red. Chisquared (I-R): " + str(chisq0/(Nr-Nk-1))
    #if verbose: print "Red. Chisquared (I-M): " + str(chisq/(Nr-Nk-1))
    
    if verbose:
        ndim=3
        for ii in range(ndim):
            for jj in range(ndim):
                figure(1); s = subplot(ndim, ndim, 1+jj + (ndim-ii-1)*ndim)
                xind = int(Nrx/(ndim+1.0))
                yind = int(Nry/(ndim+1.0))
                x0 = xx[(ii+1)*xind,(jj+1)*yind]
                y0 = yy[(ii+1)*xind,(jj+1)*yind]
                print "x0,y0>>" + str((x0,y0))
                imshow(kernel + x0*xkernel + y0*ykernel)
                s.set_ylim(s.get_ylim()[::-1])
                offset = array(centroid(kernel + x0*xkernel + y0*ykernel)) - \
                    array(centroid(ones(kernel.shape)))
                title(str((x0,y0)) + ':  ' + str(offset))
                #figure(2)
                #plot([x0],[y0], 'oc')
                #plot([x0, x0+offset[0]], [y0, y0+offset[1]], '-k')
        figure(); 
        subplot(121); imshow(r); title('Reference'); #clim([cmin, cmax]); colorbar()
        subplot(122); imshow(i); title('Current Image'); #clim([cmin, cmax]); colorbar()
        #subplot(133); imshow(m); title('M'); clim([cmin, cmax]); colorbar()

#       figure(); 
#       cmin = array([(r-i).min(), (m-i).min()]).min()
#       cmax = array([(r-i).max(), (m-i).min()]).max()
#       subplot(121); imshow(r - i); title('R - I'); clim([cmin, cmax]); colorbar();
#       subplot(122); imshow(m - i); title('M - I'); clim([cmin, cmax]); colorbar(); 

        show()

    return

#    return (m, kernel, B0, chisq)


def dia(r, i, k, w=None, verbose=False, noback=False):
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

    :EXAMPLE:
      ::

                import pylab as py
                import dia
                
                # Define reference and blurred images:
                ref = py.zeros((10,10))
                img = py.zeros((10,10))
                ref[3,3] = 1; ref[3,6] = 1;  ref[6,3] = 1
                img[2,3] = 1; img[3,2:5] = 1;  img[4,3] = 1
                img[2,6] = 1; img[3,5:8] = 1;  img[4,6] = 1
                img[5,3] = 1; img[6,2:5] = 1;  img[7,3] = 1

                # Add some noise:
                img += py.randn(10, 10)/10.

                # Define kernel basis
                kb = py.ones((3,3))    

                # Compute Difference Image Analysis:
                m, kern, bkg, chisq = dia.dia(ref, img, kb)

                # Display results:
                py.figure()
                py.subplot(231)
                py.imshow(ref)
                py.title('Reference image')
                py.subplot(232)
                py.imshow(img)
                py.title('Observed image')
                py.subplot(234)
                py.imshow(kern)
                py.title('Optimal kernel')
                py.subplot(235)
                py.imshow(m)
                py.title('Convolved Reference')
                py.subplot(236)
                py.imshow(img - m)
                py.title('Residuals')
                



    :NOTES:
        Best results are obtained with proper registration of the images.

        Also, beware of edge effects.  As a general rule, anything within
        a kernel width of the edges is suspect.

        Based on the 2D Bramich (2008) DIA algorithm


        2008-11-18 11:12 IJC: Extrapolating from my 1D algorithm

        2012-04-03 08:19 IJMC: Added example to documentation.
    """
#    """Testing Bramich's algorithm for 2D DIA."""

    from nsdata import imshow

    if noback:
        if verbose: print "Not fitting for a variable background..."

    tol = 1e-10  # tolerance for singularity

    r = array(r, copy=True)
    i = array(i, copy=True)
    k = array(k, copy=True, dtype=bool)
    if len(k.ravel())==1:
        k = k.reshape((1,1))
    elif len(k.shape)==1:
        k = k.reshape((k.shape[0], 1))

    Nrx = r.shape[0]
    Nry = r.shape[1]
    Nr  = Nrx * Nry
    Nkx = k.shape[0]
    Nky = k.shape[1]
    Nk  = k.sum()
    dx  = int(floor(Nkx/2))
    dy  = int(floor(Nky/2))
    ix  = Nrx - Nkx + 1
    iy  = Nry - Nky + 1

    pvec = find(k.ravel())
    pl = pvec/Nkx
    pm = pvec % Nky

    if w==None:
        w = ones(i.shape, dtype=float)

    #ind = arange(Nr-Nk+1, dtype=int)

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

    #b = zeros(Nk+1, dtype=float)
    #for ii in arange(Nk):
    #    b[ii] = (w[ind+dx] * i[ind+dx] * r[ind+ii]).sum()
    #b[Nk] = (w[ind+dx] * i[ind+dx]).sum()

    b = zeros(Nk+1, dtype=float)
    itemp = i[dx:dx+ix,dy:dy+iy]
    wtemp = w[dx:dx+ix,dy:dy+iy]
    if verbose: print "itemp>>" + str(itemp)
    for ii in arange(Nk):
        b[ii] = (itemp * wtemp * r[pl[ii]:pl[ii]+ix,pm[ii]:pm[ii]+iy]).sum()
    b[Nk] = (w[dx:dx+ix] * i[dx:dx+ix]).sum()
    if verbose: print "b>>" + str(b)
    
    if verbose: print "Made it through 'b' Calculation!"
    
        
    U = zeros((Nk+1,Nk+1))
    # This is optimized by only computing the upper diagonal elements...
    for ii in arange(Nk):
        l = pl[ii]
        m = pm[ii]
        if verbose: print "l, m>>" + str((l,m))
        rwtemp = r[l :l +ix,m :m +ix] * wtemp
        U[Nk, ii] = rwtemp.sum()
        U[ii, Nk] = U[Nk, ii]
        for jj in arange(ii,Nk):
            l2 = pl[jj]
            m2 = pm[jj]
            if verbose: print "l2, m2>>" + str((l2,m2))
            U[ii,jj] = ( rwtemp * r[l2:l2+ix,m2:m2+ix] ).sum()
            U[jj,ii] = U[ii,jj]
    if verbose: print "U>>" + str(U)

    U[Nk, Nk] = wtemp.sum()

    if verbose: print "U>>" + str(U)

    if noback:
        U = U[0:Nk, 0:Nk]
        b = b[0:Nk]

    detU = linalg.det(U)
    if verbose: print "det(U) is:  " + str(detU)

    if detU<tol:
        if verbose: print "Singular matrix: det(U) < tol.  Using pseudoinverse..."
        a = dot(linalg.pinv(U), b)
    else:
        a = dot(inv(U), b)
    
    if noback:
        K = a
        B0 = 0.0
    else:
        K = a[0:len(a)-1]
        B0 = a[-1]
    if verbose: print "K>>" + str(K)
    if verbose: print "B0>>" + str(B0)
    if verbose: print "find(k.ravel())>>" + str(find(k.ravel()))

    kernel = zeros(Nkx*Nky, dtype=float)
    kernel[find(k.ravel())] = K

    kernel = kernel.reshape((Nkx, Nky))
                        
    # Bramich's convolution is not Numpy-standard, so we do it by hand:
    m = rconvolve2d(r, kernel, mode='valid') + B0

    if verbose: print "itemp.shape>>" + str(itemp.shape)
    if verbose: print "wtemp.shape>>" + str(wtemp.shape)
    if verbose: print "r.shape>>" + str(r.shape)
    if verbose: print "m.shape>>" + str(m.shape)
    if verbose: print "m[dx:dx+ix,dy:dy+iy].shape>>" + str(m[dx:dx+ix,dy:dy+iy].shape)

    chisq  = ( wtemp * (itemp - m[dx:dx+ix,dy:dy+iy])**2 ).sum()
    chisq0 = ( wtemp * (itemp - r[dx:dx+ix,dy:dy+iy])**2 ).sum()

    #if verbose: print "Kernel is:  " + str(kernel)
    if verbose: print "Background: " + str(B0)
    if verbose: print "Phot. scaling: " + str(kernel.sum())
    if verbose: print "For the (" + str(Nr) + " - " + str(Nk+1) + ") = " + str(Nr-Nk-1) + " DOF:"
    if verbose: print "Red. Chisquared (I-R): " + str(chisq0/(Nr-Nk-1))
    if verbose: print "Red. Chisquared (I-M): " + str(chisq/(Nr-Nk-1))
    
    if verbose:
        figure(); 
        cmin = array([i.min(), m.min()]).min()
        cmax = array([i.max(), m.max()]).max()
        subplot(131); imshow(r); title('R'); clim([cmin, cmax]); colorbar()
        subplot(132); imshow(i); title('I'); clim([cmin, cmax]); colorbar()
        subplot(133); imshow(m); title('M'); clim([cmin, cmax]); colorbar()

        figure(); 
        cmin = array([(r-i).min(), (m-i).min()]).min()
        cmax = array([(r-i).max(), (m-i).min()]).max()
        subplot(121); imshow(r - i); title('R - I'); clim([cmin, cmax]); colorbar();
        subplot(122); imshow(m - i); title('M - I'); clim([cmin, cmax]); colorbar(); 

        figure();
        imshow(kernel); title('Kernel'); colorbar()

        show()

    return (m, kernel, B0, chisq)

def rconvolve2d(a, b, mode='valid', extend=0, verbose=False):
    """
    Compute a 2D reverse convolution in the style of Bramich 2008.

    :INPUTS:
        'a' should be larger than 'b' -- i.e., 'b' is the kernel -- and
        both should be square.

        'extend' tells how to extend the boundaries -- either
        'nearest'-neighbor or a number
        
    :NOTES:
        This is "reversed" from the canonical definition of the convolution.

        Note that I could also implement an FFT convolution algorithm, for
        increased speed.

    :SEE_ALSO: 
      :func:`dsa`, :func:`rconvolve1d`
    """
    # 2008-11-17 21:57 IJC: Created @ UCLA

    a = array(a, copy=True)
    b = array(b, copy=True)
    nax = a.shape[0]
    nbx = b.shape[0]
    nay = a.shape[1]
    nby = b.shape[1]
    nx = max(nax, nbx)
    ny = max(nay, nby)
    
    dx0 = int(floor(nbx/2))
    dy0 = int(floor(nby/2))
    dx1 = nbx - dx0
    dy1 = nby - dy0

    if verbose: print "dx0,dy0 >>" + str(dx0) + ", " + str(dy0)

    if extend=='nearest':
        X = a.ravel()[-1]
    else:
        X = extend

    if verbose: print "(nax+nbx-1)>>" + str(nax+nbx-1)
    if verbose: print "(nay+nby-1)>>" + str(nay+nby-1)
    if verbose: print "X>>" + str(X)

    a2 = zeros((nax+nbx-1,nay+nby-1), dtype=float) + X
    a2[dx0:nax+dx0,dy0:nay+dy0] = a

    c = zeros((nx, ny), dtype=float)

    for ii in range(nax):
        for jj in range(nay):
            if verbose: print "ii, jj>> " + str(ii) + ", " + str(jj)
            if verbose: print "ii, jj>> " + str(ii) + ", " + str(jj)
            #if verbose: print "a2[ii:ii+nb+1,ii:ii+nb+1]>>" + str(a2[ii:ii+nb,ii:ii+nb])
            #if verbose: print "b>>" + str(b)
            c[ii, jj] = (a2[ii:ii+nbx,jj:jj+nby] * b).sum()

    return c

def rconvolve1d(a, b, mode='valid', extend='nearest'):
    """
    Compute a 1D reverse convolution in the style of Bramich 2008.

    :INPUTS:
        'a' should be longer than 'b' -- i.e., 'b' is the kernel.

        'extend' tells how to extend the boundaries -- either
        'nearest'-neighbor or a number

    :NOTES:
      This is "reversed" from the canonical definition of the convolution.

    :SEE_ALSO:   
      :func:`dsa`
    """
    # 2008-11-14 18:26 IJC: Created
    na = len(a)
    nb = len(b)
    n = max(na, nb)
    
    dx = int(floor(nb/2))

    if extend=='nearest':
        X = a[-1]
    else:
        X = extend

    a2 = X + zeros(na+nb-1, dtype=float)
    a2[dx:dx+na] = a
    #a2 = concatenate((a, X + zeros(nb-1)))
    #c = zeros(n, dtype='float')
    
    bmat = tile(b, (n,1))
    amat = zeros((n, nb), dtype='float')
    for ii in range(na):
        amat[ii,:] = a2[range(ii,ii+nb)]

    c = sum(amat * bmat, axis=1)
        
    return c



def exam(case=1):
    if case==1:  # Single-pixel intensity variation
        r = zeros((5,5))
        i = zeros((5,5))
        r[1,1] = 1;   r[1,3] = 1;   r[3,1] = 1
        i[1,1] = 1;   i[1,3] = 3;   i[3,1] = 5
        k = ones(1)
    elif case==2: # Constant broadening and intensity variation 
        r = zeros((8,8))
        i = zeros((8,8))
        r[2,2] = 1;   r[2,5] = 1;  r[5,2] = 1
        i[1,2] = 1; i[2,1:4] = 1;  i[3,2] = 1
        i[1,5] = 4; i[2,4:7] = 4;  i[3,5] = 4
        i[4,2] = 7; i[5,1:4] = 7;  i[6,2] = 7
        k = zeros((3,3))
        k[0,1] = 1;  k[1,:] = 1; k[2,1] = 1
    elif case==3:
        r = zeros((10,10))
        i = zeros((10,10))
        r[3,3] = 1;   r[3,6] = 1;  r[6,3] = 1
        i[2,3] = 1; i[3,2:5] = 1;  i[4,3] = 1
        i[2,6] = 1; i[3,5:8] = 1;  i[4,6] = 1
        i[5,3] = 1; i[6,2:5] = 1;  i[7,3] = 1
        k = zeros((3,3))
        k[0,1] = 1;  k[1,:] = 1; k[2,1] = 1
    elif case==4:  # Barrel distortion
        r = zeros((9,9))
        i = zeros((9,9))
        r[3,3] = 1; r[3,5] = 1; r[5,3] = 1; r[5,5] = 1
        i[2,2] = 1; i[2,6] = 1; i[6,2] = 1; i[6,6] = 1
        k = ones((3,3))
    elif case==5:  # Less sparse Barrel distortion
        r = zeros((9,9))
        i = zeros((9,9))
        r[3:6,3:6] = 1;
        i[2,2] = 1; i[2,4] = 1; i[2,6] = 1; 
        i[4,2] = 1; i[4,4] = 1; i[4,6] = 1;
        i[6,2] = 1; i[6,4] = 1; i[6,6] = 1
        k = ones((3,3))
    elif case==6: # Variable broadening, no intensity variations
        r = zeros((8,8))
        i = zeros((8,8))
        r[2,2] = 1;   r[2,5] = 1;  r[5,2] = 1
        i[1,2] = 1;   i[2,1] = 1;   i[2,2] = 4; i[2,3] = 1;    i[3,2] = 1
        i[1,5] = 0.5; i[2,4] = 0.5; i[2,5] = 6; i[2,6] = 0.5;  i[3,5] = 0.5
        i[4,2] = 1.5; i[5,1] = 1.5; i[5,2] = 2; i[5,3] = 1.5;  i[6,2] = 1.5
        k = zeros((3,3))
        k[0,1] = 1;  k[1,:] = 1; k[2,1] = 1
    elif case==7: # Variable broadening, no intensity variations
        r = zeros((8,8))
        i = zeros((8,8))
        r[2,2] = 1;   r[2,5] = 1;   r[5,2] = 1
        i[1,2] = 1;   i[2,1] = 1;   i[2,2] = 6; i[2,3] = 1;      i[3,2] = 1
        i[1,5] = 1;   i[2,4] = 1;   i[2,5] = 3.5; i[2,6] = 3.5;  i[3,5] = 1
        i[4,2] = 1;   i[5,1] = 1;   i[5,2] = 1.5; i[5,3] = 1;    i[6,2] = 5.5
        k = zeros((3,3))
        k[0,1] = 1;  k[1,:] = 1; k[2,1] = 1
    elif case==8: # Variable broadening, no intensity variations
        r = zeros((8,8))
        i = zeros((8,8))
        r[2,2] = 1;   r[2,5] = 1;   r[5,2] = 1;  r[5,5] = 1
        i[1,2] = 1;   i[2,1] = 1;   i[2,2] = 8;   i[2,3] = 1;    i[3,2] = 1
        i[1,5] = 1;   i[2,4] = 1;   i[2,5] = 5.5; i[2,6] = 3.5;  i[3,5] = 1
        i[4,2] = 1;   i[5,1] = 1;   i[5,2] = 3.5; i[5,3] = 1;    i[6,2] = 5.5
        i[4,5] = 1;   i[5,4] = 1;   i[5,5] = 1;   i[5,6] = 3.5;  i[6,5] = 5.5
        k = zeros((3,3))
        k[0,1] = 1;  k[1,:] = 1; k[2,1] = 1

    dialin(r, i, k, noback=True, verbose=True)


    return

def centroid(a):
    """ Return the pixel coordinates of the centroid of an array.

    (xc, yc) = centroid(a)


    Uses the algorithm from: http://www.dfanning.com/tips/centroid.html
    """
    # 2008-11-21 11:31 IJC: Created
    a = array(a, copy=True)
    ash = a.shape
    x = arange(ash[0])
    y = arange(ash[1])
    
    xc = (x * a.sum(axis=1)).sum()/a.sum()
    yc = (y * a.sum(axis=0)).sum()/a.sum()

    return (xc,yc)

