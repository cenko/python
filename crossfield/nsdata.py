import os, sys
import numpy as np
from numpy import *
from scipy import interpolate
from warnings import warn
import pdb
import analysis as an


"""
 IJC's attempt at writing a useful Python library for use with
 astronomical data reduction.  Developed at UC Los Angeles.

 
 2008-06-30 17:13 IJC: Made write_exptime able to recursively use file lists

 2010-10-29 09:18 IJC: Updated documentation for Sphinx.  Removed pad,
 fix_quadnoise_old.

 2011-04-08 11:51 IJC: Moved amedian() to analysis.py.

 :REQUIREMENTS:
     :doc:`analysis`

     :doc:`matplotlib`
"""

class aperture:
    "IRAF aperture class."
    def __init__(self):
        self.filename = ''
        self.image = []
        self.nap = 0
        self.center = []
        self.low = []
        self.high = []

def intstr(num, numplaces=4):
    """A simple function to map an input number into a string padded with
    zeros (default 4).  Syntax is: out = intstr(6, numplaces=4) -->
    0006

    2008-05-27 17:12 IJC: Created"""

    formatstr = "%(#)0"+str(numplaces)+"d"

    return formatstr % {"#":int(num)}

def sfilelist(prefix, postfix, numlist, numplaces=4, delim=','):
    """2008-05-27 17:12 IJC: Create a delimited string of filenames based
    on a specified prefix, and a list of numeric values.  You can also
    specify the number of digits in the numeric part filenames;
    default is 4.  Default delimiter is a comma.

    :EXAMPLE: 
       ::
       
         files = sfilelist(prefix, postfix, numlist, numplaces=4, delim=',')

    """
    #2008-05-27 17:13 IJC:
    #2008-07-21 10:24 IJC: Updated to use list instead of start/end num.

    fnlist = ""
    for element in numlist:
        fnlist = fnlist + prefix + intstr(element, numplaces) + postfix + delim

    fnlist = fnlist[0:len(fnlist)-1]

    return fnlist

def filelist(prefix, postfix, numlist, numplaces=4):
    """2008-05-27 17:12 IJC: Create a list of filenames based on a
    specified prefix, and a list of numeric values.  You can also
    specify the number of digits in the filenames; default is 4.

    :EXAMPLE:
      ::
      
         files = filelist(prefix, postfix, numlist, numplaces=4)

    :SEE ALSO:  :func:`wfilelist`, :func:`sfilelist`, :func:`file2list`
    """
    #2008-05-27 17:13 IJC: 
    #2008-07-21 10:24 IJC: Updated to use list instead of start/end num.
    
    fnlist = []

    for element in numlist:
        fnlist = fnlist + [prefix + intstr(element, numplaces) + postfix]

    return fnlist

def wfilelist(prefix, postfix, numlist, numplaces=4, tempname="wfilelist_py.tmp"):
    """2008-05-27 17:12 IJC: Create an ASCII file of a list of filenames
    based on a specified prefix, starting number, and ending number.
    You can also specify the number of digits in the filenames;
    default is 4.  

    If the file already exists, it is overwritten.

    :EXAMPLE:
       ::

           filelist(prefix, postfix, numlist, numplaces=4, tempname='wfilelist_py.tmp')

    :SEE ALSO:  :func:`filelist`, :func:`sfilelist`, :func:`file2list`
    """#
    #2008-05-27 17:13 IJC:
    #2008-07-21 10:25 IJC: Updated to use list instead of start/end num

    f = open(tempname, "w")
    for element in numlist:
        strtowrite = prefix + intstr(element, numplaces) + postfix + "\n"
        f.write(strtowrite)

    f.close()

    return tempname

def strl2f(filename, strl, clobber=True, EOL='\n'):
    """Write a list of strings to the specified filename.

    :INPUTS:
        filename:  string, name of file to write to
        strl: list, to be written to specified file.

    Returns the filename

    :Note: this is only designed for single-depth lists 
       (i.e., no doubly-deep string lists).
    """
    # 2009-04-28 10:46 IJC: Created to convert filelist to wfilelist

    import os

    if clobber==True and os.path.isfile(filename):
        os.remove(filename)
        
    f = open(filename, 'w')
    for el in strl:
        f.write(el+EOL)
    f.close()
    
    return filename

def showfits(filename):
    """2008-05-28 13:15 IJC: Routine to load and show a FITS file.

    showfits('blah.fits')"""
    import pyfits
    from pylab import figure, imshow, cm

    im = pyfits.getdata(filename)
    imgsize = (im.shape)[0]
    figure()
    imshow(im, aspect='equal', cmap = cm.gray)

    return

def getval(filename, key, *ext, **kw):
    """Get a keyword's value from a header in a FITS file, or a list of
    files.

    Syntax is the same as pyfits.getval:

    @type filename: string or list

    @param filename: input FITS file name, or list of filenames

    @type key: string

    @param key: keyword name

    @param ext: The rest of the arguments are for extension specification.
       See L{getdata} for explanations/examples.

    @return: keyword value

    @rtype: string, integer, or float

    An extra keyword is path='' -- a path to prepend to the filenames.
    """
    # 2009-02-11 13:54 IJC: Created
    # 
    import pyfits
    
    defaults = dict(path='')

    for keyword in defaults:
        if (not kw.has_key(keyword)):
            kw[keyword] = defaults[keyword]
    
    if (filename.__class__ <> str) and (len(filename)>1):
        ret = []
        for file in filename:
            ret.append(pyfits.getval(kw['path']+file, key, *ext))
    else:
        ret = pyfits.getval(kw['path']+filename, key, *ext)

    return ret

def write_exptime(filename, coadds='coadds'):
    """Read 'itime' and 'coadds' from a specified file, and write into it
       an "exptime' header keyword (IRAF likes it this way).  If the
       filename does not have a '.fits' extension, write_exptime will
       attempt to add one in order to find the file.

       :EXAMPLE:
          ::

            nsdata.write_exptime('blah.fits')
            """

    # 2008-06-10 15:54 IJC: Updated to accept files with or w/out .fits extension
    # 2008-06-30 17:08 IJC: Make it able to recursively use file lists
    # 2010-01-21 13:15 IJC: Added coadds keyword
    # 2010-11-28 12:26 IJC: Added ignore_missing_end call

    from pyraf import iraf as ir
    import pyfits
    import pdb

    ir.load('noao')
    ir.load('imred')
    ir.load('ccdred')


    if os.path.isfile(filename):
        
        # Check (crudely) to see if it's a FITS file.  Otherwise, assume it's a list.
        infile = open(filename)
        if infile.read(6)=='SIMPLE':   # It's a FITS file!
            infile.close()
            coa = pyfits.getval(filename,coadds, ignore_missing_end=True)
            exptime = 1.0 * coa * pyfits.getval(filename, 'itime', \
                                                    ignore_missing_end=True)
            ir.ccdhedit(filename, 'exptime', exptime)
        else:
            infile.seek(0, 0)
            for line in infile:
                write_exptime(line.strip())

    elif os.path.isfile(filename + ".fits"):
        filename = filename + ".fits"
        coa = pyfits.getval(filename, coadds, ignore_missing_end=True)
        exptime = 1.0 * coa * pyfits.getval(filename, 'itime', \
                                                ignore_missing_end=True)
        ir.ccdhedit(filename, 'exptime', exptime)
    

    return

def dark_correct(rawpre, procpre, sdark, postfix, startnum, endnum, numplaces=4):
    """Use similar parameters as nsdata.filelist() to correct frames for
    dark current.  Frames will have a 'd' appended to their filenames.

    :EXAMPLE:
       ::

         nsdata.dark_correct('/raw/targ_', '/proc/targ_', '/proc/sflat', '', 10, 20)

    2008-06-11 16:53 IJC: Created
    """

    from pyraf import iraf as ir

    raw  = wfilelist(rawpre, postfix, startnum, endnum, numplaces=numplaces, tempname='tempraw')
    proc = wfilelist(procpre, 'd'+postfix, startnum, endnum, numplaces=numplaces, tempname='tempproc')
    rawlist = wfilelist(rawpre, postfix, startnum, endnum, numplaces=numplaces)

    for ii in range(len(rawlist)):  # write "exptime" header keyword
        write_exptime(rawlist[ii])
        
    ir.ccdproc('@tempraw', output='@tempproc', ccdtype="", fixpix="no", overscan="no",trim="no",zerocor="no",darkcor="yes",flatcor="no", dark=sdark)

    os.remove('tempraw')
    os.remove('tempproc')

    return

def find_features(spec, width=4, absorption=True):
    """Scan through a spectrum and attempt to identify features based
       on the sign of derivatives.  Return the indices of the centers
       of the features.

       :EXAMPLE:
         ::

           ind = find_features(spec, width=4, absorption=True)

       2008-06-20 10:06 IJC: Created
       """

    from pylab import diff, sign, sum, zeros

    if width<2:
        print "Width too small, increasing to 2..."
        width = 2
    elif (int(width)/2.0)<>int(int(width)/2):
        print "Width not an even number, increasing to the next-largest even number."
        width = int(width)+1

    if (absorption):
        s=-1
    else:
        s=1

    w2 = width/2
    ds = diff(spec)
    sigds = sign(ds)

    numiter = len(spec)-width+1
    ind_temp = zeros(len(spec))
    ind_iter = 0

    for ii in range(numiter):
        # If the first w/2 signs are negative and the next w/2 positive, it's a feature.
        if (sum(s*sigds[ii:(ii+w2)]>0)==w2) and (sum(s*sigds[(ii+w2):(ii+width)]<0)==w2):
            ind_temp[ind_iter] = ii+w2
            ind_iter = ind_iter+1
    

    ind = ind_temp[0:ind_iter]
    return ind

def dispeval(c, olim, xlim, shift=0, function='chebyshev', fullxlim=None):
    """Evaluate an IRAF-generated dispersion function (for now, Chebyshev
       only).  Needs the input "C" matrix from ECIDENTIFY, as well as
       the limits of orders and pixels:
          ::

             w = nsdata.dispeval(C, [32,37], [1,1024])

       It can also be used in conjunction with :func:`nsdata.getdisp`:
          ::

             fn = 'ec/ecmar21s0165s'
             d  = nsdata.getdisp(fn)
             w  = nsdata.dispeval(d[0], d[1], d[2], shift=d[3])

       :SEE ALSO: :func:`nsdata.getdisp`, :func:`nsdata.interp_spec`

       :NOTE:
          May not be correct for multi-order (echelle) dispersion solutions!

       2008-06-24 11:48 IJC"""

    # 2008-06-26 10:05 IJC: Now returns w.tranpose() for proper
    #     arrangements.  Also fixed a bug in calculating the Q
    #     polynomial coefficients... was incorrect before.

    # 2012-12-13 15:34 IJMC: Now return precisely the same wavelength
    #      scale as displayed by IRAF; may not yet be correct for
    #      multi-order (echelle) dispersion solutions!

    from pylab import sort, arange, zeros, size, dot

    if function.lower()<>'chebyshev':
        print "nsdata.dispeval only works with chebyshev functions for now"
        return -1

    olim = sort(olim[0:2])
    xlim = sort(xlim[0:2])

    o = olim[1] - arange(olim[1]-olim[0]+1)
    x = arange(xlim[1]-xlim[0]+1)+xlim[0]
    on = (2.0*o-(o.max()+o.min())) / (o.max() - o.min()) 
    xn = (2.0*x-(x.max()+x.min())) / (x.max() - x.min())
    if fullxlim is not None:
        if hasattr(fullxlim, '__iter__') and len(fullxlim)>1:
            fullx = arange(fullxlim[1]-fullxlim[0]+1) + fullxlim[0]
        else:
            fullx = arange(fullxlim) +1.5
        fullxn = (2.0*fullx-(x.max()+x.min())) / (x.max() - x.min())
        
    #pdb.set_trace()
    w = zeros([len(x), len(o)], float)

    nord = c.ndim
    np = size(c,0)
    p = zeros(np, float)
    p[0] = 1

    if nord>1:
        nq = size(c,1)
        q = zeros(nq, float)
        q[0] = 1

        for i in range(len(x)):
            p[1] = xn[i]
            for k in (arange(np-2)+2):
                p[k] = 2.0*xn[i]*p[k-1] - p[k-2]

            for j in range(len(o)):
                if nord>1:
                    q[1] = on[j]
                    for k in (arange(nq-2)+2):
                        # Debugging: print i,j,k,np, nq, size(xn), size(on)
                        q[k] = 2.0*on[j]*q[k-1] - q[k-2]
                else:
                    q = array([1])

                f = dot(dot(p, c), q)
                w[i,j] = (f + shift)/o[j]
        w = w.transpose()

    else:
        if fullxlim:
            xn = fullxn.copy()
        w = c[0] + zeros(len(xn))
        for ii in range(1, c.size):
            if ii==1:
                zi = xn
                zi1 = array([1])
            else:
                zi2 = zi1.copy()
                zi1 = zi.copy()
                zi = 2*xn*zi1 - zi2
            w += c[ii]*zi
        

    return w

def getdisp(filename, mode='echelle'):
    """Read the most recent dispersion function values from an IRAF
       dispersion database file:
         ::

           D = nsdata.getdisp('ec/ecmar21s0165s')
       
       where:
          D[0] = coefficient matrix, C_mn

          D[1] = [min_order, max_order]

          D[2] = [min_pix,   max_pix]

          D[3] = pixel shift applied to the fit (units of PIXELS!)

          D[4] = type: 1=chebyshev, 2=legendre

       It is designed for use in conjunction with NSDATA.DISPEVAL:
         ::

             fn = 'ec/ecmar21s0165s'
             d  = nsdata.getdisp(fn)
             w  = nsdata.dispeval(d[0], d[1], d[2], shift=d[3])

       :SEE ALSO: :func:`nsdata.dispeval`

       2008-06-24 14:06 IJC"""
    # 2009-11-04 11:47 IJC: Print filename of missing file

    from pylab import array

    if os.path.isfile(filename):
        infile = open(filename, 'r')
    else:
        print "File %s not found!" % filename
        return -1

    raw = infile.readlines()
    infile.close()
    raw.reverse()

    ii=0
    coef = []
    while raw[ii].find('coefficients')==-1:
        if len(raw[ii].strip())>0:
            coef.append(float(raw[ii].strip()))
        ii = ii+1

    shift = 0
    jj = 0
    while (shift==0) and (jj<10):
        if raw[ii+jj].find('shift')==1:
            shift = float(  raw[ii+jj].replace('shift', '').strip()  )
        jj=jj+1
            

    func_type = int(coef[-1])  

    # Some gymnastics to get "c" matrix in the right format:
    if mode=='echelle':
        dispcoef = coef[0:xpow*opow]
        dispcoef.reverse()
        xpow = int(coef[-2])
        opow = int(coef[-3])
        xlim = [int(coef[-5]), int(coef[-6])]
        olim = [int(coef[-7]), int(coef[-8])]
        c = array(dispcoef).reshape(opow,xpow).transpose()
    else:
        opow = [coef[-2]]
        xlim = [coef[-3], coef[-4]]
        c = array(coef[::-1][4:])
        olim = [0, 0]

    return [c, olim, xlim, shift, func_type]

def isfits(filename):
    """ Test (CRUDELY!) whether a file is FITS format or not.
        ::

          result = nsdata.isfits(filename)

        :Inputs:
           filename  -- a string representing a filename

        :Outputs:
           result -- 0 if not a FITS file; 
                     1 if the explicitly passed filename is a FITS file; 
                     2 if conditions for 1 are not met, but 
                        'filename.fits' is a FITS file.

        If 'filename' does not exist, also checks filename+'.fits'.
        """
    # 2008-07-01 11:15 IJC: Created @ UCLA

    if os.path.isfile(filename):
        existence = True
        isFITSfile = 1
    elif os.path.isfile(filename+'.fits'):
        filename = filename + '.fits'
        existence = True
        isFITSfile = 2
    else:
        existence = False
        isFITSfile = 0

    if existence:
        infile = open(filename)
        if (not infile.read(6)=='SIMPLE'):
            isFITSfile = 0
        infile.close()

    return isFITSfile

def wl_grid(w, dispersion, method='log', mkhdr=False, verbose=False):
    """Generate a linear or logarithmic wavelength map with constant dispersion.
        ::

          w_new = wl_grid(w_old, dispersion, method='log')

        :INPUTS:
           wl_old -- array; the original wavelength map.  If composed
                     of N echelle orders of length M, should be shape
                     N x M
           dispersion -- maximum desired dispersion (wavelength units
                         per pixel).  

        :OPTIONAL_INPUTS:
           method  -- str; 'log' (DEFAULT) or 'linear' 
           mkhdr   -- bool; False (DEFAULT) or True.  If True, output
                        FITS header keys for each of the N echelle orders
           verbose -- bool; False (DEFAULT) or True.

        :OUTPUTS:
           w_new   -- the new wavelength map, with constant dispersion
           hdr     -- list of dict, each containing FITS wavelength headers
           

        :EXAMPLE1:
           ::

               w = array([1, 1.1, 1.2])
               w2, h = wl_grid(w, 0.05, method='log', mkhdr=True)
           
           Which gives: w2 = CRVAL1 * CDELT1**X, where X={0,...,M-1}

        :EXAMPLE2:
           ::

               w = array([1.00, 1.10, 1.21])
               w2, h = wl_grid(w, 0.05, method='linear', mkhdr=True)
           
           Which gives: w2 = CRVAL1 + CDELT1*X, where X={0,...,M-1}

        :SEE ALSO: :func:`interp_spec`

    """
# 2008-07-02 10:05 IJC: TBD: Make it work for linear wavelength spacing.
# 2008-12-01 13:35 IJC:  Linear wavelength spacing, plus FITS headers.
    
    w = array(w).copy()
    if len(w.shape)==1:
        w = w.reshape((1, w.shape[0]))
    n_order = w.shape[0]
    if verbose:  print "n_order>>" + str(n_order)

    if method=='linear':
        n_w = int(  (  (w[:,-1] - w[:,0]) / dispersion + 1 ).max() )
        if verbose:  print "n_w>>" + str(n_w)
        grid = meshgrid(range(n_w), w[:,0])
        w_interp = grid[1] + dispersion*grid[0]

    elif method=='log':
        c = dispersion/w.max() + 1.0   # logarithmic scaling constant
        if verbose:  print "c>>" + str(c)
        n_w = int(  (  log10(w[:,-1] / w[:,0]) / log10(c)  ).max()  )   + 1
        grid = meshgrid(range(n_w), w[:,0])
        w_interp = grid[1] * c**grid[0]
    else: 
        w_interp = wl_grid(w, dispersion, method='log')
    
    if verbose:  print "w_interp.shape>>" + str(w_interp.shape)

    # Make FITS header keywords
    if mkhdr:
        headers = []
        for ii in range(n_order):
            if method=='log':
                keys = dict(CTYPE1='log', CRPIX1=1, CRVAL1=w_interp[ii,0],
                            CDELT1=c)
            elif method=='linear':
                keys = dict(CTYPE1='linear', CRPIX1=1, CRVAL1=w_interp[ii,0],
                            CDELT1=dispersion)
            headers.append(keys)
        return w_interp, headers
    else:
        return w_interp

def interp_spec(filename, w, w_interp, suffix='int', k=1, badval=0, clobber=True, verbose=False):
    """ Reinterpolate a spectrum from a FITS file with a given wavelength
        calibration into a given wavelength grid.
          ::

             result = nsdata.interp_spec(filename, wl, wl_new, k=1, badval=0, clobber=True)

        :Inputs:
           filename      -- a FITS file, a list of FITS files, or a file list
                              ('.fits' optional)

           wavelengths   -- a wavelength map of the FITS files to be reinterpolated

           wl_new        -- the new wavelength map (preferably from nsdata.wl_grid)

        :Output:
           result        -- 0 if something went wrong; 1 otherwise.

        :Example:
           d_cor  = nsdata.getdisp('ec/ecscifile')

           w_old  = nsdata.dispeval(d[0], d[1], d[2], shift=d[3])

           w_new  = nsdata.wl_grid(w_old, 0.075)

           result = nsdata.interp_spec('mar21s0161s', w_old, w_new)

        :SEE ALSO: :func:`getdisp`, :func:`dispeval`, :func:`wl_grid`

    """
    # 2008-07-01 11:27 IJC: Written @ UCLA. 
    # 2008-07-22 10:14 IJC: Made it work for list input
    # 2009-07-10 10:55 IJC: Set pyfits.writeto's 'output_verify' to 'ignore'
    import pyfits
    from pylab import find

    badResult  = 0
    goodResult = 1

    # ---- Check File ---
    if filename.__class__==list:
        returnlist = []
        for element in filename: 
            returnlist.append( interp_spec(element, w, w_interp, suffix=suffix, 
                                           badval=badval, clobber=clobber) )
        return returnlist

    if (not os.path.isfile(filename)) and (not os.path.isfile(filename+'.fits')):
        return badResult

    if (not isfits(filename)):   # assume it's a file list
        infile = open(filename)
        for line in infile:
            thisfile = line.strip()
            if os.path.isfile(thisfile):
                interp_spec(thisfile, w, dispersion, suffix=suffix, badval=badval, clobber=clobber)
            elif os.path.isfile(thisfile+'.fits'):
                interp_spec(thisfile+'.fits', w, dispersion, suffix=suffix, badval=badval, clobber=clobber)
            
    else:   # must be a FITS file!
        
        if isfits(filename)==2:
            filename = filename + '.fits'

        # ---- Get data; check it ---
        irafspectrum = pyfits.getdata(  filename)
        irafheader   = pyfits.getheader(filename)

        spec_dim = rank(irafspectrum)
        if spec_dim<3:
            n_bands = 1
            if spec_dim<2:
                n_ap = 1
            else:
                n_ap = size(irafspectrum, 0)
        else:
            n_bands = size(irafspectrum, 0)
            n_ap    = size(irafspectrum, 1)
        if verbose:
            print "n_bands>>" + str(n_bands)
            print "n_ap>>" + str(n_ap)
            print "w_interp.shape>>" + str(w_interp.shape)

        n_w = w.size/len(w)
        s_interp = zeros([n_bands+1, n_ap, size(w_interp,1)])
        s_interp[n_bands,:,:] = w_interp
        if verbose:
            print "irafspectrum.shape>>" + str(irafspectrum.shape)
            print "s_interp.shape>>" + str(s_interp.shape)
        # ---- Reshape data; interpolate it ---
        irafspectrum = irafspectrum.reshape(n_bands, n_ap, size(irafspectrum)/n_bands/n_ap)
        for i_band in range(n_bands):
            for i_ap in range(n_ap):
                spec = irafspectrum[i_band,i_ap,:]
                spline = interpolate.UnivariateSpline(w[i_ap,:], spec, s=0.0, k=k)
                s_interp[i_band,i_ap,:] = spline(w_interp[i_ap,:])
                
                out_of_range = find(w_interp[i_ap,:]>w[i_ap,:].max())
                s_interp[i_band,i_ap,out_of_range] = badval
                
        # ---  Write data to new FITS file. ---
        irafheader.update('BANDID'+str(n_bands+1), 'lambda - reinterpolated wavelengths (IJC)')
        irafheader.update('WLINTERP', 'Reinterpolated! (IJC)')
        irafheader.update('BADVAL', badval)
        pyfits.writeto(filename.replace('.fits','')+suffix+'.fits', s_interp, header=irafheader, clobber=clobber, output_verify='ignore')


    return goodResult

def wspectext(inputspec):
    """Write a FITS-file spectrum to a column-format ASCII file.
       ::

           nsdata.wspectext('filelist')   ## Not yet working
           nsdata.wspectext('target.fits')
           nsdata.wspectext(specdata)     ## Not yet working

       :Inputs:

         If 3D: The first band is the data; the next-to-last is the
                 error; the last is the wavelength

       Outputs:  TBW

       2008-07-07 10:44 IJC
    """
    # 2008-07-07 10:44 IJC: Created for use with XTELLCOR; still work to do.
    # 2008-07-18 14:23 IJC: Works for FITS files; 

    import pyfits

    if inputspec.__class__==str:  # it's a filename or file list
        filetype = isfits(inputspec)
        filename = inputspec + '.fits'*(filetype-1)
        if filetype==0:   # Filelist: recursively call each line of list
            if os.path.isfile(filename):
                infile = open(filename)
                for line in infile:
                    thisfile = line.strip()
                    wspectext(thisfile)
            else:
                warn('Filelist not found!')

        elif (filetype==1) or (filetype==2):  # FITS File
            # Initialize things
            filename_out = filename.replace('.fits','')+'.dat'

            rawspec = pyfits.getdata(  filename)
            rawhead = pyfits.getheader(filename)

            naxis = rawhead['NAXIS']
            dims  = []
            for ii in range(naxis):
                dims.append(rawhead['NAXIS'+str(ii+1)])

            if dims.__len__==2:
                rawspec.reshape(1, dims[1], dims[0])
            elif dims.__len__==1:
                rawspec.reshape(1, 1, dims[0])

            # Write a file.  The first band is the data; the
            #  next-to-last is the error; the last is the wavelength

            file = open(filename_out, 'w')
            for i_ap in range(dims[1]):
                for i_line in range(dims[0]):
                    outstr = ('%.8e' % rawspec[-1, i_ap, i_line] + '  ' +
                              '%.8e' % rawspec[1,  i_ap, i_line] + '  ' +
                              '%.8e' % rawspec[-2, i_ap, i_line] + '\n')
                    file.write(outstr)

            file.close()
            

        else:
            warn("FAIL: Unexpected error")
        

    else:   # input is a spectrum
        warn("FAIL: Raw spectrum input not yet working.")

    return

def getirafap(filename, filelist=False, verbose=False):
    """ Get various parameters from an IRAF aperture file, or a list thereof
          ::

            ap = getirafap(filename, filelist=False)
            ap = getirafap(listname, filelist=True )

    :Input: filename / listname -- a string
           filelistlist        -- set to True if input is a list of files

    :Output: ap -- IRAF 'aperture' object (or a list of these objects)
    """
    # 2008-07-18 16:33 IJC: Created for files & lists

    # Subfunctions:
    def apapply(ap, filestr, keyword):
                    
        if filestr.find(keyword)<>-1:
            filestr = filestr.replace(keyword, '').strip()
            exec( 'ap.' + keyword + '.append(map(float, filestr.split()))' )
        return ap 

    # Initialize & check inputs
    ap = aperture()
    
    keywords = ['center\t', 'low\t', 'high\t']

    
    if filelist:
        infile = open(filename)
        ap = []
        for line in infile:
            if verbose:  print line.strip()
            ap.append(getirafap(line.strip()))
        infile.close()

    elif filename.__class__==list:
        ap = []
        for file in filename:
            if verbose:  print file
            ap.append(getirafap(file))

    else:
        # Read file
        if (not os.path.isfile(str(filename))):
            warn('File ' + str(filename) + ' not found!')
            return []
        ap.filename = filename
        infile = open(filename)
        raw = infile.readlines()
        infile.close()

    # Parse values for each keyword
        for line in raw:
            temp = line.strip()
            for element in keywords:
                apapply(ap, temp, element)

    return ap

def collapse_objlist(object_list, keyword, suffix='', doarray = False):
    """ Given a LIST of identical objects, collapse the objects to form a
        list of only a single keyword.
         ::

           aplist = nsdata.getirafap('ap1')
           cenlist = nsdata.collapse_objlist(aplist, 'keyword')

        If you want to get fancy, you can add a suffix to the keyword
        (e.g., to subscript):
         ::

           newlist = nsdata.collapse_objlist(object_list, 'keyword', suffix='[0]')
        """
    # 2008-07-21 10:47 IJC: Created.
    
    # allaps  = ns.getirafap(aperture, filelist=True)
    # temp    = ns.collapse_objlist(allaps, 'center', suffix='[:,1]', doarray=True)
    # centers = array(temp).ravel().reshape(len(temp), 6)
    # py.plot(x, centers - py.mean(centers))

    if type(object_list)<>list:
        warn('Input object not a list!')
        return []

    object_keywords = dir(object_list)
    found_keyword = False
    for val in dir(object_list[0]):
        found_keyword = (found_keyword or (val==keyword))
    if (not found_keyword):
        warn('Specified keyword "' + keyword + '" not found!')
        return []


    collapsed_list = []
    mainclass = object_list[0].__class__
    for element in object_list:
        if element.__class__<>mainclass:
            warn('Input list is not homogeneous!')
            return []
        elif doarray:
            collapsed_list.append(eval('array(element.' + keyword + ')' + suffix) )
        else:
            collapsed_list.append(eval('element.' + keyword + suffix))

    return collapsed_list

def file2list(filelist, prefix='', suffix=''):
    """Convert an IRAF-like data file list of images to a Python list.
       ::

         newlist = nsdata.file2list('filelist', prefix='', suffix='')

       A prefix or suffix can also be appended to each filename.  This
       is useful if, e.g., you need an explicit file extension identifier.


       :SEE ALSO:  :func:`wfilelist`, :func:`filelist`, :func:`sfilelist`
       """

    #2008-07-25 16:30 IJC: Created

    newlist = []
    
    if filelist.__class__<>str: 
        warn('Input to file2list must be a string identifying a filename.')
    elif (not os.path.isfile(filelist)):
        warn('Input file "' + filelist + '" not found!')

    else:   # Everything is good.
        infile = open(filelist)

        for line in infile:
            filename = prefix + line.strip() + suffix
            if len(filename)>0:
                newlist.append(filename)

        infile.close()

    return newlist

def imshow(data, x=[], y=[], aspect='auto', interpolation='nearest', cmap=None, vmin=[], vmax=[]):
    """ Version of pylab's IMSHOW with my own defaults:
    ::

      imshow(data, aspect='auto', interpolation='nearest', cmap=cm.gray, vmin=[], vmax=[])

    Other IMSHOW options are default, but a new one exists: 
          x=  and y=  let you set the axes values by passing in the x and y coordinates."""
    #2008-07-25 18:30 IJC: Created to save a little bit of time and do axes.

    from pylab import arange, cm, imshow

    if cmap==None:
        cmap = cm.gray

    def getextent(data, x, y):
        """ Gets the extent of the data for plotting.  Subfunc of IMSHOW."""
        dsh = data.shape

        if len(x)==0:
            x = arange(dsh[1])
        if len(y)==0:
            y = arange(dsh[0])

        dx = 1.0* (x.max() - x.min()) / (len(x) - 1)
        xextent = [x.min() - dx/2.0, x.max() + dx/2.0]
        xextent = [x[0] - dx/2.0, x[-1] + dx/2.0]

        dy = 1.0* (y.max() - y.min()) / (len(y) - 1)
        yextent = [y.max() + dy/2.0, y.min() - dy/2.0]
        yextent = [y[-1] + dy/2.0, y[0] - dy/2.0]

        extent = xextent + yextent
        
        return extent

    def getclim(data, vmin, vmax):
        if vmin.__class__==list:
            vmin = data.min()
        if vmax.__class__==list:
            vmax = data.max()
        return [vmin, vmax]
    
    #------------- Start the actual routine -------------

    extent = getextent(data, x,y)
    clim   = getclim(data, vmin, vmax)
    imshow(data, aspect=aspect, interpolation=interpolation, cmap=cmap, 
              vmin=clim[0], vmax=clim[1], extent=extent)


def subdata(data, op='median', axis=None, returndata=False):
    """Take the mean/median along a specified direction and subtract it
       from the rest of the data.

       :EXAMPLE:
         ::

          p   = [[1,2,3], [4,5,7]]
          q1 =  nsdata.subdata(p, op='mean',   axis=0)
          q2 =  nsdata.subdata(p, op='median', axis=1)

       :Gives:
             q1:   [[-1.5, -1.5, -2], [1.5, 1.5, 2]]

             q2:   [[  -1,    0,  1], [ -1,   0, 2]]

       :KEYWORDS:
            *op:   operation to perform; either 'median' (DEFAULT) or 'mean'

            *axis: axis along which to perform 'op'; if None (DEFAULT),
                   'op' is performed on the entire data set as a whole.

            *returndata: Whether to also return the data series by which the
                          division was performed.  DEFAULT is False.

        :REQUIREMENTS:
            :doc:`analysis`

        :SEE ALSO:  :func:`divdata`
     """

    # 2008-07-25 19:42 IJC: Created, and proud of it.
    import analysis as an

    data = array(data)
    dsh = list(data.shape)

    if op=='median':
        chunk = an.amedian(data, axis=axis)
    elif op=='mean':
        chunk = mean(data, axis=axis)

    if axis<>None:
        dsh[axis] = 1
        chunk = chunk.reshape(dsh)

    newdata = data - chunk

    if returndata:
        return (newdata, chunk)
    else:
        return newdata


def divdata(data, op='median', axis=None, badval=nan, returndata=False):
    """Take the mean/median along a specified direction and divide
       the rest of the data by it.

    :EXAMPLE:
      ::

        p   = [[1,2,3], [4,5,7]]
        q1 =  nsdata.divdata(p, op='mean',   axis=0)
        q2 =  nsdata.divdata(p, op='median', axis=1)

    :Gives:
          q1:   [[-1.5, -1.5, -2], [1.5, 1.5, 2]]

          q2:   [[  -1,    0,  1], [ -1,   0, 2]]

    :KEYWORDS:
        *op:         operation to perform; either 'median' (DEFAULT),
                       'mean', or 'none' (i.e., divide by 1)

        *axis:       axis along which to perform 'op'; if None (DEFAULT),
                      'op' is performed on the entire data set as a whole.

        *badval:     value to replace any residual nan/inf values with.
                       DEFAULT is nan.  Makes two passes, pre- and
                       post-division.

        *returndata: Whether to also return the data series by which the
                      division was performed.  DEFAULT is False.

    :REQUIREMENTS:
        :doc:`analysis`

    :SEE ALSO:  :func:`subdata`
     """
    
    # 2008-07-25 19:42 IJC: Created, and proud of it.
    # 2009-10-19 22:18 IJC: There's always room for improvement.  Uses
    #                       fixval now.

    from pylab import find
    from analysis import fixval, amedian

    data = array(data)
    dsh = list(data.shape)
    nsh = list(dsh)
    nsh[axis] = 1

    #print dsh

    fixval(data, badval)

    if op=='median':
        chunk = amedian(data, axis=axis)
    elif op=='mean':
        chunk = mean(data, axis=axis)
    elif op=='none':
        chunk = ones(nsh,float)

    if axis<>None:
        chunk = chunk.reshape(nsh)

    newdata = 1.0 * data / chunk

    fixval(newdata, badval)

    if returndata:
        return (newdata, chunk)
    else:
        return newdata

def repval(data, badval, newval):
    """ Replace all occurrences of one value with another value.  This
        handles nan and inf as well.

        :EXAMPLE:             To set all 'nan' values to zero, just type:
          ::
 
            C = repval(data, nan, 0)

        """
    #2008-07-29 10:25 IJC: Created
    
    from pylab import find

    data = array(data)
    dsh = data.shape
    data = data.ravel()
    
    if isfinite(badval):
        ind = find(data==badval)
    elif isnan(badval):
        ind = find(isnan(data))
    elif isinf(badval):
        ind = find(isinf(data))
    else:
        warn('Value to replace is neither finite, nan, nor inf... error!!')
        return array([])

    data[ind] = newval
    data = data.reshape(dsh)

    return data

def mederr(data, ntrials=1000, mode='median'):
    """ Return the median or mode, and the 68.3% error on that
    quantity, using bootstrapping."""
    # 2008-07-29 16:54 IJC: Created
    # 2009-08-27 12:43 IJC: Added 'mode="mean"' option.
    data = array(data.ravel())

    if mode=='median':
        medianval = median(data)
    elif mode=="mean":
        medianval = mean(data)

    nd    = len(data)
    sigcutoff = [int(0.1586*ntrials+0.5), int(0.8414*ntrials+0.5)]  # contains 68.27% of the data -- one standard deviation.

    meds = zeros(ntrials)
    for ii in range(ntrials):
        ind = (random.rand(nd)*nd - 0.5).round()
        newdata = data[ind.tolist()]
        if mode=='median':
            meds[ii] = median(newdata)
        elif mode=='mean':
            meds[ii] = mean(newdata)

    meds.sort()
    lowsig = meds[sigcutoff[0]]
    hisig  = meds[sigcutoff[1]]

    medianstd = (lowsig - hisig)/2.0

    return (medianval, medianstd)

def quadd(arg1, arg2):
    """ Add two arguments in quadrature.
      ::

        print quadd(0.1, 0.2)       -->  0.2236
        print sqrt(0.1**2 + 0.2**2) -->  0.2236
    """
    # 2008-07-30 19:53 IJC: Created
    arg1 = array(arg1, subok=True, copy=True)
    arg2 = array(arg2, subok=True, copy=True)

    return sqrt(arg1**2 + arg2**2)

def gd2jd(datestr):
    """ Convert a string Gregorian date into a Julian date using Pylab.
        If no time is given (i.e., only a date), then noon is assumed.
        Timezones can be given, but UTC is assumed otherwise.

       :EXAMPLES:
          ::

            print gd2jd('Aug 11 2007')   #---------------> 2454324.5
            print gd2jd('Aug 11 2007, 12:00 PST')  #-----> 2454324.29167
            print gd2jd('12:00 PM, January 1, 2000')  #--> 2451545.0

       :REQUIREMENTS: :doc:`matplotlib`

       :SEE ALSO: :func:`jd2gd`
       """
# 2008-08-26 14:03 IJC: Created        
# 2010-12-08 13:00 IJC: Removed "+ 3442850" from num2julian call
# 2011-05-19 11:37 IJMC: Put the factor back in for error-catching...
    
    import matplotlib.dates as dates
    
    if datestr.__class__==str:
        d = dates.datestr2num(datestr)
        jd = dates.num2julian(d) 
        if jd<0:
            jd = dates.num2julian(d + 3442850)
            print "You are probably using an old version of Matplotlib..."
    else:
        jd = []

    return jd

def jd2gd(juldat):
    """ Convert a numerial Julian date into a Gregorian date using Pylab.
        Timezone returned will be UTC.

       :EXAMPLES:
         ::

          print jd2gd(2454324.5)  #--> 2007-08-12 00:00:00
          print jd2gd(2451545)    #--> 2000-01-01 12:00:00

       :SEE ALSO: :func:`gd2jd`"""
    # 2008-08-26 14:03 IJC: Created    
    # 2011-01-22 16:24 IJC: Removed arbitrary (?) subtraction of 3442850 from 'd'
    # 2011-10-21 14:11 IJMC: Put it back in, but with MPL version-checking.
    import matplotlib.dates as dates
    from matplotlib import __version__

    if __version__ < '1.0.0':
        print "You are probably using an old version of Matplotlib..."
        d = dates.julian2num(juldat - 3442850)
    else:
        d = dates.julian2num(juldat)
    gd = dates.num2date(d )

    return gd

def fix_quadnoise(*args, **kw):
    """Fix the 8-row coherent patterns in each quadrant of NIRSPEC using
       linear least-squares after removing outliers.

    :INPUTS:
       file -- a filename or list of FITS files.  The file suffix
             '.fits' is appended if the file cannot be found.  If this
             parameter begins with the '@' ('at') symbol, it is
             interpreted as an IRAF file list.

    :OPTIONAL_INPUTS:
       prefix  -- prefix to add to the fixed files.  

       clobber -- overwrite existing files

       verbose -- boolean flag for more output printed to the screen

    :OUTPUTS:
       none

    :EXAMPLE:
       ::

       fix_quadnoise(file, verbose=False, clobber=True)
       fix_quadnoise('mar21s0165.fits')
    """
    # 2009-11-03 09:33 IJC: Created anew; fit to background levels.
    # 2010-09-07 11:04 IJC: Substantially revised to combat occasional 1e9 pedestals
    import pyfits
    from phot import estbg
    from numpy import zeros, arange, round, hstack, dot, prod
    from numpy.linalg import pinv
    from analysis import removeoutliers

    sigma = 7
    nrow = 8

    # -------- initialize inputs --------
    if len(args)==0:
        print "No filename given.  Exiting."
        return
    else:
        file = args[0]

    defaults = dict(prefix='qfix', verbose=False, clobber=False)

    if len(kw)==0:
        verbose = False
        clobber = False
    else:
        for key in defaults:
            if (not kw.has_key(key)):
                kw[key] = defaults[key]
    
    verbose  = bool(kw['verbose'])
    prefix   = str(kw['prefix'])
    clobber = bool(kw['clobber'])

    if file.__class__==list:
        for element in file:
            if verbose: print "Python file list, file: " + str(element)
            fix_quadnoise(element, row=row, quadrant=quadrant, verbose=verbose, prefix=prefix)
        return
    elif file[0]=='@':
        f = open(file[1::])
        for line in f:
            if verbose: print "IRAF-type file list, file: " + str(line.strip())
            fix_quadnoise(line.strip(), row=row, quadrant=quadrant, verbose=verbose, prefix=prefix)
        f.close()
        return

    if (not os.path.isfile(file)):
        file = file + '.fits'
    try:
        data = pyfits.getdata(file)
        hdr  = pyfits.getheader(file)
    except:
        print "PYFITS Could not read from file '" + file + "' -- exiting."
        return

    x,y = meshgrid(arange(data.shape[1]),arange(data.shape[0]))
    basis = zeros([data.shape[0]/2,data.shape[1]/2,8],float)
    for ii in range(0,nrow):
        basis[ii:(512+ii):nrow,:,ii] += 1.0

    for ii in range(4):
        if ii==0:
            ind = ((x<512) * (y<512))
        elif ii==1:
            ind = ((x<512) * (y>=512))
        elif ii==2:
            ind = ((x>=512) * (y<512))
        elif ii==3:
            ind = ((x>=512) * (y>=512))

        rowmeans = zeros(nrow, float)
        d = data[ind].reshape(512,512)
        for jj in range(nrow):
            goodvalues, goodind = removeoutliers(d[jj::8,:], sigma, retind=True)
            rowmeans[jj] = goodvalues.mean()

        new = basis * (rowmeans.reshape(1,1,nrow) - median(rowmeans))
        data[ind] -= new.sum(2).ravel()

        print "row-means for quad %i are>>" % ii, rowmeans
        if verbose:
            print "row-means for quad %i are>>" % ii, rowmeans
            #from pylab import *



    if verbose: print "file>>" + str(file)
    if verbose: print "os.path.split(file)[0]>>" + os.path.split(file)[0]
    initpath = os.path.split(file)[0]
    if len(initpath)==0:  
        initpath = '.'
    outfn = initpath + os.sep + prefix + os.path.split(file)[1]
    if verbose: print "Writing file..." + outfn
    pyfits.writeto(outfn, data, hdr, clobber=clobber, output_verify='ignore')

    return

def preprocess(*args, **kw):
    """
    Basic processing of NIRSPEC data: set JD and HJD fields, fix
    bad-row 'quadnoise', clean cosmic rays, flatten and remove bad
    pixels.

    :INPUTS:
         input = 

         output = 

         Input and output must both be specified, and must both be
         different files.

    :OPTIONAL_INPUTS:
         qfix    = True
           qpref = '' -- string to preface all quad-fixed frames

         flat    = None -- flat field for iraf.ccdproc

         mask    = None -- bad pixel mask for iraf.ccdproc

         clobber = False -- overwrite file if input and output files are same

         verbose = False

         cleanec = False -- run nsdata.cleanec (a "poor-man's iraf.cosmicrays")
             cthreshold = 300 -- threshold for nsdata.cleanec

             csigma = 20 -- sigma threshold for nsdata.cleanec

             cwindow = 25 -- window size for nsdata.cleanec

         cleancr = False -- run iraf.cosmicrays
             rthreshold = 300 -- threshold for iraf.cosmicrays

             rratio = 5 -- fluxratio threshold for iraf.cosmicrays

         Any inputs set to 'None' disable that part of the processing.

    :EXAMPLE:
      ::

         ns.preprocess('@'+rawcal, '@'+proccal, qfix=, qpref=\
                        flat=, mask=, clobber=, verbose=)


    Note that although the syntax is similar to Pyraf tasks, only
    Python Booleans or their equivalent should be used for flags
    (i.e., use True or False instead of 'yes' or 'no'
    """
    # 2008-11-26 19:23 IJC: Created on an Amtrak train
    # 2009-11-03 11:15 IJC: Updated to use new LLS fix_quadnoise
    # 2010-08-27 13:22 IJC: Added 'fluxratio' parameter to cosmicrays
    # 2010-08-28 07:28 IJC: Added 'cleanec' function as well
    # 2010-09-06 14:48 IJC: Better integrated cleanec, and added an option thereto
    # 2010-09-08 08:56 IJC: Re-added iraf.cosmicrays, with an option

    from pyraf import iraf as ir
    ir.load('crutil')
    # Parse inputs:
    if len(args)<2:
        print "No output images specified!  Exiting..."
        return
    else:
        input = args[0]
        output = args[1]
        if input==output:
            print "Input and output files must not be the same! Exiting..."
            return

    defaults = dict(qfix=False, qpref='', flat=None, mask=None, \
                        cleanec=False, clobber=False, verbose=False, \
                        cthreshold=300, cwindow=25, csigma=20, \
                        cleancr=False, rthreshold=300, rratio=5)

    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]

    verbose = bool(kw['verbose'])
    doflat = kw['flat']<>None
    dobfix = kw['mask']<>None
    clobber = bool(kw['clobber'])

    if verbose: print "kw>>" + str(kw)
    if verbose: print "Keywords are:  "; print kw

    # Check whether inputs are lists, filelists, or files:

    if input.__class__<>output.__class__:
        print "Files or file lists must be of same type.  Exiting..."
        return
    elif input.__class__==list:
        for ii in range(len(input)):
            if verbose: print "File list, file:  " + input[ii]
            preprocess(input[ii], output[ii], **kw)
        return
    elif (input.__class__==str) and (input[0]=='@'):
        fin  = open(input[ 1::])
        fout = open(output[1::])

        ## Don't flatten, bad-pixel correct, or cleanec the first time through:

        for line in fin:
            if verbose: print "IRAF-style file list, file:  " + line.strip()
            preprocess(line.strip(), fout.readline().strip(), **kw)
        fin.close()
        fout.close()
        return
 
    # Begin processing tasks
    if kw['clobber'] and input<>output:
        ir.imdelete(output)

    ir.imcopy(input, output)

    # I don't like IRAF.setjd because it gave me incorrect values.
    # Instead, use my own setjd from astrolib.py:
    setjd(output, date='date-obs', time='UTC', jd='JD', hjd='HJD')
    
    if kw['qfix']:
        fix_quadnoise(output, prefix=kw['qpref'],clobber=clobber)
        ir.hedit(output, 'quadnois', \
                     'NIRSPEC bad row fixed by nsdata.fix_quadnoise')
    
    if doflat:

        ir.ccdproc(output, ccdtype="", fixpix="no", overscan="no",
                   trim="no", zerocor="no", darkcor="no", flatcor=doflat, 
                   flat=kw['flat'], fixfile=None,
                   minreplace=0.25, interactive="no")

    if dobfix: 
        # Make an extra bad-pixel mask from any _negative_ values, and
        # combine it with the global bad-pixel mask; necessary because
        # some negative-valued pixels manage to make it through the
        # cleaning pipeline.
        indiv_mask = output + 'imask.fits'
        cutoffmask(output, clobber=True, cutoff=[0, Inf], writeto=indiv_mask)
        ir.imcalc(kw['mask'] + "," + indiv_mask, indiv_mask, "im1||im2")
        ir.ccdproc(output, ccdtype="", fixpix=dobfix, overscan="no",
                   trim="no", zerocor="no", darkcor="no", flatcor="no", 
                   flat=None, fixfile=indiv_mask,
                   minreplace=0.25, interactive="no")

    if kw['cleancr']:
        ir.cosmicrays(output, output, threshold=kw['rthreshold'], fluxratio=kw['rratio'], \
                          npasses=5, interactive='no')
    if kw['cleanec']:
        cleanec(output, output, npasses=1, verbose=verbose, threshold=kw['cthreshold'], \
                    nsigma=kw['csigma'], window=kw['cwindow'], clobber=True)

    if verbose: print "Successfully processed '" + input + \
            "' into '" + output + "'"

    return

def setjd(filename, **kw):
    """Set JD and HJD fields in a FITS file, based on the UTC date and
    time header keywords.

    :INPUTS:
       filename : str

    :OPTIONS:
       ra=None : Right Ascension in decimal degrees.  If None, use "RA" header key

       dec=None : Declination in decimal degrees.  If None, use "DEC" header key 

    :EXAMPLE:
      ::

         dict(date='date-obs', time='UTC', epoch='equinox', jd='JD', hjd='hjd', \
            verbose=False, ra=None, dec=None)


    This does not take your position on Earth into account, so it must
    be less accurate than ~0.1 sec.
    """
    # 2010-09-07 15:14 IJC: Created

    import astrolib
    import pyfits

    defaults = dict(date='date-obs', time='UTC', epoch='equinox', jd='JD', hjd='hjd', \
                        verbose=False, ra=None, dec=None)
    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]
    
    verbose = kw['verbose']
    ra = kw['ra']
    dec = kw['dec']

    if filename.__class__==list:
        for ii in range(len(filename)):
            if verbose: print "File list, file:  " + filename[ii]
            setjd(filename[ii], **kw)
        return

    elif (filename.__class__==str) and (filename[0]=='@'):
        fin  = open(filename[ 1::])
        for line in fin:
            setjd(line.strip(), **kw)

    else: 
        try:
            hdulist = pyfits.open(filename)
        except IOError:
            filename += '.fits'
            hdulist = pyfits.open(filename)

        hdr = hdulist[0].header
        datetime = hdr[kw['date']] + ' ' + hdr[kw['time']]
        if ra is None:
            ra = hdr['ra']
        if dec is None:
            dec = hdr['dec']
        epoch = hdr[kw['epoch']]
        if epoch <> 2000:
            print "Epoch must be 2000!  Exiting..."
            return -1

        if verbose:
            print 'datetime>>', datetime
            print 'ra, dec>>', ra, dec
        jd = gd2jd(datetime)
        hjd = astrolib.helio_jd(jd - 2400000., ra, dec) + 2400000.

        hdr.update(kw['jd'], jd)
        hdr.update(kw['hjd'], hjd)
        if verbose: print 'jd, hjd, dt>>', jd, hjd, hjd-jd

        hdulist.writeto(filename, output_verify='ignore', clobber=True)
    return

def linespec(loc, ew, win, **kw):
    """
    Create a delta-function line spectrum based on a wavelength grid
    and a list of line locations and equivalent widths.

    :INPUTS:
       loc -- location of lines in the emission frame of reference

       ew  -- equivalent widths of lines, in units of wavelength grid.
               Positive values are emission lines.

       w_in -- wavelength grid in the emission frame, with values
              monotonically increasing (best if it is linearly spaced)

       All inputs should be lists or one-dimensional arrays of scalars

    :OPTIONAL_INPUTS:
       cont=None -- set continuum values in the emission frame;

       nearest=False  -- if True, use full pixels instead of partial

       verbose=False  -- if True, print out various messages

    :OUTPUTS:
      s  -- delta-function line spectrum, with a continuum level of zero
    
    :EXAMPLE: (NEEDS TO BE UPDATED!):
       ::

          w   = [2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7]
          loc = [2.1, 2.35, 2.62]
          ew  = [0.1, .02, .01]
          s = linespec(loc, ew, w)
          print s  #  --->  [0, 1, 0, 0.1, 0.1, 0, 0.08, 0.02]

    :NOTE:  This may give incorrect results for saturated lines.
    """
    # 2008-12-05 13:31 IJC: Created
    # 2008-12-10 13:30 IJC: Added continuum option, reworked code some.
    # 2008-12-12 12:33 IJC: Removed RV option

    from pylab import find

    # Check inputs
    loc = array(loc).copy().ravel()
    ew  = array(ew ).copy().ravel()
    win   = array(win  ).copy().ravel()

    defaults = dict(cont=None, nearest=False, verbose=False)
    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]
    verbose = bool(kw['verbose'])
    nearest = bool(kw['nearest'])
    contset = kw['cont']<>None

    if contset:
        cont = array(kw['cont']).copy()
        if len(cont)<>len(win):
            print "Wavelength grid and continuum must have the same length!"
            return -1
    else:
        cont = ones(win.shape)

    nlines = len(loc)
    if nlines <> len(ew):
        if verbose:  print "len(loc)>>" + str(len(loc))
        if verbose:  print "len(ew)>>" + str(len(ew))
        print "Line locations and equivalent widths must have same length!"
        return -1

    #Only use lines in the proper wavelength range
    nlineinit = len(loc)
    lind = (loc>=win.min()) * (loc<=win.max())
    loc = loc[lind]
    ew  =  ew[lind]
    nlines = len(loc)

    s = cont.copy()
    d = diff(win).mean()

    if verbose:  print "s>>" + str(s)

    for ii in range(nlines):
        lineloc = loc[ii]
        lineew  = ew[ii]
        index = (win<lineloc).sum() - 1
        if nearest:
            s[index+1] = s[index]-cont[index]*lineew/d
        elif index==len(win):
            s[index] = s[index] - cont[index]*lineew/d
        else:
            s[index] = s[index] - lineew*cont[index]* \
                (win[index+1] - lineloc)/d/d
            s[index+1] = s[index+1] - lineew*cont[index+1] * \
                (lineloc - win[index])/d/d
        
        if verbose:  
            print "(lineloc, lineew)>>" + str((lineloc, lineew))
            print "(index, d)>>" + str((index,d))

    if verbose:
        print "(nlineinit, nline)>>" + str((nlineinit, nlines))
    return s
            
def readfile(fn, cols=None):
    """ Read data from an space- or tab-delimited ASCII file."""

    # 2008-12-12 11:42 IJC: created

    f = open(fn, 'r')
    raw = f.readlines()
    f.close()

    if cols==None:
        dat = array([map(float, line.split()) for line in raw])
    else:
        dat = array([map(float, line.split()[cols]) for line in raw])

    return dat
    
def initobs(date, **kw):
    """Initialize variables for Nirspec data analysis.

    :INPUT:
       date -- a string of type YYYYMMMDD (e.g., 2008mar21 or 2008jun15a)

    :OPTIONAL_INPUT:
       remote=False -- changes processed-data directory as specified
                       within.

       interp=True -- whether to load wavelength-interpolated
                      (upsampled) spectra (if True) or the raw
                      1024-pixel spectra (if False)

    :OUTPUT:
       a tuple containing the following values, in order:
         planet -- name of planet for use in analysis.planet()

         datalist -- list of data file numbers to analyse

         _proc  -- processed data directory

         wavefilename -- filename of the wavelength solution

         starmodelfilename -- path and filename of stellar model

         planetmodelfilename -- path and filename of planet model

         aplist -- a list of the IRAF aperture filenames (for use with
                   nsdata.getirafap)

         telluric -- the FITS file spectrum of the telluric/A0V spectrum

         n_aperture -- number of echelle apertures for this setup

         filter -- NIRSPEC filter used

         prefix -- filename prefix

         calnod -- whether calibrator stars were nodded.

         rowfix -- list of four lists; which rows to fix in each of
                   four quadrants (see FIX_QUADNOISE).
"""

    # 2009-07-09 16:58 IJC: Added HD 189733b set (2008jun15b)
    # 2009-07-31 15:17 IJC: Added multi-chop runs (e.g. 2008jun15)
    # 2009-08-05 12:17 IJC: Now multi-chop runs only return dates
    # 2010-08-15 14:57 IJC: Flagged a bad file in 2008jun15b
    # 2010-08-24 15:42 IJC: Added 2010aug16 dataset
    # 2010-09-03 23:45 IJC: Added 2010sep04 dataset
    # 2010-11-10 17:02 IJC: Added 2008jul12 A/B nodding datasets
    # 2011-08-10 18:06 IJMC: Added 2011aug05 SpeX/GJ 1214b run.

    import analysis as an

    defaults = dict(remote=False, interp=True)
    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]


    if date=='2010apr22':
        planet = 'WASP-12 b' 
        prefix = 'NS.20100422.'
        framelist = [19873, 20165, 20421, 20688, 21022, 21298, 21549, 21814, 22148, 22423, 22674, 22943, 23242, 23508, 23759, 24030] #, 
        flatlist = [25878] #[10108, 10178, 10216, 10254, 10292, 10330, 10367, 10405, 10443, 10481, 10519, 10556, 10594, 10640]
        darklist = [10808, 10852, 10890, 10928, 10966]
        callist = [25127, 25212, 25285, 25369, 25499, 25583, 25657, 25744]
        datadir = 'wasp12_2010apr22/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_wasp12_2010apr22_ns_20100422'
        n_aperture = 6
        filter = 'K'
        calnod = False


    elif date=='2008mar21':
        planet = '55 Cnc e' #'55cnce'
        prefix = 'mar21s'
        #framelist = range(29,43) + range(45,50) + range(51, 161)
        framelist = range(61,161)
        flatlist = range(166,186)
        darklist = range(186,206)
        callist = range(161,166)
        datadir = '55cnc_2008mar21/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_55cnc_2008mar21_mar21s'
        n_aperture = 6
        filter = 'K'
        calnod = False
        rowfix = ['a']
    elif date=='2008jun15':
        suffix = 'bcdefg'
        ret = [date+suf for suf in suffix]
        return ret
    elif date=='2008jul12':
        suffix = 'abcd'
        ret = [date+suf for suf in suffix]
        return ret
    elif date=='2009jul29':
        suffix = 'abcdefgh'
        ret = [date+suf for suf in suffix]
        return ret

    elif date=='2011aug05':
        planet = 'GJ 1214 b' 
        prefix = 'spectra'
        framelist = range(251,462)
        flatlist = range(602,701)
        darklist = range(502,601)
        callist = range(462,502)
        datadir = 'gj1214_2011aug05/'
        ap_suffix = 'undetermined' #'database/ap_Users_ianc_proj_pcsa_data_proc_gj1214_2010aug16_aug16s'
        n_aperture = 6
        filter = 'SXD'
        calnod = False

    elif date=='2010aug16':
        planet = 'GJ 1214 b' #hd209458b'
        prefix = 'aug16s'
        framelist = range(366,491) + range(492, 549) #range(366,549)
        flatlist = range(260,310)
        darklist = range(310,360)
        callist = range(360,366)  # range(550,565)
        datadir = 'gj1214_2010aug16/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_gj1214_2010aug16_aug16s'
        n_aperture = 7
        filter = 'K'
        calnod = False

    elif date=='2010aug16b':
        planet = 'GJ 1214 b' 
        prefix = 'aug16s'
        framelist = range(500,549)
        flatlist = range(260,310)
        darklist = range(310,360)
        callist =  range(550,565)
        datadir = 'gj1214_2010aug16/'
        ap_suffix = 'ec/ap_Users_ianc_proj_pcsa_data_proc_gj1214_2010aug16_aug16s'
        n_aperture = 7
        filter = 'K'
        calnod = False

    elif date=='2010sep04':
        planet = 'GJ 1214 b' 
        prefix = 'sep04s'
        framelist = range(144,152) + range(154,182) + range(185,198) + range(199,271)
        framelist = range(185,198) + range(199,271)
        flatlist = range(52, 102)
        darklist = range(103, 123) + range(124, 144)
        callist = range(272, 280)
        datadir = 'gj1214_2010sep04/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_gj1214_2010sep04_sep04s'
        n_aperture = 8
        filter = 'H'
        calnod = False

    elif date=='2009jul29a':
        planet = 'HD 209458 b' #hd209458b'
        prefix = 'jul29s'
        #framelist = range(29,43) + range(45,50) + range(51, 161)
        framelist = range(194,207)
        flatlist = range(616,716)
        darklist = range(516,616)
        callist = range(188,194)  
        datadir = 'hd209458_2009jul29/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd209458_2009jul29_jul29s'
        n_aperture = 6
        filter = 'K'
        calnod = False
    elif date=='2009jul29b':
        planet = 'HD 209458 b' #'hd209458b'
        prefix = 'jul29s'
        framelist = range(213,235)
        flatlist = range(616,716)
        darklist = range(516,616)
        callist = range(207,213)  
        datadir = 'hd209458_2009jul29/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd209458_2009jul29_jul29s'
        n_aperture = 6
        filter = 'K'
        calnod = False
    elif date=='2009jul29c':
        planet = 'HD 209458 b'#'hd209458b'
        prefix = 'jul29s'
        framelist = range(247,289)
        flatlist = range(616,716)
        darklist = range(516,616)
        callist = range(236, 246)  
        datadir = 'hd209458_2009jul29/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd209458_2009jul29_jul29s'
        n_aperture = 6
        filter = 'K'
        calnod = False
    elif date=='2009jul29d':
        planet = 'HD 209458 b'#'hd209458b'
        prefix = 'jul29s'
        framelist = range(297,354)
        flatlist = range(616,716)
        darklist = range(516,616)
        callist = range(289,297)  
        datadir = 'hd209458_2009jul29/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd209458_2009jul29_jul29s'
        n_aperture = 6
        filter = 'K'
        calnod = False
    elif date=='2009jul29e':
        planet = 'HD 209458 b'#'hd209458b'
        prefix = 'jul29s'
        framelist = range(368,414)
        flatlist = range(616,716)
        darklist = range(516,616)
        callist = range(362,368)  
        datadir = 'hd209458_2009jul29/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd209458_2009jul29_jul29s'
        n_aperture = 6
        filter = 'K'
        calnod = False
    elif date=='2009jul29f':
        planet = 'HD 209458 b'#'hd209458b'
        prefix = 'jul29s'
        framelist = range(422,464)  # 464 has bad headers!
        flatlist = range(616,716)
        darklist = range(516,616)
        callist = range(414,422)  
        datadir = 'hd209458_2009jul29/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd209458_2009jul29_jul29s'
        n_aperture = 6
        filter = 'K'
        calnod = False
    elif date=='2009jul29g':
        planet = 'HD 209458 b'#'hd209458b'
        prefix = 'jul29s'
        framelist = range(473,506)
        flatlist = range(616,716)
        darklist = range(516,616)
        callist = range(465,473)  
        datadir = 'hd209458_2009jul29/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd209458_2009jul29_jul29s'
        n_aperture = 6
        filter = 'K'
        calnod = False
    elif date=='2009jul29h':
        planet = 'HD 209458 b'#'hd209458b'
        prefix = 'jul29s'
        framelist = []
        flatlist = range(616,716)
        darklist = range(516,616)
        callist = range(506,516)  
        datadir = 'hd209458_2009jul29/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd209458_2009jul29_jul29s'
        n_aperture = 6
        filter = 'K'
        calnod = False
    elif date=='2008jun15b':
        planet = 'HD 189733 b' #'hd189733b'
        prefix = 'jun15s'
        framelist = range(315,332)+range(333,335)
        flatlist = range(1,101)
        darklist = range(101,201)
        callist = [313,314]
        datadir = 'hd189733_2008jun15/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jun15_jun15s'
        n_aperture = 6
        filter = 'K'
        calnod = True
    elif date=='2008jun15c':
        base = initobs('2008jun15b')
        planet = base[0]
        prefix = base[16]
        framelist = range(337,352)
        flatlist = range(1,101)
        darklist = range(101,201)
        callist = [335,336]
        datadir = 'hd189733_2008jun15/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jun15_jun15s'
        n_aperture = base[14]
        filter = base[15]
        calnod = base[17]
    elif date=='2008jun15d':
        base = initobs('2008jun15b')
        planet = base[0]
        prefix = base[16]
        framelist = range(354,369)
        flatlist = range(1,101)
        darklist = range(101,201)
        callist = [352,353]
        datadir = 'hd189733_2008jun15/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jun15_jun15s'
        n_aperture = base[14]
        filter = base[15]
        calnod = base[17]
    elif date=='2008jun15e':
        base = initobs('2008jun15b')
        planet = base[0]
        prefix = base[16]
        framelist = range(371,386) + range(387,401)
        flatlist = range(1,101)
        darklist = range(101,201)
        callist = [369,370]
        datadir = 'hd189733_2008jun15/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jun15_jun15s'
        n_aperture = base[14]
        filter = base[15]
        calnod = base[17]
    elif date=='2008jun15f':
        base = initobs('2008jun15b')
        planet = base[0]
        prefix = base[16]
        framelist = range(403,443)
        flatlist = range(1,101)
        darklist = range(101,201)
        callist = [401,402]
        datadir = 'hd189733_2008jun15/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jun15_jun15s'
        n_aperture = base[14]
        filter = base[15]
        calnod = base[17]
    elif date=='2008jun15g':
        base = initobs('2008jun15b')
        planet = base[0]
        prefix = base[16]
        framelist = range(445,544)
        flatlist = range(1,101)
        darklist = range(101,201)
        callist = [443,444]
        datadir = 'hd189733_2008jun15/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jun15_jun15s'
        n_aperture = base[14]
        filter = base[15]
        calnod = base[17]
    elif date=='2008jun15a':
        planet = 'tau Boo b' #'tauboob'
        prefix = 'jun15s'
        framelist = range(205,245) + [246,247]+range(249,253) + [256] + range(263,277) #+ range(287,297) +[303,304]
        flatlist = range(1,101)
        darklist = range(101,201)
        callist = [305]
        datadir = 'tauboo_2008jun15/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_tauboo_2008jun15_jun15s'
        n_aperture = 6
        filter = 'K'
        calnod = False

    elif date=='2008jul12a':
        planet = 'HD 189733 b' #'hd189733b'
        prefix = 'jul12s'
        framelist = range(131, 162)
        flatlist = range(261,461)
        darklist = range(464, 564)
        callist = range(161,165)
        datadir = 'hd189733_2008jul12/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jul12_jul12s'
        n_aperture = 5
        filter = 'L'
        calnod = True

    elif date=='2008jul12b':
        base = initobs('2008jul12a')
        planet = base[0]
        prefix = base[16]
        framelist = range(166,202)
        flatlist = range(261,461)
        darklist = range(464,564)
        callist = range(202,206)
        datadir = 'hd189733_2008jul12/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jul12_jul12s'
        n_aperture = base[14]
        filter = base[15]
        calnod = base[17]

    elif date=='2008jul12c':
        base = initobs('2008jul12a')
        planet = base[0]
        prefix = base[16]
        framelist = range(206,230)
        flatlist = range(261,461)
        darklist = range(464,564)
        callist = [229,233]
        datadir = 'hd189733_2008jul12/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jul12_jul12s'
        n_aperture = base[14]
        filter = base[15]
        calnod = base[17]

    elif date=='2008jul12d':
        base = initobs('2008jul12a')
        planet = base[0]
        prefix = base[16]
        framelist = range(234,257)
        flatlist = range(261,461)
        darklist = range(464,564)
        callist = [257,261]
        datadir = 'hd189733_2008jul12/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jul12_jul12s'
        n_aperture = base[14]
        filter = base[15]
        calnod = base[17]

    elif date=='2008jul12A':
        base = initobs('2008jul12a')
        planet = base[0]
        prefix = base[16]
        framelist = [131, 132, 135, 136, 139, 140, 143, 144, 147, 148, 151, 152, 155, \
                         156, 159, 160, 167, 168, 171, 172, 175, 176, 179, 180, 183, 184, \
                         187, 188, 191, 192, 195, 196, 199, 200, 207, 208, 211, 212, 215, \
                         216, 219, 220, 223, 224, 227, 228, 235, 236, 239, 240, 243, 244, \
                         250, 251, 254, 255]
        flatlist = range(261,461)
        darklist = range(464,564)
        callist = [229,233]
        datadir = 'hd189733_2008jul12/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jul12_jul12s'
        n_aperture = base[14]
        filter = base[15]
        calnod = base[17]

    elif date=='2008jul12B':
        base = initobs('2008jul12a')
        planet = base[0]
        prefix = base[16]
        framelist = [133, 134, 137, 138, 141, 142, 145, 146, 149, 150, 153, 154, 157, \
                         158, 161, 166, 169, 170, 173, 174, 177, 178, 181, 182, 185, 186, \
                         189, 190, 193, 194, 197, 198, 201, 206, 209, 210, 213, 214, 217, \
                         218, 221, 222, 225, 226, 229, 234, 237, 238, 241, 242, 245, 246, \
                         247, 248, 249, 252, 253, 256]
        flatlist = range(261,461)
        darklist = range(464,564)
        callist = [229,233]
        datadir = 'hd189733_2008jul12/'
        ap_suffix = 'database/ap_Users_ianc_proj_pcsa_data_proc_hd189733_2008jul12_jul12s'
        n_aperture = base[14]
        filter = base[15]
        calnod = base[17]





    else:
        print date+' not found!!!!'
        return []

    if len(date)==9:
        meancal = 'avgcal'
    else:
        meancal = 'avgcal'+date[9::]

    date = date[0:9]

    if kw['interp']:
        postfix = 'sinttel.fits'
    else:
        postfix = 's.fits'

    wavefilename = 'winterp.fits'

    if kw['remote']:
        _proc = "/Users/ianc/atwork/proj/pcsa/data/proc/" + datadir
        _raw = "/Users/ianc/atwork/proj/pcsa/data/raw/" + date + "/spec/"
        _model = "/Users/ianc/atwork/proj/pcsa/data/model/"
    else:
        _proc = "/Users/ianc/proj/pcsa/data/proc/" + datadir
        _raw = "/Users/ianc/proj/pcsa/data/raw/" + date + "/spec/"
        _model = "/Users/ianc/proj/pcsa/data/model/"
        
    if planet=='55 Cnc e': #'55cnce':
        starmodelfilename = _model + \
            'lte5243_4.33_+0.25_55Cnc.hires.7.fits'
        planetmodelfilename = _model + \
            'lte0125-3.5.rainout.HD75732e.redist_0.50.hires.7.fits'
    elif planet=='tau Boo b': #'tauboob':
        starmodelfilename = _model + \
            'lte6309_4.30_0.28_tau_boo.hires.7.fits'
        planetmodelfilename = _model + \
            'lte0125-3.5.rainout.HD120136b.redist=0.5.hires_id.7.fits'
    elif planet=='HD 189733 b': #'hd189733b':
        starmodelfilename = _model + \
            'XXXX_Lband_189model'
        planetmodelfilename = _model + \
            'XXXX_Lband_189bmodel'
    elif planet=='HD 209458 b': #'hd209458b':
        print "WARNING: HD209458 DOES NOT HAVE CORRECT STELLAR MODEL!"
        starmodelfilename = _model + \
            'lte5243_4.33_+0.25_55Cnc.hires.7.fits'
        print "WARNING: HD209458 DOES NOT HAVE CORRECT PLANETARY MODEL!"
        planetmodelfilename = _model + \
            'lte0125-3.5.rainout.HD75732e.redist_0.50.hires.7.fits'
    elif planet=='GJ 1214 b':
        if filter=='K':
            starmodelfilename = ''
            planetmodelfilename = _model + \
                'lte0125-3.00.rainout_irrad.GJ1214b_redist=0.25.aces_rlam.fits'
            print "WARNING: GJ1214b K-band DOES NOT HAVE CORRECT STELLAR MODEL!"
            planetmodelfilename = _model + \
                'noplanetmodel'
        elif filter=='H':
            starmodelfilename = _model + 'lte3000_5.00-0.0.GJ1214_hires_Hband.7.fits'
            planetmodelfilename = _model + \
                'lte0125-3.00.rainout_irrad.GJ1214b_redist=0.25.aces_rlam.fits'
            print "WARNING: GJ1214b H-band DOES HAS A _TRANSMISSION SPECTRUM_ MODEL!"
        elif filter=='SXD':
            starmodelfilename = _model + 'NOSTELLARMODEL'
            planetmodelfilename = _model + 'lte0125-3.00.rainout_irrad.GJ1214b_redist=0.25.aces_rlam.fits'
    elif planet=='WASP-12 b':
        print "WARNING: WASP12 DOES NOT HAVE CORRECT STELLAR MODEL!"
        starmodelfilename = _model + \
            'lte5243_4.33_+0.25_55Cnc.hires.7.fits'
        print "WARNING: WASP12 DOES NOT HAVE CORRECT PLANETARY MODEL!"
        planetmodelfilename = _model + \
            'lte0125-3.5.rainout.HD75732e.redist_0.50.hires.7.fits'


    else:
        stop


    datalist = filelist(prefix, postfix, framelist)

    darkfilelist = filelist(_raw+prefix, '.fits', darklist, numplaces=4)
    flatfilelist = filelist(_raw+prefix, '.fits', flatlist, numplaces=4)
    rawcalfilelist = filelist( _raw+prefix, '.fits', callist, numplaces=4)
    proccalfilelist = filelist( _proc+prefix, 'fn', callist, numplaces=4)
    rawtargfilelist = filelist( _raw+prefix, '.fits', framelist, numplaces=4)
    proctargfilelist = filelist( _proc+prefix, 'fn', framelist, numplaces=4)
    speccalfilelist = filelist(prefix, 's', callist, numplaces=4)
    spectargfilelist = filelist(prefix, 's', framelist, numplaces=4)

    aplist     = filelist(_proc+ap_suffix, 'fn', framelist)
    aplist_cal = filelist(_proc+ap_suffix, 'fn', callist)
    
    telluric_fn = _proc + prefix + meancal + 'int.fits'
    disp = 0.075

    ret = (planet, _proc, datalist, wavefilename, 
           starmodelfilename, planetmodelfilename, aplist, telluric_fn, _raw,
           darkfilelist, flatfilelist, (rawcalfilelist, proccalfilelist),
           (rawtargfilelist, proctargfilelist), (speccalfilelist, spectargfilelist), 
           n_aperture, filter, prefix, calnod, meancal, disp, aplist_cal)
    return ret

def cleanec(input, output, **kw):
    """Clean raw echelleogram files of bad pixels.

    Similar to IRAF's 'cosmicrays' or 'crmed' tasks, but this only
    looks along the spectrum's dispersion direction.  The less
    parallel your spectra are to the pixel grid, the less well this
    may work; in this case try increasing the 'window' parameter
    somewhat.

    :INPUTS:

      file : (str or numpy.array) 
           input echellogram file to be cleaned.

    :DEFAULT_OPTIONS:
    
       dict(dispaxis=0, npasses=1, mapout=None, verbose=False, 
                        threshold=300, nsigma=15, window=25, clobber=False, 
                        hdr=None)
           """
    # 2010-08-27 13:34 IJC: CReated
    # 2010-08-30 17:01 IJC: Activated 'npasses' option

    import pyfits
    from pylab import find
    from time import time

    # Check if filein is str or array; load file if the former.  Set
    # dispersion direction via transpose.
    defaults = dict(dispaxis=0, npasses=1, mapout=None, verbose=False, \
                        threshold=300, nsigma=15, window=25, clobber=False, \
                        hdr=None)

    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]

    verbose = kw['verbose']
    hdr = kw['hdr']
    if verbose: print "CLEANEC begun, input>>", input

    # Test for type of input; iterate, if necessary:
    if input.__class__<>output.__class__:
        print "Files or file lists must be of same type.  Exiting..."
        return
    elif input.__class__==list:
        for ii in range(len(input)):
            if verbose: print "File list, file:  " + input[ii]
            cleanec(input[ii], output[ii], **kw)
        return
    elif (input.__class__==str) and (input[0]=='@'):
        fin  = open(input[ 1::], 'r')
        fout = open(output[1::], 'r')
        if verbose: 
            print "using input file %s and output file %s" % (fin.name, fout.name)
        kw_bak = kw.copy()
        for infile, outfile in zip(fin, fout):
            infile = infile.strip()
            outfile = outfile.strip()
            if verbose: print "IRAF-style file list, file %s --> %s" % \
                    (infile.strip(), outfile.strip())
            cleanec(infile, outfile, **kw)
        fin.close()
        fout.close()
        return
    elif input.__class__==str:
        try:
            ec = pyfits.getdata(input)
            hdr = pyfits.getheader(input)
        except:
            print "Couldn't open file: %s -- adding .fits extension and re-trying" % input
            try:
                input += '.fits'
                if output.find('.fit')<0:
                    output += '.fits'
                ec = pyfits.getdata(input)
                hdr = pyfits.getheader(input)
            except:
                print "PYFITS could not open file: %s" % input
                return
    else:
        try:
            ec = input.copy()
        except:
            print "Input should be a list of filenames, a string, or " + \
                "a Numpy array; it appears to be none of these."

    if kw['dispaxis']<>0:
        ec = ec.transpose()

    ncol, nrow = ec.shape

    passNumber = 0
    while passNumber<kw['npasses']:
        t0, t1, t2, t25, t27, t3 = 0,0,0,0,0,0
        allrows = arange(nrow)
        for ipix in range(ncol):
            # Set indices and extract the segment to examine
            tic = time()
            minind = max(0, ipix - kw['window']/2)
            maxind = min(ncol, ipix + kw['window']/2)
            segs = ec[:, minind:maxind]
            segind = set(arange(segs.shape[1]))
            residuals = abs(segs - median(segs,1).reshape(nrow,1))
            maxres = tile(residuals.max(1),(segs.shape[1],1)).transpose()
            t0 += time()-tic

            tic = time()
            maxind = array([nonzero(test)[0][0] for test in residuals==maxres])
            goodind = map(list, [segind.difference(set([maxval])) for maxval in maxind])
            t1 += time()-tic

            tic = time()
            segstd = std([ss[gg] for ss,gg in zip(segs,goodind)], 1)
            # Measure (a) its discrepancy and (b) stdev of the remainder
            t2 += time()-tic
            tic = time()
            discrepancy = residuals[allrows, maxind]
            sigma = discrepancy/segstd
            repind =  (sigma>kw['nsigma']) * (discrepancy>kw['threshold'])
            t25 += time()-tic

            prior_value = maxind[repind] - 1
            prior_value[prior_value<0] = 0
            latter_value = maxind[repind] + 1
            latter_value[prior_value<0] = 0

            # If pixel is sufficiently discrepant, throw it out!
            #   if maximum value is _first_, return the latter value
            #   if maximum value is _last_, return the former value
            #   otherwise, return average of former & latter
            tic = time()
            firstbadval = maxind[repind]==0
            lastbadval = maxind[repind]==(segs.shape[1]-1)
            midbadval = True - (firstbadval + lastbadval)
            maxind2 = maxind[repind][midbadval]
            repind2 = find(repind)[midbadval]
    
            # Need to apply averaging separately from endpoint replacement.
            segs[repind, maxind[repind]] = \
                firstbadval * segs[repind,1] + \
                lastbadval * segs[repind,segs.shape[1]-2] 
            segs[repind2, maxind2] += \
                0.5 * (segs[repind2, maxind2-1] + segs[repind2, maxind2+1])
            t3 += time() - tic
            if verbose: 
                if ipix/100==(1.0*ipix/100.):
                    print 'finished column %i/%i' % (ipix, ncol),
                    sys.stdout.flush()

        passNumber += 1
    
        if verbose: print ('pass %i complete: times elapsed per section>>' % 
                passNumber), t0, t1, t2, t25, t3, 

    if hdr is None:
        pyfits.writeto(output, ec, clobber=kw['clobber'], \
                           output_verify='ignore')
    else:
        hdr.update('cleanec', \
                 'echelleogram cleaned (%i passes) by nsdata.cleanec' \
                       % passNumber)
        pyfits.writeto(output, ec, clobber=kw['clobber'], header=hdr, \
                           output_verify='ignore')

    if verbose: print "CLEANEC complete, output to>>", output
    return 


#ir.crmed(f12, f12.replace('fn.fits', 'fn2.fits'), ncmed=15, nlmed=1, ncsig = 25, nlsig = 10, crmask='', median='', sigma='', residual='')

def bfixpix(data, badmask, n=4, retdat=False):
    """Replace pixels flagged as nonzero in a bad-pixel mask with the
    average of their nearest four good neighboring pixels.

    :INPUTS:
      data : numpy array (two-dimensional)

      badmask : numpy array (same shape as data)

    :OPTIONAL_INPUTS:
      n : int
        number of nearby, good pixels to average over

      retdat : bool
        If True, return an array instead of replacing-in-place and do
        _not_ modify input array `data`.  This is always True if a 1D
        array is input!

    :RETURNS: 
      another numpy array (if retdat is True)
    """
    # 2010-09-02 11:40 IJC: Created
    #2012-04-05 14:12 IJMC: Added retdat option
    # 2012-04-06 18:51 IJMC: Added a kludgey way to work for 1D inputs
    # 2012-08-09 11:39 IJMC: Now the 'n' option actually works.
    

    if data.ndim==1:
        data = np.tile(data, (3,1))
        badmask = np.tile(badmask, (3,1))
        ret = bfixpix(data, badmask, n=2, retdat=True)
        return ret[1]


    nx, ny = data.shape

    badx, bady = nonzero(badmask)
    nbad = len(badx)

    if retdat:
        data = array(data, copy=True)
    
    for ii in range(nbad):
        thisloc = badx[ii], bady[ii]
        rad = 0
        numNearbyGoodPixels = 0

        while numNearbyGoodPixels<n:
            rad += 1
            xmin = max(0, badx[ii]-rad)
            xmax = min(nx, badx[ii]+rad)
            ymin = max(0, bady[ii]-rad)
            ymax = min(ny, bady[ii]+rad)
            x = arange(nx)[xmin:xmax+1]
            y = arange(ny)[ymin:ymax+1]
            yy,xx = meshgrid(y,x)
            #print ii, rad, xmin, xmax, ymin, ymax, badmask.shape
            
            rr = abs(xx + 1j*yy) * (1. - badmask[xmin:xmax+1,ymin:ymax+1])
            numNearbyGoodPixels = (rr>0).sum()
        
        closestDistances = unique(sort(rr[rr>0])[0:n])
        numDistances = len(closestDistances)
        localSum = 0.
        localDenominator = 0.
        for jj in range(numDistances):
            localSum += data[xmin:xmax+1,ymin:ymax+1][rr==closestDistances[jj]].sum()
            localDenominator += (rr==closestDistances[jj]).sum()

        #print badx[ii], bady[ii], 1.0 * localSum / localDenominator, data[xmin:xmax+1,ymin:ymax+1]
        data[badx[ii], bady[ii]] = 1.0 * localSum / localDenominator

    if retdat:
        ret = data
    else:
        ret = None

    return ret
    
def posofend(str1, str2):
    """returns the position immediately _after_ the end of the occurence
    of str2 in str1.  If str2 is not in str1, return -1.
    """
    pos = str1.find(str2)
    if pos>-1:
        ret= pos+len(str2)
    else:
        ret = -1
    return ret

def readart(filename):
    """Read a Keck .arT telemetry file.
    """
    # 2010-09-04 01:11 IJC: Created

    from spitzer import genericObject


    def process(struct, valstr, names, types, delim=','):
        """Append a line of values to the output object.
        """
        vals = valstr.split(delim)
        for val, name, type in zip(vals, names, types):
            if type=='DBL':
                exec 'struct.%s.append(float(val.strip()))' % name in locals()
            else:
                exec 'struct.%s.append(val.strip())' % name in locals()
        return

    try:
        f = open(filename, 'r')
    except:
        print "Could not open file: %s" % filename
        return -1
    
    output = genericObject()

    hdrlines = [f.readline() for ii in range(3)]
    pos00, pos01 = posofend(hdrlines[0], 'INTERVAL ='), hdrlines[0].find(',')
    pos10, pos11 = posofend(hdrlines[0], 'NO_ELEMENTS ='), hdrlines[0].find('\n')
    interval = int(hdrlines[0][pos00:pos01])
    no_elements = int(hdrlines[0][pos10:pos11])

    entry_names = [name.strip().replace(':','_').replace('.','_') for name in hdrlines[1].replace('"', '').split(',')]
    entry_types = [name.strip().replace(':','_').replace('.','_') for name in hdrlines[2].replace('"', '').split(',')]
    entry_types = [name.strip() for name in hdrlines[2].replace('"', '').split(',')]

    output.interval = interval
    output.no_elements = no_elements
    output.keys = entry_names
    output.types = entry_types
    for name in entry_names:
        exec 'output.%s = []' % name in locals()
    
    for nextline in f:
        process(output, nextline, entry_names, entry_types)

    return output

def readnstemps(filename):
    """Read the Keck telemetry file "nirspecTemps"
    """
    # 2010-09-04 17:46 IJC: Created

    from spitzer import genericObject

    try:
        f = open(filename, 'r')
    except:
        print "Could not open file: %s" % filename
        return -1

    output = genericObject()
    line0 = f.readline()
    pos0 = posofend(line0, 'GMT')
    output.time = [line0[0:pos0]]
    vals = map(float, line0[pos0::].strip().split(' '))
    nvals = len(vals)
    valnames = ['t%i' % ii for ii in range(nvals)]
    for val, name in zip(vals, valnames):
        exec 'output.%s = [val]' % name in locals()

    for line in f:
        vals = map(float, line[pos0::].strip().split(' '))        
        output.time.append(line[0:pos0])
        for val, name in zip(vals, valnames):
            exec 'output.%s.append(val)' % name in locals()

    return output

def readMagiqlog(filename, dat=['guidestats']):
    """
    Read specified types of data from the Keck Magiq telemetry log.

    filename : str -- logfile name

    dat : list of str -- types of data to return
         guidestats -- centroid x/y, fwhm, star & skyflux
         """
    # 2010-09-05 12:07 IJC: Created

    def processGuidestats(struct, string):
        pos_t0, pos_t1 = 0, string.find('[LegacyCam]')
        pos_x0, pos_x1 = posofend(string, 'Centroid x='), string.find('y=')
        pos_y0, pos_y1 = posofend(string, 'y='), string.find('fwhm=')
        pos_f0, pos_f1 = posofend(string, 'fwhm='), string.find('star=')
        pos_st0, pos_st1 = posofend(string, 'star='), string.find('sky=')
        pos_sk0, pos_sk1 = posofend(string, 'sky='), string.find('stats=')
        
        struct.HSTtime.append(string[pos_t0:pos_t1])
        struct.cenx.append(float(  string[pos_x0:pos_x1]))
        struct.ceny.append(float(  string[pos_y0:pos_y1]))
        struct.fwhm.append(float(  string[pos_f0:pos_f1]))
        struct.star.append(float(  string[pos_st0:pos_st1]))
        struct.sky.append(float(   string[pos_sk0:pos_sk1]))

        return
    
    from spitzer import genericObject

    try:
        f = open(filename, 'r')
    except:
        print "Could not open file: %s" % filename
        return -1

    output = genericObject()

    for type in dat:
        if type=='guidestats':
            output.guide = genericObject()
            output.guide.HSTtime = []
            output.guide.cenx = []
            output.guide.ceny = []
            output.guide.fwhm = []
            output.guide.star = []
            output.guide.sky = []
        else:
            print "Unknown type of Magiq log data: %s" % type
            
    for line in f:
        for type in dat:
            if type=='guidestats':
                if line.find('CamImage')<0 or line.find('Centroid')<0 or line.find('fwhm')<0:
                    pass
                else:
                    processGuidestats(output.guide, line)
            else:
                pass

    return output

def readarm(filename):
    """Read a Keck .arM telemetry file.

    So far, tested only with envAut.arM
    """
    # 2010-09-04 01:11 IJC: Created

    from spitzer import genericObject


    def process(struct, line, delim=','):
        """Append a line of values to the output object.
        """
        vals = line.split(delim)
        keyname = vals[2].strip().replace(':','_').replace('.','_').replace('---','')
        endofkey = posofend(line, vals[2])
        endofdelim = line.find(delim, endofkey)+1
        keyval = line[endofdelim::].strip()

        try:
            keyval = float(keyval)
        except:
            pass
        
        if keyval=='***':
            pass
        else:
            #struct.HSTtime.append(vals[0].strip())
            #struct.HSTdate.append(vals[1].strip())
            HSTdatetime = '%s %s' % (vals[0].strip(), vals[1].strip())
            try:
                exec 'struct.%s.append(keyval)' % keyname in locals()
                exec "struct.%sHST.append(HSTdatetime)" % keyname in locals()
            except:
                exec 'struct.%s = [keyval]' % keyname in locals()
                exec "struct.%sHST = [HSTdatetime]" % keyname in locals()
        
        return

    try:
        f = open(filename, 'r')
    except:
        print "Could not open file: %s" % filename
        return -1
    
    output = genericObject()
    hdrlines = [f.readline() for ii in range(2)]
    pos00, pos01 = posofend(hdrlines[0], 'INTERVAL ='), hdrlines[0].find(',')
    pos10, pos11 = posofend(hdrlines[0], 'NO_ELEMENTS ='), hdrlines[0].find('\n')
    interval = int(hdrlines[0][pos00:pos01])
    no_elements = int(hdrlines[0][pos10:pos11])

    entry_names = [name.strip().replace(':','_').replace('.','_') for name in hdrlines[1].replace('"', '').split(',')]

    output.interval = interval
    output.no_elements = no_elements
    output.HSTtime = []
    output.HSTdate = []
    output.HSTdatetime = []

    for line in f:
        process(output, line)
        

    return output
        
    #    entry_names = [name.strip().replace(':','_').replace('.','_') for name in hdrlines[1].replace('"', '').split(',')]
    #entry_types = [name.strip().replace(':','_').replace('.','_') for name in hdrlines[2].replace('"', '').split(',')]
    #entry_types = [name.strip() for name in hdrlines[2].replace('"', '').split(',')]

def envMet(filename, tz=-10, planet=None, date=None):
    """Read the envMet.arT Keck telemetry file.

    :INPUT:
       filename -- str.

    :OPTIONAL_INPUTS:
       tz: offset (in hours) from GMT (HST = -10)

       planet: an.planet object to compute orbital phase & HJD (if desired)

       date: observation run date; grab time tags from this to use for
              averaging over (not yet implemented!!!)
    """
    # 2010-09-07 13:02 IJC: Created
    import pyfits
    import astrolib
    
    met = readart(filename)

    if met<>-1:
        met.HSTdatetime = ['%s %s' % (ddate, time) for ddate,time in zip(met.HSTdate, met.HSTtime)]

        met.jd = [gd2jd(datetime)-tz/24. for datetime in met.HSTdatetime]

        if planet is not None:
            met.hjd = [astrolib.helio_jd(jjj - 2.4e6, planet.ra, planet.dec) + 2.4e6 for jjj in met.jd]
            met.phase = planet.phase(met.hjd)

        if date is not None:
            obs = initobs(date)
            rawfiles = obs[12][0]
            jd_start = []
            exptime = []
            try:
                hdr = pyfits.getheader(file)
                jd_start.append(hdr['jd'])
                exptime.append(hdr['exptime'])
            except:
                jd_start = []

            # Loop over every observation timestap, averaging over each header keyword
            met.bin = met
            met.bin.jd = jd_start
#            for ii, jd0 in enumerate(jd_start):
#                thisindex = 
            
            

    return met


def darkbpmap(filelist, clipsigma=5, sigma=25, writeto=None, clobber=False, verbose=False):
    """Use dark frames to construct a bad pixel map based on unusually
       variable pixels.

    :INPUT:
        filelist: str, list, or 3D numpy array -- dark frame filenames or data
            str -- IRAF-style file list (beginning with '@' symbol)
            list -- Python-style list of strs
            numpy array -- 3D (L*M*N) stack of L dark frames

    :OPTIONS:
        clipsigma : scalar -- significance threshold for removing transient
                   (cosmic-ray-like) events

        sigma : scalar -- significance threshold for greater-than-average
                pixel variability.

        writeto : str -- filename to write output to.  If None, returns the array.

    :RETURNS:
      if writeto is None:
          a 2D boolean numpy array: True for bad pixels, False for other pixels.
      else:
          returns None
    """
    # 2010-09-08 09:56 IJC: Created
    
    import pyfits
    from analysis import stdr, meanr

    darkstack = []
    if filelist.__class__==list:
        for line in filelist:
            if verbose: print "File list, file:  " + line
            try:
                darkstack.append(pyfits.getdata(line))
            except:
                darkstack.append(pyfits.getdata(line+'.fits'))
    elif (filelist.__class__==str) and (filelist[0]=='@'):
        fin  = open(input[ 1::], 'r')
        if verbose: 
            print "using input file %s" % (fin.name)
        for line in fin:
            line = line.strip()
            if verbose: print "IRAF-style file list, file %s" % line
            try:
                darkstack.append(pyfits.getdata(line))
            except:
                darkstack.append(pyfits.getdata(line+'.fits'))
    else:
        darkstack = filelist

    darkstack = array(darkstack)

    darkstd = stdr(darkstack, axis=0, nsigma=clipsigma, verbose=verbose-1)
    edark = stdr(darkstd.ravel(), nsigma=clipsigma, verbose=verbose-1)
    mdark = meanr(darkstd.ravel(), nsigma=clipsigma, verbose=verbose-1)

    badpixelmap = (abs(darkstd-mdark) / edark) > sigma
    if writeto is None:
        ret = badpixelmap
    else:
        ret = None
        pyfits.writeto(writeto, badpixelmap.astype(int), clobber=clobber)

    return ret
    
def cutoffmask(filename, cutoff=[0, Inf], writeto=None, clobber=True):
    """Create a simple mask from a FITS frame or array, based on whether
       its values are above or below specified cutoff values.

    :INPUT:
        filename: str or numpy array --  frame filename or data
            str -- filename of FITS frame to load
            numpy array -- data frame

    :OPTIONS:
        cutoff : list -- of form [lower, higher]; values below 'lower' or
                above 'higher' will be flagged as bad.

        writeto : str -- filename to write output to.  If None, returns the array.

    :RETURNS:
      if writeto is None:
          a 2D boolean numpy array: True for bad pixels, False for other pixels.
      else:
          returns None
    """
    # 2010-09-08 10:36 IJC: Created
    
    import pyfits
    from analysis import stdr, meanr

    if (filename.__class__==str):
        try:
            data = pyfits.getdata(filename)
        except:
            data = pyfits.getdata(filename+'.fits')
    else:
        data = filename

    badpixelmap = (data < cutoff[0]) + (data>cutoff[1])

    if writeto is None:
        ret = badpixelmap
    else:
        ret = None
        pyfits.writeto(writeto, badpixelmap.astype(int), clobber=clobber)

    return ret

def subab(afiles, bfiles, outfiles, clobber=False):
    """
    Take raw A & B nod files from a set of observations and subtract
    them; output the difference images.

    :INPUTS:
       afiles -- list of A-position filenames
       
       bfiles -- list of B-position filenames

       outfiles -- list of output filename

    :NOTE:
       All input lists will be truncated to the length of the shortest
       input list.
       
       For now, the output file has the header of th A-file of the
       pair; this is not optimal and should be fixed to record both A
       & B headers!

    :EXAMPLE:
      import nsdata as ns
      inita = ns.initobs('2008jul12A')
      initb = ns.initobs('2008jul12B')
      afiles = inita[12][0][1:-4]
      bfiles = initb[12][0][:-4]
      outfiles = [(os.path.split(bfn)[0] + '/%s-%s' % \
            (os.path.split(afn)[1], os.path.split(bfn)[1])).replace('.fits', '') + \
            '.fits' for afn, bfn in zip(afiles, bfiles)]
      ns.subab(afiles, bfiles, outfiles)

    """
    # 2010-11-29 08:43 IJC: Created

    import pyfits
    import os
    
    for afn, bfn, outfn in zip(afiles, bfiles, outfiles):
        if not os.path.isfile(afn):
            print "file %s not found, skipping." % afn
            readfiles = False
        elif not os.path.isfile(bfn):
            print "file %s not found, skipping." % bfn
            readfiles = False
        else:  # File was found
            try:
                adata = pyfits.getdata(afn, ignore_missing_end=True)
                bdata = pyfits.getdata(bfn, ignore_missing_end=True)
                ahdr = pyfits.getheader(afn, ignore_missing_end=True)
                bhdr = pyfits.getheader(bfn, ignore_missing_end=True)
                readfiles = True
            except:
                readfiles = False
                print "Could not read data or header from FITS files " + \
                    "(either %s or %s)" % (afn, bfn)

        if readfiles:
            # Test for non-standard NIRSPEC keywords
            if ahdr.has_key('gain.spe'):
                ahdr.rename_key('gain.spe', 'gain_spe')
            if ahdr.has_key('freq.spe'):
                ahdr.rename_key('freq.spe', 'freq_spe')
            if bhdr.has_key('gain.spe'):
                bhdr.rename_key('gain.spe', 'gain_spe')
            if bhdr.has_key('freq.spe'):
                bhdr.rename_key('freq.spe', 'freq_spe')


            diff = adata - bdata
            pyfits.writeto(outfn, diff, ahdr, ignore_missing_end=True, clobber=clobber)
                
    return

def nAir_old(vaclam, T=288.15, P=101325.):
    """Return the index of refraction of air at a given wavelength.

    :INPUTS: 

       vaclam -- scalar or Numpy array -- Vacuum wavelength (in
                 microns) at which to calculate n
    
       T -- scalar -- temperature in Kelvin
       
       P -- scalar -- pressure in Pascals.

    :NOTES:

       This assumes a dry atmosphere with 0.03% CO2 by volume ("standard air")

    :REFERENCE: 

        81st CRC Handbook (c.f. Edlen 1966)
    """
    # 2011-03-10 14:12 IJC: Created

    print "This function is outdated; please use nAir() instead."

    sigma2 = 1. / vaclam**2

    nm1_stp = 1e-8 * (8342.13 + 2406030. / (130. - sigma2) + 15997. / (38.9 - sigma2))

    tp_factor =(P * (1. + P * (61.3 - (T - 273.15)) * 1e-10)) / (96095.4 * (1. + 0.003661 * (T - 273.15) ) )

    n = 1. + nm1_stp * tp_factor

    #print (nm1_stp * 1e8), (n - 1) * 1e8
    return n

def nAir(vaclam, T=293.15, P=1e5, fco2=0.0004, pph2o=0.):
    """Return the index of refraction of air at a given wavelength.

    :INPUTS: 

       vaclam: scalar or Numpy array
              Vacuum wavelength (in microns) at which to calculate n
    
       T : scalar
           temperature in Kelvin
       
       P : scalar
           pressure in Pascals

       fc02 : scalar
           carbon dioxide content, as a fraction of the total atmosphere

       pph2o : scalar
           water vapor partial pressure, in Pascals

    :REFERENCE: 
specfi        Boensch and Potulski, 1998 Metrologia 35 133
    """
    # 2011-10-07 15:14 IJMC: Created
    # 2012-12-05 20:47 IJMC: Explicitly added check for 'None' option inputs.

    if T is None:
        T = 293.15
    if P is None:
        P = 1e5
    if fco2 is None:
        fco2 = 0.0004
    if pph2o is None:
        pph2o = 0.0

    sigma = 1./vaclam
    sigma2 = sigma * sigma

    # (Eq. 6a)
    nm1_drystp =  1e-8 * (8091.37 + 2333983. / (130. - sigma2) + 15518. / (38.9 - sigma2))

    # Effect of CO2 (Eq. 7):
    nm1_dryco2 = nm1_drystp * (1. + 0.5327 * (fco2 - 0.0004))

    # Effect of temperature and pressure (Eq. 8):
    nm1_dry = ((nm1_dryco2 * P) * 1.0727933e-5) * \
        (1. + 1e-8 * (-2.10233 - 0.009876 * T) * P) / (0.0036610 * T)

    # Effect of H2O (Eq. 9):
    try:
        n = 1. + (nm1_dry - pph2o * (3.8020 - 0.0384 * sigma2) * 1e-10).astype(float64)
    except:
        n = 1. + (nm1_dry - pph2o * (3.8020 - 0.0384 * sigma2) * 1e-10)

    return n


