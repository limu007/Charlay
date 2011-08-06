# -*- coding: utf-8 -*-

'''provides simple or more complicated tools for general 1/2-D data analysis:
    width of distributions fitted with gaussian
    gaps in data
    kernel density estimators
    chi2 mapping
'''

def join(lists,meth=None):
    '''
    merging more lists in one
    '''
    rep=[]
    for l in lists:
        if l==None: continue
        if meth!=None:
            rep+=getattr(l,meth)
        else: rep+=l
    return rep

def uniq(ilist,prop):
    olist=[]
    ohash=[]
    for o in ilist:
        if hasattr(o,prop):
            v=getattr(o,prop)
            if not v in ohash:
                ohash.append(v)
                olist.append(o)
    return olist

def get_width(x,y,bord=5,rep=1,sig=1.,mod='gauss'):
    '''
    fits a gaussian to a histogram
    '''
    from numpy import log,polyfit,polyval,exp
    vals=y[bord:-bord]
    if mod=='gauss':vals=log(vals)
    res2=polyfit(x[bord:-bord],vals,2)
    wid=sqrt(abs(1/res2[0]/2))
    if rep==1: return wid
    mid=-res2[1]/res2[0]/2
    if rep>3: 
        vals=polyval(res2,x[bord:-bord])
        if mod=='gauss':vals=exp(vals)
        return mid,wid,y[bord:-bord]-vals
    elif rep==3: return mid,wid,sum(y[(x<mid-sig*wid)|(x>mid+sig*wid)])
    return mid,wid

from numpy import array,arange,ones,zeros,sum,sqrt,abs,float64,int32,bool,size
from math import log

loud=0

global a,b,c
def get_sigma(data,ndiv=30,min_count=3,min_bins=7,rep=1,nclip=2.,nsig=2.):
    from numpy import histogram,argmax
    global a,b,c
    vmax=data.mean()
    vstd=data.std()
    if nclip>0:# clipped
        sata=data[abs(data-vmax)<nclip*vstd]
        vmax=sata.mean()
        vstd=sata.std()*nsig
    if len(sata)/10<ndiv: 
        ndiv=max([int(len(sata)/10),10])
        print 'reducing no. of bins for histog. to %i'%ndiv
    vstp=2*vstd/ndiv
    vmin=vmax-vstd
    vmax=vmax+vstd
    c=arange(vmin,vmax,vstp)
    if loud: print 'histog. from %4.2f to %4.2f'%(vmin,vmax)
    a,b=histogram(data,c)
    if loud>0:
        if loud>2: print b,a
        else: print 'histog. sum %i [%i above]'%(sum(a),sum(a>min_count))
    if sum(a>min_count)<10: bord=2
    else: bord=3
    if sum(a>min_count)<min_bins:
        print 'too few bins to fit a gausian'
        if rep==2: return b[argmax(a[bord:-bord])+bord],(vmax-vmin)/2
        else: return vstd
    #if rep>3: return get_width(b[a>min_count],a[a>min_count],bord=bord,rep=rep),b[a>min_count],a[a>min_count]
    pars=get_width(b[a>min_count],a[a>min_count],bord=bord,rep=rep)
    
    return pars

def clipped(val,nsig=3.,nclip=1,loud=0,rep=0):
    '''clipped average'''
    from numpy import bool
    sel=ones(val.shape,bool)
    for i in range(nclip):
        avg=val[sel].mean()
        std=val[sel].std()
        if loud>0: print "selecting %f +- %f"%(avg,std*nsig)
        sel=abs(val-avg)<std*nsig
        if sum(sel)==len(sel): break # no sense of repeating
    if rep==1: return val[sel].min(),val[sel].max()
    return sel

def ints(data,sigma,loud=0,min_wid=1):
    if sigma<0: sele=data>-sigma
    else: sele=abs(data)<2*sigma
    dlen=len(data)
    idi=arange(dlen)[sele]
    if len(idi)==0: 
        if loud>0:print 'no good bins'
        return []
    neidiff=idi[1:]-idi[:-1]
    #gstep=pidi[1:]-pidi[:-1]
    gidi=arange(len(neidiff))[neidiff>1] #indexy mezer
    if len(gidi)==0: return [(idi[0],idi[-1])] # no gaps, just one peak
    if gidi[-1]==len(sele)-1:gbegs=list(idi[gidi[:-1]+1])
    else:gbegs=list(idi[gidi+1])
    gbegs.insert(0,idi[0])
    gends=list(idi[gidi]+1)
    if gidi[-1]!=len(idi)-1: gends.append(idi[-1]+1)
    gtis=zip(gbegs,gends)
    if loud>0: print '%i good bins, %i intervals'%(sum(sele),len(gtis))
    return [a for a in gtis if a[1]-a[0]>min_wid] #good intervals longer than ...

def mask_data(data,bias,sigma=None):
    '''selection of good data
    if sigma>0: selects data within range bias+-sigma
    
    '''
    if sigma>0: sele=abs(data-bias)<sigma
    elif sigma<0: sele=abs(data-bias)>-sigma
    else:
        if bias!=None and sigma==None: sele=data>bias
        else:
            sele=ones(data.shape,dtype=bool)
    return sele

def trapezoid(e,e0,e1,w0,w1=-1,sl=None):
    from numpy import r_
    if w1<0:w1=w0
    rep=((e>e0+w0/2.)*(e<e1-w1/2.)).astype('float32')
    sel=abs(e-e0)*2<w0
    if sl:
        dif=(e1-e0+(w0-w1)/2)*sl/2
        ntop=int(sum(rep))
        a,b=1-dif,1+dif
        rep[rep>0]+=r_[-dif:dif:ntop*1j]
    else: a,b=1.,1.
    rep[sel]+=r_[0:a:sum(sel)*1j]
    sel=abs(e-e1)*2<w1
    rep[sel]+=r_[b:0:sum(sel)*1j]
    rep/=sum(rep)
    return rep

global labs,cnts
global avgs,xpos
def ndstats(data,sigma=0,bias=0,xbin=None,errs=None,bin=100,frac_min=0.75,step=None,loud=1):
    '''regroups values by "bin" bins using ndimage library
    also the x-bins and errors if provided
    other parameters as in extra.stats
    '''
    from scipy import ndimage
    global avgs,xpos
    if sigma!=None: sele=mask_data(data,bias,sigma)
    else: sele=ones(len(data),dtype=bool)
    global labs,cnts
    if step!=None and xbin!=None:
        if step<0:
            from numpy import median
            step=bin*median(xbin[1:]-xbin[:-1])
        labs=((xbin-xbin[0])/step).astype(int32)
        if sum(labs<0)>0:
            print 'x-axis not rising'
            return [],[],[]
        bin=int(len(xbin)/labs[-1])
        if loud: print 'using step %.3f, grouping in average by %i bins'%(step,bin) 
    else: labs=(arange(len(data))/bin).astype(int32)
    if sigma!=None: labs[sele==False]=-1
    cnts=zeros((max(labs)+1,))
    if loud>1: print "labels [%i-%i], length %i"%(min(labs),max(labs),len(cnts))
    for l in labs[sele]:
        cnts[l]+=1
    idx=arange(len(cnts))[cnts>=bin*frac_min]
    if len(idx)<1:
        print "all bins empty"
        return None
    else:
        if loud>0: print "%i bin(s) empty"%(len(cnts)-len(idx))
    #print 'max. index %i'%(labs[-1])
    avgs=ndimage.mean(data,labs,idx)
    if errs!=None: 
        errs=sqrt(array(ndimage.mean(errs**2,labs,idx))/bin)
    else: errs=sqrt(array(ndimage.variance(data,labs,idx)))
    if xbin==None: xbin=arange(len(data))
    xpos=ndimage.mean(xbin,labs,idx)
    #print 'check: data %i -> avgs %i labs %i idx %i'%(len(data),len(avgs),len(labs),len(xpos))
    return array(avgs),errs,array(xpos)
    
def stats(data,sigma=0,bias=None,xbin=None,errs=None,bin=100,frac_min=0.75,step=None,nclip=1):
    '''divides data in groups by "bin" bins and calculates average and variance in each group
    selects data using mask_data with parameters bias,sigma (if sigma!=None)
    frac_min - min. fraction of valid bins in a new bin
    errs - rebins also an error dataset
    '''
    global avgs,xpos
    if bias!=None: sele=mask_data(data,bias,sigma)
    if step!=None and xbin!=None:
        if step<0:
            from numpy import median
            step=bin*median(xbin[1:]-xbin[:-1])
        ndiv=int(xbin[-1]-xbin[0]//step)
    else:ndiv=len(data)//bin
    xpos=[]
    avgs=[]
    evgs=[]
    stds=[]
    cnts=bin
    ibeg=0
    iend=0
    if xbin!=None: xmark=xbin[0]
    else: 
        xbin=range(len(data))
        xmark=0
    for i in range(ndiv):
        if step!=None:
            ibeg=iend
            xmark+=step
            iend=int(sum(xbin<xmark))
            if iend>=len(data): break
        else:
            ibeg=i*bin
            iend=(i+1)*bin
        if bias!=None: cnts=sum(sele[ibeg:iend])
        if cnts>bin*frac_min:
            if bias!=None: sdata=data[ibeg:iend][sele[ibeg:iend]]
            else: sdata=data[ibeg:iend]
            if xbin!=None: 
                if bias!=None: xpos.append(sum(xbin[ibeg:iend][sele[ibeg:iend]])/cnts)
                else: xpos.append((xbin[ibeg]+xbin[iend-1])/2.)
            bavg=sum(sdata)/cnts
            bvar=(sum(sdata**2)/cnts-bavg**2)/cnts
            if bias==None and sigma>0: # sigma clipping
                while nclip>0:
                    nclip-=1
                    sdata=sdata[(sdata-bavg)*(sdata-bavg)<sigma*sigma*bvar]
                    cnts=len(sdata)
                    bavg=sum(sdata)/cnts
                    bvar=(sum(sdata**2)/cnts-bavg**2)/cnts
            if errs!=None:
                if bias!=None: eavg=sum(errs[ibeg:iend][sele[ibeg:iend]]**2)/cnts**2
                else: eavg=sum(errs[ibeg:iend]**2)/cnts**2
            else:
                eavg=bvar
            avgs.append(bavg)
            evgs.append(eavg)
            stds.append(bvar)
            #if stds[-1]<0: break
    if errs!=None:    return array(avgs),sqrt(array(evgs)),array(xpos),sqrt(array(stds))
    else: return array(avgs),sqrt(array(stds)),array(xpos)
    
    #cnts=array([ )
    #sums=array([sum(data[i*bin:(i+1)*bin][sele[i*bin:(i+1)*bin]]) for i in range(ndiv)])
    #sums2=array([sum(data[i*bin:(i+1)*bin][sele[i*bin:(i+1)*bin]]**2) for i in range(ndiv)])
    #xsums=array([sum(xbin[i*bin:(i+1)*bin][sele[i*bin:(i+1)*bin]]) for i in range(ndiv)])
    #bin*=frac_min #now bin is required minimal number of cummulated bins
    #avgs=sums[cnts>bin]/cnts[cnts>bin]
    #xpos=xsums[cnts>bin]/cnts[cnts>bin]
    #errs=sqrt((sums2[cnts>bin]-sums[cnts>bin]**2/cnts[cnts>bin]))
    #errs/=cnts[cnts>bin]
    #errorbar(xpos,avgs,errs)
    
def polyfit(x,y,r,w=None,full=0):
    ''' replacement for matplotlib/numpy polyfit,
    using also weights if needed
    '''
    from numpy.linalg import lstsq
    #from numpy.linalg import inv as inverse
    global chi2
    assert type(r)==int
    if len(y)<r+1:
        print 'too few data to fit a polynom'
        if full==0: return None    
        else: return None,None,None,None
    arrmat=ones((size(y),r+1),float64)
    if (w!=None) and (size(y)==size(w)):
        arrmat[:,0]*=w
        ymat=y*w
    else:
        ymat=y
    for i in range(r):
        arrmat[:,i+1]=arrmat[:,i]*x
    pars,chi2,rank,eigens=lstsq(arrmat,ymat)
    if full==0:return pars[::-1]
    return pars[::-1],chi2,rank,eigens

def rob_polyfit(x,y,wei=2,grad=True):
    #starting with rank2
    d=y.std()
    a=d/x.std()
    if wei==-1:
        c=((x-x.mean())*(y-y.mean())).mean()*a/d**2
        return a,c
    b=y.mean()-a*x.mean()
    if wei==0: return a,b
    from scipy.optimize import fmin
    fun=lambda p:sum((x*p[0]+p[1]-y)**2)+d*wei*(p[0]-a)**2+wei*(p[1]-b)**2
    if grad:
        derfun=lambda p:array([2*(sum((x*p[0]+p[1]-y)*x)+d*wei*(p[0]-a)),2*(sum(x*p[0]+p[1]-y)+wei*(p[1]-b))])
        from scipy.optimize import fmin_bfgs as fmin
        return fmin(fun,[a,b],derfun,disp=0)
    return fmin(fun,[a,b],disp=0)

def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                        scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    import numpy as n
    import scipy.interpolate
    import scipy.ndimage
    if not a.dtype in [n.float64, n.float32]:
        a = n.cast[float](a)

    m1 = n.cast[int](minusone)
    ofs = n.cast[int](centre) * 0.5
    old = n.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
            "This routine currently only support " \
            "rebinning to the same number of dimensions."
        return None
    newdims = n.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = n.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = n.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = n.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [n.arange(i, dtype = n.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + range( ndims - 1 )
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = n.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = n.mgrid[nslices]

        newcoords_dims = range(n.rank(newcoords))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (n.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print "Congrid error: Unrecognized interpolation type.\n", \
            "Currently only \'neighbour\', \'nearest\',\'linear\',", \
            "and \'spline\' are supported."
        return None
    

global glob_xern
glob_xern=arange(-2.,2.,.1)

from math import log
global xern,yern,kern,step
log_conv=log(10)

def fill_kde(data,errs,fact=2.,form='gauss',ndiv=100,xbase=None,xlog=False,clip=False):
    '''computes kernel density estimators
    see http://www.maths.uwa.edu.au/~duongt/seminars/appbkde/slides.pdf
    '''
    from numpy import exp,iterable
    global xern,kern,step
    #xids=[argmin(data),argmax(data)]
    if clip:
        sele=clipped(data)
        data=data[sele]
        if iterable(errs)>0:errs=errs[sele]
    if xbase==None:
        if iterable(errs)>0: xbase=[(data-errs*fact).min(),(data+errs*fact).max()]
        else: xbase=[min(data)-errs*fact,max(data)+errs*fact]
        #return xbase
    else:
        sel=(data>xbase[0])*(data<xbase[1])
        data=data[sel]
        if type(errs)==list: errs=errs[sel]
        #xbase=[data[xids[0]]-errs*fact,data[xids[1]]+errs*fact]
    step=(xbase[1]-xbase[0])/ndiv
    print 'range %5.2f-%5.2f step:%5.2f'%(xbase[0],xbase[1],step)
    if iterable(errs)>0:
        #errm=(max(errs)+min(errs))/2.
        errm=errs[errs>0].mean()
        errs[errs<=0]=errm
        xern=arange(-errm*fact,errm*fact,step)
        kern=None
    else:
        if xlog:
            #xshi=xern[0]-step
            #lshi=0
            xern=arange(-2.,2.,.1)
            kern=exp(-(exp(xern*log_conv)-1.)**2/2/(errs*log_conv)*2)
            kern-=max(kern[0],kern[-1])
            kern[kern<0]=0.
            #print 'correcting shift by %f and %f'%(xshi,cent)
        else: 
            xern=arange(-errs*fact,errs*fact,step)
            kern=exp(-xern**2/2/errs*2)
        kern/=sum(kern)
        if xlog:xern-=sum(kern*xern)
        errm=0
    axis=arange(xbase[0],xbase[1],step)
    xlen=len(xern)
    vals=zeros((ndiv),dtype=float64)
    for i in range(len(data)):
        pos=int((data[i]+xern[0]-xbase[0])/step+0.5)
        nos=0
        alen=xlen
        if pos<0: #below lowest
            nos=-pos
            pos=0
            #alen-=pos
        if pos+alen>ndiv: #above highest
            alen=ndiv-pos
        if nos+alen>=xlen:
            alen=xlen-nos-1
        if alen>=0:
            if errm>0:
                if xlog:
                    kern=exp(-(exp(xern*log_conv)-1.)**2/2/(errs[i]*log_conv)*2)
                    kern-=max(kern[0],kern[-1])
                    kern[kern<0]=0.
                else: kern=exp(-xern**2/2/errs[i]*2)
                kern/=sum(kern)
            try:
                vals[pos:pos+alen]+=kern[nos:nos+alen]
            except:
                print "casting %i:+%i to %i:+%i (%i / %i)"%(nos,alen,pos,alen,len(vals),len(kern))
    return vals,axis

def fill_k2de(xdata,xerrs,ydata,yerrs,fact=2.,form='gauss',xdiv=100,xbase=None,ydiv=100,ybase=None,xlog=False,ylog=False):
    '''computes kernel density estimators in 2 dimensions
    see http://www.maths.uwa.edu.au/~duongt/seminars/appbkde/slides.pdf
    '''
    from numpy import exp,iterable
    global xern,yern,kern,step
    #xids=[argmin(data),argmax(data)]
    if xbase==None:
        if iterable(xerrs)>0:xbase=[(xdata-xerrs*fact).min(),(xdata+xerrs*fact).max()]
        else: xbase=[min(xdata)-xerrs*fact,max(xdata)+xerrs*fact]
    if ybase==None:
        if iterable(yerrs)>0:ybase=[(ydata-yerrs*fact).min(),(ydata+yerrs*fact).max()]
        else: ybase=[min(ydata)-yerrs*fact,max(ydata)+yerrs*fact]
    
    xstep=(xbase[1]-xbase[0])/xdiv
    if ydiv==0: 
        ystep=xstep
        ydiv=(ybase[1]-ybase[0])//ystep
        ybase[1]=ybase[0]+ystep*ydiv
    else: ystep=(ybase[1]-ybase[0])/ydiv
    print 'step:%5.2f x %5.2f'%(xstep,ystep)
    if iterable(xerrs)>0:
        #xerrm=(max(xerrs)+min(xerrs))/2.
        xerrm=xerrs[xerrs>0].mean()
        xerrs[xerrs<=0]=xerrm
        xern=arange(-xerrm*fact,xerrm*fact,xstep)
        kern=None
    else:
        if xlog:
            xern=exp(glob_xern*log_conv)-1.
            kern=exp(-xern**2/2/xerrs**2)
            kern-=max(kern[0],kern[-1])
            kern[kern<0]=0.
            kern/=kern.sum()
            xern-=sum(kern*xern)
        else:
            xern=arange(-xerrs*fact,xerrs*fact,xstep)
        xerrm=0
    if iterable(yerrs)>0:
        #yerrm=(max(yerrs)+min(yerrs))/2.
        yerrm=yerrs[yerrs>0].mean()
        yerrs[yerrs<=0]=yerrm
        yern=arange(-yerrm*fact,yerrm*fact,ystep)
        kern=None
    else:
        if ylog:
            yern=exp(glob_xern*log_conv)-1.
            yerrs*=log_conv
            kern=exp(-yern**2/2/yerrs**2)
            kern-=max(kern[0],kern[-1])
            kern[kern<0]=0.
            kern/=kern.sum()
            yern-=sum(kern*yern)
        else:
            yern=arange(-yerrs*fact,yerrs*fact,ystep)
        yerrm=0
    
    xlen=len(xern)
    ylen=len(yern)
    xern=xern.reshape(xlen,1)
    yern=yern.reshape(1,ylen)
    if xerrm+yerrm==0: # both axes with constant errors
        kern=exp(-xern**2/2/xerrs**2-yern**2/2/yerrs**2)
        kern/=sum(kern)
    elif xerrm==0: xvals=exp(-xern**2/2/xerrs**2)
    elif yerrm==0: yvals=exp(-yern**2/2/yerrs**2)
    else:
        print 'mean errors:%5.2f x %5.2f'%(xerrm,yerrm)
    xaxis=arange(xbase[0],xbase[1],xstep)
    yaxis=arange(ybase[0],ybase[1],ystep)
    vals=zeros((xdiv,ydiv),dtype=float64)
    for i in range(len(xdata)):
        xpos=int((xdata[i]+xern[0,0]-xbase[0])/xstep+0.5)
        ypos=int((ydata[i]+yern[0,0]-ybase[0])/ystep+0.5)
        xnos=0
        axlen=xlen
        if xpos<0: #below lowest
            xnos=-xpos
            xpos=0
            axlen-=xnos
        if xpos+axlen>xdiv: #above highest
            axlen=xdiv-xpos
        if xnos+axlen>=xlen:
            axlen=xlen-xnos-1
        
        ynos=0
        aylen=ylen
        if ypos<0: #below lowest
            ynos=-ypos
            ypos=0
            aylen-=ynos
        if ypos+aylen>ydiv: #above highest
            aylen=ydiv-ypos
        if ynos+aylen>=ylen:
            aylen=ylen-ynos-1
        
        if axlen>=0 and aylen>=0:
            if xerrm>0:
                if yerrm>0: kern=exp(-xern**2/2/xerrs[i]**2)*exp(-yern**2/2/yerrs[i]**2)
                else: kern=exp(-xern**2/2/xerrs[i]**2)*yvals
            elif yerrm>0: kern=xvals*exp(-yern**2/2/yerrs[i]**2)
            vals[xpos:xpos+axlen,ypos:ypos+aylen]+=kern[xnos:xnos+axlen,ynos:ynos+aylen]
    return vals,xaxis,yaxis

def rot_mat(p):
    from math import sin,cos,sqrt
    v=sin(p)
    w=abs(v)>1e-2 and sqrt(1-v**2) or (abs(v)>1e-4 and 1-v**2/2. or 1)
    return array([[w,-v],[v,w]])

global xarr
xarr=None
def chi2map_anal(arep,grid,shallow=.3,bord_frac=.3,blofit='gauss',rep=0,inparg=None):
    '''finds minimum of chi2map
        blofit: gauss
                origauss: oriented 2-D gaussian
    '''
    global xarr
    from numpy import exp
    smin=arep.argmin()
    spos=[smin//len(grid[0]),smin%len(grid[0])]
    amin=arep[spos[0],spos[1]]
    print 'min. value %.2f at %i,%i : %.3f,%.3f'%(amin,spos[0],spos[1],grid[0][spos[1]],grid[1][spos[0]])
    sbord=((arep[0,0]+arep[-1,-1])/2.-amin)*shallow
    if sum(arep-amin>sbord)<arep.size*bord_frac:
        print 'too flat map'
        return arep
    if blofit=='origauss':
        from numpy import dot
        e=lambda v,x,y:(v[0]*exp(-((dot(x,rot_mat(v[6]))-v[1:3])**2/v[3:5]).sum(1))+v[5]-y)
        v0=array([1.,0.,0.,1.,1.,amin,0.])
        if inparg!=None: 
            v0[-1]=inparg[-1]
            inparg=inparg[:-1]
    elif blofit=='gauss': 
        e=lambda v,x,y:(v[0]*exp(-((x-v[1:3])**2/v[3:5]).sum(1))+v[5]-y)
        v0=array([1.,0.,0.,1.,1.,amin])
    if inparg!=None: v0[1:len(inparg)+1]=inparg
    else: 
        e=lambda v,x,y:(v[0]+((x-v[1:3])**2/v[3:5]).sum(1)-y) # or maybe just paraboloid
        v0=array([amin,0.,0.,1.,1.])
    if rep==1: return e,v0
    from scipy.optimize import leastsq
    gsize=len(grid[0])*len(grid[1])
    xarr=array([[[dx,dy] for dy in grid[1]] for dx in grid[0]]).reshape(gsize,2)
    if rep==2: return e,v0,xarr,arep.ravel()
    vsol,corr,expar,comm,res=leastsq(e,v0,(xarr,arep.ravel()),full_output=2)
    if rep==-1: return vsol,corr,expar,comm,res
    return vsol[1:],corr[1:,1:]

def lin_errors(x,y):
    from math import sqrt
    n=len(x)
    sxx=((x-x.mean())**2).sum()
    syy=((y-y.mean())**2).sum()
    sxy=((x-x.mean())*(y-y.mean())).sum()
    s=sqrt((syy-sxy**2/sxx)/(n-2))
    return s*sqrt(1/n+x.mean()**2/sxx),s/sqrt(sxx)

def com_abs(data):
    return abs(data.real)+1j*abs(data.imag)

def rej_std(data,dir=0,mode='quad',frac=0.6):
    '''rejecting most divergent value
    '''
    from numpy import mean,array,sum
    duff=array(data)-mean(data,dir)
    if mode[:3]=='abs':duff=abs(duff)
    else: duff*=duff
    sele=duff<duff.max(dir)
    cor_sele=sum(sele,dir)<len(duff)*frac
    sele[:,cor_sele]=True
    if mode=='test': return duff,sele
    muff=sum(array(data)*sele,dir)/sum(sele,dir) #corrected mean
    duff=array(data)-muff
    if mode[:3]=='abs':duff=abs(duff)
    else: duff*=duff
    return sqrt(sum(duff*sele,dir)/(sum(sele,dir)-1))
    # return sqrt((duff.sum(0)-duff.max(0))/(len(data)-1))

def scan_kde(pars,kmap,xaxis,yaxis,grid=None,nstep=20,blofit=False):
    if grid==None: grid=[1.,1.]
    if type(grid[0])!=array: grid=[arange(-a,a,2*float(a)/nstep) for a in grid]
    ystep=(yaxis[-1]-yaxis[0])/(len(yaxis)-1)
    rep=[]
    cnt=[]
    xpos=arange(len(xaxis))
    for dx in grid[0]:
        yval=(pars[0]+dx)*xaxis+pars[1]
        ybnd=[yval[0],yval[-1]]
        if ybnd[0]>ybnd[1]:ybnd=ybnd[::-1]
        if ybnd[0]+grid[1][0]<yaxis[0] or ybnd[1]+grid[1][-1]>yaxis[-1]: #need to limit the range
            sele=(yval+grid[1][0]>yaxis[0])*(yval+grid[1][-1]<=yaxis[-1])
            if sum(sele)<2:
                print 'for slope %.3f: out of range'%(pars[0]+dx)
                cnt.append(0)
                rep.append([0.]*len(grid[1]))
                continue
            yval=yval[sele]
            if loud>1:
                print 'for slope %.3f: going from %.3f to %.3f: %i'%(pars[0]+dx,yval[0],yval[-1],sum(sele))
            rep.append([kmap[xpos[sele],((yval-yaxis[0]+dy)/ystep).astype(int)].sum() for dy in grid[1]])
            cnt.append(len(yval))
        else:
            rep.append([kmap[xpos,((yval-yaxis[0]+dy)/ystep).astype(int)].sum() for dy in grid[1]])
            cnt.append(len(xpos))
    arep=array(rep).transpose()
    if blofit:
        return chi2map_anal(arep,[grid[0]+pars[0],grid[1]+pars[1]],blofit=blofit)
    return arep,cnt,grid

def confidence(perc,mean,sprd,n,met_cent='median'):
    '''confidence regions
        sprd is mean((x-mean)**2)
    '''
    from scipy import stats
    from math import sqrt
    if perc>1: perc/=100.
    dist=sqrt(sprd/(n-1))*stats.t.ppf((1+perc)/2.,n-1)
    rep=[mean-dist,mean+dist]
    # confidence intervals for deviation
    dist=sqrt(n*sprd/2)
    #cent=stats.gengamma.cdf(sigma/dist,(n-1)/2,-2)
    if met_cent=='peak': cent=stats.gengamma.cdf(sqrt(2/n),(n-1)/2,-2)
    else: cent=1/2.
    dist1=dist*stats.gengamma.ppf(cent-perc/2.,(n-1)/2,-2)
    dist2=dist*stats.gengamma.ppf(cent+perc/2.,(n-1)/2,-2)
    rep+=[dist1,dist2]
    return rep
            
def corr_dist(xcor,ycor,slop,bias,perp=False,xerr=None,yerr=None,square=True):
    '''returns perpendicular distances from line given by bias (intercept.) and slope
    '''
    yarr=(ycor-bias)
    #if yerr!=None: assert all(yerr>0)
    xarr=xcor
    if perp and slop!=0.:
        if yerr!=None:
            rer2=(xerr/yerr)**2
            cosa=(yarr*slop*rer2+xarr)
            dif=yarr**2*rer2+xarr**2
            dif-=cosa**2/(1+slop**2*rer2)
            if not square: dif=sqrt(dif)
        else:
            if square: dif=(yarr-slop*xarr)**2/(1+slop**2)
            else: dif=abs(yarr-slop*xarr)/sqrt(1+slop**2)
            #print 'dif: %f'%sum((dif-dif2)**2)
        if xerr!=None: dif/=xerr
    else:
        if square:dif=(yarr-slop*xarr)**2
        else:dif=abs(yarr-slop*xarr)
        if yerr!=None: dif/=yerr
    return dif

global arep
arep=None
def corr_anal(xarr,yarr,xerr=None,yerr=None,xlog=False,ylog=False,selv=None,grid=None,perp=True,blofit=None,simp=False,
    xbaserr=1.,ybaserr=1.):
    '''gets linear regression and maps the chi2 minimum
    if perp: measures perpendicular distance to regression line (taking in account event. errors in both axes)
    
    in principle, the expected value should be the barycenter of the region with chi2 within some distance of the minimum
    not simply the best-chi2 value
    see e.g. D'Agostini: http://arxiv.org/abs/physics/0403086
    '''
    global arep
    from numpy import dot,log
    n=len(xarr)
    dsum=0
    eslop=0
    ebias=0
    if selv==None: selv=ones((n,),bool)
    if xlog: selv*=xarr>0
    if ylog: selv*=yarr>0
    n=sum(selv)
    if n<=2:
        print 'not enough data'
        return
    xcor=xarr[selv]
    ycor=yarr[selv]
    if xlog:
        xcor=log(xcor)
        if xerr!=None:
            assert all(xerr>0)
            xerr=log(xerr)
    if xerr!=None:
        assert all(xerr>0)
        xavg=sum(xcor/xerr)/sum(1./xerr)
    else: xavg=xcor.mean()
    if ylog:
        ycor=log(ycor)
        if yerr!=None:
            assert all(yerr>0)
            yerr=log(yerr)
    if yerr==None: yerr=xerr
    if yerr!=None:
        assert all(yerr>0)
        yavg=sum(ycor/yerr)/sum(1./yerr)
    else: yavg=ycor.mean()
        
    xcor-=xavg
    ycor-=yavg
    if xerr!=None: xcor2=xcor/xerr
    else: xcor2=xcor
    if yerr!=None: ycor2=ycor/yerr
    else: ycor2=ycor
    xycor=dot(xcor2,ycor2)
    yycor=dot(ycor2,ycor2)
    xxcor=dot(xcor2,xcor2)
    slop=xycor/xxcor
    bias=yavg-slop*xavg
    dsum=sqrt((yycor-xycor**2/xxcor)/(n-2))
    dbias=0
    if not simp:
        r=1
        arrmat=ones((n,r+1),float64)
        wei=zeros(n,float64)
        if xerr!=None:
            if xlog: wei+=xerr.max()-xerr+xbaserr # already in logarithmic
            else: wei+=log(xerr.max())-log(xerr)+xbaserr
        if yerr!=None:
            if xlog: wei+=log(yerr.max())-log(yerr)+ybaserr # already in logarithmic
            else: wei+=log(yerr.max())-log(yerr)+ybaserr
        if wei.sum()>0:
            wei/=wei.sum()
            arrmat[:,0]*=wei
            ycor2=ycor*wei
            if loud>2: print "weight:%s"%wei
        else: ycor2=ycor
        for i in range(r):
            arrmat[:,i+1]=arrmat[:,i]*xcor
        from numpy.linalg import lstsq
        pars,chi2,rank,eigens=lstsq(arrmat,ycor2)
        #return pars,chi2
        fbias,slop=tuple(pars)
        bias=yavg-slop*xavg
        #print 
    if dsum>0:
        eslop=dsum/sqrt(xxcor)
        ebias=dsum*sqrt(1./n+xavg**2/xxcor)
    print 'correlation (%5.2f) : %5.2f(+-%5.2f) *x + %5.2f(+-%5.2f)'%(dsum,slop,eslop,bias,ebias)
    
    if grid==None:
        if loud==2: return corr_dist(xcor,ycor,slop,bias,perp,xerr,yerr)
        else: return slop,bias,corr_dist(xcor,ycor,slop,bias,perp,xerr,yerr).std()
    if type(grid[0])==int:
        if eslop==0: eslop=0.1*slop
        if bias==0: ebias=0.1*bias
        grid=[arange(-grid[i],grid[i])*[eslop*ybaserr,ebias*xbaserr][i] for i in range(2)]
    rep=[]
    for dy in grid[1]:
        rep.append([(corr_dist(xcor,ycor,slop+dx,dbias+dy,perp,xerr,yerr)**2).sum() for dx in grid[0]])
    #arep=array(rep)
    if blofit:
        return chi2map_anal(array(rep),[grid[0]+slop,grid[1]+bias],blofit=blofit)
        #v=array([1.,0.,1.,0.,1.])
        #x=indices((6,6))-array([2.5,2.5]).reshape(2,1,1)
        #array([a.reshape(36,) for a in idis]).swapaxes(0,1)
    return [grid[0]+slop,grid[1]+bias]
