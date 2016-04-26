# -*- coding: utf-8 -*-
'''
basic tools for spectral analysis
probably more sophisticated algorithms could be extracted from scipy.signal
'''

ev2um=0.806554243769
ev2K=11604.5026029 # inv. boltzman. constant in eV
c0=2.99792e8 # lightspeed
h=4.1356673e-15 # eV s
wien=2898. # um K - Wien's law

global base
step=0.96423 # typical step in datafile
base=None
back=None
glob_fit_mode=0

###################### some low-energy asymptotes
cau_sio2=[1.68683e-05,0,1.961953e-03,0,1.4493423] # cauchy dispersion fit, for energy in eV

## low-e part of silicon - real
cau_si_low=[2.71e-02,1.95e-03]
exp_si_low=[11.79540,0.08392, -0.039056] # exp_si_low[0]*exp(e*(exp_si_low[0]*e+exp_si_low[1]))
cau_si_low2=[11.7124,0.0791917,0.5635043] # in e**2

#lor_si_low=[  7.0144,  [[56.29,   3.3303,   0.624]]]
lor_si_low=[  2.7977,  [[1.27805056e+02,   3.78364955e+00,5.49466137e-02]]] # precise from Aspnes
#lor_si_low2=[ 4.2298 ,  [[68.524,  3.8118,  0.8], [9.0712,  3.0635,   0.889],  [20.6587, 3.3763 ,   0.274]]]
lor_si_low2=[-2.61678, [[  1.79418744e+02,   5.41916478e+00,   1.00000000e-03],[  16.3505958,   3.17564,   1.21380228e-01],   [  1.38945163e+02,   4.60986142e+00,   1.00000000e-03]]]# precise from Aspnes

lor_si_hig=[ 1.9330, [[13.5987,  5.19495,  0.4829],[101.8363, 4.1194,  0.8273], [37.5503,  3.45024,    0.4094]]] # up to 5eV

herz_si=[2.65456e-2,1.23172e-1,3.14906,-2.66511e-8,5.45852e-14,0] 
# herberger - n=polyval(herz_si[:3],1/(l**2-0.028))+polyval(herz_si[3:],l**2)  pro l z (1.1-590 um = transparent)
sell_si=[9.39816e-1,11.6858,8.10461e-3]
# sellmeyer n=polyval(sell_si[:2],1/l**2)+sell_si[2]*l1**2/(l**2-l1**2) pro l1=1.1071

tau_si=[[1555.9980, 3.33024, 0.393927, 3.081100], [18.1562, 4.21923, 0.53848, 1e-4], [99.71856, 5.20215, 0.494745, 4.63455],0.58772]

def respect(x,e,dir=None):
    # spectral rebinning
    from scipy import interpolate
    tck=interpolate.splrep(e[0][::dir],e[1][::dir])
    rep=interpolate.splev(x,tck)
    if len(e)>2:
      tck=interpolate.splrep(e[0][::dir],e[2][::dir])
      rep=rep.astype('c16')+1j*interpolate.splev(x,tck)
    return rep

def rebin(x,y,nbins=100,type='inter',dir=1):
    '''converting nonuniform (eg. reciprocal) bins (wavelength) to energy
    either simple averaging or spline interpolation
    '''
    from extra import ndstats
    step=abs(x[0]-x[-1])/nbins
    if len(x)<nbins*2: type='inter'
    if type=='stats': return ndstats(y[::dir],xbin=x[::dir],sigma=None,bin=2,step=step,frac_min=0.2,loud=0)
    from scipy import interpolate
    tck=interpolate.splrep(x[::dir],y[::dir])
    nx=arange(x.min(),x.max(),step)
    return interpolate.splev(nx,tck),None,nx

from math import pi

def get_ids(xb,sub_band):
    # bin indices corresponding to given energy band in binning 'xb'
    if (xb[0]-xb[-1])*(sub_band[0]-sub_band[1])<0:
        sub_band=(sub_band[1],sub_band[0])
    if sub_band[0]>sub_band[1]:
        amin=sum(xb>sub_band[0])
        amax=sum(xb>sub_band[1])
    else:
        amin=sum(xb<sub_band[0])
        amax=sum(xb<sub_band[1])
    if amin>amax: amin,amax=amax,amin
    return amin,amax

def corr(lst,rep,shift=0,scal=1.,sel=None):
    '''correlation
    '''
    out=dict()
    if type(scal) not in (list,dict): a=scal
    if type(shift) not in (list,dict): b=shift
    for i in range(len(lst)):
        s=lst[i]
        if type(scal)==list: a=scal[i]
        elif type(scal)==dict: a=scal[s]
        if type(shift)==list: b=shift[i]
        elif type(shift)==dict: b=shift[s]
        if not s in rep: continue
        if sel!=None: out[s]=rep[s][sel[0]:sel[1]]*a+b
        else: out[s]=rep[s]*a+b
    return out
    
def dbload(name='SiO2_gl',freq=None,connection="http"):
    from numpy import loadtxt,arange,array
    if freq!=None:
        if freq[0]>freq[-1]: freq=freq[::-1]
    if connection=="http":
        try:
            import urllib.request as urllib2
        except ImportError:
            import urllib2 
        if name==None:
            query='http://physics.muni.cz/~munz/nanoCHARM/cgi-bin/base.cgi?show=mater'
            aa=urllib2.urlopen(query)
            return [a.strip() for a in aa.readlines()]
        query='http://physics.muni.cz/~munz/nanoCHARM/cgi-bin/base.cgi?show=0&sample=%s&epsil=1'%name
        if freq!=None:
            imin=freq[0]
            imax=freq[-1]
            query+='&emin=%.3f&emax=%.3f'%(imin,imax)
        aa=urllib2.urlopen(query)
        if aa==None:
            print('cannot load data for '+name)
            return
        skip=0
    else:
        query='select energy,reepsilon,imepsilon from nanoCHARM where vzorek="%s" order by energy;'%name
        import os
        try:
            aa=os.popen('''echo '%s' | ssh %s "mysql ufkl -u ufklread -p'CTEpris#'"'''%(query,connection))
        except:
            print('failed '+'''echo '%s' | ssh %s "mysql ufkl "'''%(query,connection))
            return
        skip=1
    data=array([a.strip().split() for a in aa.readlines()[skip:]]).astype(float).transpose()
    #data=loadtxt(aa,unpack=True,skiprows=skip)
    if freq!=None:
        from scipy import interpolate
        if False: #not sorted
            from numpy import argsort
            sel=(data[0]>freq[0])*(data[0]<freq[-1])
            cord=argsort(data[0][sel])
            cord=arange(data[0].size)[sel][cord]
        else:
            cord=(data[0]>freq[0])*(data[0]<freq[-1])
        if sum(cord)<3:
            print('too few points tabulated')
            return
        if sum(cord)==3:
            pixs=arange(len(cord))[cord]
            cord[pixs[0]-1]=True
            cord[pixs[-1]+1]=True
        if len(freq)==2: #only min/max values given, not rebinning
            return data[0][cord],data[1][cord]+1j*data[2][cord]
        val_re=interpolate.splev(freq,interpolate.splrep(data[0][cord],data[1][cord]))#,xb=freq[0],xe=freq[-1]))
        val_im=interpolate.splev(freq,interpolate.splrep(data[0][cord],data[2][cord]))#,xb=freq[0],xe=freq[-1]))
        return freq,val_re+1j*val_im
    return data[0],data[1]+1j*data[2]
    
def dbsave(base,data,mat='t1034_com',skup='metals_meta',matskup=None,connection="munz@mono",nper=10,passwd=None,shi=1):
    vars='material,skupina,podskupina,komentar'
    if matskup==None: matskup=mat[:mat.find('_')]
    query='insert into bessy_spectra ('+vars+') values ("'+mat+'","'+skup+'","'+matskup+'","composite spectrum");'
    vars='poradi,energy,reepsilon,imepsilon,vzorek'
    query='insert into bessy ('+vars+') values (%i,%f,%f,%f,"'+mat+'");'
    fullq=''
    import os
    if passwd==None:
        print('need a password')
        return
    for i in range(len(data)):
        fullq+=query%(i+shi,base[i],data[i].real,data[i].imag)
        if (i+1%nper)==0:
            aa=os.popen('''echo '%s' | ssh %s mysql -u hemzal -p%s ufkl'''%(fullq,connection,passwd))
            fullq=''
    if fullq!='':aa=os.popen('''echo '%s' | ssh %s mysql -u hemzal -p%s ufkl'''%(fullq,connection,passwd))


global sele,foure
sele=None
foure=None
def filter(dat,xb=None,min_sig=0,sub_band=None,pol_deg=2,min_ampl=1.,init_wid=0,enlarg=4,do_ints=True,min_wid=0,
        rep=None,clean=True,weight=False,resid=False):
    '''looks for peaks in fourier image
    erases all but the first one
    subtracts polynomial fit first (of degree pol_deg)
    init_wid: width of initial (low freq) peak to be preserved (auto-estimated if zero)
    in FFT: freq=(2*pi)*max(fft)/(t_max-t_min)*n_bin_time/n_bin_freq
    (the last ratio is usually 1. - FFT returns array of the same number of samples)
    '''
    global sele,foure
    from numpy import array,arange,polyfit,polyval,real,convolve,ones,argmax,abs
    from numpy.fft import fft,ifft
    if xb==None: xb=base
    if xb!=None:
        if sub_band!=None:
            amin,amax=get_ids(xb,sub_band)
            ib=xb[amin:amax]
            ae=dat[amin:amax]
            print('selecting bins %i:%i'%(amin,amax))
        else:
            ib=xb
            ae=dat
    else:
        ib=arange(len(dat))
        ae=dat
    ids=polyfit(ib,ae,pol_deg)
    scal=ae.mean()
    be=ae-polyval(ids,ib)
    foure=fft(be/scal)
    re=abs(foure)
    #rep=[]
    if min_ampl<0: return re
    if min_ampl==0:
        min_ampl=re.mean()*5.
        print('setting peak lower limit to %.2f'%min_ampl)
    if do_ints:
        import extra
        reg=extra.ints(re,-min_ampl,min_wid=min_wid) # intervals of Four. spec. above min_ampl
        print('found %i intervals [%i bins] - sample length %.1f'%(len(reg),sum([a[1]-a[0] for a in reg]),(ib[-1]-ib[0])))
        if rep=='ints': return reg
        if clean:
            if len(reg)<2:
                print('no frequency cleaned (max. %.2f)'%re[10:-10].max())
                if rep==None: return ae
            if reg[0][0]<=1:
                del reg[0]
                del reg[-1]
        for r in reg:
            if clean: 
                if enlarg>=0: foure[r[0]-enlarg:r[1]+enlarg]=0j
                else: # tracing
                    if enlarg<=-2: # from top
                        cent=r[0]+argmax(re[r[0]:r[1]])
                        r=(cent,cent)
                    for i in range(r[0],0,-1):
                        if re[i]<re[i-1]: break
                    for j in range(r[1],len(re)):
                        if re[j]>re[j-1]: break
                    r=(i,j)
                    foure[r[0]:r[1]]=0j
            wei=re[r[0]:r[1]].sum()
            cent=(arange(r[0],r[1])*re[r[0]:r[1]]).sum()/wei #barycenter
            if rep!=None:
                if weight: rep.append([2*pi*cent/(ib[-1]-ib[0]),wei])
                else: rep.append(2*pi*cent/(ib[-1]-ib[0]))
            if 2*r[1]<len(ae):
                print('%i bins cent %.2f [%.2f] max %.2f'%(r[1]-r[0],cent,2*pi*cent/(ib[-1]-ib[0]),re[r[0]:r[1]].max()))
                #rep.append(cent*(ib[-1]-ib[0])/len(re))
        if enlarg<-2:
            reg2=extra.ints(re,-min_ampl,min_wid=min_wid)
            print('found %i extra peaks'%len(reg2))
            reg+=reg2
            if rep!=None:
                for r in reg2:
                    wei=re[r[0]:r[1]].sum()
                    cent=(arange(r[0],r[1])*re[r[0]:r[1]]).sum()/wei #barycenter
                    if weight: rep.append([2*pi*cent/(ib[-1]-ib[0]),wei])
                    else: rep.append(2*pi*cent/(ib[-1]-ib[0]))
        #if resid: 
    else:
        sele=real(foure)>min_ampl
        if init_wid==0: init_wid=list(sele).index(False) # where ends initial peak
        sele[:init_wid]=False
        sele[-init_wid:]=False
        isum=sum(sele)
        if isum==0:
            print('no frequency cleaned')
            return ae
        if enlarg>0:
            sele=convolve(sele,ones(2*enlarg+1),'same')>0
            print('cleaning %i (prev %i) tim. bins'%(sum(sele),isum))
        foure[sele]=0j
    if rep!=None:
        return rep
    ge=ifft(foure)*scal
    return real(ge)+polyval(ids,ib)

def basel(lst,rep,bas,lim,ord=1,fit=None):
    '''baseline subtraction below peaks
    currently lorentzian/gaussian fit to peaks available
    lst: list of names
    rep: repository
    bas: x-axis baseline
    lim: localization of the peak
    '''
    amin,amax=get_ids(bas,lim)
    from numpy import polyfit,polyval,abs,exp
    bin=bas[amin:amax]
    if fit!=None:
        from scipy.optimize import leastsq
        #fun=lambda par:par[0]/((x-par[1])**2+par[2]**2)+par[3]-y
        if fit[:3]=="voi":
            from profit import profiles
            fun=lambda par:profiles(bin-par[1],[par[0],par[2],par[4],par[5]])+par[3]-val
        elif fit[:3]=="lor":fun=lambda par:par[0]/(((bin-par[1])*par[2])**2+1)+par[3]-val
        else:fun=lambda par:par[0]*exp(-((bin-par[1])/par[2])**2)+par[3]-val
        par=[1.,bas[(amin+amax)//2],200.,0.]
    print('using %i bins (%i:%i = %.2g:%.2g)'%(len(bin),amin,amax,bin[0],bin[-1]))
    qep=[]
    for s in lst:
        if s not in rep: continue
        idx=polyfit(bin,rep[s][amin:amax],ord)
        val=rep[s][amin:amax]-polyval(idx,bin)
        flus=val.std()*2.
        #ids=arange(len(val))[abs(val)<flus]
        idx2=polyfit(bin[abs(val)<flus],val[abs(val)<flus],ord)
        print('%s: using %i bins for bkg'%(s,sum(abs(val)<flus)))
        val-=polyval(idx2,bin)
        if fit!=None:
            pos=val.argmax()
            if -val.min()>val[pos]: 
                pos=val.argmin()
                print(s+':absorption')
            par[0]=val[pos]
            par[1]=bin[pos]
            out=leastsq(fun,par,full_output=2)
            qep.append((val,idx+idx2,out[0]))
        else: qep.append((val,idx+idx2))
    return qep

global par
par=None

def phasogram(x,y,frq,nphas=20,count=False):
    '''folding with frequency frq
    '''
    from scipy import ndimage
    from numpy import arange,int32,array,where
    ph=x*frq
    ph-=ph.astype(int32)
    iph=(ph*nphas).astype(int32)
    rep=array(ndimage.mean(y,iph,arange(nphas)))
    if count: 
        iph.sort()
        sumcnt=where(iph[1:]-iph[:-1]>0)
        holes=where(iph[1:]-iph[:-1]>1)
        cnt=[sumcnt[0][0]]+list(sumcnt[0][1:]-sumcnt[0][:-1])
        return rep,cnt,len(holes[0])
    return rep

def perfind(x,y,cent,wid,nfrq=20,nphas=20,max_null_frac=0):
    '''find the most probable frequency using phasogram
    (searching for the value giving the highest variation of folded data)
    '''
    rep=[]
    from numpy import array,arange
    frqis=arange(cent-wid,cent+wid,wid/nfrq)
    for frq in frqis: # find that "least uniform"
        phg=phasogram(x,y,frq,nphas)
        if max_null_frac>0:
            frac_null=sum(phg<=0.)/float(nphas)
            if frac_null>max_null_frac:
                print('frq %.2f makes %.0f%% empty bins: possible aliasing'%(frq,frac_null*100))
                rep.append([frq,phg[phg>0].mean(),-1])
            else:
                rep.append([frq,phg.mean(),phg.std()])
        else: rep.append([frq,phg.mean(),phg.std()])
    return array(rep).transpose()
    # period=33.3376

def perslide(x,y,cent,wid=10,nstep=20,nphas=20):
    '''sliding period estimation
    '''
    rep=[]
    from numpy import arange,int32,array
    if type(wid)==float: wid=sum(x>wid+x[0])
    step=(len(x)-wid)//nstep
    from scipy import ndimage
    ph=x*cent
    ph-=ph.astype(int32)
    iph=(ph*nphas).astype(int)
    for i in range(0,nstep*step,step):
        phg=array(ndimage.mean(y[i:i+wid],iph[i:i+wid],arange(nphas)))
        # find that "least uniform"
        sel=phg>0
        rep.append([i,phg[sel].mean(),phg[sel].std()])
    return array(rep).transpose()

def beats(b,bin=20,rep=0,mev=10,lim=0.7):
    '''finds positions of beat packets
    rep=3: returns positions of minima
    rep=2:
    '''
    imax=b.size//bin
    bmin=b[:imax*bin].reshape(imax,bin).min(1)
    bmax=b[:imax*bin].reshape(imax,bin).max(1)
    if rep==0:return bmax-bmin
    else: bdif=bmax-bmin
    from numpy import array,argmin,arange
    koo=array([i+argmin(bdif[i:i+mev]) for i in range(len(bdif))]) #koo is sorted by definition
    sel=arange(len(koo))[koo[1:]-koo[:-1]>0]
    dsel=sel[1:]-sel[:-1]
    mee=arange(len(sel))[dsel>int(0.7*mev)]
    pos=koo[sel[mee]+1]
    if lim>0: pos=pos[bdif[pos]<bdif.mean()*lim] #some false minima
    if rep==3: return zip(pos*bin,bdif[pos])
    if rep==2: return zip(pos*bin,dsel[dsel>int(0.7*mev)])
    return pos*bin #positions of beat minima

def ext_ran_fit(pos,x,b,extsel=5,extchi=False,extamp=False,rang=2,alim=None,loud=1):
    '''should return extrema position ev. amplitude
    extsel: half-size of analysed region     
    '''
    from numpy import polyfit,polyval
    rep=[] #positions
    arep=[] #amplitudes
    crep=[] #fit qual.
    for i in pos:
        if extsel==0:
            if i==0: rep.append([0,0])
            else: rep.append([i+extsel,i-extsel])
            continue
        if i==0: 
            rep.append(0)
            arep.append(0)
            crep.append(0)
        if i<extsel: continue
        if i+extsel>len(b): break
        idx=polyfit(b[i-extsel:i+extsel],x[i-extsel:i+extsel],rang)
        base=alim>0 and x[i-extsel:i+extsel].min() or x[i-extsel:i+extsel].max()
        xbary=(abs(x[i-extsel:i+extsel]-base)*b[i-extsel:i+extsel]).sum()/abs(x[i-extsel:i+extsel]-base).sum()
        if alim: # limit for the highest order coef.
            if (alim<0 and idx[0]>alim) or (alim>0 and idx[0]<alim):
                rep.append(xbary)#b[i])
                if loud>0: print('spectra: problems in polynom fitting around %.4f (curv. %.4f)'%(b[i],idx[0]))
                if extamp: arep.append(alim>0 and x[i-extsel:i+extsel].max() or x[i-extsel:i+extsel].min())
                if extchi: crep.append(-1)
                continue
        xpos=-idx[1]/2/idx[0]
        if xpos<b[i-extsel] or xpos>b[i+extsel-1]: 
            rep.append(xbary)
            if extamp: arep.append(alim>0 and x[i-extsel:i+extsel].max() or x[i-extsel:i+extsel].min())
            if loud>0: print('spectra: problems in polynom fitting around %.4f (%.4f)'%(b[i],xpos-b[i]))
        else: 
            rep.append(xpos)
            if extamp: arep.append(idx[2]+xpos*idx[1]/2.)
        if extchi: crep.append(sum((x[i-extsel:i+extsel]-polyval(idx,b[i-extsel:i+extsel]))**2))
    if extamp: 
            if extchi: return rep,arep,crep
            return rep,arep
    if extchi: return rep,crep
    return rep

curv_limit=[0.001,0.001]

def extrema(x,b=None,lbin=0,poly=0,ret_all=False,check_bins=True,msplit=False,retype='nd',
    polysel=None,extsel=None,ampsel=None,ampmeas=False,loud=1):
    '''gets points making an "envelope" of a function - maxima/minima of a periodic function
    
    requires rather smooth function - no noise fluctuation
        lbin:: rebinning function by bins of [lbin] points
            check_bins:: using extra.rebin with spline interpolation
    poly>0::
        fits a polynom of given order to the evelope
        ret_all:: not only fitted coefficients, but all points are returned
        polysel:: fits polynom only to reduced interval
    poly==0: 
        returns (minima,maxima)
    msplit:: split by the middle line
    extsel:: region around extrema fitted with quadrat. polyn.
    '''
    if b!=None: ob=b.copy()
    if lbin>1: #length of the bin
        nbin=len(x)//lbin
        if b!=None and check_bins: #is binning correct
            if (b[::lbin].std()>b[::lbin].mean()):
                print('binning not uniform: using spline interpolation')
                if b[0]>b[-1]:dir=-1
                else: dir=1
                try:
                    x,e,b=rebin(b,x,nbin,dir=dir,type=retype)
                except:
                    return b,x
                nbin=0
        if nbin>0:
            x=x[:nbin*lbin].reshape(nbin,lbin).mean(1)
            if b!=None:b=b[:nbin*lbin].reshape(nbin,lbin).mean(1)
            if loud>0: print('lbinning to %i bins'%len(x))
    else: lbin=1
    if b!=None and b[0]>b[-1]:
        b=b[::-1]
        x=x[::-1]
    from numpy import where
    pmins=where((x[1:-1]<x[2:])*(x[:-2]>x[1:-1]))[0]+1
    pmaxs=where((x[1:-1]>x[2:])*(x[:-2]<x[1:-1]))[0]+1
    #dif=x[1:]-x[:-1]
    #ale=where(dif[1:]*dif[:-1]<0)[0] 
    #ale - positions of all exetrema
    if len(pmins)<max(2,poly+1):
        print('only %i minima found in %i bins'%(len(pmins),len(x)))
        return [],[]
    #shi=0
    #if x[ale[0]]>x[ale[1]]: shi=1
    if loud>1: print('found %i/%i extrema'%(len(pmins),len(pmaxs)))
    
    from numpy import polyfit,polyval,arange,array
    if b==None: b=arange(len(x))
    if msplit: # split by the middle line
        mbin=len(b)//20 # should fit to some smoothed version - or simple mid-points [mbin//2::mbin]
        mline=polyfit(b[:mbin*20].reshape(20,mbin).mean(1),x[:mbin*20].reshape(20,mbin).mean(1),1) # middle line
        tops=x[pmaxs]>polyval(mline,b[pmaxs])
        bots=x[pmins]<polyval(mline,b[pmins])
        if loud>1: print('reduced to %i/%i extrema'%(sum(bots),sum(tops)))
        if msplit=='sort':
            from numpy import concatenate
            #ale=concatenate([b[pmaxs][tops],b[pmins][bots]])
            ale=concatenate([pmaxs[tops],pmins[bots]])
            isep=sum(tops)
            sale=ale.argsort() # how do maxima and minima mix
            pale=arange(len(ale))[sale<isep]
            sep=(pale[1:]-pale[:-1])>1 # where are maxima separated
            tep=[0]+list(arange(len(sep))[sep]+1)+[len(sep)+1]
            if loud>2: print('checking %i intervals'%len(tep))
            vals=x[pmaxs][tops]
            lpos=arange(len(tops))[tops]
            for i in range(len(tep)-1): 
                if tep[i+1]-tep[i]>1:
                    mpos=vals[tep[i]:tep[i+1]].argmax()
                    spos=list(lpos[tep[i]:tep[i+1]])
                    del spos[mpos]
                    tops[spos]=False # unmarking all other maxima
                    
            pale2=arange(len(ale))[sale>=isep]
            sep=(pale2[1:]-pale2[:-1])>1 # where are minima separated
            tep=[0]+list(arange(len(sep))[sep]+1)+[len(sep)+1]
            if loud>2: print('checking %i intervals'%len(tep))
            vals=x[pmins][bots]
            lpos=arange(len(bots))[bots]
            for i in range(len(tep)-1): 
                if tep[i+1]-tep[i]>1:
                    mpos=vals[tep[i]:tep[i+1]].argmin()
                    spos=list(lpos[tep[i]:tep[i+1]])
                    del spos[mpos]
                    bots[spos]=False # unmarking all other minima
            #qale=array([1]*sum(tops)+[0]*sum(bots))[sale] #check it
        atops=pmaxs[tops]
        abots=pmins[bots]
        if loud>1: print('finally found %i/%i extrema'%(sum(bots),sum(tops)))
    else:
        abots,atops=pmins,pmaxs
    if ampsel:
        ybots,ytops=[],[]
        sh=0 #shift
        if abots[0]<atops[0]: 
            sh=1
            ybots.append(x[atops[0]]-x[abots[0]])
        dh=len(atops)-len(abots)+sh
        for i in range(len(atops)-dh):
            d1=x[atops[i]]-x[abots[i+sh]]
            if sh==0 and i==0:
                ytops.append(d1)
            else:
                d2=x[atops[i]]-x[abots[i+sh-1]]
                ytops.append((d1+d2)/2.)
        if dh>0: 
            i=len(atops)-1
            ytops.append(x[atops[i]]-x[abots[i+sh-1]])
        for i in range(len(abots)-sh-1+dh):# min(len(abots),len(atops))-1):
            d1=x[atops[i]]-x[abots[i+sh]]
            d2=x[atops[i+1]]-x[abots[i+sh]]
            ybots.append((d1+d2)/2.)
        if dh==0:
            i=len(abots)-sh-1
            ybots.append(x[atops[i]]-x[abots[i+sh]])
        ytops=array(ytops)/max(ytops)
        ybots=array(ybots)/max(ybots)
        atops[ytops<ampsel]=0
        abots[ybots<ampsel]=0
        if loud>1: print('%i/%i extrema too weak'%(sum(abots==0),sum(atops==0)))
    if poly<0: return abots,atops
    if extsel!=None:
        if extsel<0: # reporting quality of the fit
            btops,mtops,ctops=ext_ran_fit(atops,x,b,-extsel,True,True,loud=loud)
            bbots,mbots,cbots=ext_ran_fit(abots,x,b,-extsel,True,True,loud=loud)
            if ampmeas: return bbots,btops,mbots,mtops,cbots,ctops
            return bbots,btops,cbots,ctops
        else: 
            btops,mtops=ext_ran_fit(atops,x,b,extsel,True,False,alim=-curv_limit[0],loud=loud)
            bbots,mbots=ext_ran_fit(abots,x,b,extsel,True,False,alim=curv_limit[1],loud=loud)            
            if ampmeas: return bbots,btops,mbots,mtops
            return bbots,btops
    if poly>0: # making polynomial fit
        if polysel:
            abots=abots[polysel[0]:polysel[1]]
            atops=atops[polysel[0]:polysel[1]]
        atops=atops[atops>0]
        abots=abots[abots>0]
        if len(abots)<poly or len(atops)<poly:
            print('too few points, cannot fit')
            return b[abots],b[atops]
        p1=polyfit(b[atops],x[atops],poly)
        p2=polyfit(b[abots],x[abots],poly)
        if ret_all: return (p2+p1)/2.,(p1-p2)/2.,b[abots],b[atops]
        else: return (p2+p1)/2.,(p1-p2)/2.
    return b[abots],b[atops]
    
global sin_fin,par_con,sin_inv_wei
sin_fin=None
par_con=[]

si_cau_mod=[ 0.0782, 0.0101] # term x^2 x^4 [ 0.03460822,  0.26845202,  3.43098771 ]

sin_inv_wei=8.5#for MIR data
sin_inv_wei=1.6#for simulated data

def sin_com(p,x,y=None,np=None,nord=2,inv=False):
    '''universal function for fitting periodic signals
    if np: phase depends like Cauchy on [x]
    inv=='frac': sin()/sin()/arctan([4])
    int=='else': 1/(3+sin(ph+[2]))
    '''
    from numpy import sin,arctan,polyval
    #q=p[:nord+1].copy()
    p[0]*=-p[0] #simple constraint, quadratic parameter should be negative
    a=polyval(p[:nord+1],x) #amplitude of periodic func.
    if (np!=None) and (len(np)>1): ph=x*(1+np[0]*x**2+np[1]*x**4)/p[nord+1]
    else: ph=x/p[nord+1]
    if inv=='frac':
        f=a*(sin(2*pi*(ph+p[nord+2]))+p[nord+5])
        f/=arctan(p[nord+4])/pi*2*sin(2*pi*(ph*p[nord+1]/np[2]+p[nord+2]+p[nord+3]))+1.1 
        # alternatively - using only positive values [0-1] (1-exp(-p[nord+4]))*0.9
        
    elif inv: 
        if (type(inv)==int) or (type(inv)==float): f=a*(inv-(inv+1)*(inv-1)/(inv+sin(2*pi*(ph+p[nord+2]))))
        else: f=a*(2-3/(2+sin(2*pi*(ph+p[nord+2]))))
    else: f=a*sin(2*pi*(ph+p[nord+2]))
    f+=polyval(p[-nord-1:],x) #baseline profile
    if y!=None: f-=y
    return f

global last_pars
last_pars=None

def fitting(x,y,p0=None,prange=None,bounds=None,loud=1,fit_mode=1,refr_mode=None,lbin=0,nprof=None,nord=2):
    '''periodic part fitted with slightly slanted sinusoid
    prange : range of periods to use..
    // added possible polynomial dependence of refraction index (oscilation frequency) on wavelength
        refr_mode='poly' : nprof gives 2 Cauchy profile parameters
    nord: order of polynom of baseline/amplitude fit
    nprof: profile of index of refraction [parameters of Cauchy]
    returns parameters (ampl.polynom + period + phase-shift + bckg. polynom [+ 2-par Cauchy profile]) and chi2
    '''
    global sin_fin,par_con,last_pars
    if fit_mode<0: fit_mode=glob_fit_mode
    from numpy import polyval,array,sign,ndarray,iterable
    if (sin_fin==None) or nprof:
        if refr_mode=='poly':
            if iterable(nprof)>0:#type(nprof)==list or type(nprof)==array or type(nprof)==ndarray: 
                sin_fin=lambda p,x,y:sin_com(p,x,y,np=nprof,nord=nord)
                #sin_fin=lambda p,x,y:((p[0]*x*x+p[1]*x+p[2])*sin(2*pi*(x/p[3]*(1+nprof[0]*x**2+nprof[1]*x**4)+p[4]))+p[5]*x*x+p[6]*x+p[7]-y)
                if loud>2: print('defined 7par function')
            else:
                sin_fin=lambda p,x,y:sin_com(p[:4+2*nord],x,y,np=p[4+2*nord:],nord=nord)#+(1-sign(p[0]))*5.
            if type(bounds)==list and len(bounds)==0:
                for i in range(len(bounds),4+2*nord):
                    bounds.append(None)
                    bounds.append([0,1000.])
                    bounds.append([0,1000.])
            #sin_fin=lambda p,x,y:((p[0]*x*x+p[1]*x+p[2])*sin(2*pi*(x/p[3]*(1+p[8]*x**2+p[9]*x**4)+p[4]))+p[5]*x*x+p[6]*x+p[7]-y)
        elif refr_mode=='frac': sin_fin=lambda p,x,y:sin_com(p,x,y,nord=nord,np=nprof,inv='frac')
        elif refr_mode=='inve': 
            if iterable(nprof)>0: 
                sin_fin=lambda p,x,y:sin_com(p[:4+2*nord],x,y,nord=nord,inv=sin_inv_wei,np=nprof)
            else: 
                if nprof!=None: sin_fin=lambda p,x,y:sin_com(p[:4+2*nord],x,y,nord=nord,inv=nprof,np=p[4+2*nord:])
                else: sin_fin=lambda p,x,y:sin_com(p[:4+2*nord],x,y,nord=nord,inv=sin_inv_wei,np=p[4+2*nord:])
        else: 
            sin_fin=lambda p,x,y:sin_com(p,x,y,nord=nord)#+(1-sign(p[0]))*5.
        #if bounds==None: bounds=[[0,1e6]] 
            #if nord==3: sin_fin=lambda p,x,y:((p[0]*x*x*x+p[1]*x*x+p[2]*x+p[3])*sin(2*pi*(x/p[4]+p[5]))+p[6]*x*x*x+p[7]*x*x+p[8]*x+p[9]-y)
            #elif nord==2: sin_fin=lambda p,x,y:((p[0]*x*x+p[1]*x+p[2])*sin(2*pi*(x/p[3]+p[4]))+p[5]*x*x+p[6]*x+p[7]-y)
            #else: sin_fin=lambda p,x,y:((p[0]*x+p[1])*sin(2*pi*(x/p[2]+p[3]))+p[4]*x+p[5]-y)
    if p0==[]: return
    if p0=='last': p0=last_pars
    elif p0==None:
        if lbin<0 and len(x)>200: lbin=20
        try:
            am,ad,ip,iq=extrema(y,x,poly=nord,all=True,loud=2,lbin=lbin,check_bins='full',msplit='sort')
        except:
            return [array([0.]*(nord+6)),10000]
        phas=0
        if am[0]>0: 
            am=-am
            phas+=0.5
        from math import sqrt
        am[0]=-sqrt(-am[0])
        p0=list(am)+[0,phas] #period + phase
        if refr_mode=='frac':p0+=[0,phas+0.5,1.] #phase and relat. intensity (<1.) of nominator
        p0+=list(ad)
        if refr_mode=='poly':p0+=[0,0]
        if loud>0:
            print('mean value %f, mean amplitude %f'%(polyval(am,x).mean(),polyval(ad,x).mean()))
        dp=(ip[1:]-ip[:-1]).mean()
        dq=(iq[1:]-iq[:-1]).mean()
        if abs(dp-dq)>(iq[1:]-iq[:-1]).std():
            print('unstable period %f vs %f'%(dp,dq))
        else: print('estimated period %f'%((dp+dq)/2))
        ppos=nord+1 # which parameter is period?
        p0[ppos]=(dp+dq)/2
        if prange!=None:
            if abs(p0[ppos])<prange[0]:p0[ppos]=prange[0]
            elif abs(p0[ppos])>prange[1]:p0[ppos]=prange[1]
            else: p0[ppos]=abs(p0[ppos])
        if loud>2: print('init pars: %s'%p0)
        #p0[4]=iq[0]/dq-.25)
    chi2=lambda p,x,y:(sin_fin(p,x,y)**2).sum()
    if bounds: # constrained fit
        from scipy.optimize import fmin_cobyla as fmin
        if type(bounds)==list:
            par_con=[]
            for i in range(len(bounds)):
                if bounds[i]!=None:
                    par_con.append(lambda p,a,b:(p[i]-bounds[i][0])*(bounds[i][1]-p[i]))
        else:
            print('dont know how to interpret bounds "%s"'%str(bounds))
    else:
        fun=chi2
        if fit_mode==1: from scipy.optimize import fmin_bfgs as fmin
        elif fit_mode==2: from scipy.optimize import fmin_ncg as fmin
        else: 
            from scipy.optimize import leastsq as fmin
            fun=sin_fin
    if bounds: 
        print('constrained call with 2 bounds')
        last_pars=fmin(sin_fin,p0,par_con,args=(x,y),consargs=None)
    elif fit_mode>0: last_pars=fmin(fun,p0,extra_der,args=(x,y),disp=loud)
    else: last_pars,ok=fmin(fun,p0,args=(x,y))
    #if ok!=1: return zoo
    return last_pars,chi2(last_pars,x,y)

def extra_fun(y,per,exten=0,p=None):
    '''some more complicated functions describing periodic spectrum/signal
    '''
    global sin_fin,par_con
    from numpy import sin
    if exten==1: 
        sin_fin=lambda p,x,y:p[0]*(sin(2*pi*(x/p[1]+p[2]))+p[6])/(sin(2*pi*(x/p[3]+p[4]))*p[7]+1)+p[5]*x-y+p[8]*x*x
    else: sin_fin=lambda p,x,y:(p[0]*(sin(2*pi*(x/p[1]+p[2]))+p[6])/(sin(2*pi*(x/p[3]+p[4]))*p[7]+1)+p[5]*x-y)
    if p==None:
        p=[0 for a in range(8)]
        p[1]=per
        p[3]=per*1.02
        p[6]=y.mean()
        p[0]=abs(y-p[6]).max()
        p[6]/=p[0]
        p[7]=0.1
    par_con+=[lambda p:(-p[7]+0.9)*(p[7]+0.9)]
    for i in [1,3]:
        par_con+=[lambda p:(-p[i]+per*0.5)*(p[i]-1.7*per)]
    if exten==1: p.append(0.)
    return p

def extra_der(p,x,y,exten=0):
    '''calculates analytic derivatives
    '''
    from numpy import sin,cos,array
    osc=sin(2*pi*(x/p[3]+p[4]))
    osc2=2*pi*cos(2*pi*(x/p[3]+p[4])) #cosine
    amp=(p[0]*x*x+p[1]*x+p[2])
    if sin_fin==None: return None
    f=sin_fin(p,x,y)
    rep=[(x*x*f*osc).sum(),(x*f*osc).sum(),(f*osc).sum()]
    famp=amp*osc2*f
    rep+=[-(famp*x/p[3]**2).sum(),famp.sum()] #period and phase
    rep+=[(x*x*f).sum(),(x*f).sum(),f.sum()] #baseline
    if len(p)>8:
        rep+=[(famp/p[3]*x**3).sum(),(famp/p[3]*x**5).sum()]
    return 2*array(rep)

def num_der(x,p0,dp0=None):
    '''calculating (relative) numerical derivatives
    returns right and left limits
    '''
    if dp0==None: dp0=p0*0.1
    rep=[]
    from numpy import arange,array
    y=sin_fin(p0,x,None)
    ap=arange(len(p0))
    rep=[[(dp0[i]>0) and (sin_fin(p0+(ap==i).astype(int)*dp0[i],x,y)).sum() or 0. for i in ap]]
    rep+=[[(dp0[i]>0) and (sin_fin(p0-(ap==i).astype(int)*dp0[i],x,y)).sum() or 0. for i in ap]]
    return array(rep)[:,dp0>0]/dp0[dp0>0]/len(x)

def num_der2(x,p0,dp0):
    #hessian matrix
    from numpy import arange,array,dot
    dp0=array(dp0)
    d0=num_der(x,p0,dp0)
    ap=arange(len(p0))
    ret=array([1/2.,-1/2.]).reshape(1,2)
    rep=array([dot(ret,num_der(x,p0+(ap==i).astype(int)*dp0[i],dp0)) for i in ap])
    return rep/dp0

def harmony(x,per,wharm=None):
    '''creating fast rise/fall oscilations by adding higher harmonics
    '''
    from numpy import array,arange,dot,sin,newaxis
    if wharm==None: wharm=array([1,.5,.2,.1,.05])
    nharm=len(wharm)
    base=sin(arange(nharm)[:,newaxis]*x[newaxis,:]/per*2*pi)
    return dot(wharm,base)

def period(x,y,per=0,lbin=20,nbin=20,npha=10,corr=True):    
    ''' probably pulses are asymetric, something like
    log(sin(x/20.+.2)+1.1)
    '''
    global par
    from numpy import polyfit,polyval,array,zeros,sin,pi
    from scipy.optimize import leastsq
    #from scipy.optimize import lbfgsb
    if par==None or len(par)<6:par=zeros(6)
    fun = lambda par:(par[0]*x+par[3])*sin(2*pi*(par[1]*x-par[2]))+par[4]*x+par[5]-y
    if per==0:
        rep=filter(y,x,min_ampl=0,rep=[])
        if len(rep)==0: 
            print('cannot pre-determine period - exiting')
            return
        else:
            if len(rep)>1:
                print('%i periods found (%s)'%(len(rep),rep))
            per=abs(rep[0])
    else: par[1]=per
    if len(x)<lbin*nbin:
        nbin=len(x)//lbin
        print('reducing to %i bins'%(nbin))
    if len(x)<lbin*nbin:
        nbin=len(x)//lbin
        print('reducing to %i bins'%(nbin))
    pos=x[:lbin*nbin].reshape(lbin,nbin).mean(1)
    pom=y[:lbin*nbin].reshape(lbin,nbin).mean(1)
    a,b=polyfit(pos,pom,1)
    par[4]=a
    par[5]=b
    pom=y[:lbin*nbin].reshape(lbin,nbin).std(1)
    a,b=polyfit(pos,pom,1)
    par[0]=a*1.4
    par[3]=b
    out=leastsq(fun,par)
    print('found period %.3f, phase %.3f'%(1./out[0][1],out[0][2]))
    if npha==0: return 1./out[0][1],out[0][2]
    pha=x*out[0][1]-out[0][2]
    pia=((pha-pha.astype(int))*npha).astype(int)
    if corr:
        z=y-polyval(par[4:6],x)
        z/=polyval(par[[0,3]],x)
    else:
        z=y
    prof=[z[pia==i].mean() for i in range(npha)]
    perr=[z[pia==i].std() for i in range(npha)]
    return array(prof),array(perr)
    
        
def crossnorm(bas1,dat1,bas2,dat2,prec=0.1,bkg1=None,bkg2=None,rat_frac=0.2,con_siz=1,loud=0,ord=0):
    '''crossnormalisation of 2 spectra 
    with some overlapping region
    if background rates are supplied - considers only bins with rates above _rat_frac_*average
    when binning is different ("prec" is precision for bin position):
        resamples to rougher binning
        if con_siz>0: sliding average
    '''
    from numpy import arange,hamming,convolve,polyfit,array
    lim2=sum(bas2>bas1[-1])
    if lim2==0:
        bas1,dat1,bas2,dat2=bas2,dat2,bas1,dat1
        lim2=sum(bas2>bas1[-1])
        if lim2==0:
            print('no overlap')
            return
    lim1=sum(bas2[0]>bas1)
    print('using %i vs. %i bins'%(lim1,lim2))
    if lim2<lim1:
        alim=arange(lim1)
        idis=[alim[abs(bas1[-lim1:]-bas2[i])<prec] for i in range(lim2)]
        idis=array([len(a)>0 and a[0] or -1 for a in idis])
        if loud>0: print('matching %i indices of %i'%(sum(idis>=0),len(idis)))
        #idis=array(idis)[:,0]
        if con_siz>0:
            n=int(lim1/lim2)+con_siz
            oof=hamming(n+2)[1:n+1]
            oof/=sum(oof)
            dat0=convolve(dat1[-lim1:],oof,'same')
            if loud>1: print('mean value in 1. set:from %.3f to %.3f'%(dat1[-lim1:].mean(),dat0.mean()))
        else:
            dat0=dat1[-lim1:]
        dat0=dat0[idis][1:-1]
        dat3=dat2[1:lim2-1]
        if bkg1!=None: bkg1=bkg1[idis][1:-1]
        if bkg2!=None: bkg2=bkg2[1:lim2-1]
    elif lim2>lim1:
        print('please put finer sampled data first')
        return
    else:
        dat0=dat1[-lim1:]
        dat3=dat2[:lim2]
        if bkg1!=None: bkg1=bkg1[-lim1:]
        if bkg2!=None: bkg2=bkg2[:lim2]
    sele=None
    if bkg1!=None:sele=bkg1>bkg1.mean()*rat_frac
    if bkg2!=None:
        if sele!=None: sele*=bkg2>bkg2.mean()*rat_frac
        else: sele=bkg2>bkg2.mean()*rat_frac
    if sele!=None:
        dat0=dat0[sele]
        dat3=dat3[sele]
        ls=sum(sele)
    else:
        ls=0
    if loud>0: print('fitting %i vs %i bins'%(len(dat0),len(dat3)))
    if ord>0:
        idx=polyfit(dat0,dat3,ord)
        print('overlay %i bins (%i):dat2=%.3f + dat1*(%.3f)'%(lim2-2,ls,idx[1],idx[0]))
        return (dat2-idx[1])/idx[0]
    if ord<0:
        rat=(dat3/dat0).mean()
        print('ratio %.3f'%rat)
        return dat2/rat
    dif=(dat3-dat0).mean()
    print('diff %.3f'%dif)
    return dat2-dif

global bord
bord=None
def crossfit(lims,bas1,dat1,bas2,dat2,bkg1=None,bkg2=None,rat_frac=0.2,ord=1,loud=1,rep=0,cross=[0,-1],step=0):
    '''cross-callibration using polynomial fits of data
    '''
    global bord
    from numpy import arange,median,polyfit,polyval,dot
    sele1=(bas1>lims[0])&(bas1<lims[1])
    if bkg1!=None: sele1&=(bkg1>(median(bkg1)*rat_frac))[:len(sele1)]
    if sum(sele1)==0:
        print('no bins matching for set1 [%.1f - %.1f]'%tuple(lims[:2]))
        return
    idx1=polyfit(bas1[sele1],dat1[sele1],ord)
    chi1=sum((dat1[sele1]-polyval(idx1,bas1[sele1]))**2)/sum(dat1[sele1]**2)
    if loud>0: print("fitting %i bins [chi %.3g]: %s"%(sum(sele1),chi1,idx1))
    if len(lims)>2: lims2=lims[2:]
    else: lims2=lims
    sele2=(bas2>lims2[0])&(bas2<lims2[1])
    if bkg2!=None: sele2&=(bkg2>(median(bkg2)*rat_frac))[:len(sele2)]
    if sum(sele2)==0:
        print('no bins matching for set2 [%.1f - %.1f]'%tuple(lims2))
        return
    idx2=polyfit(bas2[sele2],dat2[sele2],ord)
    chi2=sum((dat2[sele2]-polyval(idx2,bas2[sele2]))**2)/sum(dat2[sele2]**2)
    if loud>0: print("fitting %i bins [chi %.3g]: %s"%(sum(sele2),chi2,idx2))
    if rep==0: return idx1,idx2
    if step>0:
        x=arange(lims[cross[0]],lims[cross[1]],step)
        rat=polyval(idx1,x).sum()
        rat/=polyval(idx2,x).sum()
    else:
        # definite integral btwin lims[0] and lims[-1]
        bord=[pow(float(lims[i]),arange(1,ord+2)[::-1])/arange(1,ord+2)[::-1] for i in cross]
        rat=dot(idx1*idx2,bord[1])-dot(idx1*idx2,bord[0])
        rat/=dot(idx2**2,bord[1])-dot(idx2**2,bord[0])
    return rat
    
def interfero(theta,lamb=0.5,R=0.3,n=1.5,d=10.,unit='deg'):
    ''' lamb and d in micrometers
    '''
    from numpy import sin, cos, pi
    
    if unit=='deg': theta*=180./pi
    fi=4*pi*n*d*cos(theta)/lamb
    a=4*R*sin(fi/2.)**2
    return a/((1-R**2)+a)
    
#----------------------------------------------------------------------------    
 
def unittest(nsamp=100,size=10.,freq=[8],rep=0,loud=0,relamp=0.2):
    from numpy import linspace,exp,sin,log
    from math import pi
    aw2=linspace(-size,size,nsamp)
    do=exp(-aw2**2/2.)
    do+=.3*exp(-(aw2-6)**2/4.)
    do+=.3*exp(-(aw2+4)**2/1.)
    do2=do.copy()
    for f in freq:
        if f<0: do2+=relamp*(log(sin(-2*pi*f*aw2)+2)-log(3.)/2) # logsinus
        else: do2+=relamp*sin(2*pi*f*aw2)
    if loud>0: print('mean %.4f before, %.4f after'%(do.mean(), do2.mean()))
    if rep==1: return do2
    redo=filter(do2,aw2,min_wid=0,min_ampl=0)
    return redo
