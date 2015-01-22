# -*- coding: utf-8 -*-
'''increasingly complex library for
dielectric function modelling (lorentzian oscilators), mixing, (multi)layers:
    dielect
fitting of reflectivity on multilayers / gradual / mixing (EMA) materials:
    reflect
    plane
    gradient
spectral 
plotting of complex functions (1 or 2 panes)
'''


from math import pi

def dielect(freq,einf,osci,etauc=None,drude=None):#,n0=1.,rep=-1,ang=0,polar='s'):
    '''calculates dielectric function with lorentz oscilators
    
    rep==-1: returns dielectic function
    rep==0: returns complex reflectivity (under angle [ang] from the normal)
    rep==1: returns abs. value of reflectivity
    
    drude: amplit and inverse collision time of free carriers
    osci: list of tuples (ampl,freq,absorb)
    n0: refr. index of external material 
    '''
    from numpy import zeros
    e=zeros(len(freq),dtype='complex128')
    for o in osci:
        if len(o)<3:
            print('3rd parameter should be a list of 3-tuples')
            return
        e+=o[0]*o[1]*o[2]/(o[1]**2-freq**2-1j*o[2]*freq)
    if etauc: e*=(freq-etauc)**2
    if einf: e+=einf
    if drude!=None:  e-=drude[0]**2/freq/(freq+1j*drude[1])
    return e

def dielect_tauc(freq,einf,ampl,e0,egap,cbri,rep=0):
    '''tauc-lorentz form for single oscillator with gap
    complicated analytical form - hopefully correct
    from Jellison Modine APL 69,2137 (1996)  + erratum
    '''
    from numpy import zeros,sqrt,arctan,abs,log
    e=zeros(len(freq),dtype='complex128')
    alpha=sqrt(4*e0**2-cbri**2)
    mgam=sqrt(e0**2-cbri**2/2)
    a_ln=freq**2*(egap**2-e0**2)+egap**2*cbri**2-e0**2*(3*egap**2+e0**2)
    a_atan=(freq**2-e0**2)*(egap**2+e0**2)+egap**2*cbri**2
    ksi4=(freq**2-mgam**2)**2+alpha**2*cbri**2/4

    brac1=(egap**2+e0**2+alpha*egap)/(egap**2+e0**2-alpha*egap)
    brac2=pi-arctan((alpha+2*egap)/cbri)+arctan((alpha-2*egap)/cbri)
    brac3=pi+2*arctan(2*(mgam**2-egap**2)/alpha/cbri)
    brac5=abs(freq-egap)*(freq+egap)/sqrt((egap**2-e0**2)**2+egap**2*cbri**2)
    #ear=[a_ln/2/alpha/e0*log(brac1),-a_atan/e0**2/cbri*brac2]
    ear=[a_ln/2/alpha/e0**2*log(brac1),-a_atan/e0**2/cbri*brac2]
    #ear+=[2*egap*(freq**2-mgam**2)/alpha/cbri*brac3]
    ear+=[4*egap*(freq**2-mgam**2)/alpha/cbri*brac3]
    ear+=[-(freq**2+egap**2)/freq*log(abs(freq-egap)/(freq+egap)),2*egap*log(brac5)]
    if rep==-1: return ear
    #common factor
    e.real=sum(ear,0)*ampl*e0*cbri/pi/ksi4
    sfreq=freq[freq>egap]
    e[freq>egap]+=1j*(ampl*e0*cbri*(sfreq-egap)**2/((sfreq**2-e0**2)**2+cbri**2*sfreq**2)/sfreq)
    e+=einf
    return e

def rdielect(freq,einf,osci,drude=None,n0=1.,ang=0,polar='s',rep=1):
    # reflection from bulk with given lorentz oscillators
    e=dielect(freq,einf,osci,drude)
    return reflect(e,n0,rep,ang,polar)

global wid_explo    
con_bins=4.
cen_toler=10
wid_explo=2.
pare=[-0.34518248, 1.4077] # log(amplit.) decrease with "gaus" param.
def profiles(x,lore,moment=0,loud=0,recent=False,fwhm=False):
    '''Voigtian modelling
    "lore" gives 3 parameters of loretzian profile
    4th parameter is gaussian broadening
    '''
    from math import pi,sqrt
    from numpy import abs,median
    step=median(x[1:]-x[:-1])
    y=lore[0]/(1-x**2/lore[1]**2-1j*x*lore[2]/lore[1])
    y=abs(y)
    if len(lore)==4 and lore[3]>0.:
        gaus=lore[3]
        from numpy import convolve,arange,exp
        z=arange(-gaus*con_bins,gaus*con_bins,step)
        y=convolve(y,exp(-z**2/gaus**2/2.)/sqrt(2*pi)*step/gaus)
        p=(len(y)-len(x))
        if p%2==0: y=y[p//2:-p//2]
        else: 
            y=(y[p//2:-p//2]+y[p//2+1:-p//2+1])/2.
    #if moment<0: return y
    if moment!=0:
        s0=y.sum()
        rep=[s0]
        if recent: 
            x0=(x*y).sum()/s0
            z=x-x0
            v=x-x0
            if loud>0: print('centered %.2f '%x0)
        else:
            z=x.copy()
            v=x.copy()
        if fwhm:
            cent=len(x)//2
            hmax=y[cent-cen_toler:cent+cen_toler].max()/2.
            p1=sum(y[:cent]<hmax)
            p2=sum(y[cent:]<hmax)
            fwid=x[-p2-1]-x[p1]
            if wid_explo!=1.:
                p1=cent-(cent-p1)*wid_explo
                p2=cent-(cent-p2)*wid_explo
                
            if loud>0: print("width %.2f [%i bins]"%(fwid,len(x)-p2-p1))
            rep.append(fwid)
            y=y[p1:-p2]
            z=z[p1:-p2]
            v=v[p1:-p2]
        for i in range(2,abs(moment)):
            z*=v
            rep.append((z*y).sum()/s0)
            #wei=pow(x-x0,arange(1,moments)) 
        return rep
    return y

def comp_profile(z,rang=[-100,100]):
    '''
    see http://www.physics.muni.cz/ufkl/Publications/clanky/JH_folia84.pdf
    adapted for complex z
    first fit [ -0.34492182, 10.2065963 , 3.33501317, 3.08889307, 25.43125216,  18.70707501]
    fun=lambda p:(p[0]/(p[1]-1j*p[2]-z1)+p[3]/(p[4]-1j*p[5]-z1))
    '''
    from scipy import integrate
    from numpy import exp
    fr=lambda t:exp(-t**2)*(z.real-t)/((z.real-t)**2+z.imag**2)
    fi=lambda t:exp(-t**2)*z.imag/((z.real-t)**2+z.imag**2) #*(-1)
    w=1/pi*(1j*integrate.quad(fr,rang[0],rang[1])[0]+integrate.quad(fi,rang[0],rang[1])[0])
    return w

######### MIXING METHODS ##########################

def garnett(e_i,e_m,delt=0.1,depol=1):
    '''
    spherical/elipsoidal inclusions in matrix
    Maxwell-Garnett model
    '''
    e_f=e_m*(e_i*(1 + 2*delt) - e_m*(2*delt - 2))
    #e_f/=e_m-(e_m-e_i)*depol*(1 - delt)
    e_f/=e_m*(2+delt)+e_i*depol*(1 - delt)
    return e_f 
def brugemann(e_i,e_m,delt=0.1):
    #solve(delt*(e_i-e_f)/(2*e_f+e_i)-(delt-1)*(e_m-e_f)/(2*e_f+e_m))
    from numpy import sqrt
    def qsolve(a,b,c):
        return (-b+sqrt(b**2-4*a*c))/2/a
    return qsolve(2,e_i-2*e_m+delt*(e_m-e_i)*3,-e_m*e_i)
    e_f=e_m/2 - e_i/4 - 3*delt*(e_m+e_i)/4
    e_f+= sqrt((4+18*delt)*e_i*e_m + (1- 6*delt)*e_i**2 + (4 -12*delt)*e_m**2  + 9*delt**2*e_i**2 + 9*delt**2*e_m**2 - 18*e_i*e_m*delt**2)/4
    return e_f
def ll(e_i,e_m,delt=0.1):
    """

    :param e_i: inclusion dielectric
    :param e_m: bulk dielectric
    :param delt: fraction of inclusion
    :return:
    """
    from numpy import pow
    return pow(delt*pow(e_i,1/3.)+(1-delt)*pow(e_m,1/3.),3.)
def cpa(e_i,e_m,delt=0.1):
    '''general form
    solve((e_f-e_m)*(3*e_f-(1-delt)*(e_i-e_m))-3*e_f*delt*(e_i-e_m))
    '''
    from numpy import sqrt
    e_f=e_m/3 + e_i/6 - delt*(e_m+e_i)/3 
    e_f+=sqrt(-8*e_i*e_m + 16*delt*e_i*e_m + e_i**2 + 16*e_m**2 - 20*delt*e_m**2 + 4*delt*e_i**2 + 4*delt**2*e_i**2 + 4*delt**2*e_m**2 - 8*e_i*e_m*delt**2)/6
    return e_f

############ reflectivity & comp. of layers ##############

def reflect(e,n0=1.,rep=-1,ang=0,polar='s'):
    '''calculates reflectivity under angle ang and polarization 'r'/'s'
        angle in degrees
        if rep==0: returns complex
    '''
    from numpy import sqrt,cos,all,conj
    #from numpy import conj,sqrt,cos,ones,all
    n=sqrt(e)
    if all(ang==0): 
        #print 'normal'
        r=(n0-n)/(n0+n) # (1-n^2-k^2- 2ik)/((1+n)^2+k^2)
    else:
        cang=cos(ang*pi/180)
        canh=sqrt(1-n0**2*(1-cang**2)/n**2) # snell's law
        if rep==-2: return canh
        if polar=='s':r=(n0*cang-n*canh)/(n0*cang+n*canh) # fresnel law
        else: r=(n*cang-n0*canh)/(n0*canh+n*cang) # perpend. polarization
    if rep==0: return r
    return (r*conj(r)).real #abs(r)**2

def friter(r,sh,psi=None,ang=0,aver=False):
    '''helper function - calculating global reflectivity from 
    reflectivities of individual interfaces [r] and phase shifts at layers [sh]
    '''
    from numpy import ones,exp,conj
    if psi==None: psi=ones(len(r))
    if aver: # no coherence = no interference
        if ang==0: 
            ref=(r[0]*conj(r[0]))
            rep=(1-ref)**2*(1-n.imag**2/n.real**2)*exp(-2*sh[0])
            rep/=1-ref**2*exp(-2*sh[0])
            rep=ref*(1+rep)
        return rep
    else: 
        #phas=2*alpha*width*cos(ang)
        for i in range(len(r)-1,0,-1): #backtracking
            r[i-1]=(r[i-1]+r[i]*exp(2j*sh[i-1]*psi[i-1]))/(1+r[i-1]*r[i]*exp(2j*sh[i-1]*psi[i-1]))
        return r[0]
        #rep=(r[0]*conj(r[0])).real #is already a real number


global n,k,alpha,psi
psi=[]
def plate(freq,epsil,width,ang=0,polar='s',n0=1.,rep=1,aver=False,unit='eV'):
    '''response from (a) planparalel plate(s)
    frequency in eV (default), cm-1 or in PHz (1e15 Hz)
    width in nm (1e-9m)
    widths and epsils will be lists of values (or arrays) for each layer
    angle in degrees
    meth=1: using tensor classes

    if rep=-1: returns precalculated reflectivites on interfaces and phase shifts
            if ang>0: also returns "psi" (reflected angles)
    if rep=0: returns complex (fresnel coef.)
    '''
    global n,k,alpha,psi
    from numpy import abs,sqrt,cos,array,zeros,ndarray
    if type(epsil) not in [list,array,ndarray]: epsil=[epsil]
    if type(width) not in [list,array,ndarray]: width=[width]

    r=[] #fresnel coef. at each interface
    sh=[] #phase shift in each layer
    if len(width)==len(epsil):nend=n0 #exit layer same as input
    else:nend=0
    from spectra import c0,ev2um
    cang,canh=None,None # cosine of angles above/below the interface
    psi=[]
    i=0
    if unit=='PHz': mfrq=freq/c0*1e6
    elif unit=='cm': mfrq=1e-7*freq 
    else: mfrq=(freq*ev2um*1e-3)
    for e in epsil: # calculate fresnel coef., widths and angles at each interface
        n=sqrt(e)
        if i<len(width): sh.append(mfrq*2*pi*n*width[i])
        if ang!=0:
            if cang==None: cang=cos(ang*pi/180.)
            else: cang=canh # output angle from previous layer
            psi.append(cang)
            #ang below - check the total reflection
            canh=sqrt(1-n0**2*(1-cang**2)/n**2)
            # bang=(n0.real)**2*(1-cang**2)/(n.real)**2
            # if all(bang<1):
                # canh=sqrt(1-bang) #cos of angle below the interface
            # else: break # ERR: not correct - some wavelengths will pass
            if polar=='p': rs=(n*cang-n0*canh)/(n*cang+n0*canh) # in fact rs is rp, just to save one variable
            else: rs=(n0*cang-n*canh)/(n0*cang+n*canh)
            if ang<0: 
                rp=(n*cang-n0*canh)/(n*cang+n0*canh)
                r.append(array([rp,rs]).transpose())
            else: r.append(rs) # fresnel law
        else: 
            r.append((n0-n)/(n0+n))
            psi.append(1)
        n0=n
        i+=1
    if ang!=0: psi.append(canh)
    if nend>0: #backside scattering 
        if ang!=0: 
            cang=canh
            canh=sqrt(n0**2-n**2*(1-cang**2))
            if polar=='p': rs=(nend*cang-n*canh)/(nend*cang+n*canh)
            else:rs=(n*cang-nend*canh)/(n*cang+nend*canh)
            if ang<0:
                rp=(nend*cang-n*canh)/(nend*cang+n*canh)
                r.append(array([rp,rs]).transpose())
            else: r.append(rs)
            psi.append(canh)
        else: r.append((n0-nend)/(n0+nend))
        width[0]=-width[0]
    if rep==-1: 
        if ang==0: return r,sh,None
        else: return r,sh,psi
    elif rep==0: return friter(r,sh,psi,ang=ang,aver=aver)
    return abs(friter(r,sh,psi,ang=ang,aver=aver))**2

global delt,pl

def tensor_plate(freq,epsil,width,rep=0,meth=1,n0=1.,ang=0):
    '''using tensor algebra
    '''
    from numpy import array,ones,zeros,sqrt,conj
    from math import cos,sin
    from algebra import tensor
    sp=epsil[0].shape
    matt=tensor(array([[ones(sp),zeros(sp)],[zeros(sp),ones(sp)]],'c32'),dim=2)
    i=0
    r=[]
    sh=[]
    cang,canh=None,None
    for e in epsil: # calculate fresnel coef., widths and angles at each interface
        n=sqrt(e)
        if i<len(width): sh.append(mfrq*2*pi*n*width[i])
        if ang!=0:
            if cang==None: cang=cos(ang*pi/180)
            else: cang=canh # output angle from previous layer
            psi.append(cang)
            canh=sqrt(1-(n0.real)**2*(1-cang**2)/(n.real)**2) #cos of angle below the interface
            pl=n*cang
            delt=sh[-1]*cang
            selt,celt=sin(delt),cos(delt)
            matt*=tensor([[celt,-1j/pl*selt],[-1j*pl*selt,celt]])
        else:
            selt,celt=sin(sh[-1]),cos(sh[-1])
            pl=1
            matt*=tensor([[celt,-1j*selt],[-1j*selt,celt]])
        n0=n
        i+=1

    if meth==3: return matt
    a=(matt[0,0]+matt[0,1]*pl)*n0*cos(ang*pi/180)
    b=(matt[1,0]+matt[1,1]*pl)
    r=(a-b)/(a+b)
    return (r*conj(r)).real
    

def matter_plate(freq,epsil,width,rep=0,meth=0,polar='s',ang=0,n0=1):
    '''alternative approach
        without angular dependance yet
    '''
    from numpy import sqrt,cos,sin,exp,conj,zeros,ones
    from spectra import ev2um
    global delt,pl
    matt=None;#matrix([[1.,0],[0.,1.]],'c16')
    cang=1
    if polar=='p': pl0=n0/cang
    else: pl0=n0*cang
    o1=ones(len(freq))
    z1=zeros(len(freq))
    from algebra import tensor
    for i in range(len(width)):
        n=sqrt(epsil[i])
        pl=n*cang
        delt=2*pi*(freq*ev2um*1e-3)*pl*width[i]
        if polar=='p': pl=n/cang
        if meth==1:
            mult=tensor([[cos(delt),-1j/pl*sin(delt)],[-1j*pl*sin(delt),cos(delt)]])
        else:
            r12=(pl0-pl)/(pl0+pl)
            t12=2*pl0/(pl0+pl)
            mult=tensor([[o1,r12],[r12,o1]])*(1./t12)
            if width[i]>0:
                delt=2*pi*(freq*ev2um*1e-3)*pl*width[i]
                mult*=tensor([[exp(1j*delt),z1],[z1,exp(-1j*delt)]])
        if matt==None: matt=mult
        else: matt*=mult
        #cang=sqrt(pl0**2-pl**2*(1-cang**2))
        pl0=pl
    #now pl0 and pl are that of the incident medium and substrate, resp.
    if rep==-1: return matt
    if len(epsil)>len(width): pl=sqrt(epsil[len(width)])*cang
    a=(matt[0,0]+matt[0,1]*pl)*sqrt(epsil[0])*cang
    b=(matt[1,0]+matt[1,1]*pl)
    r=(a-b)/(a+b)
    return (r*conj(r))
'''parametrization of profiles

ampl. set to 1.
without gaussian
profile(x,[1,a,b,0])

fwhm=a+p1(a)*exp(-p2(a)*(b-a))
surf=integ(-2*fwhm,+2*fwhm)
    =q0(a)+q1(a)*exp(-q2(a)*(b-a))
    (alternatively a log-quadratic)
    
p2=0.47-0.50
p1=8.-9.

[p2,p1] pro a = [1,1.5,2.,2.5,3.,3.5]
[ 0.46965852,  5.73013826],
[ 0.48068463,  8.20836805],
[ 0.4787242 ,  8.8263208 ],
[ 0.48346287,  8.66321754],
[ 0.50532908,  8.33762317],
[ 0.55771015,  8.22010904]]

then come 2nd,4th .. moments
'''    

def gradient(epsil,grange=[0.,1.],ndiv=30,meth=garnett,freq=None,width=1000):
    '''create multilayer with gradual mixing of 2 materials (given in epsil)
    '''
    from numpy import linspace,ones
    lays=[meth(epsil[0],epsil[1],a) for a in linspace(grange[0],grange[1],ndiv)]
    ways=ones(ndiv)*width/ndiv
    if freq==None: return lays,ways
    return plate(freq,lays,list(ways)[:-1],rep=0)
    
def calc_ellips(epsil,n0=1.,rep=-1,ang=45,conv=1,to_fourier=None):
    '''ellipsometry angles from dielect. function
    '''
    from ellipse import calc_fourier
    from numpy import arctan2
    frac=reflect(epsil,n0,rep=0,ang=ang,polar='p')/reflect(epsil,n0,rep=0,ang=ang,polar='s')
    if rep==-2: return frac
    if conv: conv=180./pi
    else: conv=1
    if to_fourier!=None:
        afrac=abs(frac)
        cfrac=frac.real/afrac #cos delta
        rep=[]
        for f in to_fourier:
            rep.append(calc_fourier(afrac,cfrac,f/conv))
        return rep
    out=arctan2(abs(frac),1.)*conv,abs(arctan2(frac.imag,frac.real))*conv
    if rep==1: return out[0]+1j*out[1]
    else: return out

def calc_ellips_plate(freq,epsil,wid=[],ang=60,conv=1,rep=1,corr=None):
    '''calculates ellipsometric angles for given layers'''
    from numpy import arctan,arctan2,abs
    if conv: conv=180./pi
    if len(wid)==0: #bulk material
        zz=reflect(epsil[0],ang=ang,polar='p',rep=0)
        zz/=reflect(epsil[0],ang=ang,polar='s',rep=0)
    else:
        if ang<0:
            fr=plate(freq,epsil,wid,ang=ang,rep=-1)
            zz=friter([f[:,0] for f in fr[0]],fr[1],fr[2])/friter([f[:,1] for f in fr[0]],fr[1],fr[2])
        else:
            zz=plate(freq,epsil,wid,ang=ang,rep=0,polar='p')
            zz/=plate(freq,epsil,wid,ang=ang,rep=0,polar='s')
    if rep==-2: return zz
    out=arctan(abs(zz))*conv,-arctan2(zz.imag,zz.real)*conv
    if corr=='pos':
        out[1][out[1]<-90]+=360.
    if rep==0: return out[0]+1j*out[1]
    else: return out

global plout
plout=None
def selpoint(x,meas,epsil=[10+4j,2.5,10+2j],vari=[1,0,2],wid=[230,60.],angs=[60,70,80],rep=0):
	global plout
	from numpy import zeros,dot,array
	vepsil=zeros((len(epsil),max(vari)))
	for i in range(len(vari)): 
		if vari[i]>0:vepsil[i,vari[i]-1]=1
	fuc=lambda p,k:array([calc_ellips_plate(x[k],epsil+dot(vepsil,p),wid,rep=0,ang=a) for a in angs])
	if rep==1: return fuc
	from scipy import optimize
	plout=[zeros(max(vari))]
	out=[]
	for k in range(len(meas)):	
	    plout.append(optimize.fmin(lambda p:sum(abs(fuc(p,k)-meas[k])),plout[-1]))
	    out.append(sum(abs(fuc(plout[-1],k)-meas[k])))
	return out

global ualc
ualc=None
def multipoint(dang,wids,weig=None,repeaf=None,xpts=None,ypts=None,poly_deg=2,alog=None,aform='UV_%i',rep=0,check=False,dowid=True):
    global ualc
    from numpy import array,abs,polyfit,polyval,dot,ones,arange,iterable
    if alog!=None: #evaluate function fitted to measured values
        if iterable(alog):
            ualc=alog
        else:
            valc=[polyfit(alog.base,alog[aform%i],poly_deg) for i in dang]
            if type(xpts)==list: xpts=alog.base[xpts]
            ualc=array([polyval(a,xpts) for a in valc])
            if rep==1: return xpts,ualc
            print('fitted for polynom deg %i'%poly_deg)
            if check:
                for j in range(len(dang)):
                    print('ang %i:chi2 %.4f'%(dang[j],sum((polyval(valc[j],alog.base)-alog[aform%dang[j]])**2)))
    if weig==None: weig=ones(len(dang))
    if ypts!=None: #p0=list(ypts.real.flat)+list(ypts.imag.flat)+list(wids)
        p0=list(array([ypts.real,ypts.imag]).swapaxes(0,2).swapaxes(0,1).flat)+wids
        neps=len(ypts)
    else: 
        p0=None
        neps=len(wids)+1
    if repeaf==None: repeaf=arange(neps)
    #print 'shape:',neps,len(xpts),2
    feps=lambda p:dot(array(p[:neps*len(xpts)*2]).reshape(neps,len(xpts),2),array([1,1j]))
    if dowid:
        from ellipse import calc_ellips_plate
        wfit=lambda p:sum([sum(abs(calc_ellips_plate(xpts,feps(p)[repeaf],p[-len(wids):],ang=dang[i],rep=0,corr='pos')-ualc[i]))*weig[i] for i in range(len(dang))])
        return wfit,p0
    else:
        wfit=lambda p:sum([sum(abs(calc_ellips_plate(xpts,feps(p)[repeaf],wids,ang=dang[i],rep=0,corr='pos')-ualc[i]))*weig[i] for i in range(len(dang))])
        return wfit,p0[:-len(wids)]
    #ival=(p0!=None) and wfit(p0) or None
    # fitting only the width: dfit=lambda p:a(p0[:-2]+[p[0],p[1]])

#testing
#walc=array([profit.calc_ellips_plate(base[::5]],eps[:,::5][repeaf],wids,ang=i,rep=1,corr='pos') for i in dang])
#fres=[optimize.fmin(*list(profit.multipoint(dang,wids,xpts=base[[(i)*5]],repeaf=repeaf,ypts=eps[:,[(i-1)*5]],dowid=False,alog=walc[:,[i]]))) for i in range(1,len(walc))]

   
def stepbypoint(dang,wids,xpts=None,ypts=None):
    neps=len(wids)
    feps=lambda p:dot(array(p[:neps*len(xpts)*2]).reshape(neps,len(xpts),2),array([1,1j]))
    for i in range(1,len(xpts)):
        oldpts=feps(ypts)
        newpts=array([oldpts,oldpts])+b[-len(wids):]
        a,b,c=multipoint(dang,wids,xpts=xpts[i-1:i+1])
        b2=b[2:4]+b[6:8]
        #a2=lambda p:a(b[:2]+[p[0],p[1]]+b[4:6]+[p[2],p[3]]+b[8:])
        a2=lambda p:a(b[:2]+list(p[:2])+b[4:6]+list(p[2:4])+b[8:])
        
 

def wid_spread(freq,epsil,width,wrange=[0.9,1.1],ndiv=10,wind=0,rep=0,ang=0):
    '''consider non-uniformity of layer thickness
    with some optimization of repeated calculations
    wind: index of layer to vary
    '''
    from numpy import linspace,add
    wmod=linspace(wrange[0],wrange[1],ndiv)
    if rep<0: #adaptation for fitting procedures
          r,sh,psi=freq,epsil,width
          rep=-rep
    else: r,sh,psi=plate(freq,epsil,width,rep=-1,ang=ang)
    shtop=sh[:wind]
    shbot=sh[wind+1:]
    res=reduce(add,[friter([t.copy() for t in r],shtop+[m*sh[wind]]+shbot,psi) for m in wmod])/ndiv
    if rep==1: return res
    return (res*res.conj()).real

#def wid_fit(rep=1,opti='lbf',fit=None,fix=None,args=()):
#    out=mizer(fit,array(ipars),None,args=args,approx_grad=True,bounds=blims)
#    return out

global einf,elims,gang
einf=None
gang=0
elims=[-40.,40.]
global dims
dims=[]

def scan(fit,pars,ndiv=20,multi=False,imin=0,loud=0,ret=0):
    '''sweeping parameter space
    parameters given as values or ranges (then are divided to ndiv bins)
    profit.scan(profit.fit,[[.8,1.2,1],[1,1.2]],multi=True)
        in this case the range is divided as logspace(logarithmic bins)
    imin: first index to vary
    ret=1: find position of a minimum
    '''
    global dims
    rep=[]
    if ret>=1: dims=[]
    from math import log10
    from numpy import argmin,linspace,logspace
    for i in range(imin,len(pars)):
        if type(pars[i])==list:
            if len(pars[i])>3: p=pars[i] 
            elif len(pars[i])==3: p=logspace(log10(pars[i][0]),log10(pars[i][1]),ndiv)
            else: p=linspace(pars[i][0],pars[i][1],ndiv)
            #print 'scanning param %i [%i points]'%(i,len(p))
            if multi:
                j=0
                for d in p:
                    if loud:print('now %i:'%i+str(pars[:i]+[d]+pars[i+1:]))
                    rep.extend(scan(fit,pars[:i]+[d]+pars[i+1:],imin=i+1,loud=loud,ndiv=ndiv,ret=-j))
                    j+=1
                if ret>=0: dims.append(p)
                break
            else:
                rep.extend([fit(pars[:i]+[d]+pars[i+1:]) for d in p])
                if loud:print(str(pars[:i]+[d]+pars[i+1:]))
                if ret>=0: dims.append(p)
    if len(rep)==0: 
        try:
            rep=[fit(pars)]
            if loud:print(pars)
        except:
            print('error '+str(pars))
    else:
        if ret==2: return dims
        elif ret==1:
            print('calculating position of minima')
            pos=[]
            ipos=argmin(rep)
            for d in dims:
                print('now pos %i - min at %i'%(ipos,ipos%len(d)))
                pos.append(d[ipos%len(d)])
                ipos/=len(d)
            return pos
    return rep

def dofit(x,y,pars,lims=[[1e-2,1e6],[1,4000.],[0.001,100.]],yerr=None,einf=1.,drude=None,rep=1,opti='lbf',fix=None,args=()):
    '''fitting reflectivity/transmissivity with N-resonator model
    pars is a 2-d array (ampl,freq,absorb)
    errors not yet implemented
    fitting method either TNC, light BFGS or Cobyla ('tnc/lbf/cob')
    '''
    global gang,fit
    extrapars={}
    if opti=='tnc':
        from scipy.optimize import tnc as optimod
        mizer=optimod.fmin_tnc
        extrapars={'messages':0}
    else:
        from scipy.optimize import lbfgsb as optimod
        mizer=optimod.fmin_l_bfgs_b
    from numpy import array
    if fit==None: #no fit functions entered
        dierep=-1
        if len(y.shape)==2 and y.shape[1]==2: 
            y=y[:,0]+1j*y[:,1]
            if yerr!=None: weight=1/yerr[:,0]+1j/yerr[:,1]
            print("fitting complex numbers: ellipsometry")
            #dierep=
        elif str(y.dtype)[:7]=='complex': 
            dierep=0 # calculating in complex plane
            if yerr!=None: weight=1/yerr.real+1j/yerr.imag
            print("fitting complex numbers")
        else:
            if yerr!=None: weight=1/yerr
        if drude!=None:
            def fit(spars):
                ospars=array(spars[3:]).reshape((len(spars)-3)//3,3) #oscilator parameters
                dif=y.copy()
                if dierep==0: dif-=dielect(x,spars[0],ospars,drude=spars[1:3])
                else: dif-=rdielect(x,spars[0],ospars,drude=spars[1:3],ang=gang)
                if yerr!=None: dif*=weight
                return sum(abs(dif)**2)#dif*conj(dif))
        else:
            def fit(spars):
                ospars=array(spars[1:]).reshape((len(spars)-1)//3,3)
                dif=y.copy()
                if dierep==0: dif-=dielect(x,spars[0],ospars)
                else: dif-=rdielect(x,spars[0],ospars,ang=gang)
                if yerr!=None: dif*=weight
                return sum(abs(dif)**2)#sum(dif*conj(dif))
        if rep==-3: return fit
        #def fit(spars,args):
        #    return sum((args[1]-dielect(args[0],spars[0],array(spars[1:]).reshape((len(spars)-1)//3,3)))**2)
    if type(pars)==list: pars=array(pars)
    if len(pars.shape)==2: pars=pars.reshape(pars.shape+(1,))
    if einf==None: ipars=[1.]
    else: ipars=type(einf)==list and einf[:1] or [einf]
    if drude!=None: ipars+=list(drude)
    ipars+=list(pars[:,:,0].flat)
    if rep==-2: return ipars
    if lims!=None:
        if len(lims)==len(ipars): blims=lims
        else:
            blims=[elims]
            if drude!=None: blims+=[[0.1,100],[0.1,100.]]
            if len(lims)==len(ipars)-1: blims+=lims
            elif len(lims)==3:
                print('setting limits')
                for i in range(len(pars)):
                    for j in range(3):
                        if lims[j][0]<0: blims.append([pars[i,j,0]+a for a in lims[j]])
                        else: blims.append(lims[j])
            #if type(lims[0]==list)
        if rep==-1: return blims
        blims=array(blims)
    else:
        blims=None
    if opti=='cob':
        from scipy.optimize import fmin_cobyla as mizer
        par_con=[lambda p:(p[i]-blims[i][0])*(blims[i][1]-p[i]) for i in range(len(blims))]
        out=mizer(fit,array(ipars),par_con,args=args)
    else:
        out=mizer(fit,array(ipars),None,args=args,approx_grad=True,bounds=blims,**extrapars)
    return out


def ptbypt(meas,dang,freq,nlay=1):
    '''gets epsilon and layer thickness from ellips. measurements at different angles
    point-by-point
    '''
    dfit=lambda p,w:sum([abs(calc_ellips_plate(freq,p,w,ang=dang[i],rep=0)-meas[i])**2 for i in range(len(dang))])
    wfit=lambda q:dfit([q[i]+1j*q[i+1] for i in range(nlay+1)],q[-nlay:])
    return wfit
    # NOT FINISHED

global fit,idix
fit=None
per_conv=606.8 #conversion factors
beat_conv=1150. #430.
per_conv=[498.9,8.44e-04]

floating_norm=-1 #base level and calibration is adjusted with every fit
non_uniform=None
non_uniform_layer=0
idix=None
mod_layer=0

def multifit(x=None,y=None,epsil=None,wid=None,ang=0,mix=None,opti='lbf',yerr=None,smooth=None,dierep=0,ifit=None,nfit_lays=0,lims=None):
    '''fitting thickness of multiple layers / mixing ratios in 
    dierep: fitting complete dielectric function
    uses tabulated values of dielec. function (given in epsil)
    additional parameters to fit can produce polynomial shift of diel. fun
    (layer adjusted is specified by global "mod_layer" parameter)
    ADDed: 
        smoothing option to reduce sensitivity to fast variations
        floating_norm global parameter (in the case of smoothing)
    global setting    
        non_uniform: accounts for spread of widths over illuminated area (large diaphragm)
    '''
    global idix,fit
    #if ang>0: print 'incidence at %.1f deg '%ang
    if opti=='tnc':
        from scipy.optimize import tnc as optimod
        mizer=optimod.fmin_tnc
    else:
        from scipy.optimize import lbfgsb as optimod
        mizer=optimod.fmin_l_bfgs_b
    from numpy import array,all,conj,sqrt,ones
    if wid[0]==0: #no assumptions about initial parameters
        from spectra import fitting
        out=fitting(x,y)
        wid[0]=per_conv[0]/((out[0][3]-per_conv[1])*sqrt(epsil[0].real.mean()))
        #wid[0]=per_conv/(out[3]*sqrt(epsil[0].real.mean()))
        print('estimated principal layer width %f'%wid[0])
    if yerr!=None and  all(yerr.real>0): 
        weight=1/yerr.real
        if all(yerr.imag>0): weight+=1j/yerr.imag
    if ifit==None: #no fit functions entered
        r,sh,psi=plate(x,epsil,wid,ang=ang,rep=-1)
        from numpy import convolve,polyfit,polyval
        if smooth!=None:
            y=convolve(y,smooth)[len(smooth):-len(smooth)]
            print('data shape %i'%y.shape)
            def fit(spars,rep=0):
                global idix
                #print spars
                slate=convolve(abs(friter([t.copy() for t in r],[sh[i]*spars[i] for i in range(len(sh))],psi))**2,smooth)[len(smooth):-len(smooth)]
                #slate=convolve(plate(x,epsil,list(spars),meth=0,ang=ang,rep=dierep),smooth)[len(smooth):-len(smooth)]
                if floating_norm>0:
                    idix=polyfit(slate,y,floating_norm)
                    dif=y-polyval(idix,slate)
                else: dif=y-slate
                if yerr!=None: dif*=weight
                if rep==1: return dif
                return sum((dif*conj(dif)).real)
        else:
            def fit(spars,rep=0):
                global idix
                #print spars
                rlays=[t.copy() for t in r]
                shlays=[sh[i]*spars[i] for i in range(min(len(sh),len(spars)))]
                if len(spars)<len(sh):
                    shlays.extend(sh[len(spars):])
                elif len(spars)>len(sh):
                    prof=polyval(spars[:len(sh)-1:-1],x)
                    rlays[mod_layer]*=prof
                    shlays[mod_layer]*=prof
                if non_uniform:
                    dif=y-wid_spread(rlays,shlays,psi,[1-abs(non_uniform),1+abs(non_uniform)],wind=non_uniform_layer,ang=ang,rep=-2)
                if floating_norm>=0:
                    slate=abs(friter(rlays,shlays,psi))**2
                    if floating_norm==0:
                        from extra import rob_polyfit
                        a,c=rob_polyfit(slate,y,wei=-1)
                        if c>0.6: # minimal correlation to do renormalization
                            idix=rob_polyfit(slate,y,wei=2)
                            dif=y-polyval(idix,slate)
                        else: dif=y-slate
                    else:
                        idix=polyfit(slate,y,floating_norm)
                        dif=y-polyval(idix,slate)
                else: 
                    dif=y-abs(friter(rlays,shlays,psi))**2
                #dif=y-plate(x,epsil,list(spars),meth=0,ang=ang,rep=dierep)
                if yerr!=None: dif*=weight
                if rep==1: return dif
                return sum((dif*conj(dif)).real)
        if dierep==-2: return fit
    else: fit=ifit
    if nfit_lays==0:nfit_lays=len(wid)
    if lims: out=mizer(fit,ones(nfit_lays),None,args=(),approx_grad=True,bounds=array(lims))
    else: out=mizer(fit,ones(nfit_lays),None,args=(),approx_grad=True)
    return out,wid

def get_errors(mlog,i,nwid=None,rep=1): ###unfinished###
    '''estimation of errors of multilayer fit
    '''
    if nwid: fit=multifit(mlog.base[imin:imax],(mlog[slist[i]]*norm)[imin:imax],epsil[[0,1,0]][:,imin:imax],nwid,dierep=-2)
    rep=scan(fit,[[.99,1.01],[0.95,1.05]],ndiv=20,multi=True)
    contour(array(rep).reshape(20,20))
    from extra import chi2map_anal
#=========================================================================================

def unittest(mode=1,base=1,comp=["SiO2_gl","cSi_asp"],frange=[0.3,3.3],wid=3000):
    '''single layer-substrate reflection
        comp: component list (if only one given, use Si standard substrate)

        loading of the database:profit.unittest(base=e2,comp=["SiO2_gl","cSi_asp"],frange=None)
    '''
    from numpy import arange,loadtxt,concatenate
    from spectra import dbload
    if len(comp)<=1:
        fa2,sir,sim=loadtxt('/home/limu/Lab/si_dielfun.dat',unpack=True)
        va2=sir+1j*sim
    else:
        fa2,va2=dbload(comp[1],connection="http")
        sir,sim=va2.real,va2.imag
    if (type(base)==int) and (base==0): freq=fa2
    else:
        if type(base) in [int,float]:
            step=0.005/base
            freq=concatenate([arange(0.01,0.9,step),arange(0.9,6,step*4)])
        else:
            freq=base
        from spectra import respect
        va2=respect(freq,[fa2,sir,sim])
    fa1,va1=dbload(comp[0],freq)
    if frange!=None:
        sel=freq>frange[0]
        sel*=freq<frange[1]
    else: sel=freq>0
    if mode==1: res=plate
    elif mode==2: res=matter_plate
    return freq[sel],[va1[sel],va2[sel]],res(freq[sel],[va1[sel],va2[sel]],[wid])


def plotwo(e,idata,fig=None,mode='pd',clean=True,ang=75):
    '''you can plot according to 'mode':
    pd: psi/delta
    pd-tc: tan(psi),cos(delta)
    eps: dielectric function
    nk: refractive indices
    '''
    from numpy import sqrt,iterable
    if fig==None:
        from matplotlib.figure import Figure
        fig=Figure()
    elif clean: fig.clf()
    from matplotlib.pyplot import subplot
    if iterable(idata[0]): data=idata[0]+1j*idata[1]
    else: data=idata
    if mode[:2]!='pd':data=from_ellips(data.real,data.imag,ang=ang,unit='deg')
    if mode=='nk':data=sqrt(data)
    if mode=='pd_ct':
        from numpy import tan,cos
        data=tan(data.real*pi/180)+1j*cos(data.imag*pi/180)
    for j in range(2):
        axa=subplot(2,1,j+1)
        if clean:
            axa.set_xlabel('energy [eV]')
            if mode=='nk':axa.set_ylabel(['n','k'][j])
            elif mode=='eps':axa.set_ylabel('$\epsilon$ '+['real','imag'][j])
            elif mode=='pd_ct': axa.set_ylabel(['tan $\Psi$','cos $\Delta$'][j])
            elif mode=='pd': axa.set_ylabel(['$\Psi$','$\Delta$'][j])
        if j==0:axa.plot(e,data.real)
        else:axa.plot(e,data.imag)

def load_huml(fname='Lab/Prakt/sio2k1.out'):
    '''
    usage: mee=profit.load_huml('Lab/spectra/prakt/sio2k2.out')
    data=loadtxt('Lab/spectra/prakt/sio2k2.dpt',unpack=True)
    sele=data[0]<1700
    elims=mee[2][0]
    rep=profit.dofit(data[0][sele],data[1][sele],mee[1],mee[2][1:],opti='tnc')[0]
    '''
    from numpy import array,size
    res=open(fname).readlines()
    p=[i+1 for i in range(len(res)) if res[i].find('C.')>0]
    pars=array([float(a.split()[2]) for a in res[p[1]+1:-6]])
    pars=pars.reshape(len(pars)//3,3)
    ipars=array([float(a.split()[-3]) for a in res[p[0]+1:p[0]+1+size(pars)]]).reshape(pars.shape)
    einf=float(res[p[0]].split()[-3])
    lims=[map(float,a.split()[-2:]) for a in res[p[0]:p[0]+1+size(pars)]]
    return pars,ipars,lims

global mfit
mfit=None
def splfit(base,rep,pix,ivals,drang=0.05):
    '''looks for optimal places of spline points to interpolate rep=func(base): 
        doing shifts from initial (pix) positions by at most _drang_
        could also take into account some weights
    '''
    from scipy import interpolate,optimize
    from numpy import zeros,ones,array
    dpix=zeros(pix.shape)
    def xfit(cpix):
        global mfit
        mfit=optimize.leastsq(lambda y0:interpolate.splev(base,interpolate.splrep(pix+cpix,y0))-rep.real,ivals)
        return sum((interpolate.splev(base,interpolate.splrep(pix+cpix,mfit[0]))-rep.real)**2)
    return optimize.fmin_l_bfgs_b(xfit,dpix,approx_grad=True,bounds=array([-drang*ones(pix.shape),drang*ones(pix.shape)]).transpose())
