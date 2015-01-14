# -*- coding: utf-8 -*-
''' manipulation of ellipsometric data
conversion to/from ellipsometric measurements
    from_fourier
    from_ellips/calc_ellips
'''

from math import pi

def from_ellips(ipsi,idelt=None,n0=1.,rep=-1,ang=45,unit='deg',errs=None):
    '''calc. dielect, function from ellipsometric angles
    errs=True: with errors
    '''
    from numpy import sin,tan,cos,exp
    if idelt==None:
        idelt=ipsi.imag
        ipsi=ipsi.real
    if unit=='deg':
        ang*=pi/180
        psi=ipsi*pi/180
        delt=idelt*pi/180
        if errs!=None: 
            dpsi=errs.real*pi/180
            ddelt=errs.imag*pi/180
    else:
        psi,delt=ipsi,idelt
        if errs!=None:
            dpsi,ddelt=errs.real,errs.imag
    pom=1+sin(2*psi)*cos(delt)
    pom2=cos(2*psi)**2-(sin(2*psi)*sin(delt))**2
    #if ang<0: return pom,pom2
    epsr=(n0*sin(ang))**2*(1+(tan(ang)/pom)**2*pom2)
    epsi=(n0*sin(ang)*tan(ang)/pom)**2*sin(4*psi)*sin(delt)
    if errs==None: return epsr+1j*epsi
    rho=tan(psi)*exp(1j*delt)
    deps2drho=4*sin(ang)**2*tan(ang)**2*(1+rho)/(1-rho)**3
    deps2dpsi=-exp(1j*delt)*deps2drho/cos(psi)**2
    deps2delt=-1j*exp(1j*delt)*deps2drho*tan(psi)
    return epsr+1j*epsi,deps2dpsi*dpsi+deps2delt*ddelt
    #deps2dang=sin(2*ang)*(1+(1+1/cos(ang)**2)*tan(ang)**2*(1+rho)**2/(1-rho)**2)


def from_fourier(al,be=None,corr=0,norm=1.,conv='no',nois=None,calerr=None,loud=0,rep=0):
    ''' calculates psi/delta (in radians) from fourier coefs with angle corr (in radians)
    corr - analyser correction (in radians)
    norm - tangent of ds(polarizer+correction)
    conv - some conventions for delta
        bes: 180-delta
        tc: returns tangent psi/ cosinus delta
        deg: using degrees instead of radians
    nois: array of errors
    calerr: errors of calibration
    '''
    from numpy import cos,sin,tan,sqrt,arctan,arccos,abs,sign,iterable
    if conv.find('inv')>=0:al,be=be,al
    if be==None: #complex numbers
        if hasattr(corr,'imag') and corr.imag!=0: 
            al=al*corr
            corr=0
        be=al.imag
        al=al.real
    if conv.find('deg')>=0:
        corr*=pi/180.
        norm=tan(norm*pi/180.)
    if conv.find('bes')>=0: corr=2*(pi-corr)
    ial,ibe=al.copy(),be.copy()
    if nois!=None:
        if iterable(nois):
            dal=nois.real
            dbe=nois.imag
        else:
            dal,dbe=1/(1-al)+(1+al)/(1-al)**2,   2*sqrt((1+al)/(1-al))
            dtpsi=dal*nois/dbe
            dcdelta=1/sqrt(1-al**2)*nois+(al*be)/(1-al**2)**(3/2)*nois;
    if corr!=0: 
        #dtpsi*=corr
        al,be=(al*cos(corr)-be*sin(corr)),  (al*sin(corr)+be*cos(corr))
    if loud>0: print('rotated to %.5f, %.5f'%(al,be))
    if norm==0: return al,be
    if nois==None and rep==3: nois=True

    al[al>1]=0.999
    al[al<-1]=-0.999
    cdel=be/sqrt(1-al**2)
    cdel[cdel>1]=0.999
    cdel[cdel<-1]=-0.999
    if conv.find('tc')>=0: 
        if conv.find('pos')>=0: cdel=abs(cdel)
        return sqrt((1+al)/(1-al))*abs(norm),cdel
    delta=sign(norm)*arccos(cdel)
    if conv.find('pos')>=0 and delta.mean()<0: delta=pi+delta
    #elif conv=='bes': delta=pi-delta
    psi=arctan(sqrt((1+al)/(1-al))*abs(norm))
    if nois!=None:
        dpsi=abs(norm)/(1+tan(psi)**2)*sqrt((1-ial)/(1+ial))*dal/(1-ial)**2
        ddelt=1/sin(delta)*sqrt(dal**2*ial**2*ibe**2/(1-ial**2)**3 + dbe**2/(1-ial**2))
    if conv.find('deg')>=0:
        psi*=180./pi
        delta*=180./pi
        if nois!=None:
            dpsi*=180./pi
            ddelt*=180./pi
    if rep==3: return psi+1j*delta,dpsi+1j*ddelt
    if rep==1: return psi+1j*delta
    return psi,delta

#simple calibration test
#dall=array([profit.from_fourier(ruv[w],None,corr=a0/180.*pi,norm=tan((int(w[16:-4])-p0)/180.*pi)) for w in w1])

global calfun,alist,blist,elist,plist
elist=None
def calib_ellips(flist=None,errs=True,anal0=0,polar0=0,sep=','):
    '''trying to find calibration parameters'''
    global calfun,alist,blist,plist,elist
    from numpy import loadtxt
    from math import tan
    import os
    if flist!=None: alist,blist,plist=[],[],[]
    else: flist=[]
    for f in flist:
        try:
            e,a,b,ea,eb=loadtxt(f,skiprows=2,delimiter=sep,unpack=True)
        except:
            continue
        pars=os.path.splitext(os.path.basename(f))[0].split('_') 
        #naming convention is "sample","band","angle","polar"
        alist.append(a)
        blist.append(b)
        plist.append(float(pars[-1]))
    if len(alist)>0: print('loaded %i measurements of %i bins'%(len(alist),len(alist[0])))
    def calfun(angs,both=False,ran=[None,None]):
        from numpy import std
        psi,delt=[],[]
        norm=1.
        if len(angs)>2:norm=angs[2]
        for i in range(len(plist)):
            py,dy=from_fourier(alist[i][ran[0]:ran[1]]*norm,blist[i][ran[0]:ran[1]]*norm,(angs[0])/180.*pi,tan((plist[i]+angs[1])/180.*pi))
            psi.append(py)
            delt.append(dy)
        if elist:
            return sum(std([psi[i]/elist[i] for i in range(len(psi))],0))/norm
        if both: return (sum(std(psi,0)**2)+sum(std(delt,0)**2))/norm
        else: return sum(std(psi,0))/norm
    return calfun

# usage: ala=glob("Lab/neboj/elipso/Hopg/hop*eab")
# kop3=calib_ellips([zub[i] for i in [2,5]])
# optimize.fmin(kop3,[70,0])

#condictivity (eps.imag-i*(eps.real-1))*freq*eps_0
def calc_nk(n2_k2,nk2):
    from numpy import sqrt
    n=sqrt((sqrt(n2_k2**2+nk2**2)+n2_k2)/2.)
    k=nk2/n/2.
    return n,k

def calc_epsi(n,k):
    return (n+k*1j)**2

def calc_fourier(afrac,cfrac,ang):
    from numpy import tan
    norm=tan(ang)
    al=afrac**2-norm**2
    be=2*afrac*norm*cfrac
    return (al+1j*be)/(afrac**2+norm**2)

def calc_ellips_plate(freq,epsil,wid=[],ang=60,conv=1,rep=0,corr=None):
    '''calculates ellipsometric angles for given layers'''
    from numpy import arctan,arctan2,abs
    from profit import plate,reflect
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
    if rep==1: return out[0]+1j*out[1]
    else: return out

def ptbypt(meas,dang,freq,nlay=1):
    '''gets epsilon and layer thickness from ellips. measurements at different angles
    point-by-point
    '''
    dfit=lambda p,w:sum([abs(calc_ellips_plate(freq,p,w,ang=dang[i],rep=1)-meas[i])**2 for i in range(len(dang))])
    wfit=lambda q:dfit([q[i]+1j*q[i+1] for i in range(nlay+1)],q[-nlay:])
    return wfit
    # another example
    # moje=lambda p0:dot(weig,abs(array([profit.calc_ellips_plate(xpts[[1]],dot(r_[1,1j],array(p0[:4]).reshape(2,2))[[1,0,1]].reshape(3,1),p0[4:],ang=i,rep=1) for i in dang])[:,0]-walc[:,1]))


#-----------------high-level process

def neboj_load(name,dlot=0):
    from pylab import load
    dut=load(name,skiprows=2)
    e1=dut[0][2::4]
    e2=dut[1][2::4]#.reshape(832,4)[:,0]
    f=dut[:,:2].transpose()
    g=dut[:,2:].reshape(dut.shape[0],832,4)
    if dlot:
        [plot(e1,q[:,0]) for q in e[::dlot]]
        legend(f[1][::dlot])
    return g,f,e1,e2
    
def neboj_group(f,g,side=None,nbin=12,min_cnt=1):
    from numpy import int as Int
    from math import sqrt
    if side=='rise': imin,imax=None,f[1][65:].argmax()
    elif side=='fall': imin,imax=f[1][65:].argmax(),None
    else: imin,imax=None,None
    ids=f[1][imin:imax].copy()
    ids-=ids.min()
    ids*=nbin/ids.max()
    ids=ids.astype(Int)
    gids=[i for i in range(max(ids)+1) if sum(ids==i)>=min_cnt]
    cnts=[sum(ids==i) for i in gids]
    rep=[g[imin:imax][ids==i].mean(0) for i in gids]
    for k in range(len(rep)):
        for i in [2,3]:
            rep[k][:,i]/=sqrt(cnts[i])
    tep=[f[1][imin:imax][ids==i].mean(0) for i in gids]
    return tep,rep,cnts
    
def neboj_fit(l):
    l2=l[1][60:,0]+1j*l[1][60:,1]
    import profit
    #fun=profit.dielect(e2[60:],0.5,[[2.,2.5,.5]],rep=0)
    profit.einf=1
    drude=[3.66,0.02]
    out=profit.dofit(e2[60:],abs(l[1][60:,0]),array([[-2.,2.5,.5]]),drude=drude)
