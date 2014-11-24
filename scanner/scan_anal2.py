#!/usr/bin/env python
# -*- coding: utf-8 -*-

import spectra,profit
from numpy import load,array,sqrt,concatenate,argsort,argmin,polyfit,polyval,sign,sum,median,r_
import extra

global px,cl,dt,sel
px,cl,dt=None,None,None
sel=None

conv_per=281.0#-5.28*(wid2-1) oxide layer in microns
#conv_fun=620./per-73-0.15*(wid2-1)
conv_per=609.7#-10.70*(wid2-1)

# calibration error 0.13%
half_anal_width=5

def variab(dt,blen=50):
    '''is signal periodic?'''
    bnum=len(dt)//blen
    dmid=dt[:bnum*blen].reshape(bnum,blen).mean(1)
    dmax=dt[:bnum*blen].reshape(bnum,blen).max(1)
    dmin=dt[:bnum*blen].reshape(bnum,blen).min(1)
    return (dmax-dmin)/dmid

commchar="#"
def init(dtfile="soinv022.npy",emin=1.42,emax=2.4,ddir='data/oct10/',smth=11,pxfile="pixnew.npy",clfile="calinv.npy",dord=-1):
    
    from numpy import loadtxt
    global px,cl,dt,sel
    exrgs={'unpack':True,'comments':commchar}
    if ddir:
        if clfile: cl=loadtxt(ddir+clfile,**exrgs) if (clfile.find('npy')<0) else load(ddir+clfile)
        if cl!=None and len(cl.shape)==2:
            px=cl[0]
            cl=cl[1]
        else:
            if pxfile: px=loadtxt(ddir+pxfile,**exrgs) if (pxfile.find('npy')<0) else load(ddir+pxfile)
        if dtfile: dt=loadtxt(ddir+dtfile,**exrgs) if (dtfile.find('npy')<0) else load(ddir+dtfile)

    if smth>0:
        from numpy import hamming,convolve
        mom=hamming(smth)
        mom/=mom.sum()
        for i in range(len(dt)):
            dt[i]=convolve(dt[i],mom,'same')
        if clfile: cl=convolve(cl,mom,'same')
    if sel==None or pxfile!=None: sel=(px>emin)*(px<emax)
    if clfile: cl=cl[sel][::dord]
    dt=dt[:,sel][:,::dord]
    if pxfile: px=px[sel][::dord]
    if cl!=None: dt/=cl

def init_single(fname,emin=1.42,emax=2.4):
    global px,cl,dt,sel
    from numpy import loadtxt
    px,dx=loadtxt(fname,unpack=True,comments="%")
    #imin,imax=spectra.get_ids(px,[emin,emax])
    sel=(px>emin)*(px<emax)
    dt=[dx]

global mids,pers,dn1
mids,pers=None,None
dn1=None

global var_limit,sig_limit
var_limit=0.1
dat_rang=[0.9,1.9]
sig_limit=2
min_max=0

def ind_deriv(mids,pars=None,der_step=0.01,absorb=0.1):
    if pars==None: pars=spectra.lor_si_low
    if absorb==None:
        rep=(mids+der_step)*sqrt(profit.dielect(mids+der_step,pars[0],pars[1]))
        rep-=(mids-der_step)*sqrt(profit.dielect(mids-der_step,pars[0],pars[1]))
    else:
        rep=(mids+der_step)*sqrt(profit.dielect(mids+der_step,pars[0],pars[1]).real+1j*absorb)
        rep-=(mids-der_step)*sqrt(profit.dielect(mids-der_step,pars[0],pars[1]).real+1j*absorb)
    return rep/2/der_step
    
def anal(inp=None,j=0,rep=1,cmin=3,crang=None,vlim=0.3,amplim=None,der_step=0.01,absorb=0.1,erng=None,use_min=0):
    '''j=0/1 - hodnoty z min/max
    absorb: neni-li None, bereme realnou cast dielekt. fce + danou (konstatni) absorbci
    '''
    global mids,pers,dn1
    fail=[0]*(rep+1)
    from numpy import median
    ipx=px
    if type(inp)==int:
        i=inp
        if erng!=None:
            inp=dt[i][(px>erng[0])*(px<erng[1])]
            ipx=px[(px>erng[0])*(px<erng[1])]
        else:
            inp=dt[i]
    if crang!=None:
        ww=inp.mean()
        if (ww<crang[0]) or (ww>crang[1]): 
            fail[1]=ww
            return fail
    if vlim:
        k=variab(inp)
        if median(k)<vlim:
            print 'too low variability '
            fail[1]=median(k)
            return fail
    apos=array(spectra.extrema(inp,ipx,msplit='sort',extsel=half_anal_width,ampsel=amplim,loud=0)[j])
    if sum(apos==0)>0:
        print "%i points removed"%sum(apos==0)
        apos=apos[apos>0]
    if der_step>0: #correction by refr. index
        mids=(apos[1:]+apos[:-1])/2.
        dn1=ind_deriv(mids,der_step=der_step,absorb=absorb)
        pers=(apos[1:]-apos[:-1])*dn1.real
    else: pers=(apos[1:]-apos[:-1])
    if rep==0: return pers
    if len(pers)<cmin:return fail
    if sig_limit>0:
        psel=extra.clipped(pers,sig_limit,2) 
        pers=pers[psel]
        mids=mids[psel]
    if use_min!=0:
        from numpy import sort
        if use_min>0: out=[sort(pers)[:use_min].mean()] #minimal values
        else: out=[sort(pers)[use_min:].mean()] #maximal values
    else: out=[median(pers)]
    out+=[pers.std()]
    if rep>1:
        out.append(len(pers))
        if rep>2: out.append(polyfit(mids,pers,1)[0])
        if rep>3: out.append((pers-out[0]-out[3]*(mids-mids.mean())).std())
    return out

import sys,os
def meas_all(mask='2010_dec_1.3_0%i.npy',ddir='/home/munz/Lab/spectra/scan/dec10/',iset=30,erng=[1.2,1.6],
    pxfile='pixels.npy',clfile='calib.npy',rep=1,zig=0):
    '''make complete analysis of a 2-D dataset
    zig: odd lines are measured in opposite direction than even ones
    '''
    global sel
    sel=None
    if type(iset)==int: iset=range(iset)
    if zig>0: dird=2*(iset[0]%2)-1
    else: dird=-1
    
    init(mask%iset[0],erng[0],erng[-1],ddir=ddir,pxfile=pxfile,clfile=clfile)
    #lst=[os.path.basename(a) for a in glob('Lab/spectra/scan/dec10/20*')]
    if len(erng)>2:
        out=[[[anal(i,min_max,rep=rep,vlim=var_limit,crang=dat_rang,erng=erng[k:k+2]) for k in range(len(erng)-1)] for i in range(len(dt))][::dird]]
    else:
        out=[[anal(i,min_max,rep=rep,vlim=var_limit,crang=dat_rang) for i in range(len(dt))][::dird]]
    for j in iset[1:]:
        if zig>0: dird=2*(j%2)-1
        if not os.path.exists(ddir+mask%j): continue
        init(mask%j,erng[0],erng[-1],ddir=ddir,pxfile=None,clfile=None)
        if len(erng)>2:
            out.append([[anal(i,min_max,rep=rep,vlim=var_limit,crang=dat_rang,erng=erng[k:k+2]) for k in range(len(erng)-1)] for i in range(len(dt))][::dird])
        else:
            out.append([anal(i,min_max,rep=rep,vlim=var_limit,crang=dat_rang) for i in range(len(dt))][::dird])
        if j%5==0: 
            sys.stdout.write("finished line %i\n"%j)
            sys.stdout.flush()
    try:
        return array(out).transpose()
    except:
        return out


if __name__ == '__main__':
    import sys,os
    fname=sys.argv[1]
    exargs={}
    for a in sys.argv[2:]:
        if a[:5]=='erng=': exargs['erng']=float(a[5:])
    if fname[-4:]=='.npy':
        fdir=os.dirname(fname)
        fname=os.basename(fname)
        rep=meas_all(fname,fdir)
    else:
        init_single(fname)
        rep=anal(range(len(dt)))
    for a in rep:
        print a


global counter
counter=0
def bump(w2=310,w1=1500,vl=None,slim=[2.2,2.6],rep=1,min_dis=0.015,loud=0,fit_wid=2):
    global counter
    if vl==None: 
        oxide=polyval(spectra.cau_sio2,px)
        epsil=profit.dielect(px,spectra.lor_si_low[0],spectra.lor_si_low[1])
        vl=profit.plate(px,[epsil,oxide,epsil],[w1,w2])
    a,b=spectra.extrema(vl)
    if len(a)<len(b):a,b=b,a
    if a[0]>b[0]:
        sh=sum(b<a[0])
        b=b[sh:]
    elif a[1]<b[0]:
        sh=sum(a<b[0])
        a=a[sh-1:]
    if a[-1]<b[-1]:
        sh=sum(b>a[-1])
        b=b[:-sh]
    elif a[-2]>b[-1]:
        sh=sum(a>b[-1])
        a=a[:-sh+1]
    assert len(b)==len(a)-1
    m2=concatenate([(px[b]+px[a][:-1])/2.,(px[b]+px[a][1:])/2.])
    n2=concatenate([(vl[b]-vl[a][:-1])/2.,(vl[b]-vl[a][1:])/2.])
    s2=argsort(m2)
    n2=n2[s2]
    m2=m2[s2]
    sel=n2>min_dis
    m2=m2[sel]
    n2=n2[sel]
    if rep==0: return m2,n2
    sel=(m2>slim[0])*(m2<slim[1])
    if sum(sel)<2*fit_wid:
        return 0,0
    i1=argmin(n2[sel])
    s=(m2[sel][i1]+m2[sel][i1+1])/2.
    j1=i1>fit_wid and i1-fit_wid  or None
    j2=i1+fit_wid+1
    counter+=1
    k,l=polyfit(m2[sel][j1:j2],(n2[sel]*sign(m2[sel]-s))[j1:j2],1)
    chi2=sum((polyval([k,l],m2[sel][j1:j2])-(n2[sel]*sign(m2[sel]-s))[j1:j2])**2)
    return -l/k,chi2

################# procedures to test the thickness reconstruction ####################

'''
polynomial fit
spectra.sin_fin=None;spectra.fitting(scan_anal2.px[sel],scan_anal2.dt[23][sel],refr_mode='poly')
1.35-1.65 dec_13 15/23
(array([  3.97798673e-06,   2.87115694e+00,  -4.66096643e+00,
         1.44617886e-02,   3.05471837e+00,   3.37296768e+00,
        -1.10670957e+01,   1.03433238e+01,  -1.66418610e-02,
         4.64434510e-03]),
 1.0301616864616243)
'''

def test_prepare(wrange,w2=1000,norm=620,j=0):
    global sig_limit,dt
    import profit
    tsig=sig_limit
    sig_limit=0
    epsil=profit.dielect(px,spectra.lor_si_low[0],spectra.lor_si_low[1])
    oxide=polyval(spectra.cau_sio2,px)
    if (dt==None) or (len(wrange)!=len(dt)): dt=[profit.plate(px,[epsil,oxide,epsil],[w1,w2]) for w1 in wrange]
    out=[]
    for a in range(len(wrange)): 
        if anal(a,j,rep=2)[2]>2: out.append([mids,(norm/pers)])
    sig_limit=tsig
    return out

def test_mids(inp,erange=[1.6,1.9]):
    out=[]
    for a in inp:
        asel=(a[0]>erange[0])*(a[0]<erange[1])
        amin=a[1][asel].argmin()
        sit=polyfit(a[0][asel],a[1][asel],2)
        out.append([sum(asel),amin,a[1][asel][amin],sit[2]-(sit[1]**2/4./sit[0])])
    return array(out)

def reconstruct(x,y8,points=None,di=10,re_ord=2,n_subs=3.46,thick=2000.,im_ini=0,im_ind=-1):
    if points==None: points=range(5*di,len(x)-5*di,5*di)
    res=spectra.extrema(y8,x,poly=-1,msplit='sort',extsel=half_anal_width,loud=0)
    mis=[((array(o[1:])+array(o[:-1]))/2.,(array(o[1:])-array(o[:-1]))) for o in res]
    from extra import polyfit as rob_fit
    n_pr=array([rob_fit(m[0],1/(spectra.ev2um*m[1]*2.*thick/1000.),re_ord,clip=1,lim=1.1) for m in mis])
    n_re=polyval(n_pr.mean()/r_[re_ord+1:0.:-1.],x) #real part from oscillations
    moo=[im_ini]
    from scipy import optimize
    for i in points:
        muf=lambda z:sum((y8[i-di:i+di]-profit.plate(x[i-di:i+di],[(n_re[i-di:i+di]+1j*z)**2,n_subs**2],[thick]))**2)
        moo.append(optimize.fmin(muf,moo[im_ind])[0])
    return n_pr,moo
