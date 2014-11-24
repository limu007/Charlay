# -*- coding: utf-8 -*-
global bins,refe,mall,dall,zzall,nomal,epsil
bins,refe=None,None
zzall=None
mall,dall=None,None
nomal,epsil=None,None
import spectra
from numpy import array,hamming,loadtxt,abs,convolve,where,mean,sum,polyval,polyfit,r_
from profit import unittest,reflect,plate,dielect

def init(iban=[500,900],nfil=52,redir='/home/limu/Lab/spectra/scan/',norm='new'):
    global bins,refe,mall,dall,zzall,nomal,epsil
    import os
    if norm=='new':
        bins2,refe2=loadtxt('/home/limu/Lab/spectra/soi3/radek3/novy_bod.dat',unpack=True,comments='%')
        ## plot(bins2,scan_anal.zzall[3][3][2580:0:-1],hold=0)
        uni2=unittest(base=bins2,comp=['cSi_asp','gSiO2'])
        bins=uni2[0]
        epsil=array(uni2[1]);refe=reflect(epsil[0])
    else:
        bins,refe=loadtxt('/home/limu/Lab/TblRefSi3648.dat',unpack=True)#[:,0]
        bins=1/spectra.ev2um/bins[:2600]*1000
    from numpy import load
    zzall=[load(redir+'fline0%02i.npy'%i) for i in range(nfil)]
    # how to detect oscilations
    if os.path.exists(redir+'nomal.npy'):
        nomal=load(redir+'nomal.npy')[1][:len(bins)]
    else:
        nr=hamming(11);nr/=nr.sum()
        mall=[[abs(z[500:900:10][5:-5]-convolve(z[500:900:10],nr,'same')[5:-5]).mean() for z in zz] for zz in zzall]
        dall=[[z[500:900].mean() for z in zz] for zz in zzall]

def calib(extra=False):
    global bins,refe,mall,dall,nomal
    #calibration plate
    sele=(array(mall)<20)*(array(dall)>700)
    from scipy import ndimage
    sele2=ndimage.grey_erosion(sele,2)
    ia,ib=where(sele2>0.5)
    nomal=mean([zzall[ia[i]][ib[i]] for i in range(100)],0)
    if extra:
        #testing individual lines:
        ps=array([sum(ia<i) for i in range(30,41)])
        koral=[mean([zzall[i+30][b] for b in ib[ps[i]:ps[i+1]]],0) for i in range(10)]
        return nomal,koral
    return nomal

import sys
from numpy import zeros,polyval,sqrt
######## zz/nomal*refe[-2600:]
def fitting(rang=[0,10],res=None,erange=[1.5,1.9],plim=120,prange=[0.025,0.035]):
    '''
    prange: range of allowed periods 
    '''
    if res==None:
        res=array([ -3.6493, 12.2222, -10.0929, 0.02989, -0.5501, -0.8253, 3.0086, -2.462])

    sint=(bins>erange[0])*(bins<erange[1])
    rar=zeros((len(zzall[0]),8))
    rep=[]
    for k in range(rang[0],rang[1]):
        rar[:,:]=0
        selel=array(mall[k])>plim
        sys.stdout.write('scanning %i points in line %i\n'%(sum(selel),k))
        sys.stdout.flush()
        zz2=zzall[k]/nomal*refe[:2600]
        parad=[spectra.fitting(bins[sint],zz2[i][sint],res,prange=prange) for i in where(selel)[0]]
        rar[(array(mall[k])>120)]=array([p[0] for p in parad])
        mval=[mean(polyval(a[5:],bins[sint])) for a in rar]
        rep+=[rar.copy()]
    return rep
    
def mean_amp(resall,ib=[430,890]):
    ren=[1/3.,1/2.,1.]
    kok=array([[abs(polyval(p[:3]*ren,bins[ib[1]])*bins[ib[1]]) for p in r] for r in resall])
    kok-=array([[abs(polyval(p[:3]*ren,bins[ib[0]])*bins[ib[0]]) for p in r] for r in resall])
    return kok
    #calc. positions of minima
    disc=resall[:,:,1]**2-4*resall[:,:,0]*resall[:,:,2]
    from numpy import ma
    pos1=ma.masked_array((resall[:,:,1]-sqrt(disc))/2/resall[:,:,0])
    pos2=ma.masked_array((resall[:,:,1]+sqrt(disc))/2/resall[:,:,0])

#zz2=[scan_anal.zzall[3][5][2234:0:-1][:2000]/scan_anal.nomal[:2000]*scan_anal.refe[:2000] for z in scan_anal.zzall[3]]
#ia,ib=spectra.get_ids(scan_anal.bins,[1.8,2.6])

def rest():
	from numpy import load
	zz=[load(a) for a in alst]
	cl=load("calib.npy")
	rep=profit.unittest(base=px[:1600],comp=["SiO2_gl","cSi_asp"])

def extra_test(wid1=1300,wid2=300,step=0.002,isel=1,levabs=0.1,prec=False):
    global bins,refe,mall,dall,zzall,nomal,epsil
    if bins==None:
        bins=r_[1:3:step]
    if epsil==None:
        epsil=dielect(bins,spectra.lor_si_low2[0],spectra.lor_si_low2[1])
        mall=sqrt(epsil.real+1j*levabs).real*bins
        dall=(mall[2:]-mall[:-2])/(bins[2:]-bins[:-2])
    oxide=polyval(spectra.cau_sio2,bins)
    dd=plate(bins,[epsil.real+1j*levabs,oxide,epsil],[wid1,wid2])
    if prec: 
        aall=spectra.extrema(dd,bins,poly=-1,msplit='fine',extsel=5)
        a2=array([sum(bins<a) for a in aall[isel]])
        dp=array(aall[isel])[1:]-array(aall[isel])[:-1]
    else: 
        aall=spectra.extrema(dd)
        a2=array(aall[isel])
        dp=bins[a2][1:]-bins[a2][:-1]
    b2=(a2[1:]+a2[:-1])/2
    return bins[b2],dp,dall[b2-1]

def extra_test_range(widr=[1000,1600],wid2=300,wstep=50,isel=1,fit=1,prec=False):
    rep=[]
    import extra
    for wid1 in range(widr[0],widr[1],wstep):
        b,w,n=extra_test(wid1,wid2,isel=isel,prec=prec)
        sel=extra.clipped(w*n,nsig=2.0, nclip=2)
        rep.append([wid1,(w*n)[sel].mean(),(w*n)[sel].std()])
    rep=array(rep).transpose()
    if fit:
        idx=polyfit(1/rep[1],rep[0],fit)
        chi2=((polyval(idx,1/rep[1])-rep[0])**2).sum()
        return idx,chi2
    return rep