# -*- coding: utf-8 -*-

dat2=loadtxt('/home/limu/Lab/spectra/oxidacex2.dat',unpack=True,comments='%')
rep=profit.unittest(base=dat2[0])

ham=hamming(51)
ham/=ham.sum()

calsi=profit.reflect(rep[1][1])
smeas=convolve(dat2[1][:len(calsi)]*calsi,ham)
smeas=smeas[5:][sel]

ook=profit.multifit(rep[0][sel],smeas,array(rep[1])[:,sel],[1050],dierep=-2)


e=r_[0.1:.6:800j]
win=extra.trapezoid(e,.2,.5,.1)
tl=[profit.plate(e,[12,2.25,12],[k,600]) for k in arange(4000,4500,50)]
plot(arange(4000,4500,50),[(t*win).sum() for t in tl])



#MIR measure

e,r=loadtxt("Lab/spectra/soi3/09Sep10a25.txt",unpack=True)[:,::-1]
e2=e/spectra.ev2um/1e4
rep=profit.unittest(base=e2)
widest=array([40000,1180])
plot(rep[0],profit.plate(rep[0],[rep[1][1],rep[1][0],rep[1][1]],widest))
multifit(e2[1063:],r[1063:],[rep[1][1],rep[1][0],rep[1][1]],widest))

sel=(e2>0.4)*(e2<0.5)
multifit(e2[sel],epsil[:,sel][[1,0,1]],widest)

#estimating the longer period
rsm=convolve(r,ham)[25:-25]
em=e2[1063:3423].reshape(59,40).mean(1)
rm=rsm[1063:3423:40]
plot(em,rm)

fun=lambda p:sum((p[0]+sin(em/p[1]+p[2])*p[3]-rm)**2)
res_per=0.36096
sim_per=0.30474


#NIR with infraprobe

moo=loadtxt("Lab/spectra/soi_infra/soi_50stred.dat",unpack=True,skiprows=4)
bins=1000./spectra.ev2um/r_[moo[0][0]:moo[0][-1]:369j]

nlog=repos.obslog("Lab/spectra/soi3/mir_soi.lst")
nlog.load_data(separ=None,unit="cm")
rep=profit.unittest(base=nlog.base,comp=["SiO2_gl","cSi_asp"],frange=None)
epsil=array(rep[1])

mype=[spectra.beats(nlog[i],rep=2) for i in range(len(nlog.rlist))]
min_pos=array([[nlog.base[n[0]] for n in m] for m in mype]).mean(0)

ids=range(320,1640,300)
ia,ib=ids[0],ids[1]
eall=[spectra.fitting(nlog.base[ia:ib],nlog[i][ia:ib])[0] for i in range(len(nlog.data))]
j=0
wids=[profit.per_conv/sqrt(epsil[1][ia:ib].real.mean())/eall[j][3],profit.beat_conv/(min_pos[1:]-min_pos[:-1]).mean()/epsil[0][ia:ib].real.mean()]
gene.fit=profit.multifit(nlog.base[ia:ib],nlog[j][ia:ib],epsil[:,ia:ib][[1,0,1]],wid=wids,dierep=-2)
# or rather
gene.main(pars=[nlog.base[ia:ib],nlog[j][ia:ib],epsil[:,ia:ib][[1,0,1]]],init=wids)

wids1=[40408.853755465585, 3342.5077335855963]

######### october


wids [ 1183.26501836,   609.3411821 ]
floating norm [ 0.61294005,  0.04221493]

import spectra,profit
nm=dd[0][-3:].mean(0)
mm=hamming(21);mm/=sum(mm)
nm=convolve(nm,mm,'same')
sel=(px>1.5)*(px<2.1)
sel2=r_[:2500][sel][3::6] # kind of rebinning
spectra.fitting(px0[sel2],dd[4][12][sel].reshape(121,6).mean(1)/nm[sel2],prange=[0.08,0.12])


simulations in 2D space **************

[1500+j,10] for j in range(0,300,20)

       [ -7.01087411e-07,  -3.59042058e-04],
       [  2.57922076e-06,   5.00420813e-03],
       [  6.19637053e-06,  -1.70165112e-02],
       [ -9.41553860e-05,   1.67570377e-01]

[0.0025182675537755769,
 1.6875374820457675e-05,
 -0.00014176397196424028,
 0.1772367006704608]

[1500+j,15] for j in range(0,300,20)

       [ -1.02188955e-06,  -5.59569930e-04],
       [  3.64279791e-06,   7.63932333e-03],
       [  9.48606624e-06,  -2.55853432e-02],
       [ -9.39787728e-05,   1.67396722e-01]

[0.0029318336589490143,
 3.6908078460865978e-05,
 -0.00031067263125745696,
 0.17759278016897079]

import spectra,profit
px0=r_[1:3.55:.002]
si=profit.dielect(px0,spectra.lor_si_low2[0],spectra.lor_si_low2[1])
sio2=polyval(spectra.cau_sio2,px0)
sib=si.copy()
sib.imag/=8.

wstp=range(0,200,50)
sel2=(px>1.3)*(px<2.3)
o2=[spectra.extrema(profit.plate(px,[sib,sio2,si],[1500+j,440])[sel2],px[sel2],poly=-1,msplit=True,extsel=4) for j in wstp]
omax=array([o[1][:12] for o in o2])
omin=array([o[0][:11] for o in o2])
omid=(omin[:,1:]+omin[:,:-1])/2
simid=array([profit.dielect(o,spectra.lor_si_low2[0],spectra.lor_si_low2[1]) for o in omid])
norm=sqrt(simid).real
res=(omin[:,1:]-omin[:,:-1])*norm


import spectra,profit
from numpy import load
import extra
px=load("data/oct10/pixelnew.npy")
cl=load("data/oct10/calibinv.npy")
dt=load("data/oct10/soi_oct022.npy")
sel=(px>1.42)*(px<2.4)

cl=cl[sel][::-1]
dt=dt[:,sel][:,::-1]
px=px[sel][::-1]

def anal(i):
a1=array(spectra.extrema((dt[i]/cl),px,poly=-1,msplit='fine',extsel=5)[0])
b1=(a1[1:]+a1[:-1])/2.
dn1=sqrt((b1+0.01)*profit.dielect(b1+0.01,spectra.lor_si_low2[0],spectra.lor_si_low2[1]).real+.1j)
dn1-=sqrt((b1-0.01)*profit.dielect(b1-0.01,spectra.lor_si_low2[0],spectra.lor_si_low2[1]).real+.1j)
dn1/=0.02
c1=(a1[1:]-a1[:-1])*dn1.real
c1=c1[extra.clipped(c1,2,2)]
return c1.mean(),c1.std()

desk=[t.mean(1)>1500 for t in dd]
full=[desk[i][j] and spectra.extrema(dd[i][j][sel]/cl[sel],px[sel],poly=-1,msplit='fine',extsel=10) or ([],[]) for j in range(len(dd[i]))]
