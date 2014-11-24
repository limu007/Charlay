import pyfits
from scipy import ndimage

def shape_me(ifile="soi6um-red.fits"):
	qred=pyfits.open(ifile)[0].data

	qpos=lambda x:median(x[x>5000]) if sum(x>5000)>2 else 0
	medtst=array([[qpos(qred[i:i+50,j:j+50]) for j in range(0,750,50)] for i in range(0,500,50)])

	#rescaled map
	#medful=array([[1 if (qred[i,j]>medtst[i//50,j//50]) else 0 for j in range(0,750)] for i in range(0,500)])
	#medful[qred[0:500,0:750]<5000]=-1

	x,y=indices(medtst.shape)*50+25
	spins=interpolate.bisplrep(x[medtst>0],y[medtst>0],medtst[medtst>0],nxest=25,nyest=25)
	flat=interpolate.bisplev(r_[:500],r_[:750],spins)
	flat[qred[0:500,0:750]<5000]=4000
	qnew=qred[0:500,0:750]/flat
	qnew[qred[0:500,0:750]<5000]=0

def imanal(qnew,cmin=10,nit=1,rep=0):
	if nit<0: qsel=ndimage.binary_dilation(ndimage.binary_erosion(qnew,iterations=-nit),iterations=-nit)
	elif nit>0: qsel=ndimage.binary_erosion(ndimage.binary_dilation(qnew,iterations=nit),iterations=nit)
	qout=ndimage.label(qsel)
	if rep==1: return qout[0]
	cens=[(qout[0]==i).sum() for i in range(qout[1])]
	from numpy import r_,where
	cpos=where(r_[cens]>cmin)[0]
	print 'found %i regions -> reduced %i'%(qout[1],len(cpos))
	if rep==2: return cens,cpos
	cmax=len(cpos)
	qs=qout[0]
	j=0
	for i in range(1,qout[1]):
		if i in cpos:
			qs[qs==i]=j
			j+=1
		else:
			qs[qs==i]=cmax
	return qs


def testfile(ifile="/home/limu/Dropbox/Data/profs.pck"):
	if ifile:
		import cPickle
		wave=cPickle.load(open(ifile,"r"))
	dgrn,dred,dnir=[pyfits.open("soi3mm-"+ep+".fits")[0].data for ep in ['green','red','830']]
	wave=array([d[290:490,290:310].mean(1) for d in [dgrn,dred,dnir]])
	pnir=spectra.extrema(wave[2])[0]
	pred=spectra.extrema(wave[1])[0]

	#not flat background
def linanal():
	import spectra
	from numpy import polyfit
	pgrn=spectra.extrema(wave[0],msplit='sort')[0]
	pfit2=polyfit(pgrn,wave[0][pgrn],2) #quad. model
	y=wave[0]-polyval(pfit2,r_[:len(wave[0])])
	pgrn2=spectra.extrema(y,msplit='sort')[0]

	i=3
	pos=[pgrn2,pred,pnir]
	wlns=[532.,650.,830.]
	dir=-1
	interpt = lambda i,seq:seq[(seq>pnir[i])*(seq<pnir[i+1])]
	posk=lambda s,k,lam:(dir*(s-pnir[i])/float(pnir[i+1]-pnir[i])*830+k*830)/lam
	#sres=array([[posk(s,k,650) for s in ip] for k in range(2,20)])
	skipmax = lambda arr:arr[arange(len(arr))!=arr.argmax()]
	fracint = lambda arr:skipmax((arr-(arr+0.5).astype('int'))**2).sum() #
	gfall=[fracint(array([[posk(s,k,wlns[j]) for s in interpt(i,pos[j])] for k in range(2,20)])) for j in range(2)]

	fdif=lambda d:array([fracint((d-polyval(afit,pos[0]))/wlns[i]) for i in range(3)])
	mins=array([3373, 3764, 4388, 4897, 5342, 5816, 6409, 7003]) # nejmensi 3373, a posledni 2

	#final estimate
	afit=array([polyfit(pos[j],r_[:len(pos[j])]/2.*wlns[j],3) for j in range(3)]).mean(0) #cubic profile

def model1(pts,xsize=200,ysize=200):
	from scipy import interpolate as ip
	from numpy import mgrid,sin,indices
	gr=indices((6,6))
	#pz=sin((gr[0]**2+gr[1]**2)/2.).ravel()
	pz=exp(-((gr[0]-2)**2+(gr[1]-3)**2)/2.).ravel()
	zmodel=ip.LSQBivariateSpline(gr[0].ravel(),gr[1].ravel(),pz,r_[:6:2],r_[:6:2])

	gr2=mgrid[:4:xsize*1j,:4:ysize*1j]
	fine=zmodel.ev(gr2[0].ravel(),gr2[1].ravel())
	imshow(sin(phas*2/pi))
