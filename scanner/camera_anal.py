from numpy import r_,dot,array
def circum(p1,p2,p3):
	from numpy.linalg import inv
	m1,d1=(p2+p1),p2-p1
	m2,d2=(p3+p1),p3-p1
	rhs=r_[(m1*d1).sum(),(m2*d2).sum()] #prava strana
	lhs=2*r_[[d1,d2]]
	return dot(inv(lhs),rhs),rhs,lhs

def check_bord(zx,dmin=10.):
	from numpy import polyval,polyfit,r_
	sel=zx[:,0]>0
	ids=r_[:sum(sel)]
	while (sum(sel)>3):
		x,y=zx[sel].transpose()
		dev=abs(polyval(polyfit(y,x,2),y)-x)
		p=dev.argmax()
		if dev[p]<dmin: break
		q=ids[sel][p]
		sel[q]=False
	return sel

def get_der_edge(m,lev=1000,step=32,rep=1,do_check=True):
	from numpy import hamming,convolve
	npcs=m.shape[0]/step
	prof=m[:npcs*step].reshape(npcs,step,m.shape[1]).mean(1)
	ders=prof[:,1:]-prof[:,:-1]
	mid=m.shape[1]/2
	cv=hamming(7)
	cv/=cv.sum()
	smer=[convolve(cv,l) for l in ders]
	pos=range(step/2,m.shape[0],step)
	y=r_[pos,pos].reshape(2,len(pos)).transpose()
	pts=[(l[:mid].argmax(),l[mid:].argmin()+mid) for l in smer]
	amp=r_[[smer[i][pts[i][0]]-smer[i][pts[i][1]] for i in range(len(smer))]]
	lsel=amp<amp.max()/6.
	print sum(lsel),' lines removed'	
	z=r_[[pts,y]].swapaxes(0,2)
	z=[w[lsel==False] for w in z]
	if do_check: z=[w[check_bord(w)] for w in z]
	return z

def get_edge(m,lev=1000,step=32,rep=1,bkg=[0,40],ndsp=5,do_check=True):
	'''predpoklada fasetu vpravo/vlevo
	volny okraj vpravo/vlevo (urcen parametrem bkg) pro urceni pozadi
	'''
	from numpy import where,median
	if lev==0:
		lev=median(m[:,bkg[0]:bkg[1]])
		dsp=median(abs((m[:,bkg[0]:bkg[1]]-lev).ravel()))
		lev+=ndsp*dsp
	else:
		dsp=median(abs((m[:,bkg[0]:bkg[1]]-lev).ravel()))
	md=m>lev
	prof=md.sum(1)
	edg_pos = lambda l:where(abs(l[1:]-l[:-1]))[0][[0,-1]]
	
	npcs=prof.shape[0]/step
	pos=where(prof[:npcs*step].reshape(npcs,step).min(1)>md.shape[1]/2)[0]*step
	if rep==2: return pos
	#pos=range(100,500,50)
	if rep==3: return [(median(m[i-step/2+1:i+step/2-1,bkg[0]:bkg[1]])+ndsp*dsp) for i in pos]
	pts=[edg_pos(m[i]>(median(m[i,bkg[0]:bkg[1]])+ndsp*dsp)) for i in pos]
	y=r_[pos,pos].reshape(2,len(pos)).transpose()
	z=r_[[pts,y]].swapaxes(0,2)
	if do_check: z=[w[check_bord(w)] for w in z]
	#plot(z[0][:,0],z[0][:,1],'o')
	return z


def find_centre(z,selpt=0,side='left'):
	if side=='right': z=z[::-1]
	c1=r_[[circum(z[0][0],zx,z[1][0+selpt])[0] for zx in z[0][3:]]]
	if len(z[1])<2: return c1.mean(),c1.std() 
	c1b=r_[[circum(zx,z[0][0],z[1][-1-selpt])[0] for zx in z[0][3:]]]
	return (c1.mean(0)+c1b.mean(0))/2.,(c1.mean(0)-c1b.mean(0))/2.

def find_fasete(z,cz,rep=1,toler=5):
	from numpy import sqrt 
	s1=sqrt(((z[0]-cz)**2).sum(1))
	s2=sqrt(((z[1]-cz)**2).sum(1))
	while (toler>1):
		sel=s2<s1.mean()-toler*s1.std()
		#print s1.mean(),s1.std(),sum(sel)
		if sum(sel)>2: break
		toler-=1
	p=z[1][sel]
	if rep==2: return z[1][sel==False]
	if rep==3: return p
	d=p-p.mean(0)
	vd=d.std(0)**2
	A=(vd[1]-vd[0])/2./(d[:,1]*d[:,0]).mean()
	alp=A+sqrt(A**2+1)
	from math import atan,pi
	return 90-atan(alp)*180/pi

wkeys=['CRPIX1','CRPIX2','CDELT1']

def save_head(wlen=None):
	ff[1].header[wkeys[0]]=cz[0]
	ff[1].header[wkeys[1]]=cz[1]
	ff[1].header[wkeys[0]]=alp
	if wlen!=None: ff[1].header['WAVELENGTH']=wlen
	ff.writeto()

glob_base_disp=7
glob_fase_toler=4
def refit(z,cz=None,has_fas=True):
    from scipy import optimize
    from numpy import sqrt,concatenate,newaxis
    if cz==None: cz=find_centre(z)[0]
    if has_fas:
        rside=find_fasete(z,cz,rep=2,toler=glob_fase_toler)
        if len(rside)<3: 
            print 'only %i pts on right side of circle'%len(rside)
            if len(rside)==0: return [],[],0
        pts=concatenate([z[0],rside]).transpose()
    else:
        pts=concatenate([z[0],z[1]]).transpose()
    rad=sqrt(((pts-cz[:,newaxis])**2).sum(0).mean()) #init. values
    sqr=lambda c:(pts[0]-c[0])**2+c[3]*(pts[1]-c[1])**2-c[2]**2
    rep=optimize.leastsq(sqr,[cz[0],cz[1],rad,1.],full_output=1)[:2]
    chi2=(sqr(rep[0])**2).sum()
    if len(rep)>1 and rep[1]!=None: errs=sqrt(rep[1].diagonal())
    else: errs=[]
    return rep[0],errs,chi2

def get_scale(dd,cent,bord=0):
	''' get common region and find scale transformation
	bord cut-off (rotation relicts)
	Ex:
	for sh=[-28,162] and bord=10
	com=[dd[0][172:490,:719],dd[1][10:328,28:]]
	'''
	from numpy import r_
	sh=(cent[1]-cent[0]+0.5).astype(int)
	sz=r_[[d.shape for d in dd]] #image size

def clean1(d):
	ds=d[:250,650:].mean()
	d[250:,650:]=ds
	d[340:,600:]=ds
	d[440:,550:]=ds


def make_png(a,oname,under=None,over=1.2):
	imax=a.max()*over
	if under: 
		imin=a.min()*under
		a-=imin
		imax-=imin
	na=a.astype("float")/imax*256
	import PIL
	img=PIL.Image.new('RGB', a.shape[::-1])
	for i in range(a.shape[0]):
	  for j in range(a.shape[1]):
	    img.putpixel([j,a.shape[0]-i-1],tuple([int(na[i,j])]*3))
	img.save(oname)

def merge(inf="Lab/soi/apr14/gan_si_004_02-r1",refine=True,rep=1,msize=1000,rad=300,letlist="ab",clean=None):
	import pyfits
	from numpy import zeros,concatenate,r_,mgrid
	dd=[pyfits.open(inf+a+".fits.gz")[1].data for a in letlist]
	ina=array([[ 382.07734858,  382.7503476 ],[ 399.76129709,  195.99717147],[ 399.76129709,  195.99717147]])[:len(dd)]
	if rep==0: return dd,ina
	if clean:
		for i in range(len(dd)):
			clean(dd[i])
	if refine:
		out=[]
		for d in dd:
			z=get_edge(d,0,ndsp=glob_base_disp)
			cz=find_centre(z)[0]
			#newz=concatenate([z[0],find_fasete(z,cz,rep=2)])
			out.append(refit(z,cz)[0])
	if rep==2: return out
	from scipy import misc
	for i in range(len(dd)):
		rt=rad/out[i][2]
		nsiz=(r_[dd[i].shape]*r_[out[i][3]*rt,rt]).astype(int)
		dd[i]=misc.imresize(dd[i],nsiz)
		z=get_edge(dd[i],0,ndsp=glob_base_disp)
		cz=find_centre(z)[0]
		ang=find_fasete(z,cz)
		if abs(ang)>0.1:
			print 'rotation by *.1f deg'%ang
			dd[i]=misc.imrotate(dd[i],-ang)
			z=get_edge(dd[i],0,ndsp=glob_base_disp)
			cz=find_centre(z)[0]
		out[i]=refit(z,cz)[0]	
		ina[i]=out[i][:2]#out[i][:2]*nsiz/r_[dd[i].shape]
	print out
	if rep==3: return dd,ina
	#sel=((mgrid[:dd[0].shape[0],:dd[0].shape[1]]-a[:2,newaxis,newaxis])**2).sum(0)>a[2]**2
	dz=zeros((msize,msize))
	#mid=r_[]
	for i in range(len(dd)):
		a=ina[i]
		sel=(((mgrid[:dd[1].shape[0],:dd[1].shape[1]]-a[1::-1,newaxis,newaxis])*r_[a[3],1][:,newaxis,newaxis])**2).sum(0)>a[2]**2
		dd[i][sel]=0
		dz[a[1]:a[1]+512,a[0]:a[0]+768]=dd[i]
	nd=misc.imresize(dd[1],Out[70])
