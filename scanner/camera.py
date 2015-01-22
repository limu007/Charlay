import subprocess as sub
from matplotlib import pyplot as plt
from matplotlib import cm
from numpy import fromfile
wdir="/home/limu/Code/Bastl/miccd-0.5/example/"
import termios, fcntl
import sys,os
global dev,lev_adapt,sub_wind
dev=None
sub_wind=None
global texp
texp=1
tcon=1

def initerm():
	fd = sys.stdin.fileno()

	oldterm = termios.tcgetattr(fd)
	newattr = termios.tcgetattr(fd)
	newattr[3] = newattr[3] & ~termios.ICANON & ~termios.ECHO
	termios.tcsetattr(fd, termios.TCSANOW, newattr)

	oldflags = fcntl.fcntl(fd, fcntl.F_GETFL)
	fcntl.fcntl(fd, fcntl.F_SETFL, oldflags | os.O_NONBLOCK)
	return fd,oldterm,oldflags

#imshow(fromstring(open("sim2.dat").read(),'int16').reshape(512,768))
iname="sim2.dat"
ifile=iname

	#rep=dev.stdout.readlines()
	#print(''.join(rep))

def save_fits(ofile=None):
    import pyfits
    nl=pyfits.HDUList([pyfits.PrimaryHDU()])
    im=pyfits.ImageHDU(gcont['data'])
    nl.header["exposure"]=texp
    nl.append(im)
    if ofile==None: ofile=ifile.replace("dat","fits")
    nl.writeto(ofile)

def save_gsf(filename, imagedata,title=None,spl=1.1):
    """Write a Gwyddion GSF file.

    filename -- Name of the output file.
    imagedata -- Image data.
    xres -- Horizontal image resolution.
    yres -- Vertical image resolution.
    xreal -- Horizontal physical dimension (optional).
    yreal -- Vertical physical dimension (optional).
    xyunits -- Unit of physical dimensions (optional).
    zunits -- Unit of values (optional).
    title -- Image title (optional).

    Image data may be passed as any listable object that can be used to form
    a floating point array.array().  This includes tuples, lists, arrays,
    numpy arrays and other stuff.
    """
    import array, sys, math
    yres,xres=imagedata.shape
    xreal,yreal=None,None
    xyunits,zunits="mm",None
    
    data = array.array('f', imagedata.ravel())
    if len(data) != xres*yres:
        raise ValueError, "imagedata does not have xres*yres items"
    isinf = math.isinf
    isnan = math.isnan
    for z in data:
        if isinf(z) or isnan(z):
            raise ValueError, "GSF files may not contain NaNs and infinities"
    if sys.byteorder == 'big':
        data.byteswap()
    if xreal==None: xreal=.075/266*xres/spl
    if yreal==None: yreal=.075/266*yres
    header = ['Gwyddion Simple Field 1.0']
    header.append('XRes = %u' % xres)
    header.append('YRes = %u' % yres)
    if xreal is not None:
        header.append('XReal = %.12g' % xreal)
    if yreal is not None:
        header.append('YReal = %.12g' % yreal)
    if xyunits is not None:
        header.append('XYUnits = %s' % xyunits.encode('utf-8'))
    if zunits is not None:
        header.append('ZUnits = %s' % zunits.encode('utf-8'))
    if title is not None:
        header.append('Title = %s' % title.encode('utf-8'))

    header = ''.join(x + '\n' for x in header)
    l = len(header)
    sentinel = bytearray()
    for j in range(4 - l % 4):
        sentinel.append(0)

    fh = file(filename, 'wb')
    fh.write(header)
    fh.write(sentinel)
    fh.write(data)
    fh.close()


def sub_quad(data,bin=16,size=None):#[512,768]):
    from numpy import median,r_,newaxis,linalg,indices,ones
    if size==None: size=data.shape
    if bin>1:
        qata2=median(median(data.reshape(size[0]/bin,bin,size[1]/bin,bin),3),1)
    else:
        qata2=data
    ix,iy=indices(qata2.shape)
    mata=r_[[ones(ix.shape),ix,iy,ix**2,iy**2,ix*iy]]
    zz=linalg.lstsq(mata.reshape(6,size[0]*size[1]/bin**2).transpose(),qata2.ravel())

    jx=r_[:ix.shape[0]:1j*size[0]][:,newaxis]-0.5
    jy=r_[:ix.shape[1]:1j*size[1]][newaxis,:]-0.5
    back=zz[0][5]*(jx*jy)+zz[0][3]*(jx*jx)+zz[0][4]*(jy*jy)
    back+=zz[0][0]+zz[0][1]*jx+zz[0][2]*jy
    return back



def flatten(a,cc,rad=200,rep=0,order=2):
    from numpy import indices,linalg,ones,dot,r_
    ix,iy=indices(a.shape)
    sel=((ix-cc[0])**2+(iy-cc[1])**2)<rad**2
    jx,jy=ix[sel],iy[sel]
    mata=[ones(jx.shape),jx,jy,jx**2,jy**2,jx*jy]
    if order>=3: mata+=[jx**3,jy**3,jx**2*jy,jy**2*jx]
    zz=linalg.lstsq(r_[mata].transpose(),a[sel])
    if rep==2: return zz[0]
    mata=[ones(ix.shape),ix,iy,ix**2,iy**2,ix*iy]
    if order>=3: mata+=[ix**3,iy**3,ix**2*iy,iy**2*ix]
    qq=dot(zz[0],r_[mata].swapaxes(0,1))
    return a-qq

from threading import Thread,Event

global gcont
gcont={'data':None}
lev_adapt=0.1
debug=1

def get_lev(data,perc=0.1):
    from numpy import histogram,where
    hh=histogram(data.ravel(),100)
    sh=hh[0].cumsum()
    imin=where(sh>sh[-1]*perc)[0][0]
    imax=where(sh>sh[-1]*(1-perc))[0][0]-1
    if debug>0: print 'levels:',hh[1][imin],hh[1][imax]
    return hh[1][imin],hh[1][imax]

def implot(control):
    data=control['data']
    if sub_wind!=None:data=data[sub_wind[1]:sub_wind[3],sub_wind[0]:sub_wind[2]]
    if 'image' in control:
        control['image'].set_data(data)
        if lev_adapt>0: control['image'].set_clim(get_lev(data,lev_adapt))
    else:
        plt.clf()
        gcont['image']=plt.imshow(data,origin="lowersd")
        plt.colorbar()
    plt.show()
    plt.draw()

def adjust(cont):
    bmean=cont['back'].mean()
    dmean=cont['data'].mean()
    covar=((cont['back']-bmean)*(cont['data']-dmean)).sum()
    corel=covar/((cont['back']-bmean)**2).sum()
    cont['data']=cont['data'].astype('float32')-cont['back']*corel
    #cont['data']-=dmean

class ConThread(Thread):
    '''controlled threading 
    '''
    def __init__(self, control, exe):
        Thread.__init__(self)
        self.control=control
        self.exe=exe
    def run(self):
        self.exe(self.control)
global thplot

do_show=True
get_out=False

def shutdown():
	dev.stdin.write("TEMP 10\n")
	dev.stdin.write("END\n")
	dev.terminate()

from time import sleep

def tinput():
    global gcont,texp,lev_adapt,sub_wind
    while 1:
        inp=raw_input("cam:")
        if inp[:4]=='quit': break
        if inp[:1]=='m': 
            if os.path.exists(ifile): os.unlink(ifile)
            dev.stdin.write("MEAS\n")
            #dev.communicate()
            sleep(texp+tcon)
            if get_out:
                dev.stdout.flush()
                out=dev.stdout.read()
                print(out)
            gcont['data']=fromfile(ifile,'int16').reshape(512,768)
            if 'back' in gcont: adjust(gcont)
            implot(gcont)
            continue
        elif inp[:2]=='di': 
            if 'image' in gcont: del gcont["image"]
            continue
        elif inp[:1]=='s': 
            shutdown()
            break
        elif inp[:1]=='b': 
            if 'data' in gcont: 
                gcont['back']=sub_quad(gcont['data'])
                gcont['data']-=gcont['back']
                implot(gcont)
            continue
        elif inp[:3]=='nob':
                del  gcont['back']
        elif inp[:2]=='qt': 
            dev.stdin.write("QTEMP\n")
            continue
        larg=inp.split()
        if len(larg)<2: 
            print 'need more arguments'
            continue
        if inp[:1]=='t': 
            dev.stdin.write("TEMP "+larg[1]+"\n")
        elif inp[:3]=='pal': 
            if 'image' in gcont:
                cmap=cm.get_cmap(larg[1])
                if cmap: gcont['image'].set_cmap(cmap)
        elif inp[:1]=='e': 
            texp=int(larg[1])
            dev.stdin.write("EXP "+larg[1]+"\n")
        elif inp[:2]=='al':
            lev_adapt=float(larg[1])
        elif inp[:2]=='fi':
            save_fits(larg[1])
        elif inp[:1]=='r':
            for i in range(int(larg[1])): 
                if os.path.exists(ifile): os.unlink(ifile)
                dev.stdin.write("MEAS\n")
                sleep(texp+tcon)
                gcont['data']=fromfile(ifile,'int16').reshape(512,768)
                if 'back' in gcont: adjust(gcont)
                implot(gcont)
        if len(larg)<3: continue
        if inp[:2]=='cl':
            cols=map(float,larg[1:])
            if 'image' in gcont:
                gcont['image'].set_clim(cols[0],cols[1])
                plt.show()
        elif inp[:2]=='wi':
            size=map(int,larg[1:3])
            if len(larg)==3:
                orig=[(768-size[0])/2,(512-size[1])/2]
            else:
                orig=map(int,larg[3:5])
            sub_wind=[orig[0],orig[1],orig[0]+size[0],orig[1]+size[1]]
            print 'sub-window ',sub_wind
def inter(control):
    fd=initerm()
    try:
        while 1:
            try:
                c = sys.stdin.read(1)
                print "Got character", repr(c)
                if c == 'e': 
                	dev.stdin.write("MEAS\n")
                	#rep=dev.stdout.readlines()
                	#print(''.join(rep))h
                	data=fromfile(ifile,'int16').reshape(512,768)
                	if do_show: 
                		# if gcont['data']==None: 
                		# 	gcont['data']=data
                		# 	thplot.start()
                		# else: 
                		control['data']=data
                elif c == 's': 
                    shutdown()
                    break
            except IOError: pass
    finally:
        termios.tcsetattr(fd[0], termios.TCSAFLUSH, fd[1])
        fcntl.fcntl(fd[0], fcntl.F_SETFL, fd[2])

def proc():
    global dev,thplot
    print "running "+wdir+"miccd-test"	
    dev=sub.Popen([wdir+"miccd-test 218040304 "+iname],shell=True,stdin=sub.PIPE)#,stdout=sub.PIPE)
    if get_out:
        dev.stdout.flush()
        out=dev.stdout.read()
        print(out)
    if 'image' in gcont: del gcont["image"]
	#thplot=ConThread(gcont,inter)
    tinput()
