# -*- coding: utf-8 -*-

testonly=False 

######### CRYOSTAT #############
global ser,per
ser=None # serial port connection
per=None # persistent instrument
serspeed=19200 #38400
gloud=0

global meas_line #storage for measurements - made global for testing purposes
meas_line=[]

from scanner import labrc as rc
def unitconv(ufrom,uto,values):
    from spectra import ev2um
    conv,cinv=1.,0
    if ufrom==uto: 
        if values==None: return conv,cinv
        return values
    if ufrom==1: conv,cinv=ev2um*1e-3,1
    elif ufrom==2: conv,cinv=1/ev2um*1e-4,0
    if uto==1: 
        conv*=ev2um*1e-3
        cinv+=1
    elif uto==2: 
        if cinv%2==1: conv*=1/ev2um*1e-4
        else: conv/=1/ev2um*1e-4
    if values==None: return conv,cinv
    if cinv%2==1: return 1/(values*conv)
    return values*conv


######### SPECTROSCOPY #############
global gx,gy
gx,gy=0,0
from ctypes import *
class DarkCorrectionType(Structure):
    _fields_ = [("m_Enable", c_uint8),
                ("m_ForgetPercentage", c_uint8)]
                
#class DevStatType(Structure):
#    _fields_ = [("m_SerialId", c_uint16),
#                ("m_UserId", c_uint8),
#                ("m_Stat", c_uint8)]
                
class IdentType(Structure):
    _fields_ = [("m_SerialId", c_char*10),
                ("m_UserId", c_char*64),
                ("m_Stat", c_uint8)]

class SmoothingType(Structure):
    _fields_ = [("m_SmoothPix", c_uint16),
                ("m_SmoothModel", c_uint8)]
                
class ControlSettingsType(Structure):
    _fields_ = [("m_StrobeControl", c_uint16),
                ("m_LaserDelay", c_uint32),
                ("m_LaserWidth", c_uint32),
                ("m_LaserWaveLength", c_float),
                ("m_StoreToRam", c_uint16)]
                
class TriggerType(Structure):
    _fields_ = [("m_Mode", c_uint8),
                ("m_Source", c_uint8),
                ("m_SourceType", c_uint8)]
                
class MeasConfigType(Structure): #MeasConfigType
    _fields_ = [("m_StartPixel", c_uint16),
                ("m_StopPixel", c_uint16),
                ("m_IntegrationTime", c_float),
                ("m_IntegrationDelay", c_uint32),
                ("m_NrAverages", c_uint32),
                ("m_CorDynDark", DarkCorrectionType),
                ("m_Smoothing", SmoothingType),
                ("m_SaturationDetection", c_uint8),
                ("m_Trigger", TriggerType),
                ("m_Control", ControlSettingsType)]

class SimpleConfigType(Structure): #MeasConfigType
    _fields_ = [("m_StartPixel", c_uint16),
        ("m_StopPixel", c_uint16),
        ("m_IntegrationTime", c_uint16),
        ("m_NrAverages", c_uint32),
        ("m_CorDynDark", DarkCorrectionType),
        ("m_Smoothing", SmoothingType),
        ("Material", c_char*64)]

polycal=[ -3.46472e-10, -5.39531e-06, 3.44907e-01, 1.748616e+02] # fit of pixel2waveleng. dependence for AvaSpec-3648
global specdata
specdata=None

class specscope(object):
    '''generic class what every spectroscope should do 
    - by default returns simulated data (for testing purposes)
    '''
    dark=None
    flat=None
    last=None
    step=0.0008267
    pixtable=None
    config=None
    data={}
    pixrange=[0,-1]
    dirs=[0,0]

    def __init__(self,erange=[1.5,4.5],nbin=None):
        '''range in eV'''
        if nbin==None: nbin=int((erange[1]-erange[0])/self.step)
        else: self.step=float(erange[1]-erange[0])/nbin
        from numpy import arange
        self.pixtable=arange(erange[1],erange[0],-self.step)
        return
    def setup(self,prange=[1,1000],integ=100,aver=10,**kwargs):
        if self.config==None:
            self.config=SimpleConfigType()
        if prange!=None:
            if type(prange[0])==int: 
                if prange[1]>=len(self.pixtable):prange[1]=len(self.pixtable)-1
            else:
                prange=[(p>self.pixtable).sum() for p in prange]
                if prange[0]!=0:prange[0]-=1
                #prange[1]=(prange[1]<self.pixtable).sum()
                print('measuring from pix. %i to pix. %i'%tuple(prange))
            if len(self.data)>0: self.data={}
            self.config.m_StartPixel,self.config.m_StopPixel=tuple(prange)
        self.config.m_IntegrationTime=int(integ)
        self.config.m_NrAverages=int(aver)
        #if 'aver' in kwargs: self.config.m_NrAverages=int(kwargs['aver'])
        return
    def measure(self,npix=None,**kwargs):
        from numpy.random import normal
        from math import sqrt
        if 'noise' in kwargs: noise=float(kwargs['noise'])
        else: noise=getattr(self.config,'Noise',0.05)
        noise/=sqrt(self.config.m_NrAverages)
        rep=normal(size=self.pixtable.shape,scale=noise)
        mater=None
        if 'mater' in kwargs: mater=kwargs['mater']
        elif hasattr(self.config,'Material') and self.config.Material!='simu': mater=self.config.Material
        else: mater=rc.simu_mater
        if mater!=None:
            if type(mater)==list:
                for m in mater:
                    if not m in self.data:
                        if m in rc.eps_ext: 
                            from profit import dielect
                            self.data[m]=dielect(self.pixtable,rc.eps_ext[m][0],rc.eps_ext[m][1])
                        elif m in rc.mod_ext: 
                            from numpy import polyval
                            self.data[m]=polyval(rc.mod_ext[m],self.pixtable)**2
                        print('calculating diel. '+m)
                if len(mater)-len(rc.simu_layer) in [0,1]:
                    from profit import plate
                    wid=rc.simu_layer
                    if hasattr(rc,'simu_layer_fluc'): wid=[wid[i]*normal(1,rc.simu_layer_fluc[i]) for i in range(len(wid))]
                    rep+=plate(self.pixtable,[self.data[m] for m in mater],wid)
            else:
                if mater in rc.mod_ext: 
                    from numpy import polyval
                    self.data[mater]=polyval(rc.mod_ext[mater],self.pixtable)**2
                else:
                    from spectra import dbload
                    try:
                        inp=dbload(mater,self.pixtable)
                    except:
                        print('cannot load material '+mater)
                        return rep
                    #for i in range(len(inp)):
                    self.data[mater]=inp[1]
                if len(self.data)>0:
                    from profit import reflect
                    rep+=reflect(self.data[mater])#,self.data[mater[1]])
        rep*=self.config.m_IntegrationTime
        return rep
    def end(self):
        return
    def result(self,smooth=0,sub_dark=True,div_flat=True,maxit=0):
        '''complete measurement process
        repeats measurement several [maxit] times, until gets positive data
        smooth: applies software smoothing
        dark current and reference correction        
        '''
        from numpy import array
        from time import sleep
        data=None
        iata=None
        if maxit==0: maxit=rc.spectr_nb_trials
        for i in range(maxit):
            iata=self.measure()
            if iata==None: continue
            if len(iata)>0: 
                data=array(iata)
                if data.max()>0: break
            sleep(0.2)
        if data==None: 
                print('measurement failed')
                return # all attempts failed
        if smooth>0:
            data=iata[smooth/2-1::-1]+iata+iata[:-smooth/2-1:-1] #reflection on edges
            from numpy import convolve,hamming
            box=hamming(smooth)
            box/=box.sum()
            data=convolve(data,box,"same")[smooth/2:-smooth/2]
        if sub_dark and self.dark!=None: data-=self.dark
        if div_flat and self.flat!=None: 
                data/=self.flat
                if rc.hard_spec_min: data[data<rc.hard_spec_min]=rc.hard_spec_min
                if rc.hard_spec_max: data[data>rc.hard_spec_max]=rc.hard_spec_max
        self.last=data[:len(self.pixtable)]
        return self.last

    def acomm(self,comm,ntries=3,comm_gap=0.1,end_gap=0.05):
        from time import sleep
        for i in range(ntries):
            self.ard.write(comm+"\n")
            sleep(comm_gap)
            reply=''.join(self.ard.readlines())
            if reply.find("nvalid")<0: break
        if gloud>0: print('ARD: '+comm+' [%i]'%i)
        sleep(end_gap)
        return i
    def adrive(self,peri=10,gap=0):
        import serial
        if not hasattr(self,'ard'): 
            try:
                self.ard=serial.Serial(rc.ard_port,baudrate=rc.ard_speed,timeout=2)    
            except:
                print("cannot connect to Arduino")
                return
        if gap==0: gap=peri//2
        if rc.ard_mode==1: 
            rep=self.acomm('SING')
            if rep>2: print("Arduino init failed")
            else: print("single drive mode")
        if peri: self.acomm('PERI %i'%peri)
        if gap: self.acomm('GAP %i'%gap)
        rc.motor_lin_speed=rc.motor_lin_speed*peri/5.
        print("high %i ms low %i ms"%(peri-gap,gap))
        self.ard.readlines()

    def rate(self,ch=1,nstep=100,wait=0,ain=0,ntries=3,loud=0):
        '''controlling scanner 
        ch: defines movement axis
        adapted for different drives:
            Thorlabs rotation stage
            A-drive linear motor controlled by Avantes
                                            by Arduino
        '''
        from time import sleep
        global gx,gy
        if hasattr(self,'motor') and ch==2: # rotation stage by Thorlabs 
            gy+=nstep
            if testonly: return
            self.motor.relat(nstep)
            if wait>=0: sleep(rc.motor_rot_wait+rc.motor_rot_speed*abs(nstep)) # needs some time to reach the goal
            # for the moment we are not able to check the time of arrival 
            # should probably use GetPosition function
            return 0
        if hasattr(self,'ard') and self.ard!=None: # control by Arduino
            if wait<=0: wait=rc.comm_ard_wait
            if ch==1: gx+=nstep
            else: gy+=nstep
            if testonly: return
            dirtext=""
            if nstep<0:
                nstep=-nstep
                if self.dirs[ch-1]!=-1:
                    if ch==2: dirtext="LEFT"
                    else: dirtext="UP"
                    self.dirs[ch-1]=-1
            else:
                if self.dirs[ch-1]!=1:
                    if ch==2: dirtext="RIGHT"
                    else: dirtext="DOWN"
                    self.dirs[ch-1]=1
            if len(dirtext)>0: 
                self.acomm(dirtext)
                sleep(wait)
            step_so_far=0
            for i in range(ntries): # X attempts to communicate
                self.ard.write("STEP%i %i \r\n"%(ch,nstep-step_so_far))
                if gloud>0: print("STEP%i %i [ %i %s]"%(ch,nstep-step_so_far,self.dirs[ch-1],dirtext))
                wait=rc.motor_lin_wait+rc.motor_lin_speed*abs(nstep)
                sleep(wait)
                #sleep(wait)
                reply=''.join(self.ard.readlines())
                if not reply.find("finished")>0:
                    for j in range(5):
                        sleep(rc.motor_lin_wait)
                        reply+=''.join(self.ard.readlines())
                        if reply.find("finished")>0: break
                if reply.find("nvalid")>0 or loud>2: print(" || ".join(reply.split("\n")))
                if reply.find("nvalid")>0:
                    p0=reply.find("argv[0] =>")
                    if p0>0 or gloud>1: 
                        print("analyzing:",reply[p0+7:p0+17])
                    p1=reply.find("argv[1] =>")
                        
                #if reply.find("finished")>0: break
                else: break
                if loud>1: 
                    print("retrying\n")
                    print("["+reply+"]")
            return i
        return -1
        
    def write(self,fname,expername=None, comm=[],data=None,xdata=None,format='%f',type='ascii'):
        '''saving data in ascii format'''
        if data==None: data=self.last
        colab='#'
        if type=='matlab': colab='%%'
        file=open(fname,"w")
        if expername: file.write(colab+" %s\n"%expername)
        #if expernote!='': file.write(comm+" %s\n"%expernote)
        #file.write(comm+" INTEG %.1f AVG %.1f SMOOTH %i\n"%(self.config.m_IntegrationTime,self.mean,self.smooth))
        for li in comm: file.write(colab+li)
        from numpy import iterable
        if xdata==None:xdata=self.pixtable
        multicol=iterable(data[0]) #are data 1-D or more
        for i in range(len(data),0,-1):
            if type=='matlab':
                if multicol>0: lin=(format+rc.output_separ+"%f"+rc.output_separ+"%s\n")%(xdata[i-1],data[i-1][-1],rc.output_separ.join([format%d for d in data[i-1][:-1]]))
                else: lin=(format+rc.output_separ+"%f"+"\n")%(xdata[i-1],data[i-1])
            else:
                if multicol>0: lin=("%f"+rc.output_separ+"%s\n")%(xdata[i-1],rc.output_separ.join([format%d for d in data[i-1]]))
                else: lin=("%f"+rc.output_separ+format+"\n")%(xdata[i-1],data[i-1])
            file.write(lin)
        file.close()

    def pulse(self,peri=10,dur=100,oname='C:\\Users\\Lenovo\\Documents\\Data\\temp_line',loud=0):
        from datetime import datetime
        from time import sleep
        start=datetime.now()
        dela=peri/2.
        act=datetime.now()
        dis=(act-start).total_seconds()
        cnt=0
        while dis<dur:
            if loud:
                sys.stdout.write("Tick %f: %f\n"%(dis,dela))
                sys.stdout.flush()
            sleep(dela)
            act=datetime.now()
            dis=(act-start).total_seconds()
            if dis//peri>cnt:
                vals=self.result()
                cnt+=1
                dela=peri/2.
                if oname: self.write(oname+"_%03i.dat"%cnt,comm=[str(act)+"\r\n"])#,array([self.pixtable,self.last]))
                print("time count %i"%cnt)
                if self.parent: self.parent.graph.update_measure(instr=self,vals=vals)
            else:
                dela=(dis//peri+1)*peri-dis-0.5
                if dela<0: dela=0.1

    def scan_save(self,xstep=-500,ystep=500,xcnt=50,ycnt=50,radius=0,oname='C:\\Users\\Lenovo\\Documents\\Data\\quad_line',
                center=False,format='txt',control=None,swapaxes=False,centerfirst=False):
            ''' 2-d scanning of the sample - saving each row in separate file (given by oname)
            if radius>0: X-Y coverage of circular wafer
            if swapaxes: scanning second dir. first
            if centerfirst: just one point in first line (center)
            '''
            from numpy import array,save,zeros,savetxt
            from math import sqrt
            import sys
            from time import sleep

            global gx,gy
            gx,gy=0,0
            global meas_line
            if format=='txt': #hack for saving in text format
                def save(a,b):
                    ouf=open(a+".txt","w")
                    if len(b.shape)>1 and min(b.shape)>1:
                        ouf.write('#'+'  '.join(["pnt%i"%i for i in range(1,len(b)+1)])+'\n') #writing header
                    savetxt(ouf,b.transpose(),fmt="%8.5f")
                    ouf.close()
            #rate=newrate
            if swapaxes: xdir,ydir=2,1
            else: xdir,ydir=1,2
            if control:
                control['size']=(xcnt,ycnt)
                control['nx']=xcnt #is OK for round samples ??
                control['ny']=ycnt
                if 'stop' in control: del control['stop']
            if radius>0: # round sample
                if center: self.rate(ydir,-radius*ystep)
                if radius<ycnt-1: ycnt=radius+1
                if radius<xcnt-1: xcnt=radius+1
                xmin,xmax=-xcnt+1,xcnt
                ymin,ymax=-ycnt+1,ycnt
            else:
                xmin,xmax=0,xcnt
                ymin,ymax=0,ycnt
            for j in range(ymin,ymax):
                    meas_line=[]
                    k=0
                    for i in range(xmin,xmax):
                        if ((radius>0) and (i**2+j**2)>radius**2) or (centerfirst and i>xmin): 
                            meas_line.append(zeros(self.pixtable.shape))
                            # not measuring at this point - empty data
                            continue
                        if k!=0: 
                            self.rate(xdir,xstep)
                            #print('going '+str(xdir)+' '+str(xstep))
                        else: k=1  # first point in the line
                        if control and 'stop' in control: break
                        ############ measurement #######################
                        meas_line.append(self.result())
                        #print('just measured %i %i'%(i,j))
                        if control:
                            if 'wait' in control: sleep(control['wait']) 
                            control['x']=i
                            control['y']=j
                            control['gx']=gx
                            control['gy']=gy
                            if 'meas_evt' in control: control['meas_evt'].set() #starting everything that should come with new measurement
                            if 'queue' in control: control['queue'].put('measured %i %i'%(i,j)) #synchronization for GUI
                            if 'anal_evt' in control: # waiting for analysis to end
                                control['anal_evt'].wait()
                                control['anal_evt'].clear()
                            if 'stop' in control: break
                            if 'refer' in control and control['refer']==1:
                                if j==-ycnt+1: #save calibration
                                    if self.flat!=None: self.flat*=array(meas_line[-1])
                                    else: self.flat=array(meas_line[-1])
                                    self.last=1
                                    control['refer']==0
                            if 'anal' in control: control['anal']() #calling analysis directly, not multitasking
                            if 'meas_evt' in control: control['meas_evt'].clear()
                            if rc.single_central_pt and (j==0): break
                    if control and 'stop' in control: 
                        print('stop forced, leaving')
                        break #leaving
                    if j%2==1: 
                        meas_line=meas_line[::-1]
                        #print('line %i inverted'%j)
                    if radius<=0:
                        save(oname+"%03i"%(j+1),array(meas_line))
                    else:
                        save(oname+"%03i"%(j+ycnt),array(meas_line))
                        if j>=radius: break
                        if j<-radius: continue
                        self.rate(xdir,xstep*(int(sqrt(radius**2-(j+1)**2))-int(sqrt(radius**2-j**2))))
                    if j<ycnt-1:self.rate(ydir,ystep)
                    if not(rc.polar_stage and rc.single_central_pt and (j==0)): 
                        xstep=-xstep
            if control:
                if 'meas_evt' in control: 
                    control['stop']=2
                    control['meas_evt'].set()
                if 'queue' in control: control['queue'].put('finished')
                if 'return' in control: # return to the starting position
                    print('returning from pos %i, %i (%i, %i)'%(i,j,xmin,ymin))
                    self.rate(ydir,-ystep*(j-ymin))
                    ymod=rc.single_central_pt and 1 or 0
                    if (j-ymin)%2==ymod: self.rate(xdir,xstep*(i-xmin),wait=-1) #odd number of x-scans
                    control['x']=xmin
                    control['y']=ymin
                    xstep=-xstep
            if hasattr(self,'pixtable'): #saving calibration data
                if self.flat==None: save(oname+"calib",self.pixtable)
                else: save(oname+"calib",array([self.pixtable,self.flat]))

class avantes(specscope):
    '''developped and tested for AvaSpec 3648
    '''
    hnd=None
    device=None
    #dark=None
    #flat=None
    #last=None
    time=None
    ddev=None
    alldev=None
    def __init__(self,port=0,loud=0):
        '''port=0: USB, others are COM
        '''
        self.device = windll.as5216
        l_Port = self.device.AVS_Init(port)
        ndev=self.device.AVS_GetNrOfDevices()
        if ndev<1: return
        fdev = c_double(ndev)
        idev = (c_int*ndev)()
        types=(IdentType*ndev)() #types= (c_int*ndev*100)()
        rep=self.device.AVS_GetList(2*75*ndev,idev,types)
        self.alldev=types
        if rep==0: print('cannot get self.device list')
        elif loud==1: print('found device '+types[0].m_UserId)
        self.hnd=self.device.AVS_Activate(pointer(types[0]))

    def setup(self,prange=[1,3000],integ=50,aver=10,calib=None,smooth=None,dyndark=False,unit=0):
        inum=c_ushort()
        starting=True
        rep=self.device.AVS_GetNumPixels(self.hnd,pointer(inum))
        if rep<0 or inum.value<=10:
            print('invalid number of pixels')
            return rep
        from spectra import ev2um
        if calib==None: calib=rc.spec_pixtable_source
        if self.config==None:
            starting=True
            self.config=MeasConfigType()
        from numpy import arange,array,polyval
        if prange!=None:
            if calib==1:
                sdev = (c_double*inum.value)()
                self.device.AVS_GetLambda(self.hnd,pointer(sdev))
                self.pixtable=unitconv(1,unit,array(list(sdev)))
            else:
                self.pixtable=unitconv(1,unit,polyval(polycal,arange(inum.value)))
            print("converting pixels into %i"%(unit))
            if type(prange[0])!=int: # float - in physical units
                prange=[(p<self.pixtable).sum() for p in prange[::-1]]
                if prange[0]!=0:prange[0]-=1
                #prange[1]=(prange[1]<self.pixtable).sum()
                print('measuring from %i to %i'%tuple(prange))
            self.pixrange=prange
            self.ddev=None
            self.config.m_StartPixel=prange[0]
            self.config.m_StopPixel=prange[1]
            self.pixtable=self.pixtable[prange[0]:prange[1]]
            self.dark=None
            self.flat=None
        self.config.m_NrAverages=aver
        if False:  #short test measurement
            self.config.m_IntegrationTime=1
            mok=self.device.AVS_PrepareMeasure(self.hnd,pointer(self.config))
            if mok<0:
                print('first setup failed')
            else:
                self.measure()
                print('test measurement passed')
        self.config.m_IntegrationTime=integ
        if dyndark:
            self.config.m_CorDynDark.m_Enable=1
            self.config.m_CorDynDark.m_ForgetPercentage=rc.dyndark_forget
        else:
            self.config.m_CorDynDark.m_Enable=0
        if smooth!=None: self.config.m_Smoothing.m_SmoothPix=smooth
        mok=self.device.AVS_PrepareMeasure(self.hnd,pointer(self.config))
        if mok<0:
            print('setup failed')
            return mok
        self.measure()
        return inum.value
        
    def measure(self,npix=None,winhand=0,nmeas=1,prep_task=None):
        from time import sleep
        if prep_task!=None: eval(prep_task)
        if winhand<0 and hasattr(self,"parent"): winhand=self.parent
        zoo=self.device.AVS_Measure(self.hnd,winhand,nmeas)
        if zoo==-21:
            print('invalid state for spectra acquisition')
            return
        sleep(rc.sleepcon+rc.sleepfact*self.config.m_IntegrationTime*self.config.m_NrAverages)
        if npix==None: npix=self.config.m_StopPixel-self.config.m_StartPixel
        if self.time == None: self.time = c_int()
        if self.ddev == None: self.ddev = (c_double*(npix+1))()
        out=self.device.AVS_GetScopeData(self.hnd,pointer(self.time),pointer(self.ddev))
        if out==-8: print('no data measured')
        elif out==-4: print('bad handle')
        elif out<0: print('measurement failed [%i]'%int(out))
        else: return list(self.ddev[:npix])
            
    def revert(self,data=None,newbins=None):
        # gets data for uniform sampling
        from scipy import interpolate
        from numpy import arange
        if data==None: data=self.last
        if self.pixtable[0]<self.pixtable[-1]:
            bins=self.pixtable
        else:
            bins=self.pixtable[::-1]
            data=data[::-1]
        tck=interpolate.splrep(bins,data)
        if newbins==None: 
            step=(bins[-1]-bins[0])/(len(bins)-1)
            newbins=arange(bins[0],bins[-1]+step,step)
        return interpolate.splev(newbins,tck)
        
    def end(self):
        self.device.AVS_StopMeasure(self.hnd)
        self.device.AVS_Deactivate(self.hnd)
        self.device.AVS_Done()

    def do_all(self):
        if self.hnd==None: self.__init__()
        rep=setup(aver=50)
        if rep!=None:
            from numpy import array
            return array(self.measure())
    def rate(self,ch=1,nstep=100,wait=0,ain=0,ntries=3,loud=0):
        if super(avantes, self).rate(ch,nstep,wait,ain,ntries,loud)>=0: return
        # control by Avantes ports
        dir=0
        if type(nstep)==float: # if not converted before 
                nstep=int(nstep/rc.xstepsize)
                #print('going %i steps'%nstep)
        if ch==1: gx+=nstep
        else: 
            gy+=nstep
            ch=2
        if nstep<0:
            self.device.AVS_SetDigOut(self.hnd,ch,1)
            dir=1
            nstep=-nstep
        else: self.device.AVS_SetDigOut(self.hnd,ch,0)
        meas=[]
        if ain>0: dd=c_float()
                
        for i in range(nstep):
                self.device.AVS_SetDigOut(self.hnd,ch+3,1)
                if wait>0: sleep(wait)
                if ain>0: 
                        self.device.AVS_GetAnalogIn(self.hnd,ain,pointer(dd))
                        meas.append(dd)
                self.device.AVS_SetDigOut(self.hnd,ch+3,0)
                if wait>0: sleep(wait)
        if dir==1:
                self.device.AVS_SetDigOut(self.hnd,ch,0)
        if ain>0:
                from numpy import array
                return array(meas)

#-----------------------------------------------------------
            
class linuscope(avantes):
    tmpfile="/tmp/data.txt"
    pixorig=None 

    def __init__(self,comm="usbdir",loud=0):
        '''communicating through C++ interface on linux
        '''
        import os

        if hasattr(rc,'scope_pixtable'): pixfile=rc.scope_pixtable
            #default pixel position table (under linux, we don't receive pixel positions)
        else: 
            pixfile=self.tmpfile
        import time
        is_ok=False
        self.device = os.popen(comm,"w")
        for j in range(5):
            self.device.write("MEAS %s\n"%self.tmpfile);self.device.flush()
            for i in range(10):
                time.sleep(0.5)
                if os.path.exists(self.tmpfile): 
                    is_ok=True
                    break
            if is_ok and os.path.getsize(self.tmpfile)>1000: break
            else: is_ok=False
        if not is_ok: 
            print('init failed')
            self.pixtable=None
            return
        data=[a.strip().split() for a in open(pixfile).readlines() if a[0]!='#']   
        if len(data)<=1: #no measurement succeeded 
            self.pixtable=None
            return
        from numpy import array
        from spectra import ev2um
        self.pixorig=1000./ev2um/array([float(a[0]) for a in data])
        self.pixtable=self.pixorig

    def setup(self,prange=[1,1000],integ=100,aver=10,calib=0,smooth=None,dyndark=None,unit=0):
        if self.config==None: self.config=MeasConfigType()
        if prange!=None:
            if type(prange[0])==int: # really in pixels 
                if prange[1]>=len(self.pixorig):prange[1]=len(self.pixorig)-1
            else: # in eV?
                #prange=[(p>self.pixorig).sum() for p in prange]
                from numpy import where
                print("(got to convert "+str(prange)+")")
                prange=[where(p<self.pixorig)[0][-1] for p in prange]
            if prange[0]>prange[1]:prange=prange[::-1]
            if prange[0]>0: prange[0]-=1
            if abs(prange[1]-prange[0])==0: prange=[0,len(self.pixorig)-1]
            #prange[1]=(prange[1]<self.pixtable).sum()
            print('measuring from pix. %i to pix. %i'%tuple(prange))
            if len(self.data)>0: self.data={}
            self.pixrange=prange
            self.config.m_StartPixel=prange[0]
            self.config.m_StopPixel=prange[1]
            self.pixtable=self.pixorig[prange[0]:prange[1]] #convert to right units?
            self.pixtable=unitconv(0,unit,self.pixtable)
            self.dark=None
            self.flat=None
        if dyndark!=None:
            if dyndark==True: dyndark=50
            if type(dyndark)==int and dyndark>0:
                self.config.m_CorDynDark.m_Enable=1
                self.config.m_CorDynDark.m_ForgetPercentage=dyndark
                self.device.write("DYD %i\n"%50);self.device.flush();
            else:
                self.config.m_CorDynDark.m_Enable=0
                self.device.write("DYD %i\n"%0);self.device.flush();
        if smooth!=None: 
            self.config.m_Smoothing.m_SmoothPix=smooth
            self.device.write("SMO %i\n"%smooth);self.device.flush();
        self.config.m_IntegrationTime=int(integ)
        if aver>0: 
            self.config.m_NrAverages=aver
            self.device.write("AVG %i\n"%self.config.m_NrAverages);self.device.flush();
        if integ>0: 
            self.config.m_IntegrationTime=integ
            self.device.write("EXP %i\n"%integ);self.device.flush();
        return


    def measure(self,npix=None,**kwargs):
        import os,time
        if os.path.exists(self.tmpfile): otime=os.path.getmtime(self.tmpfile)
        else: otime=0
        self.device.write("MEAS %s\n"%self.tmpfile);self.device.flush()
        if self.config: time.sleep((self.config.m_IntegrationTime-1)/1000.)
        for i in range(30):
            if os.path.exists(self.tmpfile) and otime<os.path.getmtime(self.tmpfile): break
            time.sleep(0.3)
        if not os.path.exists(self.tmpfile): return None
        for i in range(3):
        	time.sleep(rc.linux_write_delay)
        	data=[a.strip().split() for a in open(self.tmpfile).readlines() if a[0]!='#'] 
        	if len(data)==len(self.pixorig): break   
        self.ddev=[float(a[-1]) for a in data if len(a)==2][self.pixrange[0]:self.pixrange[1]]
        return self.ddev

    def end(self):
        self.device.write("END\n");self.device.flush()
        self.device.close()
        
    def check_bord(self,sat_lev=14000):
        #border reached
        from numpy import array
        data=self.last
        if sum(data>sat_lev)>10: 
            print('saturating')
            self.config.m_IntegrationTime=int(self.config.m_IntegrationTime*0.8)

    

class linkam():
    temprof=None
    def __init__(speed=-1,gotemp=None,ramp=30,init=True,close=False,port=0):
        import serial
        global ser
        if ser==None: ser=serial.Serial(port,baudrate=serspeed,timeout=2) # com1 communication
        if init: 
            ser.write("\x3f\xa0\r")
            ser.readline()
        if gotemp: # not responding so far
            self.set_temp(speed,gotemp,ramp)
        if close: 
            ser.close()
            ser=None
        #return temp
    def get_temp(self):
        ser.write("T\r")
        out=ser.readline()
        try:
            temp=int(out[-5:-1],16)*0.1
        except:
            print('error:'+out)
            return 
        return temp
    def set_temp(self,speed=-1,gotemp=0,ramp=30):
        if speed>=0 and speed<=31:
            ser.write("P"+chr(48+speed)+"\r")
        ser.write("L1%03i\r"%int(gotemp*10))
        ser.write("R1%04i\r"%int(ramp*100))
        ser.write("Pa0\r")


    def reach_temp(self,temp,toler=1,per=2.):
        from time import sleep,time
        aa=linkam(gotemp=temp)
        at=aa.get_temp()
        ot=at
        t0=time()
        self.temprof=[]
        while abs(at-temp)>toler:
            lt=at
            sleep(per)
            aa.set_temp(gotemp=temp)
            at=aa.get_temp()
            if abs(lt-at)<toler/5.:
                print("not working")
                break
            t1=time()-t0
            if (lt-at)/per/(ot-at)*t1<0.1: 
                print("cooling too slow")
                break
            self.temprof.append((t1,at))
        return at

