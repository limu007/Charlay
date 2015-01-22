# -*- coding: utf-8 -*-

'''extra tools for spectroscope
'''

global task,anal_data
from PyQt4 import QtGui, QtCore
from spectrac import save_dir#,change_sizes
from canvas import MyStaticMplCanvas,alert

anal_data={}
from numpy import log

from . import labrc as rc
from threading import Thread,Event

class ConThread(Thread):
    '''controlled threading 
    '''
    def __init__(self, control, exe):
        Thread.__init__(self)
        self.control=control
        self.exe=exe
    def run(self):
        self.exe(self.control)

class ScanDial(QtGui.QDialog):
    parent=None
    #save_dir="C:\\Data\\SOI\\"
    def __init__(self, nx=3,ny=3,r=30,*args,**kwargs):
        from PyQt4 import uic
        QtGui.QWidget.__init__(self, *args)
        uic.loadUi(rc.inst_dir+"scanning.ui", self)
        self.setWindowTitle('Scanning sample')
        self.xsize.setValue(nx)
        self.ysize.setValue(ny)
        if 'polar' in kwargs:
            self.ystep.setSingleStep(5)
            self.xstep.setValue(rc.rad_step)
            self.ystep.setValue(rc.theta_step)
            self.ch_centerf.setChecked(1)
            self.ch_anal.setChecked(rc.scan_start_analysis)
            rc.ystepsize=1.
        else:
            if r>0: self.ch_circle.setChecked(1)
            self.xstep.setValue(rc.cart_step)
            self.ystep.setValue(rc.cart_step)
            self.xsize.setSingleStep(2)
            self.ysize.setSingleStep(2)
        self.ch_return.setChecked(rc.scan_home_ret)
        if 'oname' in kwargs: self.outname.setText(kwargs['oname'])
        elif len(self.outname.text())==0: self.buttonBox.buttons()[0].setDisabled(1)
        self.radius.setValue(r)
        self.radius.setSingleStep(10)
        if rc.scan_axes_order==1 : self.ch_swap.setChecked(1)
        self.connect(self.outname, QtCore.SIGNAL("editingFinished()"),self.can_start)
#        self.connect(self.outname, QtCore.SIGNAL("textChanged()"),self.can_start)
        if 'parent' in kwargs: self.parent=kwargs['parent']
        
    def prepare_scan(self):
        ixstep=int(self.xstep.value()/rc.xstepsize)
        iystep=int(self.ystep.value()/rc.ystepsize)
        xnum=int(self.xsize.value()*self.xstep.value()/rc.xstepsize)//ixstep
        ynum=int(self.ysize.value()*self.ystep.value()/rc.ystepsize)//iystep

        swap=(self.ch_swap.checkState()>=1)

        if not hasattr(self.parent,'control'): 
            self.parent.control={'x':0,'y':0}
        else:
            self.parent.control['x']=0
            self.parent.control['y']=0
        if swap: 
            self.parent.control['swap']=True
            self.parent.control['nx'],self.parent.control['ny']=ynum,xnum
        else: self.parent.control['nx'],self.parent.control['ny']=xnum,ynum
        self.parent.control['xstep']=ixstep*rc.xstepsize
        self.parent.control['ystep']=iystep*rc.ystepsize
        
        if self.parent.polar:
            self.parent.control['init_r']=self.radius.value()
        if self.ch_circle.checkState()>=1:
            rad=int(self.radius.value()/rc.xstepsize)//ixstep
            self.parent.control['radius']=rad
            print("round sample r=%i "%rad)
            xnum=(xnum+1)//2
            ynum=(ynum+1)//2
            center=self.ch_center.checkState()>=1
        else: 
            centerfirst=self.ch_centerf.checkState()>=1
            rad=0
        #from .labin import scan_save
        #self.parent.control['return']=0
        
        if self.ch_return.checkState()>=1: self.parent.control['return']=1
        elif 'return' in self.parent.control: del self.parent.control['return']
        if self.ch_refer.checkState()>=1: self.parent.control['refer']=1
        elif 'refer' in self.parent.control: del self.parent.control['refer']
 
        oname=str(self.outname.text())
        dpos=oname.rfind('.')
        if dpos>0:
            ext=oname[dpos+1:]
            oname=oname[:dpos]
        else: ext=None
        center=False
            
        def myscan(control={}):
            print('running scan_save(%i,%i,%i,%i,%s,center=%s,rad=%i)'%(rc.base_xdir*ixstep,iystep,xnum,ynum,oname,str(center),rad))
            if swap:
                return self.parent.instr.scan_save(iystep,rc.base_xdir*ixstep,ynum,xnum,oname=oname,center=center,radius=rad,control=control,swapaxes=True)
            else:
                return self.parent.instr.scan_save(rc.base_xdir*ixstep,iystep,xnum,ynum,oname=oname,center=center,radius=rad,control=control,swapaxes=False)
        if self.parent.multitask:
            self.parent.control['meas_evt']=Event()
            return ConThread(control=self.parent.control,exe=myscan)#target=myscan)#,Thread(target=myanal)
        else:
            return myscan
    def can_start(self):
        self.buttonBox.buttons()[0].setDisabled(0)
    def select_file(self):
        from .spectrac import save_dir
        self.fileName=self.parent.fileDialog.getSaveFileName(self,self.tr("Saving data"),
            save_dir,self.tr("Text files (*)"))
        if self.fileName.isEmpty(): return
        self.outname.setText(self.fileName)#toAscii()
        import os
        #save_dir=os.path.dirname(str(self.outname.text()))
        self.can_start()
    def make_circular(self):
        cmp=self.ch_circle.checkState()==0
        self.radius.setDisabled(cmp)
        self.label_5.setDisabled(cmp)
    def move_up(self):
        move=int(self.xstep.value()/rc.xstepsize)
        print(("moving up %i steps"%move))
        rep=self.parent.instr.rate(1,move)
        if rep>=9: alert(self,'Not able to go %i steps X-dir'%move,'Motor failure',loud=1)
    def move_down(self):
        move=int(self.xstep.value()/rc.xstepsize)
        print(("moving down %i steps"%move))
        rep=self.parent.instr.rate(1,-move)
        if rep>=9: alert(self,'Not able to go %i steps -X-dir'%move,'Motor failure',loud=1)
    def move_right(self):
        move=int(self.ystep.value()/rc.ystepsize)
        print(("moving right %i steps"%move))
        rep=self.parent.instr.rate(2,-move)
        if rep>=9: alert(self,'Not able to go %i steps Y-dir'%move,'Motor failure',loud=1)
    def move_left(self):
        move=int(self.ystep.value()/rc.ystepsize)
        print(("moving left %i steps"%move))
        rep=self.parent.instr.rate(2,move)
        if rep>=9: alert(self,'Not able to go %i steps Y-dir'%move,'Motor failure',loud=1)
        

class TimeDial(QtGui.QDialog):
    parent=None
    #save_dir="C:\\Data\\SOI\\"
    def __init__(self, *args,**kwargs):
        from PyQt4 import uic
        QtGui.QWidget.__init__(self, *args)
        uic.loadUi(rc.inst_dir+"timing.ui", self)
        self.setWindowTitle('Timing measurement')
        if 'parent' in kwargs: self.parent=kwargs['parent']
    def can_start(self):
        self.buttonBox.buttons()[0].setDisabled(0)
    def select_file(self):
        from .spectrac import save_dir
        self.fileName=self.parent.fileDialog.getSaveFileName(self,self.tr("Saving data"),
            save_dir,self.tr("Text files (*)"))
        if self.fileName.isEmpty(): return
        self.outname.setText(self.fileName)#toAscii()
        import os
        #save_dir=os.path.dirname(str(self.outname.text()))
        self.can_start()
    def pulse(self):
        self.parent.instr.pulse(peri=self.periBox.value(),dur=self.totBox.value(),oname=str(self.outname.text()))
        
from math import pi,cos,sin

class ScanMap(QtGui.QMainWindow):
    parent=None
    points=[]
    xlim=None#[-10,10]
    ylim=None#[-10,10]
    polar=False
    graph=None
    butts=[]
    mode='scat'
    text_format=None
    control=None
    last=None
    def __init__(self,graph=True,attach=False,mode='scat',**kwargs):
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Scanning")
        self.main_widget = QtGui.QWidget(self)
        l = QtGui.QVBoxLayout(self.main_widget)
        hl = QtGui.QHBoxLayout()
        if 'parent' in kwargs: 
            self.parent=kwargs['parent']
            exargs=self.parent.control
            self.control=self.parent.control
            self.polar=self.parent.polar
        else: exargs={}
        if graph:
            self.graph = MyStaticMplCanvas(self.main_widget, width=5, height=4, dpi=100, polar=self.polar)
            self.graph.change_sizes()
            self.graph.axes.set_xlabel('X [mm]')
            self.graph.axes.set_ylabel('Y [mm]')
            #if self.data:
            #    self.anal.draw()
            l.addWidget(self.graph)
        self.butts=[]
        for bid in ['Period','Fluct.','Mean','Variab.','Save','Stop']:
            self.butts.append(QtGui.QPushButton(self.tr("&"+bid)))
            #self.attButton.setIcon(ico)
        for b in self.butts: hl.addWidget(b)
        self.connect(self.butts[0], QtCore.SIGNAL("clicked()"), self.show_peri)
        self.connect(self.butts[1], QtCore.SIGNAL("clicked()"), self.show_disp)
        self.connect(self.butts[2], QtCore.SIGNAL("clicked()"), self.show_mean)
        self.connect(self.butts[3], QtCore.SIGNAL("clicked()"), self.show_vari)
        self.connect(self.butts[4], QtCore.SIGNAL("clicked()"), self.save_data)
        self.connect(self.butts[5], QtCore.SIGNAL("clicked()"), self.stop_scan)

        l.addLayout(hl)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        self.statusBar().showMessage("scanning under way")
        if self.polar:
            self.xlim=None
            if 'ny' in exargs: self.ylim=[0,(exargs['ny']-0.5)*self.parent.control['xstep']]
        else:
            if 'xlim' in kwargs: exargs['nx']=kwargs['xlim']
            if 'nx' in exargs: self.xlim=[-0.5*self.parent.control['xstep'],(exargs['nx']-0.5)*self.parent.control['xstep']]
            if 'ylim' in kwargs: exargs['ny']=kwargs['ylim']
            if 'ny' in exargs: self.ylim=[-0.5*self.parent.control['ystep'],(exargs['ny']-0.5)*self.parent.control['ystep']]
        self.clear()
        self.graph.axes.hold(1)
        if self.xlim!=None: 
            self.graph.axes.set_xlim(self.xlim)
            print('setting xlim %.2f-%.2f'%tuple(self.xlim))
        if self.ylim!=None: 
            self.graph.axes.set_ylim(self.ylim)
            print('setting ylim %.2f-%.2f'%tuple(self.ylim))
        self.mode=mode
    def clear(self):
        self.graph.axes.clear()
        self.graph.draw()        
    def add_point(self,what='peri'):
        '''draw single point'''
        self.graph.axes.hold(1)
        ix,iy=self.control['x'],self.control['y']
        val=0
        if anal_data and what in anal_data: val=anal_data[what][ix,iy]
        if what=='peri': val/=1000. #in micrometers
        px,py=self.control['gx'],self.control['gy']

        if self.polar:
            py*=pi/180.
            px*=rc.xstepsize
            #px,py=py,px
            #pr,ptheta=px,py/180.*pi
            #if 'init_r' in self.control: py+=self.control['init_r']
        else:
            px*=rc.xstepsize
            py*=rc.ystepsize
            
        if 'swap' in self.control: px,py=py,px
            #px,py=pr*cos(ptheta),pr*sin(ptheta)
            #from matplotlib.pyplot import polar
            #polar(self.parent.control['y']/180.*pi,self.parent.control['x'],'sg')
        self.last=self.graph.axes.plot(px,py,'o',color=(0.5+0.1*val,0.2,0.2),markersize=rc.disp_mark_size)[0]
        #print "plotting %.3f %.3f"%(px,py)
        if self.ylim!=None: shi=(self.ylim[1]-self.ylim[0])*rc.disp_anal_shift
        else: shi=0
        if self.text_format: 
            self.graph.axes.text(px,py+shi,self.text_format%val)
            #print 'point %i %i: %s [%s]'%(ix,iy,self.text_format%val,what)
        self.points.append([px,py])
        if self.xlim!=None: self.graph.axes.set_xlim(self.xlim)
        if self.ylim!=None: self.graph.axes.set_ylim(self.ylim)
        self.graph.mpl_connect('button_release_event',self.mouse_sel)
#        self.graph.mpl_connect('button_press_event',self.mouse_sel)
        self.graph.draw()
    def save_data(self):
        #from spectra import save
        oname = self.parent.fileDialog.getSaveFileName(self,self.tr("Output file"),
                                         self.parent.saveFileNameLabel.text(),
                                         self.tr("Text Files (*.txt);;All Files (*)"))
        if str(oname)!="": save_period(str(oname),anal_data['peri'],anal_data['peri_disp'])
    def stop_scan(self):
        if not 'stop' in self.control:
            self.control['stop']=1
            self.butts[4].setText('Run')
        else:
            del self.control['stop']
            self.butts[4].setText('Stop')
    def mouse_sel(self,event):
        self.statusBar().showMessage("pos %f %f"%(event.x,event.y))
    def image_plot(self,what='peri',fact=1,desc=None):
        from numpy import array,indices
        if len(self.points)==0: ix,iy=array(anal_data['gpos']).transpose().astype('f')
        else: ix,iy=array(self.points).transpose().astype('f')
        if not what in anal_data:
            self.graph.axes.hold(0)
            self.graph.axes.plot(ix,iy,'c',color='green')
            self.graph.axes.hold(1)
            return
        fld=anal_data[what].transpose()*fact
        if (fld>0).sum()==0: return
        vals=fld[fld>0]
        vmin=vals.min()
        vmax=vals.max()
        extent=[(2*(i%2)-1)*fld.shape[i//2]//2 for i in range(4)] #+- range
        if self.ylim!=None: shi=(self.ylim[1]-self.ylim[0])*rc.disp_anal_shift
        else: shi=0
            #self.graph.axes.set_ylim(self.ylim)
        if self.mode=='scat': 
            self.graph.axes.clear()
            print('plot %i points'%len(vals))
            #if self.xlim==None: self.xlim=[ix.min(),ix.max()]
            #if self.ylim==None: self.ylim=[iy.min(),iy.max()]
            col=(vals-vmin)/2./(vmax-vmin)
            rep=self.graph.axes.scatter(ix,iy,rc.disp_mark_size*10,color=[(a+0.5,0.6-a,0.3) for a in col])
            if self.text_format:
                for i in range(len(self.points)):
                    self.graph.axes.text(ix[i],iy[i]+shi,self.text_format%vals[i])
        else: 
            rep=self.graph.axes.imshow(fld,origin="lower",vmin=vmin,vmax=vmax)
            self.graph.figure.colorbar(rep)
        if desc: self.statusBar().showMessage(desc)
        #if self.xlim!=None: self.graph.axes.set_xlim(self.xlim)
        self.graph.draw()
        return rep
    def show_peri(self):
        self.image_plot('peri',desc='Thickness estimate [um]',fact=0.001)
    def show_mean(self):
        self.image_plot('mean',desc='Mean reflectivity')
    def show_vari(self):
        self.image_plot('vari',desc='Interference amplitude')
    def show_disp(self):
        self.image_plot('peri_disp',desc='Thickness fluctuation [nm]',fact=100.)
    def show_map(self):
        self.image_plot(None)

def unsafe_map_show(self):
    def map_mon(control=None):
        if 'scanwin' in self.gui: del self.gui['scanwin']
        window = spectrax.ScanMap(polar=self.polar,parent=self,control=control)
        self.gui['scanwin']=window
        return window.show()
    from .spectrax import ConThread
    wincon=ConThread(control=self.control,exe=map_mon)

class RangeDial(QtGui.QDialog):
    '''dialog for setting spectral ranges
    '''
    def __init__(self, axes=None, *args):
        from PyQt4 import uic
        QtGui.QDialog.__init__(self, *args)
        uic.loadUi(rc.inst_dir+"randial.ui", self)
        self.setWindowTitle('Graph range')
        if axes: self.read_axes(axes)

    def read_axes(self,axes,ystep=0.05):
        ran=axes.get_xlim()
        self.xmin.setValue(ran[0])
        self.xmax.setValue(ran[1])
        ran=axes.get_ylim()
        self.ymin.setValue(ran[0])
        self.ymax.setValue(ran[1])
        self.ymin.setSingleStep(ystep)
        self.ymin.setSingleStep(ystep)

class ScriptWindow(QtGui.QMainWindow):
    fname=None
    cursor=None
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.letters = QtGui.QTabWidget()
        self.setCentralWidget(self.letters)
        
        self.setWindowTitle(self.tr("Script Form"))
        fileMenu = QtGui.QMenu(self.tr("&File"), self)
        newAction = fileMenu.addAction(self.tr("&New..."))
        newAction.setShortcut(self.tr("Ctrl+N"))
        self.printAction = fileMenu.addAction(self.tr("&Print..."), self.printFile)
        fileMenu.addAction(self.tr("&Load..."), self.load)
        fileMenu.addAction(self.tr("&Save..."), self.save)
        self.printAction.setShortcut(self.tr("Ctrl+P"))
        self.printAction.setEnabled(False)
        quitAction = fileMenu.addAction(self.tr("E&xit"))
        quitAction.setShortcut(self.tr("Ctrl+Q"))
        
        self.menuBar().addMenu(fileMenu)
        self.fullAction = self.menuBar().addAction(self.tr("&Run..."), self.execute)
        self.connect(quitAction, QtCore.SIGNAL("triggered()"), self, QtCore.SLOT("close()"))
    def createSample(self,name='example'):
        editor = QtGui.QTextEdit()
        tabIndex = self.letters.addTab(editor, name)
        self.letters.setCurrentIndex(tabIndex)
        self.cursor = editor.textCursor()
        self.cursor.movePosition(QtGui.QTextCursor.Start)
        boldFormat = QtGui.QTextCharFormat()
        boldFormat.setFontWeight(QtGui.QFont.Bold)
        self.cursor.insertText("# memo", boldFormat)
        self.cursor.insertBlock()
        self.cursor.insertText("# write here")
        
        self.printAction.setEnabled(True)
    def printFile(self):
        document = self.letters.currentWidget().document()
        printer = QtGui.QPrinter()
        dialog = QtGui.QPrintDialog(printer, self)
        #dialog.setWindowTitle(self.tr("Print Document"))
        if dialog.exec_() != QtGui.QDialog.Accepted:
            return
        document.print_(printer)
    def content(self):
        document = self.letters.currentWidget().document()
        return str(document.toPlainText())
    def execute(self,pos=0,init=True):
        comms=[a for a in self.content().split('\n') if len(a.strip())>0 and a[0]!='#']
        if init:
            print("%i commands"%len(comms))
            if len(comms)<=0: return
            self.loops=[]
            self.poops=[]
        i=pos-1
        while i<len(comms)-1:
            i+=1
            if '}' in comms[i]:
                if len(loops)==0: 
                    print('syntax error in loop')
                    continue
                else: 
                    if self.loops[-1]>0: 
                        self.loops[-1]-=1
                        i=self.poops[-1]
                    else: # end of loop
                        self.loops.pop()
                        self.poops.pop()
                        continue
            if comms[i][:6]=='repeat':
                self.loops.append(int(comms[i][6:comms[i].rfind('{')])-1)
                self.poops.append(i+1)
                #print 'looping %i at line %i'%(loops[-1],poops[-1])
                self.mummy.incounter=1
            elif comms[i][:7]=='measure': self.mummy.measure()
            elif comms[i][:5]=='calib': self.mummy.setReference()
            elif comms[i][:6]=='gotemp':
                if self.mummy.termo==None:
                    print('cannot find termostat')
                    continue
                from .labin import reach_temp
                try:
                    temp=float(eval(comms[i][6:]))
                    reach_temp(temp)
                except:
                    print('failed to reach ',str(comms[i][6:]))
            elif comms[i][:4]=='save': 
                c=comms[i]
                fname=c[c.find('"')+1:c.rfind('"')]
                if len(fname)>1:
                    #self.mummy.saveData(fname)
                    self.mummy.setSaveFileName(fname)
            elif comms[i][:5]=='print':
                expr=comms[i][5:].strip().replace('self','self.mummy').split(',')
                print(' '.join([str(eval(e)) for e in expr]))
                     
            elif comms[i][:4]=='wait':
                from time import sleep
                try:
                    c=comms[i]
                    dt=float(c[c.find('(')+1:c.rfind(')')])
                    sleep(dt)
                except:
                    continue
                # attempt for non-blocking timer
                #from PyQt4.QtCore import QTimer
                #tim=QTimer()
                # # tim.connect()
                #tim.singleShot(1000*dt,self,self.execute(pos=i,init=False))
                # break
        if i<len(comms): return i
        print("finished")
    def load(self):
        inpName = self.mummy.fileDialog.getOpenFileName(self,self.tr("Loading script"),
                save_dir,self.tr("Text Files (*.txt);;All Files (*)"))
        if inpName.isEmpty(): return
        self.fname=str(inpName.toAscii())
        print('reading from '+self.fname)
        self.setWindowTitle("Script: "+self.fname)
        ifile=open(self.fname,"r")
        self.cursor.insertText(ifile.read())
        ifile.close()
    def save(self):
        if self.fname==None:
            inpName = self.mummy.fileDialog.getSaveFileName(self,self.tr("Saving script"),
                    save_dir,self.tr("Text Files (*.txt);;All Files (*)"))
            if inpName.isEmpty(): return
            self.fname=str(inpName.toAscii())
            self.setWindowTitle("Script: "+self.fname)
        ofile=open(self.fname,"w")
        ofile.write(self.content())
        ofile.close()
        #print "not implemented yet"
    def saveas(self):
        self.fname=None
        self.save()

class AnalyseWindow(QtGui.QMainWindow):
    anal=None
    period=''
    data=None
    mode='peri'
    def __init__(self,graph=True,attach=False,mode='peri'):
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Data analysis")
        self.main_widget = QtGui.QWidget(self)
        self.mode=mode
        self.file_menu = QtGui.QMenu('&Fourier', self)
        self.menuBar().addMenu(self.file_menu)
        l = QtGui.QVBoxLayout(self.main_widget)
        hl = QtGui.QHBoxLayout()
        hl2 = QtGui.QHBoxLayout()
        
        if mode=='peri':
            self.perButton = QtGui.QPushButton(self.tr("&Period"))
            self.connect(self.perButton, QtCore.SIGNAL("clicked()"), get_period)
            hl.addWidget(self.perButton)
            self.fitButton = QtGui.QPushButton(self.tr("&Fit interval"))
            self.connect(self.fitButton, QtCore.SIGNAL("clicked()"), get_fit)
            hl.addWidget(self.fitButton)
        if attach:
            ico=QtGui.QIcon('E:\\Documents and Settings\\Kremik\\Dokumenty\\attach.png')
            self.attButton = QtGui.QPushButton(self.tr("&Attach"))
            self.attButton.setIcon(ico)
            self.connect(self.attButton, QtCore.SIGNAL("clicked()"), get_period)
            hl.addWidget(self.attButton)
        if mode=='peri':
            self.shapeButton = QtGui.QPushButton(self.tr("&Refract.ind."))
            shapeMenu = QtGui.QMenu(self.tr("refr. index"), self)
            shapeMenu.addAction(self.tr("const"),self.set_const)
            shapeMenu.addAction(self.tr("polynom"),self.set_poly)
            shapeMenu.addAction(self.tr("lorentz"))
            self.shapeButton.setMenu(shapeMenu)
            hl.addWidget(self.shapeButton)
        if mode=='scan' and len(anal_data)>0:
            self.shapeButton = QtGui.QPushButton(self.tr("&Available"))
            shapeMenu = QtGui.QMenu(self.tr("2-D maps"), self)
            for k in list(anal_data.keys()):
                if len(getattr(anal_data[k],'shape',(0,)))==2:
                    shapeMenu.addAction(self.tr("k"),self.show_map)
            self.shapeButton.setMenu(shapeMenu)
            hl.addWidget(self.shapeButton)
            self.vmin = QtGui.QLineEdit()
            self.vmax = QtGui.QLineEdit()
            hl.addWidget(QtGui.QLabel(self.tr(" min.val. ")))
            hl.addWidget(self.vmin)
            hl.addWidget(QtGui.QLabel(self.tr(" max.val. ")))
            hl.addWidget(self.vmax)
        if graph:
            self.anal = MyStaticMplCanvas(self.main_widget, width=5, height=4, dpi=100)
            self.anal.change_sizes()
            #if self.data:
            #    self.anal.draw()
            l.addWidget(self.anal)
        l.addLayout(hl)
        if mode=='peri':
            self.perEdit1 = QtGui.QLineEdit(str(self.period))
            self.perEdit2 = QtGui.QLineEdit(str(self.period))
            self.modeLabel = QtGui.QLabel(self.tr(""))
            perLabel = QtGui.QLabel(self.tr("Period estimate: min "))
            hl2.addWidget(perLabel)
            hl2.addWidget(self.perEdit1)
            perLabel = QtGui.QLabel(self.tr(" max "))
            hl2.addWidget(perLabel)
            hl2.addWidget(self.perEdit2)
            hl2.addWidget(self.modeLabel)
            l.addLayout(hl2)
            self.parLabel = QtGui.QLabel(self.tr("Fitted parameters: <b>None</b>"))
            l.addWidget(self.parLabel)
        
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        self.statusBar().showMessage("please choose method")
        global task
        task=self
        
    def set_poly(self):
        self.modeLabel.setText('Cauchy prof.')
    def set_const(self):
        self.modeLabel.setText('const.')
    #def fileQuit(self):
    def closeEvent(self,ce):
        QtGui.QMessageBox.information(self, self.tr("Closing"), QtCore.QT_TR_NOOP("Analysis saved"))
    def get_perange(self):
        rep=[str(self.perEdit1.text()),str(self.perEdit2.text())]
        for i in [0,1]:
            if rep[i]=='': rep[i]==None
            else: rep[i]=float(rep[i])
        return rep
    def show_map():
        name='gui'
        if name in anal_data:
            vmin=fld[fld>0].min()
            vmax=fld[fld>0].max()
            extent=[(2*(i%2)-1)*fld.shape[i//2]//2 for i in range(4)] #+- range
            rep=self.anal.axes.imshow(fld,extent=extent,origin="lower",vmin=vmin,vmax=vmax)
            self.graph.figure.colorbar(rep)
def get_period():
    from .spectra import perfind
    prange=task.get_perange()
    if (prange[0]==None) or (prange[1]==None):
        task.statusBar().showMessage("<b>I need period range</b>")
        return
    task.periogram=perfind(task.base,task.data,(prange[0]+prange[1])/2.,(-prange[0]+prange[1])/2.)
    if task.anal:
        task.anal.axes.cla()
        task.anal.axes.plot(task.periogram[0],task.periogram[1])
    task.parLabel.setText('maximum at %.2f'%(task.periogram[0][task.periogram[1].argmax()]))
    
def get_fit():
    from . import spectra
    if task.data==None: rep='no data available'
    else:
        a,b=spectra.fitting(task.base,task.data)
        if task.anal:
            from numpy import polyval
            task.anal.axes.cla()
            task.anal.axes.plot(task.base,spectra.sin_fin(a))
            task.anal.axes.grid()
            task.anal.axes.set_title('residual')
            task.anal.change_sizes()
            task.anal.draw()
        rep='<b>Parameters:</b>'
        rep+='<br/><i>mean value</i>: %.2f * e<sup>2</sup> + %.2f * e + %.2f'%tuple(a[5:8])
        rep+='<br/><i>period and phase</i>: e / %.4f + %.2f'%(a[3],a[4])
        rep+='<br/><i>amplitude</i>: %.2f * e<sup>2</sup> + %.2f * e + %.2f'%tuple(a[:3])
        rep+='<br/><b>Chi<sup>2</sup>:</b> %.2f'%b
    task.parLabel.setText(rep)
    task.statusBar().showMessage("converged")

def run_external(inp,pix,cal=None,exname="MINIMA2.EXE",cfname="Minima2.txt",wdir="period",tmname="line_%i",ox=0,check=None):
    '''estimates periods using external program 
    from a file containing line scan
    '''
    import os
    from numpy import polyval,savetxt,array,zeros

    odir=os.getcwd()
    ndir=os.path.join(odir,wdir)
    if not(os.path.exists(ndir)):
        return [-1]
    if type(inp)==str:
        from numpy import loadtxt
        inp=loadtxt(inp,unpack=True)
        print('got %i meas.'%len(inp))

    from .profit import reflect
    nsi=polyval(rc.mod_ext['Si'],pix)
    siref=reflect(nsi**2)

    os.chdir(ndir)
    cflines=open(cfname).readlines()[:4]
    labs=[]
    for i in range(len(inp)):
        if check and check(inp[i])<=0: continue
        if cal!=None: inp[i]/=cal
        odata=array([pix,inp[i]*siref]).transpose()
        if pix[0]>pix[-1]: odata=odata[::-1]
        savetxt(tmname%i+".ref",odata,delimiter=",",fmt="%10.5f")
        labs.append(i)
    out=zeros((len(inp),2))
    if len(labs)==0: 
        print('no good data')
        os.chdir(odir)
        return out
    if type(ox)==int: erange=rc.eranges[ox]
    else: erange=list(ox)
    cflines+=[(tmname+', %.3f, %.3f,  %.3f, Line%i\n')%(i,erange[0],erange[1],rc.ox_cor[ox],i) for i in range(len(labs))]
    cflines.append("XXXX, XXXX\n")
    open(cfname,"w").writelines(cflines)
    if not os.path.exists(exname): 
        os.chdir(odir)
        return out
    os.popen(exname).readlines()
    for i in range(len(labs)):
        if not os.path.exists(tmname%labs[i]+".od"):
            out[labs[i],1]=-1
            continue
        res=open(tmname%labs[i]+".od").readlines() # read results
        dat=[a.split() for a in res[3:-2]]
        vals=array([float(d[2]) for d in dat if len(d)>2])
        xpos=array([float(d[1]) for d in dat if len(d)>2])
        from numpy import median,std,abs
        if type(rc.thick_range)==list and len(rc.thick_range)>1:
            sel=vals>rc.thick_range[0]
            sel*=vals<rc.thick_range[1]
        elif rc.anal_use_median:
            dmed,dstd=median(vals),std(vals)
            sel=(vals-dmed)**2<dstd**2
        else: sel=vals>0

        if sum(sel)>2: 
            drep=vals[sel].min()
            dstd=vals[sel].std()
        else:
            print('only %i values passed [%f]'%(sum(sel),vals.mean()))
            drep=vals.min()
            dstd=vals[sel].std()
        out[labs[i]]=array([drep-rc.ox_cor[ox],dstd])
    os.chdir(odir)
    return out
        #3.51906+0.0423*Posi^2+0.01401*Posi^4

def save_period(save,thick,tdisp,comment=""):
        from numpy import savetxt
        ous=open(save,"w")
        ous.write("# thickness ["+comment+"] \n")
        savetxt(ous,thick,fmt="%10.6f")
        ous.write("# variability of thickness \n")
        savetxt(ous,tdisp,fmt="%10.6f")
        ous.close()

def anal_external(inlist,pix,cal=None,ox=0,save=None,lims=[0.1,0.1,0.1]):
    '''estimates periods from 2-D scan
    '''
    from numpy import array,median
    from .scan_anal2 import variab,half_anal_width
    rep=[]
    if type(pix)==str:
        if type(inlist[0])==int:
            patt=pix.replace("calib","%03i")
            inlist=[patt%i for i in inlist]
        from numpy import loadtxt
        pix,clx=loadtxt(pix,unpack=True)
        if cal==1: cal=clx
    def meas_check(d):
        if d.mean()<lims[0]: return -1 # no reflection
        if d.std()<lims[1]: return -2 # reference
        if median(variab(d))<lims[2]: return -3 # non-periodic signal
        return 1
    for f in inlist:
        rep.append(run_external(f,pix,check=meas_check,cal=cal,ox=ox))
    arep=array(rep)
    if save: save_period(save,arep[:,:,0],arep[:,:,1],
        comment="%.2f um oxide, range %.2f-%.2f eV"%(rc.ox_wid[ox],rc.eranges[ox][0],rc.eranges[ox][1]))
    return arep[:,:,0],arep[:,:,1]
    
def disp_external(init_r=None,extfile=None,rad_scale=1000.,col_scale=1.,cut_zero=True,control=None):
    '''displaying as scatter plot
    '''
    if extfile:
        from numpy import loadtxt
        anal_data['ext_disp']=loadtxt(extfile)
        nrow=len(anal_data['ext_disp'])//2
        anal_data['ext_peri']=anal_data['ext_disp'][:nrow]
        anal_data['ext_disp']=anal_data['ext_disp'][nrow:]
    # gpos not correct
    
    rad=anal_data['ext_disp'].ravel()*rad_scale
    col=anal_data['ext_peri'].ravel()*col_scale
    
    if control and 'xstep' in control:
        from numpy import arange
        ir=arange(0.,ncol)*control['xstep']
        ithet=arange(0.,nrow)*control['ystep']
    elif 'gpos' in anal_data:
        ir,itheta=array(anal_data['gpos']).astype('float').transpose()
    else: return
    itheta*=pi/180.
    if init_r: ir=abs('ir')+init_r
    ix,iy=ir*cos(itheta),ir*sin(itheta)
    from matplotlib.pyplot import scatter

    if cut_zero:
        sel=col>0
        rad=rad[sel]
        col=col[sel]
        ix=ix[sel]
        iy=iy[sel]
    scatter(ix,iy,rad,col,hold=0)
