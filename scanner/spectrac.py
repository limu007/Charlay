#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os, random
from PyQt4 import QtGui, QtCore
try:  
    from PyQt4.QtCore import QString,QStringList
except ImportError:  
    # we are using Python3 so QString is not defined  
    QString = str
    QStringList = list
    
from numpy import arange, array, sin, pi
global scolor,sstyle,lcolor

sstyle='-,--,-.,:,.'.split(',')
scolor='bgrmkcy'
lcolor=''
#others in cm.colors.cnames
scmap=dict([('b','blue'),
('g','green'),
('r','red'),
('c','cyan'),
('m','magenta'),
('y','yellow'),
('k','black'),
('w','white')])

try:
    from scanner import labin
    from scanner import labrc as rc
except:
    from . import labin
    from . import labrc as rc

#import labin # laboratory instruments

progname = 'Avantes new cloth'
progversion = "0.97"

global aw
aw=None

units=['eV','nm','cm-1']
uquantity=['energy','wavelength','wavecount']

global save_dir#,oxide_sel
save_mode=3
#save_dir="C:\\Data\\SOI"
save_dir="/home/limu/Lab"
use_span='horizontal'
auto_init=False

try:
    from canvas import MyStaticMplCanvas,MyDynamicMplCanvas,MyNavigationToolbar,GraphWindow,checkItem,SpanSelector
    from canvas import alert,extend_alert
except:
    from scanner.canvas import MyStaticMplCanvas,MyDynamicMplCanvas,MyNavigationToolbar,GraphWindow,checkItem,SpanSelector
    from scanner.canvas import alert,extend_alert

try: 
    import spectrax as sx
except:
    from scanner import spectrax as sx
        
#############################################################################

class ApplicationWindow(GraphWindow):
    #instr=None
    termo=None
    running=False
    multitask=False
    polar=True
    mean=100
    intime=10
    deltime=0.5
    smooth=0 #orig. 11

    gui={}
    control={}
    expername=""
    expernote=""
    table=None
    randial=None
    incounter=0
    keep=0
    report=None

    def __init__(self):
        global lcolor,save_dir

        if rc.lang=='cs':
                tran=QtCore.QTranslator()
                tran.load(rc.inst_dir+"spectrac.qm")
                qApp.installTranslator(tran)
                print('setting language to '+rc.lang)

        GraphWindow.__init__(self)
        self.setWindowIconText("PySpec")
        
        self.smooth=rc.smoothing
        self.mean=rc.mean
        self.intime=rc.intime
        self.polar=rc.polar_stage
        
        self.file_menu = QtGui.QMenu(self.tr('&File'), self)
        self.file_menu.addAction(self.tr('&Experiment name'), self.setExper, QtCore.Qt.CTRL + QtCore.Qt.Key_E)
        self.file_menu.addAction(self.tr('&Load data'), self.loadData, QtCore.Qt.CTRL + QtCore.Qt.Key_L)
        self.file_menu.addAction(self.tr('&Save data'), self.setSaveFileName, QtCore.Qt.CTRL + QtCore.Qt.Key_S)
        self.file_menu.addAction(self.tr('Load &reference'), self.loadRefer, QtCore.Qt.CTRL + QtCore.Qt.Key_R)
        self.file_menu.addAction(self.tr('Save re&ference'), self.saveRefer)
        self.file_menu.addAction(self.tr('Load line scan'), self.loadMulti, QtCore.Qt.CTRL + QtCore.Qt.ALT + QtCore.Qt.Key_L)
        self.file_menu.addAction(self.tr('&Save Excel'), self.saveExcel, QtCore.Qt.CTRL + QtCore.Qt.Key_X)
        #self.file_menu.addAction(self.tr('&Upload'), self.upload)
        self.file_menu.addAction(self.tr('&Quit'), self.fileQuit, QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        if rc.scan_realtime: self.multitask=True

        self.set_menu = QtGui.QMenu(self.tr('&Setup'), self)
        self.set_menu.addAction(self.tr('&Init spectrometer'), self.initSpec, QtCore.Qt.CTRL + QtCore.Qt.Key_I)
        self.set_menu.addAction(self.tr('&Measure'), self.measure, QtCore.Qt.CTRL + QtCore.Qt.Key_M)
        self.set_menu.addAction(self.tr('&Repeat'), self.repeat, QtCore.Qt.CTRL + QtCore.Qt.Key_R)
        self.set_menu.addAction(self.tr('&Use as reference'), self.useRefer), 
        self.set_menu.addAction(self.tr('&Graph range'), self.grange, QtCore.Qt.CTRL + QtCore.Qt.Key_G)
        self.set_menu.addAction(self.tr('&Smoothing'), self.smoothing)
        self.set_menu.addAction(self.tr('Re&binning'), self.set_binning)
        self.set_menu.addAction(self.tr('S&how stack'), self.show_stack, QtCore.Qt.CTRL + QtCore.Qt.Key_K)
        self.set_menu.addAction(self.tr('&Clear'), self.clear)
        if rc.scan_extra_menu_experim:
            self.set_menu.addAction(checkItem(self.tr("Rotation stage"),self.changePolar,self,self.polar))
            self.set_menu.addAction(checkItem(self.tr("MultiTasking"),self.changeMulti,self,self.multitask))
            self.set_menu.addAction(self.tr('Init servo motors'), self.rotInit)
            self.set_menu.addAction(self.tr('Home servo motors'), self.rotHome)
        self.set_menu.addAction(self.tr('Motor speed'), self.ardSpeed)
        self.set_menu.addAction(self.tr('Scri&pted'), self.scripted)
        self.menuBar().addMenu(self.set_menu)
        
        self.tool_menu = QtGui.QMenu(self.tr('&Tools'), self)
        self.tool_menu.addAction(self.tr('&Analysis'), self.showExtra, QtCore.Qt.CTRL + QtCore.Qt.Key_A)
        self.tool_menu.addAction(self.tr('&Variability'), self.variance)
        self.tool_menu.addAction(self.tr('Tim&ing'), self.timeMeas)
            
        if rc.scan_extra_menu_experim:
            self.tool_menu.addAction(self.tr('&Scan'), self.scanSample, QtCore.Qt.CTRL + QtCore.Qt.Key_N)
            self.tool_menu.addAction(self.tr('&Thickness (single.meas)'), self.getPeri, QtCore.Qt.CTRL + QtCore.Qt.Key_T)
        if rc.scan_extra_menu_experim==2:
            self.tool_menu.addAction(self.tr('Scan thickness (external)'), self.scanPeri)
            self.tool_menu.addAction(self.tr('&Mapping (experimental)'), self.showExtraMap)
        self.menuBar().addMenu(self.tool_menu)

        self.help_menu = QtGui.QMenu(self.tr('&Help'), self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)
        self.help_menu.addAction(self.tr('&About'), self.about)
        
        frameStyle = QtGui.QFrame.Sunken | QtGui.QFrame.Panel
        self.rangeLabel = QtGui.QLabel("   ")
        self.rangeLabel.setFrameStyle(frameStyle)
        self.runLabel = QtGui.QLabel()
        self.runLabel.setFrameStyle(frameStyle)
        self.rangeButton = QtGui.QPushButton(self.tr("Set ran&ge"))
        self.termoLabel = QtGui.QLabel()

        if use_span:
            self.span = SpanSelector(self.graph.axes, self.onselect, use_span, minspan =0.01,
                     useblit=True, rectprops=dict(alpha=0.5, facecolor='#C6DEFF') )
        
        self.saveFileNameLabel = QtGui.QLabel()
        self.initButton = QtGui.QPushButton(self.tr("&Init"))
        self.saveButton = QtGui.QPushButton(self.tr("&Save"))
        self.darkButton = QtGui.QPushButton(self.tr("&Dark"))
        self.referButton = QtGui.QPushButton(self.tr("&Refer"))
        self.referLabel = QtGui.QLabel()
        self.runButton = QtGui.QPushButton(self.tr("Ru&n"))
        self.snapButton = QtGui.QPushButton(self.tr("Sna&p"))
        self.measButton = QtGui.QPushButton(self.tr("&Measure"))
        self.grkeep = QtGui.QCheckBox(self.tr("Keep graph"));
        self.grkeep.setChecked(False)
        self.incname = QtGui.QCheckBox(self.tr("Inc.name"));
        self.incname.setChecked(False)
        
        self.connect(self.initButton, QtCore.SIGNAL("clicked()"), self.initSpec)
        self.connect(self.snapButton, QtCore.SIGNAL("clicked()"), self.snap)
        self.connect(self.saveButton, QtCore.SIGNAL("clicked()"), self.saveLast)
        self.connect(self.darkButton, QtCore.SIGNAL("clicked()"), self.setDark)
        self.connect(self.referButton, QtCore.SIGNAL("clicked()"), self.setReference)
        self.connect(self.rangeButton, QtCore.SIGNAL("clicked()"), self.setRange)
        self.connect(self.runButton, QtCore.SIGNAL("clicked()"), self.run)
        self.connect(self.measButton, QtCore.SIGNAL("clicked()"), self.measure)
        
        if labin.per==None: #persistent instrument?
            self.referButton.setDisabled(1)
            self.measButton.setDisabled(1)
            self.darkButton.setDisabled(1)
            self.runButton.setDisabled(1)
            self.saveButton.setDisabled(1)
        else:
            self.instr=labin.per
            self.initButton.setDisabled(1)
            self.instr.setup(None,integ=self.intime,aver=self.mean,smooth=self.smooth)
            
        self.darkcheck = QtGui.QCheckBox(self.tr("dyn.dark"));
        self.darkcheck.setChecked(rc.use_dyndark)
        self.avgcnt = QtGui.QSpinBox()
        self.avgcnt.setValue(self.mean)
        self.avgcnt.setRange(1,1000)
        self.avgcnt.setSingleStep(10)

        self.integbox = QtGui.QSpinBox()
        self.integbox.setValue(self.intime)
        self.integbox.setRange(1,1000)
        self.integbox.setSingleStep(10)
        self.integbox.setSuffix(" ms")

        self.connect(self.grkeep, QtCore.SIGNAL("toggled(bool)"),self.set_keep)
        self.connect(self.incname, QtCore.SIGNAL("toggled(bool)"),self.incremental)
        self.connect(self.darkcheck, QtCore.SIGNAL("toggled(bool)"),self.darkening)
        self.connect(self.avgcnt, QtCore.SIGNAL("valueChanged(int)"), self.averaging)
        self.connect(self.integbox, QtCore.SIGNAL("valueChanged(int)"), self.averaging)

        h = QtGui.QHBoxLayout(self.main_widget)
        l = QtGui.QVBoxLayout()
        self.expLabel = QtGui.QLabel(self.tr("Comment:"))
        self.expEdit = QtGui.QLineEdit(self.expernote)
        self.connect(self.expEdit, QtCore.SIGNAL("editingFinished()"),self.saveNote)

        self.pbar = QtGui.QProgressBar()
        
        hl = QtGui.QHBoxLayout()
        hl2 = QtGui.QHBoxLayout()
        hl3 = QtGui.QHBoxLayout()
        hl4 = QtGui.QHBoxLayout()
        if not auto_init: hl.addWidget(self.initButton)
        hl4.addWidget(self.darkButton)
        hl4.addWidget(self.referButton)
        hl.addWidget(self.measButton)
        hl.addWidget(self.saveFileNameLabel)
        hl.addWidget(self.saveButton)
        hl.addWidget(self.grkeep)
        hl.addWidget(self.incname)
        hl3.addWidget(self.expLabel)
        hl3.addWidget(self.expEdit)
        hl3.addWidget(self.termoLabel)
        hl3.addWidget(self.pbar)
        hl2.addWidget(self.rangeButton)
        hl2.addWidget(self.rangeLabel)
        hl2.addStretch(10)
        hl2.addWidget(self.darkcheck)
        hl2.addStretch(10)
        hl2.addWidget(QtGui.QLabel(self.tr(" averaging")))
        hl2.addWidget(self.avgcnt)
        hl2.addSpacing(10)
        hl2.addWidget(QtGui.QLabel(self.tr("integration time")))
        hl2.addWidget(self.integbox)
        
        hl.addWidget(self.runButton)
        hl.addWidget(self.runLabel)
        hl.addWidget(self.snapButton)
        
        self.toolbar = MyNavigationToolbar(self.graph, self.graph)
        
        #l.addLayout(self.toolbar.layout)
        clearButton = QtGui.QPushButton(self.tr("Vymazat"))
        self.connect(clearButton, QtCore.SIGNAL("clicked()"), self.clear)

        self.vertrangeButton = QtGui.QPushButton(self.tr("Adjust &vert."))
        self.connect(self.vertrangeButton, QtCore.SIGNAL("clicked()"), self.adjust_vert)
        self.resetrangeButton = QtGui.QPushButton(self.tr("&Reset"))
        self.connect(self.resetrangeButton, QtCore.SIGNAL("clicked()"), self.reset)
        self.analrangeButton = QtGui.QPushButton(self.tr("Set &anal.range"))
        self.connect(self.analrangeButton, QtCore.SIGNAL("clicked()"), self.setAnalRange)
        self.showanal = QtGui.QCheckBox(self.tr("show"));
        self.showanal.setChecked(rc.show_default_range)
        self.connect(self.showanal, QtCore.SIGNAL("toggled(bool)"), self.analplot)
        self.showanal.setDisabled(1)
        if self.toolbar:
            self.toolbar.addWidget(clearButton)
            self.toolbar.addWidget(self.vertrangeButton)
            self.toolbar.addWidget(self.resetrangeButton)
            self.toolbar.addWidget(self.analrangeButton)
            self.toolbar.addWidget(self.showanal)

        self.oxideModes = QtGui.QComboBox()
        self.oxideModes.addItem('estim.')
        for i in range(len(rc.ox_wid)):
            self.oxideModes.addItem('oxid %.1f um'%rc.ox_wid[i])
        self.oxideModes.setCurrentIndex(rc.oxide_sel)
        self.unitModes = QtGui.QComboBox()
        for u in units:
            self.unitModes.addItem(u)
        self.connect(self.unitModes,QtCore.SIGNAL('currentIndexChanged(int)'),self.setUnit)
        self.unitModes.setCurrentIndex(0)
        unitLabel = QtGui.QLabel(self.tr("Unit"))
        hl4.addWidget(unitLabel)
        hl4.addWidget(self.unitModes)

        self.referModes = QtGui.QComboBox()
        self.referModes.addItem("Unity")
        for k in rc.ref_tables.keys():
            self.referModes.addItem(k)
        referLabel = QtGui.QLabel(self.tr("Reference sample"))
        self.connect(self.referModes,QtCore.SIGNAL('currentIndexChanged(int)'),self.setCalib)
        referLabel.setBuddy(self.referModes)
        self.referModes.setCurrentIndex(0)
        hl4.addStretch(10)
        hl4.addWidget(referLabel)
        hl4.addWidget(self.referModes)        

        if rc.scan_extra_buttons:
            self.stackButton = QtGui.QPushButton(self.tr("Sta&ck"))
            self.connect(self.stackButton, QtCore.SIGNAL("clicked()"), self.show_stack)
            hl.addWidget(self.stackButton)
            self.thickButton = QtGui.QPushButton(self.tr("&Thickness"))
            self.connect(self.thickButton, QtCore.SIGNAL("clicked()"), self.getPeri)
            hl4.addWidget(self.thickButton)

            hl4.addWidget(self.oxideModes)
            self.connect(self.oxideModes,QtCore.SIGNAL('currentIndexChanged(int)'),self.setWorkOxide)

        l.addLayout(hl)
        l.addLayout(hl2)
        l.addWidget(self.graph)
        if self.toolbar: l.addWidget(self.toolbar)
        #l.addWidget(self.span)
        l.addWidget(self.calib)
        l.addLayout(hl4)
        
        h.addLayout(l)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        if not auto_init: self.statusBar().showMessage("please press Init!")
        self.gui['toplay']=l
        self.gui['horlay']=h
        self.gui['butlay']=hl
        self.gui['butlay2']=hl2
        self.gui['commlay']=hl3
        self.gui['butlay4']=hl4
        self.fileDialog=QtGui.QFileDialog()
        self.graph.analbar=None
        
        if rc.set_maximized: self.showMaximized()
        if rc.date_dir!=None:
            from datetime import datetime
            save_dir=rc.date_dir+os.path.sep+datetime.strftime(datetime.now(),"%Y-%m-%d")
            print("saving to "+save_dir)
        #from splot import lcolor as scolor
        lcolor=list(scolor)
        lcolor.remove(rc.line_color)
        if auto_init: self.initSpec(instrname=None)
        try:
            from queue import Queue
        except:
            from Queue import Queue
        self.report=Queue()
        
    def setUnit(self,imode,mode=None,unit=None):
        '''changing x-axis units in all graphs
        modes : 
            0=electronvolts
            1=nanometers
            2=wavecount(1/cm)
        '''
        if mode==None: mode=str(self.unitModes.currentText())
        if self.xunit==imode: return
        conv,cinv=labin.unitconv(self.xunit,imode,None)
        self.xunit=imode
        xaxis=None
        if self.graph.hscale==None: self.graph.hscale=self.graph.axes.get_xlim()
        #self.graph.axes.cla()
        #self.graph.change_sizes()
        from spectra import get_ids
        if self.instr!=None:
            if xaxis==None:
                #if self.graph.hscale==None: self.graph.hscale=self.graph.axes.get_xlim()
                xaxis=self.instr.pixtable
                xids=get_ids(xaxis,self.graph.hscale)
                if cinv%2==1: xaxis=1/(xaxis*conv)
                else: xaxis*=conv
            if cinv%2==1: self.instr.pixtable=1/(self.instr.pixtable*conv)
            else: self.instr.pixtable*=conv
        for k in self.stack.keys():
            if xaxis==None:
                xaxis=self.stack[k].get_xdata()
                xids=get_ids(xaxis,self.graph.hscale)
                print('range set to %s'%str(xids))
                if cinv%2==1: xaxis=1/(xaxis*conv)
                else: xaxis*=conv
            if len(self.stack)>0:
                self.stack[k].set_xdata(xaxis)#=self.graph.axes.plot(xaxis,self.stack[k].get_ydata())[0]
        if self.graph.last!=None: self.graph.last.set_xdata(xaxis)
        if xaxis==None: 
            print('setting unit: no xaxis defined')
            return
        xids=[i>=len(xaxis) and len(xaxis)-1 or i for i in xids]
        if abs(xids[0]-xids[1])>10:
            #if (xids[0]>=len(xaxis)) or (xids[1]>=len(xaxis)):
            #    print('no xlim updating: values %i / %i overflow'%tuple(xids))
            #else:
            self.graph.hscale=(xaxis[xids[0]],xaxis[xids[1]])
            if self.graph.hscale[0]>self.graph.hscale[1]: self.graph.hscale=self.graph.hscale[::-1]
            self.graph.axes.set_xlim(*self.graph.hscale)
        else:
            print('no xlim updating: values %i / %i too close'%tuple(xids))
            self.graph.hscale=None
        if self.instr.flat!=None:
            self.calib.axes.plot(self.instr.pixtable,self.instr.flat)
            if self.graph.hscale!=None: self.calib.axes.set_xlim(*self.graph.hscale)
            self.calib.change_sizes()
            self.calib.draw()
        self.graph.xlabel=uquantity[self.xunit]+' ['+units[self.xunit]+']'
        self.graph.axes.xaxis.set_label_text(self.graph.xlabel)
        self.graph.axes.grid(1)
        self.graph.draw()


    def setCalib(self,imode,mode=None,unit=None):
        '''set absolute/relative calibration of spectra
        '''
        if mode==None: mode=str(self.referModes.currentText())
        print('setting reference table to %s'%mode)
        slist=list(self.stack.values())
        if self.graph.last!=None: slist+=[self.graph.last]
        if not (mode in rc.ref_tables): 
            if self.table!=None:
                for s in slist:
                    vals=s.get_ydata()
                    t=(len(self.table)-len(vals))/2
                    if t>=0: s.set_ydata(vals/self.table[t:t+len(vals),1])
                self.graph.draw()
            self.table=None
            return
        if self.table!=None: oldtable=self.table
        else: oldtable=None
        from numpy import loadtxt,sqrt
        import os
        if not os.path.exists(rc.ref_tables[mode]):
            print('cannot open '+rc.ref_tables[mode])
            return
        self.table=loadtxt(rc.ref_tables[mode])
        if self.table.shape[1]==3: #complex permeability
            n=sqrt(self.table[:,1]+1j*self.table[:,2])
            r=(n-1)/(n+1) #normal incidence
            self.table[:,1]=r*r.conj()
        if (unit==None) and sum(self.table[:,0]>10)>10: 
            unit='nm'
            print('assumed calibration function in nanometers')
        if unit=='nm': 
            from spectra import ev2um
            self.table[:,0]=uniconv(1,self.xunit,self.table[:,0])
        if self.instr==None: return
        else: xbase=self.instr.pixtable
        from spectra import get_ids
        ids=get_ids(self.table[:,0],xbase[[0,-1]]) #works only if pixtable in eV
        self.table=self.table[ids[0]:ids[1]:rc.ref_sample]
        if len(self.table)!=len(xbase):
            print('using %i bins'%len(self.table))
            from numpy import array
            if len(self.table)<0.5*len(xbase):
                print('need to interpolate')
                from scipy import interpolate
                if self.table[0,0]>self.table[-1,0]: dir=-1
                else: dir=1
                try:
                    tck=interpolate.splrep(self.table[::dir,0],self.table[::dir,1])
                    redim=interpolate.splev(xbase[::dir],tck)[::dir]
                except:
                    print('interpolation failed [dir %i]'%dir)
                    return mode
            else:
                from scipy import ndimage
                step=abs((xbase[:-1]-xbase[1:])).min()
                zok=(xbase/step).astype(int)
                mok=(self.table[:,0]/step).astype(int)
                redim=ndimage.mean(self.table[:,1],mok,zok)
            self.table=array([self.instr.pixtable,redim]).transpose()
        if oldtable!=None: conv=self.table/oldtable
        else: conv=self.table
        for s in slist:
            vals=s.get_ydata()
            t=(len(conv)-len(vals))/2
            if t>0: s.set_ydata(vals*conv[t:t+len(vals),1])
        self.graph.draw()
        return mode
    def setRange(self, do_conv=False):    
        suff=uquantity[self.xunit]+" ["+units[self.xunit]+"]:"
        from spectra import ev2um
        from numpy import array
        dd=array([self.emin,self.emax])
        elabel=["lower ","upper "]
        dd=labin.unitconv(0,self.xunit,dd)
        if units[self.xunit]=='nm': dd=dd[::-1]
        d, ok = QtGui.QInputDialog.getDouble(self, self.tr("get range"),self.tr(elabel[0]+suff), dd[0], 0, 10000, 2)
        if ok: dd[0]=float(d)
        d, ok = QtGui.QInputDialog.getDouble(self, self.tr("get range"),self.tr(elabel[1]+suff), dd[1], 0, 10000, 2)
        if ok:
            #if float(d)<self.emin: d=self.emin+10
            dd[1]=float(d)
            if units[self.xunit]=='nm': dd=dd[::-1]
            #	dd=[1000./d/ev2um for d in dd]
            if do_conv: dd=labin.unitconv(self.xunit,0,dd)
            self.emin,self.emax=dd
            if self.instr!=None:
                self.instr.setup([self.emin,self.emax],integ=self.intime,aver=self.mean,smooth=self.smooth,unit=self.xunit)
            self.rangeLabel.setText("%s-%s"%tuple([str(s) for s in dd])+units[self.xunit])
            #self.rangeLabel.setText(QString("%s-%s"%tuple([str(s) for s in dd])+units[self.xunit]))
        if self.instr!=None:self.instr.measure()
        self.referLabel.setText("dark: none refer: none")
        self.statusBar().showMessage("setup saved")
    
    #---------------------- saving and loading
    
    def progress():
        '''showing progressbar during measurement
        not yet working
        '''
        #t = new QTimer(this);
        connect(t, SIGNAL(timeout()), this, SLOT(perform()));
        #t->start(0);

    def snap(self,label='meas'):
        '''rememerbering last measurement'''
        if self.expername!="" :label=self.expername
        name=label+" %i"%(len(self.stack)+1)
        color=lcolor[len(self.stack)%len(lcolor)]
        print('adding to stack as '+name)
        self.stack[os.path.basename(name)]=self.graph.axes.plot(self.instr.pixtable,self.instr.last,
            label=name,color=color,lw=rc.line_width)[0]
        if hasattr(self,'stackwin'): 
            self.stackwin.additem(text=name)
            #list.insertItem(0,stack_prefix+' '+name)
            #self.stackwin.listBox.item(0).setTextColor(QtGui.QColor(scmap[color]))
        
        return
    def writeData(self,fname,data=None,xdata=None,labels=None,format='%f',type=None):
        '''saving data in ascii format'''
        if data==None: 
            if self.instr==None: 
                print('no instrument found')
                return
            data=self.instr.last
        comm=[]
        comm.append(" INTEG %.1f AVG %.1f SMOOTH %i\n"%(self.intime,self.mean,self.smooth))
        if self.termo: comm.append(" TEMP %.2f\n"%self.act_temp)
        if self.expernote: comm.apppend(self.expernote+'\n')
        if labels:
            header=' bin[%s]'%units[self.xunit]
            if type=='matlab': header+=rc.output_separ+'refer'
            header+=''.join([rc.output_separ+l for l in labels])
            comm.append(header+'\n')
        self.instr.write(fname,self.expername,comm=comm,data=data,xdata=xdata,type=type,format=format)

    def setSaveFileName(self,name=None,do_save=True):
        global save_dir
        import os
        if not os.path.exists(save_dir): os.makedirs(save_dir)
        if self.incounter:
            if name==None: name=str(self.saveFileNameLabel.text())
            name,ext=os.path.splitext(name)
            name=os.path.join(save_dir,name)
            if self.incounter>1: name=name[:name.rfind('_')]
            name+=('_%03i'%self.incounter)+ext
            self.incounter+=1
        elif name==None:
            self.fileName = self.fileDialog.getSaveFileName(self,self.tr("Saving data"),
                                             save_dir,self.tr("Matlab Files (*.dat);;Text Files (*.txt);;All Files (*)"))
            try:
                if self.fileName=='': return
            except:
                if self.fileName.isEmpty(): return
            name=str(self.fileName)#.toAscii())
        if name:
            text=os.path.basename(name)
            self.saveFileNameLabel.setText(text)
            if not self.incounter: save_dir=os.path.dirname(name)
            if not do_save: return name
            if save_mode==3:
                from numpy import array
                data=array([self.instr.flat,self.instr.last]).transpose()
            else: data=None
            if save_mode>2: type='matlab'
            else: type=None
            self.writeData(name,data=None,type=type)
            return name
            
            
    def saveLast(self,name=None):
        name=self.setSaveFileName(name,do_save=True)
        if (self.instr==None) or (name==None): 
            print('not saving')
            return
        if self.keep>0:
            color=lcolor[len(self.stack)%len(lcolor)]
            print('adding to stack as '+os.path.basename(name))
            self.stack[os.path.basename(name)]=self.graph.axes.plot(self.instr.pixtable,self.instr.last,
                label=name,color=color,lw=rc.line_width)[0]
        self.graph.draw()

    def saveRefer(self,name=None):    
        self.fileName = self.fileDialog.getSaveFileName(self,self.tr("Reference file"),
                                         self.saveFileNameLabel.text(),
                                         self.tr("Text Files (*.txt);;All Files (*)"))
        if not self.fileName.isEmpty():
            self.writeData(str(self.fileName.toAscii()),self.instr.flat)
    def saveExcel(self,name=None,head=1,calib=True):
        try:
            import pyexcel as xl
            wb=xl.Workbook()  
            sheet1=wb.add_sheet("data")
        except:
            alert(self,'XLS format not implemented')
            return            
        if name==None:
            self.fileName = self.fileDialog.getSaveFileName(self,self.tr("Saving data"),
                                     save_dir,self.tr("Excel Files (*.xls);;All Files (*)"))
            if self.fileName.isEmpty(): return
            name=str(self.fileName.toAscii())
        style=xl.StyleHelper()                #create a new style object
        style.SetBold("on") 
        if head>0:
            sheet1.write(head-1,0,'energy [eV]')
            text=self.expername
            if self.expernote: text+='\n'+self.expernote
            sheet1.write(head-1,1,text)
            if calib and self.instr.flat!=None: sheet1.write(head-1,2,'refer')
        nbin=len(self.instr.pixtable)
        for i in range(nbin):
            sheet1.write(i+head, 0, self.instr.pixtable[nbin-i-1])
            sheet1.write(i+head, 1, self.instr.last[nbin-i-1])
            if calib and self.instr.flat!=None: sheet1.write(i+head, 2, self.instr.flat[nbin-i-1])
            #sheet1.write(row,3,xl.Formula("A1*B2"))
        wb.save(name)

    def loadData(self,destin="stack",name=None,datacol=-1,type=None,min_bins=5,calib=None):    
        if name==None:
            self.refName = self.fileDialog.getOpenFileName(self,
                 self.tr("Load data"),save_dir,
                 self.tr("Text Files (*.txt);;Matlab style files (*.dat);;AvaSpec Export Files (*.ttt);;All Files (*)"))
            if self.refName.isEmpty(): return
            name=str(self.refName.toAscii())
        basecol=0
        import os
        #text=os.path.basename(name)
        #self.saveFileNameLabel.setText(text)
        asepar=rc.output_separ
        if os.path.exists(name): 
            if name[-3:]=='ttt': 
                skipped=[file.readline() for i in range(7)]
                asepar=';'
                print('treating as AvaSpec exported file')
                type='avasp'
            elif name[-3:]=='txt': 
                type='text'
                asepar=None
            elif name[-3:]=='dat': type='matlab'
            file=open(name)
            comm='#'
            print('data format used:'+type)
            if type=='matlab':
                comm='%'
                if datacol==-1: datacol=1
                basecol=0
            try:
                lines=[map(float,a.split(asepar)) for a in file.readlines() if a[0]!=comm]
            except:
                self.statusBar().showMessage("Wrong format -  aborted")
                return
            print('read %i lines of %i columns'%(len(lines),len(lines[0])))
            file.close()
            from numpy import array
            if datacol==None: # scanned line
                if (calib) or (len(lines)!=len(self.instr.pixtable)): 
                    p1=name.rfind("_")
                    p2=name.rfind(".")
                    if p1>0 and p2>p1:
                        cname=name[:p1+1]+"calib"+name[p2:]
                        if not os.path.exists(cname): 
                            self.statusBar().showMessage("Cannot find calibration "+cname)
                            return None,None
                        print('loading calibration from '+cname)
                        clines=[map(float,a.split(asepar)) for a in open(cname).readlines() if a[0]!=comm]
                        if len(clines)!=len(lines): return
                        self.instr.pixtable=array([c[0] for c in clines])
                        if len(c)>1: self.instr.flat=array([c[1] for c in clines])
                    else: 
                        self.statusBar().showMessage("Cannot guess reference file")
                        return None,None
                bins=self.instr.pixtable
                data=array(lines).transpose()
                print('loading %i spectra'%(len(data)))
                bname=os.path.basename(name)
                if bname.rfind(".")>0:
                    bname=bname[:bname.rfind(".")]
                bname+="_%02i"
                for i in range(len(data)):
                    vals=data[i]
                    color=lcolor[len(self.stack)%len(lcolor)]
                    self.stack[bname%i]=self.graph.axes.plot(bins,vals,color=color,lw=rc.line_width)[0]
                return bins,vals
            # normal mode
            bins,vals=array(lines).transpose()[[basecol,datacol]]
            if self.instr!=None:
                    sele=(bins<=self.instr.pixtable.max())*(bins>=self.instr.pixtable.min())
                    bins=bins[sele]
                    vals=vals[sele]
                    if sum(sele)<min_bins:
                        print('only %i bins within current range'%sum(sele))
                        return
            color=lcolor[len(self.stack)%len(lcolor)]
            if destin=="flat":
                if self.instr!=None:self.instr.flat=vals
                self.calib.axes.plot(bins,vals,color=color)
                self.calib.change_sizes()
                self.calib.draw()
                print('loaded reference')
            else:
                if destin=="stack":
                    self.stack[os.path.basename(name)]=self.graph.axes.plot(bins,vals,color=color,lw=rc.line_width)[0]
                    print('added to %i stack items'%len(self.stack))
                else:
                    if self.instr!=None:
                        if len(self.instr.pixtable)!=len(bins):
                            print('different base loaded')
                        else: self.instr.last=vals
                    self.graph.axes.plot(bins,vals,color=color)
                self.graph.draw()
            return bins, vals
    def loadRefer(self,type=None,name=None):
        self.loadData('flat',type=type,name=name)

    def loadMulti(self,type=None,name=None):
        self.loadData('stack',datacol=None,type=type,name=name,calib=True)
        if len(self.stack)>0: self.show_stack()
    def setExper(self):
        a,ok=QtGui.QInputDialog.getText(self, self.tr("experiment"),self.tr("name"))
        if ok: 
            self.expername=str(a)
            self.setWindowTitle("AvaSpec : "+self.expername)
        
    def saveNote(self):
        self.expernote=str(self.expEdit.text())
        return

#----------------------------------------------------------

    def getPeri(self):
        '''period estimation for a single measurement
        using either internal estimation
        '''
        if self.instr.last==None:
            alert(self,'Need measurement first',loud=3)
            return
        #global oxide_sel
        res=-1
        if self.graph.analbar!=None:
            erange=[self.instr.pixtable[i] for i in self.anrange]
            if erange[0]>erange[1]: erange=erange[::-1]
            if rc.thick_estimator=='humlicek': 
                res=sx.run_external([self.instr.last],self.instr.pixtable,ox=erange)[0]
        else:
            if rc.thick_estimator=='humlicek': 
                res=sx.run_external([self.instr.last],self.instr.pixtable,ox=rc.oxide_sel)[0]
            erange=rc.eranges[rc.oxide_sel]
        thck=0
        if res==-1:
            print('using internal thickness estimation')
            import scan_anal2
            scan_anal2.px=self.instr.pixtable#[imin:imax]
            scan_anal2.dt=[self.instr.last]
            scan_anal2.half_anal_width=rc.anal_poly_fit_width
            res=scan_anal2.anal(0,erng=erange,use_min=rc.anal_average_min_periods,amplim=rc.anal_pattern_min_ampl)
            from spectra import ev2um
            thck=1/(2*ev2um*abs(res[0]))
            res[1]*=thck/abs(res[0])
            if len(res)>4:res[4]*=thck/abs(res[0])
            res[0]=thck
            points=["mid.pos.    dist.[eV]  thick.[um]","-"*50]
            mper=1/(2*ev2um*abs(scan_anal2.pers))
            if rc.oxide_sel>=0:
                res[0]-=rc.ox_cor[rc.oxide_sel]
                mper-=rc.ox_cor[rc.oxide_sel]
            good_lim=mper.mean()*5
            ids=range(len(scan_anal2.pers))
            if scan_anal2.mids[0]>scan_anal2.mids[-1]: ids=ids[::-1]
            for i in ids:
                 if mper[i]>good_lim: continue
                 points.append("%-10.4f    %-10.4f  %-10.4f"%(scan_anal2.mids[i],abs(scan_anal2.pers[i]),mper[i]))
        if rc.show_period_estim:
            if thck>0:
                mpos=scan_anal2.mids
            #else:
                self.calib.axes.clear()
                self.calib.axes.plot(mpos,mper,'b+')
                xlim=self.calib.axes.get_xlim()
                self.calib.axes.axhline(res[0],xlim[0],xlim[1],color='r')
                self.calib.axes.figure.canvas.draw()
                self.graph.axes.set_xlim(xlim)
                self.graph.axes.figure.canvas.draw()
        text='Estimated thickness %.3f um'%res[0]
        if rc.anal_average_min_periods>0: text+=" [from %i min. vals]"%rc.anal_average_min_periods
        elif rc.anal_average_min_periods<0: text+=" [from %i max. vals]"%(-rc.anal_average_min_periods)
        extra=['(variab. %.3f um \n in range %.2f-%.2f eV)'%tuple([res[1]]+erange)]
        if len(res)>2: extra.append('trend (%.3f) corrected variab. %.3f'%list(res[-2:]))
        self.gui['result']=extend_alert(self,text,'Period analysis',extra,points)

#----------------------------------------------------------
    def zero_anal_data(self):
        from numpy import zeros
        from scanner.spectrax import anal_data
        xdim,ydim=self.control['nx'],self.control['ny']
        anal_data['mean']=zeros((xdim,ydim))#((2*self.control['nx']-1,2*self.control['ny']-1))
        anal_data['vari']=zeros((xdim,ydim))
        anal_data['peri']=zeros((xdim,ydim))
        anal_data['peri_disp']=zeros((xdim,ydim))
        anal_data['gpos']=[]

    def timeMeas(self,outname=None):
        '''timed measurements (instr.pulse)'''
        from spectrax import TimeDial
        if not hasattr(self,'timedial'):
            self.timedial=TimeDial(parent=self)
        ok=self.timedial.exec_()
        if ok: self.timedial.pulse()
            
    def scanSample(self,outname=None):
        '''scanning process using (linear and rotation) motors
        '''
        if self.instr==None: self.initSpec()
        #try:
        #    from scanner import spectrax as sx
        #except:
        #    import spectrax as sx
        if not hasattr(self,'scandial'): #
            if self.polar: self.scandial=spectrax.ScanDial(rc.rad_pts,rc.theta_pts,0,polar=True)
            else: self.scandial=spectrax.ScanDial(rc.x_pts,rc.y_pts,rc.sample_radius)
            self.scandial.parent=self
        self.control={}
        if outname: self.scandial.outname.setText(outname)
        else:
            ok=self.scandial.exec_()
            if not ok: 
                print('dialog exited unexpectedly (%s) - trying to continue'%str(ok))
                return
        rc.scan_start_analysis=self.scandial.ch_anal.checkState()>=1
        window=None
        runme=self.scandial.prepare_scan() # returns actual scanning procedure 
        if self.multitask:
            if 'scanwin' in self.gui:
                window=self.gui['scanwin']
                window.points=[]
            else:
                window = spectrax.ScanMap(polar=self.polar,parent=self)
                self.gui['scanwin']=window
                window.clear()
        self.zero_anal_data()
        if hasattr(self.instr.config,'Material') and self.instr.config.Material==b'simu':
            self.control['wait']=rc.simu_wait
        #
        if self.multitask:
            from threading import Event
            # why using timer??
            self.itimer = QtCore.QTimer(self)
            if rc.scan_start_analysis:
                self.control['anal_evt']=Event()
                QtCore.QObject.connect(self.itimer, QtCore.SIGNAL("timeout()"), self.analSpec)
                self.itimer.start(rc.anal_timer_period)
            self.control['queue']=self.report
            if window: 
                if rc.disp_anal_format: window.text_format=rc.disp_anal_format
                window.show()
            runme.start()
            if rc.scan_start_analysis: self.analSpec()
        else:
            if rc.scan_start_analysis: self.control['anal']=self.analSpec
            runme(self.control)
        # when scan finishes ############################
        if rc.scan_start_analysis:
            if self.scandial==None: 
                print('not saving')
                return
            cname=str(self.scandial.outname.text())
            if cname.find('.txt')>0: cname=cname.replace(".txt","calib.txt")
            elif cname.find('.')>0: cname=cname.replace(".","calib.")
            else: cname=cname.replace("_","_calib.txt")
            import os
            if not os.path.exists(cname): 
                print('cannot find calib. file '+cname)
                return
            if rc.thick_estimator=='munz':
                if self.anrange==None: excomm="anal.range %.2f-%.2f"%(self.instr.pixtable[self.erange[0]],self.instr.pixtable[self.erange[1]])
                else: excomm="anal.range %.2f-%.2f"%(self.instr.pixtable[self.anrange[0]],self.instr.pixtable[self.anrange[1]])
                sx.save_period(cname.replace("calib","period"),sx.anal_data['peri'],sx.anal_data['peri_disp'],comment=excomm)
            else:
                self.scanPeri(cname)
                print('finished thickness analysis using '+cname)
            if rc.show_thick_results: 
                print("showing file "+cname.replace("calib","period"))
                os.system(rc.show_text_editor+" "+cname.replace("calib","period"))
            
    def analSpec(self):
        '''analysing a single measurement
        during the scan'''
        import scan_anal2
        from numpy import median,abs
        #if 'anal_evt' in self.control: self.control['anal_evt'].clear()
        if self.instr.last==None: return # no last data measures
        if 'queue' in self.control and rc.use_queue:
            if self.report.qsize==0: return #measure not ready yet
            while self.report.qsize>0:
                mess=self.report.get()
                if 'scanwin' in self.gui: self.gui['scanwin'].statusBar().showMessage(mess)
                print()
                if mess=='finished': 
                    if hasattr(self,'itimer') and self.itimer!=None: 
                        self.itimer.stop()
                        self.itimer=None
                    return
                if mess.startswith('measured'): break
            else: 
                print('nothing measured')
                return
        from spectrax import anal_data
        if ('gpos' in anal_data) and ('gx' in self.control):
            anal_data['gpos'].append([self.control['gx'],self.control['gy']])
        if self.multitask: 
            self.graph.update_measure()
            ### merime to jeste jednou!!!?
            self.adjust_vert()
        if rc.oxide_sel>=0 and self.anrange==None: 
            from spectra import get_ids
            self.anrange=get_ids(self.instr.pixtable,rc.eranges[rc.oxide_sel])
        if self.anrange!=None: imin,imax=self.anrange
        else: imin,imax=None,None
        amean=self.instr.last[imin:imax].mean()
        scan_anal2.half_anal_width=rc.anal_poly_fit_width
        vari=scan_anal2.variab(self.instr.last[imin:imax])
        if not self.polar and 'radius' in self.control :
            i=self.control['nx']+self.control['x']-1
            j=self.control['ny']+self.control['y']-1
        else:
            i=self.control['x']#-1
            j=self.control['y']#-1
        #if 'swap' in self.control: i,j=j,i
        if i<0 or j<0: return -1
        if 'mean' in anal_data: #check sizes
            if i>=anal_data['mean'].shape[0]: return -2
            if j>=anal_data['mean'].shape[1]: return -3
            anal_data['mean'][i,j]=amean
        if 'vari' in anal_data: anal_data['vari'][i,j]=median(vari)
        print('measured: mean %.4f vari %.4f'%(amean,median(vari)))
        if 'peri' in anal_data:
            if scan_anal2.px==None or self.anrange!=None:scan_anal2.px=self.instr.pixtable[imin:imax]
            try:
                rep=scan_anal2.anal(self.instr.last[imin:imax],rep=4)
            except:
                print('period analysis failed for point %i/%i (range %i-%i)'%(i,j,imin,imax))
                return
            if abs(rep[0])>0:
                anal_data['peri'][i,j]=scan_anal2.conv_per/abs(rep[0])
                anal_data['peri_disp'][i,j]=rep[1]
            if 'peri' in self.control: # drawing in reference pane
                self.calib.axes.cla()
                self.calib.axes.plot(scan_anal2.mids,abs(scan_anal2.pers),'rs')
                self.calib.axes.set_xlim(self.graph.axes.get_xlim())
                self.calib.figure.canvas.draw()
        if 'scanwin' in self.gui: self.gui['scanwin'].add_point()
        if 'anal_evt' in self.control: self.control['anal_evt'].set()
    def scanPeri(self,cname=None):
        '''analysis of interference pattern of saved data using external program
        gets all data corresponding to given calibration file
        '''
        #global oxide_sel
        rc.oxide_sel=self.oxideModes.currentIndex()
        from spectrax import anal_data,anal_external
        if cname==None:
            self.scanName = self.fileDialog.getOpenFileName(self,
                     self.tr("Select calib data"),save_dir,
                     self.tr("Text Files (*.txt);;Matlab style files (*.dat)"))
            if self.scanName.isEmpty(): return
            cname=str(self.scanName)
        if cname.find('calib')>0:
            from glob import glob
            inlist=glob(cname.replace("calib","*"))
            inlist=[a for a in inlist if not a.find('calib.')>0]
            print('found %i data files'%len(inlist))
            save=cname.replace("calib","period")
        if len(inlist)==0: return
        rep=anal_external(inlist,cname,save=save,ox=rc.oxide_sel)
        anal_data['ext_peri']=rep[0]
        anal_data['ext_disp']=rep[1]
        #from matplotlib import figure

#-----------------------------------------------------------------

    def rconnect(self,saddr='147.251.57.145'):
        '''initial tests for remote connection 
        '''
        import conn
        conn.inverse=True;
        conn.tools={'spect':self.measure}
        try:
            self.statusBar().showMessage("Connecting to %s..."%saddr)
            conn.client(saddr)
        except:
            self.statusBar().showMessage("Connection failed...")
        else:
            self.statusBar().showMessage("Connection closed...")
        
    def upload(self,npoints=200):
        '''uploading data to database'''
        server='cuda.physics.muni.cz'
##-----------------------------------------------------------------------------------------------

    def initSpec(self,do_termo=False,show_comm=True,instrname=None,reset_instr=True):
        '''initialization process'''
        self.statusBar().showMessage("Testing spectrometer...")
        if reset_instr and self.instr!=None:
            self.instr.end()
            self.instr=None
        if instrname==None: instrname=rc.instrname
        if self.instr==None: 
        # here we create instrument instance
            if instrname=='avantes':
                if rc.scope_interface!=None:
                    self.instr=labin.linuscope(rc.scope_interface)
                    print('starting linux bypass')
                else:
                    self.instr=labin.avantes()
            elif instrname=='jaz':
                self.instr=labin.oceanjaz()
            elif instrname=='ocean':
                self.instr=labin.ocean()
            else:
                if hasattr(rc,'simu_bins'): nbin= rc.simu_bins
                else: nbin=None
                if hasattr(rc,'simu_range'): 
                    self.instr=labin.specscope(erange=rc.simu_range,nbin=nbin)
                    self.erange=rc.simu_range
                else: self.instr=labin.specscope(nbin=nbin)
        if self.instr==None: 
            alert(self,"spectrometer <b>not found</b>!")
            return
        self.instr.setup([self.emin,self.emax],integ=self.intime,aver=self.mean,smooth=self.smooth,unit=self.xunit)
        if self.instr.pixtable==None:
            alert(self,"cannot read pixel data! <b>Please restart spectroscope</b>!")
            return 
        self.rangeLabel.setText("%i-%i"%tuple(self.instr.pixrange))
        if not auto_init: alert(self,"spectrometer <b>OK</b>!")
        if do_termo:
            try:
                print("trying to init cryostat")
                self.termo=labin.linkam
                self.termoLabel.setText("temp. %.1f C",self.termo(close=False))
            except:
                self.termo=None
                print('failed')
        if show_comm: self.gui['toplay'].addLayout(self.gui['commlay']) # showing comments
        self.instr.parent=self
        self.referButton.setDisabled(0)
        self.measButton.setDisabled(0)
        self.runButton.setDisabled(0)
        self.darkButton.setDisabled(0)
        self.initButton.setDisabled(1)
        if rc.scan_extra_buttons and rc.init_ox_set: self.setWorkOxide()
        if instrname=='avantes':
            if self.polar:
                self.statusBar().showMessage("Testing motors...")
                self.rotInit()
            if rc.ard_port!=None: self.linInit()
            #rotational stage
            if rc.spec_unit:
                upos=units.index(rc.spec_unit)
                if upos>0:
                    aw.setUnit(upos) 
                    self.unitModes.setCurrentIndex(upos)

        self.statusBar().showMessage("Spectrometer %s OK!"%self.instr.name)
        
    def setDark(self):
        if self.instr!=None: self.instr.dark=self.instr.result(sub_dark=False,div_flat=False,smooth=rc.soft_smooth)
        else: return
        alert(self,"Dark <b>saved</b>!")
        self.referButton.setDisabled(0)
        self.referLabel.setText("dark %.2f"%self.instr.dark.mean())
        self.measure()
        
    def setReference(self,loud=0,vals=None):
        ''' obtaining reference (calibration) data
            checking the range of calibration measurement (should not have negative values)
        '''
        if self.instr==None:
            alert(self,"Must do <b>init</b> first!")
            return
        if vals!=None: self.instr.flat=vals
        else: self.instr.flat=self.instr.result(div_flat=False,smooth=rc.soft_smooth)
        if self.instr.flat.min()<=0: 
            nbin=sum(self.instr.flat<=0)
            text="%i reference values <b>not positive</b>!"%nbin
            print()
            #text+="<br/>please repeat calibration"
            #self.instr.flat=None
            from numpy import int,array,where
            # looking for the longest non-zero interval
            sele=(self.instr.flat>0).astype(int)
            dele=sele[1:]-sele[:-1]
            p1,p0=list(where(dele>0)[0]),list(where(dele<0)[0])
            p0.insert(0,0)
            p1.append(len(dele))
            
            if len(p0)<len(p1): p1=p1[:len(p0)]
            elif len(p0)>len(p1): p0=p0[-len(p1):]
            k=(array(p0)-array(p1)).argmax()
            prange=[p1[k]+1,p0[k]+1]
            print("widest nonzero range %i - %i"%tuple(prange))
            if prange[1]-prange[0]>10:
                conf=self.instr.config
                conf.m_StartPixel,conf.m_StopPixel=prange[0],prange[1]
                self.instr.pixtable=self.instr.pixtable[prange[0]:prange[1]]
                self.instr.flat=self.instr.flat[prange[0]:prange[1]]
                self.instr.dark=self.instr.dark[prange[0]:prange[1]]
                text+="<br/>reducing spectral range [%.2f - %.2f eV]"%(self.instr.pixtable[-1],self.instr.pixtable[0])
                alert(self,text)
            else: 
                alert(self,text+"no good calib. range found")
                self.instr.flat[self.instr.flat==0]=1. #avoid /0
        else: alert(self,self.tr("Reference <b>saved</b>!"))
        self.calib.axes.plot(self.instr.pixtable,self.instr.flat)
        text="reference avg. %.2f "%self.instr.flat.mean()
        if self.instr.dark!=None: text+="dark avg. %.2f "%self.instr.dark.mean()
        else: text+="dark: none"
        #self.calib.draw()
        #self.calib.axes.set_xlim(*self.graph.axes.get_xlim())
        self.calib.change_sizes()
        self.calib.axes.set_ylim(self.emin,self.emax)
        self.calib.draw()
        self.measure()
        self.graph.vscale=[0.8,1.2]
        #self.graph.axes.set_ylim(*self.graph.vscale)
        self.graph.draw()
        #self.graph.clear()
        self.referLabel.setText(text)
    
    # ------------- showing saved measurements

    def show_stack(self):
        if not 'stackwin' in self.gui:
            self.stackwin=StackWin()
            self.gui['stackwin']=self.stackwin
        self.stackwin.showlist(app=self)
    def incremental(self):
        if self.incname.checkState()>0: self.incounter=1
        else: self.incounter=0

    def set_keep(self):
        if self.grkeep.checkState()>0: self.keep=1
        else: self.keep=0
            
    # ----------- spectral manipulation 
    
    def useRefer(self):
        self.setReference(vals=self.instr.last)

    def repeat(self):
        d, ok = QtGui.QInputDialog.getInteger(self, self.tr("Repeat"),self.tr("nb. of repeated measurements"), self.smooth, 1, 50, 1)
        if ok: 
            rep=None
            for i in range(d):
                self.pbar.setValue(0.5)
                if rep!=None: rep+=self.measure()
                else: rep=self.measure()
                self.pbar.setValue(99.5)
            rep/=d
        self.instr.last=rep
        self.graph.update_measure(self.instr,vals=rep)

    def smoothing(self):
        #global smooth
        dial=QtGui.QInputDialog()
        d, ok = dial.getInteger(self, self.tr("smoothing"),self.tr("nb. of bins to smooth"), self.smooth, 0, 50, 1)
        if ok: 
            self.smooth=int(d)
            if self.instr!=None: self.instr.setup(None,integ=self.intime,aver=self.mean,smooth=self.smooth)
        
    def darkening(self):
        # internal dynamic darkening
        state=self.darkcheck.checkState()>0
        #mess=QtCore.QT_TR_NOOP("Setup <b>saved</b>!")
        #QtGui.QMessageBox.information(self, self.tr("Setup util."), mess)
        if self.instr!=None and hasattr(self.instr.config,'m_CorDynDark'):
            if state:
                val=self.instr.config.m_CorDynDark.m_ForgetPercentage
                d, ok = QtGui.QInputDialog.getInteger(self, self.tr("dynamic dark"),self.tr("remembered percentage"), val, 0, 100, 1)
                #self.instr.config.m_CorDynDark.m_ForgetPercentage=int(d)
            else: d=0
            self.instr.setup(None,integ=self.intime,aver=self.mean,dyndark=d,smooth=self.smooth,unit=self.xunit)

    def averaging(self):
        self.mean=self.avgcnt.value()
        self.intime=self.integbox.value()
        #mess=QtCore.QT_TR_NOOP("Setup <b>saved</b>!")
        #QtGui.QMessageBox.information(self, self.tr("Setup util."), mess)
        if self.instr!=None: self.instr.setup(None,integ=self.intime,aver=self.mean,smooth=self.smooth,unit=self.xunit)
    
    def measure(self,retrieve=True):
        if self.termo: 
            self.act_temp=self.termo(close=False)
            self.termoLabel.setText("temp. %.2f C"%self.act_temp)
        self.graph.update_measure(self.instr)
        #mess=QtCore.QT_TR_NOOP("Data <b>measured</b>!")
        #QtGui.QMessageBox.information(self, self.tr("Result"), mess)
        self.statusBar().showMessage("Data measured..", 2000)
        self.saveButton.setDisabled(0)
        if rc.auto_adjust_vertical: self.adjust_vert()
        if retrieve: return self.instr.last

    def variance(self,count=0):
        '''measures stability of the signal by comparison of repeated measurements'''
        if self.instr==None: return
        if count==0:
            a,ok=QtGui.QInputDialog.getText(self, self.tr("Variability"),self.tr("number of measurements"))
            if ok: count=int(str(a))
        from numpy import array,zeros,sqrt
        mes1=zeros(self.instr.last.shape)
        mes2=zeros(self.instr.last.shape)
        for i in range(count):
            self.graph.update_measure(self.instr)
            mes1+=self.instr.last
            mes2+=self.instr.last**2
        mes1/=count
        mes2/=count
        self.fluctdata=sqrt(abs(mes2-mes1*mes1))
        self.graph.axes.plot(self.instr.pixtable,self.fluctdata)
        self.graph.change_sizes()
        self.graph.figure.canvas.draw()
        #self.calib.change_sizes()
        return self.fluctdata
        
    def run(self):
        '''continuous measurement'''
        if self.running:
            self.running=False
            self.graph.itimer.stop()
            self.runButton.setText('Run')
            self.runButton.font().setBold(0)
        else:
            self.running=True
            if self.graph.itimer==None: self.graph.timer(self.deltime*1000)
            else: self.graph.itimer.start(self.deltime*1000)
            if self.graph.vscale==None: self.graph.vscale=self.graph.axes.get_ylim()
            self.graph.cnt=0
            self.runLabel.setText(str(self.graph.cnt))
            self.runButton.setText('Stop') #self.tr('<font style="color:red;font-weight:bold">Stop</b>'))
            self.runButton.font().setBold(1)
    
    # ---------------- external tools
    
    def showExtra(self):
        task=sx.AnalyseWindow()
        if self.instr!=None:
            task.data=self.instr.last
            if task.data==None:
                if len(self.stack)>0:
                    kk=self.stack.keys()[0]
                    task.base,task.data=array(self.stack[kk].get_xydata()).transpose()
                else:
                    print('no data available')
                    return
            else:
                task.base=self.instr.pixtable
            if self.anrange:
                task.base=task.base[self.anrange[0]:self.anrange[1]]
                task.data=task.data[self.anrange[0]:self.anrange[1]]
            if task.anal:
                task.anal.axes.cla()
                task.anal.axes.plot(task.base,task.data)
        task.parent=self
        self.analyse=task
        task.show()
        
    def showExtraMap(self):
        task=sx.AnalyseWindow(mode='scan')
        task.parent=self
        self.analyse=task
        task.show()
        
    def setAnalRange(self):
        '''setting currently selected graph zoom as range for analysis
        providing some basic characteristic in StatusBar (for peak analysis)
        '''
        bound=self.graph.hscale
        if bound==None: bound=self.graph.axes.get_xlim()

        from spectra import get_ids
        if self.instr==None:
            alert("x-axis not initialized yet")
            return
        imin,imax=get_ids(self.instr.pixtable,[bound[0],bound[1]])
        pmax=imin+self.instr.last[imin:imax].argmax()
        weights=self.instr.last[imin:imax]-self.instr.last[imin:imax].min()
        barycent=(weights*self.instr.pixtable[imin:imax]).sum()/weights.sum()
        self.rangeLabel.setText("anal.range: %.2f-%.2f eV"%tuple(bound))
        self.anrange=[imin,imax]
        self.statusBar().showMessage("max. value %.2f at %.3f eV [barycenter %.3f]"%(self.instr.last[pmax],self.instr.pixtable[pmax],barycent))
        if not self.showanal.isEnabled(): self.showanal.setEnabled(1)
        if self.showanal.isChecked(): self.analplot()
        
    def setWorkOxide(self):
        # show analysis range corresponding to given oxide
        if self.oxideModes.currentIndex()==0: return
        rc.oxide_sel=self.oxideModes.currentIndex()-1
        from spectra import get_ids
        self.anrange=get_ids(self.instr.pixtable,rc.eranges[rc.oxide_sel])
        if self.graph.analbar!=None:
            self.graph.analbar.remove()
            self.graph.analbar=None
            self.rangeLabel.setText("anal.range: %.2f-%.2f eV"%tuple(rc.eranges[rc.oxide_sel]))
        if not self.showanal.isEnabled(): self.showanal.setEnabled(1)
        if self.showanal.isChecked(): self.analplot()
        
    def getTemp(self):
        # if thermostatic module is plugged
        self.temp=self.termo(close=False)
        
    def scripted(self):
        #try: 
        #    import spectrax as sx
        #except:
        #    import scanner.spectrax as sx
        window = sx.ScriptWindow()
        window.resize(640, 480)
        window.createSample()
        self.script=window
        window.parent=self
        window.show()
        return

    def fileQuit(self):
        if self.instr!=None: 
            if labin.per==None: self.instr.end()
        if self.termo!=None: 
            try:
                self.termo(close=True)
            except:
                self.termo=None
        if hasattr(self.instr,"ard") and self.instr.ard!=None: self.instr.ard.close()
        if hasattr(self,"script"): self.script.close()
        if 'stackwin' in self.gui: self.gui['stackwin'].close()
        if 'scanwin' in self.gui:
            self.gui['scanwin'].close()
            del self.gui['scanwin']
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()

    def about(self):
        QtGui.QMessageBox.about(self, "About %s" % progname,
        u"""%(prog)s version %(version)s
        First attempt to communicate with AvaSpec-3648 machine
        Written with PyQT and matplotlib toolkit for Python
        Integrated support for Thorlabs rot. stage and Arduino mashup"""
        % {"prog": progname, "version": progversion})

    def linInit(self):
        if hasattr(self.instr,"ard") and self.instr.ard!=None: self.instr.ard.close()
        if rc.ard_port!=None: 
            self.instr.adrive(rc.motor_lin_rate)
            print("Arduino init OK")
    def rotInit(self):
        '''extra initialization for Thorlab servo rotator
        (scanner ver. 2)
        '''
        import serial
        try:
            import motor
            self.instr.motor=motor.MotorWindow(izero=rc.rot_zero_pos)
            self.instr.motor.run(0)
        except:
            print("servo init failed")
        rc.ystepsize=1.
        if self.control and 'gy' in self.control: self.control['gy']=0
        print('stage initialized')
        if self.instr.hnd: self.instr.device.AVS_SetDigOut(self.instr.hnd,2,1)
    def rotHome(self):
        self.instr.motor.run(0)
    def ardSpeed(self):
        if (not hasattr(self.instr,'ard')) or self.instr.ard==None: return
        d, ok = QtGui.QInputDialog.getInteger(self, self.tr("Speed"),self.tr("set time for 1 step in ms"), rc.motor_lin_rate, 1, 50, 1)
        e=int(d)
        if ok and e>0:
            self.instr.ard.write("PERI %i\n"%e)
            g=e//2
            if g<1: g=1
            self.instr.ard.write("GAP %i\n"%g)
            rc.motor_lin_rate=d
    def changeMulti(self):
        self.multitask=not self.multitask
    def changePolar(self):
        self.polar=not self.polar

stack_prefix="spect."

class StackWin(QtGui.QMainWindow):
    list=None
    #stack=None
    parent=None
    def __init__(self,stack=[]):
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("List of measurements")
        self.main_widget = QtGui.QWidget(self)
        self.layout = QtGui.QVBoxLayout(self.main_widget)
        self.topLabel = QtGui.QLabel(self.tr("Opened %i files:"%(len(stack))))
        self.layout.addWidget(self.topLabel)
        delButton = QtGui.QPushButton(self.tr("&Remove"))
        pltButton = QtGui.QPushButton(self.tr("&Plot"))
        savButton = QtGui.QPushButton(self.tr("&Save"))
        renButton = QtGui.QPushButton(self.tr("&Rename"))
        self.connect(delButton, QtCore.SIGNAL("clicked()"), self.delete)
        self.connect(pltButton, QtCore.SIGNAL("clicked()"), self.plotlist)
        self.connect(savButton, QtCore.SIGNAL("clicked()"), self.savelist)
        self.connect(renButton, QtCore.SIGNAL("clicked()"), self.rename)
        self.namEdit = QtGui.QLineEdit()
        hl = QtGui.QHBoxLayout()
        hl.addWidget(delButton)
        hl.addWidget(pltButton)
        hl.addWidget(savButton)
        self.layout.addLayout(hl)
        self.list = QtGui.QListWidget()
        self.list.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.layout.addWidget(self.list)
        hl = QtGui.QHBoxLayout()
        hl.addWidget(self.namEdit)
        hl.addWidget(renButton)
        self.layout.addLayout(hl)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        #self.stack=stack
    def rename(self):
        ilist=self.list.selectedItems()
        if len(ilist)==0: return
        ite=ilist[0]
        name=str(ite.text())[len(stack_prefix)+1:]
        newname=str(self.namEdit.text())
        ite.setText(stack_prefix+' '+newname)
        self.parent.stack[newname]=self.parent.stack[name]
        del self.parent.stack[name]
    def delete(self):
        ilist=self.list.selectedItems()
        name=str(ilist[0].text())[len(stack_prefix)+1:]
        print('want to remove "%s"'%name)
        for il in ilist:
            il.setHidden(1)
        self.list.removeItemWidget(ilist[0])
        if name in self.parent.stack:
            self.parent.stack[name].remove()
            del self.parent.stack[name]
    def savelist(self):
        if self.parent.stack==None: return
        ilist=[str(li.text())[len(stack_prefix)+1:] for li in self.list.selectedItems()]
        if len(ilist)==0:
            alert(self,'Select spectra to save','Message',loud=1)
            return
        slist=[]
        for s in ilist:
            if not s in self.parent.stack: 
                s=s.split()[-1]
                if not s in self.parent.stack: continue
            slist.append(s)
        name=self.parent.setSaveFileName(do_save=False)
        from numpy import array
        odata=array([self.parent.stack[s].get_ydata() for s in slist]).transpose()
        self.parent.writeData(name,odata,xdata=self.parent.stack[s].get_xdata(),labels=slist)
    def plotlist(self):
        if self.list==None or self.parent.stack==None: return
        ilist=[str(li.text())[len(stack_prefix)+1:] for li in self.list.selectedItems()]
        if len(self.parent.stack.keys())==0: return
        for k in self.parent.stack.keys():
            self.parent.stack[k].set_visible(k in ilist)
        self.parent.stack[k].figure.canvas.draw()
        if self.parent.instr and len(ilist)>0: self.parent.instr.last=self.parent.stack[ilist[-1]].get_ydata()
        #if self.app!=None: self.parent.graph.draw()
    def additem(self,id=-1,text=None):
        if text!=None:
            if id<0: id=self.list.count()
            self.list.insertItem(id,stack_prefix+' '+text)
        ite=self.list.item(id)
        if self.parent==None: return
        a=str(ite.text())[len(stack_prefix)+1:]
        if not a in self.parent.stack: a=a.split()[-1]
        if not a in self.parent.stack: return
        color=self.parent.stack[a].get_color()
        if color in scmap: color=scmap[color]
        ite.setTextColor(QtGui.QColor(color))
        #self.parent.connect(ite, QtCore.SIGNAL("itemDoubleClicked()"), self.list.openPersistentEditor(ite))
        self.topLabel.setText(self.tr("%i spect. in stack"%(len(self.parent.stack))))

    def showlist(self,type='list',app=None):
        if app==None:
            self.show()
            return
        self.parent=app
        stacklist = QStringList()
        for a in sorted(app.stack.keys()):
            stacklist.append(stack_prefix+" "+a)
        self.list.insertItems(0, stacklist)
        for i in range(len(app.stack)):
            self.additem(id=i)
        #self.list=app.listBox
        #self.stack=app.stack
        self.topLabel.setText(self.tr("%i spect. in stack"%(len(self.parent.stack))))
        self.show()
        return self
    def closeEvent(self, ce):
        if self.parent!=None and ('stackwin' in self.parent.gui):
            del self.parent.gui['stackwin']
        del self.parent.stackwin


global qApp
qApp=None
def run(range=None,init=False,refname=None,loadname=[],extra=None):
    global aw,qApp
    if qApp==None: 
        qApp = QtGui.QApplication(['myapp'])
    aw = ApplicationWindow()
    if range!=None: aw.emin,aw.emax=tuple(range)
    aw.setWindowTitle("PyAvaSpec")
    aw.show() 
    if init: 
        if extra=='sim': aw.initSpec(instrname='simula')
        else: aw.initSpec(instrname=rc.instrname)
    if refname: aw.loadRefer('real',refname)
    if len(loadname)>0:
        aw.graph.axes.cla()
        aw.graph.change_sizes()
        aw.graph.axes.hold(1)
        for l in loadname:
            xdat=aw.loadData(name=l,type='real')[0]
        aw.graph.axes.set_xlim(xdat[0],xdat[-1])
    if extra=='sim':
        aw.instr.config.Material=bytes(rc.simu_calib,'ascii')
        aw.setReference()
        aw.instr.config.Material=b'simu'
        #ite=aw.stack.keys()[0]
        #aw.graph.axes.set_xlim(aw.stack[ite].get_xdata().min(),aw.stack[ite].get_xdata().max())

if __name__=="__main__":
    qApp = QtGui.QApplication(sys.argv)
    #pylab_setup()
    aw = ApplicationWindow()
    aw.setWindowTitle("PyAvaSpec")
    if '--simu' in sys.argv: aw.initSpec(instrname='simula')
    elif '--init' in sys.argv: aw.initSpec(instrname=rc.instrname)
    aw.show()
    sys.exit(qApp.exec_())
    #qApp.exec_()
