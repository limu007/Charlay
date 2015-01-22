# -*- coding: utf-8 -*-

from PyQt4 import QtGui, QtCore
try:
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar 
except:
    from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
    from matplotlib.backends.backend_gtk import NavigationToolbar2GTK as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector 

from numpy import arange
#from . 
import labrc as rc

fluid_lims=False #can change vertical limits interactively?

def alert(self,text,title="Calib. util.",loud=0):
    if loud>=(1-rc.interact):
        mess=QtCore.QT_TR_NOOP(text)
        QtGui.QMessageBox.information(self, self.tr(title), mess)
    else:
        print(text)

def extend_alert(self,text,title="Results",extras=[],details=[]):
    box=QtGui.QMessageBox(QtGui.QMessageBox.Information,title,text)
    #title=self.tr(title)
    if len(extras)>0:box.setInformativeText("\n".join(extras))
    if len(details)>0:box.setDetailedText("\n".join(details))
    box.exec_()
    return box

#--------------------------------------------

class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    xlabel='energy [eV]'
    ylabel='callibration'
    def compute_initial_figure(self):
        return
    def __init__(self, parent=None, width=6, height=4, dpi=100, polar=False):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111,polar=polar)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)
        self.compute_initial_figure()
        FigureCanvas.__init__(self, fig)
        if True:
        #try:
            self.setParent(parent)
            #FigureCanvas.setSizePolicy(self,
            #        QtGui.QSizePolicy.Expanding,
            #        QtGui.QSizePolicy.Expanding)
            FigureCanvas.updateGeometry(self)
        #except:
        #    print("simple graphics - beware")
        self.figure.subplots_adjust(left=0.08,right=0.96)

    def clear(self,reset=True):
        self.figure.clf()
        if reset:
            self.axes=self.figure.gca()
            self.draw()
    def change_sizes(self,labsize=8):
        '''modify font sizes'''
        self.axes.xaxis.label.set_fontsize(labsize)
        for lab in self.axes.xaxis.get_ticklabels():
            lab.set_fontsize(labsize)
        self.axes.yaxis.label.set_fontsize(labsize)
        for lab in self.axes.yaxis.get_ticklabels():
            lab.set_fontsize(labsize)

class MyStaticMplCanvas(MyMplCanvas):
        """Simple canvas"""
        init=False
        def __init__(self, *args, **kwargs):
            MyMplCanvas.__init__(self, *args, **kwargs)
            if 'xlabel' in kwargs: self.xlabel=kwargs['xlabel']
            if 'ylabel' in kwargs: self.ylabel=kwargs['ylabel']
        def compute_initial_figure(self):
            self.axes.set_xlabel(self.xlabel)
            self.axes.set_ylabel(self.ylabel)
            #self.axes.grid()
            if self.init: self.draw()
            self.axes.grid(1)
            self.init=True

def pylab_setup():
    from pylab import rcParams
    params = {'backend': 'qt',
              'axes.labelsize': 8,
              'text.fontsize': 8,
              'legend.fontsize': 10,
              'xtick.labelsize': 8,
              'ytick.labelsize': 8}
    rcParams.update(params)

class MyDynamicMplCanvas(MyMplCanvas):
    """A canvas that updates itself with a new plot."""
    itimer=None
    parent=None
    vscale=None
    hscale=None
    localYMax = 0 
    last=None
    cnt=0
    analbar=None
    
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        self.mpl_connect('button_press_event', self.onclick) 
        self.ylabel=rc.spec_mode
        
    def timer(self,delay=500):
        self.itimer = QtCore.QTimer(self)
        QtCore.QObject.connect(self.itimer, QtCore.SIGNAL("timeout()"), self.update_measure)
        self.itimer.start(delay)

    def compute_initial_figure(self):
        self.axes.plot([0, 1, 2, 3], [1, 2, 0, 4], 'r')
        self.axes.set_xlabel(self.xlabel)
        self.axes.set_ylabel(self.ylabel)

    def clear(self,reset_lims=True):
        self.axes.cla()
        self.axes.grid(1)
        self.last=None
        if reset_lims:
            self.hscale=None
            self.vscale=None
        self.change_sizes()
        self.analbar=None
        self.draw()
        
    def update_measure(self,instr=None,mar=0,vals=None):
        ''' do actual measurement and plot it
        '''
        # mar: margin not drawn
        if instr==None and self.parent!=None: instr=self.parent.instr
        from time import sleep
        if instr==None:
            sleep(.5)
            return
        self.cnt+=1
        self.axes.hold(1)
        if self.last==None: 
            self.axes.cla()
            self.axes.grid(1)
            self.axes.set_ylabel(self.ylabel)
            self.axes.set_xlabel(self.xlabel)
        else: 
            if fluid_lims:
                self.hscale=self.axes.get_xlim()
                self.vscale=self.axes.get_ylim()
            if len(self.axes.lines)>0: self.last.remove()
        self.parent.runLabel.setText("<b>%i</b>"%(self.cnt))
        from numpy import arange
        # here we measure data
        if vals==None: yvals=instr.result(smooth=rc.soft_smooth)
        else: yvals=vals
        if self.parent!=None and self.parent.table!=None: 
            if len(self.parent.table)==len(yvals): yvals*=self.parent.table[:,1]
            else: 
                print('calibration table of wrong size')
                self.parent.table==None
        xvals=instr.pixtable
        if yvals==None: return
        if xvals==None: xvals=arange(yvals.shape)
        lar=min(len(xvals),len(yvals))-mar
        self.last=self.axes.plot(xvals[mar:lar],yvals[mar:lar], rc.line_color,lw=rc.line_width)[0]
        if self.hscale!=None: self.axes.set_xlim(*self.hscale)
        if self.vscale!=None: self.axes.set_ylim(*self.vscale)
        self.change_sizes()
        self.axes.grid(1)
        #MyMplCanvas.change_sizes(self)
        self.draw()
    def onclick(self,event):
        if event.ydata != None:
            self.localYMax = int(event.ydata) 

class MyNavigationToolbar(NavigationToolbar) :
    span=None
    hZoom=True
    def __init__(self , canvas , parent , direction = 'h' ) :
        #NavigationToolbar.__init__(self,parent,canvas)
        #self.layout = QVBoxLayout( self )
        self.canvas = canvas
        QtGui.QWidget.__init__( self, parent )

        if direction=='h' :    self.layout = QtGui.QHBoxLayout( self )
        else : self.layout = QtGui.QVBoxLayout( self )
        self.layout.setMargin( 2 )
        self.layout.setSpacing( 0 )
        NavigationToolbar.__init__( self, canvas, parent ) 
    def zoomToggle(self):
        #self.toolbar.zoom() #this implements the classic zoom
        if self.hZoom:
            self.hZoom = False
            if self.span: self.span.visible = False
        else:
            self.hZoom = True
            if self.span: self.span.visible = True

class GraphWindow(QtGui.QMainWindow):
    '''window containing 2 canvases
    '''
    instr=None
    stack={}
    erange=[1,3647] #depends later on spectroscope connected
    emin,emax=1.13,2.8
    xunit=0
    rebins=None
    anrange=None

    def __init__(self):
        QtGui.QMainWindow.__init__(self)

        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("application main window")
        if hasattr(rc,'meas_range'): self.emin,self.emax=rc.meas_range #only in eV?
        self.main_widget = QtGui.QWidget(self)
        self.graph = MyDynamicMplCanvas(self.main_widget, width=8, height=6, dpi=100)
        self.calib = MyStaticMplCanvas(self.main_widget, width=8, height=4, dpi=100)
        self.graph.parent=self # link back for communication
        self.graph.change_sizes()
        self.calib.change_sizes()
        return
        
    def grange(self):
        '''plot range'''
        if self.randial==None: 
            from spectrax import RangeDial
            self.randial=RangeDial(self.graph.axes)
        else:
            self.randial.read_axes(self.graph.axes)
        self.randial.exec_()
        self.graph.hscale=[self.randial.xmin.value(),self.randial.xmax.value()]
        self.graph.vscale=[self.randial.ymin.value(),self.randial.ymax.value()]
        self.graph.axes.set_xlim(*self.graph.hscale)
        self.calib.axes.set_xlim(*self.graph.hscale)
        self.graph.axes.set_ylim(*self.graph.vscale)
        self.graph.draw()
        self.calib.draw()
        
    def onselect(self, xmin, xmax):
        if self.toolbar.hZoom:
            #self.graph.axes.set_ylim(ymax = self.localYMax)
            self.graph.axes.set_xlim(xmin,xmax)
            if self.instr!=None and self.instr.flat!=None: 
                self.calib.axes.set_xlim(xmin,xmax)
                self.calib.axes.grid(1)
                self.calib.draw()
            self.graph.draw()
        self.graph.hscale=(xmin,xmax) 
        #self.graph.vscale=(xmin,xmax) 
        
    def set_binning(self):
        d, ok = QtGui.QInputDialog.getInteger(self, self.tr("binning"),self.tr("nb. of bins"), 20, 0, 1000, 10)
        if ok: 
            nbins=int(d)
            if d<0:
                self.rebins=None
                return
            step=(self.instr.pixtable[-1]-self.instr.pixtable[0])/nbins
            self.rebins=((self.instr.pixtable-self.instr.pixtable[0])/step).astype('int')
            self.rebins[-1]=self.rebins[-2]
    #ndimage.mean(y,spectrac.aw.rebins,range(19))

    def get_vscale(self,imin=None,imax=None):
        '''min/max value of all displayed graphs
        '''
        data=[f.get_ydata()[imin:imax] for f in list(self.stack.values()) if f.get_visible()]
        if len(data)==0 and self.instr.last!=None: data.append(self.instr.last[imin:imax])
        if len(data)==0:
            alert(self,"No data plotted")
            return 0,0
        #if self.graph.last: data.append(self.graph.last.get_ydata()[imin:imax])
        return [min([d.min() for d in data]),max([d.max() for d in data])]

    def reset(self,exten=0.2):
        if self.stack!=None and len(self.stack)>0: 
            self.graph.hscale=list(self.stack.values())[0].axes.dataLim.get_points()[:,0]
        else: 
            if self.xunit==1: self.graph.hscale=[self.instr.pixtable[0],self.instr.pixtable[-1]]
            else: self.graph.hscale=[self.instr.pixtable[-1],self.instr.pixtable[0]]
        self.graph.vscale=self.get_vscale()
        vspan=(self.graph.vscale[1]-self.graph.vscale[0])*exten
        self.graph.vscale=[self.graph.vscale[0]-vspan,self.graph.vscale[1]+vspan]
        self.graph.axes.set_xlim(*self.graph.hscale)
        self.calib.axes.set_xlim(*self.graph.hscale)
        self.graph.axes.set_ylim(*self.graph.vscale)
        self.graph.axes.grid(1)
        self.calib.axes.grid(1)
        self.graph.draw()
        self.calib.draw()
        
    def adjust_vert(self,exten=0,use_stack=True):
        ''' should adjust according to all graphs plotted
        '''
        if exten==0: exten=rc.vert_exten
        from spectra import get_ids
        if self.instr: 
            imin,imax=get_ids(self.instr.pixtable,self.graph.axes.get_xlim())
            #print "using range %i / %i"%(imin,imax)
        elif len(self.stack)>0:
            ite=list(self.stack.keys())[0]
            imin,imax=get_ids(self.stack[ite].get_xdata(),self.graph.axes.get_xlim())
        else: imin,imax=None,None
        self.graph.vscale=self.get_vscale(imin,imax)
        xlim=self.graph.axes.get_xlim()
        
        data=[f.get_ydata()[imin:imax] for f in list(self.stack.values())]
        if self.instr and self.instr.last!=None: data.append(self.instr.last[imin:imax])
        elif self.graph.last: data.append(self.graph.last.get_ydata()[imin:imax])
        if len(data)==0: return 
        self.graph.vscale=[min([d.min() for d in data]),max([d.max() for d in data])]
        vspan=(self.graph.vscale[1]-self.graph.vscale[0])*exten
        self.graph.vscale=[self.graph.vscale[0]-vspan,self.graph.vscale[1]+vspan]
        self.graph.axes.set_ylim(self.graph.vscale[0],self.graph.vscale[1])
        self.graph.axes.set_xlim(min(xlim),max(xlim))
        self.graph.draw()
        #self.graph.axes.grid(1)

    def analplot(self):
        '''show analytic range
        '''
        if self.anrange==None: return
        if self.graph.analbar:
        #if self.showanal.checkState()>0:
            self.showanal.setChecked(0)
            for a in self.graph.axes.patches:
                    self.graph.axes.patches.remove(a)
            #self.graph.analbar.remove()
            self.graph.analbar=None
        else:
            self.showanal.setChecked(1)
            pos=[self.instr.pixtable[self.anrange[0]],self.instr.pixtable[self.anrange[1]]]
            self.graph.analbar=self.graph.axes.axvspan(pos[0],pos[1],facecolor='g',alpha=0.2)
        self.graph.draw()
        if self.graph.hscale: self.graph.axes.set_xlim(*self.graph.hscale)
            
    def clear(self):
        self.graph.clear()
        self.stack.clear()


def checkItem(text,confunc,context,init=False):
    '''creating a checkbox with annotation
    '''
    ch_wid=QtGui.QCheckBox(context.tr(text))
    context.connect(ch_wid, QtCore.SIGNAL("clicked()"), confunc)
    mee_box=QtGui.QWidgetAction(context)
    mee_box.setDefaultWidget(ch_wid)
    ch_wid.setChecked(init)
    return mee_box

