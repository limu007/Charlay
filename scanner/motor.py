import win32com.client
ProgId="MGMOTOR.MGMotorCtrl.1"

from PyQt4 import QAxContainer,QtCore,QtGui

#import labrc as rc

class MotorWindow(QtGui.QMainWindow):
    instr=None
    stack={}
    ok=False
    actpos=0
    zeropos=0
    scale=13.33587
    def __init__(self,izero=None,sernum=83825897):
        QtGui.QMainWindow.__init__(self)

        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Motor window")

        self.main_widget = QtGui.QWidget(self)
        self.runButton = QtGui.QPushButton(self.tr("Set ran&ge"))
        self.connect(self.runButton, QtCore.SIGNAL("clicked()"), self.scan)
        hl = QtGui.QHBoxLayout()
        hl.addWidget(self.runButton)
        from PyQt4 import QAxContainer
        self.mowid=QAxContainer.QAxWidget(ProgId)
        hl.addWidget(self.mowid)
        try:
            self.mowid.setProperty('HWSerialNum',QtCore.QVariant(sernum))
        except:
            self.mowid.setProperty('HWSerialNum',sernum)
        l = QtGui.QVBoxLayout(self.main_widget)
        l.addLayout(hl)
        if izero!=None: self.zeropos=izero
        return

    def go(self,pos=5.):
        if not self.ok: 
            self.mowid.dynamicCall("StartCtrl")
            self.ok=True
        try:
            self.mowid.dynamicCall("SetAbsMovePos(int,float)",QtCore.QVariant(1),QtCore.QVariant((pos+self.zeropos)/self.scale))
            self.mowid.dynamicCall("MoveAbsolute(int,bool)",QtCore.QVariant(1),QtCore.QVariant(False))
        except:
            self.mowid.dynamicCall("SetAbsMovePos(int,float)",1,(pos+self.zeropos)/self.scale)
            self.mowid.dynamicCall("MoveAbsolute(int,bool)",1,False)
        #self.mowid.dynamicCall("MoveAbsolute",1,False)

    def relat(self,pos=5.):
        if abs(pos)<0.5: return
        if abs(pos)>500: return
        self.actpos+=pos
        if self.actpos<0: self.actpos=0 #+=360
        if self.actpos>600: self.actpos=600
        self.go(self.actpos)
        
    def scan(self,slum=10.,npt=10):
        from time import sleep
        for i in range(npt):
            self.relat(10)
            sleep(slum)

def old_init():
    '''
    streamObj = qtaxcontainer.QAxObject("ADODB.Stream")
    streamObj.dynamicCall( "Open()" );
    streamObj.dynamicCall( "WriteText(QString&)", qt.QVariant(aStringWhichContinsXmlInRightFormat) ); 
    streamObj.dynamicCall( "SetPosition(int)", qt.QVariant(0) );

    nebo import pywin
    '''
    s=QAxContainer.QAxObject(ProgId)
    s.setProperty('HWSerialNum',QtCore.QVariant(83825897))
    s.dynamicCall('SetAbsMovePos',1,50/13.33587)
    s.dynamicCall('MoveAbsolute',1,False)
    return s

global qApp,aw
qApp=None
aw=None

def run(sernum=83825897):
    global qApp,aw
    qApp = QtGui.QApplication(['myapp'])
    aw = MotorWindow(sernum)
    aw.show()
