# -*- coding: utf-8 -*-

smoothing=0
mean=20
intime=25
use_dyndark=0
soft_smooth=9
instrname='jaz'

debug=0

spec_mode='transmissivity'#'reflectivity'
spec_unit='eV'
meas_range=1.2,4.9 # in electronVolts

hard_spec_min=-0.1
hard_spec_max=20

lang='cs' #language of the interface
#inst_dir='/home/optik/code/python/scanner/'#installation directory
#inst_dir='C:\\Python26\\Lib\\site-packages\\scanner\\'
inst_dir='/home/limu/Code/Monty/scanner/'
date_dir='/tmp/data'
auto_adjust_vertical=True # after measurement

#Avantes measure delay in seconds
sleepfact=.0012
sleepcon=.22

spectr_nb_trials=3

#----------------------------
polar_stage=False
scan_axes_order=0 # for rot. stage 0= radial first, 1=rotate first
                  # for lin. stage 0= upper first, 1=lower first
base_xdir=1 # -1 for xy-scanner, +1 for rotation stage
base_ydir=1

comm_ard_wait=0.05 # delay in sec between subsequent commands to Arduino 
simu_wait=1.1
ard_port=None
ard_speed=38400
ard_mode=2 #1-single drive, 2-double
#ard_port=None #'COM3'#'/dev/ttyUSB0'

# rotation stage
motor_rot_wait=.9 #sec
motor_rot_speed=0.15 #sec/degree

# initial step sizes
theta_step=60. #degrees
rad_step=10. #mm
cart_step=10. #mm in cartesian XY grid

#default scanning field
x_pts=3
y_pts=4
theta_pts=5
rad_pts=2
single_central_pt=False
scan_home_ret=True
sample_radius=0

rot_zero_pos=124
rot_linear_zero=7100 #steps from center to right edge for rotational setup

# stepper linear motor
motor_lin_wait=0.1
motor_lin_speed=0.0005  #sec/step
motor_lin_rate=1 #ms/step

#motors
mm4step=0.0031463
xstepsize=mm4step # in mm
ystepsize=mm4step # in mm

#--------------------------------------------------
# reference materials

ref_tables={'Si':'/home/optik/code/python/data/TblRefSi3648.dat',
'Si-eps':'/home/optik/code/python/data/si.dat'}
ref_tables['Si']='C:\Documents and Settings\Kremik\Dokumenty\Python\period\TblRefSi3648.dat'
ref_tables['Si-eps']='C:\Documents and Settings\Kremik\Dokumenty\Python\period\si.dat'
ref_sample=3 # taking every ref_sample point from tables above

#--------------------------------------------------
ooihome='/opt/OmniDriver/OOI_HOME'
java_jdk='/usr/local/lib/jdk'
java_home='/opt/OmniDriver/OOI_HOME/_jvm'
java_jvm='/usr/local/lib64/libjvm.so'
java_jvm='/opt/OmniDriver/OOI_HOME/_jvm/lib/amd64/server/libjvm.so'
cool_temp=0

#--------------------------------------------------
scope_interface='/home/limu/Code/Utils/Cython/face/usbtest' #Linux workaround
linux_write_delay=0.6 #waiting for all data to arrive

eps_ext={'Si':[-2.61678, [[179.4187,5.41916,0.001],[  16.3506,3.17564,0.12138],   [138.945163, 4.60986142,0.001]]]} # dielect. function models

mod_ext={'Si':[0.01401,0,0.0423,0,3.51906]} # Si model by JH : refr. index
mod_ext['SiO2']=[1.68683e-05,0,1.961953e-03,0,1.4493423] #: 

eps_trans={'Si':"cSi_Asp",'SiO2':"SiO2_gl"}
# spectral calibration
spec_pixtable_source=1 # 0-uses stored polynom 1-uses stored pixel position
dyndark_forget=50 #percentage for dynamic dark
#--------------------------------------------------
max_aver_count = 100

scan_extra_buttons=True
scan_extra_menu_experim=True

# analysis + multitasking
use_queue=True
anal_timer_period=500 #used in itimer
#synchronization uses queue - if measurement is not ready, queue is empty 

#SOI specific parameters (see spectrax.anal_external)
oxide_sel=-1

ox_wid=[1.,0.6,0.3]
ox_cor=[0.12,0.063,0.03] # thickness correction for 1 um and 0.6mm respectively
eranges=[[1.33,1.66],[1.4,2.08],[1.46,2.28]] # ranges for period extraction in electonvolts
init_ox_set=False

scan_realtime=False # multitasking - analysis in real time
scan_start_analysis=False
thick_estimator='munz' #'munz' or 'humlicek'
thick_range=[] #[1.,3.]
anal_use_median=False
anal_average_min_periods=3
anal_poly_fit_width=5
anal_pattern_min_ampl=None

show_thick_results=False
show_text_editor="notepad"

output_separ='\t'

# simulating physical surface
simu_range=[1.1,2.8]
simu_mater="InSb"
#simu_mater=["Si","SiO2","Si"]
#simu_layer=[3800,1000]
#simu_layer_fluc=[0.01,0.01]
simu_calib="Si"

#-------------------------------------------
# interface

show_default_range=0 #show analytic. range
show_period_estim=True
disp_anal_format="%.3f"
disp_anal_size=10
disp_mark_size=10
disp_anal_shift=0.06

set_maximized=False #full screen from start

interact=0 # 1=alerts as extra boxes

line_color='r'
line_width=0.5

vert_exten=0.2

