# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 18:28:18 2013
 
@author: ramazan tekin
 
Kortikal hücreler
ref : Bazhenov M, Timofeev I, Steriade M, Sejnowski TJ (2002) Model of 
thalamocortical slow-wave sleep oscillations and transitions to activated 
States. J Neurosci 22:8691–8704. 
 
Bazhenov M, Timofeev I, Steriade M, Sejnowski TJ (1998b) Computational models 
of thalamocortical augmenting responses. J Neurosci 18:6444–6465
 
PY ve RE hücre modeli 
 
 
 
---------------------------------------------------------------
Talamik hücreler
ref : 
Destexhe,ve ark. 1996, Ionic mechanisms underlying synchronized oscillations 
and propagating waves in a model of ferret thalamic slices.
 
Destexhe, ve ark. 1994a, A model of spindle rhythmicity in the isolated thalamic reticular
 
Destexhe, ve ark. 1996b, In vivo, in vitro and computational analysis of 
dendritic calcium currents in thalamic reticular neurons.  
  
Bazhenov ve ark. 1998, Cellular and network models for intrathalamic 
augmenting responses during 10 Hz stimulation. 
 
Bazhenov ve ark. 2002, Model of Thalamocortical Slow-Wave Sleep Oscillations 
and Transitions to Activated States
 
Destexhe ve ark. 1993, Ionic Mechanisms for Intrinsic Slow Oscillations in 
Thalamic Relay Neurons. Biophysical Journal
 
TC hücre modeli 
 
fast sodyum INa ve fast potasyum IK (Traub and Miles, 1991), 
low-threshold Ca+2 IT (Huguenard & McCormick; Huguenard & Prince (1992)), 
potassium leak akım I_KL (Bazhenov ve ark. 2002),
Hiperpolarizasyon aktivasyonlu katyonik Ih (Harigawa & Irisawa,1989;McCormick & Pape, 1990;
Huguenard & McCormick 1992;Destexhe ve ark. 1993; Destexhe ve ark., 1996) 
akımları içermektedir (Destexhe ve ark., 1996;Bazhenov & ark. 2002).
 
RE hücre modeli 
 
fast sodyum INa ve fast potasyum IK (Traub & Miles, 1991), 
low-threshold Ca+2 IT (Huguenard & McCormick 1992;Huguenard & Prince 1992), 
potassium leak akım I_KL (Bazhenov ve ark. 2002),
akımları içermektedir (Bazhenov ve ark. 2002).
 
"""

from brian import *
import scipy.io as sio
import os
import random as rnd
from matplotlib.font_manager import fontManager, FontProperties
from tezmodul import *

# Genel Parametreler   
os.chdir('/home/rtekin/Desktop/scriptler/tez_scrpts/Bazhenov2002/')
#os.chdir('/media/OS/ramazan/scriptler/tez_scrpts/Bazhenov2002/')

path_data = './deneme_ubuntu/tez_data/TalamoKortikalNet/data/'
path_images = './deneme_ubuntu/tez_data/TalamoKortikalNet/images/'

if not os.path.exists(path_data):
    os.makedirs(path_data)
    print "data klasörü oluşturuldu..."
if not os.path.exists(path_images):
    os.makedirs(path_images)
    print "images klasörü oluşturuldu..."

defaultclock.dt=0.1*ms
dt = defaultclock.dt
myclock = Clock(dt=.1 * ms)
clock_rec = Clock(dt= .1 * ms)

duration = 1*1000 * ms

new = False # True # 
dispersion_rate = 0.2;
f_ext = 5 # Hz
kom_rate = 1*5; # %x sol %x sağ ve kendiysi ile %2x+1 bağ
randCon_rate = 0*0.5;
recData = False # data kayıt edilsin mi?
recPic = False # resimler kayıt edilsin mi?

Noise_sigma = 0*0.1 # std
NSigma = 0*0.1 * nA
NSynSigma = 0*0.05

wn_tau_PY = 2 * ms
wn_tau_IN = 2 * ms
wn_tau_TC = 2 * ms
wn_tau_RE = 2 * ms

wn_mu = 0*0.1
wn_tau_ampa = 2 * ms # 0.5*6*5 * ms
wn_tau_gabaa = 2 * ms # 0.5*6*5 * ms
wn_tau_gabab = 2 * ms # 0.5*250*5 * ms
wn_tau_nmda = 2 * ms # 0.5*100*5 * ms

Amax_ampa_noise = (1-new)*0.088146 + new*0.617986 #
Amax_nmda_noise = (1-new)*0.0951 + new*0.0692 #
Amax_gabaa_noise = (1-new)*0.6412 + new*0.9598 # eğer alfa=20 olursa x=0.85525 olur
Amax_gabab_noise = (1-new)*1.45e-05 + new*2.642e-4 #

ISI_thr = 8 * ms / second
bases=[3/4., 1/4.]

N=100
# PY ve IN hücre sayıları
N_PY = N; N_IN = N/4;
# TC ve RE hücre sayıları
N_TC = N/2; N_RE = N/2;

# sinaptik bağ oranları
#kom_rate = 2*5; # %x sol %x sağ ve kendiysi ile %2x+1 bağ
# Kortikal sinaptik bağ sayıları
kom_AMPA_PY_PY = (N_PY)*kom_rate/100; 
kom_AMPA_PY_IN = (N_PY)*kom_rate/100; 
kom_NMDA_PY_PY = (N_PY)*kom_rate/100; 
kom_NMDA_PY_IN = (N_PY)*kom_rate/100;
kom_GABAA_IN_PY = (N_IN)*kom_rate/100;
 
#Talamik sinaptik bağ sayıları
kom_AMPA_TC_RE = (N_TC)*kom_rate/100;
kom_GABAA_RE_TC = (N_RE)*kom_rate/100; 
kom_GABAB_RE_TC = (N_RE)*kom_rate/100; 
kom_GABAA_RE_RE = (N_RE)*kom_rate/100; 
 
#TalamoKrotikal sinaptik bağ sayıları
kom_AMPA_TC_PY = (N_TC)*kom_rate/100; 
kom_AMPA_TC_IN = (N_TC)*kom_rate/100; 
kom_AMPA_PY_TC = (N_PY)*kom_rate/100; 
kom_AMPA_PY_RE = (N_PY)*kom_rate/100; 

############################ Nernst ####################################
# Nernst parametreleri
sicaklik = 36 * celsius
Temperature = (sicaklik + 273.15 * celsius) # * kelvin
valence = 2 # valance is +1 for Na+, +1 for K+, +2 for Ca2+, -1 for Cl-, 

################################# Kortikal ################################
# PY ve IN sıcaklık ölçekleri
q10 = 2.5; tadj_T_IN = tadj_T_PY = q10 ** ((sicaklik - (24*celsius))/ (10.*celsius))
q10 = 2.3; phi_m_CaL_IN = phi_m_CaL_PY = q10 ** ((sicaklik - (23*celsius))/ (10.*celsius))
q10 = 2.3; phi_h_CaL_IN = phi_h_CaL_PY = q10 ** ((sicaklik - (23*celsius))/ (10.*celsius))
q10 = 2.3; tadj_KCa_IN = tadj_KCa_PY = q10 ** ((sicaklik - (23*celsius))/ (10.*celsius))
q10 = 2.3; tadj_Km_IN = tadj_Km_PY = q10 ** ((sicaklik - (23*celsius))/ (10.*celsius))
q10 = 2.3; tadj_Kv_IN = tadj_Kv_PY = q10 ** ((sicaklik - (23*celsius))/ (10.*celsius))
q10 = 2.3; phi_m_Na_IN = phi_m_Na_PY = q10 ** ((sicaklik - (23*celsius))/ (10.*celsius))
q10 = 2.3; phi_h_Na_IN = phi_h_Na_PY = q10 ** ((sicaklik - (23*celsius))/ (10.*celsius))
q10 = 2.7; phi_m_Nap_IN = phi_m_Nap_PY = q10 ** ((sicaklik - (22*celsius))/ (10.*celsius))
q10 = 3.0; tadj_CAN_IN = tadj_CAN_PY = q10 ** ((sicaklik - (22*celsius))/ (10.*celsius))
#----------------------------------------------------------------

# kompartmanlar arası bağlantı rezistansı
kappa_IN = kappa_PY = 10e3 * kohm; # 1/kohm = mS
# dendritik-aksosomatik alan oranı
rho_PY = 165.; # RS PY için rho=165, FS IN için rho=50 
rho_IN = 50.;

# soma alanı
diam_SOMA_PY = 10/pi * um; L_SOMA_PY = 10 * um; nseg_SOMA_PY = 1;
area_SOMA_PY = pi*(diam_SOMA_PY)*L_SOMA_PY/nseg_SOMA_PY # 1.0e-6 * cm**2;
diam_DEND_PY = 10/pi * um; L_DEND_PY = rho_PY * L_SOMA_PY; nseg_DEND_PY = 1;
area_DEND_PY = 2*pi*(diam_DEND_PY/2)*L_DEND_PY/nseg_DEND_PY #area_SOMA_PY * rho_PY;

# dendrit alanı
diam_SOMA_IN = 10/pi * um; L_SOMA_IN = 10 * um; nseg_SOMA_IN = 1;
area_SOMA_IN = pi*(diam_SOMA_IN)*L_SOMA_IN/nseg_SOMA_IN # 1.0e-6 * cm**2;
diam_DEND_IN = 10/pi * um; L_DEND_IN = rho_IN * L_SOMA_IN; nseg_DEND_IN = 1;
area_DEND_IN = 2*pi*(diam_DEND_IN/2)*L_DEND_IN/nseg_DEND_IN # area_SOMA_IN * rho_IN;

# soma spesifik eksenel özdirenç (ohm*cm)
Ra_SOMA_PY = 100 * ohm * cm
# dendrit yarı-segment eksenel rezistans http://www.neuron.yale.edu/neuron/static/docs/help/neuron/neuron/geometry.html
ri_PY=(Ra_SOMA_PY*((L_DEND_PY/nseg_DEND_PY)/2.)/(pi*(diam_DEND_PY/2.)**2)); # ohm
# dendrit eksenel rezistans
Ra_DEND_PY = Ra_SOMA_PY*kappa_PY/ri_PY; # ohm*cm
 
#kappa_SOMA_PY = L_SOMA_PY*Ra_SOMA_PY/(pi*(diam_SOMA_PY/2)**2);
#kappa_DEND_PY = L_DEND_PY*Ra_DEND_PY/(pi*(diam_DEND_PY/2)**2);

g_compart_PYd = 1.0/(kappa_PY*area_DEND_PY);
g_compart_PYs = 1.0/(kappa_PY*area_SOMA_PY);

g_compart_INd = 1.0/(kappa_IN*area_DEND_IN); 
g_compart_INs = 1.0/(kappa_IN*area_SOMA_IN);
 
Cm_SOMA_PY = Cm_DEND_PY = (0.75 * ufarad * cm ** -2)
Cm_SOMA_IN = Cm_DEND_IN = (0.75 * ufarad * cm ** -2)
 
# spesifik membran direnci
Rm_SOMA_PY = 30 * kohm * cm ** 2 # membran direnci
rm_SOMA_PY = Rm_SOMA_PY / area_SOMA_PY
#----------------------------------------------------------------

# PY reversaller
E_L_PY = -70 * mV
E_KL_PY = -95 * mV
E_Ca_CaL_PY = 140 * mV
E_KCa_PY = -90 * mV
E_Km_PY = -90 * mV
E_Kv_PY = -90 * mV
E_Na_PY = 60 * mV
E_CAN_PY = 10 * mV

#PY iletkenlikler
g_Lk_PYd = (0.033 * msiemens * cm ** -2) # 1./Rm_SOMA_PY 
g_KL_PYd = (0.0025 * msiemens * cm ** -2)  
g_CaL_PYd = (4*0.03 * msiemens * cm ** -2) 
g_KCa_PYd = (0.3 * msiemens * cm ** -2) 
g_Km_PYd = (0.01 * msiemens * cm ** -2)  
g_Na_PYd = (2*1.5 * msiemens * cm ** -2) 
g_Nap_PYd = (0.07 * msiemens * cm ** -2) 

g_CaT_PYd = (0*5e-5 * msiemens * cm ** -2)
g_Kv_PYd = (0*0.1*3.0 * msiemens * cm ** -2)
g_CAN_PYd = (0*0.01 * msiemens * cm ** -2)

g_Na_PYs = (3000 * msiemens * cm ** -2)
g_Nap_PYs = (0.07 * msiemens * cm ** -2)
g_Kv_PYs = (3*200 * msiemens * cm ** -2)

g_Km_PYs = (0*0.02 * msiemens * cm ** -2)
g_KCa_PYs = (0*0.65 * msiemens * cm ** -2)
g_CaL_PYs = (0*0.3 * msiemens * cm ** -2)
#----------------------------------------------------------------

#IN reversaller
E_L_IN = -70 * mV
E_KL_IN = -95 * mV
E_Ca_CaL_IN = 140 * mV
E_KCa_IN = -90 * mV
E_Km_IN = -90 * mV
E_Kv_IN = -90 * mV
E_Na_IN = 60 * mV

#IN iletkenlikler
g_Lk_INd = (0.033 * msiemens * cm ** -2)
g_KL_INd = (0.0025 * msiemens * cm ** -2)
g_CaL_INd = (5*0.03 * msiemens * cm ** -2)
g_KCa_INd = (0.3 * msiemens * cm ** -2)
g_Km_INd = (0.01 * msiemens * cm ** -2)
g_Na_INd = (2*1.5 * msiemens * cm ** -2)
#g_Nap_INd = (0 * msiemens * cm ** -2)

#g_Nap_INs = (0 * msiemens * cm ** -2)
g_Na_INs = (2500 * msiemens * cm ** -2)
g_Kv_INs = (3*200 * msiemens * cm ** -2)
#----------------------------------------------------------------

# Soma ve dendrit [Ca]_dış konsantrasyonu
SCa_o_IN = DCa_o_IN = SCa_o_PY = DCa_o_PY = 2.0 * mMolar #(1e-6 * mole*cm**-3) # 1 M = 1e-3 mole/cm3
# soma ve dendrit [Ca]_dış kararlı durum
SCa_i_inf_IN = DCa_i_inf_IN = SCa_i_inf_PY = DCa_i_inf_PY = 1e-4 * mMolar #2.4e-4 * mmole #(1e-6 * mole * cm**-3) # 1 M = 1e-3 mole/cm3

#A = 1e-5*(cm2/(ms*uA))*mmole# 5.18e-5*(cm2/(ms*uA))*mmole #(1e-6 * mole * cm**-3)
#A = (10*5.18e-5*mMolar*cm2/(ms*uA))

depth_cad_PY = .1 * um; taur_cad_PY = 150. * ms
depth_cad_IN = .1 * um; taur_cad_IN = 200. * ms
#----------------------------------------------------------------

# PY ve IN hücre akım sabitleri 
# kalsiyum bağımlı potasyum akım değişkenleri
Ra_KCa_IN = Ra_KCa_PY = 5*0.01 * 1/ms # maks aktivasyon oranı
Rb_KCa_IN = Rb_KCa_PY = 5*0.02 * 1/ms # maks deaktivasyon oranı
caix_KCa_IN = caix_KCa_PY = 4

# kalsiyum bağımlı nonspesifik CAN akım değişkenleri
Ra_CAN_IN = Ra_CAN_PY = 0.003 * 1/ms # maks aktivasyon oranı
Rb_CAN_IN = Rb_CAN_PY = 0.003 * 1/ms # maks deaktivasyon oranı
caix_CAN_IN = caix_CAN_PY = 8
cac_CAN_IN = cac_CAN_PY = 1e-4 * mMolar

# potasyum Km akım değişkenleri
Ra_Km_IN = Ra_Km_PY = 0.001 * 1/ms # maks aktivasyon oranı
Rb_Km_IN = Rb_Km_PY = 0.001 * 1/ms # maks deaktivasyon oranı
tha_Km_IN = tha_Km_PY = -30 * mV # v 1/2 for inf
qa_Km_IN = qa_Km_PY = 9 * mV # inf slope    
 
# hızlı potasyum Kv akım değişkenleri
Ra_Kv_IN = Ra_Kv_PY = 0.02 * 1/ms # maks aktivasyon oranı
Rb_Kv_IN = Rb_Kv_PY = 0.002 * 1/ms # maks deaktivasyon oranı
tha_Kv_IN = tha_Kv_PY = 25 * mV # v 1/2 for inf
qa_Kv_IN = qa_Kv_PY = 9 * mV # inf slope    
 
# hızlı sodyum Na akım değişkenleri
Vshift_Na_IN = Vshift_Na_PY = -10 * mV     # voltage shift (affects all)
tha_Na_IN = tha_Na_PY = -35 * mV         # v 1/2 for act        (-42)
qa_Na_IN = qa_Na_PY = 9 * mV           # act slope      
Ra_Na_IN = Ra_Na_PY = 0.182 * 1/ms       # open (v)     
Rb_Na_IN = Rb_Na_PY = 0.124 * 1/ms     # close (v)      
thi1_Na_IN = thi1_Na_PY = -50 *mV        # v 1/2 for inact  
thi2_Na_IN = thi2_Na_PY = -75 *mV        # v 1/2 for inact  
qi_Na_IN = qi_Na_PY = 5 *mV            # inact tau slope
thinf_Na_IN = thinf_Na_PY = -65 *mV      # inact inf slope  
qinf_Na_IN = qinf_Na_PY = 6.2 * mV       # inact inf slope
Rg_Na_IN = Rg_Na_PY = 0.0091 * 1/ms    # inact (v)  
Rd_Na_IN = Rd_Na_PY = 0.024 * 1/ms     # inact recov (v) 
 
# sürekli sodyum Nap akım değişkenleri
Tet_Nap_IN = Tet_Nap_PY = -42 * mV
Sig_Nap_IN = Sig_Nap_PY = 5 * mV
f_Nap_IN = f_Nap_PY = 0.02
tau_Dm_Nap_IN = tau_Dm_Nap_PY = 0.8*ms/phi_m_Nap_PY
tau_Sm_Nap_IN = tau_Sm_Nap_PY = 0.8*ms/phi_m_Nap_PY

# düşük eşikli Ca akım değişkenleri
v12m_PY=50 * mV
v12h_PY=78 * mV
vwm_PY=7.4 * mV
vwh_PY=5.0 * mV
am_PY=3 * ms
ah_PY=85 * ms
vm1_PY=25 * mV
vm2_PY=100 * mV
vh1_PY=46 * mV
vh2_PY=405 * mV
wm1_PY=20 * mV
wm2_PY=15 * mV
wh1_PY=4 * mV
wh2_PY=50 * mV

############################### Talamik ####################################
# TC sıcaklık ölçekleri
q10 = 3.55; phi_m_TC= q10 ** ((sicaklik - (24*celsius))/ (10.*celsius))
q10 = 3.0; phi_h_TC= q10 ** ((sicaklik - (24*celsius))/ (10.*celsius))
q10 = 3.0; aadj_TC= q10 ** ((sicaklik - (23.5*celsius))/ (10.*celsius))
q10 = 3.0; hadj_TC= q10 ** ((sicaklik - (36*celsius))/ (10.*celsius))
q10 = 3.0; hadj2_TC= q10 ** ((sicaklik - (35.5*celsius))/ (10.*celsius))
 
#RE sıcaklık ölçekleri
q10 = 5.0; phi_m_RE= q10 ** ((sicaklik - (24*celsius))/ (10.*celsius))
q10 = 3.0; phi_h_RE= q10 ** ((sicaklik - (24*celsius))/ (10.*celsius))
#----------------------------------------------------------------
 
# Ra = 100, nseg = 1, L=diam=96 um, area=pi*L*d
area_TC = 29000 * umetre ** 2 # yada 2.9e-4*(cm**2)
# Ra = 100, nseg = 1, L=64.86um, diam=70 um, area=pi*L*d
area_RE = 14300 * umetre ** 2 # yada 1.43e-4*(cm**2)
 
Cm_TC = (1 * ufarad * cm ** -2) 
Cm_RE = (1 * ufarad * cm ** -2) 
#----------------------------------------------------------------
 
# TC reversaller
EL_TC = -70 * mV
EKL_TC = -(95) * mV
EK_TC = -(95) * mV
ENa_TC = (50) * mV
#ECa_TC = 120 * mV
Eh_TC = -40 * mV
#----------------------------------------------------------------
#TC iletkenlikler
g_L_TC = (0.01 * msiemens * cm ** -2)
#g_KL_TC = ((0.03-0.02) * msiemens * cm ** -2)
g_Na_TC = (90. * msiemens * cm ** -2) 
g_K_TC = ((12.-2) * msiemens * cm ** -2) 
g_T_TC = ((2.2+0.8) * msiemens * cm ** -2) 
#g_h_TC = ((0.017-0.002) * msiemens * cm ** -2) 
g_A_TC = (1.0 * msiemens * cm ** -2)
#g_T2_TC = (1.75 * msiemens * cm ** -2)
 
#RE reversaller
EL_RE = -(77+13) * mV
EKL_RE = -95 * mV
EK_RE = -(100) * mV
ENa_RE = (50) * mV
#ECa_RE = 120 * mV
 
# RE iletkenlikler
g_L_RE = (0.05 * msiemens * cm ** -2) 
#g_KL_RE = ((0.005-0.002) * msiemens * cm ** -2)  
g_Na_RE = ((100.+100) * msiemens * cm ** -2) 
g_K_RE = ((10.+10) * msiemens * cm ** -2) 
g_Ts_RE = ((2.3+0.7) * msiemens * cm ** -2) 
#----------------------------------------------------------------
 
# h-akım sabitleri
nca = 4.            # ca++ iyonların KB proteinlere bağlandığı bölge sayısı. genelde 4 tür.
nexp = 1.           # Ih kanallara balanma sayısı
ginc = 2.           # Ca++ bağlanma durumuna bağlı iletkenlik ayarlama değeri. genelde 2-3 tür (Harigawa & Hirisawa, 1989).
taumin = 20          # tau min değeri
taumax = 1000        # tau max değeri
 
# Ca iyonların KB (Kalsiyum Bağlayıcı) proteinlere bağlanma zaman sabiti
k1 = 2.5e7 * (mMolar**-4)
# Ca iyonların KB proteinlerden ayrılma zaman sabiti
k2 = 0.0004
# KB proteinlerin Ih kanallara bağlanma zaman sabiti
k3 = (0.1)
# KB proteinlerin Ih kanallardan ayrılma zaman sabiti
k4 = 0.001
# KB proteinlerin yarı aktivasyonu (0.001-0.01 mM)
cac = 0.002 * mMolar #(k2/k1)**(1./nca) # CB protein yarı-aktivasyonu; yaklaşık 1 ile 10 microM.
# Ih kanallar için KB proteinlerin yarı aktivasyonu 
Pc = 0.01 #(k4/k3)**(1./nexp)
#----------------------------------------------------------------
 
# [Ca]_dış konsantrasyonu
Ca_o_RE = Ca_o_TC = 2.0 * mMolar # 1 M = 1e-3 mole/cm3
# [Ca]_dış kararlı durum
Ca_i_inf_RE = Ca_i_inf_TC = 2.4e-4 * mMolar # 1 M = 1e-3 mole/cm3
 
#A_TC = A_RE = (5.18e-5*mMolar*cm2/(ms*uA))
 
depth_cad_TC = depth_cad_RE = 1. * um
taur_cad_TC = taur_cad_RE = 5. * ms
#----------------------------------------------------------------
 
VTr_TC = (-40) * mV 
VTrK_TC = (-25) * mV 
VTr_RE = (-50) * mV
VTrK_RE = (-50) * mV
 
V_shift_TC = 2 * mV
V_shift_RE = 2 * mV
 
###################### Sinaptik parametreler ############################
# Kortikal sinaps Reversal potensiyeller
E_AMPA_PY_PY = 0 * mV
E_AMPA_PY_IN = 0 * mV
E_NMDA_PY_PY = 0 * mV
E_NMDA_PY_IN = 0 * mV
E_GABAA_IN_PY = -70 * mV
E_AMPA_EXT_PY = 0 * mV
E_AMPA_EXT_IN = 0 * mV
# Talamik sinaps Reversal potensiyeller
E_AMPA_TC_RE = 0 * mV
E_GABAA_RE_RE = -(70+10) * mV
E_GABAA_RE_TC = -85 * mV
E_GABAB_RE_TC = -95 * mV
E_AMPA_EXT_TC = 0 * mV
E_AMPA_EXT_RE = 0 * mV
# TalamoKortikal sinaps Reversal potensiyeller
E_AMPA_PY_TC = 0 * mV
E_AMPA_PY_RE = 0 * mV
E_AMPA_TC_PY = 0 * mV
E_AMPA_TC_IN = 0 * mV
#----------------------------------------------------------------
# Kortikal sinaps maks iletkenlikler
rk1 = -0.5; rk2 = -0.4;
w_AMPA_PY_PY = (0.1+0.1*rk1) * uS / area_DEND_PY / (2*kom_AMPA_PY_PY+1) # Bazhenov ve ark. 2002
w_AMPA_PY_IN = (0.05+0.05*rk2) * uS / area_DEND_IN / (2*kom_AMPA_PY_IN+1) # Bazhenov ve ark. 2002
w_NMDA_PY_PY = (0.01+0.01*rk1) * uS / area_DEND_PY / (2*kom_NMDA_PY_PY+1) # Bazhenov ve ark. 2002
w_NMDA_PY_IN = (0.008+0.008*rk2) * uS / area_DEND_IN / (2*kom_NMDA_PY_IN+1) # Bazhenov ve ark. 2002
w_GABAA_IN_PY = (0.05+0.05*rk2) * uS / area_DEND_PY / (2*kom_GABAA_IN_PY+1) # Bazhenov ve ark. 2002

w_AMPA_EXT_PY = 0*0.5*0.5 * uS / area_DEND_PY # Bazhenov ve ark. 2002
w_AMPA_EXT_IN = 0*0.5*0.5 * uS / area_DEND_IN # Bazhenov ve ark. 2002

# Talamik sinaps maks iletkenlikler
rt1 = -0.4; rt2 = -0.5;
w_GABAA_RE_TC = (0.2+0.2*rt1) * uS / area_TC / (2*kom_GABAA_RE_TC+1) # Bazhenov ve ark. 2002
w_GABAB_RE_TC = (0.04+0.04*rt1) * uS / area_TC / (2*kom_GABAB_RE_TC+1) # Bazhenov ve ark. 2002
w_AMPA_TC_RE = (0.4+0.4*rt1) * uS / area_RE / (2*kom_AMPA_TC_RE+1) # Bazhenov ve ark. 2002
w_GABAA_RE_RE = (0.2+0.2*rt2) * uS / area_RE / (2*kom_GABAA_RE_RE+1) # Bazhenov ve ark. 2002

w_AMPA_EXT_TC = 0*0.5 * uS / area_TC 
w_AMPA_EXT_RE = 0*0.5 * uS / area_RE 

# TalamoKortikal sinaps maks iletkenlikler
rtk1 = 0.3; rtk2 = -0.2;
w_AMPA_PY_TC = (0.025+0.025*rtk1) * uS / area_TC / (2*kom_AMPA_PY_TC+1) # Bazhenov ve ark. 2002
w_AMPA_PY_RE = (0.05+0.05*rtk2) * uS / area_RE / (2*kom_AMPA_PY_RE+1) # Bazhenov ve ark. 2002
w_AMPA_TC_PY = (0.1+0.1*rtk1) * uS / area_DEND_PY / (2*kom_AMPA_TC_PY+1) # Bazhenov ve ark. 2002
w_AMPA_TC_IN = (0.1+0.1*rtk2) * uS / area_DEND_IN / (2*kom_AMPA_TC_IN+1) # Bazhenov ve ark. 2002

################# Kortikal akım kinetikleri ######################
# PY CaT-akım kinetikleri
# Schaefer ve ark. 2003
mCaTinf_PY=lambda v: 1.0 / ( 1 + exp(-(v+v12m_PY)/vwm_PY))
hCaTinf_PY=lambda v: 1.0 / ( 1 + exp((v+v12h_PY)/vwh_PY))

taumCaT_PY=lambda v: (am_PY + (1.*ms)/(exp((v+vm1_PY)/wm1_PY) + exp(-(v+vm2_PY)/wm2_PY)))
tauhCaT_PY=lambda v: (ah_PY + (1.0*ms)/( exp((v+vh1_PY)/wh1_PY) + exp(-(v+vh2_PY)/wh2_PY))) 
#----------------------------------------------------------------

# Kv-akım kinetikleri
# Potassium channel, Hodgkin-Huxley style kinetics
#Kinetic rates based roughly on Sah et al. and Hamill et al. (1991)
alpham_Kv=lambda v,Ra,tha,qa: Ra * (mV**-1) * (v - tha) / (1 - exp(-(v - tha)/qa))
betam_Kv=lambda v,Rb,tha,qa: -Rb * (mV**-1) * (v - tha) / (1 - exp((v - tha)/qa))
#----------------------------------------------------------------

# Km-akım kinetikleri
# Potassium channel, Hodgkin-Huxley style kinetics Based 
# on I-M (muscarinic K channel) Slow, noninactivating
alpham_Km=lambda v,Ra,tha,qa: Ra * (mV**-1) * (v - tha) / (1 - exp(-(v - tha)/qa))
betam_Km=lambda v,Rb,tha,qa: -Rb * (mV**-1) * (v - tha) / (1 - exp((v - tha)/qa))
#----------------------------------------------------------------

###################### Talamik akım kinetikleri ############################
# TC T-akım kinetikleri
# Bazhenov ve ark. 1998a
mCainf_TC=lambda v:(1.)/(1+exp(-(v+V_shift_TC+57*mV)/(6.2*mV)))
hCainf_TC=lambda v:(1.)/(1+exp((v+V_shift_TC+81*mV)/(4*mV)))
 
# Huguenard and McCormick (1992)
#taumCa_TC=lambda v:((0.612*ms) + (1.*ms)/(exp(-(v+V_shift_TC+129.6*mV)/(16.7*mV))+\
#                exp((v+V_shift_TC+14.8*mV)/(18.2*mV))))/phi_m_TC
# yada Bazhenov ve ark. 1998a - cellular and betwork models ...
#taumCa_TC=lambda v:(0.13*ms) + (0.22*ms)/(exp((-132*mV-v)/(16.7*mV))+\
#                exp((16.8*mV+v)/(18.2*mV)))
 
# Huguenard and McCormick (1992)        
#tauhCa_TC=lambda v:(28.0*ms + (1.*ms)*exp((-22*mV-v)/(10.5*mV)))/tadj_TC if (v>=-80*mV) else \
#                (1.*ms)*exp((467*mV+v)/(66.6*mV))/tadj_TC
# yada Bazhenov ve ark. 1998a - cellular and network models ...
#tauhCa_TC=lambda v:(8.2*ms + (1.*ms)*((56.6 + 0.27*exp((115.2*mV+v)/(5*mV)))/\
#                (1+exp((86*mV+v)/(3.2*mV)))))
# yada Destexhe ve ark. 1996a Ionic mechanisms underlying
## ****** bu modelde aktivasyonun kararlı durumu bulunuyor mCa_inf *********
tauhCa_TC=lambda v:(30.8*ms + (1.*ms)*((211.4 + exp((v+V_shift_TC+113.2*mV)/(5*mV)))/\
                (1+exp((v+V_shift_TC+84*mV)/(3.2*mV)))))/phi_h_TC
#----------------------------------------------------------------
 
# RE T-akım kinetikleri
# Bazhenov ve ark. 1998a ve Destexhe ve ark. 1996b
mCainf_RE=lambda v:(1.)/(1+exp(-(v+V_shift_RE+50*mV)/(7.4*mV)))
hCainf_RE=lambda v:(1.)/(1+exp((v+V_shift_RE+78*mV)/(5*mV)))
 
# Time constants were obtained from J. Huguenard
# link : http://senselab.med.yale.edu/modeldb/showmodel.asp?model=3343&file=\DLGN_NEW\IT2.mod
taumCa_RE=lambda v:((3.0*ms) + (1.0*ms)/(exp((v+V_shift_RE+25*mV)/(10*mV))+\
                exp(-(v+V_shift_RE+100*mV)/(15*mV))))/phi_m_RE
# yada Bazhenov ve ark. 1998a - cellular and betwork models ... ile
# Destexhe ve ark. 1996b In Vivo, ln vitro, and Computational Analysis ... 
#taumCa_RE=lambda v:(1.0*ms) + (0.33*ms)/(exp((27*mV+v)/(10*mV))+\
#                exp((-102*mV-v)/(15*mV)))
 
# Time constants were obtained from J. Huguenard
# link : http://senselab.med.yale.edu/modeldb/showmodel.asp?model=3343&file=\DLGN_NEW\IT2.mod
tauhCa_RE=lambda v:((85.0*ms) + (1.0*ms)/(exp((v+V_shift_RE+46*mV)/(4*mV))+\
                exp(-(v+V_shift_RE+405*mV)/(50*mV))))/phi_h_RE
# yada Bazhenov ve ark. 1998a - cellular and betwork models ...
#tauhCa_RE=lambda v:(22.7*ms + (0.27*ms)/(exp((48*mV+v)/(4*mV))+\
#                exp((-407*mV-v)/(50*mV))))
# yada Destexhe ve ark. 1996b In Vivo, ln vitro, and Computational Analysis ... 
#tauhCa_RE=lambda v:(28.3*ms + (0.33*ms)/(exp((48*mV+v)/(4*mV))+\
#                exp((-407*mV-v)/(50*mV))))
#----------------------------------------------------------------
 
# TC A-akım
mAKinf_TC=lambda v:(1.)/(1+exp((-60*mV-v)/(8.5*mV)))
hAKinf_TC=lambda v:(1.)/(1+exp((78.*mV+v)/(6.0*mV)))
 
# Bazhenov ve ark. 1998a
#tauam_TC=lambda v:(0.1*ms) + (0.27*ms)/(exp((35.8*mV+v)/(19.7*mV))+\
#                exp((-79.7*mV-v)/(12.7*mV)))
# yada Huguenard and McCormick (1992)
taumAK_TC=lambda v:((0.37*ms) + (1.*ms)/(exp((35.82*mV+v)/(19.69*mV))+\
                exp((-79.69*mV-v)/(12.7*mV))))/aadj_TC
 
# Bazhenov ve ark. 1998a
#tauah_TC=lambda v:(0.27*ms)/(exp((46*mV+v)/(5*mV))+exp((-238*mV-v)/(37.5*mV))) if (v<-63*mV) else \
#               (1.*ms*5.1)
# yada Huguenard and McCormick (1992)
tauhAK_TC=lambda v: 1.*(v<-63*mV)*(1.*ms)/(exp((46.05*mV+v)/(5*mV))+exp((-238.4*mV-v)/(37.45*mV)))/aadj_TC \
               + 1.*(1-(v<-63*mV))*(1.*ms*19.)/aadj_TC
#----------------------------------------------------------------
""" 
#################### Harici DC girişler############################

# PY Soma harici uyartı
i_cur1_PY = 0.1/area_SOMA_PY
eqs_input_PY='''
B = x*nA : amp
x = 0 if y<100 else i_cur1_PY if y<900 else 0 : 1
dy/dt = (1)/ms : 1
'''
input_PYs=NeuronGroup(1,eqs_input_PY,clock=myclock)
 
# PY Dendrit harici uyartı
i_cur2_PY = 0.1/area_DEND_PY
eqs_input_PY='''
B = x*nA : amp
x = 0 if y<100 else i_cur2_PY if y<900 else 0 : 1
dy/dt = (1)/ms : 1
'''
input_PYd=NeuronGroup(1,eqs_input_PY,clock=myclock)
#----------------------------------------------------------------
 
# IN Soma harici uyartı
i_cur1_IN = 0.1/area_SOMA_IN
eqs_input_IN='''
B = x*nA : amp
x = 0 if y<100 else i_cur1_IN if y<900 else 0 : 1
dy/dt = (1)/ms : 1
'''
input_INs=NeuronGroup(1,eqs_input_IN,clock=myclock)
#----------------------------------------------------------------
 
# TC harici uyartı
i_cur1_TC = 0.1/area_TC
eqs_input_TC='''
B = x*nA : amp
x = 0 if y<100 else i_cur1_TC if y<1000 else 0 : 1
dy/dt = (1)/ms : 1
'''
input_TC=NeuronGroup(1,eqs_input_TC, clock=myclock)
#----------------------------------------------------------------
 
# RE harici uyartı
i_cur1_RE = 0.1/area_RE
eqs_input_RE='''
B = x*nA : amp
x = 0 if y<10 else i_cur1_RE if y<500 else 0 : 1
dy/dt = (1)/ms : 1
'''
input_RE=NeuronGroup(1,eqs_input_RE, clock=myclock)
#----------------------------------------------------------------
"""
# modelin rastsal parametreleri
w_AMPA_EXT_PY_rnd = np.random.uniform(0.25,1.0,N_PY)
w_AMPA_EXT_IN_rnd = np.random.uniform(0.25,1.0,N_IN)
w_AMPA_EXT_TC_rnd = np.random.uniform(0.25,1.0,N_TC)
w_AMPA_EXT_RE_rnd = np.random.uniform(0.25,1.0,N_RE)

Noise_sigma_tot_PY = NSigma/area_DEND_PY
Noise_sigma_tot_IN = NSigma/area_DEND_IN
Noise_sigma_tot_TC = NSigma/area_TC
Noise_sigma_tot_RE = NSigma/area_RE

######################## Hücre modelleri ##################################
 
# PY nöron modeli
# Bazhenov ve ark. (2002;1998b)
eqs_PY = Equations('''
dvd/dt = (Id_ext-I_syn-I_int_dend+g_compart_PYd*(vs-vd))/Cm_DEND_PY : volt
dvs/dt = (Is_ext-I_int_soma-g_compart_PYs*(vs-vd))/(Cm_SOMA_PY) : volt
#dvs/dt = (-vs+((vd + kappa_PY * area_SOMA_PY * g2_SOMA) / (1. + kappa_PY*area_SOMA_PY * g1_SOMA)))/ms : volt

#g1_SOMA = ((phi_m_Na_PY*g_Na_PYs*(SmNa**3)*ShNa) + tadj_Kv_PY*(g_Kv_PYs*SmKv) + \
#    (g_Nap_PYs*SmNap)+tadj_Km_PY*g_Km_PYs*mKmS) : siemens/cm2
#g2_SOMA = (Is_ext + (phi_m_Na_PY*g_Na_PYs*(SmNa**3)*ShNa)*E_Na_PY +\
#    tadj_Kv_PY*(g_Kv_PYs*SmKv)*E_Kv_PY+(g_Nap_PYs*SmNap)*E_Na_PY )+\
#    tadj_Km_PY*g_Km_PYs*mKmS*E_Km_PY : amp/cm2 

I_int_soma = Is_Na+Is_Kv+Is_Nap+Is_Km+Is_KCa+Is_CaL : amp/cm2
I_int_dend = Id_Lk+Id_KL+Id_CaL+Id_KCa+Id_Km+Id_Na+Id_Nap+Id_CaT+Id_Kv+Id_CAN : amp/cm2

dwn/dt = -(wn)/wn_tau_PY + (((2*(Noise_sigma_tot_PY)**2)/wn_tau_PY)**.5*xi) : amp/cm2

I_syn = (I_AMPA_PY_PY+I_NMDA_PY_PY+I_GABAA_IN_PY+I_AMPA_EXT_PY+I_AMPA_TC_PY)+wn : amp/cm2
#----------------------------------------------------------------
Id_Lk = g_Lk_PYd*(vd-E_L_PY) : amp/cm2
Id_KL = g_KL_PYd*(vd-E_KL_PY) : amp/cm2
Id_CaL = phi_m_CaL_PY*g_CaL_PYd*DmCaL*DmCaL*DhCaL*(vd-E_Ca_CaL_PY) : amp/cm2
Id_KCa = tadj_KCa_PY * g_KCa_PYd * DmKCa * (vd - E_KCa_PY) : amp/cm2
Id_Km = tadj_Km_PY*g_Km_PYd*DmKm*(vd-E_Km_PY) : amp/cm2
Id_Na = phi_m_Na_PY*g_Na_PYd*(DmNa**3)*DhNa*(vd-E_Na_PY) : amp/cm2
Id_Nap = g_Nap_PYd*DmNap*(vd-E_Na_PY) : amp/cm2

Id_CaT = g_CaT_PYd*(mCaT**2)*hCaT*(vd-EDCa_PY) : amp/cm2
Id_Kv = tadj_Kv_PY*g_Kv_PYd*DmKv*(vd-E_Kv_PY) : amp/cm2
Id_CAN = g_CAN_PYd * DmCAN * (vd - E_CAN_PY) : amp/cm2

Is_Na = phi_m_Na_PY*g_Na_PYs*(SmNa**3)*ShNa*(vs-E_Na_PY) : amp/cm2
Is_Nap = g_Nap_PYs*SmNap*(vs-E_Na_PY) : amp/cm2
Is_Kv = tadj_Kv_PY*g_Kv_PYs*SmKv*(vs-E_Kv_PY) : amp/cm2

Is_Km = tadj_Km_PY*g_Km_PYs*mKmS*(vs-E_Km_PY) : amp/cm2
Is_KCa = tadj_KCa_PY * g_KCa_PYs * SmKCa * (vs - E_KCa_PY) : amp/cm2
Is_CaL = phi_m_CaL_PY*g_CaL_PYs*SmCaL*SmCaL*ShCaL*(vs-E_Ca_CaL_PY) : amp/cm2
#----------------------------------------------------------------
dwn_AMPA_PY_PY/dt = -(wn_AMPA_PY_PY)/wn_tau_ampa + (Amax_ampa_noise*w_AMPA_PY_PY)*(((2*(NSynSigma)**2)/wn_tau_ampa)**.5*xi) : siemens/cm2
g_temp1 = g_AMPA_PY_PY+wn_AMPA_PY_PY : siemens/cm2
I_AMPA_PY_PY = (g_temp1>0)*g_temp1*(vd-E_AMPA_PY_PY) : amp/cm2

dwn_NMDA_PY_PY/dt = -(wn_NMDA_PY_PY)/wn_tau_nmda + (Amax_nmda_noise*w_NMDA_PY_PY)*(((2*(NSynSigma)**2)/wn_tau_nmda)**.5*xi) : siemens/cm2
g_temp2 = g_NMDA_PY_PY+wn_NMDA_PY_PY : siemens/cm2
I_NMDA_PY_PY = (g_temp2>0)*g_temp2*(vd-E_NMDA_PY_PY) : amp/cm2

dwn_GABAA_IN_PY/dt = -(wn_GABAA_IN_PY)/wn_tau_gabaa + (Amax_gabaa_noise*w_GABAA_IN_PY)*(((2*(NSynSigma)**2)/wn_tau_gabaa)**.5*xi) : siemens/cm2
g_temp3 = g_GABAA_IN_PY+wn_GABAA_IN_PY : siemens/cm2
I_GABAA_IN_PY = (g_temp3>0)*g_temp3*(vd-E_GABAA_IN_PY) : amp/cm2

dwn_AMPA_EXT_PY/dt = -(wn_AMPA_EXT_PY)/wn_tau_ampa + (Amax_ampa_noise*w_AMPA_EXT_PY)*(((2*(NSynSigma)**2)/wn_tau_ampa)**.5*xi) : siemens/cm2
g_temp4 = g_AMPA_EXT_PY+wn_AMPA_EXT_PY : siemens/cm2
I_AMPA_EXT_PY = (g_temp4>0)*g_temp4*(vd-E_AMPA_EXT_PY) : amp/cm2

dwn_AMPA_TC_PY/dt = -(wn_AMPA_TC_PY)/wn_tau_ampa + (Amax_ampa_noise*w_AMPA_TC_PY)*(((2*(NSynSigma)**2)/wn_tau_ampa)**.5*xi) : siemens/cm2
g_temp5 = g_AMPA_TC_PY+wn_AMPA_TC_PY : siemens/cm2
I_AMPA_TC_PY = (g_temp5>0)*g_temp5*(vd-E_AMPA_TC_PY) : amp/cm2

LFP = abs(I_AMPA_PY_PY)+abs(I_NMDA_PY_PY)+abs(I_GABAA_IN_PY)+ \
        abs(I_AMPA_TC_PY)+abs(I_AMPA_EXT_PY) : amp/cm2

#gurultusuz sinaptik akımlar
I_AMPA_PY_PYx = g_AMPA_PY_PYx*(vd-E_AMPA_PY_PY) : amp/cm2
I_NMDA_PY_PYx = g_NMDA_PY_PYx*(vd-E_NMDA_PY_PY) : amp/cm2
I_GABAA_IN_PYx = g_GABAA_IN_PYx*(vd-E_GABAA_IN_PY) : amp/cm2
I_AMPA_TC_PYx = g_AMPA_TC_PYx*(vd-E_AMPA_TC_PY) : amp/cm2
I_AMPA_EXT_PYx = g_AMPA_EXT_PYx*(vd-E_AMPA_EXT_PY) : amp/cm2
LFPx = abs(I_AMPA_PY_PYx)+abs(I_NMDA_PY_PYx)+abs(I_GABAA_IN_PYx)+\
        abs(I_AMPA_TC_PYx)+abs(I_AMPA_EXT_PYx) : amp/cm2
#----------------------------------------------------------------
 
# Dendritik sürekli sodyum Na akım kinetikleri
dDmNap/dt = -(DmNap - DmNap_inf)/tau_Dm_Nap_PY : 1 
DmNap_inf = f_Nap_PY/(1 + exp(-(vd - Tet_Nap_PY)/Sig_Nap_PY)) : 1

# Aksosomatik sürekli sodyum Na akım kinetikleri
dSmNap/dt = -(SmNap - SmNap_inf)/tau_Sm_Nap_PY : 1 
SmNap_inf = f_Nap_PY/(1 + exp(-(vs - Tet_Nap_PY)/Sig_Nap_PY)) : 1
#----------------------------------------------------------------
 
# Dendritik hızlı sodyum Na akım kinetikleri
dDmNa/dt = -(DmNa - DmNa_inf)/tau_DmNa : 1 
dDhNa/dt = -(DhNa - DhNa_inf)/tau_DhNa : 1 

DmNa_inf = alphaDmNa/(alphaDmNa+betaDmNa) : 1
DhNa_inf = (1/(1+exp(((vd+Vshift_Na_PY)-thinf_Na_PY)/qinf_Na_PY))) : 1
tau_DhNa = (1/(alphaDhNa+betaDhNa))/phi_h_Na_PY : second
tau_DmNa = (1/(alphaDmNa+betaDmNa))/phi_m_Na_PY : second

alphaDmNa = 1.*(fabs((vd+Vshift_Na_PY)/(tha_Na_PY))>1e-6)*((Ra_Na_PY*(mV**-1)*\
            ((vd+Vshift_Na_PY)-(tha_Na_PY)))/(1-exp(-((vd+Vshift_Na_PY)-(tha_Na_PY))/\
            qa_Na_PY))) + 1.*(1-(fabs((vd+Vshift_Na_PY)/(tha_Na_PY))>1e-6))*(Ra_Na_PY*qa_Na_PY*(mV**-1)) : Hz
#alphaDmNa = ((Ra_Na_PY*(mV**-1)*((vd+Vshift_Na_PY)-(tha_Na_PY)))/(1-exp(-((vd+Vshift_Na_PY)-(tha_Na_PY))/\
#            qa_Na_PY))) : Hz

betaDmNa = 1.*(fabs(-(vd+Vshift_Na_PY)/(-tha_Na_PY))>1e-6)*((Rb_Na_PY*(mV**-1)*\
            (-(vd+Vshift_Na_PY)-(-tha_Na_PY)))/(1-exp(-(-(vd+Vshift_Na_PY)-(-tha_Na_PY))/\
            qa_Na_PY))) + 1.*(1-(fabs(-(vd+Vshift_Na_PY)/(-tha_Na_PY))>1e-6))*(Rb_Na_PY*qa_Na_PY*(mV**-1)) : Hz
#betaDmNa = ((Rb_Na_PY*(mV**-1)*(-(vd+Vshift_Na_PY)-(-tha_Na_PY)))/(1-exp(-(-(vd+Vshift_Na_PY)-(-tha_Na_PY))/\
#            qa_Na_PY))) : Hz

alphaDhNa = 1.*(fabs((vd+Vshift_Na_PY)/(thi1_Na_PY))>1e-6)*((Rd_Na_PY*(mV**-1)*\
            ((vd+Vshift_Na_PY)-(thi1_Na_PY)))/(1-exp(-((vd+Vshift_Na_PY)-(thi1_Na_PY))/\
            qi_Na_PY))) + 1.*(1-(fabs((vd+Vshift_Na_PY)/(thi1_Na_PY))>1e-6))*(Rd_Na_PY*qi_Na_PY*(mV**-1)) : Hz
#alphaDhNa = ((Rd_Na_PY*(mV**-1)*((vd+Vshift_Na_PY)-(thi1_Na_PY)))/(1-exp(-((vd+Vshift_Na_PY)-(thi1_Na_PY))/\
#            qi_Na_PY))) : Hz

betaDhNa = 1.*(fabs(-(vd+Vshift_Na_PY)/(-thi2_Na_PY))>1e-6)*((Rg_Na_PY*(mV**-1)*\
            (-(vd+Vshift_Na_PY)-(-thi2_Na_PY)))/(1-exp(-(-(vd+Vshift_Na_PY)-(-thi2_Na_PY))/\
            qi_Na_PY))) + 1.*(1-(fabs(-(vd+Vshift_Na_PY)/(-thi2_Na_PY))>1e-6))*(Rg_Na_PY*qi_Na_PY*(mV**-1)) : Hz
#betaDhNa = ((Rg_Na_PY*(mV**-1)*(-(vd+Vshift_Na_PY)-(-thi2_Na_PY)))/(1-exp(-(-(vd+Vshift_Na_PY)-(-thi2_Na_PY))/\
#            qi_Na_PY))) : Hz

# Aksosomatik hızlı sodyum Na akım kinetikleri
dSmNa/dt = -(SmNa - SmNa_inf)/tau_SmNa : 1 
dShNa/dt = -(ShNa - ShNa_inf)/tau_ShNa : 1 

SmNa_inf = alphaSmNa/(alphaSmNa+betaSmNa) : 1
ShNa_inf = (1/(1+exp(((vs+Vshift_Na_PY)-thinf_Na_PY)/qinf_Na_PY))) : 1
tau_SmNa = (1/(alphaSmNa+betaSmNa))/phi_m_Na_PY : second
tau_ShNa = (1/(alphaShNa+betaShNa))/phi_h_Na_PY : second

alphaSmNa = 1.*(fabs((vs+Vshift_Na_PY)/(tha_Na_PY))>1e-6)*((Ra_Na_PY*(mV**-1)*\
            ((vs+Vshift_Na_PY)-(tha_Na_PY)))/(1-exp(-((vs+Vshift_Na_PY)-(tha_Na_PY))/\
            qa_Na_PY))) + 1.*(1-(fabs((vs+Vshift_Na_PY)/(tha_Na_PY))>1e-6))*(Ra_Na_PY*qa_Na_PY*(mV**-1)) : Hz
#alphaSmNa = ((Ra_Na_PY*(mV**-1)*((vs+Vshift_Na_PY)-(tha_Na_PY)))/(1-exp(-((vs+Vshift_Na_PY)-(tha_Na_PY))/\
#            qa_Na_PY))) : Hz

betaSmNa = 1.*(fabs(-(vs+Vshift_Na_PY)/(-tha_Na_PY))>1e-6)*((Rb_Na_PY*(mV**-1)*\
            (-(vs+Vshift_Na_PY)-(-tha_Na_PY)))/(1-exp(-(-(vs+Vshift_Na_PY)-(-tha_Na_PY))/\
            qa_Na_PY))) + 1.*(1-(fabs(-(vs+Vshift_Na_PY)/(-tha_Na_PY))>1e-6))*(Rb_Na_PY*qa_Na_PY*(mV**-1)) : Hz
#betaSmNa = ((Rb_Na_PY*(mV**-1)*(-(vs+Vshift_Na_PY)-(-tha_Na_PY)))/(1-exp(-(-(vs+Vshift_Na_PY)-(-tha_Na_PY))/\
#            qa_Na_PY))) : Hz

alphaShNa = 1.*(fabs((vs+Vshift_Na_PY)/(thi1_Na_PY))>1e-6)*((Rd_Na_PY*(mV**-1)*\
            ((vs+Vshift_Na_PY)-(thi1_Na_PY)))/(1-exp(-((vs+Vshift_Na_PY)-(thi1_Na_PY))/\
            qi_Na_PY))) + 1.*(1-(fabs((vs+Vshift_Na_PY)/(thi1_Na_PY))>1e-6))*(Rd_Na_PY*qi_Na_PY*(mV**-1)) : Hz
#alphaShNa = ((Rd_Na_PY*(mV**-1)*((vs+Vshift_Na_PY)-(thi1_Na_PY)))/(1-exp(-((vs+Vshift_Na_PY)-(thi1_Na_PY))/\
#            qi_Na_PY))) : Hz

betaShNa = 1.*(fabs(-(vs+Vshift_Na_PY)/(-thi2_Na_PY))>1e-6)*((Rg_Na_PY*(mV**-1)*\
            (-(vs+Vshift_Na_PY)-(-thi2_Na_PY)))/(1-exp(-(-(vs+Vshift_Na_PY)-(-thi2_Na_PY))/\
            qi_Na_PY))) + 1.*(1-(fabs(-(vs+Vshift_Na_PY)/(-thi2_Na_PY))>1e-6))*(Rg_Na_PY*qi_Na_PY*(mV**-1)) : Hz
#betaShNa = ((Rg_Na_PY*(mV**-1)*(-(vs+Vshift_Na_PY)-(-thi2_Na_PY)))/(1-exp(-(-(vs+Vshift_Na_PY)-(-thi2_Na_PY))/\
#            qi_Na_PY))) : Hz
#----------------------------------------------------------------

# I_CaT akım kinetikleri kinetikleri
dmCaT/dt=(mCaT_inf-mCaT)/tau_mCaT : 1
mCaT_inf=mCaTinf_PY(vd) : 1   
tau_mCaT=taumCaT_PY(vd) : second

dhCaT/dt=(hCaT_inf-hCaT)/tau_hCaT : 1
hCaT_inf=hCaTinf_PY(vd) : 1   
tau_hCaT=tauhCaT_PY(vd) : second
#----------------------------------------------------------------

# Dendritik Hight-threshold Ca2+ CaL-akım kinetikleri
dDmCaL/dt = -(DmCaL-DmCaL_inf)/tau_DmCaL : 1
dDhCaL/dt = -(DhCaL-DhCaL_inf)/tau_DhCaL : 1

alphaDmCaL = 0.055*(mV**-1)*(-27*mV-vd)/ \
    (exp((-27*mV-vd)/(3.8*mV))-1.)/ms : Hz
betaDmCaL = 0.94*exp((-75*mV-vd)/(17*mV))/ms : Hz
tau_DmCaL = (1/(alphaDmCaL+betaDmCaL))/phi_m_CaL_PY : second
DmCaL_inf = (alphaDmCaL/(alphaDmCaL+betaDmCaL)) : 1

alphaDhCaL = 0.000457*exp((-13*mV-vd)/(50*mV))/ms : Hz
betaDhCaL = 0.0065/(1+exp((-vd-15*mV)/(28*mV)))/ms : Hz
tau_DhCaL = (1/(alphaDhCaL+betaDhCaL))/phi_h_CaL_PY : second
DhCaL_inf = (alphaDhCaL/(alphaDhCaL+betaDhCaL)) : 1

# akso-somatik Hight-threshold Ca2+ CaL-akım kinetikleri
dSmCaL/dt = -(SmCaL-SmCaL_inf)/tau_SmCaL : 1
dShCaL/dt = -(ShCaL-ShCaL_inf)/tau_ShCaL : 1

alphaSmCaL = 0.055*(mV**-1)*(-27*mV-vs)/ \
    (exp((-27*mV-vs)/(3.8*mV))-1.)/ms : Hz
betaSmCaL = 0.94*exp((-75*mV-vs)/(17*mV))/ms : Hz
tau_SmCaL = (1/(alphaSmCaL+betaSmCaL))/phi_m_CaL_PY : second
SmCaL_inf = (alphaSmCaL/(alphaSmCaL+betaSmCaL)) : 1

alphaShCaL = 0.000457*exp((-13*mV-vs)/(50*mV))/ms : Hz
betaShCaL = 0.0065/(1+exp((-vs-15*mV)/(28*mV)))/ms : Hz
tau_ShCaL = (1/(alphaShCaL+betaShCaL))/phi_h_CaL_PY : second
ShCaL_inf = (alphaShCaL/(alphaShCaL+betaShCaL)) : 1
#----------------------------------------------------------------
 
# Ca_i dinamiği
# dendrit Ca_i dinamiği
#dDCa_i/dt = -(A*(Id_CaL))+(DCa_i_inf_PY-DCa_i)/(taur_cad_PY) : Molar #mole/cm3
dDCa_i/dt = -(Id_CaL+Id_CaT)/(2*F*depth_cad_PY)+(DCa_i_inf_PY-DCa_i)/(taur_cad_PY) : Molar #mole/cm3
#dDCa_i/dt = -(1e4*I_T/mA)/(2*96489*1*ms)+((2.4e-4)-DCa_i)/(taur_cad_PY) : 1

# soma Ca_i dinamiği
dSCa_i/dt = -(Is_CaL)/(2*F*depth_cad_PY)+(SCa_i_inf_PY-SCa_i)/(taur_cad_PY) : Molar #mole/cm3
#----------------------------------------------------------------
 
# Kalsiyum reversal potansiyel dinamiği
# dendrit ve akso-somatik Kalsiyum reversal potansiyel dinamiği
dEDCa_PY/dt = (Nernst(valence,Temperature, DCa_i,DCa_o_PY)-EDCa_PY)/(ms) : volt

dESCa_PY/dt = (Nernst(valence,Temperature, SCa_i,SCa_o_PY)-ESCa_PY)/(ms) : volt
#----------------------------------------------------------------

# Aksosomatik hızlı potasyum Kv akım kinetikleri
# Dendritik hızlı potasyum Kv akım kinetikleri
dDmKv/dt = -(DmKv - DmKv_inf)/tau_DmKv : 1 
tau_DmKv = (1/(alphaDmKv+betaDmKv))/tadj_Kv_PY : second
DmKv_inf = alphaDmKv/(alphaDmKv+betaDmKv) : 1

alphaDmKv = alpham_Kv(vd,Ra_Kv_PY,tha_Kv_PY,qa_Kv_PY) : Hz
betaDmKv = betam_Kv(vd,Rb_Kv_PY,tha_Kv_PY,qa_Kv_PY) : Hz

# Aksosomatik hızlı potasyum Kv akım kinetikleri
dSmKv/dt = -(SmKv - SmKv_inf)/tau_SmKv : 1 
tau_SmKv = (1/(alphaSmKv+betaSmKv))/tadj_Kv_PY : second
SmKv_inf = alphaSmKv/(alphaSmKv+betaSmKv) : 1

alphaSmKv = alpham_Kv(vs,Ra_Kv_PY,tha_Kv_PY,qa_Kv_PY) : Hz
betaSmKv = betam_Kv(vs,Rb_Kv_PY,tha_Kv_PY,qa_Kv_PY) : Hz
#----------------------------------------------------------------
 
# Dendritik potasyum Km-akım kinetikleri
dDmKm/dt = -(DmKm - DmKm_inf)/tau_DmKm : 1 
tau_DmKm = (1/(alphaDmKm+betaDmKm))/tadj_Km_PY : second
DmKm_inf = alphaDmKm/(alphaDmKm+betaDmKm) : 1

alphaDmKm = alpham_Km(vd,Ra_Km_PY,tha_Km_PY,qa_Km_PY) : Hz
betaDmKm = betam_Km(vd,Rb_Km_PY,tha_Km_PY,qa_Km_PY) : Hz

# akso-somatik potasyum Km-akım kinetikleri
dmKmS/dt = -(mKmS - mKmS_inf)/tau_mKmS : 1 
tau_mKmS = (1/(alphamKmS+betamKmS))/tadj_Km_PY : second
mKmS_inf = alphamKmS/(alphamKmS+betamKmS) : 1

alphamKmS = alpham_Km(vs,Ra_Km_PY,tha_Km_PY,qa_Km_PY) : Hz
betamKmS = betam_Km(vs,Rb_Km_PY,tha_Km_PY,qa_Km_PY) : Hz
#----------------------------------------------------------------
 
# Dendritik Ca-bağımlı potasyum KCa-akım kinetikleri
dDmKCa/dt = -(DmKCa-DmKCa_inf)/tau_DmKCa : 1
DmKCa_inf = (alphaDmKCa/(alphaDmKCa+betaDmKCa)) : 1
tau_DmKCa = (1/(alphaDmKCa+betaDmKCa))/tadj_KCa_PY : second

alphaDmKCa = (Ra_KCa_PY*(DCa_i/mMolar)**caix_KCa_PY) : Hz
betaDmKCa = Rb_KCa_PY : Hz

# akso-somatik Ca-bağımlı potasyum KCa-akım kinetikleri
dSmKCa/dt = -(SmKCa-SmKCa_inf)/tau_SmKCa : 1
SmKCa_inf = (alphaSmKCa/(alphaSmKCa+betaSmKCa)) : 1
tau_SmKCa = (1/(alphaSmKCa+betaSmKCa))/tadj_KCa_PY : second

alphaSmKCa = (Ra_KCa_PY*(SCa_i/mMolar)**caix_KCa_PY) : Hz
betaSmKCa = Rb_KCa_PY : Hz
#----------------------------------------------------------------

# Dendritik Ca-bağımlı nonspesifik katyonik CAN-akım kinetikleri
dDmCAN/dt = -(DmCAN-DmCAN_inf)/tau_DmCAN : 1
DmCAN_inf = (alphaDmCAN/(alphaDmCAN+betaDmCAN)) : 1
tau_DmCAN = (1/(alphaDmCAN+betaDmCAN))/tadj_CAN_PY : second

alphaDmCAN = (Ra_CAN_PY*(DCa_i/cac_CAN_PY)**caix_CAN_PY) : Hz
betaDmCAN = Rb_CAN_PY : Hz

# akso-somatik Ca-bağımlı nonspesifik katyonik CAN-akım kinetikleri
dSmCAN/dt = -(SmCAN-SmCAN_inf)/tau_SmCAN : 1
SmCAN_inf = (alphaSmCAN/(alphaSmCAN+betaSmCAN)) : 1
tau_SmCAN = (1/(alphaSmCAN+betaSmCAN))/tadj_CAN_PY : second

alphaSmCAN = (Ra_CAN_PY*(SCa_i/cac_CAN_PY)**caix_CAN_PY) : Hz
betaSmCAN = Rb_CAN_PY : Hz
#----------------------------------------------------------------

g_AMPA_PY_PY : siemens/cm2
g_NMDA_PY_PY : siemens/cm2
g_GABAA_IN_PY : siemens/cm2
g_AMPA_TC_PY : siemens/cm2
g_AMPA_EXT_PY : siemens/cm2

#gürültüsüz sinaptik iletkenlikler
g_AMPA_PY_PYx : siemens/cm2
g_NMDA_PY_PYx : siemens/cm2
g_GABAA_IN_PYx : siemens/cm2
g_AMPA_TC_PYx : siemens/cm2
g_AMPA_EXT_PYx : siemens/cm2

Id_ext : amp/cm2
Is_ext : amp/cm2
 
''')
#----------------------------------------------------------------
 
# IN nöron modeli
# Bazhenov ve ark. (2002;1998b)
eqs_IN = Equations('''
dvd/dt = (Id_ext-I_syn-I_int_dend+g_compart_INd*(vs-vd))/Cm_DEND_IN : volt
dvs/dt = (Is_ext-I_int_soma-g_compart_INs*(vs-vd))/(Cm_SOMA_IN) : volt
#dvs/dt = (-vs+((vd + kappa_IN * area_SOMA_IN * g2_SOMA) / (1. + kappa_IN*area_SOMA_IN * g1_SOMA)))/ms : volt

#g1_SOMA = ((phi_m_Na_IN*g_Na_INs*(SmNa**3)*ShNa) + tadj_Kv_IN*(g_Kv_INs*SmKv) ) : siemens/cm2
#g2_SOMA = (Is_ext + (phi_m_Na_IN*g_Na_INs*(SmNa**3)*ShNa)*E_Na_IN +tadj_Kv_IN*(g_Kv_INs*SmKv)*E_Kv_IN ) : amp/cm2 

I_int_soma = Is_Na+Is_Kv : amp/cm2
I_int_dend = Id_Lk+Id_KL+Id_CaL+Id_KCa+Id_Km+Id_Na : amp/cm2
 
dwn/dt = -(wn)/wn_tau_IN + (((2*(Noise_sigma_tot_IN)**2)/wn_tau_IN)**.5*xi) : amp/cm2

I_syn = (I_AMPA_PY_IN+I_NMDA_PY_IN+I_AMPA_EXT_IN+I_AMPA_TC_IN)+wn : amp/cm2
#----------------------------------------------------------------
 
Id_Lk = g_Lk_INd*(vd-E_L_IN) : amp/cm2
Id_KL = g_KL_INd*(vd-E_KL_IN) : amp/cm2
Id_CaL = phi_m_CaL_IN*g_CaL_INd*DmCaL*DmCaL*DhCaL*(vd-E_Ca_CaL_IN) : amp/cm2
Id_KCa = tadj_KCa_IN * g_KCa_INd * DmKCa * (vd - E_KCa_IN) : amp/cm2
Id_Km = tadj_Km_IN*g_Km_INd*DmKm*(vd-E_Km_IN) : amp/cm2
Id_Na = phi_m_Na_IN*g_Na_INd*(DmNa**3)*DhNa*(vd-E_Na_IN) : amp/cm2
#Id_Nap = g_Nap_INd*DmNap*(vd-E_Na_IN) : amp/cm2

Is_Na = phi_m_Na_IN*g_Na_INs*(SmNa**3)*ShNa*(vs-E_Na_IN) : amp/cm2
Is_Kv = tadj_Kv_IN*g_Kv_INs*SmKv*(vs-E_Kv_IN) : amp/cm2
#Is_Nap = g_Nap_INs*SmNap*(vs-E_Na_IN) : amp/cm2

dwn_AMPA_PY_IN/dt = -(wn_AMPA_PY_IN)/wn_tau_ampa + (Amax_ampa_noise*w_AMPA_PY_IN)*(((2*(NSynSigma)**2)/wn_tau_ampa)**.5*xi) : siemens/cm2
g_temp1 = g_AMPA_PY_IN+wn_AMPA_PY_IN : siemens/cm2
I_AMPA_PY_IN = (g_temp1>0)*g_temp1*(vd-E_AMPA_PY_IN) : amp/cm2

dwn_NMDA_PY_IN/dt = -(wn_NMDA_PY_IN)/wn_tau_nmda + (Amax_nmda_noise*w_NMDA_PY_IN)*(((2*(NSynSigma)**2)/wn_tau_nmda)**.5*xi) : siemens/cm2
g_temp2 = g_NMDA_PY_IN+wn_NMDA_PY_IN : siemens/cm2
I_NMDA_PY_IN = (g_temp2>0)*g_temp2*(vd-E_NMDA_PY_IN) : amp/cm2

dwn_AMPA_EXT_IN/dt = -(wn_AMPA_EXT_IN)/wn_tau_ampa + (Amax_ampa_noise*w_AMPA_EXT_IN)*(((2*(NSynSigma)**2)/wn_tau_ampa)**.5*xi) : siemens/cm2
g_temp3 = g_AMPA_EXT_IN+wn_AMPA_EXT_IN : siemens/cm2
I_AMPA_EXT_IN = (g_temp3>0)*g_temp3*(vd-E_AMPA_EXT_IN) : amp/cm2

dwn_AMPA_TC_IN/dt = -(wn_AMPA_TC_IN)/wn_tau_ampa + (Amax_ampa_noise*w_AMPA_TC_IN)*(((2*(NSynSigma)**2)/wn_tau_ampa)**.5*xi) : siemens/cm2
g_temp4 = g_AMPA_TC_IN+wn_AMPA_TC_IN : siemens/cm2
I_AMPA_TC_IN = (g_temp4>0)*g_temp4*(vd-E_AMPA_TC_IN) : amp/cm2

LFP = abs(I_AMPA_PY_IN)+abs(I_NMDA_PY_IN)+abs(I_AMPA_TC_IN)+abs(I_AMPA_EXT_IN) : amp/cm2

# gurultusuz sinaptik akımlar
I_AMPA_PY_INx = g_AMPA_PY_INx*(vd-E_AMPA_PY_IN) : amp/cm2
I_NMDA_PY_INx = g_NMDA_PY_INx*(vd-E_NMDA_PY_IN) : amp/cm2
I_AMPA_TC_INx = g_AMPA_TC_INx*(vd-E_AMPA_TC_IN) : amp/cm2
I_AMPA_EXT_INx = g_AMPA_EXT_INx*(vd-E_AMPA_EXT_IN) : amp/cm2
LFPx = abs(I_AMPA_PY_INx)+abs(I_NMDA_PY_INx)+abs(I_AMPA_TC_INx)+abs(I_AMPA_EXT_INx) : amp/cm2
#----------------------------------------------------------------

# Dendritik hızlı sodyum Na akım kinetikleri
dDmNa/dt = -(DmNa - DmNa_inf)/tau_DmNa : 1 
dDhNa/dt = -(DhNa - DhNa_inf)/tau_DhNa : 1 

DmNa_inf = alphaDmNa/(alphaDmNa+betaDmNa) : 1
DhNa_inf = (1/(1+exp(((vd+Vshift_Na_IN)-thinf_Na_IN)/qinf_Na_IN))) : 1
tau_DhNa = (1/(alphaDhNa+betaDhNa))/phi_h_Na_IN : second
tau_DmNa = (1/(alphaDmNa+betaDmNa))/phi_m_Na_IN : second

alphaDmNa = 1.*(fabs((vd+Vshift_Na_IN)/(tha_Na_IN))>1e-6)*((Ra_Na_IN*(mV**-1)*\
            ((vd+Vshift_Na_IN)-(tha_Na_IN)))/(1-exp(-((vd+Vshift_Na_IN)-(tha_Na_IN))/\
            qa_Na_IN))) + 1.*(1-(fabs((vd+Vshift_Na_IN)/(tha_Na_IN))>1e-6))*(Ra_Na_IN*qa_Na_IN*(mV**-1)) : Hz
#alphaDmNa = ((Ra_Na_IN*(mV**-1)*((vd+Vshift_Na_IN)-(tha_Na_IN)))/(1-exp(-((vd+Vshift_Na_IN)-(tha_Na_IN))/\
#            qa_Na_IN))) : Hz

betaDmNa = 1.*(fabs(-(vd+Vshift_Na_IN)/(-tha_Na_IN))>1e-6)*((Rb_Na_IN*(mV**-1)*\
            (-(vd+Vshift_Na_IN)-(-tha_Na_IN)))/(1-exp(-(-(vd+Vshift_Na_IN)-(-tha_Na_IN))/\
            qa_Na_IN))) + 1.*(1-(fabs(-(vd+Vshift_Na_IN)/(-tha_Na_IN))>1e-6))*(Rb_Na_IN*qa_Na_IN*(mV**-1)) : Hz
#betaDmNa = ((Rb_Na_IN*(mV**-1)*(-(vd+Vshift_Na_IN)-(-tha_Na_IN)))/(1-exp(-(-(vd+Vshift_Na_IN)-(-tha_Na_IN))/\
#            qa_Na_IN))) : Hz

alphaDhNa = 1.*(fabs((vd+Vshift_Na_IN)/(thi1_Na_IN))>1e-6)*((Rd_Na_IN*(mV**-1)*\
            ((vd+Vshift_Na_IN)-(thi1_Na_IN)))/(1-exp(-((vd+Vshift_Na_IN)-(thi1_Na_IN))/\
            qi_Na_IN))) + 1.*(1-(fabs((vd+Vshift_Na_IN)/(thi1_Na_IN))>1e-6))*(Rd_Na_IN*qi_Na_IN*(mV**-1)) : Hz
#alphaDhNa = ((Rd_Na_IN*(mV**-1)*((vd+Vshift_Na_IN)-(thi1_Na_IN)))/(1-exp(-((vd+Vshift_Na_IN)-(thi1_Na_IN))/\
#            qi_Na_IN))) : Hz

betaDhNa = 1.*(fabs(-(vd+Vshift_Na_IN)/(-thi2_Na_IN))>1e-6)*((Rg_Na_IN*(mV**-1)*\
            (-(vd+Vshift_Na_IN)-(-thi2_Na_IN)))/(1-exp(-(-(vd+Vshift_Na_IN)-(-thi2_Na_IN))/\
            qi_Na_IN))) + 1.*(1-(fabs(-(vd+Vshift_Na_IN)/(-thi2_Na_IN))>1e-6))*(Rg_Na_IN*qi_Na_IN*(mV**-1)) : Hz
#betaDhNa = ((Rg_Na_IN*(mV**-1)*(-(vd+Vshift_Na_IN)-(-thi2_Na_IN)))/(1-exp(-(-(vd+Vshift_Na_IN)-(-thi2_Na_IN))/\
#            qi_Na_IN))) : Hz

# Aksosomatik hızlı sodyum Na akım kinetikleri
dSmNa/dt = -(SmNa - SmNa_inf)/tau_SmNa : 1 
dShNa/dt = -(ShNa - ShNa_inf)/tau_ShNa : 1 

SmNa_inf = alphaSmNa/(alphaSmNa+betaSmNa) : 1
ShNa_inf = (1/(1+exp(((vs+Vshift_Na_IN)-thinf_Na_IN)/qinf_Na_IN))) : 1
tau_SmNa = (1/(alphaSmNa+betaSmNa))/phi_m_Na_IN : second
tau_ShNa = (1/(alphaShNa+betaShNa))/phi_h_Na_IN : second

alphaSmNa = 1.*(fabs((vs+Vshift_Na_IN)/(tha_Na_IN))>1e-6)*((Ra_Na_IN*(mV**-1)*\
            ((vs+Vshift_Na_IN)-(tha_Na_IN)))/(1-exp(-((vs+Vshift_Na_IN)-(tha_Na_IN))/\
            qa_Na_IN))) + 1.*(1-(fabs((vs+Vshift_Na_IN)/(tha_Na_IN))>1e-6))*(Ra_Na_IN*qa_Na_IN*(mV**-1)) : Hz
#alphaSmNa = ((Ra_Na_IN*(mV**-1)*((vs+Vshift_Na_IN)-(tha_Na_IN)))/(1-exp(-((vs+Vshift_Na_IN)-(tha_Na_IN))/\
#            qa_Na_IN))) : Hz

betaSmNa = 1.*(fabs(-(vs+Vshift_Na_IN)/(-tha_Na_IN))>1e-6)*((Rb_Na_IN*(mV**-1)*\
            (-(vs+Vshift_Na_IN)-(-tha_Na_IN)))/(1-exp(-(-(vs+Vshift_Na_IN)-(-tha_Na_IN))/\
            qa_Na_IN))) + 1.*(1-(fabs(-(vs+Vshift_Na_IN)/(-tha_Na_IN))>1e-6))*(Rb_Na_IN*qa_Na_IN*(mV**-1)) : Hz
#betaSmNa = ((Rb_Na_IN*(mV**-1)*(-(vs+Vshift_Na_IN)-(-tha_Na_IN)))/(1-exp(-(-(vs+Vshift_Na_IN)-(-tha_Na_IN))/\
#            qa_Na_IN))) : Hz

alphaShNa = 1.*(fabs((vs+Vshift_Na_IN)/(thi1_Na_IN))>1e-6)*((Rd_Na_IN*(mV**-1)*\
            ((vs+Vshift_Na_IN)-(thi1_Na_IN)))/(1-exp(-((vs+Vshift_Na_IN)-(thi1_Na_IN))/\
            qi_Na_IN))) + 1.*(1-(fabs((vs+Vshift_Na_IN)/(thi1_Na_IN))>1e-6))*(Rd_Na_IN*qi_Na_IN*(mV**-1)) : Hz
#alphaShNa = ((Rd_Na_IN*(mV**-1)*((vs+Vshift_Na_IN)-(thi1_Na_IN)))/(1-exp(-((vs+Vshift_Na_IN)-(thi1_Na_IN))/\
#            qi_Na_IN))) : Hz

betaShNa = 1.*(fabs(-(vs+Vshift_Na_IN)/(-thi2_Na_IN))>1e-6)*((Rg_Na_IN*(mV**-1)*\
            (-(vs+Vshift_Na_IN)-(-thi2_Na_IN)))/(1-exp(-(-(vs+Vshift_Na_IN)-(-thi2_Na_IN))/\
            qi_Na_IN))) + 1.*(1-(fabs(-(vs+Vshift_Na_IN)/(-thi2_Na_IN))>1e-6))*(Rg_Na_IN*qi_Na_IN*(mV**-1)) : Hz
#betaShNa = ((Rg_Na_IN*(mV**-1)*(-(vs+Vshift_Na_IN)-(-thi2_Na_IN)))/(1-exp(-(-(vs+Vshift_Na_IN)-(-thi2_Na_IN))/\
#            qi_Na_IN))) : Hz
#----------------------------------------------------------------
 
# Aksosomatik hızlı potasyum Kv akım kinetikleri
dSmKv/dt = -(SmKv - SmKv_inf)/tau_SmKv : 1 
SmKv_inf = alphaSmKv/(alphaSmKv+betaSmKv) : 1
tau_SmKv = (1/(alphaSmKv+betaSmKv))/tadj_Kv_IN : second

alphaSmKv = alpham_Kv(vs,Ra_Kv_IN,tha_Kv_IN,qa_Kv_IN) : Hz
betaSmKv = betam_Kv(vs,Rb_Kv_IN,tha_Kv_IN,qa_Kv_IN) : Hz
#----------------------------------------------------------------
 
# Dendritik potasyum Km-akım kinetikleri
dDmKm/dt = -(DmKm - DmKm_inf)/tau_DmKm : 1 
DmKm_inf = alphaDmKm/(alphaDmKm+betaDmKm) : 1
tau_DmKm = (1/(alphaDmKm+betaDmKm))/tadj_Km_IN : second

alphaDmKm = alpham_Km(vd,Ra_Km_IN,tha_Km_IN,qa_Km_IN) : Hz
betaDmKm = betam_Km(vd,Rb_Km_IN,tha_Km_IN,qa_Km_IN) : Hz
#----------------------------------------------------------------
 
# Dendritik Ca-bağımlı potasyum KCa-akım kinetikleri
dDmKCa/dt = -(DmKCa-DmKCa_inf)/tau_DmKCa : 1
tau_DmKCa = (1/(alphaDmKCa+betaDmKCa))/tadj_KCa_IN : second
DmKCa_inf = (alphaDmKCa/(alphaDmKCa+betaDmKCa)) : 1

alphaDmKCa = (Ra_KCa_IN*(DCa_i/mMolar)**caix_KCa_IN) : Hz
betaDmKCa = Rb_KCa_IN : Hz
#----------------------------------------------------------------
 
# Dendritik Hight-threshold Ca2+ CaL-akım kinetikleri
dDmCaL/dt = -(DmCaL-DmCaL_inf)/tau_DmCaL : 1
dDhCaL/dt = -(DhCaL-DhCaL_inf)/tau_DhCaL : 1

alphaDmCaL = 0.055*(mV**-1)*(-27*mV-vd)/ \
    (exp((-27*mV-vd)/(3.8*mV))-1.)/ms : Hz
betaDmCaL = 0.94*exp((-75*mV-vd)/(17*mV))/ms : Hz
tau_DmCaL = (1/(alphaDmCaL+betaDmCaL))/phi_m_CaL_IN : second
DmCaL_inf = (alphaDmCaL/(alphaDmCaL+betaDmCaL)) : 1

alphaDhCaL = 0.000457*exp((-13*mV-vd)/(50*mV))/ms : Hz
betaDhCaL = 0.0065/(1+exp((-vd-15*mV)/(28*mV)))/ms : Hz
tau_DhCaL = (1/(alphaDhCaL+betaDhCaL))/phi_h_CaL_IN : second
DhCaL_inf = (alphaDhCaL/(alphaDhCaL+betaDhCaL)) : 1
#----------------------------------------------------------------
 
# dendrit Ca_i dinamiği
#dDCa_i/dt = -(A*(Id_CaL))+(DCa_i_inf_IN-DCa_i)/(taur_cad_IN) : Molar #mole/cm3
dDCa_i/dt = -(Id_CaL)/(2*F*depth_cad_IN)+(DCa_i_inf_IN-DCa_i)/(taur_cad_IN) : Molar #mole/cm3
#dDCa_i/dt = -(1e4*I_T/mA)/(2*96489*1*ms)+((2.4e-4)-DCa_i)/(taur_cad_IN) : 1
#----------------------------------------------------------------
 
# dendrit Kalsiyum reversal potansiyel dinamiği
dEDCa_IN/dt = (Nernst(valence,Temperature, DCa_i,DCa_o_IN)-EDCa_IN)/(ms) : volt
#----------------------------------------------------------------
 
g_AMPA_PY_IN : siemens/cm2
g_NMDA_PY_IN : siemens/cm2
g_AMPA_TC_IN : siemens/cm2
g_AMPA_EXT_IN : siemens/cm2

#gürültüsüz sinaptik akımlar
g_AMPA_PY_INx : siemens/cm2
g_NMDA_PY_INx : siemens/cm2
g_AMPA_TC_INx : siemens/cm2
g_AMPA_EXT_INx : siemens/cm2

Id_ext : amp/cm2
Is_ext : amp/cm2
 
''')
#----------------------------------------------------------------
 
# TC nöron modeli
# Bazhenov ve ark. (2002;1998a)
eqs_TC = Equations('''
dv/dt = (I-I_syn-I_int)/Cm_TC : volt
 
I_int = I_L+I_KL+I_Na+I_K+I_T+I_h+I_A : amp/cm2

dwn/dt = -(wn)/wn_tau_TC + (((2*(Noise_sigma_tot_TC)**2)/wn_tau_TC)**.5*xi) : amp/cm2

I_syn = (I_GABAA_RE_TC+I_GABAB_RE_TC+I_AMPA_EXT_TC+I_AMPA_PY_TC)+wn : amp/cm2
 
dwn_GABAA_RE_TC/dt = -(wn_GABAA_RE_TC)/wn_tau_gabaa + (Amax_gabaa_noise*w_GABAA_RE_TC)*(((2*(NSynSigma)**2)/wn_tau_gabaa)**.5*xi) : siemens/cm2
g_temp1 = g_GABAA_RE_TC+wn_GABAA_RE_TC : siemens/cm2
I_GABAA_RE_TC = (g_temp1>0)*g_temp1*(v-E_GABAA_RE_TC) : amp/cm2

dwn_GABAB_RE_TC/dt = -(wn_GABAB_RE_TC)/wn_tau_gabab + (Amax_gabab_noise*w_GABAB_RE_TC)*(((2*(NSynSigma)**2)/wn_tau_gabab)**.5*xi) : siemens/cm2
g_temp2 = g_GABAB_RE_TC+wn_GABAB_RE_TC : siemens/cm2
I_GABAB_RE_TC = (g_temp2>0)*g_temp2*(v-E_GABAB_RE_TC) : amp/cm2

dwn_AMPA_EXT_TC/dt = -(wn_AMPA_EXT_TC)/wn_tau_ampa + (Amax_ampa_noise*w_AMPA_EXT_TC)*(((2*(NSynSigma)**2)/wn_tau_ampa)**.5*xi) : siemens/cm2
g_temp3 = g_AMPA_EXT_TC+wn_AMPA_EXT_TC : siemens/cm2
I_AMPA_EXT_TC = (g_temp3>0)*g_temp3*(v-E_AMPA_EXT_TC) : amp/cm2

dwn_AMPA_PY_TC/dt = -(wn_AMPA_PY_TC)/wn_tau_ampa + (Amax_ampa_noise*w_AMPA_PY_TC)*(((2*(NSynSigma)**2)/wn_tau_ampa)**.5*xi) : siemens/cm2
g_temp4 = g_AMPA_PY_TC+wn_AMPA_PY_TC : siemens/cm2
I_AMPA_PY_TC = (g_temp4>0)*g_temp4*(v-E_AMPA_PY_TC) : amp/cm2

LFP = abs(I_GABAA_RE_TC)+abs(I_GABAB_RE_TC)+abs(I_AMPA_PY_TC)+abs(I_AMPA_EXT_TC) : amp/cm2

# Gürültüsüz sinaptik akımlar
I_GABAA_RE_TCx = g_GABAA_RE_TCx*(v-E_GABAA_RE_TC) : amp/cm2
I_GABAB_RE_TCx = g_GABAB_RE_TCx*(v-E_GABAB_RE_TC) : amp/cm2
I_AMPA_PY_TCx = g_AMPA_PY_TCx*(v-E_AMPA_PY_TC) : amp/cm2
I_AMPA_EXT_TCx = g_AMPA_EXT_TCx*(v-E_AMPA_EXT_TC) : amp/cm2
LFPx = abs(I_GABAA_RE_TCx)+abs(I_GABAB_RE_TCx)+abs(I_AMPA_PY_TCx)+abs(I_AMPA_EXT_TCx) : amp/cm2

I_L = g_L_TC*(v-EL_TC) : amp/cm2
I_KL = gg_KL_TC*(v-EKL_TC) : amp/cm2
I_Na = g_Na_TC*(m*m*m)*h*(v-ENa_TC) : amp/cm2
I_K = g_K_TC*(n*n*n*n)*(v-EK_TC) : amp/cm2
I_A = g_A_TC*(mAK**4)*hAK*(v-EK_TC) : amp/cm2
#----------------------------------------------------------------
 
# A-akım kinetikleri 
dmAK/dt=(mAK_inf-mAK)/tau_mAK : 1
mAK_inf=mAKinf_TC(v) : 1   
tau_mAK=taumAK_TC(v) : second
 
dhAK/dt=(hAK_inf-hAK)/tau_hAK : 1
hAK_inf=hAKinf_TC(v) : 1   
tau_hAK=tauhAK_TC(v) : second
#----------------------------------------------------------------
 
# I_h akım kinetik ifadeleri (Destexhe ve ark. 1996a)
I_h = gg_h_TC*mh*(v-Eh_TC) : amp/cm2
mh = (o1+ginc*o2) :1
dc1/dt = (betamh*(1-c1-o2) - alphamh*c1)/ms : 1
dp1/dt = (k1ca*(1-p1) - k2*p1)/ms : 1
do2/dt = (k3p*(1-c1-o2) - k4*o2)/ms : 1
k1ca = (k2 * ((Ca_i/cac)**nca)) : 1
k3p = (k4 * (p1/Pc)**nexp) : 1
p0 = (1-p1) : 1
#c1 = (1-o1) : 1
o1 = (1-c1-o2) : 1
betamh = (1-mh_inf)/tau_mh : 1
alphamh = mh_inf/tau_mh : 1
mh_inf = 1./(1+exp((v+75*mV)/(5.5*mV))) : 1
tau_mh = (taumin + (taumax)/(exp((v+71.5*mV)/(14.2*mV))+exp(-(v+89*mV)/(11.6*mV))))/hadj_TC : 1
#----------------------------------------------------------------
 
# Ca_i dinamiği
#dCa_i/dt = -(A_TC*(I_T))+(Ca_i_inf_TC-Ca_i)/(taur_cad_TC) : Molar
dCa_i/dt = -(I_T)/(2*F*depth_cad_TC)+(Ca_i_inf_TC-Ca_i)/(taur_cad_TC) : Molar
#dCa_i/dt = -mMolar*(1e4*I_T/(mA*2.9e-4))/(2*96489*1*ms)+(Ca_i_inf_TC-Ca_i)/(taur_cad_TC) : Molar
#----------------------------------------------------------------
 
# Kalsiyum reversal potansiyel dinamiği
dECa_TC/dt = (Nernst(valence,Temperature, Ca_i,Ca_o_TC)-ECa_TC)/(ms) : volt
#ECa_TC = (Nernst(valence,Temperature, Ca_i,Ca_o_TC)) : volt
#ECa_TC = R * (Temperature) / (F * valence) * log(Ca_o_TC / Ca_i) : volt
#----------------------------------------------------------------
 
# I_T akım ve kinetikleri kinetikleri
I_T = g_T_TC*(mCa**2)*hCa*(v-ECa_TC) : amp/cm2
#dmCa/dt=(mCa_inf-mCa)/tau_mCa : 1
dmCa/dt=(mCa_inf-mCa)/ms : 1
mCa_inf=mCainf_TC(v) : 1   
#tau_mCa=taumCa_TC(v) : second
dhCa/dt=(hCa_inf-hCa)/tau_hCa : 1
hCa_inf=hCainf_TC(v) : 1   
tau_hCa=tauhCa_TC(v) : second
#----------------------------------------------------------------
 
# fast sodyum ve potasyum akım kinetikleri
dm/dt = alpham*(1-m)-betam*m : 1
dn/dt = alphan*(1-n)-betan*n : 1
dh/dt = alphah*(1-h)-betah*h : 1
 
alpham = 0.32*(mV**-1)*(13*mV-v+VTr_TC)/ \
    (exp((13*mV-v+VTr_TC)/(4*mV))-1.)/ms : Hz
 
betam = 0.28*(mV**-1)*(v-VTr_TC-40*mV)/ \
    (exp((v-VTr_TC-40*mV)/(5*mV))-1)/ms : Hz
 
alphah = 0.128*exp((17*mV-v+VTr_TC)/(18*mV))/ms : Hz
betah = 4./(1+exp((40*mV-v+VTr_TC)/(5*mV)))/ms : Hz
 
alphan = 0.032*(mV**-1)*(15*mV-v+VTrK_TC)/ \
    (exp((15*mV-v+VTrK_TC)/(5*mV))-1.)/ms : Hz
betan = .5*exp((10*mV-v+VTrK_TC)/(40*mV))/ms : Hz
#----------------------------------------------------------------
 
g_GABAA_RE_TC : siemens/cm2
g_GABAB_RE_TC : siemens/cm2
g_AMPA_PY_TC : siemens/cm2
g_AMPA_EXT_TC : siemens/cm2

#gürültüsüz sinaptik akımlar
g_GABAA_RE_TCx : siemens/cm2
g_GABAB_RE_TCx : siemens/cm2
g_AMPA_PY_TCx : siemens/cm2
g_AMPA_EXT_TCx : siemens/cm2

gg_KL_TC : siemens/cm2
gg_h_TC : siemens/cm2
 
I : amp/cm2
 
''')
#----------------------------------------------------------------
 
# RE nöron modeli
eqs_RE = Equations('''
dv/dt = (I-I_syn-I_int)/Cm_RE : volt
 
I_int = I_L+I_KL+I_Na+I_K+I_Ts : amp/cm2

dwn/dt = -(wn)/wn_tau_RE + (((2*(Noise_sigma_tot_RE)**2)/wn_tau_RE)**.5*xi) : amp/cm2

I_syn = (I_AMPA_TC_RE+I_GABAA_RE_RE+I_AMPA_EXT_RE+I_AMPA_PY_RE)+wn : amp/cm2
 
dwn_AMPA_TC_RE/dt = -(wn_AMPA_TC_RE)/wn_tau_ampa + (Amax_ampa_noise*w_AMPA_TC_RE)*(((2*(NSynSigma)**2)/wn_tau_ampa)**.5*xi) : siemens/cm2
g_temp1 = g_AMPA_TC_RE+wn_AMPA_TC_RE : siemens/cm2
I_AMPA_TC_RE = (g_temp1>0)*g_temp1*(v-E_AMPA_TC_RE) : amp/cm2

dwn_GABAA_RE_RE/dt = -(wn_GABAA_RE_RE)/wn_tau_gabaa + (Amax_gabaa_noise*w_GABAA_RE_RE)*(((2*(NSynSigma)**2)/wn_tau_gabaa)**.5*xi) : siemens/cm2
g_temp2 = g_GABAA_RE_RE+wn_GABAA_RE_RE : siemens/cm2
I_GABAA_RE_RE = (g_temp2>0)*g_temp2*(v-E_GABAA_RE_RE) : amp/cm2

dwn_AMPA_EXT_RE/dt = -(wn_AMPA_EXT_RE)/wn_tau_ampa + (Amax_ampa_noise*w_AMPA_EXT_RE)*(((2*(NSynSigma)**2)/wn_tau_ampa)**.5*xi) : siemens/cm2
g_temp3 = g_AMPA_EXT_RE+wn_AMPA_EXT_RE : siemens/cm2
I_AMPA_EXT_RE = (g_temp3>0)*g_temp3*(v-E_AMPA_EXT_RE) : amp/cm2

dwn_AMPA_PY_RE/dt = -(wn_AMPA_PY_RE)/wn_tau_ampa + (Amax_ampa_noise*w_AMPA_PY_RE)*(((2*(NSynSigma)**2)/wn_tau_ampa)**.5*xi) : siemens/cm2
g_temp4 = g_AMPA_PY_RE+wn_AMPA_PY_RE : siemens/cm2
I_AMPA_PY_RE = (g_temp4>0)*g_temp4*(v-E_AMPA_PY_RE) : amp/cm2

LFP = abs(I_AMPA_TC_RE)+abs(I_GABAA_RE_RE)+abs(I_AMPA_PY_RE)+abs(I_AMPA_EXT_RE) : amp/cm2

# Gürültüsüz sinaptik akımlar
I_AMPA_TC_REx = g_AMPA_TC_REx*(v-E_AMPA_TC_RE) : amp/cm2
I_GABAA_RE_REx = g_GABAA_RE_REx*(v-E_GABAA_RE_RE) : amp/cm2
I_AMPA_PY_REx = g_AMPA_PY_REx*(v-E_AMPA_PY_RE) : amp/cm2
I_AMPA_EXT_REx = g_AMPA_EXT_REx*(v-E_AMPA_EXT_RE) : amp/cm2
LFPx = abs(I_AMPA_TC_REx)+abs(I_GABAA_RE_REx)+abs(I_AMPA_PY_REx)+abs(I_AMPA_EXT_REx) : amp/cm2

I_L = g_L_RE*(v-EL_RE) : amp/cm2
I_KL = gg_KL_RE*(v-EKL_RE) : amp/cm2
I_Na = g_Na_RE*(m*m*m)*h*(v-ENa_RE) : amp/cm2
I_K = g_K_RE*(n*n*n*n)*(v-EK_RE) : amp/cm2
#----------------------------------------------------------------
 
# Ca_i dinamiği
#dCa_i/dt = -(A_RE*I_Ts)+(Ca_i_inf_RE-Ca_i)/(taur_cad_RE) : Molar
dCa_i/dt = -(I_Ts)/(2*F*depth_cad_RE)+(Ca_i_inf_RE-Ca_i)/(taur_cad_RE) : Molar
#dCa_i/dt = -(1e4*I_Ts/mA)/(2*96489*1*ms)+((2.4e-4)-Ca_i)/(taur_cad_RE) : 1
#----------------------------------------------------------------
 
# Kalsiyum reversal potansiyel dinamiği
dECa_RE/dt = (Nernst(valence,Temperature, Ca_i,Ca_o_RE)-ECa_RE)/(ms) : volt
#----------------------------------------------------------------
 
# I_Ts akım kinetikleri
I_Ts = g_Ts_RE*(mCa**2)*hCa*(v-ECa_RE) : amp/cm2
dmCa/dt=(mCa_inf-mCa)/tau_mCa : 1
mCa_inf=mCainf_RE(v) : 1   
tau_mCa=taumCa_RE(v) : second
dhCa/dt=(hCa_inf-hCa)/tau_hCa : 1
hCa_inf=hCainf_RE(v) : 1   
tau_hCa=tauhCa_RE(v) : second
#----------------------------------------------------------------
 
# fast sodyum ve potasyum akım kinetikleri
dm/dt = alpham*(1-m)-betam*m : 1
dn/dt = alphan*(1-n)-betan*n : 1
dh/dt = alphah*(1-h)-betah*h : 1
 
alpham = 0.32*(mV**-1)*(13*mV-v+VTr_RE)/ \
    (exp((13*mV-v+VTr_RE)/(4*mV))-1.)/ms : Hz
 
betam = 0.28*(mV**-1)*(v-VTr_RE-40*mV)/ \
    (exp((v-VTr_RE-40*mV)/(5*mV))-1)/ms : Hz
 
alphah = 0.128*exp((17*mV-v+VTr_RE)/(18*mV))/ms : Hz
betah = 4./(1+exp((40*mV-v+VTr_RE)/(5*mV)))/ms : Hz
 
alphan = 0.032*(mV**-1)*(15*mV-v+VTrK_RE)/ \
    (exp((15*mV-v+VTrK_RE)/(5*mV))-1.)/ms : Hz
betan = .5*exp((10*mV-v+VTrK_RE)/(40*mV))/ms : Hz
#----------------------------------------------------------------
 
g_AMPA_TC_RE : siemens/cm2
g_GABAA_RE_RE : siemens/cm2
g_AMPA_PY_RE : siemens/cm2
g_AMPA_EXT_RE : siemens/cm2

# Gürültüsüz sinaptik iletkenlikler
g_AMPA_TC_REx : siemens/cm2
g_GABAA_RE_REx : siemens/cm2
g_AMPA_PY_REx : siemens/cm2
g_AMPA_EXT_REx : siemens/cm2

gg_KL_RE : siemens/cm2
 
I : amp/cm2
 
''')
#----------------------------------------------------------------
 
####################### Nöral gruplar #######################################
 
#Kortikal grup
PY = NeuronGroup(N_PY, model=eqs_PY, threshold=Threshold(0 * mV,state='vs'), refractory=3*ms, #threshold= 0 * mV,
    #threshold=EmpiricalThreshold(threshold= 0 * mV, refractory=1 * ms, clock=myclock),
    implicit=True, freeze=True, 
    #method='exponential_Euler', 
    method='linear', 
    clock=myclock)
 
IN = NeuronGroup(N_IN, model=eqs_IN, threshold=Threshold(0 * mV,state='vs'), refractory=2*ms, #threshold= 0 * mV,
    #threshold=EmpiricalThreshold(threshold= 0 * mV, refractory=1 * ms, clock=myclock),
    implicit=True, freeze=True, 
    #method='exponential_Euler', 
    method='linear', 
    clock=myclock)
 
#Talamik grup
TC = NeuronGroup(N_TC, model=eqs_TC, threshold= 0 * mV,refractory=3*ms,
    #threshold=EmpiricalThreshold(threshold= 0 * mV, refractory=1 * ms, clock=myclock),
    implicit=True, freeze=True, 
    #method='exponential_Euler', 
    method='linear', 
    clock=myclock)
 
RE = NeuronGroup(N_RE, model=eqs_RE, threshold= 0 * mV,refractory=2*ms,
    #threshold=EmpiricalThreshold(threshold= 0 * mV, refractory=1 * ms, clock=myclock),
    implicit=True, freeze=True, 
    #method='exponential_Euler', 
    method='linear', 
    clock=myclock)
#----------------------------------------------------------------
 
print "Bağlantılar oluşturuluyor..."
 
######################### TC spike giriş ###################################
# TC hücre için spike girişi
"""
Nspikes = 10;T=15*ms;T0=T*Nspikes;
TC_input = SpikeGeneratorGroup(1,[(0,n*T+T0) for n in range(Nspikes)])
C_TC_input=Connection(TC_input,TC,'g_AMPA_TC_TC',weight=0.5*uS)
"""
# PoissonInput birbirinden bağımsız giriş sayısı kadar poison girişi state değişkenine
# weight ağırlığa göre ekler
# deta için http://briansimulator.org/docs/inputs.html?highlight=poissoninput
# C_TC_input = PoissonInput(TC, bgs, rate=20*Hz, weight=0.5*uS, state='g_AMPA_TC_TC')
# poison dağılıma göre 
#TC_input = PoissonGroup(1,rates=20*Hz)
# Gaussian dağılıma göre
#TC_input=PulsePacket(t=100*ms,n=10,sigma=30*ms)
#C_TC_input=Connection(TC_input,TC,'g_AMPA_TC_TC',weight=0.5*uS)
# elle ayarlanmış düzenli spike dizisi
 
#spiketimes = [(0,1100*ms), (0,1200*ms), (0,1300*ms), (0,1400*ms), (0,1500*ms),\
#              (0,1600*ms), (0,1700*ms), (0,1800*ms), (0,1900*ms)]

#f_ext = 10 # Hz
#dispersion_rate = 0.2; # spike'ların frekans etrafında dağılı oranı

# Kortikal harici giriş
spiketimes_PY = spikegen(N_PY,f_ext,dispersion_rate,duration/second)
PY_input = SpikeGeneratorGroup(N_PY,spiketimes_PY, clock=myclock)

#spiketimes = [(0,(x+1)*(1./f_ext)*second) for x in range(f_ext)]
#PY_input = SpikeGeneratorGroup(1,spiketimes, clock=myclock)

spiketimes_IN = spikegen(N_IN,f_ext,dispersion_rate,duration/second)
IN_input = SpikeGeneratorGroup(N_IN,spiketimes_IN, clock=myclock)

PY_one_input = SpikeGeneratorGroup(1,[(0,100*ms)], clock=myclock)
IN_one_input = SpikeGeneratorGroup(1,[(0,100*ms)], clock=myclock)
PY_poisson_input=PoissonGroup(N_PY,rates=linspace(0*Hz,10*Hz,N_PY), clock=myclock)

# Talamik harici giriş
spiketimes_TC = spikegen(N_TC,f_ext,dispersion_rate,duration/second)
TC_input = SpikeGeneratorGroup(N_TC,spiketimes_TC, clock=myclock)

#spiketimes = [(0,(x+1)*(1./f_ext)*second) for x in range(f_ext)]
#TC_input = SpikeGeneratorGroup(1,spiketimes, clock=myclock)
spiketimes_RE = spikegen(N_RE,f_ext,dispersion_rate,duration/second)
RE_input = SpikeGeneratorGroup(N_RE,spiketimes_RE, clock=myclock)

RE_one_input = SpikeGeneratorGroup(1,[(0,100*ms)], clock=myclock)
TC_one_input = SpikeGeneratorGroup(1,[(0,100*ms)], clock=myclock)
TC_poisson_input=PoissonGroup(N_TC,rates=linspace(0*Hz,10*Hz,N_TC), clock=myclock)

#################### Small-world parametreleri ########################

#randCon_rate = 0*p_val2[1];
#----------------------------------------------------------------
 
#----------------------------------------------------------------
#################### sinaptik gürültü parametreleri ########################

#new=False

#----------------------------------------------------------------
 
############################ AMPA Sinaps ####################################
# AMPA sinaptik sabitler
Cmax_AMPA   = ((1-new)*0.5+new*1.0) * mMolar            #  maks. trasmiter konsantrasyonu
Cdur_AMPA   = ((1-new)*0.3+new*1.0) * ms                #  trasmiter süresi (rising phase)
Alpha_AMPA  = ((1-new)*0.94+new*1.1) / (ms * mMolar)    #  bağlanma (binding) oranı
Beta_AMPA   = ((1-new)*0.18+new*0.19) / ms              #  serbestlik (unbinding) oranı

Rinf_AMPA = Cmax_AMPA*Alpha_AMPA / (Cmax_AMPA*Alpha_AMPA + Beta_AMPA)
Rtau_AMPA = 1. / ((Alpha_AMPA * Cmax_AMPA) + Beta_AMPA)
 
# AMPA sinaps modeli : 
# Destexhe ve ark. 1998b - 
# Bazhenov ve ark. 1998a Cellular and network models ...
 
AMPA_PY_PY = Synapses(PY, PY,
             model='''w : siemens/cm2
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      #T = (Cmax_AMPA)/(1+exp(-(vd_pre-2*mV)/(5*mV))) : Molar #mole/cm3
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_AMPA,dt)*Cmax_AMPA : Molar #mole/cm3
                      #dx_ampa/dt = heaviside3(t_post,lr1,Cdur_AMPA,dt)*(Rinf_AMPA - x_ampa)/(Rtau_AMPA) - (1-heaviside3(t_post,lr1,Cdur_AMPA,dt))*x_ampa*Beta_AMPA : 1
                      dwn/dt = -(wn-wn_mu)/wn_tau_ampa + (((2*(Noise_sigma)**2)/wn_tau_ampa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_ampa_noise : 1
                      gampa1 = w*(x_ampa + wn*Amax_ampa_noise) : siemens/cm2
                      gampa = (gampa1>0)*gampa1 : siemens/cm2
                      gampax = w*x_ampa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear'
                      method='exponential_Euler'
                      )
 
AMPA_TC_PY = Synapses(TC, PY,
             model='''w : siemens/cm2
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      #T = (Cmax_AMPA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : Molar #mole/cm3
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_AMPA,dt)*Cmax_AMPA : Molar #mole/cm3
                      #dx_ampa/dt = heaviside3(t_post,lr1,Cdur_AMPA,dt)*(Rinf_AMPA - x_ampa)/(Rtau_AMPA) - (1-heaviside3(t_post,lr1,Cdur_AMPA,dt))*x_ampa*Beta_AMPA : 1
                      dwn/dt = -(wn-wn_mu)/wn_tau_ampa + (((2*(Noise_sigma)**2)/wn_tau_ampa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_ampa_noise : 1
                      gampa1 = w*(x_ampa + wn*Amax_ampa_noise) : siemens/cm2
                      gampa = (gampa1>0)*gampa1 : siemens/cm2
                      gampax = w*x_ampa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear'
                      method='exponential_Euler'
                      )

AMPA_PY_IN = Synapses(PY, IN,
             model='''w : siemens/cm2
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      #T = (Cmax_AMPA)/(1+exp(-(vd_pre-2*mV)/(5*mV))) : Molar #mole/cm3
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_AMPA,dt)*Cmax_AMPA : Molar #mole/cm3
                      #dx_ampa/dt = heaviside3(t_post,lr1,Cdur_AMPA,dt)*(Rinf_AMPA - x_ampa)/(Rtau_AMPA) - (1-heaviside3(t_post,lr1,Cdur_AMPA,dt))*x_ampa*Beta_AMPA : 1
                      dwn/dt = -(wn-wn_mu)/wn_tau_ampa + (((2*(Noise_sigma)**2)/wn_tau_ampa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_ampa_noise : 1
                      gampa1 = w*(x_ampa + wn*Amax_ampa_noise) : siemens/cm2
                      gampa = (gampa1>0)*gampa1 : siemens/cm2
                      gampax = w*x_ampa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear'
                      method='exponential_Euler'
                      ) 

AMPA_TC_IN = Synapses(TC, IN,
             model='''w : siemens/cm2
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      #T = (Cmax_AMPA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : Molar #mole/cm3
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_AMPA,dt)*Cmax_AMPA : Molar #mole/cm3
                      #dx_ampa/dt = heaviside3(t_post,lr1,Cdur_AMPA,dt)*(Rinf_AMPA - x_ampa)/(Rtau_AMPA) - (1-heaviside3(t_post,lr1,Cdur_AMPA,dt))*x_ampa*Beta_AMPA : 1
                      dwn/dt = -(wn-wn_mu)/wn_tau_ampa + (((2*(Noise_sigma)**2)/wn_tau_ampa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_ampa_noise : 1
                      gampa1 = w*(x_ampa + wn*Amax_ampa_noise) : siemens/cm2
                      gampa = (gampa1>0)*gampa1 : siemens/cm2
                      gampax = w*x_ampa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear'
                      method='exponential_Euler'
                      )

#----------------------------------------------------------------
 
AMPA_PY_TC = Synapses(PY, TC,
             model='''w : siemens/cm2
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      #T = (Cmax_AMPA)/(1+exp(-(vd_pre-2*mV)/(5*mV))) : Molar #mole/cm3
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_AMPA,dt)*Cmax_AMPA : Molar #mole/cm3
                      #dx_ampa/dt = heaviside3(t_post,lr1,Cdur_AMPA,dt)*(Rinf_AMPA - x_ampa)/(Rtau_AMPA) - (1-heaviside3(t_post,lr1,Cdur_AMPA,dt))*x_ampa*Beta_AMPA : 1
                      dwn/dt = -(wn-wn_mu)/wn_tau_ampa + (((2*(Noise_sigma)**2)/wn_tau_ampa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_ampa_noise : 1
                      gampa1 = w*(x_ampa + wn*Amax_ampa_noise) : siemens/cm2
                      gampa = (gampa1>0)*gampa1 : siemens/cm2
                      gampax = w*x_ampa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear'
                      method='exponential_Euler'
                      )
 
AMPA_TC_RE = Synapses(TC, RE,
             model='''w : siemens/cm2
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      #T = (Cmax_AMPA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : Molar #mole/cm3
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_AMPA,dt)*Cmax_AMPA : Molar #mole/cm3
                      #dx_ampa/dt = heaviside3(t_post,lr1,Cdur_AMPA,dt)*(Rinf_AMPA - x_ampa)/(Rtau_AMPA) - (1-heaviside3(t_post,lr1,Cdur_AMPA,dt))*x_ampa*Beta_AMPA : 1
                      dwn/dt = -(wn-wn_mu)/wn_tau_ampa + (((2*(Noise_sigma)**2)/wn_tau_ampa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_ampa_noise : 1
                      gampa1 = w*(x_ampa + wn*Amax_ampa_noise) : siemens/cm2
                      gampa = (gampa1>0)*gampa1 : siemens/cm2
                      gampax = w*x_ampa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear',
                      method='exponential_Euler',
                      )
 
AMPA_PY_RE = Synapses(PY, RE,
             model='''w : siemens/cm2
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      #T = (Cmax_AMPA)/(1+exp(-(vd_pre-2*mV)/(5*mV))) : Molar #mole/cm3
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_AMPA,dt)*Cmax_AMPA : Molar #mole/cm3
                      #dx_ampa/dt = heaviside3(t_post,lr1,Cdur_AMPA,dt)*(Rinf_AMPA - x_ampa)/(Rtau_AMPA) - (1-heaviside3(t_post,lr1,Cdur_AMPA,dt))*x_ampa*Beta_AMPA : 1
                      dwn/dt = -(wn-wn_mu)/wn_tau_ampa + (((2*(Noise_sigma)**2)/wn_tau_ampa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_ampa_noise : 1
                      gampa1 = w*(x_ampa + wn*Amax_ampa_noise) : siemens/cm2
                      gampa = (gampa1>0)*gampa1 : siemens/cm2
                      gampax = w*x_ampa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear'
                      method='exponential_Euler'
                      )

AMPA_EXT_PY = Synapses(PY_input, PY,
             model='''w : siemens/cm2
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      #T = (Cmax_AMPA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : Molar #mole/cm3
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_AMPA,dt)*Cmax_AMPA : Molar #mole/cm3
                      #dx_ampa/dt = heaviside3(t_post,lr1,Cdur_AMPA,dt)*(Rinf_AMPA - x_ampa)/(Rtau_AMPA) - (1-heaviside3(t_post,lr1,Cdur_AMPA,dt))*x_ampa*Beta_AMPA : 1
                      dwn/dt = -(wn-wn_mu)/wn_tau_ampa + (((2*(Noise_sigma)**2)/wn_tau_ampa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_ampa_noise : 1
                      gampa1 = w*(x_ampa + 0*wn*Amax_ampa_noise) : siemens/cm2
                      gampa = (gampa1>0)*gampa1 : siemens/cm2
                      gampax = w*x_ampa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear'
                      method='exponential_Euler'
                      )
#----------------------------------------------------------------

AMPA_EXT_IN = Synapses(IN_input, IN,
             model='''w : siemens/cm2
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      #T = (Cmax_AMPA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : Molar #mole/cm3
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_AMPA,dt)*Cmax_AMPA : Molar #mole/cm3
                      #dx_ampa/dt = heaviside3(t_post,lr1,Cdur_AMPA,dt)*(Rinf_AMPA - x_ampa)/(Rtau_AMPA) - (1-heaviside3(t_post,lr1,Cdur_AMPA,dt))*x_ampa*Beta_AMPA : 1
                      dwn/dt = -(wn-wn_mu)/wn_tau_ampa + (((2*(Noise_sigma)**2)/wn_tau_ampa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_ampa_noise : 1
                      gampa1 = w*(x_ampa + 0*wn*Amax_ampa_noise) : siemens/cm2
                      gampa = (gampa1>0)*gampa1 : siemens/cm2
                      gampax = w*x_ampa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear'
                      method='exponential_Euler'
                      )

AMPA_EXT_TC = Synapses(TC_input, TC,
             model='''w : siemens/cm2
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      #T = (Tmax_AMPA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : Molar
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_AMPA,dt)*Cmax_AMPA : Molar
                      dwn/dt = -(wn-wn_mu)/wn_tau_ampa + (((2*(Noise_sigma)**2)/wn_tau_ampa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_ampa_noise : 1
                      gampa1 = w*(x_ampa + 0*wn*Amax_ampa_noise) : siemens/cm2
                      gampa = (gampa1>0)*gampa1 : siemens/cm2
                      gampax = w*x_ampa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear',
                      method='exponential_Euler',
                      )

AMPA_EXT_RE = Synapses(RE_input, RE,
             model='''w : siemens/cm2
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      #T = (Tmax_AMPA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : Molar
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_AMPA,dt)*Cmax_AMPA : Molar
                      dwn/dt = -(wn-wn_mu)/wn_tau_ampa + (((2*(Noise_sigma)**2)/wn_tau_ampa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_ampa_noise : 1
                      gampa1 = w*(x_ampa + 0*wn*Amax_ampa_noise) : siemens/cm2
                      gampa = (gampa1>0)*gampa1 : siemens/cm2
                      gampax = w*x_ampa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear',
                      method='exponential_Euler',
                      )

#----------------------------------------------------------------

############################ NMDA Sinaps ############################
# NMDA sinaptik sabitler
Cmax_NMDA   = ((1-new)*0.5+new*1.0) * mMolar        #  maks. trasmiter konsantrasyonu
Cdur_NMDA   = ((1-new)*0.3+new*1.0) * ms            #  trasmiter süresi (rising phase)
Alpha_NMDA  = ((1-new)*1.0+new*0.072) / (ms*mMolar) #  bağlanma (binding) oranı
Beta_NMDA   = ((1-new)*0.0067+new*0.0066) / ms      #  serbestlik (unbinding) oranı

mg_NMDA = 0*2 * mMolar                 # [Mg]_o = 1 - 2 mM

Rinf = Cmax_NMDA*Alpha_NMDA / (Cmax_NMDA*Alpha_NMDA + Beta_NMDA)
Rtau = 1. / ((Alpha_NMDA * Cmax_NMDA) + Beta_NMDA)
 
 
# NMDA sinaps modeli : 
# Destexhe ve ark. 1998b - 
# Bazhenov ve ark. 1998a Cellular and network models ...

NMDA_PY_PY = Synapses(PY, PY,
             model='''w : siemens/cm2
                      dx_nmda/dt = (Alpha_NMDA*T)*(1-x_nmda)-Beta_NMDA*x_nmda : 1
                      #T = (Cmax_NMDA)/(1+exp(-(vd_pre-2*mV)/(5*mV))) : Molar #mole/cm3
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_NMDA,dt)*Cmax_NMDA : Molar #mole/cm3
                      #dx_nmda/dt = heaviside3(t_post,lr1,Cdur_NMDA,dt)*(Rinf - x_nmda)/(Rtau) - (1-heaviside3(t_post,lr1,Cdur_NMDA,dt))*x_nmda*Beta_NMDA : 1                      
                      mgblock = (1.0 / (1.0 + exp(-0.062 * (vd_post/mV)) * (mg_NMDA / (3.57 * mMolar)))) : 1
                      #mgblock = 1/(1+exp(-(vd_post - (-25*mV))/(12.5*mV))) : 1        
                      dwn/dt = -(wn-wn_mu)/wn_tau_nmda + (((2*(Noise_sigma)**2)/wn_tau_nmda)**.5*xi) : 1 #siemens/cm2
                      r_nmda = x_nmda*mgblock : 1
                      #wnn = (wn>0)*wn * Amax_nmda_noise : 1
                      gnmda1 = w*(r_nmda + wn*Amax_nmda_noise) : siemens/cm2
                      gnmda = (gnmda1>0)*gnmda1 : siemens/cm2
                      gnmdax = w*r_nmda : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t',
                      method='exponential_Euler'
                      #method='linear'
                      )
 
NMDA_PY_IN = Synapses(PY, IN,
             model='''w : siemens/cm2
                      dx_nmda/dt = (Alpha_NMDA*T)*(1-x_nmda)-Beta_NMDA*x_nmda : 1
                      #T = (Cmax_NMDA)/(1+exp(-(vd_pre-2*mV)/(5*mV))) : Molar #mole/cm3
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_NMDA,dt)*Cmax_NMDA : Molar #mole/cm3
                      #dx_nmda/dt = heaviside3(t_post,lr1,Cdur_NMDA,dt)*(Rinf - x_nmda)/(Rtau) - (1-heaviside3(t_post,lr1,Cdur_NMDA,dt))*x_nmda*Beta_NMDA : 1                      
                      mgblock = (1.0 / (1.0 + exp(-0.062 * (vd_post/mV)) * (mg_NMDA / (3.57 * mMolar)))) : 1
                      #mgblock = 1/(1+exp(-(vd_post - (-25*mV))/(12.5*mV))) : 1        
                      dwn/dt = -(wn-wn_mu)/wn_tau_nmda + (((2*(Noise_sigma)**2)/wn_tau_nmda)**.5*xi) : 1 #siemens/cm2
                      r_nmda = x_nmda*mgblock : 1
                      #wnn = (wn>0)*wn * Amax_nmda_noise : 1
                      gnmda1 = w*(r_nmda + wn*Amax_nmda_noise) : siemens/cm2
                      gnmda = (gnmda1>0)*gnmda1 : siemens/cm2
                      gnmdax = w*r_nmda : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t',
                      method='exponential_Euler'
                      #method='linear'
                      )
 
#----------------------------------------------------------------
 
############################ GABAA Sinaps ############################
# GABAA sinaptik sabitler
Cmax_GABAA   = ((1-new)*0.5+new*1.0) * mMolar           #  maks. trasmiter konsantrasyonu
Cdur_GABAA   = ((1-new)*0.3+new*1.0) * ms               #  trasmiter süresi (rising phase)
Alpha_GABAA  = ((1-new)*10.5+new*5.0) / (ms * mMolar)   #  bağlanma (binding) oranı
Beta_GABAA   = ((1-new)*0.166+new*0.18) / ms            #  serbestlik (unbinding) oranı
 
Rinf_GABAA = Cmax_GABAA*Alpha_GABAA / (Cmax_GABAA*Alpha_GABAA + Beta_GABAA)
Rtau_GABAA = 1. / ((Alpha_GABAA * Cmax_GABAA) + Beta_GABAA)
 
# GABAA sinaps modeli : 
# Destexhe ve ark. 1998b - 
# Bazhenov ve ark. 1998a Cellular and network models ...
GABAA_IN_PY = Synapses(IN, PY,
             model='''w : siemens/cm2
                      dx_gabaa/dt = (Alpha_GABAA*T)*(1-x_gabaa)-Beta_GABAA*x_gabaa : 1
                      #T = (Cmax_GABAA)/(1+exp(-(vd_pre-2*mV)/(5*mV))) : Molar
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_GABAA,dt)*Cmax_GABAA : Molar
                      #dx_gabaa/dt = heaviside3(t_post,lr1,Cdur_GABAA,dt)*(Rinf_GABAA - x_gabaa)/(Rtau_GABAA) - (1-heaviside3(t_post,lr1,Cdur_GABAA,dt))*x_gabaa*Beta_GABAA : 1                                            
                      dwn/dt = -(wn-wn_mu)/wn_tau_gabaa + (((2*(Noise_sigma)**2)/wn_tau_gabaa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_gabaa_noise : 1
                      ggabaa1 = w*(x_gabaa + wn*Amax_gabaa_noise) : siemens/cm2
                      ggabaa = (ggabaa1>0)*ggabaa1 : siemens/cm2
                      ggabaax = w*x_gabaa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear'
                      method='exponential_Euler'
                      )
#----------------------------------------------------------------
                       
GABAA_RE_RE = Synapses(RE, RE,
              model='''w : siemens/cm2
                      dx_gabaa/dt = (Alpha_GABAA*T)*(1-x_gabaa)-Beta_GABAA*x_gabaa : 1
                      #T = (Cmax_GABAA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : Molar
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_GABAA,dt)*Cmax_GABAA : Molar
                      #dx_gabaa/dt = heaviside3(t_post,lr1,Cdur_GABAA,dt)*(Rinf_GABAA - x_gabaa)/(Rtau_GABAA) - (1-heaviside3(t_post,lr1,Cdur_GABAA,dt))*x_gabaa*Beta_GABAA : 1                                            
                      dwn/dt = -(wn-wn_mu)/wn_tau_gabaa + (((2*(Noise_sigma)**2)/wn_tau_gabaa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_gabaa_noise : 1
                      ggabaa1 = w*(x_gabaa + wn*Amax_gabaa_noise) : siemens/cm2
                      ggabaa = (ggabaa1>0)*ggabaa1 : siemens/cm2
                      ggabaax = w*x_gabaa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear', 
                      method='exponential_Euler'
                      )

GABAA_RE_TC = Synapses(RE, TC,
              model='''w : siemens/cm2
                      dx_gabaa/dt = (Alpha_GABAA*T)*(1-x_gabaa)-Beta_GABAA*x_gabaa : 1
                      #T = (Cmax_GABAA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : Molar
                      T = heaviside3(t_post+(1-new)*dt,lr1,Cdur_GABAA,dt)*Cmax_GABAA : Molar
                      #dx_gabaa/dt = heaviside3(t_post,lr1,Cdur_GABAA,dt)*(Rinf_GABAA - x_gabaa)/(Rtau_GABAA) - (1-heaviside3(t_post,lr1,Cdur_GABAA,dt))*x_gabaa*Beta_GABAA : 1                                            
                      dwn/dt = -(wn-wn_mu)/wn_tau_gabaa + (((2*(Noise_sigma)**2)/wn_tau_gabaa)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_gabaa_noise : 1
                      ggabaa1 = w*(x_gabaa + wn*Amax_gabaa_noise) : siemens/cm2
                      ggabaa = (ggabaa1>0)*ggabaa1 : siemens/cm2
                      ggabaax = w*x_gabaa : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear', 
                      method='exponential_Euler'
                      )

#----------------------------------------------------------------

############################ GABAB Sinaps ############################
# GABAB sinaptik sabitler
Cmax_GABAB   = ((1-new)*0.5+new*1.0) * mMolar           #  maks. trasmiter konsantrasyonu
Cdur_GABAB   = ((1-new)*0.3+new*1.0) * ms               #  trasmiter süresi (rising phase)

K1_GABAB    = ((1-new)*0.52+new*0.09) / (ms * mMolar)   #  reseptörlere bağlanma oranı
K2_GABAB    = ((1-new)*0.0013+new*0.0012) / ms          #  reseptörlerin serbest kalma oranı
K3_GABAB    = ((1-new)*0.098+new*0.18) / ms             #  G-protein oluşma oranı
K4_GABAB    = ((1-new)*0.033+new*0.034) / ms            #  G-protein yok olma oranı
KD_GABAB    = 100                                       #  K+ kanalı ayrışma sabiti (uMolar**4)

n_GABAB     = 4             # G-proteininin K+ kanalına bağlanma nokta sayısı
 
# AMPA sinaps modeli : 
# Destexhe ve ark. 1998b - 
# Bazhenov ve ark. 1998a Cellular and network models ...
GABAB_RE_TC = Synapses(RE, TC,
             model='''w : siemens/cm2
                      dx_gabab/dt = (K1_GABAB*T)*(1-x_gabab)-K2_GABAB*x_gabab : 1
                      #T = (Cmax_GABAB)/(1+exp(-(v_pre-2*mV)/(5*mV))) : Molar
                      T = heaviside3(t_post+0*(1-new)*dt,lr1,Cdur_GABAB+0*new*dt,dt)*Cmax_GABAB : Molar
                      ds_gabab/dt = (K3_GABAB*x_gabab)-(K4_GABAB*s_gabab) : 1
                      dwn/dt = -(wn-wn_mu)/wn_tau_gabab + (((2*(Noise_sigma)**2)/wn_tau_gabab)**.5*xi) : 1 #siemens/cm2
                      #wnn = (wn>0)*wn * Amax_gabab_noise : 1
                      r_gabab = (s_gabab**n_GABAB)/((s_gabab**n_GABAB)+KD_GABAB) : 1
                      ggabab1 = w*(r_gabab + wn*Amax_gabab_noise) : siemens/cm2
                      ggabab = (ggabab1>0)*ggabab1 : siemens/cm2
                      ggababx = w*r_gabab : siemens/cm2
                      lr1 : second
                      ''', pre='lr1 = t', 
                      #method='linear', # linear olduğunda x_gabab'i bir dt eksik hesaplıyor
                      method='exponential_Euler'
                      )
 
####################### Sinaptik bağlantılar ##############################
""""
# full bağlantı yapısı
# PY->PY AMPA bağ
PY.g_AMPA_PY_PY = AMPA_PY_PY.gampa
AMPA_PY_PY[:,:]=True
AMPA_PY_PY.w=w_AMPA_PY_PY/N_PY
 
# PY->IN AMPA bağ
IN.g_AMPA_PY_IN = AMPA_PY_IN.gampa
AMPA_PY_IN[:,:]=True
AMPA_PY_IN.w=w_AMPA_PY_IN/N_IN
 
# PY->PY NMDA bağ
PY.g_NMDA_PY_PY = NMDA_PY_PY.gnmda
NMDA_PY_PY[:,:]=True
NMDA_PY_PY.w=w_NMDA_PY_PY/N_PY
 
# PY->IN NMDA bağ
IN.g_NMDA_PY_IN = NMDA_PY_IN.gnmda
NMDA_PY_IN[:,:]=True
NMDA_PY_IN.w=w_NMDA_PY_IN/N_IN
 
# IN->PY GABAA bağ
PY.g_GABAA_IN_PY = GABAA_IN_PY.ggabaa
GABAA_IN_PY[:,:]=True
GABAA_IN_PY.w=w_GABAA_IN_PY/N_PY
#----------------------------------------------------------------
 
# RE->TC GABAA bağ
TC.g_GABAA_RE_TC = GABAA_RE_TC.ggabaa
GABAA_RE_TC[:,:]=True
GABAA_RE_TC.w=w_GABAA_RE_TC/N_RE
 
# RE->TC GABAB bağ
TC.g_GABAB_RE_TC = GABAB_RE_TC.ggabab
GABAB_RE_TC[:,:]=True
GABAB_RE_TC.w=w_GABAB_RE_TC/N_RE
 
# RE->RE GABAA bağ
RE.g_GABAA_RE_RE = GABAA_RE_RE.ggabaa
GABAA_RE_RE[:,:]=True #'i!=j'#True
GABAA_RE_RE.w=w_GABAA_RE_RE/N_RE
 
# TC->RE AMPA bağ
RE.g_AMPA_TC_RE = AMPA_TC_RE.gampa
AMPA_TC_RE[:,:]=True
AMPA_TC_RE.w=w_AMPA_TC_RE/N_TC
#----------------------------------------------------------------
# PY->TC AMPA bağ
TC.g_AMPA_PY_TC = AMPA_PY_TC.gampa
AMPA_PY_TC[:,:]=True
AMPA_PY_TC.w=w_AMPA_PY_TC/N_PY
 
# PY->RE AMPA bağ
RE.g_AMPA_PY_RE = AMPA_PY_RE.gampa
AMPA_PY_RE[:,:]=True
AMPA_PY_RE.w=w_AMPA_PY_RE/N_PY
 
# TC->PY AMPA bağ
PY.g_AMPA_TC_PY = AMPA_TC_PY.gampa
AMPA_TC_PY[:,:]=True
AMPA_TC_PY.w=w_AMPA_TC_PY/N_TC
 
# TC->IN AMPA bağ
IN.g_AMPA_TC_IN = AMPA_TC_IN.gampa
AMPA_TC_IN[:,:]=True
AMPA_TC_IN.w=w_AMPA_TC_IN/N_TC
"""
################ ÖNEMLİ NOT #######################
# small-world rewiring işlemi sırasımda bazı nöronların 
# belirli sinaptik türüne ait hiçbir bağlantısı kalmayabilir.
# bu durumda simülasyon hata veriyor (IndexError: index out of bounds). 
# tekrar çalıştırıldığında bu durumun oluşmadığı rastsallıkla sorunsuz 
# çalısırsa daha uygun olur. her nöronun her sinaptik bağ türünden  
# bağlantısı olması daha uygundur.
#----------------------------------------------------------------
 
# Kortikal sinaptik bağlantı indexleri
randCon_AMPA_PY_PY = randCon_rate;randCon_AMPA_PY_IN = randCon_rate;
randCon_NMDA_PY_PY = randCon_rate;randCon_NMDA_PY_IN = randCon_rate;
randCon_GABAA_IN_PY = randCon_rate;
 
index_AMPA_PY_PY = arange(0,N_PY); index_AMPA_PY_IN = arange(0,N_PY);
index_NMDA_PY_PY = arange(0,N_PY); index_NMDA_PY_IN = arange(0,N_PY);
index_GABAA_IN_PY = arange(0,N_IN);

bindex_AMPA_PY_PY = (ones((1,N_PY),bool)); bindex_AMPA_PY_IN = (ones((1,N_PY),bool))
bindex_NMDA_PY_PY = (ones((1,N_PY),bool)); bindex_NMDA_PY_IN = (ones((1,N_PY),bool))
bindex_GABAA_IN_PY = (ones((1,N_IN),bool));

bindex_AMPA_PY_PY[0,arange(kom_AMPA_PY_PY+1,N_PY-kom_AMPA_PY_PY)] = False
bindex_NMDA_PY_PY[0,arange(kom_NMDA_PY_PY+1,N_PY-kom_NMDA_PY_PY)] = False
#bindex_NMDA_PY_PY[0,0]=False
bindex_AMPA_PY_IN[0,arange(kom_AMPA_PY_IN+1,N_PY-kom_AMPA_PY_IN)] = False
bindex_NMDA_PY_IN[0,arange(kom_NMDA_PY_IN+1,N_PY-kom_NMDA_PY_IN)] = False
bindex_GABAA_IN_PY[0,arange(kom_GABAA_IN_PY+1,N_IN-kom_GABAA_IN_PY)] = False
#----------------------------------------------------------------
 
# Talamik sinaptik bağlantı indexleri
randCon_GABAA_RE_TC = randCon_rate;randCon_GABAA_RE_RE = randCon_rate;
randCon_GABAB_RE_TC = randCon_rate;randCon_AMPA_TC_RE = randCon_rate;
 
index_GABAA_RE_TC = arange(0,N_RE); index_GABAB_RE_TC = arange(0,N_RE);
index_GABAA_RE_RE = arange(0,N_RE); index_AMPA_TC_RE = arange(0,N_TC);
 
bindex_GABAA_RE_TC = (ones((1,N_RE),bool)); bindex_GABAB_RE_TC = (ones((1,N_RE),bool))
bindex_GABAA_RE_RE = (ones((1,N_RE),bool)); bindex_AMPA_TC_RE = (ones((1,N_TC),bool))
 
bindex_GABAA_RE_TC[0,arange(kom_GABAA_RE_TC+1,N_RE-kom_GABAA_RE_TC)] = False
bindex_GABAA_RE_RE[0,arange(kom_GABAA_RE_RE+1,N_RE-kom_GABAA_RE_RE)] = False
#bindex_GABAA_RE_RE[0,0]=False
bindex_GABAB_RE_TC[0,arange(kom_GABAB_RE_TC+1,N_RE-kom_GABAB_RE_TC)] = False
bindex_AMPA_TC_RE[0,arange(kom_AMPA_TC_RE+1,N_TC-kom_AMPA_TC_RE)] = False
#----------------------------------------------------------------
 
# TalamoKortikal sinaptik bağlantı indexleri
randCon_AMPA_TC_PY = randCon_rate;randCon_AMPA_TC_IN = randCon_rate;
randCon_AMPA_PY_TC = randCon_rate;randCon_AMPA_PY_RE = randCon_rate;
 
index_AMPA_TC_PY = arange(0,N_TC); index_AMPA_TC_IN = arange(0,N_TC);
index_AMPA_PY_TC = arange(0,N_PY); index_AMPA_PY_RE = arange(0,N_PY);
 
bindex_AMPA_TC_PY = (ones((1,N_TC),bool)); bindex_AMPA_TC_IN = (ones((1,N_TC),bool))
bindex_AMPA_PY_TC = (ones((1,N_PY),bool)); bindex_AMPA_PY_RE = (ones((1,N_PY),bool))
 
bindex_AMPA_TC_PY[0,arange(kom_AMPA_TC_PY+1,N_TC-kom_AMPA_TC_PY)] = False
bindex_AMPA_TC_IN[0,arange(kom_AMPA_TC_IN+1,N_TC-kom_AMPA_TC_IN)] = False
bindex_AMPA_PY_TC[0,arange(kom_AMPA_PY_TC+1,N_PY-kom_AMPA_PY_TC)] = False
bindex_AMPA_PY_RE[0,arange(kom_AMPA_PY_RE+1,N_PY-kom_AMPA_PY_RE)] = False
#----------------------------------------------------------------
 
# PY->PY AMPA bağlantılar 
PY.g_AMPA_PY_PY = AMPA_PY_PY.gampa; PY.g_AMPA_PY_PYx = AMPA_PY_PY.gampax; indxs_AMPA_PY_PY=[]
for i in range(N_PY):
    indx = array(shift_right(list(bindex_AMPA_PY_PY[0]),(i*N_PY/N_PY)))
    indxs_AMPA_PY_PY.append(indx)

aindx = array(indxs_AMPA_PY_PY).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_AMPA_PY_PY))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_AMPA_PY_PY=(aindx.reshape(N_PY,N_PY)).tolist()
 
for i in range(N_PY):
    indexs = [x for x in range(len(indxs_AMPA_PY_PY[i])) if indxs_AMPA_PY_PY[i][x]==True] 
    if len(indexs): AMPA_PY_PY[index_AMPA_PY_PY[indexs],i]=True
AMPA_PY_PY.w=w_AMPA_PY_PY
#----------------------------------------------------------------
 
# PY->IN AMPA bağ 
IN.g_AMPA_PY_IN = AMPA_PY_IN.gampa; IN.g_AMPA_PY_INx = AMPA_PY_IN.gampax; indxs_AMPA_PY_IN=[]
for i in range(N_IN):
    indx = array(shift_right(list(bindex_AMPA_PY_IN[0]),(i*N_PY/N_IN)))
    indxs_AMPA_PY_IN.append(indx)
 
aindx = array(indxs_AMPA_PY_IN).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_AMPA_PY_IN))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_AMPA_PY_IN=(aindx.reshape(N_IN,N_PY)).tolist()
 
for i in range(N_IN):
    indexs = [x for x in range(len(indxs_AMPA_PY_IN[i])) if indxs_AMPA_PY_IN[i][x]==True] 
    if len(indexs): AMPA_PY_IN[index_AMPA_PY_IN[indexs],i]=True 
AMPA_PY_IN.w=w_AMPA_PY_IN
#----------------------------------------------------------------
 
# PY->PY NMDA bağ 
PY.g_NMDA_PY_PY = NMDA_PY_PY.gnmda; PY.g_NMDA_PY_PYx = NMDA_PY_PY.gnmdax; indxs_NMDA_PY_PY=[]
for i in range(N_PY):
    indx = array(shift_right(list(bindex_NMDA_PY_PY[0]),(i*N_PY/N_PY)))
    indxs_NMDA_PY_PY.append(indx)
 
aindx = array(indxs_NMDA_PY_PY).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_NMDA_PY_PY))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False] # if ((aindx[x]==False) & (x%(N_RE+1)!=0))]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_NMDA_PY_PY=(aindx.reshape(N_PY,N_PY)).tolist()
 
for i in range(N_PY):
    indexs = [x for x in range(len(indxs_NMDA_PY_PY[i])) if indxs_NMDA_PY_PY[i][x]==True] 
    if len(indexs): NMDA_PY_PY[index_NMDA_PY_PY[indexs],i]=True
NMDA_PY_PY.w=w_NMDA_PY_PY
#----------------------------------------------------------------
 
# PY->IN NMDA bağ 
IN.g_NMDA_PY_IN = NMDA_PY_IN.gnmda; IN.g_NMDA_PY_INx = NMDA_PY_IN.gnmdax; indxs_NMDA_PY_IN=[]
for i in range(N_IN):
    indx = array(shift_right(list(bindex_NMDA_PY_IN[0]),(i*N_PY/N_IN)))
    indxs_NMDA_PY_IN.append(indx)
 
aindx = array(indxs_NMDA_PY_IN).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_NMDA_PY_IN))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_NMDA_PY_IN=(aindx.reshape(N_IN,N_PY)).tolist()
 
for i in range(N_IN):
    indexs = [x for x in range(len(indxs_NMDA_PY_IN[i])) if indxs_NMDA_PY_IN[i][x]==True] 
    if len(indexs): NMDA_PY_IN[index_NMDA_PY_IN[indexs],i]=True 
NMDA_PY_IN.w=w_NMDA_PY_IN
#----------------------------------------------------------------
 
# IN->PY GABAA bağ 
PY.g_GABAA_IN_PY = GABAA_IN_PY.ggabaa; PY.g_GABAA_IN_PYx = GABAA_IN_PY.ggabaax; indxs_GABAA_IN_PY=[]
for i in range(N_PY):
    indx = array(shift_right(list(bindex_GABAA_IN_PY[0]),(i*N_IN/N_PY)))
    indxs_GABAA_IN_PY.append(indx)

aindx = array(indxs_GABAA_IN_PY).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_GABAA_IN_PY))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_GABAA_IN_PY=(aindx.reshape(N_PY,N_IN)).tolist()
 
for i in range(N_PY):
    indexs = [x for x in range(len(indxs_GABAA_IN_PY[i])) if indxs_GABAA_IN_PY[i][x]==True]
    if len(indexs): GABAA_IN_PY[index_GABAA_IN_PY[indexs],i]=True
GABAA_IN_PY.w=w_GABAA_IN_PY
#----------------------------------------------------------------
 
# RE->TC GABAA bağlantılar 
TC.g_GABAA_RE_TC = GABAA_RE_TC.ggabaa; TC.g_GABAA_RE_TCx = GABAA_RE_TC.ggabaax; indxs_GABAA_RE_TC=[]
for i in range(N_TC):
    indx = array(shift_right(list(bindex_GABAA_RE_TC[0]),(i*N_RE/N_TC)))
    indxs_GABAA_RE_TC.append(indx)
 
aindx = array(indxs_GABAA_RE_TC).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_GABAA_RE_TC))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_GABAA_RE_TC=(aindx.reshape(N_TC,N_RE)).tolist()
 
for i in range(N_TC):
    indexs = [x for x in range(len(indxs_GABAA_RE_TC[i])) if indxs_GABAA_RE_TC[i][x]==True] 
    if len(indexs): GABAA_RE_TC[index_GABAA_RE_TC[indexs],i]=True
GABAA_RE_TC.w=w_GABAA_RE_TC
#----------------------------------------------------------------

# RE->TC GABAB bağ 
TC.g_GABAB_RE_TC = GABAB_RE_TC.ggabab; TC.g_GABAB_RE_TCx = GABAB_RE_TC.ggababx; indxs_GABAB_RE_TC=[]
for i in range(N_TC):
    indx = array(shift_right(list(bindex_GABAB_RE_TC[0]),(i*N_RE/N_TC)))
    indxs_GABAB_RE_TC.append(indx)
 
aindx = array(indxs_GABAB_RE_TC).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_GABAB_RE_TC))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_GABAB_RE_TC=(aindx.reshape(N_TC,N_RE)).tolist()
 
for i in range(N_TC):
    indexs = [x for x in range(len(indxs_GABAB_RE_TC[i])) if indxs_GABAB_RE_TC[i][x]==True]
    if len(indexs): GABAB_RE_TC[index_GABAB_RE_TC[indexs],i]=True
GABAB_RE_TC.w=w_GABAB_RE_TC
#----------------------------------------------------------------

# RE->RE GABAA bağ 
RE.g_GABAA_RE_RE = GABAA_RE_RE.ggabaa; RE.g_GABAA_RE_REx = GABAA_RE_RE.ggabaax; indxs_GABAA_RE_RE=[]
for i in range(N_RE):
    indx = array(shift_right(list(bindex_GABAA_RE_RE[0]),(i*N_RE/N_RE)))
    indxs_GABAA_RE_RE.append(indx)
 
aindx = array(indxs_GABAA_RE_RE).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_GABAA_RE_RE))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False] # if ((aindx[x]==False) & (x%(N_RE+1)!=0))]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_GABAA_RE_RE=(aindx.reshape(N_RE,N_RE)).tolist()
 
for i in range(N_RE):
    indexs = [x for x in range(len(indxs_GABAA_RE_RE[i])) if indxs_GABAA_RE_RE[i][x]==True] 
    if len(indexs): GABAA_RE_RE[index_GABAA_RE_RE[indexs],i]=True
GABAA_RE_RE.w=w_GABAA_RE_RE
#----------------------------------------------------------------

# TC->RE AMPA bağ 
RE.g_AMPA_TC_RE = AMPA_TC_RE.gampa; RE.g_AMPA_TC_REx = AMPA_TC_RE.gampax; indxs_AMPA_TC_RE=[]
for i in range(N_RE):
    indx = array(shift_right(list(bindex_AMPA_TC_RE[0]),(i*N_TC/N_RE)))
    indxs_AMPA_TC_RE.append(indx)
     
aindx = array(indxs_AMPA_TC_RE).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_AMPA_TC_RE))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_AMPA_TC_RE=(aindx.reshape(N_RE,N_TC)).tolist()
 
for i in range(N_RE):
    indexs = [x for x in range(len(indxs_AMPA_TC_RE[i])) if indxs_AMPA_TC_RE[i][x]==True]
    if len(indexs): AMPA_TC_RE[index_AMPA_TC_RE[indexs],i]=True
AMPA_TC_RE.w=w_AMPA_TC_RE
#----------------------------------------------------------------

# PY->TC AMPA bağlantılar 
TC.g_AMPA_PY_TC = AMPA_PY_TC.gampa; TC.g_AMPA_PY_TCx = AMPA_PY_TC.gampax; indxs_AMPA_PY_TC=[]
for i in range(N_TC):
    indx = array(shift_right(list(bindex_AMPA_PY_TC[0]),(i*N_PY/N_TC)))
    indxs_AMPA_PY_TC.append(indx)
 
aindx = array(indxs_AMPA_PY_TC).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_AMPA_PY_TC))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_AMPA_PY_TC=(aindx.reshape(N_TC,N_PY)).tolist()
 
for i in range(N_TC):
    indexs = [x for x in range(len(indxs_AMPA_PY_TC[i])) if indxs_AMPA_PY_TC[i][x]==True] 
    if len(indexs): AMPA_PY_TC[index_AMPA_PY_TC[indexs],i]=True
AMPA_PY_TC.w=w_AMPA_PY_TC
#----------------------------------------------------------------
 
# PY->RE AMPA bağlantılar 
RE.g_AMPA_PY_RE = AMPA_PY_RE.gampa; RE.g_AMPA_PY_REx = AMPA_PY_RE.gampax; indxs_AMPA_PY_RE=[]
for i in range(N_RE):
    indx = array(shift_right(list(bindex_AMPA_PY_RE[0]),(i*N_PY/N_RE)))
    indxs_AMPA_PY_RE.append(indx)
 
aindx = array(indxs_AMPA_PY_RE).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_AMPA_PY_RE))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_AMPA_PY_RE=(aindx.reshape(N_RE,N_PY)).tolist()
 
for i in range(N_RE):
    indexs = [x for x in range(len(indxs_AMPA_PY_RE[i])) if indxs_AMPA_PY_RE[i][x]==True] 
    if len(indexs): AMPA_PY_RE[index_AMPA_PY_RE[indexs],i]=True 
AMPA_PY_RE.w=w_AMPA_PY_RE
#----------------------------------------------------------------
 
# TC->PY AMPA bağlantılar 
PY.g_AMPA_TC_PY = AMPA_TC_PY.gampa; PY.g_AMPA_TC_PYx = AMPA_TC_PY.gampax; indxs_AMPA_TC_PY=[]
for i in range(N_PY):
    indx = array(shift_right(list(bindex_AMPA_TC_PY[0]),(i*N_TC/N_PY)))
    indxs_AMPA_TC_PY.append(indx)
 
aindx = array(indxs_AMPA_TC_PY).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_AMPA_TC_PY))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_AMPA_TC_PY=(aindx.reshape(N_PY,N_TC)).tolist()
 
for i in range(N_PY):
    indexs = [x for x in range(len(indxs_AMPA_TC_PY[i])) if indxs_AMPA_TC_PY[i][x]==True] 
    if len(indexs): AMPA_TC_PY[index_AMPA_TC_PY[indexs],i]=True
AMPA_TC_PY.w=w_AMPA_TC_PY
#----------------------------------------------------------------
 
# TC->IN AMPA bağlantılar 
IN.g_AMPA_TC_IN = AMPA_TC_IN.gampa; IN.g_AMPA_TC_INx = AMPA_TC_IN.gampax; indxs_AMPA_TC_IN=[]
for i in range(N_IN):
    indx = array(shift_right(list(bindex_AMPA_TC_IN[0]),(i*N_TC/N_IN)))
    indxs_AMPA_TC_IN.append(indx)
 
aindx = array(indxs_AMPA_TC_IN).flatten('C')
connects = [x for x in range(len(aindx)) if aindx[x]==True]
rand_connects = rnd.sample(connects, int(len(connects)*randCon_AMPA_TC_IN))
aindx[rand_connects] = ~aindx[rand_connects]
unconnects = [x for x in range(len(aindx)) if aindx[x]==False]
randindx = rnd.sample(unconnects,len(rand_connects))
aindx[randindx] = ~aindx[randindx]
indxs_AMPA_TC_IN=(aindx.reshape(N_IN,N_TC)).tolist()

for i in range(N_IN):
    indexs = [x for x in range(len(indxs_AMPA_TC_IN[i])) if indxs_AMPA_TC_IN[i][x]==True] 
    if len(indexs): AMPA_TC_IN[index_AMPA_TC_IN[indexs],i]=True 
AMPA_TC_IN.w=w_AMPA_TC_IN
#----------------------------------------------------------------

#################### harici hirişler #####################
#Kortikal
#w_AMPA_EXT_PY_rnd = np.random.uniform(0.25,1.0,N_PY)
PY.g_AMPA_EXT_PY = AMPA_EXT_PY.gampa; PY.g_AMPA_EXT_PYx = AMPA_EXT_PY.gampax
AMPA_EXT_PY[:,:]='i==j' # giriş spike matriksini hedef ağa ile bire bir eşleştiriyoruz
AMPA_EXT_PY.w = w_AMPA_EXT_PY*w_AMPA_EXT_PY_rnd # 'exp(-abs(j-i-N_PY/2)*.1)*w_AMPA_EXT_PY' # 

IN.g_AMPA_EXT_IN = AMPA_EXT_IN.gampa; IN.g_AMPA_EXT_INx = AMPA_EXT_IN.gampax
AMPA_EXT_IN[:,:]='i==j' #True
AMPA_EXT_IN.w = w_AMPA_EXT_IN*w_AMPA_EXT_IN_rnd #'0*exp(-abs(j-i-N_PY/2)*.2)*w_AMPA_EXT_IN' #
#----------------------------------------------------------------
#Talamik
TC.g_AMPA_EXT_TC = AMPA_EXT_TC.gampa; TC.g_AMPA_EXT_TCx = AMPA_EXT_TC.gampax;
AMPA_EXT_TC[:,:]='i==j' # giriş spike matriksini hedef ağa ile bire bir eşleştiriyoruz
AMPA_EXT_TC.w = w_AMPA_EXT_TC*w_AMPA_EXT_TC_rnd # 'exp(-abs(j-i-N_TC/2)*.1)*w_AMPA_EXT_TC' # 

RE.g_AMPA_EXT_RE = AMPA_EXT_RE.gampa; RE.g_AMPA_EXT_REx = AMPA_EXT_RE.gampax;
AMPA_EXT_RE[:,:]='i==j' #True
AMPA_EXT_RE.w = w_AMPA_EXT_RE*w_AMPA_EXT_RE_rnd #'0*exp(-abs(j-i-N_RE/2)*.2)*w_AMPA_EXT_RE' #

################## Simulas ön hazırlıklar ##############################
print "Gerekli hazirliklar yapılıyor..."
 
pref_str='TalamoKortikalNet_'+'sig_'+('%.3f'%(Noise_sigma)).replace('.','_')+ \
    '_p_'+('%.2f'%(randCon_rate)).replace('.','_')+'_k_'+('%.0f'%(kom_rate)).replace('.','_')+ \
    '_dr_'+('%.2f'%(dispersion_rate)).replace('.','_')+'_t_'+repr(int(duration/ms))
 
#----------------------------------------------------------------
 
# PY nöron hazırlık
#PY.Is_ext=linked_var(input_PYs,'B')
#PY.Id_ext=linked_var(input_PYd,'B')
 
PY.DCa_i = 1e-4 * mMolar;PY.SCa_i = 1e-4 * mMolar
PY.EDCa_PY = 120 * mV; PY.ESCa_PY = 120 * mV
PY.vd = -70*mV # np.random.normal(-70,0.1,N_PY)*mV 
PY.vs = -70*mV # np.random.normal(-70,0.1,N_PY)*mV 
"""
x0 = INaCX_init(-70.); PY.DmNa = x0[0]; PY.DhNa = x0[1]; #PY.SmNa = x0[0]; PY.ShNa = x0[1];
x0 = ICaL_init(-70.,phi_m_CaL_PY,phi_h_CaL_PY); PY.DmCaL = x0[0];PY.DhCaL = x0[1];
x0 = IKCa_init(1e-4); PY.DmKCa = x0;
x0 = IKm_init(-70); PY.DmKm = x0;
x0 = INap_init(-70,phi_m_Nap_PY); PY.DmNap = x0; PY.SmNap = x0;
"""
# IN nöron hazırlık
#IN.Is_ext=linked_var(input_INs,'B')
 
IN.DCa_i = 1e-4 * mMolar #(1e-6 * mole * cm**-3) # 1 M = 1e-3 mole/cm3 
IN.EDCa_IN = 120 * mV
IN.vd = -68*mV # np.random.normal(-68,0.1,N_IN)*mV 
IN.vs = -68*mV # np.random.normal(-68,0.1,N_IN)*mV 
"""
x0 = INaCX_init(-68.); IN.DmNa = x0[0]; IN.DhNa = x0[1]; #IN.SmNa = x0[0]; IN.ShNa = x0[1];
x0 = ICaL_init(-68.,phi_m_CaL_IN,phi_h_CaL_IN); IN.DmCaL = x0[0]; IN.DhCaL = x0[1];
x0 = IKCa_init(1e-4); IN.DmKCa = x0;
x0 = IKm_init(-68); IN.DmKm = x0;
"""
#----------------------------------------------------------------
 
# TC nöron hazırlık
g_h_TC_vals = np.random.uniform(0.25,1.25,N_TC)
g_KL_TC_vals = np.random.uniform(0.55,0.6,N_TC)

TC.gg_h_TC = g_h_TC_vals * (0.02 * msiemens * cm ** -2)
TC.gg_KL_TC = g_KL_TC_vals * (0.02 * msiemens * cm ** -2)
#TC.gg_h_TC = array([0.005, 0.02]) * (msiemens * cm ** -2) 
#TC.gg_KL_TC = array([0.005, 0.004]) * usiemens / area_TC
TC.o2=0; TC.p1=0; TC.c1=1;
 
#TC.I=linked_var(input_TC,'B')
TC.Ca_i = 2.4e-4 * mMolar # 1 M = 1e-3 mole/cm3  
TC.ECa_TC = 120 * mV #+ (rand(len(TC)) * 10 - 10) * mV
TC.v = -68*mV # + (randn(len(TC)) * 5 - 5) * mV
"""
x0 = INaKTK_init(-68.); TC.m = x0[0]; TC.h = x0[1]; TC.n = x0[2];
x0 = IT_init(-68.,phi_m_TC,phi_h_TC,36,1e-4, 2); TC.mCa = x0[0]; TC.hCa = x0[1]; TC.ECa_TC = x0[2]*mV;
x0 = Ih_init(-68.,1e-4,hadj_TC); TC.p1 = x0[0]; TC.c1 = x0[1]; TC.o2 = x0[2];
x0 = IA_init(-68.); TC.am = x0[0]; TC.ah = x0[1]
"""

# RE nöron hazırlık
g_KL_RE_vals = np.random.uniform(0.5,1.5,N_RE)
RE.gg_KL_RE = g_KL_RE_vals * (0.002 * msiemens * cm ** -2) 

#RE.I=linked_var(input_RE,'B')
RE.Ca_i = 2.4e-4 * mMolar # 1 M = 1e-3 mole/cm3  
RE.ECa_RE = 120 * mV #+ (rand(len(RE)) * 10 - 10) * mV
RE.v = -61*mV # np.random.normal(-61,0.05,N_RE)*mV # + (randn(len(RE)) * 5 - 5) * mV
"""
x0 = INaKTK_init(-61.); RE.m = x0[0]; RE.h = x0[1]; RE.n = x0[2];
x0 = ITs_init(-61.,phi_m_RE,phi_h_RE,36,1e-4, 2); RE.mCa = x0[0]; RE.hCa = x0[1]; RE.ECa_RE = x0[2]*mV;
"""
#----------------------------------------------------------------
 
########################### Kayıtlar ###################################
# PY hücre Akım veya kinetik kayıtları
M_PY = SpikeMonitor(PY)
frate_PY = PopulationRateMonitor(PY)
MVd_PY = StateMonitor(PY, 'vd', record=True, clock=clock_rec)
MVs_PY = StateMonitor(PY, 'vs', record=True, clock=clock_rec)
LFP_PY = StateMonitor(PY, 'LFP', record=True, clock=clock_rec)
LFPx_PY = StateMonitor(PY, 'LFPx', record=True, clock=clock_rec)

M_PY_input = SpikeMonitor(PY_input)
 
MI_AMPA_PY_PY = StateMonitor(PY, 'I_AMPA_PY_PY', record=True, clock=clock_rec)
MI_NMDA_PY_PY = StateMonitor(PY, 'I_NMDA_PY_PY', record=True, clock=clock_rec)
MI_GABAA_IN_PY = StateMonitor(PY, 'I_GABAA_IN_PY', record=True, clock=clock_rec)
Mg_AMPA_PY_PY = StateMonitor(PY, 'g_AMPA_PY_PY', record=True, clock=clock_rec)
Mg_NMDA_PY_PY = StateMonitor(PY, 'g_NMDA_PY_PY', record=True, clock=clock_rec)
Mg_GABAA_IN_PY = StateMonitor(PY, 'g_GABAA_IN_PY', record=True, clock=clock_rec)

# IN hücre Akım veya kinetik kayıtları
M_IN = SpikeMonitor(IN)
frate_IN = PopulationRateMonitor(IN)
MVd_IN = StateMonitor(IN, 'vd', record=True, clock=clock_rec)
MVs_IN = StateMonitor(IN, 'vs', record=True, clock=clock_rec)

LFP_IN = StateMonitor(IN, 'LFP', record=True, clock=clock_rec)
LFPx_IN = StateMonitor(IN, 'LFPx', record=True, clock=clock_rec)

M_IN_input = SpikeMonitor(IN_input)

"""
MI_IN=StateMonitor(IN,'I',record=True)
MI_AMPA_PY_IN = StateMonitor(IN, 'I_AMPA_PY_IN', record=True, clock=clock_rec)
MI_NMDA_PY_IN = StateMonitor(IN, 'I_NMDA_PY_IN', record=True, clock=clock_rec)
Mg_AMPA_PY_IN = StateMonitor(IN, 'g_AMPA_PY_IN', record=True, clock=clock_rec)
Mg_NMDA_PY_IN = StateMonitor(IN, 'g_NMDA_PY_IN', record=True, clock=clock_rec)

MDCa_i_PY = StateMonitor(PY, 'DCa_i', record=True, clock=clock_rec)
MDCa_i_IN = StateMonitor(IN, 'DCa_i', record=True, clock=clock_rec)
 
MEDCa_PY = StateMonitor(PY, 'EDCa_PY', record=True, clock=clock_rec)
MEDCa_IN = StateMonitor(IN, 'EDCa_IN', record=True, clock=clock_rec)
"""
#----------------------------------------------------------------
 
# TC hücre Akım veya kinetik kayıtları 
#M_TC_rate = PopulationRateMonitor(TC)
M_TC = SpikeMonitor(TC)
frate_TC = PopulationRateMonitor(TC)
MV_TC = StateMonitor(TC, 'v', record=True, clock=clock_rec)

LFP_TC = StateMonitor(TC, 'LFP', record=True, clock=clock_rec)
LFPx_TC = StateMonitor(TC, 'LFPx', record=True, clock=clock_rec)

M_TC_input = SpikeMonitor(TC_input)

"""
MI_TC = StateMonitor(TC,'I', record=True, clock=clock_rec)
MIL_TC = StateMonitor(TC,'I_L', record=True, clock=clock_rec)
MIKL_TC = StateMonitor(TC,'I_KL', record=True, clock=clock_rec)
Mm_TC = StateMonitor(TC, 'm', record=True, clock=clock_rec)
Mh_TC = StateMonitor(TC, 'h', record=True, clock=clock_rec)
Mn_TC = StateMonitor(TC, 'n', record=True, clock=clock_rec)
"""
 
#MIT_TC = StateMonitor(TC,'I_T', record=True, clock=clock_rec)
"""
MmCa_TC = StateMonitor(TC,'mCa', record=True, clock=clock_rec)
MhCa_TC = StateMonitor(TC,'hCa', record=True, clock=clock_rec)
MmCa_inf_TC = StateMonitor(TC,'mCa_inf', record=True, clock=clock_rec)
MhCa_inf_TC = StateMonitor(TC,'hCa_inf', record=True, clock=clock_rec)
Mtau_hCa_TC = StateMonitor(TC,'tau_hCa', record=True, clock=clock_rec)
"""
 
#MIh_TC = StateMonitor(TC,'I_h', record=True, clock=clock_rec)
"""
Mmh_inf_TC = StateMonitor(TC,'mh_inf', record=True, clock=clock_rec)
Mtau_mh_TC = StateMonitor(TC,'tau_mh', record=True, clock=clock_rec)
Malphamh_TC = StateMonitor(TC,'alphamh', record=True, clock=clock_rec)
Mbetamh_TC = StateMonitor(TC,'betamh', record=True, clock=clock_rec)
Mk1ca_TC = StateMonitor(TC,'k1ca', record=True, clock=clock_rec)
Mk3p_TC = StateMonitor(TC,'k3p', record=True, clock=clock_rec)
Mc1_TC = StateMonitor(TC,'c1', record=True, clock=clock_rec)
Mp0_TC = StateMonitor(TC,'p0', record=True, clock=clock_rec)
Mo1_TC = StateMonitor(TC,'o1', record=True, clock=clock_rec)
Mo2_TC = StateMonitor(TC,'o2', record=True, clock=clock_rec)
"""
 
#MI_GABAA_RE_TC = StateMonitor(TC, 'I_GABAA_RE_TC', record=True, clock=clock_rec)
#MI_GABAB_RE_TC = StateMonitor(TC, 'I_GABAB_RE_TC', record=True, clock=clock_rec)
 
#----------------------------------------------------------------
 
# RE hücre Akım veya kinetik kayıtları
M_RE = SpikeMonitor(RE)
frate_RE = PopulationRateMonitor(RE)
MV_RE = StateMonitor(RE, 'v', record=True, clock=clock_rec)

LFP_RE = StateMonitor(RE, 'LFP', record=True, clock=clock_rec)
LFPx_RE = StateMonitor(RE, 'LFPx', record=True, clock=clock_rec)

M_RE_input = SpikeMonitor(RE_input)

#MI_RE=StateMonitor(RE,'I',record=True, clock=clock_rec)
 
#MIT_RE = StateMonitor(RE,'I_Ts', record=True, clock=clock_rec)
"""
MmCa_RE = StateMonitor(RE,'mCa', record=True, clock=clock_rec)
MhCa_RE = StateMonitor(RE,'hCa', record=True, clock=clock_rec)
MmCa_inf_RE = StateMonitor(RE,'mCa_inf', record=True, clock=clock_rec)
MhCa_inf_RE = StateMonitor(RE,'hCa_inf', record=True, clock=clock_rec)
Mtau_mCa_RE = StateMonitor(RE,'tau_mCa', record=True, clock=clock_rec)
Mtau_hCa_RE = StateMonitor(RE,'tau_hCa', record=True, clock=clock_rec)
"""
 
#MI_AMPA_TC_RE = StateMonitor(RE, 'I_AMPA_TC_RE', record=True, clock=clock_rec)
#MI_GABAA_RE_RE = StateMonitor(RE, 'I_GABAA_RE_RE', record=True, clock=clock_rec)
 
"""
MS_AMPA_TC_RE_x = StateMonitor(AMPA_TC_RE,'x_ampa',record=True, clock=clock_rec)
MS_GABAA_RE_RE_x = StateMonitor(GABAA_RE_RE,'x_gabaa',record=True, clock=clock_rec)
MS_GABAA_RE_TC_x = StateMonitor(GABAA_RE_TC,'x_gabaa',record=True, clock=clock_rec)
MS_GABAB_RE_TC_x = StateMonitor(GABAB_RE_TC,'x_gabab',record=True, clock=clock_rec)
 
MS_AMPA_TC_RE_g = StateMonitor(AMPA_TC_RE,'gampa',record=True, clock=clock_rec)
MS_GABAA_RE_RE_g = StateMonitor(GABAA_RE_RE,'ggabaa',record=True, clock=clock_rec)
MS_GABAA_RE_TC_g = StateMonitor(GABAA_RE_TC,'ggabaa',record=True, clock=clock_rec)
MS_GABAB_RE_TC_g = StateMonitor(GABAB_RE_TC,'ggabab',record=True, clock=clock_rec)
 
MS_AMPA_TC_RE_T = StateMonitor(AMPA_TC_RE,'T',record=True, clock=clock_rec)
MS_GABAA_RE_RE_T = StateMonitor(GABAA_RE_RE,'T',record=True, clock=clock_rec)
MS_GABAA_RE_TC_T = StateMonitor(GABAA_RE_TC,'T',record=True, clock=clock_rec)
MS_GABAB_RE_TC_T = StateMonitor(GABAB_RE_TC,'T',record=True, clock=clock_rec)
"""
"""
MCa_i_TC = StateMonitor(TC, 'Ca_i', record=True, clock=clock_rec)
MCa_i_RE = StateMonitor(RE, 'Ca_i', record=True, clock=clock_rec)
 
MECa_TC = StateMonitor(TC, 'ECa_TC', record=True, clock=clock_rec)
MECa_RE = StateMonitor(RE, 'ECa_RE', record=True, clock=clock_rec)
"""
####################### simülasyon ######################

print "Simulasyon başlatılıyor..."

run(duration,report='text')

######################### son ############################

print('spike sayilari: PY=%s / IN=%s / TC=%s / RE=%s' \
        %(str(M_PY.nspikes),str(M_IN.nspikes),str(M_TC.nspikes),str(M_RE.nspikes)))

####################### Analizler #######################
PY_spikes=M_PY.spiketimes
IN_spikes=M_IN.spiketimes
TC_spikes=M_TC.spiketimes
RE_spikes=M_RE.spiketimes
"""
M_PY_CVs, M_PY_FRs, M_PY_VSs, M_PY_ISIs = [], [], [], [];
for i in range(N_PY):
    #M_PY_CVs.append(statistics.CV(M_PY[i])) # CV_and_CC_calc ile aynı
    #M_PY_FRs.append(firing_rate(M_PY[i])) # network_mean_frequency ile aynı
    #M_PY_VSs.append(vector_strength(M_PY[i],1*second/f_ext)) # population_VS_Calc ile aynı
    M_PY_ISIs.extend(myISI(PY_spikes[i],1000))

M_IN_CVs, M_IN_FRs, M_IN_VSs, M_IN_ISIs = [], [], [], []; 
for i in range(N_IN):
    #M_IN_CVs.append(statistics.CV(M_IN[i])) # CV_and_CC_calc ile aynı
    #M_IN_FRs.append(firing_rate(M_IN[i])) # network_mean_frequency ile aynı
    #M_IN_VSs.append(vector_strength(M_IN[i],1*second/f_ext)) # population_VS_Calc ile aynı
    M_IN_ISIs.extend(myISI(IN_spikes[i],1000))

M_TC_CVs, M_TC_FRs, M_TC_VSs, M_TC_ISIs = [], [], [], [];
for i in range(N_TC):
    #M_TC_CVs.append(statistics.CV(M_TC[i])) # CV_and_CC_calc ile aynı
    #M_TC_FRs.append(firing_rate(M_TC[i])) # network_mean_frequency ile aynı
    #M_TC_VSs.append(vector_strength(M_TC[i],1*second/f_ext)) # population_VS_Calc ile aynı
    M_TC_ISIs.extend(myISI(TC_spikes[i],1000))

M_RE_CVs, M_RE_FRs, M_RE_VSs, M_RE_ISIs = [], [], [], []; 
for i in range(N_RE):
    #M_RE_CVs.append(statistics.CV(M_RE[i])) # CV_and_CC_calc ile aynı
    #M_RE_FRs.append(firing_rate(M_RE[i])) # network_mean_frequency ile aynı
    #M_RE_VSs.append(vector_strength(M_RE[i],1*second/f_ext)) # population_VS_Calc ile aynı
    M_RE_ISIs.extend(myISI(RE_spikes[i],1000))
"""

############################ LFP ############################

# toplam ve ort. PY voltaj ve LFP
total_LFP_PY = sum(LFP_PY.values, axis=0); total_LFPx_PY = sum(LFPx_PY.values, axis=0)
mean_LFP_PY = mean(LFP_PY.values, axis=0); mean_LFPx_PY = mean(LFPx_PY.values, axis=0)
total_membrane_PY = sum(MVs_PY.values, axis=0); mean_membrane_PY = mean(MVs_PY.values, axis=0)
# toplam ve ort. IN voltaj ve LFP
total_LFP_IN = sum(LFP_IN.values, axis=0); total_LFPx_IN = sum(LFPx_IN.values, axis=0)
mean_LFP_IN = mean(LFP_IN.values, axis=0); mean_LFPx_IN = mean(LFPx_IN.values, axis=0)
total_membrane_IN = sum(MVs_IN.values, axis=0); mean_membrane_IN = mean(MVs_IN.values, axis=0)
# toplam ve ort. TC voltaj ve LFP
total_LFP_TC = sum(LFP_TC.values, axis=0); total_LFPx_TC = sum(LFPx_TC.values, axis=0)
mean_LFP_TC = mean(LFP_TC.values, axis=0); mean_LFPx_TC = mean(LFPx_TC.values, axis=0)
total_membrane_TC = sum(MV_TC.values, axis=0); mean_membrane_TC = mean(MV_TC.values, axis=0)
# toplam ve ort. RE voltaj ve LFP
total_LFP_RE = sum(LFP_RE.values, axis=0); total_LFPx_RE = sum(LFPx_RE.values, axis=0)
mean_LFP_RE = mean(LFP_RE.values, axis=0); mean_LFPx_RE = mean(LFPx_RE.values, axis=0)
total_membrane_RE = sum(MV_RE.values, axis=0); mean_membrane_RE = mean(MV_RE.values, axis=0)

"""
#----------------------------------------------------------------------------

# empedans hesaplanıyor (Ohm)

R = 105.0*um;                    # kaynak yarıçapı (um)
rext = 500.0*um;                 # kayıt elektrodunun uzaklığı (um)
rmax = 200*R;               # max integrasyon mesafesi (um)
dr = rmax/1000;             # integrasyon adımı (um)

#sigmaR = 1.56*siemens/meter;  			# kaynak yakınında mutlak ietkenlik değeri (S/m)
sigmaR = 0.3*siemens/meter;
#sigmaR = sigmaR * 1e-6;                # um birim cinsinden iletkenlik

Lambda = 5*500*um;       		   # iletkenlik zayıflama uzay sabiti (um)
sigma1 = sigmaR/sigmaR;             	# normalize edilmiş max iletkenlik değeri
sigma2 = 0.2;       		          # normalize edilmiş min iletkenlik değeri
                                        # (makroskopik asimptotik değer)
sigma2 = (3.5e-9*siemens/meter)/sigmaR;  # normalize edilmiş min iletkenlik değeri
                                        # (membranda min değer)
epsilon1 = (1e-10*farad/meter)/sigmaR;  # sabit normalize edilmiş dielektrik sabite
#epsilon1 = 0.0001*second;                # heuristik

# herbir PY nöronun sinaptik akımların LFP etkisi ayrı ayrı değerlendiriliyor
dt = defaultclock.dt/second # 0.1ms = 0.0001 sec
Vexts_PY = [];
for i in range(N_PY):
    rexts = (100+rand()*1000)*um;
    x=LFP_PY.values[i]*area_DEND_PY;
    x, Zr, Zi, f, Y, Ir, Ii, Vext = calc_lfp(x,dt,rexts,rmax,dr, \
        R,sigma1,sigma2,Lambda,epsilon1,sigmaR);
    Vexts_PY.append(Vext)
    print (i)

Vexts_PY = squeeze(array(Vexts_PY)[:,:len(LFP_PY.times):])

figure();
subplot(211); plot(LFP_PY.times/ms, LFP_PY.values[0]/(nA/area_DEND_PY))
ylabel('PY LFP (nA)')
subplot(212); plot(LFP_PY.times/ms,Vexts_PY[0]/mV)
ylabel('PY Vexts (mV)')
show();

# herbir IN nöronun sinaptik akımların LFP etkisi ayrı ayrı değerlendiriliyor
dt = defaultclock.dt/second # 0.1ms = 0.0001 sec
Vexts_IN = [];
for i in range(N_IN):
    rexts = (100+rand()*1000)*um;
    x=LFP_IN.values[i]*area_DEND_IN;
    x, Zr, Zi, f, Y, Ir, Ii, Vext = calc_lfp(x,dt,rexts,rmax,dr, \
        R,sigma1,sigma2,Lambda,epsilon1,sigmaR);
    Vexts_IN.append(Vext)
    print (i)

Vexts_IN = squeeze(array(Vexts_IN)[:,:len(LFP_IN.times):])

# toplam ve ort. LFP'ler
total_Vexts_PY = sum(Vexts_PY, axis=0); mean_Vexts_PY = mean(Vexts_PY, axis=0);
total_Vexts_IN = sum(Vexts_IN, axis=0); mean_Vexts_IN = mean(Vexts_IN, axis=0);

#################### Spektral Güç Dağılımı (PSD) ############################
npow1=512 
npow2=int(1*second/dt)

figure()
subplot(331); plot(LFP_PY.times/ms,total_Vexts_PY/mV); ylabel('Voltaj (mV)'); title('PY total Vexts (mV)')
subplot(334); Pxx_PY, freqs_PY, bins, im = specgram(detrend_linear(total_Vexts_PY/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); ylabel('frequency (Hz)')
subplot(337); PxxTV_PY, freqsTV_PY = psd(detrend_linear(total_Vexts_PY/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(332); plot(LFP_PY.times/ms,mean_Vexts_PY/mV); title('PY mean Vexts (mV)')
subplot(335); Pxx_PY, freqs_PY, bins, im = specgram(detrend_linear(mean_Vexts_PY/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); 
subplot(338); PxxMV_PY, freqsMV_PY = psd(detrend_linear(mean_Vexts_PY/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt,label=''); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(333); plot(MVs_PY.times/ms,total_membrane_PY/mV); title('total membran (mV)')
subplot(336); Pxx_PY, freqs_PY, bins, im = specgram(detrend_linear(total_membrane_PY/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)');
subplot(339); PxxTM_PY, freqsTM_PY = psd(detrend_linear(total_membrane_PY/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') #ylabel(r'$PSD (mV^2/Hz)$')

figure()
subplot(331); plot(LFP_IN.times/ms,total_Vexts_IN/mV); ylabel('Voltaj (mV)'); title('IN total Vexts (mV)')
subplot(334); Pxx_IN, freqs_IN, bins, im = specgram(detrend_linear(total_Vexts_IN/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); ylabel('frequency (Hz)')
subplot(337); PxxTV_IN, freqsTV_IN = psd(detrend_linear(total_Vexts_IN/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(332); plot(LFP_IN.times/ms,mean_Vexts_IN/mV); title('IN mean Vexts (mV)')
subplot(335); Pxx_IN, freqs_IN, bins, im = specgram(detrend_linear(mean_Vexts_IN/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); 
subplot(338); PxxMV_IN, freqsMV_IN = psd(detrend_linear(mean_Vexts_IN/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt,label=''); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(333); plot(MVs_IN.times/ms,total_membrane_IN/mV); title('total membran (mV)')
subplot(336); Pxx_IN, freqs_IN, bins, im = specgram(detrend_linear(total_membrane_IN/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)');
subplot(339); PxxTM_IN, freqsTM_IN = psd(detrend_linear(total_membrane_IN/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') #ylabel(r'$PSD (mV^2/Hz)$')

show();

##################### PY için stokastik rezonasn (SR) #######################
# toplam memb. voltaja göre
maxindTM_PY, maxvalTM_PY = max(enumerate(PxxTM_PY), key=operator.itemgetter(1))
base1=int(maxindTM_PY*bases[0]); base2=int(maxindTM_PY*bases[1]);
if base1==base2:base1=maxindTM_PY
snrTM_PY1,bindxleft,bindexright = stochastic_resonance(PxxTM_PY,base1,base2) 
snrTM_PY2,bindxleft,bindexright = stochastic_resonance(PxxTM_PY)

# toplam Vexts voltaja göre
maxindTV_PY, maxvalTV_PY = max(enumerate(PxxTV_PY), key=operator.itemgetter(1))
base1=int(maxindTV_PY*bases[0]); base2=int(maxindTV_PY*bases[1]);
if base1==base2:base1=maxindTV_PY
snrTV_PY1,bindxleft,bindexright = stochastic_resonance(PxxTV_PY,base1,base2) 
snrTV_PY2,bindxleft,bindexright = stochastic_resonance(PxxTV_PY)

# ort. Vexts voltaja göre
maxindMV_PY, maxvalMV_PY = max(enumerate(PxxMV_PY), key=operator.itemgetter(1))
base1=int(maxindMV_PY*bases[0]); base2=int(maxindMV_PY*bases[1]);
if base1==base2:base1=maxindMV_PY
snrMV_PY1,bindxleft,bindexright = stochastic_resonance(PxxMV_PY,base1,base2) 
snrMV_PY2,bindxleft,bindexright = stochastic_resonance(PxxMV_PY)

##################### IN için stokastik rezonasn (SR) #######################
# toplam memb. voltaja göre
maxindTM_IN, maxvalTM_IN = max(enumerate(PxxTM_IN), key=operator.itemgetter(1))
base1=int(maxindTM_IN*bases[0]); base2=int(maxindTM_IN*bases[1]);
if base1==base2:base1=maxindTM_IN
snrTM_IN1,bindxleft,bindexright = stochastic_resonance(PxxTM_IN,base1,base2) 
snrTM_IN2,bindxleft,bindexright = stochastic_resonance(PxxTM_IN)

# toplam Vexts voltaja göre
maxindTV_IN, maxvalTV_IN = max(enumerate(PxxTV_IN), key=operator.itemgetter(1))
base1=int(maxindTV_IN*bases[0]); base2=int(maxindTV_IN*bases[1]);
if base1==base2:base1=maxindTV_IN
snrTV_IN1,bindxleft,bindexright = stochastic_resonance(PxxTV_IN,base1,base2) 
snrTV_IN2,bindxleft,bindexright = stochastic_resonance(PxxTV_IN)

# ort. Vexts voltaja göre
maxindMV_IN, maxvalMV_IN = max(enumerate(PxxMV_IN), key=operator.itemgetter(1))
base1=int(maxindMV_IN*bases[0]); base2=int(maxindMV_IN*bases[1]);
if base1==base2:base1=maxindMV_IN
snrMV_IN1,bindxleft,bindexright = stochastic_resonance(PxxMV_IN,base1,base2) 
snrMV_IN2,bindxleft,bindexright = stochastic_resonance(PxxMV_IN)

####################### PY için koherens rezonans (CR) ####################
# toplam memb. voltaja göre
maxindTM_PY, maxvalTM_PY = max(enumerate(PxxTM_PY), key=operator.itemgetter(1))
base1=int(maxindTM_PY*bases[0]); base2=int(maxindTM_PY*bases[1]);
if base1==base2:base1=maxindTM_PY
CFactor_betaTM_PY1 = coherence_resonance(PxxTM_PY,freqsTM_PY,base1,base2); # koherens1
CFactor_betaTM_PY2 = coherence_resonance(PxxTM_PY,freqsTM_PY); # koherens2

# toplam Vexts voltaja göre
maxindTV_PY, maxvalTV_PY = max(enumerate(PxxTV_PY), key=operator.itemgetter(1))
base1=int(maxindTV_PY*bases[0]); base2=int(maxindTV_PY*bases[1]);
if base1==base2:base1=maxindTV_PY
CFactor_betaTV_PY1 = coherence_resonance(PxxTV_PY,freqsTV_PY,base1,base2); # koherens1
CFactor_betaTV_PY2 = coherence_resonance(PxxTV_PY,freqsTV_PY); # koherens2

# ort. Vexts voltaja göre
maxindMV_PY, maxvalMV_PY = max(enumerate(PxxMV_PY), key=operator.itemgetter(1))
base1=int(maxindMV_PY*bases[0]); base2=int(maxindMV_PY*bases[1]);
if base1==base2:base1=maxindMV_PY
CFactor_betaMV_PY1 = coherence_resonance(PxxMV_PY,freqsMV_PY,base1,base2); # koherens1
CFactor_betaMV_PY2 = coherence_resonance(PxxMV_PY,freqsMV_PY); # koherens2

####################### IN için koherens rezonans (CR) ####################
# toplam memb. voltaja göre
maxindTM_IN, maxvalTM_IN = max(enumerate(PxxTM_IN), key=operator.itemgetter(1))
base1=int(maxindTM_IN*bases[0]); base2=int(maxindTM_IN*bases[1]);
if base1==base2:base1=maxindTM_IN
CFactor_betaTM_IN1 = coherence_resonance(PxxTM_IN,freqsTM_IN,base1,base2); # koherens1
CFactor_betaTM_IN2 = coherence_resonance(PxxTM_IN,freqsTM_IN); # koherens2

# toplam Vexts voltaja göre
maxindTV_IN, maxvalTV_IN = max(enumerate(PxxTV_IN), key=operator.itemgetter(1))
base1=int(maxindTV_IN*bases[0]); base2=int(maxindTV_IN*bases[1])
if base1==base2:base1=maxindTV_IN
CFactor_betaTV_IN1 = coherence_resonance(PxxTV_IN,freqsTV_IN,base1,base2); # koherens1
CFactor_betaTV_IN2 = coherence_resonance(PxxTV_IN,freqsTV_IN); # koherens2

# ort. Vexts voltaja göre
maxindMV_IN, maxvalMV_IN = max(enumerate(PxxMV_IN), key=operator.itemgetter(1))
base1=int(maxindMV_IN*bases[0]); base2=int(maxindMV_IN*bases[1]);
if base1==base2:base1=maxindMV_IN
CFactor_betaMV_IN1 = coherence_resonance(PxxMV_IN,freqsMV_IN,base1,base2); # koherens1
CFactor_betaMV_IN2 = coherence_resonance(PxxMV_IN,freqsMV_IN); # koherens2

"""

#################### Spektral Güç Dağılımı (PSD) ############################
#dt=0.0001 # 0.1ms = 0.0001 sec
#nextpow2=32768 # here 2000ms of length
npow1=512 
npow2=int(duration/dt)

# psd(...) komutu ile Welch ortlama periodogram yöntemiyle power spectral density hesaplanıyor
############################ PY PSD ####################################
figure()
subplot(331); plot(LFP_PY.times/ms,total_LFP_PY/(nA/area_DEND_PY)); ylabel(r'$Amplitude (nA)$'); title('PY total LFP (nA)')
subplot(334); PxxTP_PY, freqsTP_PY, bins, im = specgram(detrend_linear(total_LFP_PY/(nA/area_DEND_PY)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); ylabel('frequency (Hz)')
subplot(337); PxxTP_PY, freqsTP_PY = psd(detrend_linear(total_LFP_PY/(nA/area_DEND_PY)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(332); plot(LFP_PY.times/ms,mean_LFP_PY/(nA/area_DEND_PY)); title('PY mean LFP (nA)')
subplot(335); PxxMP_PY, freqsMP_PY, bins, im = specgram(detrend_linear(mean_LFP_PY/(nA/area_DEND_PY)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); 
subplot(338); PxxMP_PY, freqsMP_PY = psd(detrend_linear(mean_LFP_PY/(nA/area_DEND_PY)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt,label=''); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(333); plot(MVs_PY.times/ms,total_membrane_PY/mV); title('PY total membran (mV)')
subplot(336); PxxTM_PY, freqsTM_PY, bins, im = specgram(detrend_linear(total_membrane_PY/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)');
subplot(339); PxxTM_PY, freqsTM_PY = psd(detrend_linear(total_membrane_PY/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') #ylabel(r'$PSD (mV^2/Hz)$')

figure()
subplot(331); plot(LFPx_PY.times/ms,total_LFPx_PY/(nA/area_DEND_PY)); ylabel(r'$Amplitude (nA)$'); title('PY total LFPx (nA)')
subplot(334); PxxTPx_PY, freqsTPx_PY, bins, im = specgram(detrend_linear(total_LFPx_PY/(nA/area_DEND_PY)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); ylabel('frequency (Hz)')
subplot(337); PxxTPx_PY, freqsTPx_PY = psd(detrend_linear(total_LFPx_PY/(nA/area_DEND_PY)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(332); plot(LFPx_PY.times/ms,mean_LFPx_PY/(nA/area_DEND_PY)); title('PY mean LFPx (nA)')
subplot(335); PxxMPx_PY, freqsMPx_PY, bins, im = specgram(detrend_linear(mean_LFPx_PY/(nA/area_DEND_PY)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); 
subplot(338); PxxMPx_PY, freqsMPx_PY = psd(detrend_linear(mean_LFPx_PY/(nA/area_DEND_PY)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt,label=''); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(333); plot(MVs_PY.times/ms,mean_membrane_PY/mV); title('PY mean membran (mV)')
subplot(336); PxxTM_PY, freqsTM_PY, bins, im = specgram(detrend_linear(mean_membrane_PY/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)');
subplot(339); PxxTM_PY, freqsTM_PY = psd(detrend_linear(mean_membrane_PY/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (mV^2/Hz)$')

############################ IN PSD ####################################
figure()
subplot(331); plot(LFP_IN.times/ms,total_LFP_IN/(nA/area_DEND_IN)); ylabel(r'$Amplitude (nA)$'); title('IN total LFP (nA)')
subplot(334); PxxTP_IN, freqsTP_IN, bins, im = specgram(detrend_linear(total_LFP_IN/(nA/area_DEND_IN)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); ylabel('frequency (Hz)')
subplot(337); PxxTP_IN, freqsTP_IN = psd(detrend_linear(total_LFP_IN/(nA/area_DEND_IN)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(332); plot(LFP_IN.times/ms,mean_LFP_IN/(nA/area_DEND_IN)); title('IN mean LFP (nA)')
subplot(335); PxxMP_IN, freqsMP_IN, bins, im = specgram(detrend_linear(mean_LFP_IN/(nA/area_DEND_IN)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); 
subplot(338); PxxMP_IN, freqsMP_IN = psd(detrend_linear(mean_LFP_IN/(nA/area_DEND_IN)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt,label=''); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(333); plot(MVs_IN.times/ms,total_membrane_IN/mV); title('IN total membran (mV)')
subplot(336); PxxTM_IN, freqsTM_IN, bins, im = specgram(detrend_linear(total_membrane_IN/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)');
subplot(339); PxxTM_IN, freqsTM_IN = psd(detrend_linear(total_membrane_IN/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') #ylabel(r'$PSD (mV^2/Hz)$')

############################ TC PSD ####################################
figure()
subplot(331); plot(LFP_TC.times/ms,total_LFP_TC/(nA/area_TC)); ylabel(r'$Amplitude (nA)$'); title('TC total LFP (nA)')
subplot(334); PxxTP_TC, freqsTP_TC, bins, im = specgram(detrend_linear(total_LFP_TC/(nA/area_TC)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); ylabel('frequency (Hz)')
subplot(337); PxxTP_TC, freqsTP_TC = psd(detrend_linear(total_LFP_TC/(nA/area_TC)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(332); plot(LFP_TC.times/ms,mean_LFP_TC/(nA/area_TC)); title('TC mean LFP (nA)')
subplot(335); PxxMP_TC, freqsMP_TC, bins, im = specgram(detrend_linear(mean_LFP_TC/(nA/area_TC)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); 
subplot(338); PxxMP_TC, freqsMP_TC = psd(detrend_linear(mean_LFP_TC/(nA/area_TC)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt,label=''); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(333); plot(MV_TC.times/ms,total_membrane_TC/mV); title('TC total membran (mV)')
subplot(336); PxxTM_TC, freqsTM_TC, bins, im = specgram(detrend_linear(total_membrane_TC/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)');
subplot(339); PxxTM_TC, freqsTM_TC = psd(detrend_linear(total_membrane_TC/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') #ylabel(r'$PSD (mV^2/Hz)$')

figure()
subplot(331); plot(LFPx_TC.times/ms,total_LFPx_TC/(nA/area_TC)); ylabel(r'$Amplitude (nA)$'); title('TC total LFPx (nA)')
subplot(334); PxxTPx_TC, freqsTPx_TC, bins, im = specgram(detrend_linear(total_LFPx_TC/(nA/area_TC)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); ylabel('frequency (Hz)')
subplot(337); PxxTPx_TC, freqsTPx_TC = psd(detrend_linear(total_LFPx_TC/(nA/area_TC)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(332); plot(LFPx_TC.times/ms,mean_LFPx_TC/(nA/area_TC)); title('TC mean LFPx (nA)')
subplot(335); PxxMPx_TC, freqsMPx_TC, bins, im = specgram(detrend_linear(mean_LFPx_TC/(nA/area_TC)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); 
subplot(338); PxxMPx_TC, freqsMPx_TC = psd(detrend_linear(mean_LFPx_TC/(nA/area_TC)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt,label=''); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(333); plot(MV_TC.times/ms,mean_membrane_TC/mV); title('TC mean membran (mV)')
subplot(336); PxxTM_TC, freqsTM_TC, bins, im = specgram(detrend_linear(mean_membrane_TC/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)');
subplot(339); PxxTM_TC, freqsTM_TC = psd(detrend_linear(mean_membrane_TC/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (mV^2/Hz)$')

############################ RE PSD ####################################
figure()
subplot(331); plot(LFP_RE.times/ms,total_LFP_RE/(nA/area_RE)); ylabel(r'$Amplitude (nA)$'); title('RE total LFP (nA)')
subplot(334); PxxTP_RE, freqsTP_RE, bins, im = specgram(detrend_linear(total_LFP_RE/(nA/area_RE)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); ylabel('frequency (Hz)')
subplot(337); PxxTP_RE, freqsTP_RE = psd(detrend_linear(total_LFP_RE/(nA/area_RE)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(332); plot(LFP_RE.times/ms,mean_LFP_RE/(nA/area_RE)); title('RE mean LFP (nA)')
subplot(335); PxxMP_RE, freqsMP_RE, bins, im = specgram(detrend_linear(mean_LFP_RE/(nA/area_RE)), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)'); 
subplot(338); PxxMP_RE, freqsMP_RE = psd(detrend_linear(mean_LFP_RE/(nA/area_RE)), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt,label=''); ylabel(r'$PSD (dB/Hz)$') # ylabel(r'$PSD (\mu A^2/Hz)$')

subplot(333); plot(MV_RE.times/ms,total_membrane_RE/mV); title('RE total membran (mV)')
subplot(336); PxxTM_RE, freqsTM_RE, bins, im = specgram(detrend_linear(total_membrane_RE/mV), window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('time (ms)');
subplot(339); PxxTM_RE, freqsTM_RE = psd(detrend_linear(total_membrane_RE/mV), window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); ylabel(r'$PSD (dB/Hz)$') #ylabel(r'$PSD (mV^2/Hz)$')

# gürültü ve sistemin toplam LFP değeri arasında cross-covariance
# PY için
total_LFP_PY1=detrend_linear(total_LFP_PY/nA);total_LFP_PY2=detrend_linear((total_LFP_PY-total_LFPx_PY)/nA);
covTP_PY, lagTP_PY = xcov(total_LFP_PY1,total_LFP_PY2, maxlags=1000, norm='coeff')
mCovind_total_LFP_PY, mCovval_total_LFP_PY = max(enumerate(covTP_PY), key=operator.itemgetter(1))
# IN için
total_LFP_IN1=detrend_linear(total_LFP_IN/nA);total_LFP_IN2=detrend_linear((total_LFP_IN-total_LFPx_IN)/nA);
covTP_IN, lagTP_IN = xcov(total_LFP_IN1,total_LFP_IN2, maxlags=1000, norm='coeff')
mCovind_total_LFP_IN, mCovval_total_LFP_IN = max(enumerate(covTP_IN), key=operator.itemgetter(1))
# TC için
total_LFP_TC1=detrend_linear(total_LFP_TC/nA);total_LFP_TC2=detrend_linear((total_LFP_TC-total_LFPx_TC)/nA);
covTP_TC, lagTP_TC = xcov(total_LFP_TC1,total_LFP_TC2, maxlags=1000, norm='coeff')
mCovind_total_LFP_TC, mCovval_total_LFP_TC = max(enumerate(covTP_TC), key=operator.itemgetter(1))
# RE için
total_LFP_RE1=detrend_linear(total_LFP_RE/nA);total_LFP_RE2=detrend_linear((total_LFP_RE-total_LFPx_RE)/nA);
covTP_RE, lagTP_RE = xcov(total_LFP_RE1,total_LFP_RE2, maxlags=1000, norm='coeff')
mCovind_total_LFP_RE, mCovval_total_LFP_RE = max(enumerate(covTP_RE), key=operator.itemgetter(1))

figure(); 
subplot(221); plot(lagTP_PY,covTP_PY);title('PY Capraz-Kovaryans')
subplot(222); plot(lagTP_IN,covTP_IN);title('IN Capraz-Kovaryans')
subplot(223); plot(lagTP_TC,covTP_TC);title('TC Capraz-Kovaryans')
subplot(224); plot(lagTP_RE,covTP_RE);title('RE Capraz-Kovaryans')
show()

# SNR'i; SNR ∝ (ϵ∆U/D)exp(-∆U/D) denklemine göre çözmek için Stacey ve ark. 2009
# D = mean(mean_LFPx_PY/nA)**2 # Noise_intensity 
##################### PY için stokastik rezonans (SR) #######################
# toplam memb. voltaja göre
maxindTM_PY, maxvalTM_PY = max(enumerate(PxxTM_PY), key=operator.itemgetter(1))
base1=int(maxindTM_PY*bases[0]); base2=int(maxindTM_PY*bases[1]);
if base1==base2:base1=maxindTM_PY
snrTM_PY1,bindxleft,bindexright = stochastic_resonance(PxxTM_PY,base1,base2) 
snrTM_PY2,bindxleft,bindexright = stochastic_resonance(PxxTM_PY)
# toplam LFP akıma göre
maxindTP_PY, maxvalTP_PY = max(enumerate(PxxTP_PY), key=operator.itemgetter(1))
base1=int(maxindTP_PY*bases[0]); base2=int(maxindTP_PY*bases[1]);
if base1==base2:base1=maxindTP_PY
snrTP_PY1,bindxleft,bindexright = stochastic_resonance(PxxTP_PY,base1,base2) 
snrTP_PY2,bindxleft,bindexright = stochastic_resonance(PxxTP_PY)
# ort. LFP akıma göre
maxindMP_PY, maxvalMP_PY = max(enumerate(PxxMP_PY), key=operator.itemgetter(1))
base1=int(maxindMP_PY*bases[0]); base2=int(maxindMP_PY*bases[1]);
if base1==base2:base1=maxindMP_PY
snrMP_PY1,bindxleft,bindexright = stochastic_resonance(PxxMP_PY,base1,base2) 
snrMP_PY2,bindxleft,bindexright = stochastic_resonance(PxxMP_PY)

##################### IN için stokastik rezonans (SR) #######################
# toplam memb. voltaja göre
maxindTM_IN, maxvalTM_IN = max(enumerate(PxxTM_IN), key=operator.itemgetter(1))
base1=int(maxindTM_IN*bases[0]); base2=int(maxindTM_IN*bases[1]);
if base1==base2:base1=maxindTM_IN
snrTM_IN1,bindxleft,bindexright = stochastic_resonance(PxxTM_IN,base1,base2) 
snrTM_IN2,bindxleft,bindexright = stochastic_resonance(PxxTM_IN)
# toplam LFP akıma göre
maxindTP_IN, maxvalTP_IN = max(enumerate(PxxTP_IN), key=operator.itemgetter(1))
base1=int(maxindTP_IN*bases[0]); base2=int(maxindTP_IN*bases[1]);
if base1==base2:base1=maxindTP_IN
snrTP_IN1,bindxleft,bindexright = stochastic_resonance(PxxTP_IN,base1,base2) 
snrTP_IN2,bindxleft,bindexright = stochastic_resonance(PxxTP_IN)
# ort. LFP akıma göre
maxindMP_IN, maxvalMP_IN = max(enumerate(PxxMP_IN), key=operator.itemgetter(1))
base1=int(maxindMP_IN*bases[0]); base2=int(maxindMP_IN*bases[1]);
if base1==base2:base1=maxindMP_IN
snrMP_IN1,bindxleft,bindexright = stochastic_resonance(PxxMP_IN,base1,base2) 
snrMP_IN2,bindxleft,bindexright = stochastic_resonance(PxxMP_IN)

##################### TC için stokastik rezonasn (SR) #######################
# toplam memb. voltaja göre
maxindTM_TC, maxvalTM_TC = max(enumerate(PxxTM_TC), key=operator.itemgetter(1))
base1=int(maxindTM_TC*bases[0]); base2=int(maxindTM_TC*bases[1]);
if base1==base2:base1=maxindTM_TC
snrTM_TC1,bindxleft,bindexright = stochastic_resonance(PxxTM_TC,base1,base2) 
snrTM_TC2,bindxleft,bindexright = stochastic_resonance(PxxTM_TC)
# toplam LFP akıma göre
maxindTP_TC, maxvalTP_TC = max(enumerate(PxxTP_TC), key=operator.itemgetter(1))
base1=int(maxindTP_TC*bases[0]); base2=int(maxindTP_TC*bases[1]);
if base1==base2:base1=maxindTP_TC
snrTP_TC1,bindxleft,bindexright = stochastic_resonance(PxxTP_TC,base1,base2) 
snrTP_TC2,bindxleft,bindexright = stochastic_resonance(PxxTP_TC)
# ort. LFP akıma göre
maxindMP_TC, maxvalMP_TC = max(enumerate(PxxMP_TC), key=operator.itemgetter(1))
base1=int(maxindMP_TC*bases[0]); base2=int(maxindMP_TC*bases[1]);
if base1==base2:base1=maxindMP_TC
snrMP_TC1,bindxleft,bindexright = stochastic_resonance(PxxMP_TC,base1,base2) 
snrMP_TC2,bindxleft,bindexright = stochastic_resonance(PxxMP_TC)

##################### RE için stokastik rezonasn (SR) #######################
# toplam memb. voltaja göre
maxindTM_RE, maxvalTM_RE = max(enumerate(PxxTM_RE), key=operator.itemgetter(1))
base1=int(maxindTM_RE*bases[0]); base2=int(maxindTM_RE*bases[1]);
if base1==base2:base1=maxindTM_RE
snrTM_RE1,bindxleft,bindexright = stochastic_resonance(PxxTM_RE,base1,base2) 
snrTM_RE2,bindxleft,bindexright = stochastic_resonance(PxxTM_RE)
# toplam LFP akıma göre
maxindTP_RE, maxvalTP_RE = max(enumerate(PxxTP_RE), key=operator.itemgetter(1))
base1=int(maxindTP_RE*bases[0]); base2=int(maxindTP_RE*bases[1]);
if base1==base2:base1=maxindTP_RE
snrTP_RE1,bindxleft,bindexright = stochastic_resonance(PxxTP_RE,base1,base2) 
snrTP_RE2,bindxleft,bindexright = stochastic_resonance(PxxTP_RE)
# ort. LFP akıma göre
maxindMP_RE, maxvalMP_RE = max(enumerate(PxxMP_RE), key=operator.itemgetter(1))
base1=int(maxindMP_RE*bases[0]); base2=int(maxindMP_RE*bases[1]);
if base1==base2:base1=maxindMP_RE
snrMP_RE1,bindxleft,bindexright = stochastic_resonance(PxxMP_RE,base1,base2) 
snrMP_RE2,bindxleft,bindexright = stochastic_resonance(PxxMP_RE)

####################### PY için koherens rezonans (CR) ####################
# toplam memb. voltaja göre
maxindTM_PY, maxvalTM_PY = max(enumerate(PxxTM_PY), key=operator.itemgetter(1))
base1=int(maxindTM_PY*bases[0]); base2=int(maxindTM_PY*bases[1]);
if base1==base2:base1=maxindTM_PY
CFactor_betaTM_PY1 = coherence_resonance(PxxTM_PY,freqsTM_PY,base1,base2); # koherens1
CFactor_betaTM_PY2 = coherence_resonance(PxxTM_PY,freqsTM_PY); # koherens2
# toplam LFP akıma göre
maxindTP_PY, maxvalTP_PY = max(enumerate(PxxTP_PY), key=operator.itemgetter(1))
base1=int(maxindTP_PY*bases[0]); base2=int(maxindTP_PY*bases[1]);
if base1==base2:base1=maxindTP_PY
CFactor_betaTP_PY1 = coherence_resonance(PxxTP_PY,freqsTP_PY,base1,base2); # koherens1
CFactor_betaTP_PY2 = coherence_resonance(PxxTP_PY,freqsTP_PY); # koherens2
# ort. LFP akıma göre
maxindMP_PY, maxvalMP_PY = max(enumerate(PxxMP_PY), key=operator.itemgetter(1))
base1=int(maxindMP_PY*bases[0]); base2=int(maxindMP_PY*bases[1]);
if base1==base2:base1=maxindMP_PY
CFactor_betaMP_PY1 = coherence_resonance(PxxMP_PY,freqsMP_PY,base1,base2); # koherens1
CFactor_betaMP_PY2 = coherence_resonance(PxxMP_PY,freqsMP_PY); # koherens2

####################### IN için koherens rezonans (CR) ####################
# toplam memb. voltaja göre
maxindTM_IN, maxvalTM_IN = max(enumerate(PxxTM_IN), key=operator.itemgetter(1))
base1=int(maxindTM_IN*bases[0]); base2=int(maxindTM_IN*bases[1]);
if base1==base2:base1=maxindTM_IN
CFactor_betaTM_IN1 = coherence_resonance(PxxTM_IN,freqsTM_IN,base1,base2); # koherens1
CFactor_betaTM_IN2 = coherence_resonance(PxxTM_IN,freqsTM_IN); # koherens2
# toplam LFP akıma göre
maxindTP_IN, maxvalTP_IN = max(enumerate(PxxTP_IN), key=operator.itemgetter(1))
base1=int(maxindTP_IN*bases[0]); base2=int(maxindTP_IN*bases[1]);
if base1==base2:base1=maxindTP_IN
CFactor_betaTP_IN1 = coherence_resonance(PxxTP_IN,freqsTP_IN,base1,base2); # koherens1
CFactor_betaTP_IN2 = coherence_resonance(PxxTP_IN,freqsTP_IN); # koherens2
# ort. LFP akıma göre
maxindMP_IN, maxvalMP_IN = max(enumerate(PxxMP_IN), key=operator.itemgetter(1))
base1=int(maxindMP_IN*bases[0]); base2=int(maxindMP_IN*bases[1]);
if base1==base2:base1=maxindMP_IN
CFactor_betaMP_IN1 = coherence_resonance(PxxMP_IN,freqsMP_IN,base1,base2); # koherens1
CFactor_betaMP_IN2 = coherence_resonance(PxxMP_IN,freqsMP_IN); # koherens2

####################### TC için koherens rezonans (CR) ####################
# toplam memb. voltaja göre
maxindTM_TC, maxvalTM_TC = max(enumerate(PxxTM_TC), key=operator.itemgetter(1))
base1=int(maxindTM_TC*bases[0]); base2=int(maxindTM_TC*bases[1]);
if base1==base2:base1=maxindTM_TC
CFactor_betaTM_TC1 = coherence_resonance(PxxTM_TC,freqsTM_TC,base1,base2); # koherens1
CFactor_betaTM_TC2 = coherence_resonance(PxxTM_TC,freqsTM_TC); # koherens2
# toplam LFP akıma göre
maxindTP_TC, maxvalTP_TC = max(enumerate(PxxTP_TC), key=operator.itemgetter(1))
base1=int(maxindTP_TC*bases[0]); base2=int(maxindTP_TC*bases[1]);
if base1==base2:base1=maxindTP_TC
CFactor_betaTP_TC1 = coherence_resonance(PxxTP_TC,freqsTP_TC,base1,base2); # koherens1
CFactor_betaTP_TC2 = coherence_resonance(PxxTP_TC,freqsTP_TC); # koherens2
# ort. LFP akıma göre
maxindMP_TC, maxvalMP_TC = max(enumerate(PxxMP_TC), key=operator.itemgetter(1))
base1=int(maxindMP_TC*bases[0]); base2=int(maxindMP_TC*bases[1]);
if base1==base2:base1=maxindMP_TC
CFactor_betaMP_TC1 = coherence_resonance(PxxMP_TC,freqsMP_TC,base1,base2); # koherens1
CFactor_betaMP_TC2 = coherence_resonance(PxxMP_TC,freqsMP_TC); # koherens2

####################### RE için koherens rezonans (CR) ####################
# toplam memb. voltaja göre
maxindTM_RE, maxvalTM_RE = max(enumerate(PxxTM_RE), key=operator.itemgetter(1))
base1=int(maxindTM_RE*bases[0]); base2=int(maxindTM_RE*bases[1]);
if base1==base2:base1=maxindTM_RE
CFactor_betaTM_RE1 = coherence_resonance(PxxTM_RE,freqsTM_RE,base1,base2); # koherens1
CFactor_betaTM_RE2 = coherence_resonance(PxxTM_RE,freqsTM_RE); # koherens2
# toplam LFP akıma göre
maxindTP_RE, maxvalTP_RE = max(enumerate(PxxTP_RE), key=operator.itemgetter(1))
base1=int(maxindTP_RE*bases[0]); base2=int(maxindTP_RE*bases[1]);
if base1==base2:base1=maxindTP_RE
CFactor_betaTP_RE1 = coherence_resonance(PxxTP_RE,freqsTP_RE,base1,base2); # koherens1
CFactor_betaTP_RE2 = coherence_resonance(PxxTP_RE,freqsTP_RE); # koherens2
# ort. LFP akıma göre
maxindMP_RE, maxvalMP_RE = max(enumerate(PxxMP_RE), key=operator.itemgetter(1))
base1=int(maxindMP_RE*bases[0]); base2=int(maxindMP_RE*bases[1]);
if base1==base2:base1=maxindMP_RE
CFactor_betaMP_RE1 = coherence_resonance(PxxMP_RE,freqsMP_RE,base1,base2); # koherens1
CFactor_betaMP_RE2 = coherence_resonance(PxxMP_RE,freqsMP_RE); # koherens2
#---------------------------------- Son ----------------------------------

# her hücrenin ve ağın ortalama frekansı değeri
netmf_PY, netfreq_PY = network_mean_frequency(PY_spikes)
netmf_IN, netfreq_IN = network_mean_frequency(IN_spikes)
netmf_TC, netfreq_TC = network_mean_frequency(TC_spikes)
netmf_RE, netfreq_RE = network_mean_frequency(RE_spikes)

# ağın ve hücrelerin birbirleriyle ortlama faz uyumu
netmphase_PY, netphases_PY = net_mean_phase_coherence(PY_spikes,True)
netmphase_IN, netphases_IN = net_mean_phase_coherence(IN_spikes,True)
netmphase_TC, netphases_TC = net_mean_phase_coherence(TC_spikes,True)
netmphase_RE, netphases_RE = net_mean_phase_coherence(RE_spikes,True)

# ağın senkronize burst ölçütü
sync_burst_PY = Synchronous_bursting(PY_spikes)
sync_burst_IN = Synchronous_bursting(IN_spikes)
sync_burst_TC = Synchronous_bursting(TC_spikes)
sync_burst_RE = Synchronous_bursting(RE_spikes)

# iki ayrı uygulama arasındaki benzerlik ölçüsü: parametrik uzaklık
# burada örnek olsun diye aynı uygulama değeri bir sabitle çarpılmıştır
param_dist_PY = parametric_distance([netmf_PY,netmf_PY*1.5],[netmphase_PY,netmphase_PY*1.5],[sync_burst_PY,sync_burst_PY*1.5])
param_dist_IN = parametric_distance([netmf_IN,netmf_IN*1.5],[netmphase_IN,netmphase_IN*1.5],[sync_burst_IN,sync_burst_IN*1.5])
param_dist_TC = parametric_distance([netmf_TC,netmf_TC*1.5],[netmphase_TC,netmphase_TC*1.5],[sync_burst_TC,sync_burst_TC*1.5])
param_dist_RE = parametric_distance([netmf_RE,netmf_RE*1.5],[netmphase_RE,netmphase_RE*1.5],[sync_burst_RE,sync_burst_RE*1.5])

# Ağ aktivite durumu
spkCount_PY,passCount_PY = net_activity(PY_spikes)
spkCount_IN,passCount_IN = net_activity(IN_spikes)
spkCount_TC,passCount_TC = net_activity(TC_spikes)
spkCount_RE,passCount_RE = net_activity(RE_spikes)

figure()
subplot(221); pcolor(array(netphases_PY), cmap='jet',alpha=1.0);colorbar();
title('PY - Faz uyumu (PC=%s)' %(str(round(netmphase_PY,2))))
subplot(222); pcolor(array(netphases_IN), cmap='jet',alpha=1.0);colorbar();
title('IN - Faz uyumu (PC=%s)' %(str(round(netmphase_IN,2))))
subplot(223); pcolor(array(netphases_TC), cmap='jet',alpha=1.0);colorbar();
title('TC - Faz uyumu (PC=%s)' %(str(round(netmphase_TC,2))))
subplot(224); pcolor(array(netphases_RE), cmap='jet',alpha=1.0);colorbar();
title('RE - Faz uyumu (PC=%s)' %(str(round(netmphase_RE,2))))

figure(); title('Spike frekans')
plot(netfreq_PY,'bo'); plot(netfreq_IN,'ks');
plot(netfreq_TC,'g*'); plot(netfreq_RE,'r.');
xlabel('Noron ID'); ylabel('#spike'); legend(['PY','IN','TC','RE'],numpoints = 1) # scatter için scatterpoints = 1
show()

######################################################################
# ISI ile ilişkili coefﬁcient of variation - CV ve coefficient of correlation-CC
CVs_PY,CCs_PY,mCVs_PY,mCCs_PY = CV_and_CC_calc(PY_spikes,ISI_thr)
CVs_IN,CCs_IN,mCVs_IN,mCCs_IN = CV_and_CC_calc(IN_spikes,ISI_thr)
CVs_TC,CCs_TC,mCVs_TC,mCCs_TC = CV_and_CC_calc(TC_spikes,ISI_thr)
CVs_RE,CCs_RE,mCVs_RE,mCCs_RE = CV_and_CC_calc(RE_spikes,ISI_thr)

###################### popülasyonun VS parametreleri #####################
period_PY= 1./freqsTP_PY[maxindTP_PY] #1./f_ext;
period_IN= period_PY # 1./freqsTM_IN[maxindTM_IN] #1./f_ext;
period_TC= 1./freqsTP_TC[maxindTP_TC] #1./f_ext;
period_RE= period_TC # 1./freqsTM_RE[maxindTM_RE] #1./f_ext;

bins=linspace(0,2*pi,13); # 0:pi/6:2*pi, artım=pi/6=30
VS_PY,VSTh_PY,VSTheta_PY,VSX_PY,VSY_PY,VSThetaM_PY,mVS_PY=population_VS_Calc(PY_spikes,period_PY,bins)
VS_IN,VSTh_IN,VSTheta_IN,VSX_IN,VSY_IN,VSThetaM_IN,mVS_IN=population_VS_Calc(IN_spikes,period_IN,bins)
VS_TC,VSTh_TC,VSTheta_TC,VSX_TC,VSY_TC,VSThetaM_TC,mVS_TC=population_VS_Calc(TC_spikes,period_TC,bins)
VS_RE,VSTh_RE,VSTheta_RE,VSX_RE,VSY_RE,VSThetaM_RE,mVS_RE=population_VS_Calc(RE_spikes,period_RE,bins)

################## gerilim bağımlı senkronizasyon ####################
dt=defaultclock.dt;
basla=0*ms; bitir=duration-dt;
senk_PY = Voltaj_dep_senk_calc(MVs_PY.values,int(basla/dt),int(bitir/dt),smoothing=True,winlength=30,window='flat');
senk_IN = Voltaj_dep_senk_calc(MVs_IN.values,int(basla/dt),int(bitir/dt),smoothing=True,winlength=30,window='flat');
senk_TC = Voltaj_dep_senk_calc(MV_TC.values,int(basla/dt),int(bitir/dt),smoothing=True,winlength=30,window='flat');
senk_RE = Voltaj_dep_senk_calc(MV_RE.values,int(basla/dt),int(bitir/dt),smoothing=True,winlength=30,window='flat');

################## popülasyon spike ateşleme uyumluluğu (coherence) ########
dt=defaultclock.dt; bin=1*ms;
basla=0*ms; bitir=duration;

CC_PY,CC_map_PY,spike_map_PY = pop_firing_coherence(basla, bitir, bin, dt, PY_spikes)
CC_IN,CC_map_IN,spike_map_IN = pop_firing_coherence(basla, bitir, bin, dt, IN_spikes)
CC_TC,CC_map_TC,spike_map_TC = pop_firing_coherence(basla, bitir, bin, dt, TC_spikes)
CC_RE,CC_map_RE,spike_map_RE = pop_firing_coherence(basla, bitir, bin, dt, RE_spikes)

figure()
subplot(121); pcolor(CC_map_PY, cmap='jet',alpha=1.0);colorbar();
title('PY - CC map (CC=%s)' %(str(round(CC_PY,2))))
subplot(122); pcolor(CC_map_IN, cmap='jet',alpha=1.0);colorbar();
title('IN - CC map (CC=%s)' %(str(round(CC_IN,2))))

figure()
subplot(121); pcolor(CC_map_TC, cmap='jet',alpha=1.0);colorbar();
title('TC - CC map (CC=%s)' %(str(round(CC_TC,2))))
subplot(122); pcolor(CC_map_RE, cmap='jet',alpha=1.0);colorbar();
title('RE - CC map (CC=%s)' %(str(round(CC_RE,2))))

show()

##################### benzetimin genel analiz raporu ##################
AR_PY = {}; AR_IN = {}; AR_TC = {}; AR_RE = {};

AR_PY['p']=randCon_rate
AR_PY['Sig']=Noise_sigma
AR_PY['spkCount']=spkCount_PY
AR_PY['passCount']=passCount_PY
AR_PY['VSync']=round(senk_PY,3)
AR_PY['SpikeCoh']=round(CC_PY,3)
AR_PY['meanVS']=round(mVS_PY,3)
AR_PY['meanCV']=round(mCVs_PY,3)
AR_PY['meanCC']=round(mCCs_PY,3)
AR_PY['meanFreq']=round(netmf_PY,3)
AR_PY['meanPhase']=round(netmphase_PY,3)
AR_PY['SyncBurst']=round(sync_burst_PY,3)
AR_PY['StocRes']=round(snrTP_PY1,3)
AR_PY['CohRes']=round(CFactor_betaTP_PY1,3)

# PY raporu
infstr="";
infstr="Toplam Spike sayisi [PY]=%s\n\
Pasif Hucre Sayisi [PY]=%s\n\
Voltaj Senk. [PY]=%s\n\
Populasyon spike uyumu [PY]=%s\n\
Ortalama VS [PY]=%s\n\
Ortalama CV [PY]=%s\n\
Ortalama CC [PY]=%s\n\
Ort. frekans[PY]=%s\n\
Ort. faz uyumu[PY]=%s\n\
senk. burst[PY]=%s\n\
SR(SNR) toplam memb.[PY]=%s\n\
CR faktor toplam memb.[PY]=%s\n" %(str(AR_PY['spkCount']),str(AR_PY['passCount']),\
str(AR_PY['VSync']), str(AR_PY['SpikeCoh']),str(AR_PY['meanVS']),\
str(AR_PY['meanCV']),str(AR_PY['meanCC']),str(AR_PY['meanFreq']),\
str(AR_PY['meanPhase']),str(AR_PY['SyncBurst']),\
str(AR_PY['StocRes']),str(AR_PY['CohRes']))

figure()
text(0.5, 0.1, squeeze(infstr), color = 'b', visible = True, linespacing = 2, weight = 'bold', size = 'large')
xlim(0,10);ylim(0,10);

# IN raporu
AR_IN['p']=randCon_rate
AR_IN['Sig']=Noise_sigma
AR_IN['spkCount']=spkCount_IN
AR_IN['passCount']=passCount_IN
AR_IN['VSync']=round(senk_IN,3)
AR_IN['SpikeCoh']=round(CC_IN,3)
AR_IN['meanVS']=round(mVS_IN,3)
AR_IN['meanCV']=round(mCVs_IN,3)
AR_IN['meanCC']=round(mCCs_IN,3)
AR_IN['meanFreq']=round(netmf_IN,3)
AR_IN['meanPhase']=round(netmphase_IN,3)
AR_IN['SyncBurst']=round(sync_burst_IN,3)
AR_IN['StocRes']=round(snrTP_IN1,3)
AR_IN['CohRes']=round(CFactor_betaTP_IN1,3)

infstr="";
infstr="Toplam Spike sayisi [IN]=%s\n\
Pasif Hucre Sayisi [IN]=%s\n\
Voltaj Senk. [IN]=%s\n\
Populasyon spike uyumu [IN]=%s\n\
Ortalama VS [IN]=%s\n\
Ortalama CV [IN]=%s\n\
Ortalama CC [IN]=%s\n\
Ort. frekans[IN]=%s\n\
Ort. faz uyumu[IN]=%s\n\
senk. burst[IN]=%s\n\
SR(SNR) toplam memb.[IN]=%s\n\
CR faktor toplam memb.[IN]=%s\n" %(str(AR_IN['spkCount']),str(AR_IN['passCount']),\
str(AR_IN['VSync']), str(AR_IN['SpikeCoh']),str(AR_IN['meanVS']),\
str(AR_IN['meanCV']),str(AR_IN['meanCC']),str(AR_IN['meanFreq']),\
str(AR_IN['meanPhase']),str(AR_IN['SyncBurst']),\
str(AR_IN['StocRes']),str(AR_IN['CohRes']))

figure()
text(0.5, 0.1, squeeze(infstr), color = 'b', visible = True, linespacing = 2, weight = 'bold', size = 'large')
xlim(0,10);ylim(0,10);

AR_TC['p']=randCon_rate
AR_TC['Sig']=Noise_sigma
AR_TC['spkCount']=spkCount_TC
AR_TC['passCount']=passCount_TC
AR_TC['VSync']=round(senk_TC,3)
AR_TC['SpikeCoh']=round(CC_TC,3)
AR_TC['meanVS']=round(mVS_TC,3)
AR_TC['meanCV']=round(mCVs_TC,3)
AR_TC['meanCC']=round(mCCs_RE,3)
AR_TC['meanFreq']=round(netmf_TC,3)
AR_TC['meanPhase']=round(netmphase_TC,3)
AR_TC['SyncBurst']=round(sync_burst_TC,3)
AR_TC['StocRes']=round(snrTP_TC1,3)
AR_TC['CohRes']=round(CFactor_betaTP_TC1,3)

# TC raporu
infstr="";
infstr="Toplam Spike sayisi [TC]=%s\n\
Pasif Hucre Sayisi [TC]=%s\n\
Voltaj Senk. [TC]=%s\n\
Populasyon spike uyumu [TC]=%s\n\
Ortalama VS [TC]=%s\n\
Ortalama CV [TC]=%s\n\
Ortalama CC [TC]=%s\n\
Ort. frekans[TC]=%s\n\
Ort. faz uyumu[TC]=%s\n\
senk. burst[TC]=%s\n\
SR(SNR) toplam memb.[TC]=%s\n\
CR faktor toplam memb.[TC]=%s\n" %(str(AR_TC['spkCount']),str(AR_TC['passCount']),\
str(AR_TC['VSync']), str(AR_TC['SpikeCoh']),str(AR_TC['meanVS']),\
str(AR_TC['meanCV']),str(AR_TC['meanCC']),str(AR_TC['meanFreq']),\
str(AR_TC['meanPhase']),str(AR_TC['SyncBurst']),\
str(AR_TC['StocRes']),str(AR_TC['CohRes']))

figure()
text(0.5, 0.1, squeeze(infstr), color = 'b', visible = True, linespacing = 2, weight = 'bold', size = 'large')
xlim(0,10);ylim(0,10);

# RE raporu
AR_RE['p']=randCon_rate
AR_RE['Sig']=Noise_sigma
AR_RE['spkCount']=spkCount_RE
AR_RE['passCount']=passCount_RE
AR_RE['VSync']=round(senk_RE,3)
AR_RE['SpikeCoh']=round(CC_RE,3)
AR_RE['meanVS']=round(mVS_RE,3)
AR_RE['meanCV']=round(mCVs_RE,3)
AR_RE['meanCC']=round(mCCs_RE,3)
AR_RE['meanFreq']=round(netmf_RE,3)
AR_RE['meanPhase']=round(netmphase_RE,3)
AR_RE['SyncBurst']=round(sync_burst_RE,3)
AR_RE['StocRes']=round(snrTP_RE1,3)
AR_RE['CohRes']=round(CFactor_betaTP_RE1,3)

infstr="";
infstr="Toplam Spike sayisi [RE]=%s\n\
Pasif Hucre Sayisi [RE]=%s\n\
Voltaj Senk. [RE]=%s\n\
Populasyon spike uyumu [RE]=%s\n\
Ortalama VS [RE]=%s\n\
Ortalama CV [RE]=%s\n\
Ortalama CC [RE]=%s\n\
Ort. frekans[RE]=%s\n\
Ort. faz uyumu[RE]=%s\n\
senk. burst[RE]=%s\n\
SR(SNR) toplam memb.[RE]=%s\n\
CR faktor toplam memb.[RE]=%s\n" %(str(AR_RE['spkCount']),str(AR_RE['passCount']),\
str(AR_RE['VSync']), str(AR_RE['SpikeCoh']),str(AR_RE['meanVS']),\
str(AR_RE['meanCV']),str(AR_RE['meanCC']),str(AR_RE['meanFreq']),\
str(AR_RE['meanPhase']),str(AR_RE['SyncBurst']),\
str(AR_RE['StocRes']),str(AR_RE['CohRes']))

figure()
text(0.5, 0.1, squeeze(infstr), color = 'b', visible = True, linespacing = 2, weight = 'bold', size = 'large')
xlim(0,10);ylim(0,10);

show()

################################### Kayıt ################################
print "Veriler kaydediliyor..."

#recData = True

if recData: sio.savemat(path_data+pref_str,\
    mdict={'MV_PY_times':MVs_PY.times,'MV_IN_times':MVs_IN.times,\
           'MV_PY':MVs_PY.values,'MV_IN':MVs_IN.values,\
           'M_PY_spikes':M_PY.spikes,'M_IN_spikes':M_IN.spikes,\
           'M_PY_input_spikes':M_PY_input.spikes,'M_IN_input_spikes':M_IN_input.spikes,\
           'M_PY_spiketimes':M_PY.spiketimes.values(),'M_IN_spiketimes':M_IN.spiketimes.values(),\
           'M_PY_input_spiketimes':M_PY_input.spiketimes.values(),'M_IN_input_spiketimes':M_IN_input.spiketimes.values(),\
           'LFP_PY':LFP_PY.values,'LFP_IN':LFP_IN.values,\
           'LFPx_PY':LFPx_PY.values,'LFP_IN':LFPx_IN.values,\
           'M_PY_rateSR2ms':frate_PY.smooth_rate(2 * ms),'M_IN_rateSR2ms':frate_IN.smooth_rate(2 * ms),\
           
           'w_AMPA_PY_PY':w_AMPA_PY_PY,'w_AMPA_PY_IN':w_AMPA_PY_IN,\
           'w_NMDA_PY_PY':w_NMDA_PY_PY,'w_NMDA_PY_IN':w_NMDA_PY_IN,\
           'w_GABAA_IN_PY':w_GABAA_IN_PY,'w_AMPA_EXT_PY':w_AMPA_EXT_PY,\
           'w_AMPA_EXT_IN':w_AMPA_EXT_IN,\
           
           'MV_TC_times':MV_TC.times,'MV_RE_times':MV_RE.times,\
           'MV_TC':MV_TC.values,'MV_RE':MV_RE.values,\
           'M_TC_spikes':M_TC.spikes,'M_RE_spikes':M_RE.spikes,\
           'M_TC_input_spikes':M_TC_input.spikes,'M_RE_input_spikes':M_RE_input.spikes,\
           'M_TC_spiketimes':M_TC.spiketimes.values(),'M_RE_spiketimes':M_RE.spiketimes.values(),\
           'M_TC_input_spiketimes':M_TC_input.spiketimes.values(),'M_RE_input_spiketimes':M_RE_input.spiketimes.values(),\
           'LFP_TC':LFP_TC.values,'LFP_RE':LFP_RE.values,\
           'LFPx_TC':LFPx_TC.values,'LFP_RE':LFPx_RE.values,
           'M_TC_rateSR2ms':frate_TC.smooth_rate(2 * ms),'M_RE_rateSR2ms':frate_RE.smooth_rate(2 * ms),\
           
           'w_GABAA_RE_TC':w_GABAA_RE_TC,'w_GABAB_RE_TC':w_GABAB_RE_TC,\
           'w_GABAA_RE_RE':w_GABAA_RE_RE,'w_AMPA_TC_RE':w_AMPA_TC_RE,\
           'w_AMPA_EXT_TC':w_AMPA_EXT_TC,'w_AMPA_EXT_RE':w_AMPA_EXT_RE,\
           
           'w_AMPA_PY_TC':w_AMPA_PY_TC,'w_AMPA_PY_RE':w_AMPA_PY_RE,\
           'w_AMPA_TC_PY':w_AMPA_TC_PY,'w_AMPA_TC_IN':w_AMPA_TC_IN,\
           
           'N_PY':N_PY,'N_IN':N_IN,'N_TC':N_TC,'N_RE':N_RE,
           'duration':duration,'dt':dt, \
           'f_ext':f_ext,'dispersion_rate':dispersion_rate,\
           'Noise_sigma':Noise_sigma,'randCon_rate':randCon_rate,'kom_rate':kom_rate,\
           }, appendmat=True)
###########################################################################
############################# grafik ayarları #############################
#recPic = False
size1=(12,8); size2=(8,8); size3=(16,8); size4=(16,12); size5=(8,16);
dpi300=300; figformat1='.png'

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12}

matplotlib.rc('font', **font)
fs1=26; fs2=18; fs3=16

plt_i1 = 0; plt_i2PY = N_PY/2; plt_i2IN = N_IN/2; 
plt_i2TC = N_TC/2; plt_i2RE = N_RE/2;  

legend_font= FontProperties(weight='normal',size=22)
#----------------------------------------------------------------

######################### grafikler #######################################
###### PY ve IN hücrelerin zamana ve nörona göre spike frekansları #######
fig=figure(figsize=size3); 
# PY
subplot(431); plot(frate_PY.times / ms, frate_PY.smooth_rate(2 * ms),lw=1);
ylabel('PY #Spike')
subplot(432); plot(netfreq_PY,lw=1);
xlabel('PY ID'); ylabel('Frekans (Hz)');
subplot(433); maxHist = ceil(max(netfreq_PY))
a=hist(netfreq_PY,[i*maxHist/30. for i in range(31)],normed=0,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
xlabel('Frekans (Hz)'); ylabel('PY Hucre #');
# IN
subplot(434); plot(frate_IN.times / ms, frate_IN.smooth_rate(2 * ms),lw=1)
ylabel('IN #Spike');xlabel('Zaman (ms)')
subplot(435); plot(netfreq_IN,lw=1);
xlabel('IN ID'); ylabel('Frekans (Hz)');
subplot(436); maxHist = ceil(max(netfreq_IN))
a=hist(netfreq_IN,[i*maxHist/30. for i in range(31)],normed=0,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
xlabel('Frekans (Hz)'); ylabel('IN Hucre #');
# TC
subplot(437); plot(frate_TC.times / ms, frate_TC.smooth_rate(2 * ms),lw=1);
ylabel('TC #Spike')
subplot(438); plot(netfreq_TC,lw=1);
xlabel('TC ID'); ylabel('Frekans (Hz)');
subplot(439); maxHist = ceil(max(netfreq_TC))
a=hist(netfreq_TC,[i*maxHist/30. for i in range(31)],normed=0,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
xlabel('Frekans (Hz)'); ylabel('TC Hucre #');
# RE
subplot(4,3,10); plot(frate_RE.times / ms, frate_RE.smooth_rate(2 * ms),lw=1)
ylabel('RE #Spike');xlabel('Zaman (ms)')
subplot(4,3,11); plot(netfreq_RE,lw=1);
xlabel('RE ID'); ylabel('Frekans (Hz)');
subplot(4,3,12); maxHist = ceil(max(netfreq_RE))
a=hist(netfreq_RE,[i*maxHist/30. for i in range(31)],normed=0,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
xlabel('Frekans (Hz)'); ylabel('RE Hucre #');

if recPic : fig.savefig(path_images+ pref_str +'_(SpikeRateAnaliz)' + figformat1,dpi=dpi300)
###########################################################################

####################### VS Parametre Grafikleri ###########################
# herbir PY ve IN hücrenin VS değeri
fig=figure(figsize=size3)
# PY
subplot(241); #title('kortikal VS''ler')
plot(VS_PY,lw=1); xlabel('Noron ID'); ylabel('PY'); ylim(0,1)
subplot(242); maxHist = ceil(max(VS_PY))
a=hist(VS_PY,[i*maxHist/30. for i in range(31)],normed=0,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
ylabel('Frekans'); xlabel('VS'); #xlim(0,150)
# IN
subplot(243); plot(VS_IN,lw=1); xlabel('Noron ID'); ylabel('IN'); ylim(0,1)
subplot(244); maxHist = ceil(max(VS_IN))
a=hist(VS_IN,[i*maxHist/30. for i in range(31)],normed=0,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
ylabel('Frekans'); xlabel('VS'); #xlim(0,150)
# TC
subplot(245); #title('Talamik VS''ler')
plot(VS_TC,lw=1); xlabel('Noron ID'); ylabel('TC'); ylim(0,1)
subplot(246); maxHist = ceil(max(VS_TC))
a=hist(VS_TC,[i*maxHist/30. for i in range(31)],normed=0,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
ylabel('Frekans'); xlabel('VS'); #xlim(0,150)
# RE
subplot(247); plot(VS_RE,lw=1); xlabel('Noron ID'); ylabel('RE'); ylim(0,1)
subplot(248); maxHist = ceil(max(VS_RE))
a=hist(VS_RE,[i*maxHist/30. for i in range(31)],normed=0,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
ylabel('Frekans'); xlabel('VS'); #xlim(0,150)

if recPic : fig.savefig(path_images+ pref_str +'_(VSAnaliz)' + figformat1,dpi=dpi300)
#----------------------------------------------------------------

# radil VS gösterim
fig=figure(figsize=size4)
ara=10;
# PY
subplot(431,polar=True); 
title('PY'); vlen_PY=np.round(sqrt(array(list(flatten(VSX_PY))[::ara])**2+array(list(flatten(VSY_PY))[::ara])**2))
[x,y] = pol2cart(array(list(flatten(VSTheta_PY))[::ara]),vlen_PY); ax1=compass(x,y);
subplot(432); maxHist = 2*pi
a=hist(array(list(flatten(VSTheta_PY))[::]),[i*maxHist/30. for i in range(31)],normed=1,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
ylabel('Oran'); xlabel(r'$\theta$'); 
subplot(433); 
im=imshow(VSThetaM_PY, aspect='auto',origin='lower', \
    vmin=array(VSThetaM_PY).min(), vmax=array(VSThetaM_PY).max(), \
    extent=(0, 2*pi, 0, N_PY-1))
xlabel(r'$\theta$'); ylabel('PY ID'); 
cax = fig.add_axes([0.91, 0.73, 0.02, 0.172]); fig.colorbar(im,cax=cax)
# IN
subplot(434,polar=True); 
title('IN'); vlen_IN=np.round(sqrt(array(list(flatten(VSX_IN))[::ara])**2+array(list(flatten(VSY_IN))[::ara])**2))
[x,y] = pol2cart(array(list(flatten(VSTheta_IN))[::ara]),vlen_IN); ax2=compass(x,y);
subplot(435); maxHist = 2*pi
a=hist(array(list(flatten(VSTheta_IN))[::]),[i*maxHist/30 for i in range(31)],normed=1,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
ylabel('Frekans'); xlabel(r'$\theta$'); 
subplot(436); 
im=imshow(VSThetaM_IN, aspect='auto',origin='lower', \
    vmin=array(VSThetaM_IN).min(), vmax=array(VSThetaM_IN).max(), \
    extent=(0, 2*pi, 0, N_IN-1))
xlabel(r'$\theta$'); ylabel('IN ID'); 
cax = fig.add_axes([0.91, 0.52, 0.02, 0.172]); fig.colorbar(im,cax=cax)

# TC
subplot(437,polar=True); 
title('TC'); vlen_TC=np.round(sqrt(array(list(flatten(VSX_TC))[::ara])**2+array(list(flatten(VSY_TC))[::ara])**2))
[x,y] = pol2cart(array(list(flatten(VSTheta_TC))[::ara]),vlen_TC); ax1=compass(x,y);
subplot(438); maxHist = 2*pi
a=hist(array(list(flatten(VSTheta_TC))[::]),[i*maxHist/30. for i in range(31)],normed=1,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
ylabel('Oran'); xlabel(r'$\theta$'); 
subplot(439); 
im=imshow(VSThetaM_TC, aspect='auto',origin='lower', \
    vmin=array(VSThetaM_TC).min(), vmax=array(VSThetaM_TC).max(), \
    extent=(0, 2*pi, 0, N_TC-1))
xlabel(r'$\theta$'); ylabel('TC ID'); 
cax = fig.add_axes([0.91, 0.31, 0.02, 0.172]); fig.colorbar(im,cax=cax)
# RE
subplot(4,3,10,polar=True); 
title('RE'); vlen_RE=np.round(sqrt(array(list(flatten(VSX_RE))[::ara])**2+array(list(flatten(VSY_RE))[::ara])**2))
[x,y] = pol2cart(array(list(flatten(VSTheta_RE))[::ara]),vlen_RE); ax2=compass(x,y);
subplot(4,3,11); maxHist = 2*pi
a=hist(array(list(flatten(VSTheta_RE))[::]),[i*maxHist/30 for i in range(31)],normed=1,alpha=0.6,facecolor='g')
tick_params(axis='both', which='major')
ylabel('Frekans'); xlabel(r'$\theta$'); 
subplot(4,3,12); 
im=imshow(VSThetaM_RE, aspect='auto',origin='lower', \
    vmin=array(VSThetaM_RE).min(), vmax=array(VSThetaM_RE).max(), \
    extent=(0, 2*pi, 0, N_RE-1))
xlabel(r'$\theta$'); ylabel('RE ID'); 
cax = fig.add_axes([0.91, 0.104, 0.02, 0.172]); fig.colorbar(im,cax=cax)
if recPic : fig.savefig(path_images+ pref_str +'_(VSPhaseAnaliz)' + figformat1,dpi=dpi300)
###########################################################################

######################## PY ve IN hücre gerilimleri #######################
fig=figure(figsize=size3);
ax1=subplot(221); plot(MVd_PY.times/ms, MVd_PY[plt_i1]/mV)
title('[PY-Dendtritik - '+ repr(plt_i1+1) + '] - ' + repr(len(M_PY.spiketimes[plt_i1])));setp( ax1.get_xticklabels(), visible=False)
ax2=subplot(222); plot(MVd_PY.times/ms, MVd_PY[plt_i2PY]/mV)
title('[PY-Dendtritik - '+ repr(plt_i2PY+1) + '] - ' + repr(len(M_PY.spiketimes[plt_i2PY])));setp( ax2.get_xticklabels(), visible=False)
ax3=subplot(223); plot(MVs_PY.times/ms, MVs_PY[plt_i1]/mV)
title('[PY-Somatik - '+ repr(plt_i1+1) + '] - ' + repr(len(M_PY.spiketimes[plt_i1])))
ax4=subplot(224); plot(MVs_PY.times/ms, MVs_PY[plt_i2PY]/mV)
title('[PY-Somatik - '+ repr(plt_i2PY+1) + '] - ' + repr(len(M_PY.spiketimes[plt_i2PY])))
if recPic : fig.savefig(path_images+ pref_str +'_(VmPY)' + figformat1,dpi=dpi300)

fig=figure(figsize=size3);
ax1=subplot(221); plot(MVd_IN.times/ms, MVd_IN[plt_i1]/mV)
title('[IN-Dendtritik - '+ repr(plt_i1+1) + '] - ' + repr(len(M_IN.spiketimes[plt_i1])));setp( ax1.get_xticklabels(), visible=False)
ax2=subplot(222); plot(MVd_IN.times/ms, MVd_IN[plt_i2IN]/mV)
title('[IN-Dendtritik - '+ repr(plt_i2IN+1) + '] - ' + repr(len(M_IN.spiketimes[plt_i2IN])));setp( ax2.get_xticklabels(), visible=False)
ax3=subplot(223); plot(MVs_IN.times/ms, MVs_IN[plt_i1]/mV)
title('[IN-Somatik - '+ repr(plt_i1+1) + '] - ' + repr(len(M_IN.spiketimes[plt_i1])))
ax4=subplot(224); plot(MVs_IN.times/ms, MVs_IN[plt_i2IN]/mV)
title('[IN-Somatik - '+ repr(plt_i2IN+1) + '] - ' + repr(len(M_IN.spiketimes[plt_i2IN])))
if recPic : fig.savefig(path_images+ pref_str +'_(VmIN)' + figformat1,dpi=dpi300)

######################## TC ve RE hücre gerilimleri #######################
fig=figure(figsize=size2);
ax1=subplot(211); plot(MV_TC.times/ms, MV_TC[plt_i1]/mV)
title('[TC - '+ repr(plt_i1+1) + '] - ' + repr(len(M_TC.spiketimes[plt_i1])));
ax2=subplot(212); plot(MV_TC.times/ms, MV_TC[plt_i2TC]/mV)
title('[TC - '+ repr(plt_i2TC+1) + '] - ' + repr(len(M_TC.spiketimes[plt_i2TC])));
if recPic : fig.savefig(path_images+ pref_str +'_(VmTC)' + figformat1,dpi=dpi300)

fig=figure(figsize=size2);
ax1=subplot(211); plot(MV_RE.times/ms, MV_RE[plt_i1]/mV)
title('[RE - '+ repr(plt_i1+1) + '] - ' + repr(len(M_RE.spiketimes[plt_i1])))
ax2=subplot(212); plot(MV_RE.times/ms, MV_RE[plt_i2RE]/mV)
title('[RE - '+ repr(plt_i2RE+1) + '] - ' + repr(len(M_RE.spiketimes[plt_i2RE])))
if recPic : fig.savefig(path_images+ pref_str +'_(VmRE)' + figformat1,dpi=dpi300)
########################################################################

######################### raster plot ##################################
fig=figure(figsize=size1);
#my_raster_plot(M_IN, M_PY, M_PY_input,colorlist=['k.','b.','r.'],xlabe='Zaman (ms)',ylabel='Grup',yticks = ['IN', 'PY', 'EXT'])
#my_raster_plot(M_RE, M_TC, M_TC_input,colorlist=['k.','b.','r.'],xlabe='Zaman (ms)',ylabel='Grup',yticks = ['RE', 'TC', 'EXT'])
my_raster_plot(M_RE,M_TC,M_IN,M_PY,colorlist=['k.','b.','r.','g.'],xlabe='Zaman (ms)',ylabel='Grup',yticks = ['RE','TC','IN','PY'])
if recPic : fig.savefig(path_images+ pref_str +'_(raster)' + figformat1,dpi=dpi300)
#----------------------------------------------------------------

######################### Vm color plot full ########################
fig=figure(figsize=size1)
# PY
ax1=subplot(411); 
im=imshow(MVs_PY[:,:]/mV,  aspect='auto',origin='lower',vmin=(MVs_PY.getvalues()/mV).min(), vmax=(MVs_PY.getvalues()/mV).max(),\
        extent=(0, duration/ms, 0, N_PY-1))
ylabel('PYs - Noron ID'); title('Soma Potansiyel (mV)')
setp(ax1.get_xticklabels(), visible=False)
cax = fig.add_axes([0.91, 0.73, 0.02, 0.172])
fig.colorbar(im,cax=cax)
# IN
ax2=subplot(412); 
im=imshow(MVs_IN[:,:]/mV, aspect='auto',origin='lower',vmin=(MVs_IN.getvalues()/mV).min(), vmax=(MVs_IN.getvalues()/mV).max(),\
        extent=(0, duration/ms, 0, N_IN-1))
ylabel('IN - Noron ID'); xlabel('Zaman (ms)')
cax = fig.add_axes([0.91, 0.52, 0.02, 0.172])
fig.colorbar(im,cax=cax)
# TC
ax3=subplot(413); 
im=imshow(MV_TC[:,:]/mV, aspect='auto',origin='lower',vmin=(MV_TC.getvalues()/mV).min(), vmax=(MV_TC.getvalues()/mV).max(),\
        extent=(0, duration/ms, 0, N_TC-1))
ylabel('TC - Noron ID'); title('Potansiyel (mV)')
setp(ax1.get_xticklabels(), visible=False)
cax = fig.add_axes([0.91, 0.31, 0.02, 0.172])
fig.colorbar(im,cax=cax)
# RE
ax4=subplot(414); 
im=imshow(MV_RE[:,:]/mV, aspect='auto',origin='lower',vmin=(MV_RE.getvalues()/mV).min(), vmax=(MV_RE.getvalues()/mV).max(),\
        extent=(0, duration/ms, 0, N_RE-1))
ylabel('RE - Noron ID'); xlabel('Zaman (ms)')
cax = fig.add_axes([0.91, 0.104, 0.02, 0.172])
fig.colorbar(im,cax=cax)
if recPic : fig.savefig(path_images+ pref_str +'_(VmColorFull)' + figformat1,dpi=dpi300)
#----------------------------------------------------------------
#Vm color plot 700ms 
fig=figure(figsize=size1)
# PY
ax1=subplot(411); 
im=imshow(MVs_PY[:,:]/mV, aspect='auto',origin='lower',vmin=(MVs_PY.getvalues()/mV).min(), vmax=(MVs_PY.getvalues()/mV).max(),\
        extent=(0, duration/ms, 0, N_PY-1))
ylabel('PY - Noron ID'); title('Soma Potansiyel (mV)')
xlim((50*ms)/ms,(1000*ms)/ms);setp(ax1.get_xticklabels(), visible=False)
cax = fig.add_axes([0.91, 0.73, 0.02, 0.172])
fig.colorbar(im,cax=cax)
# IN
ax2=subplot(412); 
im=imshow(MVs_IN[:,:]/mV, aspect='auto',origin='lower',vmin=(MVs_IN.getvalues()/mV).min(), vmax=(MVs_IN.getvalues()/mV).max(),\
        extent=(0, duration/ms, 0, N_IN-1))
ylabel('IN - Noron ID'); xlabel('Zaman (ms)')
xlim((50*ms)/ms,(1000*ms)/ms)
cax = fig.add_axes([0.91, 0.52, 0.02, 0.172])
fig.colorbar(im,cax=cax)
# TC
ax1=subplot(413); 
im=imshow(MV_TC[:,:]/mV, aspect='auto',origin='lower',vmin=(MV_TC.getvalues()/mV).min(), vmax=(MV_TC.getvalues()/mV).max(),\
        extent=(0, duration/ms, 0, N_TC-1))
ylabel('TC - Noron ID'); title('Potansiyel (mV)')
xlim((50*ms)/ms,(1000*ms)/ms);setp(ax1.get_xticklabels(), visible=False)
cax = fig.add_axes([0.91, 0.31, 0.02, 0.172])
fig.colorbar(im,cax=cax)
# RE
ax2=subplot(414); 
im=imshow(MV_RE[:,:]/mV, aspect='auto',origin='lower',vmin=(MV_RE.getvalues()/mV).min(), vmax=(MV_RE.getvalues()/mV).max(),\
        extent=(0, duration/ms, 0, N_RE-1))
ylabel('RE - Noron ID'); xlabel('Zaman (ms)')
xlim((50*ms)/ms,(1000*ms)/ms)
cax = fig.add_axes([0.91, 0.104, 0.02, 0.172])
fig.colorbar(im,cax=cax)
if recPic : fig.savefig(path_images+ pref_str +'_(VmColorSection)' + figformat1,dpi=dpi300)
#----------------------------------------------------------------

# Kortikal sinaptik bağlar
fig=figure(figsize=size1)
ax1=subplot(321)
imshow(AMPA_PY_PY.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('PY');ylabel('PY');title('AMPA PY-PY');
ax2=subplot(322)
imshow(AMPA_PY_IN.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('IN');ylabel('PY');title('AMPA PY-IN');
ax3=subplot(323)
imshow(NMDA_PY_PY.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('PY');ylabel('PY');title('NMDA PY-PY')
ax4=subplot(324)
imshow(NMDA_PY_IN.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('IN');ylabel('PY');title('NMDA_PY_IN');
ax5=subplot(325)
imshow(GABAA_IN_PY.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('PY');ylabel('IN');title('GABAA IN-PY');
if recPic : fig.savefig(path_images+ pref_str +'_(SynMatrixKortikal)' + figformat1,dpi=dpi300)

# Talamik sinaptik bağlar
fig=figure(figsize=size2)
ax1=subplot(221)
imshow(GABAA_RE_TC.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('TC');ylabel('RE');title('GABAA RE-TC');
ax2=subplot(222)
imshow(GABAB_RE_TC.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('TC');ylabel('RE');title('GABAB RE-TC');
setp(ax2.get_xticklabels(), visible=False);
ax3=subplot(223)
imshow(GABAA_RE_RE.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('RE');ylabel('RE');title('GABAA RE-RE')
ax4=subplot(224)
imshow(AMPA_TC_RE.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('RE');ylabel('TC');title('AMPA TC-RE');
if recPic : fig.savefig(path_images+ pref_str +'_(SynMatrixTalamik)' + figformat1,dpi=dpi300)

# Talamik sinaptik bağlar
fig=figure(figsize=size2)
ax1=subplot(221)
imshow(AMPA_PY_TC.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('TC');ylabel('PY');title('AMPA PY-TC');
ax2=subplot(222)
imshow(AMPA_PY_RE.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('RE');ylabel('PY');title('AMPA PY-RE');
setp(ax2.get_xticklabels(), visible=False);
ax3=subplot(223)
imshow(AMPA_TC_PY.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('PY');ylabel('TC');title('AMPA TC-PY')
ax4=subplot(224)
imshow(AMPA_TC_IN.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('IN');ylabel('TC');title('AMPA TC-IN');
if recPic : fig.savefig(path_images+ pref_str +'_(SynMatrixTalamoKortikal)' + figformat1,dpi=dpi300)
#----------------------------------------------------------------

show()