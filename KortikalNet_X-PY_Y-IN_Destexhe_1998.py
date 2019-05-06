# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 18:28:18 2013

@author: ramazan tekin

PY hücre modeli 
ref : Destexhe A, Contreras D, Steriade M (1998a) Mechanisms underlying the 
synchronizing action of corticothalamic feedback through inhibition of 
thalamic relay cells. J Neurophysiol 79:999–1016. 

fast sodyum INa ve fast potasyum IK (Traub and Miles, 1991), 
Yavaş deaktive olmayan K akım (McCormick ve ark. 1993)
low-threshold Ca+2 IT (Destexhe ve ark.  LTS cells in cerebral cortex (2001)), 
high-threshold Ca+2 IL (Reuveni ve ark. 1993;Pospischil ve ark. 2008)
akımları içermektedir (Destexhe ve ark. (1998a,2001);Pospischil ve ark. 2008).

IN hücre modeli 
ref : Destexhe A, Contreras D, Steriade M (1998a) Mechanisms underlying the 
synchronizing action of corticothalamic feedback through inhibition of 
thalamic relay cells. J Neurophysiol 79:999–1016. 

fast sodyum INa ve fast potasyum IK (Traub and Miles, 1991), 
akımları içermektedir (Destexhe A, Contreras D, Steriade M (1998a)).

"""

from brian import *

# Parameters
#defaultclock.dt=0.05*ms
duration = 500 * ms # 4*1000 * ms

# Nernst parameters
R=8.31451 * joule/(mole*kelvin)
F=96.485 * kcoulomb/mole
Temperature = (36 + 273.15) * kelvin
valence = 2 # valance is +1 for Na+, +1 for K+, +2 for Ca2+, -1 for Cl-, 

Nernst = lambda valence, Temperature, InConc, OutConc:\
    R * (Temperature) / (F * valence) * \
    log(OutConc / InConc) 

# PY sıcaklık ölçekleri
q10 = 2.5; tadj_PY= q10 ** ((36. - 24.)/ 10.)

# PY ve IN hücre sayıları
N_PY = 2*50; N_IN = 2*50;

# L=d=96 um, area=pi*L*d
area_PY = 29000 * umetre ** 2 # yada 2.9e-4*(cm**2)
# L=64.86um ve d=70 um, area=pi*L*d
area_IN = 14300 * umetre ** 2 # yada 1.43e-4*(cm**2)

Cm_PY = (1 * ufarad * cm ** -2) * area_PY
Cm_IN = (1 * ufarad * cm ** -2) * area_IN

# PY reversaller
EL_PY = -70 * mV
EKL_PY = -100 * mV
EK_PY = -100 * mV
ENa_PY = 50 * mV

#PY iletkenlikler
g_Lk_PY = (0.01 * msiemens * cm ** -2) * area_PY # RS için 0.1 mS/cm2
g_KL_PY = (0.005 * msiemens * cm ** -2) * area_PY 
g_Na_PY = (50. * msiemens * cm ** -2) * area_PY
g_K_PY = (5. * msiemens * cm ** -2) * area_PY
g_M_PY = (0.03 * msiemens * cm ** -2) * area_PY # RS için 0.07 mS
g_T_PY = (0.8 * msiemens * cm ** -2) * area_PY # 0.8 mS Destexhe ve ark. 2001

g_L_PY = (0.2 * msiemens * cm ** -2) * area_PY
P_T_PY= (0.2e-3 * cm/second) * area_PY

#IN reversaller
EL_IN = -70 * mV
EKL_IN = -70 * mV
EK_IN = -100 * mV
ENa_IN = 50 * mV

# IN iletkenlikler
g_Lk_IN = (0.15 * msiemens * cm ** -2) * area_IN
g_Na_IN = (50. * msiemens * cm ** -2) * area_IN
g_K_IN = (10. * msiemens * cm ** -2) * area_IN

# [Ca]_dış konsantrasyonu
Ca_o = 2.0 * (1e-6 * mole*cm**-3) # 1 M = 1e-3 mole/cm3
# [Ca]_dış kararlı durum
Ca_i_inf = 2.4e-4 * (1e-6 * mole * cm**-3) # 1 M = 1e-3 mole/cm3

A = (5.18e-5*(1e-6 * mole * cm**-3)*cm2/(ms*uA))

depth	= .1 * um

# M-akım zaman sabiti (tau_p) tepe değeri
tau_max_PY = 1000 * ms

VT_PY = (-55) * mV
VT_IN = (-55) * mV

# zaman sabitleri
tau_AMPA_PY_PY = 10 * ms
tau_AMPA_PY_IN = 10 * ms
tau_GABAA_IN_IN = 10 * ms
tau_GABAA_IN_PY = 10 * ms
tau_GABAB_IN_IN = 10 * ms
tau_GABAB_IN_PY = 10 * ms

# sinaptik Reversal potensiyeller
E_AMPA_PY_PY = 0 * mV
E_GABAA_IN_PY = -80 * mV
E_GABAB_IN_PY = -80 * mV
E_AMPA_PY_IN = 0 * mV
E_GABAA_IN_IN = -80 * mV
E_GABAB_IN_IN = -80 * mV

w_AMPA_PY_PY = 0.6 * uS # 0-0.9 uS Destexhe ve ark. 1998a
w_GABAA_IN_PY = 0.15 * uS # 0.09-0.2 uS Destexhe ve ark. 1998a
w_GABAB_IN_PY = 0.03 * uS # 0-0.2 uS Destexhe ve ark. 1998a
w_AMPA_PY_IN = 0.2 * uS # 0.1-0.4 uS Destexhe ve ark. 1998a
#w_GABAA_IN_IN = 0.2 * uS 
#w_GABAB_IN_IN = 0.06 * uS


################################################################
# PY T-akım kinetikleri
# (Destexhe ve ark. 1998c) 
mCainf_PY=lambda v:(1.)/(1+exp((-56*mV-v)/(6.2*mV)))
hCainf_PY=lambda v:(1.)/(1+exp((80*mV+v)/(4*mV)))

# Huguenard and McCormick (1992)
#taumCa_PY=lambda v:((0.612*ms) + (1.*ms)/(exp((-131*mV-v)/(16.7*mV))+\
#                exp((15.8*mV+v)/(18.2*mV))))/tadj_PY
# yada Bazhenov ve ark. 1998a - cellular and network models ...
taumCa_PY=lambda v:(0.13*ms) + (0.22*ms)/(exp((-132*mV-v)/(16.7*mV))+\
                exp((16.8*mV+v)/(18.2*mV)))

# Huguenard and McCormick (1992)        
#tauhCa_PY=lambda v:(28.0*ms + (1.*ms)*exp((-21*mV-v)/(10.5*mV)))/tadj_PY if (v>=-81*mV) else \
#                (1.*ms)*exp((466*mV+v)/(66.6*mV))/tadj_PY
# yada Bazhenov ve ark. 1998a - cellular and network models ...
tauhCa_PY=lambda v:(8.2*ms + (1.*ms)*((56.6 + 0.27*exp((115.2*mV+v)/(5*mV)))/\
                (1+exp((86*mV+v)/(3.2*mV)))))
# yada Destexhe ve ark. 1996a Ionic mechanisms underlying
## ****** bu modelde aktivasyonun kararlı durumu bulunuyor mCa_inf *********
#tauhCa_PY=lambda v:(30.8*ms + (1.*ms)*((211.4 + exp((115.2*mV+v)/(5*mV)))/\
#                (1+exp((86*mV+v)/(3.2*mV)))))/tadj_TC

################################################################

################################################################
# M-akım
# McCormick ve ark. (1993)
# http://senselab.med.yale.edu/modeldb/showmodel.asp?model=3817&file=\cortex\IM.mod
pinf_PY=lambda v:(1.)/(1+exp((-35*mV-v)/(10.0*mV)))
# McCormick ve ark. (1993)
taup_PY=lambda v:(tau_max_PY)/(3.3*exp((35.0*mV+v)/(20.0*mV))+exp((-35.0*mV-v)/(20.0*mV)))

################################################################

################################################################
# PY harici uyartı
i_cur1_PY = 0*0.7; i_cur2_PY = 0*-0.1
eqs_input_PY='''
B = x*nA : amp
x = 0 if y<50 else i_cur1_PY if y<150 else 0 if y<200 else i_cur2_PY if y<220 else 0 : 1
dy/dt = (1)/ms : 1
'''
input_PY=NeuronGroup(1,eqs_input_PY)
################################################################

################################################################
# IN harici uyartı
i_cur1_IN = 0*0.1; i_cur2_IN = 0*-0.1
eqs_input_IN='''
B = x*nA : amp
x = 0 if y<10 else i_cur1_IN if y<50 else 0 if y<200 else i_cur2_IN if y<220 else 0 : 1
dy/dt = (1)/ms : 1
'''
input_IN=NeuronGroup(1,eqs_input_IN)
################################################################

vtrap=lambda x,y : (y*(1-x/y/2)) if abs(x/y)<1e-6 else (x/(Exp(x/y)-1))
Exp=lambda x : 0 if x<-100 else exp(x)

def GHK(v, Cai, Cao) : 
    z = 2*(F)*v/(R*T)
    ghk=2*(F)*(Cai*(-z/(exp(-z) - 1)) - Cao*z/(exp(z) - 1))
    return ghk


print "Model oluşturuluyor..."

# PY nöron modeli
# Destexhe A, Contreras D, Steriade M (1998a)
eqs_PY = Equations('''
dv/dt = (I+I_syn-I_Lk-0*I_KL-I_Na-I_K-I_T-I_M)/Cm_PY : volt

I_syn = I_AMPA_PY_PY+I_GABAA_IN_PY+I_GABAB_IN_PY : amp
I_AMPA_PY_PY = g_AMPA_PY_PY*(E_AMPA_PY_PY-v) : amp
I_GABAA_IN_PY = g_GABAA_IN_PY*(E_GABAA_IN_PY-v) : amp
I_GABAB_IN_PY = g_GABAB_IN_PY*(E_GABAB_IN_PY-v) : amp

I_Lk = g_Lk_PY*(v-EL_PY) : amp
I_KL = g_KL_PY*(v-EKL_PY) : amp
I_Na = g_Na_PY*(m*m*m)*h*(v-ENa_PY) : amp
I_K = g_K_PY*(n*n*n*n)*(v-EK_PY) : amp
I_T = g_T_PY*(mCa**2)*hCa*(v-ECa_PY) : amp
I_M = g_M_PY*p*(v-EK_PY) : amp

#############################
# M-akım kinetikleri 
dp/dt=(p_inf-p)/tau_p : 1
p_inf=pinf_PY(v) : 1   
tau_p=taup_PY(v) : second

#############################

#############################
# Ca_i dinamiği
dCa_i/dt = -(A*(I_T/area_PY))+(Ca_i_inf-Ca_i)/(5*ms) : mole/cm3
#dCa_i/dt = -(10*I_T/area_PY)/(2*F*depth)+(Ca_i_inf-Ca_i)/(5*ms) : mole/cm3
#dCa_i/dt = -(1e4*I_T/mA)/(2*96489*1*ms)+((2.4e-4)-Ca_i)/(5*ms) : 1
#############################

#############################
# Kalsiyum reversal potansiyel dinamiği
dECa_PY/dt = (Nernst(valence,Temperature, Ca_i,Ca_o)-ECa_PY)/(ms) : volt
#############################

#############################
# I_T akım kinetikleri
dmCa/dt=(mCa_inf-mCa)/tau_mCa : 1
mCa_inf=mCainf_PY(v) : 1   
tau_mCa=taumCa_PY(v) : second
dhCa/dt=(hCa_inf-hCa)/tau_hCa : 1
hCa_inf=hCainf_PY(v) : 1   
tau_hCa=tauhCa_PY(v) : second
#############################

#############################
# fast sodyum ve potasyum akım kinetikleri
dm/dt = alpham*(1-m)-betam*m : 1
dn/dt = alphan*(1-n)-betan*n : 1
dh/dt = alphah*(1-h)-betah*h : 1

alpham = 0.32*(mV**-1)*(13*mV-v+VT_PY)/ \
    (exp((13*mV-v+VT_PY)/(4*mV))-1.)/ms : Hz

betam = 0.28*(mV**-1)*(v-VT_PY-40*mV)/ \
    (exp((v-VT_PY-40*mV)/(5*mV))-1)/ms : Hz

alphah = 0.128*exp((17*mV-v+VT_PY)/(18*mV))/ms : Hz
betah = 4./(1+exp((40*mV-v+VT_PY)/(5*mV)))/ms : Hz

alphan = 0.032*(mV**-1)*(15*mV-v+VT_PY)/ \
    (exp((15*mV-v+VT_PY)/(5*mV))-1.)/ms : Hz
betan = .5*exp((10*mV-v+VT_PY)/(40*mV))/ms : Hz
#############################

g_AMPA_PY_PY : siemens
g_GABAA_IN_PY : siemens
g_GABAB_IN_PY : siemens

I : amp
''')

# IN nöron modeli
# Destexhe A, Contreras D, Steriade M (1998a)
eqs_IN = Equations('''
dv/dt = (I+I_syn-I_Lk-I_Na-I_K)/Cm_IN : volt

I_syn = I_AMPA_PY_IN+I_GABAA_IN_IN + I_GABAB_IN_IN : amp
I_AMPA_PY_IN = g_AMPA_PY_IN*(E_AMPA_PY_IN-v) : amp
I_GABAA_IN_IN = g_GABAA_IN_IN*(E_GABAA_IN_IN-v) : amp
I_GABAB_IN_IN = g_GABAB_IN_IN*(E_GABAB_IN_IN-v) : amp

I_Lk = g_Lk_IN*(v-EL_IN) : amp
I_Na = g_Na_IN*(m*m*m)*h*(v-ENa_IN) : amp
I_K = g_K_IN*(n*n*n*n)*(v-EK_IN) : amp

#############################
# fast sodyum ve potasyum akım kinetikleri
dm/dt = alpham*(1-m)-betam*m : 1
dn/dt = alphan*(1-n)-betan*n : 1
dh/dt = alphah*(1-h)-betah*h : 1

alpham = 0.32*(mV**-1)*(13*mV-v+VT_IN)/ \
    (exp((13*mV-v+VT_IN)/(4*mV))-1.)/ms : Hz

betam = 0.28*(mV**-1)*(v-VT_IN-40*mV)/ \
    (exp((v-VT_IN-40*mV)/(5*mV))-1)/ms : Hz

alphah = 0.128*exp((17*mV-v+VT_IN)/(18*mV))/ms : Hz
betah = 4./(1+exp((40*mV-v+VT_IN)/(5*mV)))/ms : Hz

alphan = 0.032*(mV**-1)*(15*mV-v+VT_IN)/ \
    (exp((15*mV-v+VT_IN)/(5*mV))-1.)/ms : Hz
betan = .5*exp((10*mV-v+VT_IN)/(40*mV))/ms : Hz
#############################

dg_GABAA_IN_IN/dt = -g_GABAA_IN_IN*(1./tau_GABAA_IN_IN) : siemens
dg_GABAB_IN_IN/dt = -g_GABAB_IN_IN*(1./tau_GABAB_IN_IN) : siemens
g_AMPA_PY_IN : siemens

I : amp
''')


PY = NeuronGroup(N_PY, model=eqs_PY,
    threshold=EmpiricalThreshold(threshold= -0 * mV, refractory=3 * ms),
    implicit=True, freeze=True)

IN = NeuronGroup(N_PY, model=eqs_IN,
    threshold=EmpiricalThreshold(threshold= -0 * mV, refractory=3 * ms),
    implicit=True, freeze=True)

print "Bağlantılar oluşturuluyor..."

######################### AMPA Sinaps ######################################
# AMPA sinaptik sabitler
Tmax_AMPA   = 0.5 * (1e-6 * mole * cm**-3)               #  maks. trasmiter konsantrasyonu
Cdur_AMPA   = 0.3 * ms                  #  trasmiter süresi (rising phase)
Alpha_AMPA  = 0.94 / (ms * (1e-6 * mole * cm**-3))       #  bağlanma (binding) oranı
Beta_AMPA   = 0.18 / ms                 #  serbestlik (unbinding) oranı

# AMPA sinaps modeli : 
# Destexhe ve ark. 1998b - 
# Bazhenov ve ark. 1998a Cellular and network models ...

AMPA_PY_PY = Synapses(PY, PY,
             model='''w : siemens
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      T = (Tmax_AMPA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : mole/cm3
                      gampa = w*x_ampa : siemens
                      ''')

AMPA_PY_IN = Synapses(PY, IN,
             model='''w : siemens
                      dx_ampa/dt = (Alpha_AMPA*T)*(1-x_ampa)-Beta_AMPA*x_ampa : 1
                      T = (Tmax_AMPA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : mole/cm3
                      gampa = w*x_ampa : siemens
                      ''')
######################### AMPA Son #########################################

######################### NMDA Sinaps #####################################
# GABAA sinaptik sabitler
Tmax_NMDA   = 1.0 * mmole           #  maks. trasmiter konsantrasyonu
Cdur_NMDA   = 1 * ms                #  trasmiter süresi (rising phase)
Alpha_NMDA  = 0.072 / (ms * mmole)    #  bağlanma (binding) oranı
Beta_NMDA   = 0.0066 / ms             #  serbestlik (unbinding) oranı

mg_NMDA = 1 * mmole                 # [Mg]_o = 1 - 2 mM

# AMPA sinaps modeli : 
# Destexhe ve ark. 1998b - 
# Bazhenov ve ark. 1998a Cellular and network models ...
S2 = Synapses(PY, IN,
             model='''w : siemens
                      dx_nmda/dt = (Alpha_NMDA*T)*(1-x_nmda)-Beta_NMDA*x_nmda : 1
                      T = (Tmax_NMDA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : mole
                      mgblock = (1.0 / (1.0 + exp(-0.062 * (v_post/mV)) * (mg_NMDA / (3.57 * mmole)))) : 1
                      gnmda = w*x_nmda*mgblock : siemens
                      ''')

######################### NMDA Son ########################################

######################### GABAA Sinaps #####################################
# GABAA sinaptik sabitler
Tmax_GABAA   = 0.5 * (1e-6 * mole * cm**-3)           #  maks. trasmiter konsantrasyonu
Cdur_GABAA   = 0.3 * ms                #  trasmiter süresi (rising phase)
Alpha_GABAA  = 20.0 / (ms * (1e-6 * mole * cm**-3))    #  bağlanma (binding) oranı
Beta_GABAA   = 0.162 / ms             #  serbestlik (unbinding) oranı

# GABAA sinaps modeli : 
# Destexhe ve ark. 1998b - 
# Bazhenov ve ark. 1998a Cellular and network models ...

GABAA_IN_PY = Synapses(IN, PY,
              model='''w : siemens
                      dx_gabaa/dt = (Alpha_GABAA*T)*(1-x_gabaa)-Beta_GABAA*x_gabaa : 1
                      T = (Tmax_GABAA)/(1+exp(-(v_pre-2*mV)/(5*mV))) : mole/cm3
                      ggabaa = w*x_gabaa : siemens
                      ''')

######################### GABAA Son ########################################

######################### GABAB Sinaps #####################################
# GABAB sinaptik sabitler
Tmax_GABAB   = 0.5 * (1e-6 * mole * cm**-3)           #  maks. trasmiter konsantrasyonu
K1_GABAB    = 0.5 / (ms * (1e-6 * mole * cm**-3))     #  reseptörlere bağlanma oranı
K2_GABAB    = 0.0012 / ms            #  reseptörlerin serbest kalma oranı
K3_GABAB    = 0.18 / ms             #  G-protein oluşma oranı
K4_GABAB    = 0.034 / ms             #  G-protein yok olma oranı
KD_GABAB    = 100                    #  K+ kanalı ayrışma sabiti (umole**4)

n_GABAB     = 4             # G-proteininin K+ kanalına bağlanma nokta sayısı

# AMPA sinaps modeli : 
# Destexhe ve ark. 1998b - 
# Bazhenov ve ark. 1998a Cellular and network models ...
GABAB_IN_PY = Synapses(IN, PY,
              model='''w : siemens
                      dx_gabab/dt = (K1_GABAB*T)*(1-x_gabab)-K2_GABAB*x_gabab : 1
                      ds_gabab/dt = (K3_GABAB*x_gabab)-(K4_GABAB*s_gabab) : 1
                      T = (Tmax_GABAB)/(1+exp(-(v_pre-2*mV)/(5*mV))) : mole/cm3
                      ggabab = w*(s_gabab**n_GABAB)/((s_gabab**n_GABAB)+KD_GABAB) : siemens
                      ''')
######################### GABAB Son ########################################

def shift_left(seq, n):
    n = n % len(seq)
    return seq[n:] + seq[:n]

def shift_right(seq, n):
    n = n % len(seq)
    return seq[-n:] + seq[:-n]

kom_PY = N_PY/10; kom_IN = N_IN/10 # 10% sağ, 10% sol toplam %20+1
indexPY = arange(0,N_PY); indexIN = arange(0,N_IN);
bindexPY = (ones((1,N_PY),bool)); bindexIN = (ones((1,N_IN),bool))
bindexPY[0,arange(kom_PY+1,N_PY-kom_PY)] = False
bindexIN[0,arange(kom_IN+1,N_IN-kom_IN)] = False
bindex2PY = bindexPY.copy(); bindex2PY[0,0] = False
bindex2IN = bindexIN.copy(); bindex2IN[0,0] = False

# PY->PY AMPA bağ
PY.g_AMPA_PY_PY = AMPA_PY_PY.gampa
for i in range(N_PY):
    AMPA_PY_PY[indexPY[array(shift_right(list(bindex2PY[0]),i))],i]=True
AMPA_PY_PY.w=w_AMPA_PY_PY/N_PY

# PY->IN AMPA bağ
IN.g_AMPA_PY_IN = AMPA_PY_IN.gampa
for i in range(N_IN):
    AMPA_PY_IN[indexPY[array(shift_right(list(bindexPY[0]),i))],i]=True
AMPA_PY_IN.w=w_AMPA_PY_IN/N_PY

# IN->PY GABAA bağ
PY.g_GABAA_IN_PY = GABAA_IN_PY.ggabaa
for i in range(N_PY):
    GABAA_IN_PY[indexIN[array(shift_right(list(bindexIN[0]),i))],i]=True
GABAA_IN_PY.w=w_GABAA_IN_PY/N_IN

# IN->PY GABAB bağ
PY.g_GABAB_IN_PY = GABAB_IN_PY.ggabab
for i in range(N_PY):
    GABAB_IN_PY[indexIN[array(shift_right(list(bindexIN[0]),i))],i]=True
GABAB_IN_PY.w=w_GABAB_IN_PY/N_IN


######################### PY spike giriş ###################################
# PY hücre için spike girişi
"""
Nspikes = 10;T=15*ms;T0=T*Nspikes;
PY_input = SpikeGeneratorGroup(1,[(0,n*T+T0) for n in range(Nspikes)])
C_PY_input=Connection(PY_input,PY,'g_AMPA_PY_PY',weight=0.5*uS)
"""
# PoissonInput birbirinden bağımsız giriş sayısı kadar poison girişi state değişkenine
# weight ağırlığa göre ekler
# deta için http://briansimulator.org/docs/inputs.html?highlight=poissoninput
# C_PY_input = PoissonInput(PY, bgs, rate=20*Hz, weight=0.5*uS, state='g_AMPA_PY_PY')
# poison dağılıma göre 
#PY_input = PoissonGroup(1,rates=20*Hz)
# Gaussian dağılıma göre
#PY_input=PulsePacket(t=100*ms,n=10,sigma=30*ms)
#C_PY_input=Connection(PY_input,PY,'g_AMPA_PY_PY',weight=0.5*uS)
# elle ayarlanmış düzenli spike dizisi
spiketimes = [(0,100*ms), (0,200*ms), (0,300*ms), (0,400*ms), (0,500*ms),\
              (0,600*ms), (0,700*ms), (0,800*ms), (0,900*ms)]
#PY_input = SpikeGeneratorGroup(1,spiketimes)
#IN_input = SpikeGeneratorGroup(1,spiketimes)

#C_PY_input=Connection(PY_input,PY,'g_AMPA_PY_PY',weight=lambda i,j:exp(-abs(j-i-N_PY/2)*.2)*5*uS)
#C_IN_input=Connection(IN_input,IN,'g_AMPA_PY_IN',weight=lambda i,j:exp(-abs(j-i-N_IN/2)*.2)*5*uS)

######################### PY spike giriş son ###############################


print "Gerekli hazırlıklar yapılıyor..."

# PY nöron hazırlık
#PY[1].g_T_PY = (0.01 * msiemens * cm ** -2) * area_PY

#PY.I=linked_var(input_PY,'B')
PY.Ca_i = 2.4e-4 * (1e-6 * mole * cm**-3) # 1 M = 1e-3 mole/cm3 
PY.ECa_PY = 120 * mV
PY.v = -80*mV

# IN nöron hazırlık
#IN.I=linked_var(input_IN,'B')
IN.v = -70*mV

# PY hücre Akım veya kinetik kayıtları
M_PY = SpikeMonitor(PY)
MV_PY = StateMonitor(PY, 'v', record=True)
MI_PY = StateMonitor(PY,'I', record=True)
MI_AMPA_PY_PY = StateMonitor(PY, 'I_AMPA_PY_PY', record=True)
MI_GABAA_IN_PY = StateMonitor(PY, 'I_GABAA_IN_PY', record=True)
MI_GABAB_IN_PY = StateMonitor(PY, 'I_GABAB_IN_PY', record=True)
Mg_AMPA_PY_PY = StateMonitor(PY, 'g_AMPA_PY_PY', record=True)
Mg_GABAA_IN_PY = StateMonitor(PY, 'g_GABAA_IN_PY', record=True)
Mg_GABAB_IN_PY = StateMonitor(PY, 'g_GABAB_IN_PY', record=True)


# IN hücre Akım veya kinetik kayıtları
M_IN = SpikeMonitor(IN)
MV_IN = StateMonitor(IN, 'v', record=True)
MI_IN=StateMonitor(IN,'I',record=True)
MI_AMPA_PY_IN = StateMonitor(IN, 'I_AMPA_PY_IN', record=True)
MI_GABAA_IN_IN = StateMonitor(IN, 'I_GABAA_IN_IN', record=True)
MI_GABAB_IN_IN = StateMonitor(IN, 'I_GABAB_IN_IN', record=True)
Mg_AMPA_PY_IN = StateMonitor(IN, 'g_AMPA_PY_IN', record=True)
Mg_GABAA_IN_IN = StateMonitor(IN, 'g_GABAA_IN_IN', record=True)
Mg_GABAB_IN_IN = StateMonitor(IN, 'g_GABAB_IN_IN', record=True)

MCa_i_PY = StateMonitor(PY, 'Ca_i', record=True)

MECa_PY = StateMonitor(PY, 'ECa_PY', record=True)

print "Simulasyon başlatılıyor..."

######################## Simulasyon zamanı ##########################

#run(duration,report='text')

I_ext = 2*5.0 * nA

spasif = 100*ms
saktif = 100*ms

PY.I = 0; run(spasif,report='text')
ssure = spasif

I_m = []

while (ssure<duration):
    print repr(i+1)
    irnd = I_ext*rand(N_PY)
    I_m.append(tile(irnd,(saktif/defaultclock.dt,1)).T)
    PY.I = irnd; run(saktif,report='text')
    irnd = 0*rand(N_PY)
    I_m.append(tile(irnd,(spasif/defaultclock.dt,1)).T)
    PY.I = irnd; run(spasif,report='text')
    ssure = ssure+spasif+saktif
    
irnd = 0*rand(N_PY)
I_m.append(tile(irnd,(spasif/defaultclock.dt,1)).T)
PY.I = irnd; run(spasif,report='text')

#####################################################################

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12}

matplotlib.rc('font', **font)

plt_i1 = 0; plt_i2 = N_PY/2+1; 
# PY ve IN hücre gerilimleri
figure()
ax1=subplot(221); plot(MV_PY.times/ms, MV_PY[plt_i1]/mV)
title('[PY - '+ repr(plt_i1+1) + '] - ' + repr(len(M_PY.spiketimes[plt_i1])));setp( ax1.get_xticklabels(), visible=False)
ax2=subplot(222); plot(MV_PY.times/ms, MV_PY[plt_i2]/mV)
title('[PY - '+ repr(plt_i2+1) + '] - ' + repr(len(M_PY.spiketimes[plt_i2])));setp( ax2.get_xticklabels(), visible=False)
ax3=subplot(223); plot(MV_IN.times/ms, MV_IN[plt_i1]/mV)
title('[IN - '+ repr(plt_i1+1) + '] - ' + repr(len(M_IN.spiketimes[plt_i1])))
ax4=subplot(224); plot(MV_IN.times/ms, MV_IN[plt_i2]/mV)
title('[IN - '+ repr(plt_i2+1) + '] - ' + repr(len(M_IN.spiketimes[plt_i2])))

# harici uyartı akımları
figure()
ax1=subplot(221); plot(MI_PY.times/ms, MI_PY[plt_i1]/nA,'r');setp( ax1.get_xticklabels(), visible=False)
title('[PY - '+ repr(plt_i1+1) + '] - harici akim (nA)')
ax2=subplot(222); plot(MI_PY.times/ms, MI_PY[plt_i2]/nA,'r');setp( ax2.get_xticklabels(), visible=False)
title('[PY - '+ repr(plt_i2+1) + '] - harici akim (nA)')
ax3=subplot(223); plot(MI_IN.times/ms, MI_IN[plt_i1]/nA,'r')
title('[IN - '+ repr(plt_i1+1) + '] - harici akim (nA)')
ax4=subplot(224); plot(MI_IN.times/ms, MI_IN[plt_i2]/nA,'r')
title('[IN - '+ repr(plt_i2+1) + '] - harici akim (nA)')

# Sinaptik akımlar
# PY'e ait sinaptik akımlar
figure()
ax1=subplot(321); plot(MI_AMPA_PY_PY.times/ms, MI_AMPA_PY_PY[plt_i1]/nA,'r'); 
title('PY - '+ repr(plt_i1+1) + ' - AMPA sinaptik akim (nA)');setp(ax1.get_xticklabels(), visible=False)
ax2=subplot(322); plot(MI_AMPA_PY_PY.times/ms, MI_AMPA_PY_PY[plt_i2]/nA,'r'); 
title('PY - '+ repr(plt_i2+1) + ' - AMPA sinaptik akim (nA)');setp(ax2.get_xticklabels(), visible=False)
ax3=subplot(323); plot(MI_GABAA_IN_PY.times/ms, MI_GABAA_IN_PY[plt_i1]/nA,'r')
title('PY - '+ repr(plt_i1+1) + ' - GABAA sinaptik akim (nA)');setp(ax3.get_xticklabels(), visible=False)
ax4=subplot(324); plot(MI_GABAA_IN_PY.times/ms, MI_GABAA_IN_PY[plt_i2]/nA,'r')
title('PY - '+ repr(plt_i2+1) + ' - GABAA sinaptik akim (nA)');setp(ax4.get_xticklabels(), visible=False)
ax5=subplot(325); plot(MI_GABAB_IN_PY.times/ms, MI_GABAB_IN_PY[plt_i1]/nA,'r')
title('PY - '+ repr(plt_i1+1) + ' - GABAB sinaptik akim (nA)')
ax6=subplot(326); plot(MI_GABAB_IN_PY.times/ms, MI_GABAB_IN_PY[plt_i2]/nA,'r')
title('PY - '+ repr(plt_i2+1) + ' - GABAB sinaptik akim (nA)')

# IN'e ait sinaptik akımlar
figure()
ax1=subplot(321); plot(MI_AMPA_PY_IN.times/ms, MI_AMPA_PY_IN[plt_i1]/nA,'r'); 
title('IN - '+ repr(plt_i1+1) + ' - AMPA sinaptik akim (nA)');setp(ax1.get_xticklabels(), visible=False)
ax2=subplot(322); plot(MI_AMPA_PY_IN.times/ms, MI_AMPA_PY_IN[plt_i2]/nA,'r'); 
title('IN - '+ repr(plt_i2+1) + ' - AMPA sinaptik akim (nA)');setp(ax2.get_xticklabels(), visible=False)
ax3=subplot(323); plot(MI_GABAA_IN_IN.times/ms, MI_GABAA_IN_IN[plt_i1]/nA,'r')
title('IN - '+ repr(plt_i1+1) + ' - GABAA sinaptik akim (nA)');setp(ax3.get_xticklabels(), visible=False)
ax4=subplot(324); plot(MI_GABAA_IN_IN.times/ms, MI_GABAA_IN_IN[plt_i2]/nA,'r')
title('IN - '+ repr(plt_i2+1) + ' - GABAA sinaptik akim (nA)');setp(ax4.get_xticklabels(), visible=False)
ax5=subplot(325); plot(MI_GABAB_IN_IN.times/ms, MI_GABAB_IN_IN[plt_i1]/nA,'r')
title('IN - '+ repr(plt_i1+1) + ' - GABAB sinaptik akim (nA)')
ax6=subplot(326); plot(MI_GABAB_IN_IN.times/ms, MI_GABAB_IN_IN[plt_i2]/nA,'r')
title('IN - '+ repr(plt_i2+1) + ' - GABAB sinaptik akim (nA)')

# Sinaptik iletkenlikler
# PY'e ait sinaptik iletkenlikler
figure()
ax1=subplot(321); plot(Mg_AMPA_PY_PY.times/ms, Mg_AMPA_PY_PY[plt_i1]/nS,'r'); 
title('PY - '+ repr(plt_i1+1) + ' - AMPA sinaptik iletkenlik (nS)');setp(ax1.get_xticklabels(), visible=False)
ax2=subplot(322); plot(Mg_AMPA_PY_PY.times/ms, Mg_AMPA_PY_PY[plt_i2]/nS,'r'); 
title('PY - '+ repr(plt_i2+1) + ' - AMPA sinaptik iletkenlik (nS)');setp(ax2.get_xticklabels(), visible=False)
ax3=subplot(323); plot(Mg_GABAA_IN_PY.times/ms, Mg_GABAA_IN_PY[plt_i1]/nS,'r')
title('PY - '+ repr(plt_i1+1) + ' - GABAA sinaptik iletkenlik (nS)');setp(ax3.get_xticklabels(), visible=False)
ax4=subplot(324); plot(Mg_GABAA_IN_PY.times/ms, Mg_GABAA_IN_PY[plt_i2]/nS,'r')
title('PY - '+ repr(plt_i2+1) + ' - GABAA sinaptik iletkenlik (nS)');setp(ax4.get_xticklabels(), visible=False)
ax5=subplot(325); plot(Mg_GABAB_IN_PY.times/ms, Mg_GABAB_IN_PY[plt_i1]/nS,'r')
title('PY - '+ repr(plt_i1+1) + ' - GABAB sinaptik iletkenlik (nS)')
ax6=subplot(326); plot(Mg_GABAB_IN_PY.times/ms, Mg_GABAB_IN_PY[plt_i2]/nS,'r')
title('PY - '+ repr(plt_i2+1) + ' - GABAB sinaptik iletkenlik (nS)')

# IN'ye ait sinaptik iletkenlikler
figure()
ax1=subplot(321); plot(Mg_AMPA_PY_IN.times/ms, Mg_AMPA_PY_IN[plt_i1]/nS,'r'); 
title('IN - '+ repr(plt_i1+1) + ' - AMPA sinaptik iletkenlik (nS)');setp(ax1.get_xticklabels(), visible=False)
ax2=subplot(322); plot(Mg_AMPA_PY_IN.times/ms, Mg_AMPA_PY_IN[plt_i2]/nS,'r'); 
title('IN - '+ repr(plt_i2+1) + ' - AMPA sinaptik iletkenlik (nS)');setp(ax2.get_xticklabels(), visible=False)
ax3=subplot(323); plot(Mg_GABAA_IN_IN.times/ms, Mg_GABAA_IN_IN[plt_i1]/nS,'r')
title('IN - '+ repr(plt_i1+1) + ' - GABAA sinaptik iletkenlik (nS)');setp(ax3.get_xticklabels(), visible=False)
ax4=subplot(324); plot(Mg_GABAA_IN_IN.times/ms, Mg_GABAA_IN_IN[plt_i2]/nS,'r')
title('IN - '+ repr(plt_i2+1) + ' - GABAA sinaptik iletkenlik (nS)');setp(ax4.get_xticklabels(), visible=False)
ax5=subplot(325); plot(Mg_GABAB_IN_IN.times/ms, Mg_GABAB_IN_IN[plt_i1]/nS,'r')
title('IN - '+ repr(plt_i1+1) + ' - GABAB sinaptik iletkenlik (nS)')
ax6=subplot(326); plot(Mg_GABAB_IN_IN.times/ms, Mg_GABAB_IN_IN[plt_i2]/nS,'r')
title('IN - '+ repr(plt_i2+1) + ' - GABAB sinaptik iletkenlik (nS)')

# Hücre içi kalsiyum konsantrasyon değişimi
# PY'nin Ca değişimi
figure()
ax1=subplot(211);plot(MCa_i_PY.times/ms, MCa_i_PY[plt_i1]/(1e-6 * mole * cm**-3),'r')
title('PY - '+ repr(plt_i1+1) + ' - Ca_i degisimi (mmolar)');setp(ax1.get_xticklabels(), visible=False)
ax2=subplot(212);plot(MCa_i_PY.times/ms, MCa_i_PY[plt_i2]/(1e-6 * mole * cm**-3),'r')
title('PY - '+ repr(plt_i2+1) + ' - Ca_i degisimi (mmolar)')

# Hücre içi kalsiyum konsantrasyon değişimi
# PY'nin ECa değişimi
figure()
ax1=subplot(211);plot(MECa_PY.times/ms, MECa_PY[plt_i1]/mV,'r')
title('PY - '+ repr(plt_i1+1) + ' - ECa degisimi (mV)');setp(ax1.get_xticklabels(), visible=False)
ax2=subplot(212);plot(MECa_PY.times/ms, MECa_PY[plt_i2]/mV,'b')
title('PY - '+ repr(plt_i2+1) + ' - ECa degisimi (mV)')

figure()
raster_plot(M_PY,M_IN)

fig=figure()
ax1=subplot(211); 
im=imshow(MV_PY[:,:]/mV, aspect='auto',origin='lower',vmin=-100, vmax=50,\
        extent=(0, duration/ms, 0, N_PY-1))
ylabel('PY - Noron ID'); title('Potansiyel (mV)')
setp(ax1.get_xticklabels(), visible=False)
subplot(212); 
im=imshow(MV_IN[:,:]/mV, aspect='auto',origin='lower',vmin=-100, vmax=50,\
        extent=(0, duration/ms, 0, N_IN-1))
ylabel('IN - Noron ID'); xlabel('Zaman (ms)')
cax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
fig.colorbar(im,cax=cax)

fig=figure()
ax1=subplot(211); 
im=imshow(MV_PY[:,:]/mV, aspect='auto',origin='lower',vmin=-80, vmax=50,\
        extent=(0, duration/ms, 0, N_PY-1))
ylabel('PY - Noron ID'); title('Potansiyel (mV)')
xlim((50*ms)/ms,(1000*ms)/ms);setp(ax1.get_xticklabels(), visible=False)
subplot(212); 
im=imshow(MV_IN[:,:]/mV, aspect='auto',origin='lower',vmin=-80, vmax=50,\
        extent=(0, duration/ms, 0, N_IN-1))
ylabel('IN - Noron ID'); xlabel('Zaman (ms)')
xlim((50*ms)/ms,(1000*ms)/ms)
cax = fig.add_axes([0.91, 0.1, 0.03, 0.8])
fig.colorbar(im,cax=cax)

figure()
ax1=subplot(221)
imshow(AMPA_PY_PY.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('PY');ylabel('PY');title('AMPA PY-PY');setp(ax1.get_xticklabels(), visible=False)
ax2=subplot(222)
imshow(AMPA_PY_IN.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('PY');ylabel('IN');title('AMPA PY-IN');setp(ax2.get_xticklabels(), visible=False)
setp(ax2.get_yticklabels(), visible=False)
ax3=subplot(223)
imshow(GABAA_IN_PY.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('IN');ylabel('PY');title('GABAA IN-PY')
ax4=subplot(224)
imshow(GABAB_IN_PY.w.to_matrix()/nS, aspect='auto',origin='lower',interpolation='nearest')
xlabel('IN');ylabel('PY');title('GABAB IN-PY');setp(ax4.get_yticklabels(), visible=False)

# harici giriş rand olduğunda
a=array(I_m).reshape(N_PY,size(I_m[:])/N_PY)
fig=figure()
im=imshow(a[:,:]/nA, aspect='auto',origin='lower',vmin=0, vmax=I_ext/nA,\
        extent=(0, duration/ms, 0, N_PY-1))
xlabel('Zaman (ms)'); ylabel('PY - Noron ID'); title('Harici akim (nA)')
colorbar(im)

show()
