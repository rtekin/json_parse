# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 00:47:44 2014

@author: rtekin
"""
from brian import *
import scipy.io as sio
import scipy.ndimage as ndimage
from scipy.signal import convolve2d
from scipy.signal import kaiserord, lfilter, firwin #, freqz, find_peaks_cwt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import colorConverter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as clrm
from collections import Iterable
import operator
import random as rnd
from detect_peaks import detect_peaks

#from math import *

is_beetween=lambda x,min,max: 0 if x<min else x if x<max else 0

################# Molar birim tanımları ##############################
Molar = Unit.create(Dimension(mole=1, m=-3), "Molar", "Mol","k")
aMolar = Unit.create_scaled_unit(Molar, "a")
cMolar = Unit.create_scaled_unit(Molar, "c")
ZMolar = Unit.create_scaled_unit(Molar, "Z")
PMolar = Unit.create_scaled_unit(Molar, "P")
dMolar = Unit.create_scaled_unit(Molar, "d")
GMolar = Unit.create_scaled_unit(Molar, "G")
fMolar = Unit.create_scaled_unit(Molar, "f")
hMolar = Unit.create_scaled_unit(Molar, "h")
daMolar = Unit.create_scaled_unit(Molar, "da")
mMolar = Unit.create_scaled_unit(Molar, "m")
nMolar = Unit.create_scaled_unit(Molar, "n")
pMolar = Unit.create_scaled_unit(Molar, "p")
uMolar = Unit.create_scaled_unit(Molar, "u")
TMolar = Unit.create_scaled_unit(Molar, "T")
yMolar = Unit.create_scaled_unit(Molar, "y")
EMolar = Unit.create_scaled_unit(Molar, "E")
zMolar = Unit.create_scaled_unit(Molar, "z")
MMolar = Unit.create_scaled_unit(Molar, "M")
kMolar = Unit.create_scaled_unit(Molar, "k")
YMolar = Unit.create_scaled_unit(Molar, "Y")

##################### Uygulama sonuçlarını yükler ####################
def load_mat_data(fname, itemLst):
    '''
    itemLst de belirtilen uygulama sonuçlarını *.mat dosyadan yükler
    '''
    mat_ws=sio.loadmat(fname,squeeze_me=True)
    AppResults = {};
    for item in itemLst:
        AppResults[item]=mat_ws[item]
    return AppResults

def dict_inames(args,ext):
    '''
    dict nesnesine ait key listesini oluşturur.
    '''
    inames =['MV_%s_times','MV_%s','M_%s_spikes','M_%s_spiketimes', \
            'M_%s_input_spikes','M_%s_input_spiketimes','LFP_%s','LFPx_%s', \
            'M_%s_rateSR2ms','N_%s']
    inames_comm = ['duration','inert_time','dt','f_ext','dispersion_rate','Noise_sigma',\
                    'randCon_rate','kom_rate']
    dict_inames = [];
    for i in range(len(args)):
        dict_inames.extend([s%(args[i]) for s in inames])
    if ext: dict_inames.extend(['N_EXT'])
    dict_inames.extend(inames_comm)
    
    return dict_inames

def release_mem(AppRes,gstr):
    '''
    istenen key değerlerini siler.
    '''
    inames =['MV_%s_times','MV_%s','M_%s_spikes','M_%s_spiketimes', \
            'M_%s_input_spikes','M_%s_input_spiketimes','LFP_%s','LFPx_%s', \
            'M_%s_rateSR2ms']
    
    inames_ext = ['duration','dt','f_ext','dispersion_rate','Noise_sigma',\
                    'randCon_rate','kom_rate']
    
    for i in range(len(inames)):
        AppRes[inames[i]%(gstr)]=array([]);  
    return AppRes

def moving_average_2d(data, window):
    """Moving average on two-dimensional data.
    """
    # Makes sure that the window function is normalized.
    window /= window.sum()
    # Makes sure data array is a numpy array or masked array.
    if type(data).__name__ not in ['ndarray', 'MaskedArray']:
        data = np.asarray(data)

    # The output array has the same dimensions as the input data 
    # (mode='same') and symmetrical boundary conditions are assumed
    # (boundary='symm').
    return convolve2d(data, window, mode='same', boundary='symm')

def smooth(x,window_len=11,window='hanning'):
    if x.ndim != 1:
        raise ValueError, "sadece 1D dizilere uygulanır."

    if x.size < window_len:
        raise ValueError, "giris vektoru pencere boyutundan buyuk olmalı."

    if window_len<3:
        return x
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "pencere 'flat', 'hanning', 'hamming', 'bartlett', 'blackman' biri olmalı..."
    
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

def smooth2(wave,sample_rate,width=100.,cutoff=110.):
    
    #------------------------------------------------
    # FIR filtre
    #------------------------------------------------
    
    #  Nyquist 
    nyq_rate = sample_rate / 2.0
    
    # bandpass frekans 
    width = width/nyq_rate
    
    # stop band titreşim,dB.
    ripple_db = 60.0
    
    # Kaiser FIR filtre.
    N, beta = kaiserord(ripple_db, width)
    
    #  cutoff frekans.
    cutoff_hz = cutoff/nyq_rate
    
    # lowpass Kaiser window  FIR filtre.
    taps = firwin(N, cutoff_hz/nyq_rate, window=('kaiser', beta))
    
    return lfilter(taps, 1.0, wave)

def downsample_1d(myarr,factor,estimator=mean):
    xs = myarr.shape[0]
    crarr = myarr[:xs-(xs % int(factor))]
    dsarr = estimator( np.concatenate([[crarr[i::factor] 
        for i in range(factor)] ]),axis=0)
    return dsarr

############################ Nernst ####################################
# Nernst parametreleri
R=8.31451 * joule/(mole*kelvin)
F=96.485 * kcoulomb/mole # coulomb=amp*second

Nernst = lambda valence, Temperature, InConc, OutConc:\
    R * (Temperature) / (F * valence) * \
    log(OutConc / InConc) 

############################ heaviside ####################################
heaviside1 = lambda x, lrls, Cdur, dt: 1 if (lrls>0 and x<=Cdur+dt) else 0

def heaviside3(x, lrls, Cdur,dt):
    x = np.array(x)
    if x.shape != ():
        y = np.zeros(x.shape)
        y[(np.around(x-lrls,decimals=9) <= (Cdur)) & (lrls>0)] = 1
        #y[((x-lrls-dt) == 0) & (lrls>0)] = 0.5
    else: # special case for 0d array (a number)
        #if ((x-lrls-dt) == 0) & (lrls>0): y = 0.5
        if ((x-lrls) <= Cdur) & (lrls>0): y = 1
        else: y = 0
    return y

###################### salınım tepe frekanslar ###############################
bandlist = {'slow': (0.1, 1), 'delta': (1, 4), 'theta': (4, 6),  \
                'spindle': (6, 15), 'beta': (15, 25), 'gamma': (25, 80),\
                'ripples': (80, 200)} # Hz
                
def peakdetOscBands(psdData, freqs, selbands = None, delta = None):
    """
    spektral veri içinde istenen bandlara ait tepe noktaları tespit eder.
    """
    peaks = {};
    
    if not type(selbands) is list: selbands=['slow','spindle','gamma']
    if not isscalar(delta): delta = freqs[1]-freqs[0]
    for i in range(len(selbands)):
        if selbands[i] in bandlist.keys():
            mn,mx = bandlist['%s'%(selbands[i])]
            D=psdData[np.bitwise_and(freqs>=(mn-delta),freqs<(mx+delta))]
            ind = find(psdData==max(D))[0];
            peaks['%s'%(selbands[i])]=[ind, psdData[ind], freqs[ind]];
    return peaks

def peakdetOscBands2(psdData, freqs, selbands = None, delta = 0, dist = 1.0):
    """
    spektral veri içinde istenen bandlara ait tepe noktaları tespit eder.
    """
    peaks = {};
    
    if not type(selbands) is list: selbands=['slow','spindle','gamma']
    if not isscalar(delta): delta = freqs[1]-freqs[0]
    for i in range(len(selbands)):
        if selbands[i] in bandlist.keys():
            fprec = diff(freqs)[0]
            mn,mx=((array(bandlist['%s'%(selbands[i])]))/fprec+1).astype(int)
            mn,mx=(int(mn-mn*delta),int(mx+mx*delta))
            
            x_temp=psdData[mn-1:mx]
            peak=find(psdData==max(x_temp[detect_peaks(x_temp,valley=False)]))[-1]
            if peak > mx:
                peak_ind=nan; peak_val=nan; peak_freq=nan;
            else:
                peak_ind=peak
                peak_val=psdData[peak]
                peak_freq=freqs[peak]
            
            if ~np.isnan(peak_ind):
                left_margin=peak_ind-int(dist/fprec)
                if left_margin<0 : left_margin=0
                x_temp=psdData[left_margin:peak_ind]
                dpeaks1=detect_peaks(x_temp, valley=True)
                x_temp=psdData[peak_ind:peak_ind+int(dist/fprec)]
                dpeaks2=detect_peaks(x_temp, valley=True)
                if (dpeaks1.size and dpeaks2.size):
                    deep1=find(psdData==min(x_temp[dpeaks1]))[-1]
                    deep2=find(psdData==min(x_temp[dpeaks2]))[-1]
                else:
                    peak_ind=nan; peak_val=nan; peak_freq=nan;
                    deep1=nan; deep2=nan
            else:
                deep1=nan; deep2=nan
            
            peaks['%s'%(selbands[i])]=[peak_ind, peak_val, peak_freq, deep1, deep2];
    return peaks

def check_in_range(val,selbands,fprec):
    for i in range(len(selbands)):
        x=val<((array(bandlist['%s'%(selbands[i])]))/fprec).astype(int)[1]
        if x: 
            break
    return x

def detect_peak_band(freq_val):
    r=None
    for band in bandlist.keys():
        if (freq_val>=bandlist['%s'%(band)][0]) and (freq_val<bandlist['%s'%(band)][1]):
            r=band
    return r            

def peakdetOscBands3(psdData, freqs, max_peak_ind=None, showplt=False):
    """
    spektral veri içinde istenen bandlara ait tepe noktaları tespit eder.
    """
    fprec = diff(freqs)[0]
    
    #if not type(selbands) is list: selbands=['slow','spindle','gamma']
    if max_peak_ind==None: max_peak_ind=int(len(psdData)*fprec)-1
    
    peaks = {};
    for bnd in bandlist.keys():
        peaks['%s'%(bnd)]=[-1, 0, 0];
    
    peak=0; peak_list=[];
    while(True):
        bnd=detect_peak_band(peak*fprec)
        if bnd!=None:
            mn=(bandlist['%s'%(bnd)][-1]+1)/fprec; mx=-1
            peak_list.append(peak)
            peaks['%s'%(bnd)]=[peak, psdData[peak], freqs[peak]];
        else:
            mn=peak+1; mx=-1
            
        x_temp=psdData[mn:mx:]
        peak=find(psdData==(x_temp[detect_peaks(x_temp,valley=False)]).max())[-1]
        if peak>(max_peak_ind/fprec): break
        
    if showplt:
        fig=figure()
        plot(freqs,10*log10(psdData)); xlim(0,max_peak_ind); 
        for a in peaks.keys():
            pp=peaks['%s'%(a)][0]
            if (pp>0):
                annotate('%s:%s'%(a,round(freqs[pp],2)), xy=(freqs[pp],10*log10(psdData[pp])), \
                        xytext = (20, 30),textcoords = 'offset points',\
                        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', color='k'))
                plot(freqs[pp],10*log10(psdData[pp]),'r.',mew=2,ms=10)
        fig.show()
    
    return peaks

#################### iyonik akım kinetik hazırlık ##########################

def IA_init (v):
    m0 = 1.0 / (1+exp(-(v+60)/8.5));
    h0 = 1.0/(1+exp((v+78)/6.));
    return m0,h0

def Ih_init (v,cai,Phi):
     ginc = 2.0; cac = 0.0015; pc = 0.01; k2 = 0.0004; k4 = 0.001; Shift = 0;
     nca = 4; nexp = 1; taum = 20;
     h_inf = 1./(1 + exp((v + 75 - Shift)/5.5));
     tau_s = (taum + 1000. / (exp((v + 71.5 - Shift)/14.2) + 
                          exp(-(v + 89 - Shift)/11.6))) / Phi;
     alpha = h_inf/tau_s;
     beta = (1 - h_inf)/tau_s;
     p10 = 1/(1 + pow((cac/cai),nca));
     o10 = 1/(1 + beta/alpha + pow((p10/pc),nexp));
     o20 = pow((p10/pc),nexp) * o10;
     c10 = (1-o10)
     return p10,c10,o20

def IT_init (v,Phi_m,Phi_h,Cels,cai,cao):
    m0 = 1. / (1+exp(-(v+59)/6.2));
    h0 = 1. / (1+exp((v+83)/4.0));
    ratio = cao/cai;
    eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489)* log(ratio);
    return m0, h0, eca0


def ITs_init (v,Phi_m,Phi_h,Cels,cai,cao):
    Shift = 2;
    m0 = 1/(1 + exp(-(v + Shift + 50)/7.4));
    h0 = 1/(1 + exp((v + Shift + 78)/5.));
    ratio = cao/cai;
    eca0 = 1000*8.31441*(273.15 + Cels)/(2*96489)* log(ratio);
    return m0, h0, eca0

def INaKTK_init (v):
    Vtr = -50.; VtrK = -50.;
    v2 = v - Vtr;
    v2K = v - VtrK;
    Alpha1 = 0.32*(13 - v2)/(exp((13 - v2)/4.) - 1);
    Beta1 = 0.28*(v2 - 40)/(exp((v2 - 40)/5.) - 1);
    m0 = Alpha1/(Alpha1 + Beta1);

    Alpha2 = 0.128*exp((17 - v2)/18.);
    Beta2 = 4./(exp((40 - v2)/5.) + 1);
    h0 = Alpha2/(Alpha2 + Beta2);

    Alpha3 = 0.032*(15 - v2K)/(exp((15 - v2K)/5.) - 1);
    Beta3 = 0.5*exp((10 - v2K)/40.);
    n0 = Alpha3/(Alpha3 + Beta3);
    return m0, h0, n0

def ICaL_init (v,Phi_m,Phi_h):
    Shift = 0;
    vm = v + Shift;
    a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1);
    b = 0.94*exp((-75-vm)/17.);
    m0 = a/(a+b)/Phi_m;
    a = 0.000457*exp((-13-vm)/50.);
    b = 0.0065/(exp((-vm-15)/28.) + 1);
    h0 = a/(a+b)/Phi_h;
    return m0, h0

def IKCa_init (cai):
    Ra = 0.01; Rb = 0.02;
    a = Ra * cai;  # caix = 1
    b = Rb;
    m0 = a/(a+b);
    return m0

def IKm_init (v):
    tha = -30; qa = 9;Ra = 0.001;Rb = 0.001;
    a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
    b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
    m0 = a/(a+b);
    return m0

def IKv_init (v):
    tha = 25; qa = 9; Ra = 0.02; Rb = 0.002; 
    a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
    b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
    m0 = a/(a+b);
    return m0


def trap0 (v, th, a, q) :
    if 1 :#(fabs(v/th) > 1.0e-6):
        return ( a * (v - th) / (1 - exp(-(v - th)/q)) );
    else :
        return (a * q );

def INaCX_init (v):
    Shift = -10; tha = -35; qa = 9; Ra = 0.182; 
    Rb = 0.124; thi1 = -50; thi2 = -75; qi = 5; thinf = -65; qinf = 6.2;
    Rg = 0.0091; Rd = 0.024;
    vm = v + Shift;
    a = trap0(vm,tha,Ra,qa);
    b = trap0(-vm,-tha,Rb,qa);
    m0 = a/(a+b);

    a = trap0(vm,thi1,Rd,qi);
    b = trap0(-vm,-thi2,Rg,qi);
    h0 = 1/(1+exp((vm-thinf)/qinf));
    return m0, h0

def INap_init (v,Phi_m):
    Tet = -42; Sig = 5; f = 0.02;
    tau_m = 0.8/Phi_m;
    m0 = f/(1 + exp(-(v - Tet)/Sig));
    return m0

# 
def calc_lfp(x, dt,rext,rmax,dr,R,sigma1,sigma2,Lambda,epsilon,sigmaR):
    x = np.asarray(x, dtype=float)
    N = x.shape[0]
    if np.log2(N) % 1 > 0:
        x=resize(x,nextpow2(len(x)),0)

    npoints = len(x);          	# nb points should be an exponent of 2
    npt = npoints/2;            # nb of points in frequency 
    df = 1/(npoints*dt);     # freq interval in Hz
    f_max = df*npt;              # max frequency in Hz
    [Zr, Zi, freq] = impedance(f_max,df,rext,rmax,dr,R,sigma1,sigma2,Lambda,epsilon,sigmaR);
    #n = length(x);
    fft_x=realft(x,1); #fft.fft(x);	
    #fft_x2=transform_radix2(x,1)
    #Y = Y([1:npoints/2]);
    Ir, Ii = [], []
    #Ir = list(real(fft_x[0:npt:])/npt); Ii = list(-imag(fft_x[0:npt])/npt);
    #Ir.append(real(fft_x[npt])/npt);Ii[0]=0;Ii.append(0);
    Ir=list(fft_x[0::2]/npt); Ii=list(fft_x[1::2]/npt);
    Ir.insert(npt,fft_x[1]/npt); Ii[0]=0;Ii.insert(npt,0);
    Ir=array(Ir);Ii=array(Ii);
    #Ir = real(fft_x); Ii = imag(fft_x);
    Vext = calc_extracellular(Zr, Zi, Ir, Ii, npt);
    return x, Zr, Zi, freq, fft_x, Ir, Ii, Vext

#----------------------------------------------------------------------------
#  Empedans hesabı

def nextpow2(i):
    n = 2
    while n < i: n = n * 2
    return n

def DFT_slow(x):
    # 1D dizi (x) için ayrık Fourier dönüşümü
    x = np.asarray(x, dtype=float)
    N = x.shape[0]
    n = np.arange(N)
    k = n.reshape((N, 1))
    M = np.exp(-2j * np.pi * k * n / N)
    return np.dot(M, x)

def FFT(x):
    #Rekürsif 1D Cooley-Tukey FFT
    x = np.asarray(x, dtype=float)
    N = x.shape[0]
    if N % 2 > 0:
        x.resize((nextpow2(len(x)),1))
        N = x.shape[0]
        #raise ValueError("size of x must be a power of 2")
    if N <= 32:  # this cutoff should be optimized
        return DFT_slow(x)
    else:
        X_even = FFT(x[::2])
        X_odd = FFT(x[1::2])
        factor = np.exp(-2j * np.pi * np.arange(N) / N)
        return np.concatenate([X_even + factor[:N / 2] * X_odd,
                               X_even + factor[N / 2:] * X_odd])

####################### get spikes in range ###########################
def get_spikes_inrange(spikes,begin=0,end=inf):
    """
    belir bir zaman aralığında oluşmuş spike dizilerini getirir.
    """
    #spks = spikes.copy()
    for i in range(len(spikes)):
        spks=array(spikes[i],ndmin=1)
        if spks.size:
            spikes[i]=spks[np.bitwise_and(spks>begin, spks<=end)];
    return spikes

####################### get data in range ###########################
def get_data_inrange(data,begin=0,end=inf):
    """
    belir bir zaman aralığına ait veriyi dündürür.
    """
    if isinf(end): return data[begin::]
    else: return data[begin:end:]

############################# raster plot #################################
def _take_options(myopts, givenopts):
    """
    varsayılan parametre değerlerini girilen parametre değerleri ile değiştirir.
    """
    for k in myopts.keys():
        if k in givenopts:
            myopts[k] = givenopts.pop(k)

def my_raster_plot(*monitors, **additionalplotoptions):
    myopts = {"title":"", "xlabel":"Time (ms)", 
                  "yticks":"", "ylabel":"Neuron number", "showplot":False,
                  "showgrouplines":False, "spacebetweengroups":0.0,
                  "grouplinecol":"k", 'newfigure':False, 
                  "colorlist":['r.','g.','b.','y.','c.'],
                  'showlast':None}
    _take_options(myopts, additionalplotoptions)
    spacebetween = myopts['spacebetweengroups']
    if myopts['newfigure']:
        figure();
    for i, m in enumerate(monitors):
        mspikes = m.spikes
        if len(mspikes):
            sn, st = array(mspikes).T
        else:
            sn, st = array([]), array([])
        st /= ms
        if len(monitors) > 1:
            sn = i + ((1. - spacebetween) / float(len(m.source))) * sn
        plot(st,sn,myopts['colorlist'][i])
    if myopts['yticks'] != "":
        ytick=[(2*x+1)/2. for x in range(len(myopts['yticks']))]
        yticks(ytick, myopts['yticks'])
    ylabel(myopts['ylabel'])
    xlabel(myopts['xlabel'])
    title(myopts["title"])

def my_raster_plot2(spikes,source,fsize=12, **additionalplotoptions):
    rc("font", size=fsize)    
    myopts = {"title":"", "xlabel":"Time (second)", 
                  "yticks":"", "ylabel":"Neuron number", "showplot":False,
                  "showgrouplines":False, "spacebetweengroups":0.0,
                  "grouplinecol":"k", 'newfigure':False, 
                  "colorlist":['r.','g.','b.','y.','c.'],
                  'showlast':None}
    _take_options(myopts, additionalplotoptions)
    spacebetween = myopts['spacebetweengroups']
    if myopts['newfigure']:
        figure();
    indx = spacebetween; ytick=[]
    for i in range(len(spikes)):
        mspikes = spikes[i]
        sn, st = array(mspikes).T
        st /= second
        if len(spikes) > 1:
            #sn = i + ((1. - spacebetween) / float(source[i])) * sn
            sn = indx + sn
        ytick.append(indx+float(source[i])/2.)
        indx += source[i] + spacebetween
        plot(st,sn,myopts['colorlist'][i])
    if myopts['yticks']!="":
        #ytick=[(2*x+1)/2. for x in range(len(myopts['yticks']))]
        yticks(ytick, myopts['yticks'])
    ylabel(myopts['ylabel'])
    xlabel(myopts['xlabel'])
    title(myopts["title"])

###################### my Plots ###########################

def my_scatter_plot(app, Vals, barw, barcolors, msize, linecolor, lwidth, fsize=12):
    rc("font", size=fsize)
    minVals = np.min(Vals, axis=1); maxVals = np.max(Vals, axis=1); 
    meanVals = np.mean(Vals, axis=1); stdVals = np.std(Vals, axis=1);
    bar(app,meanVals-stdVals,barw,color=barcolors[0])
    bar(app,2*stdVals,barw,color=barcolors[1],bottom=meanVals-stdVals)
    bar(app,maxVals-meanVals-stdVals,barw,color=barcolors[2],bottom=meanVals+stdVals)
    plot(app+barw/2.,meanVals,marker='s',color=linecolor,lw=lwidth, \
    markerfacecolor=[0.55, 0.52, 0.9], \
    markeredgecolor=[0.15, 0.08, 0.03],ms=msize,mew=3)
    xticks(app+barw/2.,[str(x) for x in app]); xlim(min(app)-2*barw,max(app)+2*barw)

def my_adjust_plot2(app,Vals,fsize=16,lwidth=2,msize=8,mewidth=1.5, alpha=1, \
        showgrid=True,glwidth=2,showAnnText=False,AnnTextstep=1,marker='o',leg=None):
    markers      = ['o','s','<','>','d','*','^']
    colors       = ['b','g','k','r','c','m','y']
    markerecolor = ['c','r','g','m','y','b','k']
    markerfcolor = ['r','b','y','k','g','c','m'] # [[rand(),rand(),rand()] for x in range(len(markers))]
    
    if showgrid: grid(color='r', linestyle='--', linewidth=glwidth, alpha=0.5)
    rc("font", size=fsize)
    for i in range(len(Vals)):
        plot(app,Vals[i],linestyle='-',color=colors[i],lw=lwidth,marker=markers[i],ms=msize,\
        mew=mewidth,markerfacecolor=markerfcolor[i],markeredgecolor=markerecolor[i],alpha=alpha,\
        label='%s'%(leg[i]))
    if leg!=None: legend(loc='best',numpoints=1, fancybox=False, \
            framealpha=0.5, prop={'size':fsize})
    tick_params(bottom='off',top='off',left='off',right='off')    
    if showAnnText:
        for i, txt in enumerate(Vals[0::AnnTextstep]):
            annotate(txt, (app[0::AnnTextstep][i],Vals[0::AnnTextstep][i]))
    xmin = array(app).min(); xmax = array(app).max(); xmean = array(app).mean()
    ymin = array(Vals).min(); ymax = array(Vals).max(); ymean = array(Vals).mean(); 
    xlim(xmin-abs(xmean)*0.3,xmax+abs(xmean)*0.3)
    ylim(ymin-abs(ymean)*0.3,ymax+abs(ymean)*0.3)
    
    #xlim(np.min(app)-abs(np.mean(app))*0.2,np.max(app)+abs(np.mean(app))*0.2)
    #ylim(np.min(Vals)-abs(np.mean(Vals))*0.2,np.max(Vals)+abs(np.mean(Vals))*0.2)
    ticklabel_format(style='sci', axis='y', scilimits=(-3,3))

def my_adjust_plot1(app,Vals,fsize=16,lwidth=3,msize=10,mewidth=2, alpha=0.7, \
        showgrid=True,glwidth=2,showAnnText=False,AnnTextstep=1,marker='o', \
        linepoint=None, color='b',xminmax=None,yminmax=None):

    if showgrid: grid(color='r', linestyle='--', linewidth=glwidth, alpha=0.5)
    rc("font", size=fsize)
    plot(app,Vals,linestyle='-',color=color,lw=lwidth,marker=marker,ms=msize,\
        mew=mewidth,markerfacecolor=[0.9, 0.2, 0.2],markeredgecolor=[0.01, 0.01, 0.01],alpha=alpha)
    tick_params(bottom='off',top='off',left='off',right='off')    
    if showAnnText:
        for i, txt in enumerate(Vals[0::AnnTextstep]):
            annotate(txt, (app[0::AnnTextstep][i],Vals[0::AnnTextstep][i]))
    if xminmax==None: xminmax = (min(app), max(app)) #xmean = (xmax-xmin)/2 # mean(app)
    if yminmax==None: yminmax = (min(Vals), max(Vals)) # ymean = (ymax-ymin)/2 # mean(Vals); 
    
    xlim(xminmax[0]-abs(mean(xminmax))*0.1, xminmax[1]+abs(mean(xminmax))*0.1)
    ylim(yminmax[0]-abs(mean(yminmax))*0.1, yminmax[1]+abs(mean(yminmax))*0.1)
        
    if linepoint!=None:
        plt.plot([xminmax[0], xminmax[1]], [linepoint, linepoint], 'r-', lw=4, alpha=0.7)
    #xlim(np.min(app)-abs(np.mean(app))*0.2,np.max(app)+abs(np.mean(app))*0.2)
    #ylim(np.min(Vals)-abs(np.mean(Vals))*0.2,np.max(Vals)+abs(np.mean(Vals))*0.2)
    ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
             
def my_plot_trace3D(times,traces,fig,lines=[],fsize=16,color='b',viewpos=(0,0),gridon=True,axisonoff="on"):
    rc("font", size=fsize)
    ax = Axes3D(fig)
    ax.grid(gridon)
    indexs=arange(0,len(traces),1)
    for i in indexs:
        ax.plot(times,traces[i],zs=i,color=color);
    for line,pos in lines:
        ax.plot(line*ones(len(indexs)),pos*ones(len(indexs)),indexs,zdir='z',color='r',lw=2);
    ax.set_ylim3d(-100,1000)
    ax.view_init(elev=viewpos[0], azim=viewpos[1])
    ax.axis(axisonoff)
    #ndimage.rotate(ax, 90)

cc = lambda arg,alp: colorConverter.to_rgba(arg, alpha=alp)

def my_plot2_3D(times,data,fig,ax=None,yticklabels=None,lines=[],fsize=12,cmap=None,viewpos=None,gridon=True,axisonoff="on",alpha=1,lw=1):
    rc("font", size=fsize)
    if ax==None: ax = Axes3D(fig)
    ax.grid(gridon)
    if viewpos!=None: ax.view_init(elev=viewpos[0], azim=viewpos[1])
    if cmap==None: 
        color = [cc('r',alpha), cc('g',alpha), cc('y',alpha), cc('b',alpha)]
    else:
        color = [cmap(float(x)/len(data)) for x in range(len(data))]
    ax.axis(axisonoff)
    indexs = arange(0,len(data),1)
    for i, arr in enumerate(data):
        xs = times
        ys = arr
        ax.plot(xs, ys, zs=i, zdir='y', color=color, alpha=alpha,lw=lw);
    for line,pos in lines:
        ax.plot(line*ones(len(indexs)),pos*ones(len(indexs)),indexs,zdir='z',color='r',lw=lw);
    if yticklabels!=None:
        ax.set_yticks(indexs)
        ax.set_yticklabels(yticklabels,rotation=-15, verticalalignment='baseline', horizontalalignment='left')
    ax.set_xlim3d(xs.min(), xs.max())
    ax.set_ylim3d(indexs[0], indexs[-1])
    ax.set_zlim3d(ys.min(), ys.max())

def my_poly_3D(xs,data,fig,ax=None,yticklabels=None,lines=[],fsize=12, \
            cmap=None,viewpos=None,gridon=True,axisonoff="on",alpha=1,lw=1):
    rc("font", size=fsize)
    if ax==None: ax = Axes3D(fig)
    ax.grid(gridon)
    if viewpos!=None: ax.view_init(elev=viewpos[0], azim=viewpos[1])
    
    if cmap==None: 
        color = [cc('r',alpha), cc('g',alpha), cc('y',alpha), cc('b',alpha)] 
    else:
        color = [cmap(float(x)/len(data)) for x in range(len(data))]    

    ax.axis(axisonoff)
    indexs = arange(0,len(data),1)
    verts = [];
    for i, arr in enumerate(data):
        xs = xs
        ys = arr
        ys[0], ys[-1] = 0, 0
        #color.append(cmap(float(i)/len(data)))
        verts.append(list(zip(xs, ys)))
        
        # ax.plot(xs, ys, zs=i, zdir='y', color=color, alpha=alpha,lw=lw);
    poly = PolyCollection(verts, color = color)
    poly.set_alpha(alpha)
    ax.add_collection3d(poly, zs=indexs, zdir='y')
    
    for line,pos in lines:
        ax.plot(line*ones(len(indexs)),pos*ones(len(indexs)),indexs,zdir='z',color='r',lw=lw);
    if yticklabels!=None:
        ax.set_yticks(indexs)
        ax.set_yticklabels(yticklabels,rotation=-15, verticalalignment='baseline', horizontalalignment='left')
    
    ax.set_xlim3d(xs.min(), xs.max())
    ax.set_ylim3d(indexs[0], indexs[-1])
    ax.set_zlim3d(ys.min(), ys.max())

def my_plot_3D(xs,data,fig,ax=None,yticklabels=None,lines=[],fsize=12, \
        cmap=None,viewpos=None,gridon=True,axisonoff="on",alpha=1,lw=1):
    rc("font", size=fsize)
    if ax==None: ax = Axes3D(fig)
    ax.grid(gridon)
    if viewpos!=None: ax.view_init(elev=viewpos[0], azim=viewpos[1])
    if cmap==None: 
        color = [cc('r',alpha), cc('g',alpha), cc('y',alpha), cc('b',alpha)]
    else:
        color = [cmap(float(x)/len(data)) for x in range(len(data))]
    ax.axis(axisonoff)
    indexs = arange(0,len(data),1)
    verts = []; 
    for i, arr in enumerate(data):
        xs = xs
        ys = arr
        #color.append(cmap(float(i)/len(data)))
        verts.append(list(zip(xs, ys)))
        
        # ax.plot(xs, ys, zs=i, zdir='y', color=color, alpha=alpha,lw=lw);
    line = LineCollection(verts, color = color, lw=lw)
    line.set_alpha(alpha)
    ax.add_collection3d(line, zs=indexs, zdir='y')
    
    for line,pos in lines:
        ax.plot(line*ones(len(indexs)),pos*ones(len(indexs)),indexs,zdir='z',color='r',lw=lw);
    if yticklabels!=None:
        ax.set_yticks(indexs)
        ax.set_yticklabels(yticklabels,rotation=-15, verticalalignment='baseline', horizontalalignment='left')
    
    ax.set_xlim3d(xs.min(), xs.max())
    ax.set_ylim3d(indexs[0], indexs[-1])
    ax.set_zlim3d(ys.min(), ys.max())

def my_hist_3D(data,fig,bins=10,ax=None,yticklabels=None,dy=0.01,lines=[],fsize=12, \
                cmap=None,viewpos=None,gridon=True,axisonoff="on",alpha=1):
    rc("font", size=fsize)
    if ax==None: ax = Axes3D(fig)
    ax.grid(gridon)
    if viewpos!=None: ax.view_init(elev=viewpos[0], azim=viewpos[1])
    if cmap==None: 
        color = [cc('r',alpha), cc('g',alpha), cc('y',alpha), cc('b',alpha)] 
    else:
        color = [cmap(float(x)/len(data)) for x in range(len(data))]    
    ax.axis(axisonoff)
    indexs = arange(0,len(data),1)
    for i, arr in enumerate(data):
        hist, bin_edges = np.histogram(arr, bins = bins)
        x = bin_edges[:-1]
        y = i*np.ones_like(hist)
        z = np.zeros_like(hist)
        dx = np.diff(bin_edges)
        dy = dy
        dz = hist
        #color = cmap(float(i)/len(data))
        ax.bar3d(x, y, z, dx, dy, dz, color = color[i], alpha = alpha)
    
    for line,pos in lines:
        ax.plot(line*ones(len(indexs)),pos*ones(len(indexs)),indexs,zdir='z',color='r',lw=1);
    
    if yticklabels!=None:
        ax.set_yticks(indexs)
        ax.set_yticklabels(yticklabels,rotation=-15, verticalalignment='baseline', horizontalalignment='left')
    #ax.set_xlim3d(bin_edges.min(), bin_edges.max())
    #ax.set_ylim3d(indexs[0], indexs[-1])
    #ax.set_zlim3d(ys.min(), ys.max())

def my_pcolor_plot(Vals,fsize=16,cmap=plt.cm.jet,zminmax=None):
    rc("font", size=fsize)
    
    if zminmax==None: zminmax = (min(Vals), max(Vals))
    
    PC=pcolor(Vals, cmap='jet',alpha=1.0, vmin=zminmax[0], vmax=zminmax[1])
    
    levels=np.round(np.linspace(zminmax[0],zminmax[1],5),2)
    
    CB = colorbar(PC)
    CB.set_ticks(levels)
    #CB.set_label('(values)',size=fsize)
    CB.ax.tick_params(labelsize=fsize-2)

def my_adjust_hist(Vals,bin=10,fsize=24,norm=False,fcolor='g',maxHist=None,\
        showgrid=False,ecolor='k',alpha=0.7,rwidth=0.75,lwidth=2, linepoint=None, \
        xminmax=None, yminmax=None):
    
    rc("font", size=fsize)
    if showgrid: grid(color='r', linestyle='--', linewidth=1, alpha=0.5)
    if maxHist is None: maxHist = ceil(max(Vals))
    r, bins, patches=hist(Vals,bin,normed=norm,align='mid',\
        color=fcolor,edgecolor=ecolor,alpha=alpha, rwidth=rwidth,linewidth=lwidth)

    if xminmax==None: xminmax = (min(Vals), max(Vals)) #xmean = (xmax-xmin)/2 # mean(app)
    if yminmax==None: yminmax = (min(r), max(r)) # ymean = (ymax-ymin)/2 # mean(Vals); 
    
    xlim(xminmax[0]-(mean(xminmax))*0.1, xminmax[1]+(mean(xminmax))*0.1)
    ylim(yminmax[0]-(mean(yminmax))*0.1, yminmax[1]+(mean(yminmax))*0.1)
    
    if linepoint!=None:
        plt.plot([linepoint, linepoint], [yminmax[0], yminmax[1]], 'r-', lw=4, alpha=0.7)
    
    #xlim(np.min(bins)-abs(np.mean(bins))*0.2,np.max(bins)+abs(np.mean(bins))*0.2)
    #ylim(0,np.max(r)+abs(np.mean(r))*0.1)
    return r, bins, patches

def my_adjust_contourf(app1,app2,Vals,fsize=16,cmap=plt.cm.jet,alpha=1,drwline=False,norm=None,levels=None):
    rc("font", size=fsize)
    origin = 'lower'
    #origin = 'upper'
    #origin = 'image'
    X, Y = np.meshgrid(app1, app2);
    Z=Vals;
    #mn = int(np.floor(Z.min())); mx = int(np.floor(Z.max())); md = (mn+mx)/2;
    #v = np.linspace(mn-mn*0.5, mx, 3, endpoint=False)    
    
    #V=linspace(Vals.min()+1e-6,Vals.max(),20)
    if levels ==None:
        CS = contourf(X, Y, Z, alpha=alpha, cmap=cmap,origin=origin,norm=norm)
    else:
        CS = contourf(X, Y, Z, levels, alpha=alpha, cmap=cmap,origin=origin,norm=norm,extend='both');
        if (levels.max()<Z.max() and levels.min()>Z.min()): CS.extend='both'
        elif levels.max()<Z.max(): CS.extend='max'
        elif levels.min()>Z.min(): CS.extend='min'
        else : CS.extend='neither'
        CS.cmap.set_under('navy')
        CS.cmap.set_over('darkred')

    #CS.set_clim(vmin=0*mn-0*mn*0.5,vmax=mx)
    if (drwline and Vals.min()!=Vals.max()): 
        contour(X, Y, Z, colors='k',origin=origin)
    tick_params(bottom='off',top='off',left='off',right='off')
    cbar=colorbar(CS)#,format='%.3g');
    #cbar = mpl.colorbar.ColorbarBase(CS,extend='both')
    
    #cbar.set_ticks(linspace(mn,mx,5))#([mn, md, mx])
    #cbar.set_ticklabels(linspace(mn,mx,5))#([mn, md, mx])
    
    #cbar.formatter.set_useOffset(False)
    cbar.formatter.set_scientific(True)
    cbar.formatter.set_powerlimits((-3,3))
    cbar.update_ticks() 

def my_adjust_contourf2(app1,app2,Vals,fsize=16,cmap=plt.cm.jet,alpha=1,drwline=False,norm=None,levels=None):
    rc("font", size=fsize)    
    origin = 'lower'
    #origin = 'upper'
    #origin = 'image'
    X, Y = np.meshgrid(app1, app2);
    Z=Vals;
    mn = int(np.floor(Z.min()))
    mx = int(np.floor(Z.max()))
    md = (mx-mn)/2
    ls=9; la = 5
    #v = np.linspace(mn-mn*0.5, mx, 3, endpoint=False)    
    
    #V=linspace(Vals.min()+1e-6,Vals.max(),20)
    if levels ==None:
        CS = contourf(X, Y, Z, alpha=alpha, cmap=cmap,origin=origin,norm=norm)
    else:
        CS = contourf(X, Y, Z, alpha=alpha, cmap=cmap,origin=origin,norm=norm,extend='both');
        CS.cmap.set_under('navy')
        CS.cmap.set_over('maroon')

        if levels=='both':  CS.extend='both'; CS.levels = arange(mn+md/la,mx-md/la,int((mx-mn)/ls))
        if levels=='over':  CS.extend='max'; CS.levels = arange(mn,mx-md/la,int((mx-mn)/ls))
        if levels=='under':  CS.extend='min'; CS.levels = arange(mn+md/la,mx,int((mx-mn)/ls))

        """
        if (levels.max()<Z.max() and levels.min()>Z.min()): CS.extend='both'
        elif levels.max()<Z.max(): CS.extend='max'
        elif levels.min()>Z.min(): CS.extend='min'
        else: CS.extend='neither'
        """
    #CS.set_clim(vmin=0*mn-0*mn*0.5,vmax=mx)
    if (drwline and Vals.min()!=Vals.max()): 
        contour(X, Y, Z, colors='k',origin=origin)
    tick_params(bottom='off',top='off',left='off',right='off')
    cbar=colorbar(CS)#,format='%.3g');
    #cbar = mpl.colorbar.ColorbarBase(CS,extend='both')
    
    #cbar.set_ticks(linspace(mn,mx,5))#([mn, md, mx])
    #cbar.set_ticklabels(linspace(mn,mx,5))#([mn, md, mx])
    
    cbar.formatter.set_useOffset(False)
    cbar.formatter.set_scientific(True)
    cbar.formatter.set_powerlimits((-3,3))
    cbar.update_ticks() 

def my_adjust_surf3D(app1,app2,Vals,fsize=16,cmap=plt.cm.jet,alpha=1):
    rc("font", size=fsize)
    ax=gca(projection='3d')
    
    X, Y = np.meshgrid(app1, app2);
    Z=Vals;
    surf = ax.plot_surface(X, Y, Z, alpha=alpha, cmap=cmap, rstride=1, cstride=1, \
        linewidth=0, antialiased=False)
    ax.tick_params(bottom='off',top='off',left='off',right='off')
    cbar=colorbar(surf,shrink=0.5,aspect=5);
    cbar.formatter.set_scientific(True) 
    cbar.formatter.set_powerlimits((-3,3)) 
    cbar.update_ticks() 

def my_adjust_polar(X,Y,Theta,ara=1,fsize=16,showgrid=True,color='b',\
        lwidth=1,glwidth=1):
    vlen=np.round(sqrt(array(list(flatten(X))[::ara])**2+array(list(flatten(Y))[::ara])**2))
    [x,y] = pol2cart(array(list(flatten(Theta))[::ara]),vlen); 
    h=compass(x,y,fsize=fsize,showgrid=showgrid,color=color,lwidth=lwidth,\
        glwidth=glwidth);
    return h

def my_adjust_imshow(Vals,extent,fsize=16,cmap=plt.cm.jet,showgrid=False, \
    showcolorbar=True,glwidth=1,interpolation='bilinear',zminmax=None):
    rc("font", size=fsize)
    
    if zminmax==None: zminmax = (Vals.min(), Vals.max())

    if showgrid: grid(color='r', linestyle='--', linewidth=glwidth)
    im=imshow(Vals,cmap=cmap, aspect='auto',origin='lower',vmin=zminmax[0], \
    vmax=zminmax[1],extent=extent,interpolation=interpolation)
    if showcolorbar: CB = colorbar(im)

############################################################################

def plotSpectrum(y,Fs):
    """
    Plots a Single-Sided Amplitude Spectrum of y(t)
    """
    n = len(y) # length of the signal
    k = arange(n)
    T = float(n)/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range
    
    Y = fft(y)/n # fft computing and normalization
    Y = Y[range(n/2)]
 
    plot(frq,abs(Y),'r') # plotting the spectrum
    xlabel('Freq (Hz)')
    ylabel('|Y(freq)|')
    return frq, Y



def impedance(f_max,df,rext,rmax,dr,R,sigma1,sigma2,Lambda,epsilon,sigmaR):
    # calculate impedance for the whole range of frequencies
    sigmatab = [];
    
    # tabulate all values of sigma in a table "sigmatab",
    # which avoids calculating them several times
    # siz = f_max/df + 1;
    # sigmatab = (float *) malloc(sizeof(float) * siz);
    r = arange(rext,rmax,dr)
    sigmatab = sigma2 + (sigma1-sigma2) * exp(-(r-R)/Lambda)
    
    # calculate the impedances for each frequency
    sigR = sigma1;
    epsR = epsilon;
    Zr = []; Zi = [];
    freq = arange(0,f_max+df,df)
    for f in freq:          # loop on frequencies
        w = 2*pi*f;
        w2 = w*w;
        ReZ=0;
        ImZ=0;
        #r = arange(rext,rmax,dr)
        sig = sigmatab
        eps = epsilon;
        den = r*r * (sig*sig + w2 * eps*eps);
        ReZ = sum((sig*sigR + w2*eps*epsR) / den);
        ImZ = sum((sig*epsR - sigR*eps) / den);
        Zr.append(dr/(4*pi*sigmaR) * ReZ);	        # impedance (UNITS: Ohm)
        Zi.append(w * dr/(4*pi*sigmaR) * ImZ);
    return Zr, Zi, freq

#
#  Procedure to calculate FFT of source current
#

def runfft(curr):
    n = len(curr)              # length of the signal
    RFT=fft.fft(curr)/n		# fft computing and normalization
    #RFT = RFT[range(n/2)]
    return RFT

#
#  Procedure to calculate extracellular voltage at a given distance
#  (assumes that Ir, Ii contains the Re/Im values of current; UNITS nA-s)
#  (assumes that Zr, Zi contains the Re/Im values of impedance; UNITS Ohm)
#
def calc_extracellular(Zr, Zi, Ir, Ii, npt):
    # compute the product (Zr,Zi)*(Ir,Ii)
    
    Vr = Zr*Ir - Zi*Ii # store in Vr, Vi (w-freq component of Vext)
    Vi = Zr*Ii + Zi*Ir # (UNITS: nV-s)
    
    Vfreq=list(zeros((1,2*npt)));
    
    Vfreq[0][0:2*npt:2]=Vr[0:-1:]
    Vfreq[0][1:(2*npt+1):2]=Vi[0:-1:];
    Vfreq[0][1]=Vr[npt];Vfreq[0][0]=0;
    
    Vfreq=array(Vfreq).T
    #RFT = Vr + Vi*1j     # store into RFT vector for inverse FFT
    
    #Vext = fft.ifft(Vfreq)             # compute inverse FFT and store into Vext
    Vext = realft(Vfreq,-1);     
    return Vext             # Vext is the extracellular voltage (UNITS nV)


def realft(x, isign):
    """
    -------------------------------------------------------------------------------------------
    USES four1
    Calculates the Fourier transform of a set of n real-valued data points. Replaces this data 
    (which is stored in array data(1:n)) by the positive frequency half of its complex Fourier 
    transform. The real-valued first and last components of the complex transform are returned 
    as elements data(1) and data(2), respectively. n must be a power of 2. This routine 
    also calculates the inverse transform of a complex data array if it is the transform of real 
    data. (Result in this case must be multiplied by 2/n.) 
    
    * Reference:  "Numerical Recipes By W.H. Press, B. P. Flannery,    *
    *              S.A. Teukolsky and W.T. Vetterling, Cambridge       *
    *              University Press, 1986" [BIBLI 08].                 *
    
    -------------------------------------------------------------------------------------------
    """
    x = np.asarray(x, dtype=float)
    N = x.shape[0]
    if np.log2(N) % 1 > 0:
        x=resize(x,nextpow2(len(x)),0)

    data=x.copy();
    
    PI=4.0*math.atan(1.0);
    n=len(data);
    theta=PI/(n/2);                      #Initialize the recurrence
    c1 = 0.5;
    if (isign == 1):
        c2=-0.5;
        four1(data,n/2,1);                 #The forward transform is here
    else:
        c2=0.5;                            #Otherwise set up for an inverse transform
        theta=-theta;
    
    wpr=-2.0*sin(0.5*theta)*sin(0.5*theta);
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    n2p3=n+3;
    for i in range(2,int(n/4)+1):            #Case i=1 done separately below
        i1=2*i-1-1;
        i2=i1+1;
        i3=n2p3-i2-1-1;
        i4=i3+1;
        wrs=wr;
        wis=wi;
        h1r=c1*(data[i1]+data[i3]);         #The two separate transforms are separated out of data
        h1i=c1*(data[i2]-data[i4]);
        
        h2r=-c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wrs*h2r-wis*h2i;       #Here they are recombined to form the true transform
        
                                            #of the original real data
        data[i2]=h1i+wrs*h2i+wis*h2r;
        data[i3]=h1r-wrs*h2r+wis*h2i;
        data[i4]=-h1i+wrs*h2i+wis*h2r;
        wtemp=wr;                           #The recurrence
        wr=wr*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    
    if (isign == 1):
        h1r=data[1-1];
        data[1-1]=h1r+data[2-1];
        data[2-1]=h1r-data[2-1];            #Squeeze the first and last data together to get
    else:                  #them all within the original array
        h1r=data[1-1].copy();
        data[1-1]=c1*(h1r+data[2-1]);
        data[2-1]=c1*(h1r-data[2-1]);
        four1(data,n/2,-1);             #This is the inverse transform for the case isign-1
    return data
    
def four1(data, nn, isign):
    """
    ------------------------------------------------------------------------------------------- 
    Replaces data(1:2*nn) by its discrete Fourier transform, if isign is input as 1; or replaces 
    data(1:2*nn) by nn times its inverse discrete Fourier transform, if isign is input as -1. 
    data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn 
    MUST be an integer power of 2 (this is not checked for here). 
    --------------------------------------------------------------------------------------------
    """
    n=2*nn;
    j=1*0;
    i=1*0;
    while (i<=n):    #This is the bit-reversal section of the routine
        if (j > i):
            tempr=data[j].copy();        #Exchange the two complex numbers
            tempi=data[j+1].copy();
            data[j]=data[i].copy();
            data[j+1]=data[i+1].copy();
            data[i]=tempr;
            data[i+1]=tempi;
        m=nn;
        while ((m >=2) and (j >= m)):
            j=j-m;
            m=m/2;
        j=j+m;
        i = i+2;
    mmax=2;                     #Here begins the Danielson-Lanczos section of the routine
    while(n > mmax):            #Outer loop executed log2 nn times
        istep=2*mmax;
        theta=6.28318530717959/(isign*mmax);   #Initialize for the trigonometric recurrence
        wpr=-2.0*sin(0.5*theta)*sin(0.5*theta);
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        m=1*0;
        while (m<mmax):             #Here are the two nested inner loops
            i=m;
            while (i<n):
                j=i+mmax;                      #This is the Danielson-Lanczos formula:
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i]=data[i]+tempr;
                data[i+1]=data[i+1]+tempi;
                i = i+istep;
            wtemp=wr;  #Trigonometric recurrence
            wr=wr*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
            m = m+2;
        mmax=istep;

def plotSpectrum(y,Fs):
    
    """
    Plots a Single-Sided Amplitude Spectrum of y(t)
    """
    n = len(y) # length of the signal
    k = arange(n)
    T = n/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(n/2)] # one side frequency range
    
    Y = fft(y)/n # fft computing and normalization
    Y = Y[range(n/2)]
     
    plot(frq,abs(Y),'r') # plotting the spectrum
    xlabel('Freq (Hz)')
    ylabel('|Y(freq)|')

def resize(l, newsize, filling=None):
    res = list(l)
    if newsize > len(res):
        res.extend([filling for x in xrange(len(res), newsize)])
    else:
        del res[newsize:]
    return array(res)

def myISI(Spikes,un):
    ISIs = [] 
    #for i in range(N):
    if len(Spikes):
        ISIs.append(un*Spikes[0])
        ISIs.extend(diff(un*Spikes))
        ISIs = [ x for x in ISIs if not math.isnan(x) ]
    return ISIs

def pol2cart(th,r):
    '''
    kutupsal dan kartezyen koordinatlara dönüşüm.
    '''
    x = r*cos(th)
    y = r*sin(th)
    return x, y

def cart2pol(x,y):
    '''
     kartezyen den kutupsal koordinatlara dönüşüm.
    '''
    th = arctan2(y+0,x+0);
    r = hypot(x,y);
    return th, r

def compass(x,y,fsize=12,showgrid=True,color='b',lwidth=1,glwidth=1):
    if showgrid: grid(color='r', linestyle='--', linewidth=glwidth, alpha=0.5)
    rc("font", size=fsize)
    xx = array([0.0, 1.0, 0.80, 1.0, 0.8]).reshape(-1,1)
    yy = array([0.0, 0.0, 0.08, 0.0, -0.08]).reshape(-1,1)
    arrow = xx + 1j*yy
    x = x[:]; y = y[:];
    z = array(x + y*1j).reshape(-1,1);z=z.T
    a = arrow * z;
    [th,r] = cart2pol(a.real,a.imag)
    h=polar(th,r,c=color,lw=lwidth);
    return h
    
def vector_strength2(spikes, period):
    '''
    spike dizisinin vector strength değerini hesaplar.
    '''
    VS=np.nan; theta=np.nan; X=[]; Y=[]; thetas=[];
    
    if spikes.size: #len(np.array(spikes,ndmin=1)):  
        # return abs(mean(exp(array(spikes) * 1j * 2 * pi / period)))
        # N = len(spikes)
        X = (cos(spikes * 2 * pi / period))
        Y = (sin(spikes * 2 * pi / period))
        thetas=mod(spikes,period)*(2*pi/period)
        #thetas=(spikes * 2 * pi) / period;
        theta = mean(thetas)
        
        #f = sqrt(mean(X)**2+mean(Y)**2)
        #VS = abs(mean(exp(spikes * 1j * 2 * pi / period)))
        VS = abs(sum(exp(spikes * 1j * 2 * pi / period)))
    return VS, theta, X, Y, thetas
    

def population_VS_Calc2(spikes, period):
    '''
    populasyonun VS ve ilgili tüm parametrelerini hesaplar.
    refs: Ashida, G., ve Carr, C. E. (2009). Effect of sampling frequency on the measurement of phase-locked action potentials
    Carr, C. E., ve Friedman, M. A. (1999). Evolution of time coding systems
    
    0<=VS<=1; 0 asenkron, 1 senkron
    '''
    
    N=len(spikes)

    if type(spikes) is dict:
        spikesvec=list(flatten_list(list(flatten_dict(spikes))))
    else:
        spikes = [np.array(x,ndmin=1) for x in spikes]
        spikesvec=list(flatten_list(spikes))
    
    if len(spikesvec):
        spikesvec.sort()
        [VS,Th,X,Y,thetas] = vector_strength2(array(spikesvec),period)
    else:
        VS=0
    
    return VS,Th,thetas,X,Y
        
def population_VS_Calc(spikes, period, bins,N=None):
    '''
    populasyonun VS ve ilgili tüm parametrelerini hesaplar.
    refs: Ashida, G., ve Carr, C. E. (2009). Effect of sampling frequency on the measurement of phase-locked action potentials
    Carr, C. E., ve Friedman, M. A. (1999). Evolution of time coding systems
    
    0<=VS<=1; 0 asenkron, 1 senkron
    '''
    if N==None: N=len(spikes)
    
    VSs = zeros(N); VSsTh = [NaN for x in range(N)]; 
    VSsTheta = []; 
    VSsX = [[] for x in range(N)]; VSsY = [[] for x in range(N)]; 
    VSsThetaM =  [[] for x in range(N)];
    
    if not len(spikes):
        return VSs,VSsTh,VSsTheta,VSsX,VSsY,squeeze(VSsThetaM),0
        
    for i in range(N):
        [VS,Th,X,Y,thetas] = vector_strength2(spikes[i],period)
        if isnan(VS): VS=0
        VSs[i]=(VS); VSsTh[i]=(Th); VSsTheta.append(thetas)
        VSsX[i].append(X); VSsY[i].append(Y);
        a=histogram(thetas,bins,normed=0)
        VSsThetaM[i].append(a[0])
    # tüm ağın VS değeri ağın tüm spikeları birleştirilerek hesaplanıyor
    meanVS = population_VS_Calc2(spikes,period)[0]
    #meanVS = mean(VSs) #sum([ x for x in VSs if not math.isnan(x)])/len(spikes)
    
    # --->> sil ! # zayıf spike aktivitesi yanıltıcı olduğundan VS frekans ile carpıldı
    #netmf = network_mean_frequency(spikes,N)[0] 
    
    return VSs,VSsTh,VSsTheta,VSsX,VSsY,squeeze(VSsThetaM),meanVS

def flatten_list(l):
    '''
    list birleştirir.
    '''
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, basestring):
            for sub in flatten_list(el):
                yield sub
        else:
            yield el

def flatten_dict(d):
    '''
    dict birleştirir.
    '''
    for k,v in d.items():
        if isinstance(v, dict):
            for item in flatten_dict(v):
                yield [k]+item
        else:
            yield v

def shift_left(seq, n):
    n = n % len(seq)
    return seq[n:] + seq[:n]

def shift_right(seq, n):
    n = n % len(seq)
    return seq[-n:] + seq[:-n]

def spikegen(N,f,d,time):
    spiketimes = [];
    for i in range(N):
        spiketimes.extend(array([(i,((x+1)*1./f + d*(1./f)*(rnd.randrange(-100,100)/100.))*second) for x in range(int(f*time))]))
    return spiketimes

def spikegen_SWS(N,f,ff,d,time):
    spiketimes = [];
    for i in range(N):
        for x in arange(0,time,1/f):
            spkst = np.random.normal(x,d,ff)
            spiketimes.extend(array([(i,a*second) for a in sort(spkst[spkst>0])]))
    return spiketimes

def get_data_from_dict(dict_data, key):
    """
    kutuphane dizisinden istenilen anahtarın değerini alır.
    """
    res = [];
    for i in range(len(dict_data)):
        if key in dict_data[i]:
            res.append(dict_data[i][key])
    return res

def get_minmax_from_dict(dict_data, key):
    """
    kutuphane dizisinden istenilen anahtarın minmax değerini alır.
    """
    res = [];
    for i in range(len(dict_data)):
        if key in dict_data[i]:
            res.append(dict_data[i][key])
    res=list(flatten(res))
    return min(res),max(res)

def get_minmax(data, minmax=None):
    """
    diziden istenilen anahtarın minmax değerini alır.
    """
    # res=np.squeeze(list(flatten(data)))
    #res=(list(flatten(data)))
        
    if minmax==None:
        return data.min(),data.max()
    mnmx=array(minmax)
    if data.min()<minmax[0]:
        mnmx[0]=data.min()
    if data.max()>minmax[1]:
        mnmx[1]=data.max()
    return mnmx[0],mnmx[1]

def getAppindx(p_Val,sig_Val,selIndx):
    """
    p_Val ve sig_Val değerinin App liste (AppRes_Lst) indeksini bulur
    """
    mCol = list(np.round(p_Val,decimals=6)).index(selIndx[0]);
    mRow = list(np.round(sig_Val,decimals=6)).index(selIndx[1]);
    return len(p_Val)*mRow+mCol

################### popülasyon firing koherens ####################
def pop_firing_coherence(basla, bitir, bin, dt, spikes, selfcorr=False):
    '''
    nöron çiftleri arasında tutarlılık/eşzamanlılık, sıfır zaman gecikmesi ve 
    belirli bir zaman penceresine göre spike ateşleme zamanları arasındaki 
    çapraz-korelasyonu ile ölçülmektedir. 
    refs: Wang, X. J., ve Buzsáki, G. (1996). Gamma oscillation by synaptic inhibition in a hippocampal interneuronal network model
    Sun, X. J., Lei, J. Z., Perc, M., Lu, Q. S., ve Lv, S. J. (2011). Effects of channel noise on firing coherence of small-world Hodgkin-Huxley neuronal networks
    
    0<=CC<=1; 0 asenkron, 1 senkron
    '''
    CC = 0 #np.nan; 
    CC_map=[]; spike_map=[];
    if len(spikes):
        N=len(spikes)
        spike_map = calc_spikes(basla, bitir, bin, dt, spikes)
        CC_map = norm_cross_corr(array(spike_map), selfcorr)
        CC=sum(CC_map)/(N*(N-1)); # /(N*(N-1)/2);
    return CC,CC_map,spike_map

def calc_spikes(basla, bitir, bin, dt, spikes):
    # Ağın spike haritası
    s=[];
    N=len(spikes);
    s=[[] for x in range(N)]
    for i in arange(basla,bitir,bin):
        min=i; max=(i+bin/second);
        for x in arange(N):
            #MM=spiketimes[x];
            #a=1*(spikes[x]>=min); b=1*(spikes[x]<=max);
            s[x].append(1*bool(sum(np.bitwise_and(spikes[x]>=min, spikes[x]<=max))))
            #s[x].append(sum(np.bitwise_and(spikes[x]>=min, spikes[x]<=max)))
            if (s[x][-1]>1):
                raise ValueError('myApp:argChk bin=%s - uygun bir değer seçilmeli (n=%s / t=%s).' % (str(bin),str(x),str(i)))
    return s

def norm_cross_corr(M, selfcorr=False):
    # popülasyon eş zamanlılık (coherence)
    # her nöronun diğer nöronlarla spike ilişkisi (normalize edilmiş)
    N=shape(M)[0];
    K=zeros((N,N));
    for x in arange(N):
        for y in arange(N):
            if ((x!=y) or selfcorr) and (sum(M[x,:])*sum(M[y,:])!=0): # kendiyle korelasyonu hesaplanmasın
                K[x,y]=sum(M[x,:]*M[y,:])/sqrt(sum(M[x,:])*sum(M[y,:]));
    return K
####################################

def Voltaj_dep_senk_calc(MV,basla=None,bitir=None,smoothing=False,winlength=0,window='flat',fs=None):
    '''
    bu senkronizasyon ölçütü, gerilimin zaman bağlı dalgalanma aktivitesini esas almaktadır.
    refs: Brunel, N., ve Hansel, D. (2006). How noise affects the synchronization 
    properties of recurrent networks of inhibitory neurons. Neural Computation, 18(5), 1066-1110
    Guo, D., Wang, Q., ve Perc, M. (2012). Complex synchronous behavior in interneuronal 
    networks with delayed inhibitory and fast electrical synapses    
    
    0<=S<=1; 0 asenkron, 1 senkron
    '''
    MV_temp=MV[:,basla:bitir:];
    if smoothing: 
        if fs==None:
            MVm=array([smooth(MV_temp[x],winlength,window) for x in range(len(MV_temp))])
        else:
            MVm=array([smooth2(MV_temp[x],fs,100,110) for x in range(len(MV_temp))])
        
    #if smoothing: MVm=array([downsample_1d(MV_temp[x],winlength) for x in range(len(MV_temp))])
    
    Vt_pop = mean(MVm,axis=0);
    sigV = mean(Vt_pop**2)-mean(Vt_pop)**2;
    sigVi = (mean(MVm**2,axis=1)-(mean(MVm,axis=1))**2);
    #sigVi = array([mean(MVm[x]**2)-(mean(MVm[x]))**2 for x in range(len(MVm))])
    S = sqrt(sigV/mean(sigVi));
    return S

def CV_and_CC_calc2(spikes, ISI_thr=0):
    '''
    ISI ile ilişkili coefﬁcient of variation - CV 
    ve coefficient of correlation-CC
    refs: Tuckwell, HC (1988) Introduction to Theoretical Neurobiology. Cambridge, UK:Cambridge University Press
    Wang, X. J. (1998). Calcium coding and adaptive temporal computation in cortical pyramidal neurons

    Burst spike'lar için uyarlanmıştır    
    
    CV_all_net=1 düzenli, diğer düzensiz;
    CC_all_net = 1 düzenli, diğer düzensiz 
    '''

    N = len(spikes);

    CV_all_net = 1; CC_all_net=1;
    if type(spikes) is dict:
        spikesvec=array(list(flatten_list(list(flatten_dict(spikes)))))
    else:
        spikes = [np.array(x,ndmin=1) for x in spikes]
        spikesvec=array(list(flatten_list(spikes)))
    
    if len(spikesvec):
        spikesvec.sort()
        ISI_all_net=diff(spikesvec);
        std_isi = std(ISI_all_net); mean_isi = mean(ISI_all_net);
        # ağın tüm spikeları birarada olduğundan CV -> sqrt(N)'e yaklaşır
        # CV_all_net=0 için senkron, değer büyüdükçe asenkron olur
        CV_all_net = ((std_isi / mean_isi)-1)/(sqrt(N));
        if len(ISI_all_net)>=2:
            CC=[];
            for j in range(len(ISI_all_net)-1):
                CC.append((ISI_all_net[j+1]-mean_isi)*(ISI_all_net[j]-mean_isi))
            CC_all_net=(1-abs(mean(CC)/(std_isi**2)))
    
    # --->> sil ! # zayıf spike aktivitesi yanıltıcı olduğundan CV,CC frekans ile carpıldı
    # netmf = network_mean_frequency(spikes,N)[0] 
    
    return CV_all_net, CC_all_net

def CV_and_CC_calc(spikes, ISI_thr=0, N=None):
    '''
    ISI ile ilişkili coefﬁcient of variation - CV 
    ve coefficient of correlation-CC
    refs: Tuckwell, HC (1988) Introduction to Theoretical Neurobiology. Cambridge, UK:Cambridge University Press
    Wang, X. J. (1998). Calcium coding and adaptive temporal computation in cortical pyramidal neurons

    Burst spike'lar için uyarlanmıştır    
    
    CV=0 düzenli, CV>=1 düzensiz;
    -1<=CC<=1; -1, 1 asenkron 0'a yakın senkron 
    '''
    if N==None: N = len(spikes);
    
    CVs=ones(N); CCs=ones(N);
    
    if not len(spikes):
        return CVs,CCs,1,1
        
    for i in range(N):
        if spikes[i].size: #len(spikes[i]): #
            spk=[]; spk.append(0); spk.extend(np.array(spikes[i],ndmin=1))
            #spk = np.array(spikes[i],ndmin=1)
            ISIs = diff(spk) # interspike intervals
            if ISI_thr: ISIs=[x for x in ISIs if x>=ISI_thr]
            std_isi = std(ISIs)+1e-25; mean_isi = mean(ISIs)
            CVs[i]=(std_isi/mean_isi);
            if len(ISIs)>=2:
                CC=[];
                for j in range(len(ISIs)-1):
                    CC.append((ISIs[j+1]-mean_isi)*(ISIs[j]-mean_isi))
                CCs[i]=(mean(abs(array(CC)))/(std_isi**2))
                #CCs[i]=abs(mean(CC)/(std_isi**2))
        #disp(i)
                
    meanCVs = sum([x for x in CVs if not math.isnan(x)])/N
    meanCCs = sum([x for x in CCs if not math.isnan(x)])/N
    return CVs,CCs,meanCVs,meanCCs

def xcorr(x, y=None, maxlags=None, norm='biased'):
    # numpy.correlate Cross-correlation
    res=np.nan; lags=np.empty([])
    if len(x):
        N = len(x)
        if y == None:
            y = x
        assert len(x) == len(y), 'x ve y aynı uzunlukta olmalı'
        assert maxlags <= N, 'maxlags veri uzunluğundan küçük olmalı'
        
        if maxlags == None:
            maxlags = N-1
            lags = arange(0, 2*N-1)
        else:
            assert maxlags < N
            lags = arange(N-maxlags-1, N+maxlags)
                  
        res = np.correlate(x, y, mode='full')
        
        if norm == 'biased':
            Nf = float(N)
            res = res[lags] / float(N)
        elif norm == 'unbiased':
            res = res[lags] / (float(N)-abs(arange(-N+1, N)))[lags]
        elif norm == 'coeff':        
            Nf = float(N)
            rms = rms_flat(x) * rms_flat(y)+1e-20 # sıfır ortalama sorunu
            res = res[lags] / rms / Nf
        else:
            res = res[lags]
    
        lags = arange(-maxlags, maxlags+1)        
    return res, lags

def xcov(x, y=None, maxlags=None, norm='biased'):
    #  Cross-covariance 
    mx = len(x)
    if y == None:
        res, lags = xcorr(x-ones(mx)*mean(x),maxlags=maxlags,norm=norm)
    else:
        my = len(y)
        res, lags = xcorr(x-ones(mx)*mean(x),y-ones(my)*mean(y),maxlags=maxlags,norm=norm)
    
    return res, lags

def stochastic_resonance(Pxx,base1=None,base2=None,maxind=None,lrange=None,rrange=None):
    '''
    Stokastik Rezonans (SR), periyodik bir giriş için sinyal-gürültü oranı (SNR) ile ifade edilir.
    SNR hesaplanırken sinyalin spektral güç dağılımı (PSD) kullanılmıştır.
    refs: Stacey, W. C., ve Durand, D. M. (2000). Stochastic resonance improves signal detection in hippocampal CA1 neurons
    Stacey, W. C., ve Durand, D. M. (2002). Noise and coupling affect signal detection and bursting in a simulated physiological neural network
    
    snr=PSD(w_p)/mean(PSD(w_p±base1 to w_p±base2)) ref: Stacey, W. C., Krieger, A., ve Litt, B. (2011). 
    '''
    if maxind==None: maxind, maxval = max(enumerate(Pxx), key=operator.itemgetter(1))
    #peak_freq = freqs[maxind] # tepe frekans
    if maxind<=0: 
        return 0,[],[]
    
    if (lrange!=None and rrange!=None):
        base_indxs_left = lrange
        base_indxs_right = rrange
        snr = Pxx[maxind]/mean(Pxx[np.concatenate((base_indxs_left,base_indxs_right), axis=1).astype(int)])
    elif (base1!=None and base2!=None):
        base_indxs_left = baseline_indexs(maxind,base1,base2,-1)
        base_indxs_right = baseline_indexs(maxind,base1,base2,1)     
        snr = Pxx[maxind]/mean(Pxx[np.concatenate((base_indxs_left,base_indxs_right), axis=1).astype(int)])
    else:
        base_indxs_left = []; # baseline_indexs(peak_freq,peak_freq,0,-1);
        base_indxs_right = []; # baseline_indexs(peak_freq,peak_freq,0,1);
        snr = Pxx[maxind]/mean(Pxx)
    
    return snr,base_indxs_left,base_indxs_right

def coherence_resonance(Pxx,Freqs,base1=None,base2=None,maxind=None,lrange=None,rrange=None):
    '''
    Koherens Rezonans (CR), sinyalin tepe frekans değerine bağlı sinyal-gürültü oranı (SNR) ile ifade edilir.
    SNR hesaplanırken sinyalin spektral güç dağılımı (PSD) kullanılmıştır.
    refs: Wang, Y., Chik, D. T., ve Wang, Z. D. (2000). Coherence resonance and noise-induced synchronization in globally coupled Hodgkin-Huxley neurons
    Stacey, W. C., Lazarewicz, M. T., ve Litt, B. (2009). Synaptic noise and physiological coupling generate high-frequency oscillations in a hippocampal computational model
    
    snr=PSD(w_p)/mean(PSD(w_p±base1 to w_p±base2)) ref: Stacey, W. C., Krieger, A., ve Litt, B. (2011). 
    '''
    if maxind==None: maxind, maxval = max(enumerate(Pxx), key=operator.itemgetter(1))
    #peak_freq = maxind # tepe frekans
    if maxind==0: 
        return 1
    
    snr,bindxleft,bindxright = stochastic_resonance(Pxx,base1,base2,maxind=maxind,lrange=lrange,rrange=rrange) 
    h=snr;
    
    """
    if (base1!=None and base2!=None):
        half_peak_height_left = mean([Pxx[maxind],mean(Pxx[bindxleft])])
        half_peak_height_right = mean([Pxx[maxind],mean(Pxx[bindxright])])
    """
    # base ile tepe arası fazla veri olmadığından yukarıdaki gibi yarı yüksekliği 
    # geçen ilk sol ve sağ PSD frekanslar arası farkı hesaplayacak frekans değeri yok
    # bu yüzden koherans aşağıdaki gibi her iki yanda ilk minimum PSD değerin olduğu frekansla 
    # tepe frekansın toplamının yarısı alınıp, ardından bu iki yandan elde edilen frekans 
    # değerleri arasındaki fark genişlik olarak alınmıştır.
    
    w_half_peak_left = (Freqs[maxind]+Freqs[baseline_min(Pxx,maxind,-1)])/2.
    w_half_peak_right = (Freqs[baseline_min(Pxx,maxind,1)]+Freqs[maxind])/2.
    
    #### yarı tepe yükseklik frekans genişliği (genişlik küçükse(dar) daha sivri tepe)
    # eğer tepenin en yakın iki yanı toplanacaksa:
    #w_half_peak = 2*(w_half_peak_right if w_half_peak_right<w_half_peak_left else w_half_peak_left)
    # yada farklı olsalarda ikisinin doğrudan toplamı bize genişliği verir.
    w_half_peak = w_half_peak_right - w_half_peak_left;

    factor = h*Freqs[maxind]/w_half_peak; # koherens Rezonans
    return factor


def baseline_min(Pxx, maxind, direction, bound=None):
    # koherans için tepe frekansın sağ ve sol min frekanslarını hesaplar
    if bound==None:
        upperlimit = len(Pxx)-1
        lowerlimit = 0
    else:
        upperlimit = maxind+bound
        lowerlimit = maxind-bound
    
    pos=maxind;
    if direction==1:    # direction=1 sağ min
        while(pos < upperlimit):
            if Pxx[pos]<Pxx[pos+1]:
                break
            pos=pos+1;                
        
    else:               # direction=-1 sol min
        while(pos > lowerlimit):
            if Pxx[pos]<Pxx[pos-1]:
                break
            pos=pos-1;
    return pos

def baseline_indexs(maxind, base1, base2, direction):
    indxs=[];
    if direction==-1:    # direction=-1 sol base indeksler
        if (maxind-base1>=0) and (maxind-base2>=0): 
            indxs.extend([x for x in range(maxind-base1,maxind-base2)]);
        else:
            # bir fark oluşsun diye 2'ye böldüm
            indxs.extend([x for x in range(0,maxind)]); 
    else:               # direction=1 sağ base indeksler
        indxs.extend([x for x in range(maxind+base2,maxind+base1)]);
    return indxs

def network_mean_frequency(spikes,N=None):
    '''
    ağın ortalama spike frekansı. 
    ref: Bogaard, A., Parent, J., Zochowski, M., ve Booth, V. (2009). Interaction of cellular and ...
    Stacey, W. C., Krieger, A., ve Litt, B. (2011). Network recruitment to coherent oscillations in ...
    '''
    if N==None:  N=len(spikes);
    net_mean_freq = zeros(N);
    
    for i in range(len(spikes)):
        if len(np.array(spikes[i],ndmin=1)):
            spk=[]; spk.append(0); spk.extend(np.array(spikes[i],ndmin=1))
            tau_spk = mean(diff(spk));
            net_mean_freq[i]=(1./tau_spk)
    return mean(net_mean_freq), net_mean_freq

def net_mean_phase_coherence(spikes, N=None, selfcoh=False, calc_phase=False):
    '''
    ortalama faz uyumu. 
    ref: Bogaard, A., Parent, J., Zochowski, M., ve Booth, V. (2009). Interaction of cellular and ...
    Stacey, W. C., Krieger, A., ve Litt, B. (2011). Network recruitment to coherent oscillations in ...
    
    0<=mean_phase_coh<=1; faz-kilitlenme, 0 zayıf, 1 güçlü
    '''
    if N==None: N=len(spikes)
    
    net_mean_phase_coh = list([[]*a for a in range(N)]);
    
    for i in range(len(net_mean_phase_coh)):
            net_mean_phase_coh[i]=array(zeros(N))
    
    ###################### geçici silinecek ####################
    if not (calc_phase and len(spikes)):
        return 0, net_mean_phase_coh # mean(net_mean_phase_coh), squeeze(net_mean_phase_coh)
    ############################################################
    for x in range(N):
        for y in range(N):
            if (x!=y or selfcoh):
                if bool(len(np.array(spikes[x],ndmin=1))):
                    phase_diff = [];
                    spk1=[]; spk1.append(0); spk1.extend(np.array(spikes[x],ndmin=1)); spk1=array(spk1)
                    #spk1=np.array(spikes[x],ndmin=1)
                    for i in range(len(spk1)-1):
                        pdiff=0;
                        tx1=round(spk1[i],6);
                        tx2=round(spk1[i+1],6);
                        spk2=np.array(spikes[y],ndmin=1)
                        ty=[a for a,b in enumerate(spk2) if (round(b,6)>tx1 and round(b,6)<=tx2)];
                        if ty!=[]:
                            if ty[-1]!=0:
                                pdiff=(exp(1j*2*pi*(round(spk2[ty[-1]-1],6)-tx1)/(tx2-tx1)))
                        phase_diff.append(pdiff)
                    net_mean_phase_coh[x][y]=(abs(mean(array(phase_diff))))
    return mean(net_mean_phase_coh), squeeze(net_mean_phase_coh)


def net_mean_phase_coherence2(spikes, selfcoh=False):
    '''
    ortalama faz uyumu. 
    ref: Bogaard, A., Parent, J., Zochowski, M., ve Booth, V. (2009). Interaction of cellular and ...
    Stacey, W. C., Krieger, A., ve Litt, B. (2011). Network recruitment to coherent oscillations in ...
    
    0<=mean_phase_coh<=1; faz-kilitlenme, 0 zayıf, 1 güçlü
    '''
    N=len(spikes);
    net_mean_phase_coh = list([[]*a for a in range(N)]);
    for x in range(N):
        if bool(len(np.array(spikes[x],ndmin=1))):
            spk=[]; spk.append(0); spk.extend(np.array(spikes[x],ndmin=1))
            for y in range(N):
                if (x!=y or selfcoh):
                    phase_diff = [];
                    for i in range(len(spk)-1):
                        tx1=round(spk[i],6);
                        tx2=round(spk[i+1],6);
                        ty=[a for a in np.array(spikes[y],ndmin=1) if (tx1<=round(a,6)<=tx2)];
                        if ty!=[]:
                            phase_diff.append(exp(1j*2*pi*(round(ty[0],6)-tx1)/(tx2-tx1)))
                        else:
                            phase_diff.append(0+1j*0)
                    net_mean_phase_coh[x].append(abs(mean(array(phase_diff))))
        else:
            net_mean_phase_coh[x].append(abs(mean(array(0+1j*0))))
    return mean(list(flatten(net_mean_phase_coh))), net_mean_phase_coh

def range_spiketimes(spktime, spiketimes):
    for i in arange(len(spiketimes)-1):
        if (round(spiketimes[i],6)<=spktime<=round(spiketimes[i+1],6)):
            return [spiketimes[i],spiketimes[i+1]]

def net_mean_phase_coherence3(spikes, selfcoh=False):
    '''
    ortalama faz uyumu. 
    ref: Bogaard, A., Parent, J., Zochowski, M., ve Booth, V. (2009). Interaction of cellular and ...
    Stacey, W. C., Krieger, A., ve Litt, B. (2011). Network recruitment to coherent oscillations in ...
    
    0<=mean_phase_coh<=1; faz-kilitlenme, 0 zayıf, 1 güçlü
    '''
    N=len(spikes);
    net_mean_phase_coh = list([[]*a for a in range(N)]);
    for x in range(N):
        for y in range(N):
            phase_diff = [];
            if (x!=y or selfcoh):
                if bool(len(np.array(spikes[x],ndmin=1))):
                    spk=[]; spk.append(0); spk.extend(np.array(spikes[y],ndmin=1))
                    for i in range(len(np.array(spikes[x],ndmin=1))):
                        tx=round(spikes[x][i],6);
                        ty=range_spiketimes(tx,spk)
                        if ty!=None:
                            phase_diff.append(exp(1j*2*pi*(tx-ty[0])/(ty[1]-ty[0])))
                        else:
                            phase_diff.append(0+1j*0)
                else:
                    phase_diff.append(0+1j*0)
                net_mean_phase_coh[x].append(abs(mean(array(phase_diff))))
    return mean(net_mean_phase_coh), net_mean_phase_coh

def Synchronous_bursting(spikes):
    '''
    ağın senkronize spike ölçütü.
    tüm ağın sıralı spike dizisinde ISI kullanılmaktadır.
    ref: Bogaard, A., Parent, J., Zochowski, M., ve Booth, V. (2009). Interaction of cellular and ...
    Stacey, W. C., Krieger, A., ve Litt, B. (2011). Network recruitment to coherent oscillations in ...
    
    0-1 aralığında değer alır. ancak bu aralığın üstünde ve altında değerlerede sahip olabilir.
    B=1 ise güçlü burst senkronizasyon daha küçük değerler ise asenkron aktivite olduğunu gösterir.
    '''
    
    N=len(spikes);
    if type(spikes) is dict:
        spikesvec=list(flatten_list(list(flatten_dict(spikes))))
    else:
        spikes = [np.array(x,ndmin=1) for x in spikes]
        spikesvec=list(flatten_list(spikes))
    if len(spikesvec):
        spikesvec.sort()
        ISI_all_net=diff(spikesvec);
        CV_all_net = std(ISI_all_net) / mean(ISI_all_net);
        sync_burst_factor_B = (CV_all_net-1) / (sqrt(N)-1)
    else:
        sync_burst_factor_B=0
    
    return sync_burst_factor_B

def parametric_distance(Fs, Rs, Bs):
    '''
    iki benzetim uygulaması arasında parametrik uzaklık ölçüsü.
    uygulamalar arası benzerliği F,R ve B ölçütleri açısından değerlendirir.
    ref: Bogaard, A., Parent, J., Zochowski, M., ve Booth, V. (2009). Interaction of cellular and ...
    
    D küçükse 1 ve 2 uygulamaları arasındaki benzerlik büyük, büyükse benzerlik yoktur.   
    '''
    Fs_sq = ((Fs[1]-Fs[0])/(Fs[1]+Fs[0]))**2;
    Rs_sq = ((Rs[1]-Rs[0])/(Rs[1]+Rs[0]))**2;
    Bs_sq = ((Bs[1]-Bs[0])/(Bs[1]+Bs[0]))**2;
    
    return (Fs_sq+Rs_sq+Bs_sq)**0.5

def net_activity(spikes):
    '''
    Toplam spike sayısı ve hiç spike ateşlemeyen hücre sayısını verir.
    Ağın spike aktivitesi ile ilgili genel bilgi verir.
    '''
    spikeCount = 0;
    passiveCount = 0;
    for i in range(len(spikes)):
        spikes[i]=np.array(spikes[i])
        if spikes[i].size:#len(spikes[i]):
            spikeCount+=spikes[i].size#len(spikes[i]);
        else:
            passiveCount+=1;
    return spikeCount,passiveCount

