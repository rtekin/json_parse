# -*- coding: utf-8 -*-
"""
Created on Fri Apr  4 13:48:53 2014

@author: rtekin
"""

from brian import *
import os
import fnmatch
import textwrap
import time
from collections import OrderedDict
from matplotlib.font_manager import fontManager, FontProperties
from tezmodul import *
from detect_peaks import detect_peaks
import shutil

# Kortikal hücre kompartmanlar arası bağlantı rezistansı
kappa_IN = kappa_PY = 10.0e3 * kohm; # 1/kohm = mS
# dendritik-aksosomatik alan oranı
rho_PY = 165; # RS PY için rho=165, FS IN için rho=50 
rho_IN = 50; 

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

# Talamik hücre
# Ra = 100, nseg = 1, L=diam=96 um, area=pi*L*d
area_TC = 29000 * umetre ** 2 # yada 2.9e-4*(cm**2)
# Ra = 100, nseg = 1, L=64.86um, diam=70 um, area=pi*L*d
area_RE = 14300 * umetre ** 2 # yada 1.43e-4*(cm**2)

AppIndx=2;
isLong=True; # False #

AppInf = {}

if AppIndx==0: # Kortikal
    AppInf['AppIndx'] = AppIndx
    AppInf['AppName'] = 'KortikalNet'+isLong*'L'
    AppInf['AppID'] = '8s_t3' # '8s_t4' # 
    AppInf['Groups'] = ['PY','IN']
    AppInf['colorlist']=['k.','b.']
    AppInf['GrpArea'] = [area_DEND_PY,area_DEND_IN]
    #AppInf['trialSel'] = [(0.0,0.0),(0.0,-0.3),(-0.1,-0.3),(0.1,0.2),(0.1,0.5),(0.3,0.1),(0.1,-0.5),(0.4,0.0),(0.3,0.4)]
    #AppInf['trialSel'] = [(0.5,0.02),(0.5,0.05),(0.9,0.06),(0.8,0.09),(0.3,0.05),(0.7,0.04),(0.6,0.05),(0.4,0.03),(0.7,0.1),(0.9,0.06),(0.1,0.04),(0.1,0.01),(0.1,0.08),(0.6,0.02),(0.7,0.07)]  # '8s_t3'
    AppInf['trialSel'] = [(0.9,0.08),(0.9,0.07),(0.9,0.05),(0.9,0.04),(0.9,0.03),(0.7,0.08),(0.7,0.06),(0.7,0.01),(0.6,0.06),(0.4,0.09),(0.4,0.08),(0.4,0.04),(0.1,0.03),(0.7,0.03),(0.7,0.01)]  # '8s_t3'
    AppInf['trialSeld'] = [(0.1,0.03),(0.7,0.03),(0.7,0.01)]
    AppInf['sel_p_indx']=[1,6,7,9]
    AppInf['sel_sig_indx']=[1,3,4,6]
    AppInf['bases']=[25, 10] # [9/10., 5/10.]
    AppInf['SW']=[0.1, 5.]
    AppInf['calc_phase']=True
    AppInf['bands']= ['slow','delta','theta','spindle','beta','gamma','ripples'];
    AppInf['Ext_input']=True
    AppInf['NFFT']= 2**20 # 2**16 # 
    #AppInf['norm']= True
    AppInf['winlength']= (10+5) # 10 # voltaj bağımlı
    AppInf['tbin']= (10+5)*ms # 10*ms # populasyon spike koherens
    AppInf['psdylim1']=(60,-80)
    AppInf['psdylim2']=(60,-60)

elif AppIndx==1: # Talamik
    AppInf['AppIndx'] = AppIndx
    AppInf['AppName'] = 'TalamikNet'+isLong*'L'
    AppInf['AppID'] = '20s_t4' # '20s_t5' # 
    AppInf['Groups'] = ['TC','RE']
    AppInf['colorlist']=['k.','b.']
    AppInf['GrpArea'] = [area_TC,area_RE]
    #AppInf['trialSel'] = [(0.0,0.0),(0.0,-0.3),(-0.1,-0.3),(0.1,0.2),(0.1,0.5),(0.3,0.1),(0.1,-0.5),(0.4,0.0),(0.3,0.4)]
    #AppInf['trialSel'] = [(0.1,0.03),(0.7,0.03),(0.7,0.01),(0.5,0.01),(0.9,0.01),(0.1,0.08),(0.7,0.08),(0.2,0.07),(0.8,0.07),(0.8,0.04),(0.1,0.01),(0.9,0.08),(0.9,0.04),(0.3,0.01),(0.9,0.01)] # '20s_t4'
    AppInf['trialSel'] = [(0.5,0.06),(0.2,0.06),(0.7,0.06),(0.9,0.02),(0.9,0.03),(0.9,0.05),(0.2,0.05),(0.5,0.02),(0.3,0.05),(0.1,0.05),(0.5,0.05),(0.1,0.03),(0.5,0.03),(0.3,0.02),(0.4,0.02)] # '20s_t4'
    AppInf['trialSeld'] = [(0.2,0.05),(0.9,0.05),(0.9,0.03)]    
    AppInf['sel_p_indx']=[1,3,7,9]
    AppInf['sel_sig_indx']=[1,5,8,9]
    AppInf['bases']=[25, 10] # [9/10., 5/10.] # 
    AppInf['SW']=[0.1, 5.]
    AppInf['calc_phase']=True
    AppInf['bands']=['delta','spindle','beta','gamma'];
    AppInf['Ext_input']=False
    AppInf['NFFT']= 2**20 # 2**17 #  
    #AppInf['norm']= True
    AppInf['winlength']= (30) # 30+20 # voltaj bağımlı
    AppInf['tbin']= (5*20)*ms # (15+15)*ms # populasyon spike koherens
    AppInf['psdylim1']=(60,-80)
    AppInf['psdylim2']=(60,-60)
    
else: # Talamokortikal
    AppInf['AppIndx'] = AppIndx
    AppInf['AppName'] = 'TalamoKortikalNet'+isLong*'L'
    AppInf['AppID'] = '20s_t5' # '20s_t4' # 
    AppInf['Groups'] = ['PY','IN','TC','RE']
    AppInf['colorlist']=['k.','b.','r.','g.']
    AppInf['GrpArea'] = [area_DEND_PY,area_DEND_IN,area_TC,area_RE]
    #AppInf['trialSel'] = [(0.0,0.0),(0.4,-0.4),(0.2,0.0),(0.4,0.4),(0.2,0.3),(0.3,-0.1),(0.2,0.0),(0.0,-0.2),(-0.4,0.4)]
    #AppInf['trialSel'] = [(0.1,0.06),(0.5,0.06),(0.7,0.06),(0.5,0.01),(0.5,0.03),(0.6,0.03),(0.6,0.07),(0.6,0.08),(0.1,0.03),(0.3,0.03),(0.5,0.07),(0.8,0.07),(0.8,0.08),(0.8,0.09),(0.8,0.1)] # '20s_t4'
    AppInf['trialSel'] = [(0.1,0.01),(0.1,0.03),(0.9,0.03),(0.9,0.01),(0.1,0.08),(0.1,0.04),(0.8,0.04),(0.8,0.08),(0.2,0.02),(0.2,0.05),(0.9,0.02),(0.9,0.05),(0.4,0.01),(0.4,0.09),(0.9,0.09)] # '20s_t5'
    #AppInf['trialSeld'] = [(0.1,0.06),(0.1,0.03),(0.5,0.06),(0.5,0.03)] # '20s_t4'
    AppInf['trialSeld'] = [(0.1,0.08),(0.1,0.04),(0.8,0.08),(0.8,0.04)] # '20s_t5'
    AppInf['sel_p_indx']=[1,3,7,9]
    AppInf['sel_sig_indx']=[1,3,5,9]
    AppInf['bases']=[25, 10] # [9/10., 5/10.]
    AppInf['SW']=[0.1, 5.]
    AppInf['calc_phase']=True
    AppInf['bands']= ['slow','delta','theta','spindle','beta','gamma','ripples'];
    AppInf['Ext_input']=False
    AppInf['NFFT']= 2**20 # 2**17 #  
    #AppInf['norm']= True
    AppInf['winlength']=(10+20) # voltaj bağımlı
    AppInf['tbin']=(100)*ms # populasyon spike koherens
    AppInf['psdylim1']=(60,-80)
    AppInf['psdylim2']=(60,-60)
    
#main_dir = './deneme_ubuntu/tez_data/%s' %(AppInf['AppName']) 
#main_dir = '/media/rtekin/TOSHIBA EXT/tez_data/%s_20s_t4' %(AppInf['AppName'])
main_dir = '/media/rtekin/TOSHIBA EXT/tez_data/%s_%s' %(AppInf['AppName'],AppInf['AppID'])
AppInf['path_data'] = '%s/data/' %(main_dir)
AppInf['path_images'] = '%s/images/' %(main_dir)
AppInf['AppStr'] = 'SWEandNoE'

os.chdir('/home/rtekin/Desktop/scriptler/tez_scrpts/Bazhenov2002/')
#os.chdir('/media/rtekin/OS/ramazan/scriptler/tez_scrpts/Bazhenov2002/')


if not os.path.exists(AppInf['path_data']):
    raise ValueError('myApp:path=%s - mevcut değil yada herhangi bir kayıt bulunmuyor...' % (AppInf['path_data'] ))

AppRes_Lst = []
itemLst=dict_inames(AppInf['Groups'],AppInf['Ext_input'])

def get_keys(x):
    pos1=x.find("_sig_"); pos2=x.find("_p_"); pos3=x.find("_k_")
    sigx=float((x[pos1+5:pos2]).replace('_','.'))
    px=float((x[pos2+3:pos3]).replace('_','.'))
    return (sigx,px)
    
filtre='*sig*p*.mat'
filelist=sorted([f for f in os.listdir(AppInf['path_data']) if fnmatch.fnmatchcase(f,filtre)],key=get_keys)

if filelist[0].find("nS")>=0: AppInf['AppStr'] = AppInf['AppStr']+'-nS';
#ind=[0,5,10,55,60,65,110,115,120]; filelist = [ filelist[i] for i in ind]
for i in range(len(filelist)):
    
    tic=time.clock();
    fname = AppInf['path_data']+filelist[i]
    AppRes=load_mat_data(fname,itemLst)
    for j in range(len(AppInf['Groups'])):
        #j=1
        gstr = AppInf['Groups'][j];
        N=AppRes['N_%s'%(gstr)];
        dt=double(AppRes['dt'])*second;
        duration = double(AppRes['duration'])*second
        inert_time = double(AppRes['inert_time'])*second
        cond=(round(AppRes['randCon_rate'],6),round(AppRes['Noise_sigma'],6)) in AppInf['trialSel']

        AppRes['M_%s_spiketimes'%(gstr)]=get_spikes_inrange(AppRes['M_%s_spiketimes'%(gstr)],begin=inert_time/second); #
        if not i: AppInf['%s_MnMxspktimes'%(gstr)]=get_minmax(np.array(list(flatten(AppRes['M_%s_spiketimes'%(gstr)])),ndmin=1))
        if cond: AppInf['%s_MnMxspktimes'%(gstr)]=get_minmax(np.array(list(flatten(AppRes['M_%s_spiketimes'%(gstr)])),ndmin=1), AppInf['%s_MnMxspktimes'%(gstr)])
        
        AppRes['M_%s_rateSR2ms'%(gstr)]=get_data_inrange(AppRes['M_%s_rateSR2ms'%(gstr)],begin=(inert_time-dt)/dt); #
        if not i: AppInf['%s_MnMxrateSR2ms'%(gstr)]=get_minmax(AppRes['M_%s_rateSR2ms'%(gstr)])       
        if cond: AppInf['%s_MnMxrateSR2ms'%(gstr)]=get_minmax(AppRes['M_%s_rateSR2ms'%(gstr)], AppInf['%s_MnMxrateSR2ms'%(gstr)])
        
        # Ağ aktivite durumu
        spkCount,passCount = net_activity(AppRes['M_%s_spiketimes'%(gstr)])
        AppRes['%s_spkCount'%(gstr)]=spkCount;
        if not i: AppInf['%s_MnMxspkCount'%(gstr)]=get_minmax(np.array(AppRes['%s_spkCount'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxspkCount'%(gstr)]=get_minmax(np.array(AppRes['%s_spkCount'%(gstr)],ndmin=1), AppInf['%s_MnMxspkCount'%(gstr)])
        
        AppRes['%s_passCount'%(gstr)]=passCount
        if not i: AppInf['%s_MnMxpassCount'%(gstr)]=get_minmax(np.array(AppRes['%s_passCount'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxpassCount'%(gstr)]=get_minmax(np.array(AppRes['%s_passCount'%(gstr)],ndmin=1), AppInf['%s_MnMxpassCount'%(gstr)])
        
        ############################ LFP ############################
        # voltaj ve LFP
        AppRes['%s_total_LFP'%(gstr)]=sum(AppRes['LFP_%s'%(gstr)], axis=0); 
        if not i: AppInf['%s_MnMxtotLFP'%(gstr)]=get_minmax(np.array(detrend(AppRes['%s_total_LFP'%(gstr)]/(mA)),ndmin=1))
        if cond: AppInf['%s_MnMxtotLFP'%(gstr)]=get_minmax(np.array(detrend(AppRes['%s_total_LFP'%(gstr)]/(mA)),ndmin=1), AppInf['%s_MnMxtotLFP'%(gstr)])
        
        AppRes['%s_total_LFPx'%(gstr)]=sum(AppRes['LFPx_%s'%(gstr)], axis=0)
        if not i: AppInf['%s_MnMxtotLFPx'%(gstr)]=get_minmax(np.array(detrend(AppRes['%s_total_LFPx'%(gstr)]/(mA)),ndmin=1))
        if cond: AppInf['%s_MnMxtotLFPx'%(gstr)]=get_minmax(np.array(detrend(AppRes['%s_total_LFPx'%(gstr)]/(mA)),ndmin=1), AppInf['%s_MnMxtotLFPx'%(gstr)])
        
        AppRes['%s_total_membrane'%(gstr)]=sum(AppRes['MV_%s'%(gstr)], axis=0);
        if not i: AppInf['%s_MnMxtotmemb'%(gstr)]=get_minmax(np.array(detrend(AppRes['%s_total_membrane'%(gstr)]/mV),ndmin=1))
        if cond: AppInf['%s_MnMxtotmemb'%(gstr)]=get_minmax(np.array(detrend(AppRes['%s_total_membrane'%(gstr)]/mV),ndmin=1), AppInf['%s_MnMxtotmemb'%(gstr)])
        
        #if not i: AppInf['%s_MnMxmV'%(gstr)]=get_minmax(np.array(AppRes['MV_%s'%(gstr)]/mV,ndmin=1))
        #AppInf['%s_MnMxmV'%(gstr)]=get_minmax(np.array(AppRes['MV_%s'%(gstr)]/mV,ndmin=1), AppInf['%s_MnMxmV'%(gstr)])
        
        # stokastik rezonans (SR) 
        # toplam memb. voltaja göre
        npow1=512;
        npow2=AppInf['NFFT'] # int((duration-inert_time)/dt);
        fprec=float(1*second/dt/AppInf['NFFT'])
        
        PxxTM, freqsTM = psd(detrend(AppRes['%s_total_membrane'%(gstr)]/mV), \
            window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=int(1*second/dt))
        AppRes['%s_PxxTM'%(gstr)]=PxxTM; AppRes['%s_freqsTM'%(gstr)]=freqsTM;
        
        maxindTM, maxvalTM = max(enumerate(PxxTM), key=operator.itemgetter(1))
        #base1=int(maxindTM*AppInf['bases'][0]); base2=int(maxindTM*AppInf['bases'][1]);
        base1=AppInf['bases'][0]; base2=AppInf['bases'][1];
        if base1==base2:base1=maxindTM
        snrTM1,bindxleft,bindexright = stochastic_resonance(PxxTM,base1,base2) 
        snrTM2,bindxleft,bindexright = stochastic_resonance(PxxTM)
        AppRes['%s_snrTM1'%(gstr)]=round(snrTM1,3)
        AppRes['%s_snrTM2'%(gstr)]=round(snrTM2,3)
        AppRes['%s_maxindTM'%(gstr)]=round(freqsTM[maxindTM],3)
        AppRes['%s_maxvalTM'%(gstr)]=round(10*log10(maxvalTM),3) 
        AppRes['%s_mainfreqTM'%(gstr)]=is_beetween(round(freqsTM[maxindTM],3),AppInf['SW'][0],AppInf['SW'][1])

        if not i: AppInf['%s_MnMxsnrTM1'%(gstr)]=get_minmax(np.array(AppRes['%s_snrTM1'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxsnrTM1'%(gstr)]=get_minmax(np.array(AppRes['%s_snrTM1'%(gstr)],ndmin=1), AppInf['%s_MnMxsnrTM1'%(gstr)])
        
        # toplam LFP akıma göre
        PxxTP, freqsTP = psd(detrend(AppRes['%s_total_LFP'%(gstr)]/(mA)), \
            window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt);
        AppRes['%s_PxxTP'%(gstr)]=PxxTP; AppRes['%s_freqsTP'%(gstr)]=freqsTP;
        
        maxindTP, maxvalTP = max(enumerate(PxxTP), key=operator.itemgetter(1))
        #base1=int(maxindTP*AppInf['bases'][0]); base2=int(maxindTP*AppInf['bases'][1]);
        base1=AppInf['bases'][0]; base2=AppInf['bases'][1];
        if base1==base2:base1=maxindTP
        """        
        bandwidth=float(diff(bandlist['slow']))/2
        lrange=range(int(bandlist['slow'][0]/fprec),int((freqsTP[maxindTP]-0.5*bandwidth)/fprec))
        rrange=range(int((freqsTP[maxindTP]+0.5*bandwidth)/fprec),int(bandlist['slow'][1]/fprec))
        """        
        left_margin=maxindTP-int(2*50/fprec); # 2*25
        if left_margin<0 : left_margin=0
        lrange=range(left_margin,maxindTP-int(1/fprec)) # 1
        rrange=range(maxindTP+int(1/fprec),maxindTP+int(2*50/fprec))
        
        snrTP1,bindxleft,bindexright = stochastic_resonance(PxxTP,base1,base2,lrange=lrange,rrange=rrange)
        snrTP2,bindxleft,bindexright = stochastic_resonance(PxxTP)
        AppRes['%s_snrTP1'%(gstr)]=round(snrTP1,3)
        AppRes['%s_snrTP2'%(gstr)]=round(snrTP2,3)
        AppRes['%s_maxindTP'%(gstr)]=round(freqsTP[maxindTP],3)
        AppRes['%s_maxvalTP'%(gstr)]=round(10*log10(maxvalTP),3)
        AppRes['%s_mainfreqTP'%(gstr)]=is_beetween(round(freqsTP[maxindTP],3),AppInf['SW'][0],AppInf['SW'][1])
        
        if not i: AppInf['%s_MnMxsnrTP1'%(gstr)]=get_minmax(np.array(AppRes['%s_snrTP1'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxsnrTP1'%(gstr)]=get_minmax(np.array(AppRes['%s_snrTP1'%(gstr)],ndmin=1), AppInf['%s_MnMxsnrTP1'%(gstr)])
        
        # koherens rezonans (CR) 
        # toplam memb. voltaja göre
        maxindTM, maxvalTM = max(enumerate(PxxTM), key=operator.itemgetter(1))
        #base1=int(maxindTM*AppInf['bases'][0]); base2=int(maxindTM*AppInf['bases'][1]);
        base1=AppInf['bases'][0]; base2=AppInf['bases'][1];
        if base1==base2:base1=maxindTM
        CFactor_betaTM1 = coherence_resonance(PxxTM,freqsTM,base1,base2); # koherens1
        CFactor_betaTM2 = coherence_resonance(PxxTM,freqsTM); # koherens2
        AppRes['%s_CFactor_betaTM1'%(gstr)]=round(CFactor_betaTM1,3)
        AppRes['%s_CFactor_betaTM2'%(gstr)]=round(CFactor_betaTM2,3)
        
        if not i: AppInf['%s_MnMxbetaTM1'%(gstr)]=get_minmax(np.array(AppRes['%s_CFactor_betaTM1'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxbetaTM1'%(gstr)]=get_minmax(np.array(AppRes['%s_CFactor_betaTM1'%(gstr)],ndmin=1), AppInf['%s_MnMxbetaTM1'%(gstr)])
        
        # toplam LFP akıma göre
        maxindTP, maxvalTP = max(enumerate(PxxTP), key=operator.itemgetter(1))
        #base1=int(maxindTP*AppInf['bases'][0]); base2=int(maxindTP*AppInf['bases'][1]);
        base1=AppInf['bases'][0]; base2=AppInf['bases'][1];
        if base1==base2:base1=maxindTP
        #lrange=range(int(bandlist['slow'][0]/fprec),int(freqsTP[maxindTP]/fprec))
        #rrange=range(int((freqsTP[maxindTP]+fprec)/fprec),int(bandlist['slow'][1]/fprec))
        
        CFactor_betaTP1 = coherence_resonance(PxxTP,freqsTP,base1,base2,lrange=lrange,rrange=rrange); # koherens1
        CFactor_betaTP2 = coherence_resonance(PxxTP,freqsTP); # koherens2
        AppRes['%s_CFactor_betaTP1'%(gstr)]=round(CFactor_betaTP1,3)
        AppRes['%s_CFactor_betaTP2'%(gstr)]=round(CFactor_betaTP2,3)
        
        if not i: AppInf['%s_MnMxbetaTP1'%(gstr)]=get_minmax(np.array(AppRes['%s_CFactor_betaTP1'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxbetaTP1'%(gstr)]=get_minmax(np.array(AppRes['%s_CFactor_betaTP1'%(gstr)],ndmin=1), AppInf['%s_MnMxbetaTP1'%(gstr)])
        
        # toplam LFP akıma göre farklı frekans bandları için SNR ve CR
        #peaks=peakdetOscBands(PxxTP,freqsTP,AppInf['bands'],delta = 0)
        #peaks=peakdetOscBands2(PxxTP,freqsTP,AppInf['bands'],delta = 0, dist=1.)
        peaks=peakdetOscBands3(PxxTP,freqsTP,max_peak_ind=200,showplt=False)
        """        
        plot(freqsTP,10*log10(PxxTP)); xlim(0,50)
        for x in range(len(AppInf['bands'])):
            peakind=peaks['%s'%(AppInf['bands'][x])][0];
            if ~np.isnan(peakind):
                plot(freqsTP[peakind],10*log10(PxxTP[peakind]),'r+',ms=10,mew=2)
                annotate('%s'%(round(freqsTP[peakind],2)), xy=(freqsTP[peakind],10*log10(PxxTP[peakind])), xytext = (10, 20),textcoords = 'offset points',\
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', color='k'))
        show();
        """
        for x in range(len(AppInf['bands'])):
            maxindbndTP=peaks['%s'%(AppInf['bands'][x])][0];
            maxvalbndTP=peaks['%s'%(AppInf['bands'][x])][1];
            
            base1=AppInf['bases'][0]; base2=AppInf['bases'][1];
            if base1==base2:base1=maxindbndTP
            
            #lmarginbndTP=peaks['%s'%(AppInf['bands'][x])][3];
            #rmarginbndTP=peaks['%s'%(AppInf['bands'][x])][4];
            
            if maxindbndTP==-1:
                snrTP=0
                CFactor_betaTP=0
            else:
                mn,mx=bandlist['%s'%(AppInf['bands'][x])]
                bandwidth=float(mx-mn)
                
                #left_margin=maxindbndTP-int(2*50/fprec); # 2*50
                #if left_margin<0 : left_margin=0
                
                lrange=range(int(mn/fprec),maxindbndTP-int((2)/fprec)) # 1*2
                rrange=range(maxindbndTP+int((2)/fprec),maxindbndTP+int(2*100/fprec)) 
                
                snrTP,bindxleft,bindexright = stochastic_resonance(PxxTP,base1,base2,maxind=maxindbndTP,lrange=lrange,rrange=rrange)
                CFactor_betaTP = coherence_resonance(PxxTP,freqsTP,base1,base2,maxind=maxindbndTP,lrange=lrange,rrange=rrange); # koherens1
            
            """
            left_margin=maxindbndTP-int(50/fprec); # 2*2
            if left_margin<0 : left_margin=0
            lrange=range(left_margin,maxindbndTP-int(5/fprec)) # 1*2
            rrange=range(maxindbndTP+int(5/fprec),maxindbndTP+int(50/fprec))
            
            #bandwidth=float(diff(bandlist['%s'%(AppInf['bands'][x])]))/2
            #lrange=range(int(bandlist['%s'%(AppInf['bands'][x])][0]/fprec),int((freqsTP[maxindbndTP]-0.5*bandwidth)/fprec))
            #rrange=range(int((freqsTP[maxindbndTP]+0.5*bandwidth)/fprec),int(bandlist['%s'%(AppInf['bands'][x])][1]/fprec))
            
            snrTP,bindxleft,bindexright = stochastic_resonance(PxxTP,base1,base2,maxind=maxindbndTP,lrange=lrange,rrange=rrange)
            CFactor_betaTP = coherence_resonance(PxxTP,freqsTP,base1,base2,maxind=maxindbndTP,lrange=lrange,rrange=rrange); # koherens1
            """
            AppRes['%s_snrTP(%s)'%(gstr,AppInf['bands'][x])]=round(snrTP,3)
            AppRes['%s_CFactor_betaTP(%s)'%(gstr,AppInf['bands'][x])]=round(CFactor_betaTP,3)
            AppRes['%s_maxvalTP(%s)'%(gstr,AppInf['bands'][x])]=round(10*log10(maxvalbndTP),3)
        #################################################################
        
        # gerilim bağımlı senkronizasyon
        basla=0; bitir=int((duration-dt)/dt)
        VSync = Voltaj_dep_senk_calc(AppRes['MV_%s'%(gstr)],basla,bitir, \
                smoothing=True,winlength=AppInf['winlength'],window='hamming') # ,fs=int(1*second/dt)); 
                
        #VSync = Voltaj_dep_senk_calc(AppRes['MV_%s'%(gstr)],basla=basla,bitir=bitir, \
        #        smoothing=True,winlength=AppInf['winlength'],window='hamming',fs=int(1*second/dt));
        AppRes['%s_VSync'%(gstr)]=round(VSync,3)
        
        if not i: AppInf['%s_MnMxVSync'%(gstr)]=get_minmax(np.array(AppRes['%s_VSync'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxVSync'%(gstr)]=get_minmax(np.array(AppRes['%s_VSync'%(gstr)],ndmin=1), AppInf['%s_MnMxVSync'%(gstr)])
        
        # populasyon spike koherens        
        basla=inert_time; bitir=duration; bin=AppInf['tbin']; # 3 ms
        CC,CC_map,spike_map = pop_firing_coherence(basla, bitir, bin, dt, \
                AppRes['M_%s_spiketimes'%(gstr)])
        AppRes['%s_SpikeCoh'%(gstr)]=round(CC,3)
        AppRes['%s_CC_map'%(gstr)]=CC_map
        AppRes['%s_spike_map'%(gstr)]=spike_map
        
        if not i: AppInf['%s_MnMxCCmap'%(gstr)]=get_minmax(np.array(AppRes['%s_CC_map'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxCCmap'%(gstr)]=get_minmax(np.array(AppRes['%s_CC_map'%(gstr)],ndmin=1), AppInf['%s_MnMxCCmap'%(gstr)])
        if not i: AppInf['%s_MnMxSpkCoh'%(gstr)]=get_minmax(np.array(AppRes['%s_SpikeCoh'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxSpkCoh'%(gstr)]=get_minmax(np.array(AppRes['%s_SpikeCoh'%(gstr)],ndmin=1), AppInf['%s_MnMxSpkCoh'%(gstr)])
        
        # coefﬁcient of variation - CV ve coefficient of correlation-CC
        ISI_thr = 0*8*ms/second
        CVs,CCs,mCVs,mCCs = CV_and_CC_calc(AppRes['M_%s_spiketimes'%(gstr)],ISI_thr,N)
        #mCVs,mCCs = CV_and_CC_calc2(AppRes['M_%s_spiketimes'%(gstr)],ISI_thr) # toplu aktiviteye göre uyarladım
        AppRes['%s_meanCV'%(gstr)]=round(mCVs,3); AppRes['%s_CVs'%(gstr)]=CVs
        AppRes['%s_meanCC'%(gstr)]=round(mCCs,3); AppRes['%s_CCs'%(gstr)]=CCs

        # her hücrenin ve ağın ortalama frekansı değeri
        netmf, netfreq = network_mean_frequency(AppRes['M_%s_spiketimes'%(gstr)],N)
        AppRes['%s_meanFreq'%(gstr)]=round(netmf,3); AppRes['%s_netfreq'%(gstr)]=netfreq
        
        if not i: AppInf['%s_MnMxnetFreq'%(gstr)]=get_minmax(np.array(AppRes['%s_netfreq'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxnetFreq'%(gstr)]=get_minmax(np.array(AppRes['%s_netfreq'%(gstr)],ndmin=1), AppInf['%s_MnMxnetFreq'%(gstr)])
        
        if not i: AppInf['%s_MnMxmFreq'%(gstr)]=get_minmax(np.array(AppRes['%s_meanFreq'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxmFreq'%(gstr)]=get_minmax(np.array(AppRes['%s_meanFreq'%(gstr)],ndmin=1), AppInf['%s_MnMxmFreq'%(gstr)])
        
        # ağın ve hücrelerin birbirleriyle ortalama faz uyumu
        netmphase, netphases = net_mean_phase_coherence(AppRes['M_%s_spiketimes'%(gstr)],  \
                selfcoh=True, calc_phase=AppInf['calc_phase'])
        AppRes['%s_meanPhase'%(gstr)]=round(netmphase,3); AppRes['%s_netphases'%(gstr)]=netphases
        
        if not i: AppInf['%s_MnMxnetPhase'%(gstr)]=get_minmax(np.array(AppRes['%s_netphases'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxnetPhase'%(gstr)]=get_minmax(np.array(AppRes['%s_netphases'%(gstr)],ndmin=1), AppInf['%s_MnMxnetPhase'%(gstr)])
        if not i: AppInf['%s_MnMxmPhase'%(gstr)]=get_minmax(np.array(AppRes['%s_meanPhase'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxmPhase'%(gstr)]=get_minmax(np.array(AppRes['%s_meanPhase'%(gstr)],ndmin=1), AppInf['%s_MnMxmPhase'%(gstr)])
        
        # ağın senkronize burst ölçütü
        sync_burst = Synchronous_bursting(AppRes['M_%s_spiketimes'%(gstr)])
        AppRes['%s_SyncBurst'%(gstr)]=round(sync_burst,3)
        
        if not i: AppInf['%s_MnMxSynBrst'%(gstr)]=get_minmax(np.array(AppRes['%s_SyncBurst'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxSynBrst'%(gstr)]=get_minmax(np.array(AppRes['%s_SyncBurst'%(gstr)],ndmin=1), AppInf['%s_MnMxSynBrst'%(gstr)])
        
        # VS parametreleri
        period=1./freqsTP[maxindTP] #1./double(AppRes['f_ext']);
        bins=linspace(0,2*pi,13); # 0:pi/6:2*pi, artım=pi/6=30
        VS,VSTh,VSTheta,VSX,VSY,VSThetaM,mVS=population_VS_Calc(AppRes['M_%s_spiketimes'%(gstr)],period,bins,N)
        AppRes['%s_VS'%(gstr)]=VS; AppRes['%s_meanVS'%(gstr)]=round(mVS,3)
        AppRes['%s_VSTh'%(gstr)]=VSTh; AppRes['%s_VSTheta'%(gstr)]=VSTheta
        AppRes['%s_VSX'%(gstr)]=VSX; AppRes['%s_VSY'%(gstr)]=VSY;
        AppRes['%s_VSThetaM'%(gstr)]=VSThetaM;

        if not i: AppInf['%s_MnMxVS'%(gstr)]=get_minmax(np.array(AppRes['%s_VS'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxVS'%(gstr)]=get_minmax(np.array(AppRes['%s_VS'%(gstr)],ndmin=1), AppInf['%s_MnMxVS'%(gstr)])
        
        if not i: AppInf['%s_MnMxmVS'%(gstr)]=get_minmax(np.array(AppRes['%s_meanVS'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxmVS'%(gstr)]=get_minmax(np.array(AppRes['%s_meanVS'%(gstr)],ndmin=1), AppInf['%s_MnMxmVS'%(gstr)])
        
        if not i: AppInf['%s_MnMxVStheta'%(gstr)]=get_minmax(np.array(AppRes['%s_VSThetaM'%(gstr)],ndmin=1))
        if cond: AppInf['%s_MnMxVStheta'%(gstr)]=get_minmax(np.array(AppRes['%s_VSThetaM'%(gstr)],ndmin=1), AppInf['%s_MnMxVStheta'%(gstr)])
        
        if not i: AppInf['%s_MnMxtotVStheta'%(gstr)]=get_minmax(sum(AppRes['%s_VSThetaM'%(gstr)],axis=0))
        if cond: AppInf['%s_MnMxtotVStheta'%(gstr)]=get_minmax(sum(AppRes['%s_VSThetaM'%(gstr)],axis=0), AppInf['%s_MnMxtotVStheta'%(gstr)])
        
        
        if not cond:
            release_mem(AppRes,gstr)
    AppRes_Lst.append(AppRes)
    toc=time.clock();
    print '['+str(i+1)+']'+fname+' %s sn de işlendi...'%(str(toc-tic))

# Uygulama ayarları
recPic = True # resimler kayıt edilsin mi?
recScript = False

alpha = 0.9; cmap=plt.cm.jet;
barcolors = ['#000080','#f0e68c','#cd5c5c'];
ATShow=False; ATStep = 5; fsize=16; lwidth=3; msize=10; mewidth=2;

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 20}

matplotlib.rc('font', **font)
titfs=30; lblfs=36; ticksfs=30;

dpi1=150; figformat1='.png';

pref_str=AppInf['AppName']+'('+AppInf['AppStr']+')'+'_k_'+('%.0f'%(AppRes_Lst[0]['kom_rate'])).replace('.','_')+ \
    '_dr_'+('%.2f'%(AppRes_Lst[0]['dispersion_rate'])).replace('.','_')+'_t_'+'%.0f'%(int(AppRes_Lst[0]['duration'])*second/ms)

fig_titles=['Spike sayisi','Pasif hucre sayisi','Gerilim bag. senk.', \
    'Spike uyumu ($\kappa$)','VS','CV','CC', \
    'frekans (Hz)','Faz uyumu ($F$)','Burst senkronizasyon ($B$)', \
    'Stokastik Rezonans TM','Koherens rezonans TM', \
    'Network ana frekans (TM) (Hz)',\
    'Max frekans degerleri (TM) (Hz)','Max PSD Degerleri (TM) (dB/Hz)', \
    'Stokastik Rezonans LFP','Koherens rezonans LFP', \
    'Network ana frekans (TP) (Hz)', \
    'Max frekans degerleri (TP)','Max PSD Degerleri (TP) (dB/Hz)'];

if not os.path.exists(AppInf['path_images']):
    os.makedirs(AppInf['path_images'])
    print "Klasör oluşturuldu..."

if recScript:
    shutil.copy(sys.argv[0],AppInf['path_images'])
    print "Script: "+ sys.argv[0] +" Kopyalandı..."

mxm = ceil(sqrt(len(AppRes_Lst)))
myshape=(mxm,mxm)

p = reshape(np.round(get_data_from_dict(AppRes_Lst,'randCon_rate'),1),myshape);
sig = reshape(np.round(get_data_from_dict(AppRes_Lst,'Noise_sigma'),2),myshape);
p_Val = OrderedDict.fromkeys(p.flatten()).keys(); #p_Val.sort();
sig_Val = OrderedDict.fromkeys(sig.flatten()).keys(); #sig_Val.sort()

lengrp=len(AppInf['Groups'])

# uygulamaların aktivite durumları ve ortalama frekans, 
fig=figure(figsize=(30,8*lengrp));
for i in range(lengrp):
    gstr = AppInf['Groups'][i];
    AppInf['%s_spkCount'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_spkCount'%(gstr)),myshape).T;
    AppInf['%s_passCount'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_passCount'%(gstr)),myshape).T;
    AppInf['%s_meanFreq'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_meanFreq'%(gstr)),myshape).T;
    ax=subplot(lengrp,3,3*i+1);
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_spkCount'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[0]).strip(),fontsize=titfs);
    ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
    ax.tick_params(axis='y', labelleft='on')
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs); 
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
        ax.tick_params(axis='x', labelbottom='on')
    ax=subplot(lengrp,3,3*i+2); 
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_passCount'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs);
    ax.tick_params(axis='both', labelbottom='off', labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[1]).strip(),fontsize=titfs);
    #ax.tick_params(axis='y', labelleft='on')
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
        ax.tick_params(axis='x', labelbottom='on')
    ax=subplot(lengrp,3,3*i+3);
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_meanFreq'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs); 
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[7]).strip(),fontsize=titfs);
    #ax.tick_params(axis='y', labelleft='on')
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs); 
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
        ax.tick_params(axis='x', labelbottom='on')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_SpkCnt-PassCnt-Freq)' + figformat1,dpi=dpi1, bbox_inches='tight')

# Voltaj senk. ve Spike uyumu
fig=figure(figsize=(20,8*lengrp));
for i in range(lengrp):
    gstr = AppInf['Groups'][i];
    AppInf['%s_VSync'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_VSync'%(gstr)),myshape).T;
    AppInf['%s_SpikeCoh'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_SpikeCoh'%(gstr)),myshape).T;
    ax=subplot(lengrp,2,2*i+1);
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_VSync'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs); 
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[2]).strip(),fontsize=titfs);
    ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
    ax.tick_params(axis='y', labelleft='on')
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs); 
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
        ax.tick_params(axis='x', labelbottom='on')
    ax=subplot(lengrp,2,2*i+2); 
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_SpikeCoh'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[3]).strip(),fontsize=titfs);
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_VmSync-PopCoh)' + figformat1,dpi=dpi1, bbox_inches='tight')

# ortalama VS, CV ve CC
fig=figure(figsize=(30,8*lengrp));
for i in range(lengrp):
    gstr = AppInf['Groups'][i];
    AppInf['%s_meanVS'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_meanVS'%(gstr)),myshape).T;
    AppInf['%s_meanCV'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_meanCV'%(gstr)),myshape).T;
    AppInf['%s_meanCC'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_meanCC'%(gstr)),myshape).T;
    
    ax=subplot(lengrp,3,3*i+1);
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_meanVS'%(gstr)]/AppInf['%s_meanVS'%(gstr)].max(), \
        cmap=cmap,alpha=alpha,fsize=ticksfs);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[4]).strip(),fontsize=titfs);
    ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
    ax.tick_params(axis='y', labelleft='on')
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs); 
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
    ax=subplot(lengrp,3,3*i+2); 
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_meanCV'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[5]).strip(),fontsize=titfs);
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
    ax=subplot(lengrp,3,3*i+3); 
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_meanCC'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[6]).strip(),fontsize=titfs);
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_VS-CV-CC)' + figformat1,dpi=dpi1, bbox_inches='tight')

# faz uyum ve Burst senk.

levels_phase=[None for a in range(lengrp)];
levels_burst=[None for a in range(lengrp)];
"""
levels_phase=[]; 
levels_phase.append(arange(0.14,0.18,0.005)) #(None) # 
levels_phase.append(arange(0.045,0.063,0.0025)) #(None) # 

levels_burst=[]; 
levels_burst.append(arange(0.05,0.25,0.025));  #(None) # 
levels_burst.append(arange(0.5,3.3,0.35));  #(None) # 
"""
fig=figure(figsize=(20,8*lengrp));
for i in range(lengrp):
    gstr = AppInf['Groups'][i];
    AppInf['%s_meanPhase'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_meanPhase'%(gstr)),myshape).T;
    AppInf['%s_SyncBurst'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_SyncBurst'%(gstr)),myshape).T;
    ax=subplot(lengrp,2,2*i+1); 
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_meanPhase'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs,levels=levels_phase[i]);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[8]).strip(),fontsize=titfs);
    ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
    ax.tick_params(axis='y', labelleft='on')
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')    
    ax=subplot(lengrp,2,2*i+2); 
    if not i: 
        title(textwrap.dedent(fig_titles[9]).strip(),fontsize=titfs);
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_SyncBurst'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs,levels=levels_burst[i]);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_PhsCoh-BrstSync)' + figformat1,dpi=dpi1, bbox_inches='tight')

# Stokastik ve Koherens rezonans
# toplam memb. voltaj'a göre
fig=figure(figsize=(30,8*lengrp));
for i in range(lengrp):
    gstr = AppInf['Groups'][i];
    AppInf['%s_snrTM1'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_snrTM1'%(gstr)),myshape).T;
    AppInf['%s_CFactor_betaTM1'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_CFactor_betaTM1'%(gstr)),myshape).T;
    AppInf['%s_mainfreqTM'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_mainfreqTM'%(gstr)),myshape).T;
    ax=subplot(lengrp,3,3*i+1);
    #my_adjust_contourf(sig_Val,p_Val,AppInf['%s_snrTM1'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs);
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_snrTM1'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs); 
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[10]).strip(),fontsize=titfs);
    ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
    ax.tick_params(axis='y', labelleft='on')
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
    ax=subplot(lengrp,3,3*i+2); 
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_CFactor_betaTM1'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[11]).strip(),fontsize=titfs);
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')    
    ax=subplot(lengrp,3,3*i+3); 
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_mainfreqTM'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[12]).strip(),fontsize=titfs);
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_StocRes1-CohRes1-mainfreq)' + figformat1,dpi=dpi1, bbox_inches='tight')

# tepe PSD ve freq
fig=figure(figsize=(20,8*lengrp));
for i in range(lengrp):
    gstr = AppInf['Groups'][i];
    AppInf['%s_maxindTM'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_maxindTM'%(gstr)),myshape).T;
    AppInf['%s_maxvalTM'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_maxvalTM'%(gstr)),myshape).T;
    ax=subplot(lengrp,2,2*i+1);
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_maxindTM'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs); 
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
    ax.tick_params(axis='y', labelleft='on')
    if not i: 
        title(textwrap.dedent(fig_titles[13]).strip(),fontsize=titfs);
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')    
    ax=subplot(lengrp,2,2*i+2); 
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_maxvalTM'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[14]).strip(),fontsize=titfs);
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_maxindTM-maxvalTM)' + figformat1,dpi=dpi1, bbox_inches='tight')

# toplam LFP ye göre
"""
levels_SR=[None for a in range(lengrp)];
levels_CR=[None for a in range(lengrp)];
levels_mainfreqTP=[None for a in range(lengrp)];
"""
levels_SR=[];
levels_SR.append(arange(0,1501,200)) # 
levels_SR.append(arange(0,1751,250)) # 
levels_SR.append(arange(0,1001,140)) # 
levels_SR.append(arange(0,651,90)) 

levels_CR=[]; 
levels_CR.append(arange(0,3001,400)) # 
levels_CR.append(arange(0,3501,500)) # 
levels_CR.append(arange(0,2101,300)) # 
levels_CR.append(arange(0,1501,200)) # 

levels_mainfreqTP=[]; 
levels_mainfreqTP.append(arange(0.0,1.1,0.1)) # 
levels_mainfreqTP.append(arange(0.0,1.1,0.1)) # 
levels_mainfreqTP.append(arange(0.15,0.51,0.05)) # 
levels_mainfreqTP.append(arange(0.15,0.51,0.05)) # 

fig=figure(figsize=(30,8*lengrp));
for i in range(lengrp):
    gstr = AppInf['Groups'][i];
    AppInf['%s_snrTP1'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_snrTP1'%(gstr)),myshape).T;
    AppInf['%s_CFactor_betaTP1'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_CFactor_betaTP1'%(gstr)),myshape).T;
    AppInf['%s_mainfreqTP'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_mainfreqTP'%(gstr)),myshape).T;
    
    ax=subplot(lengrp,3,3*i+1);
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_snrTP1'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs,levels=levels_SR[i])
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    #my_adjust_contourf(sig_Val,p_Val,AppInf['%s_snrTP1'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs,levels=arange(0,121,15))
    if not i: 
        title(textwrap.dedent(fig_titles[15]).strip(),fontsize=titfs);
    ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
    ax.tick_params(axis='y', labelleft='on')
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
    ax=subplot(lengrp,3,3*i+2); 
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_CFactor_betaTP1'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs,levels=levels_CR[i]);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    #my_adjust_contourf(sig_Val,p_Val,AppInf['%s_CFactor_betaTP1'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs,levels=arange(0,180,20));
    if not i: 
        title(textwrap.dedent(fig_titles[16]).strip(),fontsize=titfs);
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')    
    ax=subplot(lengrp,3,3*i+3); 
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_mainfreqTP'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs,levels=levels_mainfreqTP[i]);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[17]).strip(),fontsize=titfs);
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_StocRes2-CohRes2-mainfreq)' + figformat1,dpi=dpi1, bbox_inches='tight')

# tepe PSD ve freq

levels_mxindTP=[None for a in range(lengrp)];
"""
levels_mxindTP=[]; 
levels_mxindTP.append(arange(0.15,0.51,0.05)) # 
levels_mxindTP.append(arange(0.15,0.51,0.05)) # 
levels_mxindTP.append(arange(0.15,0.51,0.05)) # 
levels_mxindTP.append(arange(0.15,0.51,0.05)) # 
"""
fig=figure(figsize=(20,8*lengrp));
for i in range(lengrp):
    gstr = AppInf['Groups'][i];
    AppInf['%s_maxindTP'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_maxindTP'%(gstr)),myshape).T;
    AppInf['%s_maxvalTP'%(gstr)] = reshape(get_data_from_dict(AppRes_Lst,'%s_maxvalTP'%(gstr)),myshape).T;
    ax=subplot(lengrp,2,2*i+1);
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_maxindTP'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs,levels=levels_mxindTP[i]); 
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[18]).strip(),fontsize=titfs);
    ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
    ax.tick_params(axis='y', labelleft='on')
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')    
    ax=subplot(lengrp,2,2*i+2); 
    my_adjust_contourf(sig_Val,p_Val,AppInf['%s_maxvalTP'%(gstr)],cmap=cmap,alpha=alpha,fsize=ticksfs);
    ax.tick_params(axis='both', labelbottom='off',labelleft='off')
    if not i: 
        title(textwrap.dedent(fig_titles[19]).strip(),fontsize=titfs);
    if (i==(lengrp-1)): 
        xlabel(r'$\sigma_g$',fontsize=lblfs);
        ax.tick_params(axis='x', labelbottom='on')
        ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_maxindTP-maxvalTP)' + figformat1,dpi=dpi1, bbox_inches='tight')

# farklı frekans bandlarının tepe PSD değerleri
# AppInf['bands']=AppInf['bands'][3:]

lenbnd = len(AppInf['bands'])
fig=figure(figsize=(10*lenbnd,8*lengrp));
suptitle(textwrap.dedent(fig_titles[19]).strip(),fontsize=titfs);
for bndindx in range(lenbnd):
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        bnd = AppInf['bands'][bndindx]
        AppInf['%s_maxvalTP(%s)'%(gstr,bnd)] = reshape(get_data_from_dict(AppRes_Lst,'%s_maxvalTP(%s)'%(gstr,bnd)),myshape).T;
        ax=subplot2grid((lengrp,lenbnd),(i,bndindx));
        my_adjust_contourf(sig_Val,p_Val,AppInf['%s_maxvalTP(%s)'%(gstr,bnd)],cmap=cmap,alpha=alpha,fsize=ticksfs); 
        ax.tick_params(axis='both', labelbottom='off',labelleft='off')
        if not i: 
            title('(%s)'%(bnd),fontsize=titfs);
        if (i==(lengrp-1)): 
            xlabel(r'$\sigma_g$',fontsize=lblfs);
            ax.tick_params(axis='x', labelbottom='on')
            ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
        if not bndindx: 
            ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
            ax.tick_params(axis='y', labelleft='on')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_freqbands_maxvalTP)' + figformat1,dpi=dpi1, bbox_inches='tight')

# p ve sig'e bağlı en baskın frekans bandı şeçimi
#bands=AppInf['bands'][0:]
AppInf['bands']=['delta', 'spindle', 'beta', 'gamma']
fig=figure(figsize=(10,8*lengrp));
for i in range(lengrp):
    maxvalTP_PSD=[]
    for ind, bnd in enumerate(AppInf['bands']):
        gstr = AppInf['Groups'][i];
        maxvalTP_PSD.append(AppInf['%s_maxvalTP(%s)'%(gstr,bnd)])
    maxvalTP_PSD=np.array(maxvalTP_PSD)
    bbb=zeros((11,11))
    for a in range(11):
        for b in range(11):
            aaa=list(maxvalTP_PSD[:,a,b])
            bbb[a,b]=aaa.index(max(aaa))
    ax=subplot(lengrp,1,i+1);
    if not i: 
        title("Frekans max. PSD",fontsize=titfs);
    my_adjust_contourf(sig_Val,p_Val,bbb,cmap=cmap,alpha=alpha,fsize=ticksfs,levels=array(range(len(AppInf['bands'])-1))); 
    ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
    ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_freqbands_maxvalTP_PSD)' + figformat1,dpi=dpi1, bbox_inches='tight')

# farklı frekans bandları: toplam LFP ye göre SNR

levels_SRfbands=[[None]*lengrp for a in range(len(AppInf['bands']))];
"""
levels_SRfbands=[[] for a in range(len(AppInf['bands']))]; 
levels_SRfbands[0].append(arange(0,561,80)); 
levels_SRfbands[0].append(arange(0,561,80))
levels_SRfbands[0].append(arange(-1,1,1)); 
levels_SRfbands[0].append(arange(-1,1,1))

levels_SRfbands[1].append(arange(100,521,60)); 
levels_SRfbands[1].append(arange(80,501,60))
levels_SRfbands[1].append(arange(300,1876,225)); 
levels_SRfbands[1].append(arange(200,1351,150))

levels_SRfbands[2].append(arange(20,63,6)); 
levels_SRfbands[2].append(arange(20,50,4))
levels_SRfbands[2].append(arange(60,105,6)); 
levels_SRfbands[2].append(arange(30,60,4))

levels_SRfbands[3].append(arange(12,48,5)); 
levels_SRfbands[3].append(arange(15,37,3))
levels_SRfbands[3].append(arange(20,71,7)); 
levels_SRfbands[3].append(arange(12,34,3))
"""
# AppInf['bands']=AppInf['bands'][3:]
lenbnd = len(AppInf['bands'])
fig=figure(figsize=(10*lenbnd,8*lengrp));
suptitle(textwrap.dedent(fig_titles[15]).strip(),fontsize=titfs);
for bndindx in range(lenbnd):
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        bnd = AppInf['bands'][bndindx]
        AppInf['%s_snrTP(%s)'%(gstr,bnd)] = reshape(get_data_from_dict(AppRes_Lst,'%s_snrTP(%s)'%(gstr,bnd)),myshape).T;
        ax=subplot2grid((lengrp,lenbnd),(i,bndindx));
        my_adjust_contourf(sig_Val,p_Val,AppInf['%s_snrTP(%s)'%(gstr,bnd)], \
                cmap=cmap,alpha=alpha,fsize=ticksfs,levels=levels_SRfbands[bndindx][i]); 
        ax.tick_params(axis='both', labelbottom='off',labelleft='off')
        if not i: 
            title('(%s)'%(bnd),fontsize=titfs);
        if (i==(lengrp-1)): 
            xlabel(r'$\sigma_g$',fontsize=lblfs);
            ax.tick_params(axis='x', labelbottom='on')
            ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
        if not bndindx: 
            ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
            ax.tick_params(axis='y', labelleft='on')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_freqbands_StocRes)' + figformat1,dpi=dpi1, bbox_inches='tight')

# farklı frekans bandları: toplam LFP ye göre CohRes

levels_CRfbands=[[None]*lengrp for a in range(len(AppInf['bands']))];
"""
# talamokortikal ağ için ayarlı
levels_CRfbands=[[] for a in range(len(AppInf['bands']))]; 
levels_CRfbands[0].append(arange(0,2e3,2e2)); levels_CRfbands[0].append(arange(0,3e3,3e2))
levels_CRfbands[1].append(arange(0,1e4,1e3)); levels_CRfbands[1].append(arange(0,1e4,1e3))
levels_CRfbands[2].append(arange(0,4e4,4e3)); levels_CRfbands[2].append(arange(0,3e4,2.5e3))
levels_CRfbands[3].append(arange(0,7.5e4,7.5e3)); levels_CRfbands[3].append(arange(0,4e4,4e3))
"""
"""
# talamik ağ için ayarlı
levels_CRfbands=[[] for a in range(len(AppInf['bands']))]; 
levels_CRfbands[0].append(arange(-1,1,1)); levels_CRfbands[0].append(arange(-1,1,1))
levels_CRfbands[1].append(arange(0,2.4e5,3e4)); levels_CRfbands[1].append(arange(0,1.5e5,2e4))
levels_CRfbands[2].append(arange(0,4e4,5e3)); levels_CRfbands[2].append(arange(0,1.5e4,2e3))
levels_CRfbands[3].append(arange(0,3e4,4e3)); levels_CRfbands[3].append(arange(0,1.5e4,2e3))
"""
lenbnd = len(AppInf['bands'])
fig=figure(figsize=(10*lenbnd,8*lengrp));
suptitle(textwrap.dedent(fig_titles[16]).strip(),fontsize=titfs);
for bndindx in range(lenbnd):
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        bnd = AppInf['bands'][bndindx]
        AppInf['%s_CFactor_betaTP(%s)'%(gstr,bnd)] = reshape(get_data_from_dict(AppRes_Lst,'%s_CFactor_betaTP(%s)'%(gstr,bnd)),myshape).T;
        ax=subplot2grid((lengrp,lenbnd),(i,bndindx));
        my_adjust_contourf(sig_Val,p_Val,AppInf['%s_CFactor_betaTP(%s)'%(gstr,bnd)],cmap=cmap,alpha=alpha,fsize=ticksfs,levels=levels_CRfbands[bndindx][i]); 
        ax.tick_params(axis='both', labelbottom='off',labelleft='off')
        if not i: 
            title('(%s)'%(bnd),fontsize=titfs);
        if (i==(lengrp-1)): 
            xlabel(r'$\sigma_g$',fontsize=lblfs); 
            ax.tick_params(axis='x', labelbottom='on')
            ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
        if not bndindx: 
            ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
            ax.tick_params(axis='y', labelleft='on')
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_freqbands_CohRes)' + figformat1,dpi=dpi1, bbox_inches='tight')

"""
# toplam LFP ye göre (farklı frekans bandları)
for bnd in AppInf['bands']:
    fig=figure(figsize=(16,6*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        AppInf['%s_snrTP(%s)'%(gstr,bnd)] = reshape(get_data_from_dict(AppRes_Lst,'%s_snrTP(%s)'%(gstr,bnd)),myshape).T;
        AppInf['%s_CFactor_betaTP(%s)'%(gstr,bnd)] = reshape(get_data_from_dict(AppRes_Lst,'%s_CFactor_betaTP(%s)'%(gstr,bnd)),myshape).T;
        subplot2grid((lengrp,2),(i,0));
        my_adjust_contourf(sig_Val,p_Val,AppInf['%s_snrTP(%s)'%(gstr,bnd)],cmap=cmap,alpha=alpha,fsize=ticksfs); 
        ax.tick_params(axis='both', labelbottom='off',labelleft='off')        
        if not i: 
            title(textwrap.dedent(fig_titles[15]).strip()+'(%s)'%(bnd),fontsize=titfs);
        ylabel('%s\n$p$'%(gstr),fontsize=lblfs);
        ax.tick_params(axis='y', labelleft='on')
        if (i==(lengrp-1)): 
            xlabel(r'$\sigma_g$',fontsize=lblfs);
            ax.tick_params(axis='x', labelbottom='on')
            ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
        subplot2grid((lengrp,2),(i,1)); 
        my_adjust_contourf(sig_Val,p_Val,AppInf['%s_CFactor_betaTP(%s)'%(gstr,bnd)],cmap=cmap,alpha=alpha,fsize=ticksfs);
        ax.tick_params(axis='both', labelbottom='off',labelleft='off')
        if not i: 
            title(textwrap.dedent(fig_titles[16]).strip()+'(%s)'%(bnd),fontsize=titfs);
        if (i==(lengrp-1)): 
            xlabel(r'$\sigma_g$',fontsize=lblfs);
            ax.tick_params(axis='x', labelbottom='on')
            ax.set_xticklabels(ax.get_xticks().tolist(),rotation='vertical')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_StocRes-CohRes(%s))'%(bnd) + figformat1,dpi=dpi1, bbox_inches='tight')
"""
if not recPic : show()


####### sabit bir p/sig değerine bağlı farklı SW/gürültü etkisi #######
#sbtprm = 0;
figKat=5; nofticks=5;
titfs=28; lblfs=30; ticksfs=28;

for sbtprm in range(2):
    if not sbtprm: # SW parametresi (p) sabit
        app=sig_Val; p_indx=AppInf['sel_p_indx']; sig_indx=range(len(app)); 
        idq_str= 'Sig'; idq_lab=r'$\sigma_g$'; inf_str = '' # '\n(p=%s)'%str(array(p_Val)[p_indx])
    else: # gürültü parametresi (sig) sabit
        app=p_Val; p_indx=range(len(app)); sig_indx=AppInf['sel_sig_indx']; 
        idq_str='p'; idq_lab=r'$p$'; inf_str = '' # '\n(Sig=%s)'%str(array(sig_Val)[sig_indx])
    
    # uygulamaların aktivite durumları
    fig=figure(figsize=(30,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        
        ax=subplot(lengrp,3,3*i+1);
        if not i: title(textwrap.dedent(fig_titles[0]+inf_str).strip(),fontsize=titfs)
        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_spkCount'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_spkCount'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);        
        ax.tick_params(labelbottom='off')
        #ax.locator_params(axis='y', nbins = nofticks)
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);

        ax=subplot(lengrp,3,3*i+2); 
        if not i: title(textwrap.dedent(fig_titles[1]+inf_str).strip(),fontsize=titfs);
        vals = []; legends=[]     
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_passCount'%(gstr)][p_indx[a],sig_indx])
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_passCount'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);       
        ax.tick_params(labelbottom='off')
        #ax.locator_params(axis='y', nbins = nofticks)
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        
        ax=subplot(lengrp,3,3*i+3);
        if not i: title(textwrap.dedent(fig_titles[7]+inf_str).strip(),fontsize=titfs);
        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_meanFreq'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_meanFreq'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);        
        ax.tick_params(labelbottom='off')
        #ax.locator_params(axis='y', nbins = nofticks)
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        #ylabel('%s'%(gstr),fontsize=lblfs);
        
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_SpkCnt-PassCnt-Freq)'+'-'+idq_str+ figformat1,dpi=dpi1, bbox_inches='tight')
    
    # Voltaj senk. ve Spike uyumu
    fig=figure(figsize=(20,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        ax=subplot(lengrp,2,2*i+1);
        if not i: title(textwrap.dedent(fig_titles[2]+inf_str).strip(),fontsize=titfs);

        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_VSync'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_VSync'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        ax=subplot(lengrp,2,2*i+2); 
        if not i: title(textwrap.dedent(fig_titles[3]+inf_str).strip(),fontsize=titfs);

        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_SpikeCoh'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_SpikeCoh'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);        
        
        tick_params(labelbottom='off')    
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_VmSync-PopCoh)'+'-'+idq_str+ figformat1,dpi=dpi1, bbox_inches='tight')
    
    # ortalama VS, CV ve CC
    fig=figure(figsize=(30,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        subplot(lengrp,3,3*i+1);
        if not i: title(textwrap.dedent(fig_titles[4]+inf_str).strip(),fontsize=titfs);

        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_meanVS'%(gstr)][p_indx[a],sig_indx]/AppInf['%s_meanVS'%(gstr)].max());
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_meanVS'%(gstr)][p_indx,sig_indx[a]]/AppInf['%s_meanVS'%(gstr)].max())
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        subplot(lengrp,3,3*i+2);
        if not i: title(textwrap.dedent(fig_titles[5]+inf_str).strip(),fontsize=titfs);

        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_meanCV'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_meanCV'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        
        tick_params(labelbottom='off')    
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        subplot(lengrp,3,3*i+3); 
        if not i: title(textwrap.dedent(fig_titles[6]+inf_str).strip(),fontsize=titfs);

        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_meanCC'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_meanCC'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        
        tick_params(labelbottom='off')    
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_VS-CV-CC)'+'-'+idq_str+ figformat1,dpi=dpi1, bbox_inches='tight')
    
    # ortalama frekans, faz uyum ve Burst senk.
    fig=figure(figsize=(20,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        subplot(lengrp,2,2*i+1); 
        if not i: title(textwrap.dedent(fig_titles[8]+inf_str).strip(),fontsize=titfs);

        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_meanPhase'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_meanPhase'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        subplot(lengrp,2,2*i+2); 
        if not i: title(textwrap.dedent(fig_titles[9]+inf_str).strip(),fontsize=titfs);

        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_SyncBurst'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_SyncBurst'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_PhsCoh-BrstSync)'+'-'+idq_str + figformat1,dpi=dpi1, bbox_inches='tight')
    
    # Stokastik ve Koherens rezonans
    # toplam memb. voltaj'a göre
    fig=figure(figsize=(20,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        subplot(lengrp,2,2*i+1);
        if not i: title(textwrap.dedent(fig_titles[10]+inf_str).strip(),fontsize=titfs);

        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_snrTM1'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_snrTM1'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        subplot(lengrp,2,2*i+2); 
        if not i: title(textwrap.dedent(fig_titles[11]+inf_str).strip(),fontsize=titfs);

        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_CFactor_betaTM1'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_CFactor_betaTM1'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_StocRes1-CohRes1)'+'-'+idq_str + figformat1,dpi=dpi1, bbox_inches='tight')
    
    # toplam LFP ye göre
    fig=figure(figsize=(20,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        subplot(lengrp,2,2*i+1);
        if not i: title(textwrap.dedent(fig_titles[15]+inf_str).strip(),fontsize=titfs);

        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_snrTP1'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_snrTP1'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        subplot(lengrp,2,2*i+2); 
        if not i: title(textwrap.dedent(fig_titles[16]+inf_str).strip(),fontsize=titfs);

        vals = []; legends=[]
        if not sbtprm:
            for a in range(len(p_indx)):
                vals.append(AppInf['%s_CFactor_betaTP1'%(gstr)][p_indx[a],sig_indx]);
                legends.append('$p=%s$'%(p_Val[p_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        else:
            for a in range(len(sig_indx)):
                vals.append(AppInf['%s_CFactor_betaTP1'%(gstr)][p_indx,sig_indx[a]])
                legends.append('$\sigma_g=%s$'%(sig_Val[sig_indx[a]]))
            my_adjust_plot2(app,vals,leg=legends,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_StocRes2-CohRes2)'+'-'+idq_str + figformat1,dpi=dpi1, bbox_inches='tight')

if not recPic : show()
#show()
matplotlib.pylab.close('all')

####### p/sig değerine bağlı  SW/gürültü etkisi istatistikleri #######
for sbtprm in range(2):
    if not sbtprm: # SW parametresi (p) sabit
        app=sig_Val; p_indx=AppInf['sel_p_indx']; sig_indx=range(len(app)); 
        idq_str= 'Sig'; idq_lab=r'$\sigma_g$'; inf_str = '' # '\n(p=%s)'%str(array(p_Val)[p_indx])
    else: # gürültü parametresi (sig) sabit
        app=p_Val; p_indx=range(len(app)); sig_indx=AppInf['sel_sig_indx']; 
        idq_str= 'p'; idq_lab=r'$p$'; inf_str = '' # '\n(Sig=%s)'%str(array(sig_Val)[sig_indx])
    
    # uygulamaların aktivite durumları
    fig=figure(figsize=(30,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        
        ax=subplot(lengrp,3,3*i+1);
        if not i: title(textwrap.dedent(fig_titles[0]+inf_str).strip(),fontsize=titfs)
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_spkCount'%(gstr)],axis=0)
            #err=np.std(AppInf['%s_spkCount'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_spkCount'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

        ax=subplot(lengrp,3,3*i+2); 
        if not i: title(textwrap.dedent(fig_titles[1]+inf_str).strip(),fontsize=titfs);

        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_passCount'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_passCount'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);       
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        
        ax=subplot(lengrp,3,3*i+3);
        if not i: title(textwrap.dedent(fig_titles[7]+inf_str).strip(),fontsize=titfs);
        vals = []; legends=[]
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_meanFreq'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_meanFreq'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        #ylabel('%s'%(gstr),fontsize=lblfs)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_SpkCnt-PassCnt-Freq-Stat)'+'-'+idq_str+ figformat1,dpi=dpi1, bbox_inches='tight')

    # Voltaj senk. ve Spike uyumu
    fig=figure(figsize=(20,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        subplot(lengrp,2,2*i+1);
        if not i: title(textwrap.dedent(fig_titles[2]+inf_str).strip(),fontsize=titfs);
        
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_VSync'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_VSync'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        subplot(lengrp,2,2*i+2); 
        if not i: title(textwrap.dedent(fig_titles[3]+inf_str).strip(),fontsize=titfs);

        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_SpikeCoh'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_SpikeCoh'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);        
        
        tick_params(labelbottom='off')    
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_VmSync-PopCoh-stat)'+'-'+idq_str+ figformat1,dpi=dpi1, bbox_inches='tight')
    
    # ortalama VS, CV ve CC
    fig=figure(figsize=(30,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        subplot(lengrp,3,3*i+1);
        if not i: title(textwrap.dedent(fig_titles[4]+inf_str).strip(),fontsize=titfs);
        
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_meanVS'%(gstr)]/AppInf['%s_meanVS'%(gstr)].max(),axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_meanVS'%(gstr)]/AppInf['%s_meanVS'%(gstr)].max(),axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        subplot(lengrp,3,3*i+2);
        if not i: title(textwrap.dedent(fig_titles[5]+inf_str).strip(),fontsize=titfs);
        
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_meanCV'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_meanCV'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        
        tick_params(labelbottom='off')    
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        subplot(lengrp,3,3*i+3); 
        if not i: title(textwrap.dedent(fig_titles[6]+inf_str).strip(),fontsize=titfs);
        
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_meanCC'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_meanCC'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        
        tick_params(labelbottom='off')    
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_VS-CV-CC-stat)'+'-'+idq_str+ figformat1,dpi=dpi1, bbox_inches='tight')
    
    # ortalama frekans, faz uyum ve Burst senk.
    fig=figure(figsize=(20,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        subplot(lengrp,2,2*i+1); 
        if not i: title(textwrap.dedent(fig_titles[8]+inf_str).strip(),fontsize=titfs);
        
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_meanPhase'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_meanPhase'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        subplot(lengrp,2,2*i+2); 
        if not i: title(textwrap.dedent(fig_titles[9]+inf_str).strip(),fontsize=titfs);
        
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_SyncBurst'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_SyncBurst'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_PhsCoh-BrstSync-stat)'+'-'+idq_str + figformat1,dpi=dpi1, bbox_inches='tight')
    
    # Stokastik ve Koherens rezonans
    # toplam memb. voltaj'a göre
    fig=figure(figsize=(20,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        subplot(lengrp,2,2*i+1);
        if not i: title(textwrap.dedent(fig_titles[10]+inf_str).strip(),fontsize=titfs);
        
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_snrTM1'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_snrTM1'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        subplot(lengrp,2,2*i+2); 
        if not i: title(textwrap.dedent(fig_titles[11]+inf_str).strip(),fontsize=titfs);
        
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_CFactor_betaTM1'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_CFactor_betaTM1'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_StocRes1-CohRes1-stat)'+'-'+idq_str + figformat1,dpi=dpi1, bbox_inches='tight')
    
    # toplam LFP ye göre
    fig=figure(figsize=(20,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i];
        ax=subplot(lengrp,2,2*i+1);
        if not i: title(textwrap.dedent(fig_titles[15]+inf_str).strip(),fontsize=titfs);
        
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_snrTP1'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_snrTP1'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        
        ax=subplot(lengrp,2,2*i+2); 
        if not i: title(textwrap.dedent(fig_titles[16]+inf_str).strip(),fontsize=titfs);
        
        if not sbtprm:
            # sigma'ya göre değişen ortalama
            vals=np.mean(AppInf['%s_CFactor_betaTP1'%(gstr)],axis=0)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        else:
            # p'ye göre değişen ortalama
            vals=np.mean(AppInf['%s_CFactor_betaTP1'%(gstr)],axis=1)
            my_adjust_plot1(app,vals,fsize=ticksfs);
        
        tick_params(labelbottom='off')
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        
        if (i==(lengrp-1)): 
            xlabel(idq_lab,fontsize=lblfs); tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Apps_StocRes2-CohRes2)'+'-'+idq_str + figformat1,dpi=dpi1, bbox_inches='tight')

matplotlib.pylab.close('all')


############### seçilen p ve sig değerleri için grafikler ##################
TrialTitles = ["Zaman Frekans\n($p=%s, sig=%s$)","Uzay Frekans (Hz)\n($p=%s, sig=%s$)",\
                "Uzay Frekans Dagılım\n($p=%s, sig=%s$)","VS\n($p=%s, sig=%s$)",\
                "VS Dagılım\n($p=%s, sig=%s$)", \
                "Polar VS\n($p=%s, sig=%s$)","Spike Faz Dagılım\n($p=%s, sig=%s$)",\
                "Noron Faz Oruntu\n($p=%s, sig=%s$)",\
                "Gerilim Degisimi (mV)\n($p=%s, sig=%s$)","Raster Grafik\n($p=%s, sig=%s$)",\
                "Toplam Membran (TM)\n($p=%s, sig=%s$)","TM PSD\n($p=%s, sig=%s$)","TM Spektogram\n($p=%s, sig=%s$)",
                "Toplam LFP\n($p=%s, sig=%s$)","LFP PSD\n($p=%s, sig=%s$)","LFP Spektogram\n($p=%s, sig=%s$)",
                "Vm Oruntu\n($p=%s, sig=%s$)","Vm Oruntu($%s ms$)\n($p=%s, sig=%s$)", \
                "CC Harita\n($p=%s, sig=%s$)","Spike Harita\n($p=%s, sig=%s$)",\
                "Faz Uyumu ($R=%s$)\n($p=%s, sig=%s$)"]

trialSel=AppInf['trialSel']; # (pval,sigval)
ticksfs = 30; lblfs = 36; titfs=36; textfs = 30
lwidth=2; nofticks=4


# uzaya ve zamana göre spike frekansları
for x in range(len(trialSel)):
    tsel=trialSel[x]; 
    trial_idq='_(%s - %s)'%(str(tsel[0]).replace('.','_'),str(tsel[1]).replace('.','_'))
    fig=figure(figsize=(30,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i]; appindx = getAppindx(p_Val,sig_Val,tsel)
        ax=subplot(lengrp,3,3*i+1); 
        ax.locator_params(axis='y', nbins = nofticks)
        if not i: title(textwrap.dedent(TrialTitles[0]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        rate_time=AppRes_Lst[appindx]['MV_%s_times'%(gstr)]/second;
        rates=AppRes_Lst[appindx]['M_%s_rateSR2ms'%(gstr)];
        ratemx=AppInf['%s_MnMxrateSR2ms'%(gstr)][1]
        my_adjust_plot1(rate_time,\
            rates,fsize=ticksfs, \
            lwidth=lwidth,marker='',glwidth=1,yminmax=AppInf['%s_MnMxrateSR2ms'%(gstr)]);
        #tick_params(labelbottom='off')
        text(max(rate_time),ratemx-ratemx*0.1, \
            'Ort.=%s - STD=%s'%(round(mean(rates),2),round(std(rates),2)),va='bottom',ha='right')
        if (i==(lengrp-1)):
            xlabel('Zaman (saniye)',fontsize=lblfs); #tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        ax=subplot(lengrp,3,3*i+2);
        ax.locator_params(axis='y', nbins = nofticks)
        if not i: title(textwrap.dedent(TrialTitles[1]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        N=int(AppRes_Lst[appindx]['N_%s'%(gstr)]);
        netf=AppRes_Lst[appindx]['%s_netfreq'%(gstr)]
        netfmx=AppInf['%s_MnMxnetFreq'%(gstr)][1]
        my_adjust_plot1(arange(N),netf,fsize=ticksfs, \
            lwidth=lwidth,marker='',glwidth=1,yminmax=AppInf['%s_MnMxnetFreq'%(gstr)],linepoint=mean(netf));
        #tick_params(labelbottom='off')
        text(N,netfmx-netfmx*0.1, \
            'Ort.=%s - STD=%s'%(round(mean(netf),2),round(std(netf),2)),va='bottom',ha='right')
        if (i==(lengrp-1)): 
            xlabel('Hucre ID',fontsize=lblfs); #tick_params(labelbottom='on')
        #xminmax=get_minmax_from_dict(AppRes_Lst,'N_%s'%(gstr))
        #yminmax=get_minmax_from_dict(AppRes_Lst,'%s_netfreq'%(gstr))
        ax=subplot(lengrp,3,3*i+3);
        ax.locator_params(axis='y', nbins = nofticks)
        if not i: title(textwrap.dedent(TrialTitles[2]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        r=my_adjust_hist(netf,bin=10,fsize=ticksfs,\
        xminmax=AppInf['%s_MnMxnetFreq'%(gstr)],yminmax=(0,N*3/5),linepoint=mean(netf))
        #tick_params(labelbottom='off')
        #text(AppInf['%s_MnMxnetFreq'%(gstr)][1],int(N*3/5-(N*3/5)*0.1), \
        #    'STD=%s'%(round(std(r[0]),2)),va='bottom',ha='right')
        if (i==(lengrp-1)): 
            xlabel('Frekans (Hz)',fontsize=lblfs); #tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(SpikeRateAnaliz)'+trial_idq + figformat1,dpi=dpi1, bbox_inches='tight')

# herbir hücrenin VS değeri
for x in range(len(trialSel)):
    tsel=trialSel[x]; 
    trial_idq='_(%s - %s)'%(str(tsel[0]).replace('.','_'),str(tsel[1]).replace('.','_'))
    fig=figure(figsize=(20,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i]; appindx = getAppindx(p_Val,sig_Val,tsel)
        ax=subplot(lengrp,2,2*i+1);
        ax.locator_params(axis='y', nbins = nofticks)
        if not i: title(textwrap.dedent(TrialTitles[3]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        N=int(AppRes_Lst[appindx]['N_%s'%(gstr)]);
        VS_vals = AppRes_Lst[appindx]['%s_VS'%(gstr)]/AppInf['%s_MnMxVS'%(gstr)][1]
        VSmx=AppInf['%s_MnMxnetFreq'%(gstr)][1]
        my_adjust_plot1(arange(N),VS_vals,fsize=ticksfs,lwidth=lwidth, \
            marker='',glwidth=1,yminmax=(0,1),linepoint=mean(VS_vals)) 
        #tick_params(labelbottom='off')
        text(N,1-1*0.1, \
            'Ort.=%s - STD=%s'%(round(mean(VS_vals),2),round(std(VS_vals),2)),va='bottom',ha='right')
        if (i==(lengrp-1)): 
            xlabel('Hucre ID',fontsize=lblfs); #tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        ax=subplot(lengrp,2,2*i+2); 
        ax.locator_params(axis='y', nbins = nofticks)
        if not i: title(textwrap.dedent(TrialTitles[4]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        r=my_adjust_hist(VS_vals,bin=10,fsize=ticksfs, \
        xminmax=(0,1),yminmax=(0,AppRes_Lst[appindx]['N_%s'%(gstr)]*1/2),linepoint=mean(VS_vals))
        #tick_params(labelbottom='off')
        #text(1,N*1/2-(N*1/2)*0.1, \
        #    'STD=%s'%(round(std(r[0]),2)),va='bottom',ha='right')
        if (i==(lengrp-1)): 
            xlabel('VS',fontsize=lblfs); #tick_params(labelbottom='on')
        #ylabel('Sayi',fontsize=lblfs);
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(VSAnaliz)'+trial_idq + figformat1,dpi=dpi1, bbox_inches='tight')

# radial VS gösterim
#ara=10;
polar=False;

figKat=6; nofticks=5;
titfs=32; lblfs=30; ticksfs=30;

angle_ticks = linspace(0,2*pi,7)
angle_ticks_lbl = [r"$0$", r"$\frac{\pi}{3}$", r"$\frac{2\pi}{3}$", r"$\pi$", r"$\frac{4\pi}{3}$", r"$\frac{5\pi}{3}$", r"$2\pi$"]
angle_ticks_lbl_polar = [r"$0$", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$", r"$\frac{3\pi}{4}$", r"$\pi$", r"$\frac{5\pi}{4}$", r"$\frac{3\pi}{2}$",  r"$\frac{7\pi}{4}$"]

for x in range(len(trialSel)):
    tsel=trialSel[x]; 
    trial_idq='_(%s - %s)'%(str(tsel[0]).replace('.','_'),str(tsel[1]).replace('.','_'))
    fig=figure(figsize=(30,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i]; appindx = getAppindx(p_Val,sig_Val,tsel)
        
        ax=subplot(lengrp,3,3*i+1,polar=True); 
        ara=AppRes_Lst[appindx]['%s_spkCount'%(gstr)]/100+1; # 1:100 atlama
        if not i: title(textwrap.dedent(TrialTitles[5]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        X=AppRes_Lst[appindx]['%s_VSX'%(gstr)];
        Y=AppRes_Lst[appindx]['%s_VSY'%(gstr)];
        Theta=AppRes_Lst[appindx]['%s_VSTheta'%(gstr)];
        h=my_adjust_polar(X,Y,Theta,ara=ara,fsize=ticksfs)
        ax.set_xticks(arange(0,2*pi,pi/4)); ax.set_xticklabels(angle_ticks_lbl_polar, fontsize=ticksfs+10)
        #ylabel('%s'%(gstr),fontsize=lblfs);
        text(pi, 1.5, '%s'%(gstr), fontsize=lblfs)
        #ax.tick_params(axis='x', size=lblfs)
        
        ax=subplot(lengrp,3,3*i+2,polar=polar); 
        ax.locator_params(axis='y', nbins = nofticks)
        if not i: title(textwrap.dedent(TrialTitles[6]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        bins=linspace(0,2*pi,13);
        my_adjust_hist(array(list(flatten(Theta))[::]),bin=bins,fsize=ticksfs,maxHist = 2*pi,\
            yminmax=AppInf['%s_MnMxtotVStheta'%(gstr)])
        ax.set_xticks(angle_ticks); ax.set_xticklabels(angle_ticks_lbl,fontsize=ticksfs+10)
        #ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        ax.tick_params(axis='y', labelsize=ticksfs)
        ax.tick_params(axis='x', labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel(r'$\theta$',fontsize=lblfs);
            ax.tick_params(axis='x', labelbottom='on')
        #ylabel('Sayi',fontsize=lblfs);
        
        ax=subplot(lengrp,3,3*i+3);
        if not i: title(textwrap.dedent(TrialTitles[7]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        N=int(AppRes_Lst[appindx]['N_%s'%(gstr)])
        my_adjust_imshow(AppRes_Lst[appindx]['%s_VSThetaM'%(gstr)],extent=(0, 2*pi, 0, N-1),zminmax=AppInf['%s_MnMxVStheta'%(gstr)])
        #ylabel('%s'%(gstr),fontsize=lblfs);
        ax.set_xticks(angle_ticks); ax.set_xticklabels(angle_ticks_lbl,fontsize=ticksfs+10)
        ax.tick_params(axis='y', labelsize=ticksfs)
        ax.tick_params(axis='x', labelbottom='off')
        if (i==(lengrp-1)):
            xlabel(r'$\theta$',fontsize=lblfs);
            ax.tick_params(axis='x', labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(VSPhaseAnaliz)'+trial_idq+ figformat1,dpi=dpi1, bbox_inches='tight')

##################################### 3D çizimler ###############################
# 3d gerilim değişimleri ile uzaya ve zamana göre spike frekansları
# AppInf['trialSeld'] = [(0.1,0.08),(0.1,0.04),(0.8,0.08),(0.8,0.04)]
cmap=None
alpha = 0.7
fsize=24; titfs3d = 30
lw=2;

fig=figure(figsize=(32,figKat*lengrp));

for i in range(lengrp):
    yticklabels=[]
    data1=[]; data2=[]; data3=[]; data4=[]; 
    gstr = AppInf['Groups'][i];
    
    for x in range(len(AppInf['trialSeld'])):
        tsel=AppInf['trialSeld'][x]; 
        trial_idq=r'$(%s - %s)$'%(tsel[0],tsel[1])
        yticklabels.append(trial_idq)
        appindx = getAppindx(p_Val,sig_Val,tsel)
        N=int(AppRes_Lst[appindx]['N_%s'%(gstr)])
        
        xs4=AppRes_Lst[appindx]['MV_%s_times'%(gstr)]/second
        data4.append(AppRes_Lst[appindx]['MV_%s'%(gstr)][N*3/4-1+6]/mV)  
        
        xs1=AppRes_Lst[appindx]['MV_%s_times'%(gstr)]/second
        data1.append(AppRes_Lst[appindx]['M_%s_rateSR2ms'%(gstr)])   
        
        xs2=arange(N)
        data2.append(AppRes_Lst[appindx]['%s_netfreq'%(gstr)])  
            
    ax = fig.add_subplot(lengrp,4,4*i+1, projection='3d')
    ax.locator_params(axis='z', nbins = nofticks)
    if not i: title("Gerilim Degisimi",fontsize=titfs3d);    
    my_plot_3D(xs4,data4,fig,ax=ax,yticklabels=yticklabels,alpha=alpha,cmap=cmap,fsize=fsize)
    ax.set_xlabel('Zaman (s)',va='center', ha='left',fontsize=fsize); 
    ax.set_zlabel('mV', rotation=90,va='center', ha='left',fontsize=fsize)
    ax.text2D(0.05, 0.8, gstr, transform=ax.transAxes,fontsize=titfs3d)
    
    ax = fig.add_subplot(lengrp,4,4*i+2, projection='3d')
    ax.locator_params(axis='z', nbins = nofticks)
    if not i: title("Zaman Frekans",fontsize=titfs3d);
    my_plot_3D(xs1,data1,fig,ax=ax,yticklabels=yticklabels,alpha=alpha,cmap=cmap,fsize=fsize)
    ax.set_xlabel('Zaman (s)',va='center', ha='left',fontsize=fsize); 
    ax.set_zlabel('Hz', rotation=90,va='center', ha='left',fontsize=fsize)
    #ax.text2D(0.05, 0.8, gstr, transform=ax.transAxes)
    
    ax = fig.add_subplot(lengrp,4,4*i+3, projection='3d')
    ax.locator_params(axis='z', nbins = nofticks)
    if not i: title("Uzay Frekans (Hz)",fontsize=titfs3d);
    my_poly_3D(xs2,data2,fig,ax=ax,yticklabels=yticklabels,alpha=alpha,cmap=cmap,fsize=fsize)
    ax.set_xlabel('ID',va='center', ha='left',fontsize=fsize); 
    ax.set_zlabel('Sayi', rotation=90,va='center', ha='left',fontsize=fsize)
    #ax.text2D(0.05, 0.8, gstr, transform=ax.transAxes)
    
    ax = fig.add_subplot(lengrp,4,4*i+4, projection='3d')
    ax.locator_params(axis='z', nbins = nofticks)
    if not i: title("Uzay Frekans Dag.",fontsize=titfs3d);
    my_hist_3D(data2,fig,ax=ax,yticklabels=yticklabels,alpha=alpha,cmap=cmap,fsize=fsize)
    ax.set_xlabel('Hz',va='center', ha='left',fontsize=fsize); 
    ax.set_zlabel('Sayi', rotation=90,va='center', ha='left',fontsize=fsize)
    #ax.text2D(0.05, 0.8, gstr, transform=ax.transAxes)
fig.tight_layout()
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Vm_SpikeRateAnaliz_3D)' + figformat1,dpi=dpi1) #, bbox_inches='tight')

# 3d VS analizler
fig=figure(figsize=(24,figKat*lengrp));

for i in range(lengrp):
    yticklabels=[]
    data1=[]; data2=[]; data3=[]; data4=[]; 
    gstr = AppInf['Groups'][i];
    
    for x in range(len(AppInf['trialSeld'])):
        tsel=AppInf['trialSeld'][x]; 
        trial_idq=r'$(%s - %s)$'%(tsel[0],tsel[1])
        yticklabels.append(trial_idq)
        appindx = getAppindx(p_Val,sig_Val,tsel)
        N=int(AppRes_Lst[appindx]['N_%s'%(gstr)])
        
        xs1=arange(N)
        data1.append(AppRes_Lst[appindx]['%s_VS'%(gstr)]/AppInf['%s_MnMxVS'%(gstr)][1])  
        
        Theta = AppRes_Lst[appindx]['%s_VSTheta'%(gstr)]
        data2.append(array(list(flatten(Theta))[::]))
            
    ax = fig.add_subplot(lengrp,3,3*i+1, projection='3d')
    ax.locator_params(axis='z', nbins = nofticks)
    if not i: title("VS Degerleri",fontsize=titfs3d);    
    my_poly_3D(xs1,data1,fig,ax=ax,yticklabels=yticklabels,alpha=alpha,cmap=cmap,fsize=fsize,lw=lw)
    ax.set_xlabel('ID',va='center', ha='left',fontsize=fsize); 
    ax.set_zlabel(r'$VS_{ort.}$', rotation=90,va='center', ha='left',fontsize=fsize)
    ax.text2D(0.05, 0.8, gstr, transform=ax.transAxes,fontsize=titfs3d)
    
    ax = fig.add_subplot(lengrp,3,3*i+2, projection='3d')
    ax.locator_params(axis='z', nbins = nofticks)
    if not i: title("VS Dagılımları",fontsize=titfs3d);
    my_hist_3D(data1,fig,ax=ax,yticklabels=yticklabels,alpha=alpha,cmap=cmap,fsize=fsize)
    ax.set_xlabel(r'$VS_{ort.}$',va='center', ha='left',fontsize=fsize); 
    ax.set_zlabel('Sayi', rotation=90,va='center', ha='left',fontsize=fsize)
    #ax.text2D(0.05, 0.8, gstr, transform=ax.transAxes)
    
    ax = fig.add_subplot(lengrp,3,3*i+3, projection='3d')
    ax.locator_params(axis='z', nbins = nofticks)
    if not i: title("Faz Dagılımları",fontsize=titfs3d);
    my_hist_3D(data2,fig,ax=ax,yticklabels=yticklabels,alpha=alpha,cmap=cmap,fsize=fsize)
    ax.set_xlabel(r'$\theta$',va='center', ha='left',fontsize=fsize); 
    ax.set_zlabel('Sayi', rotation=90,va='center', ha='left',fontsize=fsize)
    ax.set_xticks(angle_ticks); ax.set_xticklabels(angle_ticks_lbl,fontsize=ticksfs)
    #ax.text2D(0.05, 0.8, gstr, transform=ax.transAxes)
fig.tight_layout()
if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(VSAnaliz_3D)' + figformat1,dpi=dpi1) #, bbox_inches='tight')

matplotlib.pylab.close('all')

"""
fig= figure()
ax = fig.add_subplot(1,1,1, projection='3d')
my_hist_3D(data3,fig,ax=ax,yticklabels=yticklabels,alpha=alpha,cmap=cmap,fsize=fsize)

fig= figure()
ax = fig.add_subplot(1,1,1, projection='3d')
my_poly_3D(xs1,data1,fig,ax=ax,yticklabels=yticklabels,alpha=alpha,cmap=cmap,fsize=fsize,lw=lw)

show()
"""
###############################################################################

# hücre gerilimleri 
for x in range(len(trialSel)):
    tsel=trialSel[x]; 
    trial_idq='_(%s - %s)'%(str(tsel[0]).replace('.','_'),str(tsel[1]).replace('.','_'))
    fig=figure(figsize=(20,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i]; appindx = getAppindx(p_Val,sig_Val,tsel)
        ax=subplot2grid((lengrp,1),(i,0));
        ax.locator_params(axis='y', nbins = nofticks)
        N=int(AppRes_Lst[appindx]['N_%s'%(gstr)])
        if not i: title(textwrap.dedent(TrialTitles[8]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs+6);
        my_adjust_plot1(AppRes_Lst[appindx]['MV_%s_times'%(gstr)]/second, \
            AppRes_Lst[appindx]['MV_%s'%(gstr)][N*3/4-1]/mV,fsize=ticksfs+6, \
            lwidth=lwidth,marker='',showgrid=False) # yminmax=AppInf['%s_MnMxmV'%(gstr)]);
        tick_params(labelbottom='off')
        if (i==(lengrp-1)):
            xlabel('Zaman (saniye)',fontsize=lblfs+6); tick_params(labelbottom='on')
        ylabel('%s' %(gstr),fontsize=lblfs+6);
    fig.tight_layout()
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Vm)'+trial_idq + figformat1,dpi=dpi1) # , bbox_inches='tight')

# raster plot
for x in range(len(trialSel)):
    raster_data=[];N_data=[];
    tsel=trialSel[x]; appindx = getAppindx(p_Val,sig_Val,tsel);
    trial_idq='_(%s - %s)'%(str(tsel[0]).replace('.','_'),str(tsel[1]).replace('.','_'))
    colorlst=list(np.copy(AppInf['colorlist'][::-1]));yticks = list(np.copy(AppInf['Groups'][::-1]));
    fig=figure(figsize=(20,figKat*lengrp)); # figure(figsize=(16+2,figKat*lengrp*0.75));
    if AppInf['Ext_input']:
        colorlst.extend(['r.']); yticks.extend(['EXT'])
        raster_data.append(AppRes_Lst[appindx]['M_%s_%s_spikes'%(AppInf['Groups'][0],'input')])
        N_data.append(int(AppRes_Lst[appindx]['N_%s'%('EXT')]))
    for i in range(lengrp):
        gstr = AppInf['Groups'][i]; 
        raster_data.append(AppRes_Lst[appindx]['M_%s_spikes'%(gstr)])
        N_data.append(int(AppRes_Lst[appindx]['N_%s'%(gstr)]))
    #subplot2grid((lengrp,2),(0,1),rowspan=lengrp);
    rast_title=textwrap.dedent(TrialTitles[9]%(str(tsel[0]),str(tsel[1]))).strip()
    my_raster_plot2(raster_data[::-1],N_data[::-1],colorlist=colorlst,fsize=ticksfs+10, \
        title=rast_title,xlabel='Zaman (saniye)',ylabel='Grup',spacebetweengroups=2.,yticks = yticks)
    fig.tight_layout()
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Raster)'+trial_idq + figformat1,dpi=dpi1) #, bbox_inches='tight')

# show()
matplotlib.pylab.close('all')
"""
# LFP ve PSD
# total memb.

for x in range(len(trialSel)):
    tsel=trialSel[x];
    trial_idq='_(%s - %s)'%(str(tsel[0]).replace('.','_'),str(tsel[1]).replace('.','_'))
    fig=figure(figsize=(24+3,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i]; appindx = getAppindx(p_Val,sig_Val,tsel);
        npow1=512 ; dt=double(AppRes_Lst[appindx]['dt'])*second; 
        duration=double(AppRes_Lst[appindx]['duration'])*second
        inert_time=double(AppRes_Lst[appindx]['inert_time'])*second
        npow2=AppInf['NFFT'] # int((duration-inert_time)/dt)
        totMemb = detrend(AppRes_Lst[appindx]['%s_total_membrane'%(gstr)]/mV)
        
        subplot2grid((lengrp,3),(i,0));
        if not i: title(textwrap.dedent(TrialTitles[10]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        my_adjust_plot1(AppRes_Lst[appindx]['MV_%s_times'%(gstr)]/second,totMemb, \
        fsize=ticksfs,lwidth=lwidth,marker='',showgrid=False, yminmax=AppInf['%s_MnMxtotmemb'%(gstr)]);
        tick_params(labelbottom='off')
        if (i==(lengrp-1)):
            xlabel('Zaman (saniye)',fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s\nmV'%(gstr),fontsize=lblfs);
        
        totMemb = detrend(AppRes_Lst[appindx]['%s_total_membrane'%(gstr)]/mV)
        subplot2grid((lengrp,3),(i,1));
        if not i: title(textwrap.dedent(TrialTitles[11]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        PxxTM, freqsTM = psd(totMemb, window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt); xlabel('')
        #PxxTM=AppRes_Lst[appindx]['%s_PxxTM'%(gstr)]; freqsTM=AppRes_Lst[appindx]['%s_freqsTM'%(gstr)]
        #my_adjust_plot1(freqsTM,10*log10(PxxTM),lwidth=lwidth,marker='',showgrid=False,fsize=ticksfs);
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel('Frekans (Hz)',fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('dB/Hz',fontsize=lblfs) #ylabel(r'$PSD (mV^2/Hz)$',fontsize=lblfs)
        xlim(0,100); #maxPxxTM=10*log10(max(PxxTM)); ylim(-20,maxPxxTM+maxPxxTM*0.1);

        subplot2grid((lengrp,3),(i,2));
        if not i: title(textwrap.dedent(TrialTitles[12]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        PxxTM, freqsTM, bins, im = specgram(totMemb, window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('')
        tick_params(labelbottom='off')
        ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
        if (i==(lengrp-1)): 
            xlabel('Zaman (ms)',fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('Frekans (Hz)',fontsize=lblfs) #ylabel(r'$PSD (mV^2/Hz)$',fontsize=lblfs)
        
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(totMemb_PSD)'+trial_idq + figformat1,dpi=dpi1, bbox_inches='tight')
"""
# total LFP
values=[]
for x in range(len(trialSel)):
    tsel=trialSel[x];
    trial_idq='_(%s - %s)'%(str(tsel[0]).replace('.','_'),str(tsel[1]).replace('.','_'))
    fig=figure(figsize=(30,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i]; appindx = getAppindx(p_Val,sig_Val,tsel);
        npow1=512 ; dt=double(AppRes_Lst[appindx]['dt'])*second; 
        duration=double(AppRes_Lst[appindx]['duration'])*second
        inert_time=double(AppRes_Lst[appindx]['inert_time'])*second
        npow2=AppInf['NFFT'] # int((duration-inert_time)/dt)
        totLFP = detrend(AppRes_Lst[appindx]['%s_total_LFP'%(gstr)]/mA) #(nA))
        
        ax=subplot2grid((lengrp,4-1),(i,0));
        ax.locator_params(nbins = nofticks)

        if not i: title(textwrap.dedent(TrialTitles[13]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        my_adjust_plot1(AppRes_Lst[appindx]['MV_%s_times'%(gstr)]/second,totLFP, \
        fsize=ticksfs,lwidth=lwidth,marker='',showgrid=False, yminmax=AppInf['%s_MnMxtotLFP'%(gstr)]);
        tick_params(labelbottom='off')
        if (i==(lengrp-1)):
            xlabel('Zaman (saniye)',fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s\nmV'%(gstr),fontsize=lblfs);
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
        
        ax=subplot2grid((lengrp,4-1),(i,1));
        ax.locator_params(nbins = nofticks)
        
        if not i: title(textwrap.dedent(TrialTitles[14]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        PxxTP, freqsTP = psd(totLFP, window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt,lw=lwidth); xlabel('')
        #PxxTP=AppRes_Lst[appindx]['%s_PxxTM'%(gstr)]; freqsTP=AppRes_Lst[appindx]['%s_freqsTM'%(gstr)]
        #my_adjust_plot1(freqsTP,10*log10(PxxTP),lwidth=lwidth,marker='',showgrid=False,fsize=ticksfs);
        ax.tick_params(labelbottom='off'); ax.grid(False)
        if (i==(lengrp-1)): 
            xlabel('Frekans (Hz)',fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('') #ylabel(r'$PSD (dB/Hz)$',fontsize=lblfs) # 
        xlim(0,100); #maxPxxTM=10*log10(max(PxxTP)); 
        ylim(AppInf['psdylim1'][1],AppInf['psdylim1'][0]);
        text(5,AppInf['psdylim1'][0]-20,'SR=%s, CR=%s'%(int(AppRes_Lst[appindx]['%s_snrTP1'%(gstr)]), \
                int(AppRes_Lst[appindx]['%s_CFactor_betaTP1'%(gstr)])), fontsize=textfs)
        #text(5,AppInf['psdylim1'][0]-20,'SR=%s'%(int(AppRes_Lst[appindx]['%s_snrTP1'%(gstr)])), fontsize=textfs)
        #text(5,AppInf['psdylim1'][0]-20,'CR=%s'%(int(AppRes_Lst[appindx]['%s_CFactor_betaTP1'%(gstr)])), fontsize=textfs)
        """
        peaks=peakdetOscBands3(PxxTP,freqsTP,max_peak_ind=100)
        for a in peaks.keys():
            pp=peaks['%s'%(a)][0]
            if (pp>0):
                annotate(r'$\%s:%s$'%(a,round(freqsTP[pp],2)), xy=(freqsTP[pp],10*log10(PxxTP[pp])), \
                        xytext = (20, 30),textcoords = 'offset points', fontsize=textfs,\
                        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5', color='k'))
                plot(freqsTP[pp],10*log10(PxxTP[pp]),'r.',mew=2,ms=10)
        """
        
        ax=subplot2grid((lengrp,4-1),(i,2));
        ax.locator_params(nbins = nofticks)
        if not i: title(textwrap.dedent(TrialTitles[14]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        PxxTP, freqsTP = psd(totLFP, window=hamming(npow2), noverlap=npow2/2, NFFT=npow2, Fs=1*second/dt,lw=lwidth); xlabel('')
        #PxxTP=AppRes_Lst[appindx]['%s_PxxTM'%(gstr)]; freqsTP=AppRes_Lst[appindx]['%s_freqsTM'%(gstr)]
        #my_adjust_plot1(freqsTP,10*log10(PxxTP),lwidth=lwidth,marker='',showgrid=False,fsize=ticksfs);
        tick_params(labelbottom='off'); ax.grid(False)
        if (i==(lengrp-1)):
            xlabel('Frekans (Hz)',fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('') #ylabel(r'$PSD (dB/Hz)$',fontsize=lblfs) # ylabel('dB/Hz',fontsize=lblfs) # 
        xlim(0,25); #maxPxxTM=10*log10(max(PxxTP)); 
        ylim(AppInf['psdylim2'][1],AppInf['psdylim2'][0]);
        maxindTP, maxvalTP = max(enumerate(PxxTP), key=operator.itemgetter(1))
        fr=freqsTP[maxindTP]
        ps=10*log10(maxvalTP)
        annotate('%s Hz - %s dB/Hz'%(round(fr,2),round(ps,2)),xy=(fr,ps), fontsize=textfs, \
            xycoords='data',xytext=(0.2, 0.85), textcoords='axes fraction', \
            arrowprops=dict(facecolor='black', shrink=0.05))
        """
        subplot2grid((lengrp,4),(i,3));
        if not i: title(textwrap.dedent(TrialTitles[15]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        PxxTP, freqsTP, bins, im = specgram(totLFP, window=hamming(npow1), noverlap=npow1/2, NFFT=npow1, Fs=1*second/dt); xlabel('')
        tick_params(labelbottom='off')
        ticklabel_format(style='sci', axis='y', scilimits=(-3,3))
        if (i==(lengrp-1)): 
            xlabel('Zaman (ms)',fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('Frekans (Hz)',fontsize=lblfs) #ylabel(r'$PSD (mV^2/Hz)$',fontsize=lblfs)
        """
        values.append("p=%s sig=%s R=%s B=%s"%(tsel[0],tsel[1],AppRes_Lst[appindx]['%s_meanPhase'%(gstr)],AppRes_Lst[appindx]['%s_SyncBurst'%(gstr)]))
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(totLFP_PSD)'+trial_idq + figformat1,dpi=dpi1, bbox_inches='tight')


"""
# Vm color plot full 
for x in range(len(trialSel)):
    tsel=trialSel[x]; 
    trial_idq='_(%s - %s)'%(str(tsel[0]).replace('.','_'),str(tsel[1]).replace('.','_'))
    fig=figure(figsize=(16,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i]; appindx = getAppindx(p_Val,sig_Val,tsel)
        subplot(lengrp,2,2*i+1);
        if not i: title(textwrap.dedent(TrialTitles[16]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        N=int(AppRes_Lst[appindx]['N_%s'%(gstr)])
        duration=(int(AppRes_Lst[appindx]['duration'])*second/ms)
        inert_time=double(AppRes_Lst[appindx]['inert_time'])*second/ms
        my_adjust_imshow(AppRes_Lst[appindx]['MV_%s'%(gstr)]/mV,extent=(inert_time, duration, 0, N-1))
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel('Zaman (ms)',fontsize=lblfs); tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        subplot(lengrp,2,2*i+2);
        kduration=duration/4# 5000*ms
        if not i: title(textwrap.dedent(TrialTitles[17]%(str(int(kduration/ms)),str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        my_adjust_imshow(AppRes_Lst[appindx]['MV_%s'%(gstr)]/mV,extent=(inert_time, duration, 0, N-1))
        xlim(inert_time+duration/10,inert_time+kduration);
        tick_params(labelbottom='off')
        if (i==(lengrp-1)): 
            xlabel('Zaman (ms)',fontsize=lblfs); tick_params(labelbottom='on')
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(VmColorFull)'+trial_idq + figformat1,dpi=dpi1, bbox_inches='tight')
"""
# popülasyon spike ateşleme uyumluluğu (coherence)
for x in range(len(trialSel)):
    tsel=trialSel[x]; 
    trial_idq='_(%s - %s)'%(str(tsel[0]).replace('.','_'),str(tsel[1]).replace('.','_'))
    fig=figure(figsize=(20,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i]; appindx = getAppindx(p_Val,sig_Val,tsel)
        ax=subplot(lengrp,2,2*i+1);
        CC=AppRes_Lst[appindx]['%s_SpikeCoh'%(gstr)]
        if not i: title(textwrap.dedent(TrialTitles[18]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        my_pcolor_plot(AppRes_Lst[appindx]['%s_CC_map'%(gstr)], zminmax=(0,1),fsize=ticksfs) #zminmax=AppInf['%s_MnMxCCmap'%(gstr)])
        xlabel('%s'%(gstr),fontsize=lblfs); ylabel('%s (CC=%s)'%(gstr,str(CC)),fontsize=lblfs);
        #ax.set_xticklabels(ax.get_xticks().tolist(),fontsize=ticksfs)
        
        ax=subplot(lengrp,2,2*i+2);
        if not i: title(textwrap.dedent(TrialTitles[19]%(str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        N=int(AppRes_Lst[appindx]['N_%s'%(gstr)])
        my_adjust_imshow(array(AppRes_Lst[appindx]['%s_spike_map'%(gstr)]),fsize=ticksfs,\
            extent=(basla/bin, (bitir-basla)/bin, 0, N-1),interpolation='nearest',\
            cmap=plt.cm.binary,showcolorbar=False)
        #tick_params(labelbottom='off')
        ax.locator_params(axis='x', nbins = nofticks)
        if (i==(lengrp-1)): 
            xlabel(r'Zaman $(t/{\Delta}t)$',fontsize=lblfs);
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(CC_map-spike_map)'+trial_idq + figformat1,dpi=dpi1, bbox_inches='tight')

# faz uyumu
for x in range(len(trialSel)):
    tsel=trialSel[x]; 
    trial_idq='_(%s - %s)'%(str(tsel[0]).replace('.','_'),str(tsel[1]).replace('.','_'))
    fig=figure(figsize=(10,figKat*lengrp));
    for i in range(lengrp):
        gstr = AppInf['Groups'][i]; appindx = getAppindx(p_Val,sig_Val,tsel)
        subplot(lengrp,1,i+1);
        netmphase=round(AppRes_Lst[appindx]['%s_meanPhase'%(gstr)],2)
        if not i: title(textwrap.dedent(TrialTitles[20]%(str(netmphase),str(tsel[0]),str(tsel[1]))).strip(),fontsize=titfs);
        my_pcolor_plot(array(AppRes_Lst[appindx]['%s_netphases'%(gstr)]),zminmax=(0,1),fsize=ticksfs)#AppInf['%s_MnMxnetPhase'%(gstr)])
        #tick_params(labelbottom='off')
        #if (i==(lengrp-1)): 
        xlabel('%s'%(gstr),fontsize=lblfs); #tick_params(labelbottom='on')
        ylabel('%s'%(gstr),fontsize=lblfs);
        
    if recPic : fig.savefig(AppInf['path_images']+ pref_str +'_(Phase_Coh)'+trial_idq + figformat1,dpi=dpi1, bbox_inches='tight')

# show()
matplotlib.pylab.close('all')

# tüm bandlara ait SR, CR değerleri
lenbnd = len(AppInf['bands']);
TrialbandSpects=[]
for x in range(len(AppInf['trialSeld'])):
    bandSpects = {}
    tsel=trialSel[x];
    for bndindx in range(lenbnd):
        for i in range(lengrp):
            gstr = AppInf['Groups'][i];
            bnd = AppInf['bands'][bndindx]
            appindx = getAppindx(p_Val,sig_Val,tsel);
            SR=int(AppRes_Lst[appindx]['%s_snrTP(%s)'%(gstr,bnd)])
            CR=int(AppRes_Lst[appindx]['%s_CFactor_betaTP(%s)'%(gstr,bnd)])
            
            
            bandSpects['%s'%(bnd)]=[SR, CR]
    TrialbandSpects.append(bandSpects)
