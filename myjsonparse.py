import pickle, json
import matplotlib.pyplot as plt
import math
from collections import Iterable
import numpy as np
from scipy import signal
import matplotlib.mlab as mlab


try:
    basestring
except NameError:
    basestring = str

############################# raster plot #################################
def _take_options(myopts, givenopts):
    """
    varsayılan parametre değerlerini girilen parametre değerleri ile değiştirir.
    """
    for k in myopts.keys():
        if k in givenopts:
            myopts[k] = givenopts.pop(k)

def my_raster_plot(spikes,cellid,npopCell,popColors, **additionalplotoptions):
	myopts = {"title":"", "xlabel":"Time (ms)",
			  "yticks":"", "ylabel":"Neuron id", "showplot":False,
			  "showgrouplines":False, "spacebetweengroups":0.0,
			  "grouplinecol":"k", 'newfigure':False,
			  'showlast':None, 'invert_yaxis':False}
	_take_options(myopts, additionalplotoptions)
	spacebetween = myopts['spacebetweengroups']
	if myopts['newfigure']:	figure()
	ytick=[]
	for lbl, arr in npopCell.items():
		crng=npopCell[lbl]
		cindx=[i for i, e in enumerate(spkid) if e in range(crng[0],crng[1])]
		plt.scatter(spikes[cindx],cellid[cindx],c=popColors[lbl], marker='.')
	if myopts['invert_yaxis']:
		plt.gca().invert_yaxis()
	if myopts['yticks']!="":
		#ytick=[(2*x+1)/2. for x in range(len(myopts['yticks']))]
		yticks(ytick, myopts['yticks'])
	plt.ylabel(myopts['ylabel'])
	plt.xlabel(myopts['xlabel'])
	plt.title(myopts["title"])

def pol2cart(th,r):
    '''
    kutupsal dan kartezyen koordinatlara dönüşüm.
    '''
    x = r*np.cos(th)
    y = r*np.sin(th)
    return x, y

def cart2pol(x,y):
    '''
     kartezyen den kutupsal koordinatlara dönüşüm.
    '''
    th = np.arctan2(y+0,x+0);
    r = np.hypot(x,y);
    return th, r

def compass(x,y,fsize=12,showgrid=True,color='b',lwidth=1,glwidth=1):
    if showgrid: plt.grid(color='r', linestyle='--', linewidth=glwidth, alpha=0.5)
    plt.rc("font", size=fsize)
    xx = np.array([0.0, 1.0, 0.80, 1.0, 0.8]).reshape(-1,1)
    yy = np.array([0.0, 0.0, 0.08, 0.0, -0.08]).reshape(-1,1)
    arrow = xx + 1j*yy
    x = x[:]; y = y[:];
    z = np.array(x + y*1j).reshape(-1,1);z=z.T
    a = arrow * z;
    [th,r] = cart2pol(a.real,a.imag)
    h=plt.polar(th,r,c=color,lw=lwidth);
    return h

def my_adjust_polar(X,Y,Theta,ara=1,fsize=16,showgrid=True,color='b',\
        lwidth=1,glwidth=1):
    vlen=np.round(np.sqrt(np.array(list(flatten_list(X))[::ara])**2+np.array(list(flatten_list(Y))[::ara])**2))
    [x,y] = pol2cart(np.array(list(flatten_list(Theta))[::ara]),vlen); 
    h=compass(x,y,fsize=fsize,showgrid=showgrid,color=color,lwidth=lwidth,glwidth=glwidth);
    return h

def my_adjust_hist(Vals,bin=10,fsize=24,norm=False,fcolor='g',maxHist=None,\
        showgrid=False,ecolor='k',alpha=0.7,rwidth=0.75,lwidth=2, linepoint=None, \
        xminmax=None, yminmax=None):
    
    plt.rc("font", size=fsize)
    if showgrid: plt.grid(color='r', linestyle='--', linewidth=1, alpha=0.5)
    if maxHist is None: maxHist = np.ceil(max(Vals))
    r, bins, patches=plt.hist(Vals,bin,normed=norm,align='mid',\
        color=fcolor,edgecolor=ecolor,alpha=alpha, rwidth=rwidth,linewidth=lwidth)
    
    if len(Vals):
        if xminmax==None: xminmax = (min(Vals), max(Vals)) #xmean = (xmax-xmin)/2 # mean(app)
        if yminmax==None: yminmax = (min(r), max(r)) # ymean = (ymax-ymin)/2 # mean(Vals); 
        
        plt.xlim(xminmax[0]-(np.mean(xminmax))*0.1, xminmax[1]+(np.mean(xminmax))*0.1)
        plt.ylim(yminmax[0]-(np.mean(yminmax))*0.1, yminmax[1]+(np.mean(yminmax))*0.1)
        
        if linepoint!=None:
            plt.plot([linepoint, linepoint], [yminmax[0], yminmax[1]], 'r-', lw=4, alpha=0.7)
    
    #xlim(np.min(bins)-abs(np.mean(bins))*0.2,np.max(bins)+abs(np.mean(bins))*0.2)
    #ylim(0,np.max(r)+abs(np.mean(r))*0.1)
    return r, bins, patches

def my_adjust_imshow(Vals,extent,fsize=16,cmap=plt.cm.jet,showgrid=False, \
    showcolorbar=True,glwidth=1,interpolation='bilinear',zminmax=None):
    plt.rc("font", size=fsize)
    
    if zminmax==None: zminmax = (Vals.min(), Vals.max())

    if showgrid: plt.grid(color='r', linestyle='--', linewidth=glwidth)
    im=plt.imshow(Vals,cmap=cmap, aspect='auto',origin='lower',vmin=zminmax[0], vmax=zminmax[1],extent=extent,interpolation=interpolation)
    if showcolorbar: 
        CB = plt.colorbar(im)
        CB.ax.tick_params(labelsize=fsize) 

def my_adjust_contourf(app1,app2,Vals,fsize=16,cmap=plt.cm.jet,alpha=1,drwline=False,norm=None,levels=None):
    plt.rc("font", size=fsize)
    origin = 'lower'
    X, Y = np.meshgrid(app1, app2);
    Z=Vals;
    
    if levels ==None:
        CS = plt.contourf(X, Y, Z, alpha=alpha, cmap=cmap,origin=origin,norm=norm)
    else:
        CS = plt.contourf(X, Y, Z, levels, alpha=alpha, cmap=cmap,origin=origin,norm=norm,extend='both');
        if (levels.max()<Z.max() and levels.min()>Z.min()): CS.extend='both'
        elif levels.max()<Z.max(): CS.extend='max'
        elif levels.min()>Z.min(): CS.extend='min'
        else : CS.extend='neither'
        CS.cmap.set_under('navy')
        CS.cmap.set_over('darkred')

    #CS.set_clim(vmin=0*mn-0*mn*0.5,vmax=mx)
    if (drwline and Vals.min()!=Vals.max()): 
        plt.contour(X, Y, Z, colors='k',origin=origin)
    cbar=plt.colorbar(CS)#,format='%.3g');
    #cbar = mpl.colorbar.ColorbarBase(CS,extend='both')
    
    #cbar.set_ticks(linspace(mn,mx,5))#([mn, md, mx])
    #cbar.set_ticklabels(linspace(mn,mx,5))#([mn, md, mx])
    
    #cbar.formatter.set_useOffset(False)
    cbar.formatter.set_scientific(True)
    cbar.formatter.set_powerlimits((-3,3))
    cbar.update_ticks() 

def vector_strength2(spikes, freq):
    '''
    spike dizisinin vector strength değerini hesaplar.
    '''
    VS=np.nan; theta=np.nan; X=[]; Y=[]; thetas=[];
    
    if spikes.size: #len(np.array(spikes,ndmin=1)):  
        # return abs(mean(exp(np.array(spikes) * 1j * 2 * pi / period)))
        # N = len(spikes)
        X = (np.cos(spikes * 2 * np.pi * freq))
        Y = (np.sin(spikes * 2 * np.pi * freq))
        thetas=np.mod(spikes,1.0/freq)*(2*np.pi * freq)
        #thetas=(spikes * 2 * np.pi) / period;
        theta = np.mean(thetas)
        
        VS = np.sqrt(np.mean(X)**2+np.mean(Y)**2)
        #VS = abs(mean(np.exp(spikes * 1j * 2 * np.pi *freq)))
        #VS = abs(sum(np.exp(spikes * 1j * 2 * np.pi *freq)))
    return VS, theta, X, Y, thetas
    

def population_VS_Calc2(spikes, freq):
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
        [VS,Th,X,Y,thetas] = vector_strength2(np.array(spikesvec),freq)
    else:
        VS=0
    
    return VS,Th,thetas,X,Y
        
def population_VS_Calc(spikes, freq, bins,N=None):
    '''
    populasyonun VS ve ilgili tüm parametrelerini hesaplar.
    refs: Ashida, G., ve Carr, C. E. (2009). Effect of sampling frequency on the measurement of phase-locked action potentials
    Carr, C. E., ve Friedman, M. A. (1999). Evolution of time coding systems
    
    0<=VS<=1; 0 asenkron, 1 senkron
    '''
    if N==None: N=len(spikes)
    
    VSs = np.zeros(N); VSsTh = [np.NaN for x in range(N)]; 
    VSsTheta = []; 
    VSsX = [[] for x in range(N)]; VSsY = [[] for x in range(N)]; 
    VSsThetaM =  [[] for x in range(N)];
    
    if not(len(spikes) and freq):
        return VSs,VSsTh,VSsTheta,VSsX,VSsY,np.squeeze(VSsThetaM),0
    
    #period=1.0/freq
    for i in range(N):
        [VS,Th,X,Y,thetas] = vector_strength2(np.array(spikes[i]),freq)
        if np.isnan(VS): VS=0
        VSs[i]=(VS); VSsTh[i]=(Th); VSsTheta.append(thetas)
        VSsX[i].append(X); VSsY[i].append(Y);
        a=np.histogram(thetas,bins,normed=False)
        VSsThetaM[i].append(a[0])
    # tüm ağın VS değeri ağın tüm spikeları birleştirilerek hesaplanıyor
    meanVS = population_VS_Calc2(spikes,freq)[0]
    #meanVS = mean(VSs) #sum([ x for x in VSs if not math.isnan(x)])/len(spikes)
    
    # --->> sil ! # zayıf spike aktivitesi yanıltıcı olduğundan VS frekans ile carpıldı
    #netmf = network_mean_frequency(spikes,N)[0] 
    
    return VSs,VSsTh,VSsTheta,VSsX,VSsY,np.squeeze(VSsThetaM),meanVS

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


def butter_highpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y


fname = 'v53_batch12_0_0_0' # 'v53_batch12_1_0_0' # 'v53_batch12_2_0_0' # 'M1detailed' # 

with open(fname+'.json', 'r') as f:
    M1detailed_dict = json.load(f)
    print(fname+' file is loaded')

#for M1detailed in M1detailed_dict
#print(M1detailed_dict["simData"]["spkt"])
#print(M1detailed_dict["simData"].keys())

spkt=M1detailed_dict["simData"]["spkt"]
spkid=M1detailed_dict["simData"]["spkid"]
spikesTimes=np.array(spkt)
cellid=np.array(spkid)
duration=M1detailed_dict["simConfig"]["duration"] # in ms
rec_dt=M1detailed_dict["simConfig"]["recordStep"] # in ms
Fs=1000.0/rec_dt
modelLFPfreqs=None

avgRate=M1detailed_dict["simData"]["avgRate"]
popRates=M1detailed_dict["simData"]["popRates"]

figformat1='.png';dpi1=150;
alpha = 0.9; cmap=plt.cm.jet;

if 'LFP' in M1detailed_dict["simData"].keys():
    LFP=M1detailed_dict["simData"]["LFP"]
    print('LFP exist: OK')
    modelLFPfreqs=[]
    nsmooth=50
    fig=plt.figure()
    for i in range(len(LFP[0])):
        LFP1=np.array(LFP)[...,i]
        LFP1 = butter_highpass_filter(LFP1,3,Fs)
        win=np.hanning(len(LFP1))
        pxx, freq=mlab.psd(((LFP1-LFP1.mean())*win)/len(LFP1), NFFT=len(LFP1), Fs=Fs, pad_to=len(LFP1),window=mlab.window_none)
        psdVal=10*np.log10(pxx)
        weights = np.ones(nsmooth, dtype=float) / nsmooth
        psd_s = np.convolve(10*np.log10(pxx), weights, mode='valid')
        freqs_s = np.convolve(freq, weights, mode='valid')
        
        plt.plot(freqs_s,psd_s)
        m=max(psdVal)
        maxindxs=[i for i, j in enumerate(psdVal) if j == m]
        modelLFPfreqs.append(freq[maxindxs[-1]])
    plt.xlim(0,100)
    #plt.show()
    fig.savefig('LFP_PSDs_' + fname + figformat1,dpi=dpi1, bbox_inches='tight')

"""
fig=plt.figure(10,8)
ax=plt.subplot(3,1,1)
plot(LFP1)

LFP1 = butter_highpass_filter(LFP1,3,Fs)
ax=plt.subplot(3,1,2)
plot(LFP1)

ax=plt.subplot(3,1,3)
win=np.hanning(len(LFP1))
pxx, freq=plt.psd(((LFP1-LFP1.mean())*win)/len(LFP1), NFFT=len(LFP1), Fs=Fs, pad_to=len(LFP1),window=mlab.window_none)
plt.xlim(0,200)
plt.show()

plt.plot(freq,10*log10(pxx))
plt.xlim(0,200)
plt.show()
"""

with open('../cells/popColors.pkl', 'rb') as fileObj: popColors = pickle.load(fileObj)['popColors']

layer = {'2': [0.12,0.31], '4': [0.31,0.42], '5A': [0.42,0.52], '45A':[0.31,0.52], '5B': [0.52,0.77], '6': [0.77,1.0], 'long': [2.0,3.0]}  # normalized layer boundaries

"""
from netpyne import sim
sim.load('v53_batch12_1_0_0.json', instantiate=False)
sim.loadAll('v53_batch12_1_0_0.json')
sim.loadNetParams('v53_batch12_1_0_0.json')
sim.loadNet('v53_batch12_1_0_0.json')
sim.net.createPops()   
sim.net.createCells()  

sim.net.pops['SOM6'].cellGids

npopCell2={}
for lbl in sim.net.pops.keys():
    npopCell2[lbl]=[np.min(sim.net.pops[lbl].cellGids),np.max(sim.net.pops[lbl].cellGids)]

print(npopCell2)

{'IT2': [0, 1729], 'SOM2': [1730, 1856], 'PV2': [1857, 2113], 'IT4': [2114, 2880], 'IT5A': [2881, 3524], 
 'SOM5A': [3525, 3580], 'PV5A': [3581, 3694], 'IT5B': [3695, 5129], 'PT5B': [5130, 6564], 'SOM5B': [6565, 6815], 
 'PV5B': [6816, 7326], 'IT6': [7327, 8562], 'CT6': [8563, 9798], 'SOM6': [9799, 9888], 'PV6': [9889, 10072], 'TPO': [10073, 10073]}

{'IT2': [0, 1730], 'SOM2': [1730, 1857], 'PV2': [1857, 2114], 'IT4': [2114, 2881], 'IT5A': [2881, 3525], 
 'SOM5A': [3525, 3581], 'PV5A': [3581, 3695], 'IT5B': [3695, 5130], 'PT5B': [5130, 6565], 'SOM5B': [6565, 6816], 
 'PV5B': [6816, 7327], 'IT6': [7327, 8563], 'CT6': [8563, 9799], 'SOM6': [9799, 9889], 'PV6': [9889, 10073], 'TPO': [10073, 11073]}
"""

nIT2=0+M1detailed_dict["net"]["params"]["popParams"]["IT2"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*density[('M1','E')][0]*(layer['2'][1]-layer['2'][0]))
nSOM2=nIT2+M1detailed_dict["net"]["params"]["popParams"]["SOM2"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*density[('M1','SOM')][0]*(layer['2'][1]-layer['2'][0]))
nPV2=nSOM2+M1detailed_dict["net"]["params"]["popParams"]["PV2"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*density[('M1','PV')][0]*(layer['2'][1]-layer['2'][0]))
nIT4=nPV2+M1detailed_dict["net"]["params"]["popParams"]["IT4"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*density[('M1','E')][1]*(layer['4'][1]-layer['4'][0]))
nIT5A=nIT4+M1detailed_dict["net"]["params"]["popParams"]["IT5A"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*density[('M1','E')][2]*(layer['5A'][1]-layer['5A'][0]))
nSOM5A=nIT5A+M1detailed_dict["net"]["params"]["popParams"]["SOM5A"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*sum(density[('M1','SOM')][1:3])/2*(layer['45A'][1]-layer['45A'][0]))
nPV5A=nSOM5A+M1detailed_dict["net"]["params"]["popParams"]["PV5A"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*sum(density[('M1','PV')][1:3])/2*(layer['45A'][1]-layer['45A'][0]))
nIT5B=nPV5A+M1detailed_dict["net"]["params"]["popParams"]["IT5B"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*0.5*density[('M1','E')][3]*(layer['5B'][1]-layer['5B'][0]))
nPT5B=nIT5B+M1detailed_dict["net"]["params"]["popParams"]["PT5B"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*0.5*density[('M1','E')][3]*(layer['5B'][1]-layer['5B'][0]))
nSOM5B=nPT5B+M1detailed_dict["net"]["params"]["popParams"]["SOM5B"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*density[('M1','SOM')][3]*(layer['5B'][1]-layer['5B'][0]))
nPV5B=nSOM5B+M1detailed_dict["net"]["params"]["popParams"]["PV5B"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*density[('M1','PV')][3]*(layer['5B'][1]-layer['5B'][0]))
nIT6=nPV5B+M1detailed_dict["net"]["params"]["popParams"]["IT6"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*0.5*density[('M1','E')][4]*(layer['6'][1]-layer['6'][0]))
nCT6=nIT6+M1detailed_dict["net"]["params"]["popParams"]["CT6"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*0.5*density[('M1','E')][4]*(layer['6'][1]-layer['6'][0]))
nSOM6=nCT6+M1detailed_dict["net"]["params"]["popParams"]["SOM6"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*density[('M1','SOM')][4]*(layer['6'][1]-layer['6'][0]))
nPV6=nSOM6+M1detailed_dict["net"]["params"]["popParams"]["PV6"]["numCells"] # math.floor((22/7)*1.35*0.2*0.2*density[('M1','PV')][4]*(layer['6'][1]-layer['6'][0]))

nTPO=nPV6+M1detailed_dict["net"]["params"]["popParams"]["TPO"]["numCells"] 

npopCell={'IT2': [0, nIT2], 
		  'SOM2': [nIT2, nSOM2], 
		  'PV2': [nSOM2, nPV2], 
          'IT4': [nPV2, nIT4],
		  'IT5A': [nIT4, nIT5A], 
		  'SOM5A': [nIT5A, nSOM5A], 
		  'PV5A': [nSOM5A, nPV5A], 
		  'IT5B': [nPV5A, nIT5B], 
		  'PT5B': [nIT5B, nPT5B], 
		  'SOM5B': [nPT5B, nSOM5B], 
		  'PV5B': [nSOM5B, nPV5B], 
		  'IT6': [nPV5B, nIT6], 
		  'CT6': [nIT6, nCT6], 
		  'SOM6': [nCT6, nSOM6], 
		  'PV6': [nSOM6, nPV6]
         }

# Vector Strength
bins=np.linspace(0,2*np.pi,13)

figKat=len(npopCell); nofticks=5;lengrp=4;
titfs=12; lblfs=12; ticksfs=12; lwidth=2;

angle_ticks = np.linspace(0,2*np.pi,13)
angle_ticks_lbl = [r"$0$", r"$\frac{\pi}{6}$", r"$\frac{\pi}{3}$", r"$\frac{\pi}{2}$", r"$\frac{2\pi}{3}$", r"$\frac{5\pi}{6}$", r"$\pi$", r"$\frac{7\pi}{6}$", \
                   r"$\frac{4\pi}{3}$", r"$\frac{3\pi}{2}$", r"$\frac{5\pi}{3}$", r"$\frac{11\pi}{6}$", r"$2\pi$"]
angle_ticks_lbl_polar = [r"$0$", r"$\frac{\pi}{4}$", r"$\frac{\pi}{2}$", r"$\frac{3\pi}{4}$", r"$\pi$", r"$\frac{5\pi}{4}$", r"$\frac{3\pi}{2}$",  r"$\frac{7\pi}{4}$"]

popVSDict={};popmeanVSDict={}
freqStart=1;freqEnd=50;freqStep=1;
cidx=0;fig=plt.figure(figsize=(lengrp*6,figKat*3)); # each subplot size: width=6 inch and height=3 inch
for lbl, arr in npopCell.items():
    print('iter %d: %s'%(cidx,lbl))
    if modelLFPfreqs!=None: 
        cfrequency=max(set(modelLFPfreqs), key = modelLFPfreqs.count) # most frequent element in list
    else:
        cfrequency=popRates[lbl]
    cfrequency=10
    crng=npopCell[lbl]
    cindx=[i for i, e in enumerate(spkid) if e in range(crng[0],crng[1])]
    cellNum=crng[1]-crng[0]
    a=spikesTimes[cindx]
    b=cellid[cindx]
    spkSeries = [[] for _ in range(cellNum)]
    for i,idx in enumerate(b): 
        spkSeries[int(idx)-crng[0]].append(a[i])
    spkSeries2=[spkSeries[i] for i in range(cellNum) if len(spkSeries[i])]
    
    popVSDict[lbl]=list();popmeanVSDict[lbl]=list();
    for i in np.arange(freqStart,freqEnd,freqStep):
        cVS,cVSTh,cVSTheta,cVSX,cVSY,cVSThetaM,cmVS=population_VS_Calc(spkSeries2,i,bins)
        popVSDict[lbl].append(cVS)
        popmeanVSDict[lbl].append(np.mean(cVS))
    
    cVS,cVSTh,cVSTheta,cVSX,cVSY,cVSThetaM,cmVS=population_VS_Calc(spkSeries2,cfrequency,bins)
    if(cVSThetaM.ndim==1): 
        cVSThetaM=np.reshape(cVSThetaM,(1,len(cVSThetaM)))
    
    ax=plt.subplot(figKat,lengrp,lengrp*cidx+1,polar=True, projection='polar'); 
    ara=int(len(cindx)/100)+1; # 1:100 atlama
    X=cVSX;
    Y=cVSY;
    Theta=cVSTheta;
    nSpikes=len(list(flatten_list(X)))
    rnSpikes=len(np.array(list(flatten_list(X))[::ara]))
    h=my_adjust_polar(X,Y,Theta,ara=ara,fsize=ticksfs)
    ax.set_xticks(np.arange(0,2*np.pi,np.pi/6)); 
    #ax.set_xticklabels(angle_ticks_lbl_polar, fontsize=ticksfs+10)
    ax.tick_params(axis='y', labelsize=ticksfs)
    ax.tick_params(axis='x', labelsize=ticksfs)
    plt.text(np.pi, 0.275+bool(nSpikes)*(2.925), '%s\n# of cell:%d\n# of Active cell:%d\n#Spikes:%d of %d\nf:%.2fHz\nflatten mean VS:%.2f\npopulation mean VS:%.2f'\
         %(lbl,cellNum,len(spkSeries2),rnSpikes,nSpikes,cfrequency,cmVS,cVS.mean()), fontsize=lblfs)
    
    ax=plt.subplot(figKat,lengrp,lengrp*cidx+2,polar=False); 
    ax.locator_params(axis='y', nbins = nofticks)
    yminmax=(None,None)
    if len(list(flatten_list(list(sum(cVSThetaM))))):
        yminmax=(sum(cVSThetaM).min(),sum(cVSThetaM).max())
    my_adjust_hist(np.array(list(flatten_list(Theta))[::]),bin=bins,fsize=ticksfs,maxHist = 2*np.pi,yminmax=yminmax)
    ax.set_xticks(angle_ticks); 
    ax.set_xticklabels(angle_ticks_lbl,fontsize=ticksfs+10)
    #ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax.tick_params(axis='y', labelsize=ticksfs)
    ax.tick_params(axis='x', labelsize=ticksfs)
    #plt.xlabel(r'$\theta$',fontsize=lblfs)
    ax.tick_params(axis='x', labelbottom='on')
    plt.ylabel('#of spikes',fontsize=lblfs)
    if (cidx==(figKat-1)):
        plt.xlabel(r'$\theta$',fontsize=lblfs)

    ax=plt.subplot(figKat,lengrp,lengrp*cidx+3);
    N=len(cVSThetaM)
    zminmax=(None,None)
    if len(list(flatten_list(list(cVSThetaM)))):
        zminmax=(cVSThetaM.min(),cVSThetaM.max())
    my_adjust_imshow(cVSThetaM,fsize=ticksfs,extent=(0, 2*np.pi, 0, N-1),zminmax=zminmax)
    ax.set_xticks(angle_ticks); 
    ax.set_xticklabels(angle_ticks_lbl,fontsize=ticksfs+10)
    ax.tick_params(axis='y', labelsize=ticksfs)
    ax.tick_params(axis='x', labelsize=ticksfs)
    ax.tick_params(axis='x', labelbottom='on')
    plt.ylabel('cell id',fontsize=lblfs)
    if (cidx==(figKat-1)):
        plt.xlabel(r'$\theta$',fontsize=lblfs)
        
    ax=plt.subplot(figKat,lengrp,lengrp*cidx+4);
    ax.locator_params(axis='y', nbins = nofticks)
    VS_vals = cVS
    r=my_adjust_hist(VS_vals,bin=10,fsize=ticksfs,xminmax=(0,1),yminmax=(0,len(VS_vals)*1/2),linepoint=np.mean(VS_vals))
    ax.tick_params(axis='y', labelsize=ticksfs)
    ax.tick_params(axis='x', labelsize=ticksfs)
    plt.ylabel('frequency(# of val)',fontsize=lblfs)
    if (cidx==(figKat-1)):
        plt.xlabel('VS',fontsize=lblfs); 
    
    cidx=cidx+1
fig.savefig('VS_Analysis_%dHz_'%(cfrequency) + fname + figformat1,dpi=dpi1, bbox_inches='tight')

figRow=np.ceil(np.sqrt(len(popVSDict.keys())))
figCol=np.ceil(np.sqrt(len(popVSDict.keys())))
cidx=0;fig=plt.figure(figsize=(figCol*5,figRow*3)); # each subplot size: width=5 inch and height=3 inch
for lbl in popVSDict.keys():
    cidx=cidx+1
    f2VSval=np.array(popVSDict[lbl]).T
    N=len(f2VSval)
    
    ax=plt.subplot(figRow,figCol,cidx)
    zminmax=(None,None)
    if len(list(flatten_list(list(f2VSval)))):
        zminmax=(f2VSval.min(),f2VSval.max())
    
    my_adjust_imshow(f2VSval,fsize=ticksfs,extent=(freqStart, freqEnd, 0, N-1),zminmax=zminmax)
    
    plt.xlabel('frequency (Hz)',fontsize=lblfs)
    if cidx%figRow==1: plt.ylabel('cell id',fontsize=lblfs)
    plt.title(lbl,fontsize=titfs)
    
fig.savefig('VS_Analysis_Freq_%d-%d-%d_'%(freqStart,freqEnd,freqStep) + fname + figformat1,dpi=dpi1, bbox_inches='tight')


fig=plt.figure(figsize=(4,3)); # each subplot size: width=4 inch and height=3 inch
ax=plt.subplot(1,1,1)
popKeyList=list(popmeanVSDict.keys())
f2PopmeanVSval= np.array([popmeanVSDict[lbl] for lbl in popmeanVSDict.keys()])
N=len(f2PopmeanVSval)
zminmax=(None,None)
if len(list(flatten_list(list(f2PopmeanVSval)))):
    zminmax=(f2PopmeanVSval.min(),f2PopmeanVSval.max())

my_adjust_imshow(f2PopmeanVSval,fsize=ticksfs,extent=(freqStart, freqEnd, 0, N-1),zminmax=zminmax)

ax.set_yticks(np.arange(0,len(popKeyList),1)); 
ax.set_yticklabels(popKeyList,fontsize=ticksfs)

plt.xlabel('frequency (Hz)',fontsize=lblfs)
plt.ylabel('Population Name',fontsize=lblfs)
plt.title('mean populations VSs',fontsize=titfs)

fig.savefig('VS_Analysis_meanPopVS_Freq_%d-%d-%d_'%(freqStart,freqEnd,freqStep) + fname + figformat1,dpi=dpi1, bbox_inches='tight')


"""
fig=plt.figure(figsize=(24,3));
lbl= 'IT2' # 'PV5A' # 'PV6' # 'SOM6' # 'CT6' # 'IT6' # 'PV5B' # 'SOM5B' # 'PT5B' # 'IT5B' # 'SOM5A' # 'IT5A' # 'PV2' # 'SOM2' # 
if modelLFPfreqs!=None:
    cfrequency=max(set(modelLFPfreqs), key = modelLFPfreqs.count) # most frequent element in list
else:
    cfrequency=popRates[lbl]

crng=npopCell[lbl]
cindx=[i for i, e in enumerate(spkid) if e in range(crng[0],crng[1])]
cellNum=crng[1]-crng[0]
a=spikesTimes[cindx]
b=cellid[cindx]
spkSeries = [[] for _ in range(cellNum)]
for i,idx in enumerate(b):
	#print(i,idx-crng[0])
	spkSeries[int(idx)-crng[0]].append(a[i])
spkSeries2=[spkSeries[i] for i in range(cellNum) if len(spkSeries[i])]
cVS,cVSTh,cVSTheta,cVSX,cVSY,cVSThetaM,cmVS=population_VS_Calc(spkSeries2,cfrequency,bins)
if(cVSThetaM.ndim==1): 
    cVSThetaM=np.reshape(cVSThetaM,(1,len(cVSThetaM)))

ax=plt.subplot(figKat*0+1,lengrp,lengrp*cidx+1,polar=True, projection='polar'); 
ara=int(len(cindx)/100)+1; # 1:100 atlama
X=cVSX;
Y=cVSY;
Theta=cVSTheta;
nSpikes=len(list(flatten_list(X)))
rnSpikes=len(np.array(list(flatten_list(X))[::ara]))
h=my_adjust_polar(X,Y,Theta,ara=ara,fsize=ticksfs)
ax.set_xticks(np.arange(0,2*np.pi,np.pi/6)); 
ax.set_xticklabels(angle_ticks_lbl_polar, fontsize=ticksfs+10)
ax.tick_params(axis='y', labelsize=ticksfs)
ax.tick_params(axis='x', labelsize=ticksfs)
plt.text(np.pi, 0.275+bool(nSpikes)*(2.925), '%s\n# of cell:%d\n# of Active cell:%d\n#Spikes:%d of %d\nf:%.2fHz\nflatten mean VS:%.2f\npopulation mean VS:%.2f'\
         %(lbl,cellNum,len(spkSeries2),rnSpikes,nSpikes,cfrequency,cmVS,cVS.mean()), fontsize=lblfs)

ax=plt.subplot(figKat*0+1,lengrp,lengrp*cidx+2,polar=False); 
ax.locator_params(axis='y', nbins = nofticks)
yminmax=(None,None)
if len(list(flatten_list(list(sum(cVSThetaM))))):
    yminmax=(sum(cVSThetaM).min(),sum(cVSThetaM).max())
my_adjust_hist(np.array(list(flatten_list(Theta))[::]),bin=bins,fsize=ticksfs,maxHist = 2*np.pi,yminmax=yminmax)
ax.set_xticks(angle_ticks); 
ax.set_xticklabels(angle_ticks_lbl,fontsize=ticksfs+10)
#ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
ax.tick_params(axis='y', labelsize=ticksfs)
ax.tick_params(axis='x', labelsize=ticksfs)
ax.tick_params(axis='x', labelbottom='off')
#plt.xlabel(r'$\theta$',fontsize=lblfs)
ax.tick_params(axis='x', labelbottom='on')
plt.ylabel('#of spikes',fontsize=lblfs)

ax=plt.subplot(figKat*0+1,lengrp,lengrp*cidx+3,polar=False);
N=len(cVSThetaM)
zminmax=(None,None)
if len(list(flatten_list(list(cVSThetaM)))):
    zminmax=(cVSThetaM.min(),cVSThetaM.max())
my_adjust_imshow(cVSThetaM,fsize=ticksfs,extent=(0, 2*np.pi, 0, N-1),zminmax=zminmax)
ax.set_xticks(angle_ticks); 
ax.set_xticklabels(angle_ticks_lbl,fontsize=ticksfs+10)
ax.tick_params(axis='y', labelsize=ticksfs)
ax.tick_params(axis='x', labelsize=ticksfs)
#ax.tick_params(axis='x', labelbottom='off')
#plt.xlabel(r'$\theta$',fontsize=lblfs);
#ax.tick_params(axis='x', labelbottom='on')
plt.ylabel('cell id',fontsize=lblfs)

ax=plt.subplot(figKat*0+1,lengrp,lengrp*cidx+4,polar=False);
ax.locator_params(axis='y', nbins = nofticks)
VS_vals = cVS
r=my_adjust_hist(VS_vals,bin=10,fsize=ticksfs,xminmax=(0,1),yminmax=(0,len(VS_vals)*1/2),linepoint=np.mean(VS_vals))
ax.tick_params(axis='y', labelsize=ticksfs)
ax.tick_params(axis='x', labelsize=ticksfs)
plt.ylabel('cell id',fontsize=lblfs)
plt.show()
"""

"""
# raster plot
for lbl, arr in npopCell.items():
	#print(lbl,arr[0],popColors[lbl])
	crng=npopCell[lbl]
	cindx=[i for i, e in enumerate(spkid) if e in range(crng[0],crng[1])]
	plt.scatter(spikesTimes[cindx],cellid[cindx],c=popColors[lbl], marker='.')
plt.gca().invert_yaxis()
plt.ylabel('ylabel')
plt.xlabel('xlabel')
plt.title('title')
plt.show()


fig=plt.figure()
my_raster_plot(spikesTimes,cellid,npopCell,popColors,invert_yaxis=True)
plt.show()

"""


