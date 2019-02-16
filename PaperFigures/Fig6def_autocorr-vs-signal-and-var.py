# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 15:21:43 2016

@author: tweetybaldwin
"""

#Create a synthetic AR1 time series.
#Length (365 days)*(100 years) = 36500
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from copy import copy
#%%
# stdv/mean warming = 0:2.5, autocorrelation = 0:1 
size = 36500#36500 #make 36500 in the end

def syntemp(scale,phi):
    e = np.random.normal(loc = 0, scale = scale, size = size) #innovation (random gaussian variable)
    #phi autoregressive parameter (equivalent to r1 autocorrelation)
    x = np.zeros(size) #synthetic temp variable to be populated
    mu = 0 #mean of synthetic variable
    x[0] = mu + e[0]
    for t in range(1,size):
        x[t] = mu + phi*(x[t-1]-mu) + e[t]
    return x
                
#%%
#Create collection of these time series
steps = 30
scales = np.linspace(0.01,1.6,steps)
phis = np.linspace(0,1,steps)
c = np.zeros((size,steps,steps))
for i in range(0,steps):
    for j in range(0,steps):
       c[:,i,j] = syntemp(scales[i],phis[j])

#Shift mean to create global warming AR1 time series w
warming = 1
w = c + warming

#Calculate 90th percentile of control c
thres = np.percentile(c,90,axis = 0)

#Calculate exceedances (when c/w is greater than threshold).
c_ex = np.zeros((size,steps,steps))
w_ex = np.zeros((size,steps,steps))

for i in range(0,size):
    c_ex[i,...] = c[i,...]>thres
    w_ex[i,...] = w[i,...]>thres
c_ex2 = np.zeros((size,steps,steps))
w_ex2 = np.zeros((size,steps,steps))
c_ex2[c_ex>0]=c[c_ex>0]
w_ex2[w_ex>0]=w[w_ex>0]
    
#%%
#ehfs=c_ex2
###########DEFINTION-el3-bmax1-e2l1-dm4####################
# Inputs
bmax = 2 #maximum length of a break
e2length = 1 #required length of heat waves that can compound
elength = 6 #required length of first heat wave(>=e2length)
dmin = 7 #required minimum number of hot days for compound event

def identify_hw(ehfs):
    """identify_hw locates heatwaves from EHF and returns duration indicators hw, chw, ahw.
    """
     # Agregate consecutive days with EHF>0, and consecutive days where EHF<0 (not heat wave)
    events = (ehfs!=0.).astype(np.int)
    breaks = (events-1)/(-1) 
    evente = events.copy()
    breake = breaks.copy()
    #event and break durations on first day of each:
    for i in xrange(ehfs.shape[0]-2,-1,-1):
         events[i,events[i,...]>0] = events[i+1,events[i,...]>0]+1 
         breaks[i,breaks[i,...]>0] = breaks[i+1,breaks[i,...]>0]+1
    #event and break durations on last day of each:
    for i in xrange(1,ehfs.shape[0],1):
        #evente[i,evente[i,...]>0] = evente[i-1,evente[i,...]>0]+1 #event duration on last day of event
        breake[i,breake[i,...]>0] = breake[i-1,breake[i,...]>0]+1 #break duration on last day of break
    #Identify when heatwaves and breaks start and end with duration
    #Note that doesn't work for events that start on the first day (maybe should correct later)
    diffes = np.zeros(ehfs.shape)
    diffbs = np.zeros(ehfs.shape)
    diffbe = np.zeros(ehfs.shape)
    diffes[1:,...] = np.diff(events, axis=0).astype(np.int)
    diffbs[1:,...] = np.diff(breaks, axis=0).astype(np.int)
    diffbe[0:ehfs.shape[0]-1,...] = np.diff(breake, axis=0).astype(np.int)
    #Edit ends so breaks/events can end on last day, and start on first day.
    diffes[0,...]=events[0,...]
    diffbs[0,...]=breaks[0,...]
    diffbe[ehfs.shape[0]-1,...]=-1*breake[ehfs.shape[0]-1,...]
    endes = np.zeros(ehfs.shape,dtype=np.int)
    endbs = np.zeros(ehfs.shape,dtype=np.int)
    endbe = np.zeros(ehfs.shape,dtype=np.int)
    endes[diffes>0] = events[diffes>0]
    endbs[diffbs>0] = breaks[diffbs>0]
    endbe[diffbe<0] = breake[diffbe<0]
    #del breake, breaks, diffes, diffbs, diffbe
    ##################
    #HW: Find non-compound heat waves of minimum length elength
    hw = np.zeros(ehfs.shape).astype(np.int)
    hw[endes>=elength] = endes[endes>=elength]
    ##################
    #CHW: Find compound heat waves with first event >= elength, compounded events >= e2length, and total number of hot days >= dmin
    #AHW: Find compound and non-compound heat waves
    #Place in array break durations before and after events at same elements as event starts
    lasttime=ehfs.shape[0]-1
    breakbefore = np.roll(endbe,1, axis=0)
    breakafter=np.zeros(ehfs.shape).astype(np.int)
    for i in xrange(0,lasttime,1):
        xy=endes[i,...]>0 #only evaluate for first days of events
        eduration=endes[i,xy]
        bsi=i+eduration
        bsi[bsi>lasttime]=0 #first day chw always=0, so use this index
        breakafter[i,xy]=endbs[bsi,xy]
    #Populate CHW output array with only events that might compound (break before or after <=bmax, length >= e2length), and AHW output array with only events >= e2length
    chw = np.zeros(ehfs.shape).astype(np.int)
    bora = np.logical_or(breakafter<=bmax,breakbefore<=bmax)
    e2 = endes>=e2length
    boraande2 = np.logical_and(bora,e2)
    chw[boraande2] = endes[boraande2]
    #del bora, boraande2, e2
    #Add events with breaks<=bmax to ones following in time, iterated backwards in time so accumulate
    #chw counts only event (hot) days, while endlength counts days that are hot and days of breaks
    chw=np.append(chw,np.zeros((1,steps,steps)),axis=0) #Append additional day with value 0 to help with next loop.
    endlength=chw.copy()
    for i in xrange(lasttime,-1,-1):
        xy=np.logical_and(chw[i,...]>0,breakafter[i,...]<=bmax) #only evaluation for first days of events part of compound events, where the break after is less than max break length
        eduration=endes[i,xy]    
        bafter=breakafter[i,xy]
        eaddi=i+eduration+bafter
        eaddi[eaddi>lasttime]=size #appended last day always=0, so use this index
        chw[i,xy]=eduration+chw[eaddi,xy] #no. of event days in compound heat wave
        endlength[i,xy]=eduration+bafter+endlength[eaddi,xy] #no. of days (breaks and events) in compound heat wave
    chw=np.delete(chw,size,0) #Remove additional day with value 0.
    endlength=np.delete(endlength,size,0) #Remove additional day with value 0.
    chw[endes<elength]=0 #Turn compound events that can't be starts (original event length < elength) to zero
    chw[chw==endes] = 0 #Gets rid of heat waves at end that aren't compound
    chw[chw<dmin] = 0 #Get rid of heat waves that don't meet a minimum duration requirement
    #combine hw and chw to get preliminary ahw
    ahw=chw+hw
    chwhw = np.logical_and(chw>0,hw>0)
    ahw[chwhw]=chw[chwhw] #ensures don't double count hw that are summed into chw already
    del chwhw
    #Turn compound events that are not at beginning of event to zero for chw and ahw
    for i in xrange(0,lasttime,1):
        for index, j in np.ndenumerate(chw[i,...]):
                if chw[i,index[0],index[1]]>0: #add ,index[1] after  index[0] if 3-D array
                    length=int(endlength[i,index[0],index[1]]) #add ,index[1] after  index[0] if 3-D array
                    elasti=i+length-1
                    chw[(i+1):elasti,index[0],index[1]]=0 #add ,index[1] after  index[0] if 3-D array
                    ahw[(i+1):elasti,index[0],index[1]]=0
    #secondary ahw (events that compound on to other events)
    ahw2 = chw-hw
    ahw2[ahw2<0] = 0
    #primary ahw (compound starting events, and events with no compounding)
    ahw1 = ahw-ahw2
    
    return hw, chw, ahw, ahw1, ahw2
#%%    
c_hw, c_chw, c_ahw, c_ahw1, c_ahw2 = identify_hw(c_ex2)
w_hw, w_chw, w_ahw, w_ahw1, w_ahw2 = identify_hw(w_ex2)

#%%
#Calculate total number of exceedances
c_tot = np.sum(c_ex,axis=0)
w_tot = np.sum(w_ex,axis=0)
d_tot = w_tot - c_tot

c_hw_tot = np.sum(c_hw,axis=0).astype(np.float)
w_hw_tot = np.sum(w_hw,axis=0).astype(np.float)
d_hw_tot = w_hw_tot - c_hw_tot

c_chw_tot = np.sum(c_chw,axis=0).astype(np.float)
w_chw_tot = np.sum(w_chw,axis=0).astype(np.float)
d_chw_tot = w_chw_tot - c_chw_tot

c_ahw_tot = np.sum(c_ahw,axis=0).astype(np.float)
w_ahw_tot = np.sum(w_ahw,axis=0).astype(np.float)
d_ahw_tot = w_ahw_tot - c_ahw_tot

c_ahw2_tot = np.sum(c_ahw2,axis=0).astype(np.float)
w_ahw2_tot = np.sum(w_ahw2,axis=0).astype(np.float)

c_cp_tot = c_ahw2_tot/c_ahw_tot
w_cp_tot = w_ahw2_tot/w_ahw_tot
d_cp_tot = w_cp_tot-c_cp_tot
#%%
ndh=d_hw_tot/c_hw_tot
ndch=d_chw_tot/((c_chw_tot+w_chw_tot)/2)
#%%

#Shaded plot
plt.rcParams.update({'font.size': 12})
plt.figure(1)
dcp_ma = np.ma.masked_where(d_cp_tot==np.inf,d_cp_tot)
dcp_ma = np.ma.masked_invalid(dcp_ma)
cmap=copy(cm.jet)
cmap.set_bad('0.8',1.)
plt.pcolormesh(phis,np.square(scales)/warming,dcp_ma*100,cmap=cmap,vmin=0,vmax=100)
#plt.figtext(0.2,0.8,'mean=%.2f'% np.mean(dcp_ma*100),fontsize=18,color='white')
cb = plt.colorbar()
plt.ylim((0,2.5))
cb.set_label('2XCO2-Control Proportion\nof Compounded Days [%]')
plt.ylabel("[Standard Deviation]\n/[Mean Warming Signal]")
plt.xlabel("Lag 1-Day Autocorrelation")
#plt.title("90th 3214")
plt.savefig('3114syn.png')

##Shaded plot
#plt.figure(1)
#ndhma = np.ma.masked_where(ndh==np.inf,ndh)
#plt.pcolor(phis,np.square(scales)/warming,ndhma,cmap=cm.jet)
#plt.text(0.1,8.5,'mean=%.2f'% np.mean(ndh),fontsize=16,color='white')
#cb = plt.colorbar()
#plt.ylim((0,10))
#cb.set_label(r'$\Delta$HWF/HWF')
#plt.ylabel("[Standard Deviation]/[Mean Warming Signal]")
#plt.xlabel("Lag 1-Day Autocorrelation")
#plt.title("Non-compound 3")
##%%
#plt.figure(2)
#ndchma = np.ma.masked_where(ndch==np.inf,ndch)
#ndchma2 = np.ma.masked_where(np.isnan(ndchma),ndchma)
#plt.pcolor(phis,np.square(scales)/warming,ndchma2,cmap=cm.jet)
#plt.text(0.5,8.5,'mean=%.2f'% np.mean(ndchma2),fontsize=16,color='white')
#cb = plt.colorbar()
#plt.ylim((0,10))
#cb.set_label(r'$\Delta$HWF/HWF')
#plt.ylabel("[Standard Dev.]/[Mean Warming Signal]")
#plt.xlabel("Lag 1-Day Autocorrelation")
#plt.title("Compound 3-1")

##%%
##Shaded plot
#plt.figure(1)
#plt.contourf(phis,np.square(scales)/warming,d_tot/c_tot)
#cb = plt.colorbar()
#cb.set_label(r'$\Delta$HWF/HWF')
#plt.ylabel("[Standard Dev.]/[Mean Warming Signal]")
#plt.xlabel("Lag 1-Day Autocorrelation")
#plt.title("Non-compound, 1 day min. event")

##%%
##Plot extremes of time series
##Same Y axes
#plt.figure(2)
#plt.subplot(221)
#plt.plot(c[steps-1,0])
#plt.title("High Var., Low Autocorr.")
#plt.ylim((-120,120))
#
#plt.subplot(222)
#plt.plot(c[steps-1,steps-1])
#plt.title("High Var., High Autocorr.")
#plt.ylim((-120,120))
#
#plt.subplot(223)
#plt.plot(c[0,0])
#plt.xlabel("Low Var., Low Autocorr.",fontsize=12)
#plt.ylim((-120,120))
#
#plt.subplot(224)
#plt.plot(c[0,steps-1])
#plt.xlabel("Low Var., High Autocorr.",fontsize=12)
#plt.ylim((-120,120))
#
#plt.suptitle("Same Y Axes", fontsize=12)
#plt.show()
#
##Different Y axes
#plt.figure(2)
#plt.subplot(221)
#plt.plot(c[steps-1,0])
#plt.title("High Var., Low Autocorr.")
#
#plt.subplot(222)
#plt.plot(c[steps-1,steps-1])
#plt.title("High Var., High Autocorr.")
#
#plt.subplot(223)
#plt.plot(c[0,0])
#plt.xlabel("Low Var., Low Autocorr.",fontsize=12)
#
#plt.subplot(224)
#plt.plot(c[0,steps-1])
#plt.xlabel("Low Var., High Autocorr.",fontsize=12)
#
#plt.suptitle("Different Y Axes", fontsize=12)
#plt.show()


#%%





