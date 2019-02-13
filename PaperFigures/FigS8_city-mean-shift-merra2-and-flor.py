import os
os.chdir('/Users/tweetybaldwin/Dropbox/1. RESEARCH/STEP/City Analysis/')

import warnings
warnings.filterwarnings('ignore')
import sys
import matplotlib.pyplot as plt
try:
    modulename = 'pandas'
    import pandas as pd
    modulename = 'numpy'
    import numpy as np
    modulename = 'datetime'
    import datetime as dt
    modulename = 'math'
    import math
    modulename = 'qtiler'
    import qtiler
    modulename = 'netCDF4'
    from netCDF4 import MFDataset, MFTime, Dataset, num2date
    modulename = 'netcdftime'
    import netcdftime
    modulename = 'optparse'
    from optparse import OptionParser
    modulename = 'distutils.version'
    from distutils.version import LooseVersion
except ImportError:
    print modulename, " is missing. Please install missing packages."
    sys.exit(2)
if LooseVersion(np.__version__) < LooseVersion('1.8.0'):
    print "Please install numpy version 1.8.0 or higher."
    sys.exit(2)
#%% HEAT WAVE DEFINITION SCRIPT
         
def hw(tmaxfile,tmaxvname,tminfile,tminvname,bp,hemisphere,sfx,sfn):
    season = 'summer'
    timevname = 'TIME'
    pcntl = 90
    qtilemethod = 'climpact'
    
    bpstart = int(bp[:4])
    bpend = int(bp[5:9])
    
    # Load time data
    try:
        tmaxnc = MFDataset(tmaxfile, 'r')
    except IndexError:
        tmaxnc = Dataset(tmaxfile, 'r')
    nctime = tmaxnc.variables[timevname]
    try:
        nctime = MFTime(nctime)
    except AttributeError:
        pass
    calendar = nctime.calendar
    if not calendar:
        print 'Unrecognized calendar. Using gregcorian.'
        calendar = 'gregorian'
    elif calendar=='360_day':
        daysinyear = 360
        # 360 day season start and end indices
        SHS = (301,451)
        SHW = (121,271)
        dayone = num2date(nctime[0], nctime.units,
                calendar=calendar)
        daylast = num2date(nctime[-1], nctime.units,
                calendar=calendar)
        class calendar360():
            def __init__(self,sdate,edate):
                self.year = np.repeat(range(sdate.year,edate.year+1), 360, 0)
                nyears = len(xrange(sdate.year,edate.year+1))
                self.month = np.tile(np.repeat(range(1,12+1), 30, 0), nyears)
                self.day = np.tile(np.tile(range(1,30+1), 12), nyears)
                if (sdate.day!=1)|(edate.month!=1):
                    sdoyi = (sdate.month-1)*30+sdate.day-1
                    self.year = self.year[sdoyi:]
                    self.month = self.month[sdoyi:]
                    self.day = self.day[sdoyi:]
                if (edate.day!=30)|(edate.month!=12):
                    edoyi = (12-edate.month)*30+(30-edate.day)
                    self.year = self.year[:-edoyi]
                    self.month = self.month[:-edoyi]
                    self.day = self.day[:-edoyi]
        dates = calendar360(dayone, daylast)
        shorten = 0
        if (daylast.day!=30)|(daylast.month!=12):
            shorten = 30*(13-daylast.month) - daylast.day
    else:
        daysinyear = 365
        # 365 day season start and end indices
        SHS = (304,455)
        SHW = (121,274)
        if tmaxnc.variables[timevname].units=='day as %Y%m%d.%f':
            st = str(int(nctime[0]))
            nd = str(int(nctime[-1]))
            dayone = dt.datetime(int(st[:4]), int(st[4:6]), int(st[6:]))
            daylast = dt.datetime(int(nd[:4]), int(nd[4:6]), int(nd[6:]))
        else:
            dayone = num2date(nctime[0], nctime.units,
                    calendar=calendar)
            daylast = num2date(nctime[-1], nctime.units,
                    calendar=calendar)
        dates = pd.period_range(str(dayone), str(daylast))
        shorten = 0
        if (daylast.day!=30)|(daylast.month!=12):
            endofdata = dt.datetime(2000, daylast.month, daylast.day)
            shorten = dt.datetime(2000, 12, 31) - endofdata
            shorten = shorten.days
    
    # Load base period data
    try:
        tminnc = MFDataset(tminfile, 'r')
    except IndexError:
        tminnc = Dataset(tminfile, 'r')
    try:
        tmaxnc = MFDataset(tmaxfile, 'r')
    except IndexError:
        tmaxnc = Dataset(tmaxfile, 'r')
    vname = tmaxvname
    bptime = tmaxnc.variables[timevname]
    try:
        bptime = MFTime(bptime)
    except AttributeError:
        pass
    if tmaxnc.variables[timevname].units=='day as %Y%m%d.%f':
        st = str(int(bptime[0]))
        nd = str(int(bptime[-1]))
        bpdayone = dt.datetime(int(st[:4]), int(st[4:6]), int(st[6:]))
        bpdaylast = dt.datetime(int(nd[:4]), int(nd[4:6]), int(nd[6:]))
    else:
        bpdayone = num2date(bptime[0], bptime.units, calendar=calendar)
        bpdaylast = num2date(bptime[-1], bptime.units, calendar=calendar)
    if calendar=='360_day': bpdates = calendar360(bpdayone, bpdaylast)
    else: 
        bpdates = pd.period_range(str(bpdayone), str(bpdaylast))
        dates_base = bpdates[(bpstart<=bpdates.year)&(bpdates.year<=bpend)]
    tmax = tmaxnc.variables[vname][(bpstart<=bpdates.year)&(bpdates.year<=bpend)]
    if len(tmax.shape)==4: tmax = tmax.squeeze()
    original_shape = tmax.shape
    tmax -= 273.15
    vname = tminvname
    tmin = tminnc.variables[vname][(bpstart<=bpdates.year)&(bpdates.year<=bpend)]
    if len(tmin.shape)==4: tmin = tmin.squeeze()
    tmin -= 273.15
    tave_base = (tmax + tmin)/2.
    
    # Remove leap days in gregorian calendars
    if (calendar=='gregorian')|(calendar=='proleptic_gregorian')|\
                (calendar=='standard'):
        tave_base = tave_base[(dates_base.month!=2)|(dates_base.day!=29),...]
        tmax = tmax[(dates_base.month!=2)|(dates_base.day!=29),...]
        tmin = tmin[(dates_base.month!=2)|(dates_base.day!=29),...]
        del dates_base
    
    # Caclulate 90th percentile
    tpct = np.ones(((daysinyear,)+tave_base.shape[1:]))*np.nan
    txpct = tpct.copy()
    tnpct = tpct.copy()
    window = np.zeros(daysinyear,dtype=np.bool)
    wsize = 15.
    window[-np.floor(wsize/2.):] = 1
    window[:np.ceil(wsize/2.)] = 1
    window = np.tile(window,bpend+1-bpstart)
    if qtilemethod=='python':
        percentile = np.percentile
        parameter = 0
    elif qtilemethod=='zhang':
        percentile = qtiler.quantile_zhang
        parameter = False
    elif qtilemethod=='matlab':
        percentile = qtiler.quantile_R
        parameter = 5
    elif qtilemethod=='climpact':
        percentile = qtiler.quantile_climpact
        parameter = False
    for day in xrange(daysinyear):
        tpct[day,...] = percentile(tave_base[window,...], pcntl, parameter)
        txpct[day,...] = percentile(tmax[window,...], pcntl, parameter)
        tnpct[day,...] = percentile(tmin[window,...], pcntl, parameter)
        window = np.roll(window,1)
    del tave_base
    del window
    
    # Load all data
    try:
        tminnc = MFDataset(tminfile, 'r')
    except IndexError:
        tminnc = Dataset(tminfile, 'r')
    try:
        tmaxnc = MFDataset(tmaxfile, 'r')
    except IndexError:
        tmaxnc = Dataset(tmaxfile, 'r')
    tmax = tmaxnc.variables[tmaxvname][:]
    if len(tmax.shape)==4: tmax = tmax.squeeze()
    tmin = tminnc.variables[tminvname][:]
    if len(tmin.shape)==4: tmin = tmin.squeeze()
    tmax -= 273.15
    tmin -= 273.15
    tave = (tmax + tmin)/2.
    
    # Remove leap days from tave
    if (calendar=='gregorian')|(calendar=='proleptic_gregorian')|\
                (calendar=='standard'):
        tave = tave[(dates.month!=2)|(dates.day!=29),...]
        tmax = tmax[(dates.month!=2)|(dates.day!=29),...]
        tmin = tmin[(dates.month!=2)|(dates.day!=29),...]
    
    # Remove incomplete starting year
    first_year = dayone.year
    if (dayone.month!=1)|(dayone.day!=1):
        first_year = dayone.year+1
        start = np.argmax(dates.year==first_year)
        tave = tave[start:,...]
        tmax = tmax[start:,...]
        tmin = tmin[start:,...]
    
    # Calculate EHF
    EHF = np.ones(tave.shape)*np.nan
    for i in xrange(32,tave.shape[0]):
        EHIaccl = tave[i-2:i+1,...].sum(axis=0)/3. - \
                tave[i-32:i-2,...].sum(axis=0)/30.
        EHIsig = tave[i-2:i+1,...].sum(axis=0)/3. - \
                tpct[i-daysinyear*int((i+1)/daysinyear),...]
        EHF[i,...] = np.maximum(EHIaccl,1.)*EHIsig
    EHF[EHF<0] = 0
    
    ###########DEFINTION-el3-bmax1-e2l1-dm4####################
    # Inputs
    bmax = 1 #maximum length of a break
    e2length = 1 #required length of heat waves that can compound
    elength = 3 #required length of first heat wave(>=e2length)
    dmin = 4 #required minimum number of hot days for compound event
    
    def identify_hw(ehfs):
        """identify_hw locates heatwaves from EHF and returns duration indicators hw, chw, ahw.
        """
         # Agregate consecutive days with EHF>0, and consecutive days where EHF<0 (not heat wave)
        events = (ehfs>0.).astype(np.int)
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
        chw = np.append(chw,np.expand_dims(np.zeros(ehfs[0,...].shape),axis=0),0) #Append additional day
        endlength=chw.copy()
        size = ehfs.shape[0]
        for i in xrange(lasttime,-1,-1):
            xy=np.logical_and(chw[i,...]>0,breakafter[i,...]<=bmax) #only evaluation for first days of events part of compound events, where the break after is less than max break length
            eduration=endes[i,xy]    
            bafter=breakafter[i,xy]
            eaddi=i+eduration+bafter
            eaddi[eaddi>lasttime]=size #first day endss always=0, so use this index
            chw[i,xy]=eduration+chw[eaddi,xy] #no. of event days in compound heat wave
            endlength[i,xy]=eduration+bafter+endlength[eaddi,xy] #no. of days (breaks and evnets) in compound heat wave
        chw = np.delete(chw,size,0) #Remove additional day with value 0
        endlength = np.delete(endlength,size,0) #Remove additional day with value 0
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
                    if chw[i,index[0]]>0: #add ,index[1] after  index[0] if 3-D array
                        length=int(endlength[i,index[0]]) #add ,index[1] after  index[0] if 3-D array
                        elasti=i+length-1
                        chw[(i+1):elasti,index[0]]=0 #add ,index[1] after  index[0] if 3-D array
                        ahw[(i+1):elasti,index[0]]=0
        #secondary ahw (events that compound on to other events)
        ahw2 = chw-hw
        ahw2[ahw2<0] = 0
        #primary ahw (compound starting events, and events with no compounding)
        ahw1 = ahw-ahw2
        return hw, chw, ahw, ahw1, ahw2    
    
    # Shifting Tmin and Tmax
    tmax = tmax + sfx
    tmin = tmin + sfn
    
    # Tx90pc exceedences
    txexceed = np.ones(tmax.shape)*np.nan
    tnexceed = txexceed.copy()
    for i in xrange(0,tmax.shape[0]):
        idoy = i-daysinyear*int((i+1)/daysinyear)
        txexceed[i,...] = tmax[i,...]>txpct[idoy,...]
        tnexceed[i,...] = tmin[i,...]>tnpct[idoy,...]
    txexceed[txexceed>0] = tmax[txexceed>0]
    tnexceed[tnexceed>0] = tmin[tnexceed>0]
    txexceed = np.squeeze(txexceed)
    tnexceed = np.squeeze(tnexceed)
    
    nyears = len(range(first_year,daylast.year+1))
    
    def hw_aspects(EHF, season, hemisphere):
        """hw_aspects takes EHF values or temp 90pct exceedences identifies
        heatwaves and calculates seasonal aspects.
        """
        # Select indices depending on calendar season and hemisphere
        if season=='summer':
            if hemisphere=='south':
                startday = SHS[0]
                endday = SHS[1]
            else:
                startday = SHW[0]
                endday = SHW[1]
        elif season=='winter':
            if hemisphere=='south':
                startday = SHW[0]
                endday = SHW[1]
            else:
                startday = SHS[0]
                endday = SHS[1]
        # Initialize arrays
        HWN = np.ones(((nyears,)+(EHF.shape[1],)))*np.nan
        HWF = HWN.copy()
        HWD = HWN.copy()
        CHWN = HWN.copy()
        CHWF = HWN.copy()
        CHWD = HWN.copy()
        AHWN = HWN.copy()
        AHWF = HWN.copy()
        AHWD = HWN.copy()
        AHW1N = HWN.copy()
        AHW1F = HWN.copy()
        AHW1D = HWN.copy()
        AHW2F = HWN.copy()
        AHW2D = HWN.copy()
        # Loop over years
        for iyear, year in enumerate(xrange(first_year,daylast.year)):
            if (year==daylast.year): continue # Incomplete yr
            # Select this years season
            allowance = 0 # For including heatwave days after the end of the season
            ifrom = startday + daysinyear*iyear
            ito = endday + daysinyear*iyear + allowance
            EHF_i = EHF[ifrom:ito,...]
            hw_i, chw_i, ahw_i, ahw1_i, ahw2_i = identify_hw(EHF_i)
            # Remove events that start after the end of the season and before start
            EHF_i = EHF_i[0:,...]
            hw_i = hw_i[0:]  #hw_i = hw_i[0:-allowance]
            chw_i = chw_i[0:] #hw_i = hw_i[0:-allowance]
            ahw_i = ahw_i[0:] #hw_i = hw_i[0:-allowance]
            ahw1_i = ahw1_i[0:] #hw_i = hw_i[0:-allowance]
            ahw2_i = ahw2_i[0:] #hw_i = hw_i[0:-allowance]
            # Calculate metrics
            HWN[iyear,...] = (hw_i>0).sum(axis=0)
            HWF[iyear,...] = hw_i.sum(axis=0)
            HWD[iyear,...] = hw_i.max(axis=0)
            CHWN[iyear,...] = (chw_i>0).sum(axis=0)
            CHWF[iyear,...] = chw_i.sum(axis=0)
            CHWD[iyear,...] = chw_i.max(axis=0)
            AHWN[iyear,...] = (ahw_i>0).sum(axis=0)
            AHWF[iyear,...] = ahw_i.sum(axis=0)
            AHWD[iyear,...] = ahw_i.max(axis=0)
            AHW1N[iyear,...] = (ahw1_i>0).sum(axis=0)
            AHW1F[iyear,...] = ahw1_i.sum(axis=0)
            AHW1D[iyear,...] = ahw1_i.max(axis=0)
            AHW2F[iyear,...] = ahw2_i.sum(axis=0)
            AHW2D[iyear,...] = ahw2_i.max(axis=0)
        return HWN, HWF, HWD, CHWN, CHWF, CHWD, AHWN, AHWF, AHWD, AHW1N, AHW1F, AHW1D, AHW2F, AHW2D
    
    
    HWN, HWF, HWD, CHWN, CHWF, CHWD, AHWN, AHWF, AHWD, AHW1N, AHW1F, AHW1D, AHW2F, AHW2D = hw_aspects(tnexceed,season,hemisphere)

    return AHWF, AHW2F

#%% CALCULATION OF HEAT WAVES
    
sfn = np.arange(0,6.25,0.25)
sfx = np.arange(0,6.25,0.25)

#MERRA2 Chicago
tmaxfile = 'merra2_chicago.nc'#'merra2_manaus.nc'
tmaxvname = 'T2MMAX_FLORGRID' #'T2MMAX_FLORGRID'#'T_REF_MAX'
tminfile = 'merra2_chicago.nc'#'merra2_manaus.nc'
tminvname = 'T2MMIN_FLORGRID'
bp = '1981-2010'
hemisphere = 'north'
AHWF_mc, AHW2F_mc = hw(tmaxfile,tmaxvname,tminfile,tminvname,bp,hemisphere,sfx,sfn)

#MERRA2 Manaus
tmaxfile = 'merra2_manaus.nc'#'merra2_manaus.nc'
tmaxvname = 'T2MMAX_FLORGRID' #'T2MMAX_FLORGRID'#'T_REF_MAX'
tminfile = 'merra2_manaus.nc'#'merra2_manaus.nc'
tminvname = 'T2MMIN_FLORGRID'
bp = '1981-2010'
hemisphere = 'south'
AHWF_mm, AHW2F_mm = hw(tmaxfile,tmaxvname,tminfile,tminvname,bp,hemisphere,sfx,sfn)

#FLOR Chicago
tmaxfile = 'flor_control_chicago.nc'
tmaxvname = 'T_REF_MAX' 
tminfile = 'flor_control_chicago.nc'
tminvname = 'T_REF_MIN'
bp = '0401-0430'
hemisphere = 'north'
AHWF_fc, AHW2F_fc = hw(tmaxfile,tmaxvname,tminfile,tminvname,bp,hemisphere,sfx,sfn)

#FLOR Manaus
tmaxfile = 'flor_control_manaus.nc'
tmaxvname = 'T_REF_MAX' 
tminfile = 'flor_control_manaus.nc'
tminvname = 'T_REF_MIN'
bp = '0401-0430'
hemisphere = 'south'
AHWF_fm, AHW2F_fm = hw(tmaxfile,tmaxvname,tminfile,tminvname,bp,hemisphere,sfx,sfn)

#%% PLOTTING

#set up plot
#f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)# sharex='col', sharey='row')
fig = plt.figure(figsize=(9,6))

#MERRA2 Chicago
ax1 = plt.subplot(2,2,1)
AHWF = AHWF_mc
AHW2F = AHW2F_mc
ax1.plot(sfn, np.mean(AHWF[:34,:],axis=0), 'b-', label= "All Heat Wave Days")
ax1.plot(sfn, np.mean(AHW2F[:34,:],axis=0),'b--', label="Compounded Days")
#ax1.set_xlabel('Warming [K]')
ax1.set_ylabel('Days', color='b')
ax1.tick_params('y', colors='b')
ax1b = ax1.twinx()
dcp = np.mean(AHW2F[:34,:],axis=0)/np.mean(AHWF[:34,:],axis=0)*100
ax1b.plot(sfn, dcp, 'r--', label = "Compound Proportion")
ax1b.plot(sfn, np.mean(AHWF[:34,:]/153*100,axis=0), 'r.', label = "Summer in Heat Waves")
#ax1b.set_ylabel('Percent', color='r')
ax1b.tick_params('y', colors='r')
ax1.set_ylim((0,153))
ax1b.set_ylim((0,100))

#MERRA2 Manaus
ax2 = plt.subplot(2,2,2)
AHWF = AHWF_mm
AHW2F = AHW2F_mm
ax2.plot(sfn, np.mean(AHWF[:34,:],axis=0), 'b-', label= "All Heat Wave Days")
ax2.plot(sfn, np.mean(AHW2F[:34,:],axis=0),'b--', label="Compounded Days")
#ax2.set_xlabel('Warming [K]')
#ax2.set_ylabel('Days', color='b')
ax2.tick_params('y', colors='b')
ax2b = ax2.twinx()
dcp = np.mean(AHW2F[:34,:],axis=0)/np.mean(AHWF[:34,:],axis=0)*100
ax2b.plot(sfn, dcp, 'r--', label = "Compound Proportion")
ax2b.plot(sfn, np.mean(AHWF[:34,:]/153*100,axis=0), 'r.', label = "Summer in Heat Waves")
ax2b.set_ylabel('Percent', color='r')
ax2b.tick_params('y', colors='r')
ax2.set_ylim((0,153))
ax2b.set_ylim((0,100))

#FLOR Chicago
ax3 = plt.subplot(2,2,3)
AHWF = AHWF_fc
AHW2F = AHW2F_fc
ax3.plot(sfn, np.mean(AHWF[:34,:],axis=0), 'b-', label= "All Heat Wave Days")
ax3.plot(sfn, np.mean(AHW2F[:34,:],axis=0),'b--', label="Compounded Days")
ax3.set_xlabel('Warming [K]')
ax3.set_ylabel('Days', color='b')
ax3.tick_params('y', colors='b')
ax3b = ax3.twinx()
dcp = np.mean(AHW2F[:34,:],axis=0)/np.mean(AHWF[:34,:],axis=0)*100
ax3b.plot(sfn, dcp, 'r--', label = "Compound Proportion")
ax3b.plot(sfn, np.mean(AHWF[:34,:]/153*100,axis=0), 'r.', label = "Summer in Heat Waves")
#ax3b.set_ylabel('Percent', color='r')
ax3b.tick_params('y', colors='r')
ax3.set_ylim((0,153))
ax3b.set_ylim((0,100))

#FLOR Manaus
ax4 = plt.subplot(2,2,4)
AHWF = AHWF_fm
AHW2F = AHW2F_fm
ax4.plot(sfn, np.mean(AHWF[:34,:],axis=0), 'b-', label= "All Heat Wave Days")
ax4.plot(sfn, np.mean(AHW2F[:34,:],axis=0),'b--', label="Compounded Days")
ax4.set_xlabel('Warming [K]')
#ax4.set_ylabel('Days', color='b')
ax4.tick_params('y', colors='b')
ax4b = ax4.twinx()
dcp = np.mean(AHW2F[:34,:],axis=0)/np.mean(AHWF[:34,:],axis=0)*100
ax4b.plot(sfn, dcp, 'r--', label = "Compound Proportion")
ax4b.plot(sfn, np.mean(AHWF[:34,:]/153*100,axis=0), 'r.', label = "Summer in Heat Waves")
ax4b.set_ylabel('Percent', color='r')
ax4b.tick_params('y', colors='r')
ax4.set_ylim((0,153))
ax4b.set_ylim((0,100))

plt.tight_layout(pad=3,h_pad=1,w_pad=1)

# Shrink current axis by 20%
#box = ax1.get_position()
#ax1.set_position([box.x0, box.y0, box.width, box.height])
#ax2.set_position([box.x0, box.y0, box.width, box.height])

# Put a legend to the right of the current axis
handles1, labels1 = ax1.get_legend_handles_labels()
handles1b, labels1b = ax1b.get_legend_handles_labels()
handles = np.append(handles1,handles1b)
labels = np.append(labels1,labels1b)
ax1.legend(handles, labels,loc='upper right',fontsize=9)

#handles, labels = ax1b.get_legend_handles_labels()
#ax1b.legend(handles, labels,loc='center left')

ax1.text(-0.23,0.5,'MERRA2',transform=ax1.transAxes,fontsize=16,rotation='vertical',ha='center',va='center')
ax3.text(-0.23,0.5,'FLOR 1990 Control',transform=ax3.transAxes,fontsize=16,rotation='vertical',ha='center',va='center')
ax1.text(0.5,1.1,'Chicago (41.8$^{\circ}$N,87.6$^{\circ}$W)',transform=ax1.transAxes,fontsize=16,rotation='horizontal',ha='center',va='center')
ax2.text(0.5,1.1,'Manaus (3.1$^{\circ}$S, 60.0$^{\circ}$W)',transform=ax2.transAxes,fontsize=16,rotation='horizontal',ha='center',va='center')

ax1.text(0.03,0.75,'std=3.47',transform=ax1.transAxes,fontsize=14,rotation='horizontal',ha='left',va='top')
ax2.text(0.03,0.75,'std=0.83',transform=ax2.transAxes,fontsize=14,rotation='horizontal',ha='left',va='top')
ax3.text(0.03,0.75,'std=3.80',transform=ax3.transAxes,fontsize=14,rotation='horizontal',ha='left',va='top')
ax4.text(0.03,0.75,'std=1.25',transform=ax4.transAxes,fontsize=14,rotation='horizontal',ha='left',va='top')

ax1.text(-0.12,132, "a", fontsize=20, bbox=dict(boxstyle='square',
            fc='0.8', alpha=1.0), zorder=100)
ax2.text(-0.12,132, "b", fontsize=20, bbox=dict(boxstyle='square',
            fc='0.8', alpha=1.0), zorder=100)
ax3.text(-0.12,132, "c", fontsize=20, bbox=dict(boxstyle='square',
            fc='0.8', alpha=1.0), zorder=100)
ax4.text(-0.12,132, "d", fontsize=20, bbox=dict(boxstyle='square',
            fc='0.8', alpha=1.0), zorder=100)

fig.savefig('cityshift.pdf')

plt.show()
