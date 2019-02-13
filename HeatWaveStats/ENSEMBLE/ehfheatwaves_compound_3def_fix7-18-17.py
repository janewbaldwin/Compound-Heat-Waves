import warnings
warnings.filterwarnings('ignore')
import sys
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

# Parse command line arguments
usage = "usage: %prog -x <FILE> -n <FILE> -m <FILE> [options]"
parser = OptionParser(usage=usage)
parser.add_option('-x', '--tmax', dest='tmaxfile', 
        help='file containing tmax', metavar='FILE')
parser.add_option('--vnamex', dest='tmaxvname', default='tasmax',
        help='tmax variable name', metavar='STR')
parser.add_option('-n', '--tmin', dest='tminfile',
        help='file containing tmin', metavar='FILE')
parser.add_option('--vnamen', dest='tminvname', default='tasmin',
        help='tmin variable name', metavar='STR')
parser.add_option('-d', '--model', dest='model', default='FLOR',
        help='name of model that supplies data (ie FLOR, NCEP, ERAInt). Default is FLOR.',
        metavar='STR')
parser.add_option('-e', '--ensrn', dest='ensrn',
        help='ensemble realization number (ie 1-30).',
        metavar='STR')
parser.add_option('--bpfx', dest='bpfx',
                help=('Indicates a future simulation, specifying a tmax file '
                'containing the historical base period to be used'),
                metavar='FILE')
parser.add_option('--bpfn', dest='bpfn',
                help=('Indicates a future simulation, specifying a tmin file '
                'containing the historical base period to be used'),
                metavar='FILE')
parser.add_option('-m', '--mask', dest='maskfile',
        help='file containing land-sea mask', metavar='FILE')
parser.add_option('--vnamem', dest='maskvname', default='sftlf',
        help='mask variable name', metavar='STR')
parser.add_option('--vnamet', dest='timevname', default='time',
                help='time variable name', metavar='STR')
parser.add_option('-s', '--season', dest='season', default='summer',
        help='austal season for annual metrics. Defaults to austral summer',
        metavar='STR')
parser.add_option('-p', dest='pcntl', type='float', default=90,
        help='the percentile to use for thresholds. Defaults to 90',
        metavar='INT')
parser.add_option('--base', dest='bp', default='1961-1990',
        help='base period to calculate thresholds. Default 1961-1990',
        metavar='YYYY-YYYY')
parser.add_option('-q', '--qmethod', dest='qtilemethod', default='climpact',
        help='quantile interpolation method. Default is climpact', 
        metavar='STR')
parser.add_option('--t90pc', action="store_true", dest='t90pc',
                help='Calculate tx90pc and tn90pc heatwaves')
(options, args) = parser.parse_args()
if not options.tmaxfile or not options.tminfile:
    print "Please specify tmax and tmin files."
    sys.exit(2)
if not options.maskfile:
    print ("You didn't specify a land-sea mask. It's faster if you do,"
        "so this might take a while.")
if len(options.bp)!=9:
    print "Incorect base period format."
    sys.exit(2)
else:
    bpstart = int(options.bp[:4])
    bpend = int(options.bp[5:9])
# Percentile
pcntl = options.pcntl
# climpact/python/matlab
qtilemethod = options.qtilemethod
# season (winter/summer)
season = options.season
if (season!='summer')&(season!='winter'):
    print "Use either summer or winter."
    sys.exit(2)

# Load time data
try:
    tmaxnc = MFDataset(options.tmaxfile, 'r')
except IndexError:
    tmaxnc = Dataset(options.tmaxfile, 'r')
nctime = tmaxnc.variables[options.timevname]
try:
    nctime = MFTime(nctime)
except AttributeError:
    pass
calendar = nctime.calendar
if not calendar:
    print 'Unrecognized calendar. Using gregorian.'
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
    if tmaxnc.variables[options.timevname].units=='day as %Y%m%d.%f':
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

# Load land-sea mask
if options.maskfile:
    masknc = Dataset(options.maskfile, 'r')
    vname = options.maskvname
    mask = masknc.variables[vname][:]
    if mask.max()>1: mask = mask>50
    mask = mask.astype(np.bool)
    masknc.close()

# Load base period data
if options.bpfn:
    try:
        tminnc = MFDataset(options.bpfn, 'r')
    except IndexError:
        tminnc = Dataset(options.bpfn, 'r')
else:
    try:
        tminnc = MFDataset(options.tminfile, 'r')
    except IndexError:
        tminnc = Dataset(options.tminfile, 'r')
if options.bpfx:
    try:
        tmaxnc = MFDataset(options.bpfx, 'r')
    except IndexError:
        tmaxnc = Dataset(options.bpfx, 'r')
else:
    try:
        tmaxnc = MFDataset(options.tmaxfile, 'r')
    except IndexError:
        tmaxnc = Dataset(options.tmaxfile, 'r')
vname = options.tmaxvname
bptime = tmaxnc.variables[options.timevname]
try:
    bptime = MFTime(bptime)
except AttributeError:
    pass
if tmaxnc.variables[options.timevname].units=='day as %Y%m%d.%f':
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
if options.maskfile:
    mask=np.squeeze(mask) #added by Jane on 12/12/15
    tmax = tmax[:,mask]
if tmaxnc.variables[vname].units=='K': tmax -= 273.15
vname = options.tminvname
tmin = tminnc.variables[vname][(bpstart<=bpdates.year)&(bpdates.year<=bpend)]
if len(tmin.shape)==4: tmin = tmin.squeeze()
if options.maskfile:
    tmin = tmin[:,mask]
if tminnc.variables[vname].units=='K': tmin -= 273.15
tave_base = (tmax + tmin)/2.
if not options.t90pc:
    del tmin
    del tmax

# Remove leap days in gregorian calendars
if (calendar=='gregorian')|(calendar=='proleptic_gregorian')|\
            (calendar=='standard'):
    tave_base = tave_base[(dates_base.month!=2)|(dates_base.day!=29),...]
    if options.t90pc:
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
    if options.t90pc:
        txpct[day,...] = percentile(tmax[window,...], pcntl, parameter)
        tnpct[day,...] = percentile(tmin[window,...], pcntl, parameter)
    window = np.roll(window,1)
del tave_base
del window

# Load all data
try:
    tminnc = MFDataset(options.tminfile, 'r')
except IndexError:
    tminnc = Dataset(options.tminfile, 'r')
try:
    tmaxnc = MFDataset(options.tmaxfile, 'r')
except IndexError:
    tmaxnc = Dataset(options.tmaxfile, 'r')
tmax = tmaxnc.variables[options.tmaxvname][:]
if len(tmax.shape)==4: tmax = tmax.squeeze()
tmin = tminnc.variables[options.tminvname][:]
if len(tmin.shape)==4: tmin = tmin.squeeze()
if options.maskfile:
    tmax = tmax[:,mask]
    tmin = tmin[:,mask]
if tmaxnc.variables[options.tmaxvname].units=='K': tmax -= 273.15
if tminnc.variables[options.tminvname].units=='K': tmin -= 273.15
tave = (tmax + tmin)/2.
if not options.t90pc:
    del tmax
    del tmin

# Remove leap days from tave
if (calendar=='gregorian')|(calendar=='proleptic_gregorian')|\
            (calendar=='standard'):
    tave = tave[(dates.month!=2)|(dates.day!=29),...]
    if options.t90pc:
        tmax = tmax[(dates.month!=2)|(dates.day!=29),...]
        tmin = tmin[(dates.month!=2)|(dates.day!=29),...]

# Remove incomplete starting year
first_year = dayone.year
if (dayone.month!=1)|(dayone.day!=1):
    first_year = dayone.year+1
    start = np.argmax(dates.year==first_year)
    tave = tave[start:,...]
    if options.t90pc:
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

# Tx90pc exceedences
if options.t90pc:
    txexceed = np.ones(tmax.shape)*np.nan
    tnexceed = txexceed.copy()
    for i in xrange(0,tmax.shape[0]):
        idoy = i-daysinyear*int((i+1)/daysinyear)
        txexceed[i,...] = tmax[i,...]>txpct[idoy,...]
        tnexceed[i,...] = tmin[i,...]>tnpct[idoy,...]
    txexceed[txexceed>0] = tmax[txexceed>0]
    tnexceed[tnexceed>0] = tmin[tnexceed>0]

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

# Calculate metrics year by year
def split_hemispheres(EHF):
    """split_hemispheres splits the input data by hemispheres, and glues them
    back together after heatwave calculations.
    """
    if south:
        if options.maskfile:
            EHF_s = EHF[:,:(mask[lats<=0]>0).sum()]
        else:
            EHF_s = EHF[:,lats<=0,...]
        # Reshape to 2D
        space = EHF_s.shape[1:]
        if len(space)>1:
            EHF_s = EHF_s.reshape(EHF_s.shape[0],space[0]*space[1])
        # Southern hemisphere aspects
        HWN_s, HWF_s, HWD_s, CHWN_s, CHWF_s, CHWD_s, AHWN_s, AHWF_s, AHWD_s, AHW1N_s, AHW1F_s, AHW1D_s, AHW2F_s, AHW2D_s  = \
                hw_aspects(EHF_s, season, 'south')
        del EHF_s
    if north:
        if options.maskfile:
            EHF_n = EHF[:,(mask[lats<=0]>0).sum():]
        else:
            EHF_n = EHF[:,lats>0,...]
        # Reshape to 2D
        space = EHF_n.shape[1:]
        if len(space)>1:
            EHF_n = EHF_n.reshape(EHF_n.shape[0],space[0]*space[1])
        # Northern hemisphere aspects
        HWN_n, HWF_n, HWD_n, CHWN_n, CHWF_n, CHWD_n, AHWN_n, AHWF_n, AHWD_n, AHW1N_n, AHW1F_n, AHW1D_n, AHW2F_n, AHW2D_n = \
                hw_aspects(EHF_n, season, 'north')
        del EHF_n
    # Glue hemispheres back together
    if north and south:
        HWN = np.append(HWN_s, HWN_n, axis=1)
        HWF = np.append(HWF_s, HWF_n, axis=1)
        HWD = np.append(HWD_s, HWD_n, axis=1)
        CHWN = np.append(CHWN_s, CHWN_n, axis=1)
	CHWF = np.append(CHWF_s, CHWF_n, axis=1)
	CHWD = np.append(CHWD_s, CHWD_n, axis=1)	
	AHWN = np.append(AHWN_s, AHWN_n, axis=1)
	AHWF = np.append(AHWF_s, AHWF_n, axis=1)
	AHWD = np.append(AHWD_s, AHWD_n, axis=1)
	AHW1N = np.append(AHW1N_s, AHW1N_n, axis=1)
	AHW1F = np.append(AHW1F_s, AHW1F_n, axis=1)
	AHW1D = np.append(AHW1D_s, AHW1D_n, axis=1)
	AHW2F = np.append(AHW2F_s, AHW2F_n, axis=1)
	AHW2D = np.append(AHW2D_s, AHW2D_n, axis=1)
    elif north:
        HWN = HWN_n
        HWF = HWF_n
        HWD = HWD_n
        CHWN = CHWN_n
        CHWF = CHWF_n
        CHWD = CHWD_n
        AHWN = AHWN_n
        AHWF = AHWF_n
        AHWD = AHWD_n
        AHW1N = AHW1N_n
        AHW1F = AHW1F_n
        AHW1D = AHW1D_n
        AHW2F = AHW2F_n
        AHW2D = AHW2D_n
    elif south:
        HWN = HWN_s
        HWF = HWF_s
        HWD = HWD_s
        CHWN = CHWN_s
        CHWF = CHWF_s
        CHWD = CHWD_s
        AHWN = AHWN_s
        AHWF = AHWF_s
        AHWD = AHWD_s
        AHW1N = AHW1N_s
        AHW1F = AHW1F_s
        AHW1D = AHW1D_s
        AHW2F = AHW2F_s
        AHW2D = AHW2D_s
    return HWN, HWF, HWD, CHWN, CHWF, CHWD, AHWN, AHWF, AHWD, AHW1N, AHW1F, AHW1D, AHW2F, AHW2D


# Split by latitude
try:
    lats = tmaxnc.variables['lat'][:]
except KeyError:
    lats = tmaxnc.variables['latitude'][:]
north = (lats>0).any()
south = (lats<=0).any()
HWN_EHF, HWF_EHF, HWD_EHF, CHWN_EHF, CHWF_EHF, CHWD_EHF, AHWN_EHF, AHWF_EHF, AHWD_EHF, AHW1N_EHF, AHW1F_EHF, AHW1D_EHF, AHW2F_EHF, AHW2D_EHF = \
        split_hemispheres(EHF)
if options.t90pc:
    HWN_tx, HWF_tx, HWD_tx, CHWN_tx, CHWF_tx, CHWD_tx, AHWN_tx, AHWF_tx, AHWD_tx, AHW1N_tx, AHW1F_tx, AHW1D_tx, AHW2F_tx, AHW2D_tx = \
            split_hemispheres(txexceed)
    HWN_tn, HWF_tn, HWD_tn, CHWN_tn, CHWF_tn, CHWD_tn, AHWN_tn, AHWF_tn, AHWD_tn, AHW1N_tn, AHW1F_tn, AHW1D_tn, AHW2F_tn, AHW2D_tn = \
            split_hemispheres(tnexceed)

# Save to netCDF
try:
    experiment = str(elength).rstrip('0').rstrip('.')+str(bmax).rstrip('0').rstrip('.')+str(e2length).rstrip('0').rstrip('.')+str(dmin).rstrip('0').rstrip('.')
    model = options.model
    realization = options.ensrn
    if realization==u'N/A': realization = u'r1i1p1'
except AttributeError:
    experiment = ''
    model = ''
    realization = ''
try:
    space = (tmaxnc.dimensions['lat'].__len__(),tmaxnc.dimensions['lon'].__len__())
    lonname = 'lon'
    latname = 'lat'
except KeyError:
    lonname = 'longitude'
    latname = 'latitude'
    space = (tmaxnc.dimensions['latitude'].__len__(),tmaxnc.dimensions['longitude'].__len__())

def save_yearly(HWN, HWF, HWD, CHWN, CHWF, CHWD, AHWN, AHWF, AHWD, AHW1N, AHW1F, AHW1D, AHW2F, AHW2D, definition, basedat):
    yearlyout = Dataset('%s_heatwaves_%s_r%s_%s_yearly_%s.nc'%(definition,model,realization,experiment,season), 'w')
    yearlyout.createDimension('time', len(range(first_year,
            daylast.year+1)))
    yearlyout.createDimension('lon', tmaxnc.dimensions[lonname].__len__())
    yearlyout.createDimension('lat', tmaxnc.dimensions[latname].__len__())
    yearlyout.createDimension('day', daysinyear)
    setattr(yearlyout, "author", "Tammas Loughran")
    setattr(yearlyout, "contact", "t.loughran@student.unsw.edu.au")
    setattr(yearlyout, "source", "https://github.com/tammasloughran/ehfheatwaves")
    setattr(yearlyout, "date", dt.datetime.today().strftime('%Y-%m-%d'))
    setattr(yearlyout, "script", sys.argv[0])
    if model:
        setattr(yearlyout, "model_id", model)
        setattr(yearlyout, "experiment", experiment)
        setattr(yearlyout, "realization", "%s"%(realization))
    setattr(yearlyout, "period", "%s-%s"%(str(first_year),str(daylast.year)))
    setattr(yearlyout, "base_period", "%s-%s"%(str(bpstart),str(bpend)))
    setattr(yearlyout, "percentile", "%sth"%(str(pcntl)))
    setattr(yearlyout, "definition", definition)
    setattr(yearlyout, "frequency", "yearly")
    setattr(yearlyout, "season", season)
    setattr(yearlyout, "season_note", ("The year of a season is the year it starts"
            "in. SH summer: Nov-Mar. NH summer: May-Sep."))
    try:
        file = open('version', 'r')
        commit = file.read()[:-2]
    except IOError:
        commit = "Unknown. Check date for latest version."
    setattr(yearlyout, "git_commit", commit)
    setattr(yearlyout, "tmax_file", options.tmaxfile)
    setattr(yearlyout, "tmin_file", options.tminfile)
    if options.maskfile:
        setattr(yearlyout, "mask_file", options.maskfile)
    setattr(yearlyout, "quantile_method", options.qtilemethod)
    otime = yearlyout.createVariable('time', 'f8', 'time', 
            fill_value=-999.99)
    setattr(otime, 'units', 'year')
    olat = yearlyout.createVariable('lat', 'f8', 'lat')
    setattr(olat, 'standard_name', 'latitude')
    setattr(olat, 'long_name', 'Latitude')
    setattr(olat, 'units', 'degrees_north')
    setattr(olat, 'axis', 'Y')
    olon = yearlyout.createVariable('lon', 'f8', 'lon')
    setattr(olon, 'standard_name', 'longiitude')
    setattr(olon, 'long_name', 'Longitude')
    setattr(olon, 'units', 'degrees_east')
    setattr(olon, 'axis', 'X')
    otpct = yearlyout.createVariable('%s%spct'%(basedat,str(pcntl).rstrip('0').rstrip('.')), 'f8', 
	    ('day','lat','lon'), fill_value=-999.99)
    setattr(otpct, 'long_name', '%s %sth percentile'%(basedat,str(pcntl).rstrip('0').rstrip('.')))
    setattr(otpct, 'units', 'degC')
    setattr(otpct, 'description', 
            '%s %sth percentile of %s-%s'%(basedat,str(pcntl).rstrip('0').rstrip('.'),str(bpstart),str(bpend)))
    HWNout = yearlyout.createVariable('HWN_%s'%(definition), 'f8', ('time', 'lat', 'lon'), 
            fill_value=-999.99)
    setattr(HWNout, 'long_name', 'Non-compound Heatwave Number')
    setattr(HWNout, 'units','heatwaves')
    setattr(HWNout, 'description', 'Number of heatwaves per year')
    HWFout = yearlyout.createVariable('HWF_%s'%(definition), 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(HWFout, 'long_name','Non-compound Heatwave Frequency')
    setattr(HWFout, 'units', 'days')
    setattr(HWFout, 'description', 'Proportion of non-compound heatwave days per season')
    HWDout = yearlyout.createVariable('HWD_%s'%(definition), 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(HWDout, 'long_name', 'Non-compound Heatwave Duration')
    setattr(HWDout, 'units', 'days')
    setattr(HWDout, 'description', 'Duration of the longest non-compound heatwave per year')
    CHWNout = yearlyout.createVariable('CHWN_%s'%(definition), 'f8', ('time', 'lat', 'lon'), 
            fill_value=-999.99)
    setattr(CHWNout, 'long_name', 'Compound Heatwave Number')
    setattr(CHWNout, 'units','heatwaves')
    setattr(CHWNout, 'description', 'Number of compound heatwaves per year')
    CHWFout = yearlyout.createVariable('CHWF_%s'%(definition), 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(CHWFout, 'long_name','Compound Heatwave Frequency')
    setattr(CHWFout, 'units', 'days')
    setattr(CHWFout, 'description', 'Proportion of compound heatwave days per season')
    CHWDout = yearlyout.createVariable('CHWD_%s'%(definition), 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(CHWDout, 'long_name', 'Compound Heatwave Duration')
    setattr(CHWDout, 'units', 'days')
    setattr(CHWDout, 'description', 'Duration of the longest compound heatwave per year')
    AHWNout = yearlyout.createVariable('AHWN_%s'%(definition), 'f8', ('time', 'lat', 'lon'), 
            fill_value=-999.99)
    setattr(AHWNout, 'long_name', 'All Heatwave Number')
    setattr(AHWNout, 'units','heatwaves')
    setattr(AHWNout, 'description', 'Number of compound and non-compound heatwaves per year')
    AHWFout = yearlyout.createVariable('AHWF_%s'%(definition), 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(AHWFout, 'long_name','All Heatwave Frequency')
    setattr(AHWFout, 'units', 'days')
    setattr(AHWFout, 'description', 'Proportion of compound and non-compound heatwave days per season')
    AHWDout = yearlyout.createVariable('AHWD_%s'%(definition), 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(AHWDout, 'long_name', 'All Heatwave Duration')
    setattr(AHWDout, 'units', 'days')
    setattr(AHWDout, 'description', 'Duration of the longest of any heatwave (compound or non-compound) per year')
    AHW1Nout = yearlyout.createVariable('AHW1N_%s'%(definition), 'f8', ('time', 'lat', 'lon'), 
            fill_value=-999.99)
    setattr(AHW1Nout, 'long_name', 'Primary Heatwave Number')
    setattr(AHW1Nout, 'units','heatwaves')
    setattr(AHW1Nout, 'description', 'Number of heatwave starting events per year')
    AHW1Fout = yearlyout.createVariable('AHW1F_%s'%(definition), 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(AHW1Fout, 'long_name','Primary Heatwave Frequency')
    setattr(AHW1Fout, 'units', 'days')
    setattr(AHW1Fout, 'description', 'Proportion of heatwave days in starting heatwave per season')
    AHW1Dout = yearlyout.createVariable('AHW1D_%s'%(definition), 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(AHW1Dout, 'long_name', 'Primary Heatwave Duration')
    setattr(AHW1Dout, 'units', 'days')
    setattr(AHW1Dout, 'description', 'Duration of the longest of any starting heatwave (compound or non-compound) per year')
    AHW2Fout = yearlyout.createVariable('AHW2F_%s'%(definition), 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(AHW2Fout, 'long_name','Secondary Heatwave Frequency')
    setattr(AHW2Fout, 'units', 'days')
    setattr(AHW2Fout, 'description', 'Proportion of secondary heatwave days per season')
    AHW2Dout = yearlyout.createVariable('AHW2D_%s'%(definition), 'f8', ('time','lat','lon'), 
            fill_value=-999.99)
    setattr(AHW2Dout, 'long_name', 'Secondary Heatwave Duration')
    setattr(AHW2Dout, 'units', 'days')
    setattr(AHW2Dout, 'description', 'Duration of the longest number of secondary heat wave days for an individual heat wave (compound or non-compound) per year')
    otime[:] = range(first_year, daylast.year+1)
    olat[:] = tmaxnc.variables[latname][:]
    olon[:] = tmaxnc.variables[lonname][:]
    dummy_array = np.ones((daysinyear,)+original_shape[1:])*np.nan
    if options.maskfile:
        dummy_array[:,mask] = eval(basedat+'pct')
        dummy_array[np.isnan(dummy_array)] = -999.99
        otpct[:] = dummy_array.copy()
        dummy_array = np.ones((nyears,)+original_shape[1:])*np.nan
        dummy_array[:,mask] = HWN
        dummy_array[np.isnan(dummy_array)] = -999.99
        HWNout[:] = dummy_array.copy()
        dummy_array[:,mask] = HWF
        dummy_array[np.isnan(dummy_array)] = -999.99
        HWFout[:] = dummy_array.copy()
        dummy_array[:,mask] = HWD
        dummy_array[np.isnan(dummy_array)] = -999.99
        HWDout[:] = dummy_array.copy()
        dummy_array[:,mask] = CHWN
        dummy_array[np.isnan(dummy_array)] = -999.99
        CHWNout[:] = dummy_array.copy()
        dummy_array[:,mask] = CHWF
        dummy_array[np.isnan(dummy_array)] = -999.99
        CHWFout[:] = dummy_array.copy()
        dummy_array[:,mask] = CHWD
        dummy_array[np.isnan(dummy_array)] = -999.99
        CHWDout[:] = dummy_array.copy()
        dummy_array[:,mask] = AHWN
        dummy_array[np.isnan(dummy_array)] = -999.99
        AHWNout[:] = dummy_array.copy()
        dummy_array[:,mask] = AHWF
        dummy_array[np.isnan(dummy_array)] = -999.99
        AHWFout[:] = dummy_array.copy()
        dummy_array[:,mask] = AHWD
        dummy_array[np.isnan(dummy_array)] = -999.99
        AHWDout[:] = dummy_array.copy() 
        dummy_array[:,mask] = AHW1N
        dummy_array[np.isnan(dummy_array)] = -999.99
        AHW1Nout[:] = dummy_array.copy()
        dummy_array[:,mask] = AHW1F
        dummy_array[np.isnan(dummy_array)] = -999.99
        AHW1Fout[:] = dummy_array.copy()
        dummy_array[:,mask] = AHW1D
        dummy_array[np.isnan(dummy_array)] = -999.99
        AHW1Dout[:] = dummy_array.copy()
        dummy_array[:,mask] = AHW2F
        dummy_array[np.isnan(dummy_array)] = -999.99
        AHW2Fout[:] = dummy_array.copy()
        dummy_array[:,mask] = AHW2D
        dummy_array[np.isnan(dummy_array)] = -999.99
        AHW2Dout[:] = dummy_array.copy()
    else:
        otpct[:] = eval(basedat+'pct')
        HWNout[:] = HWN.reshape((nyears,)+space)
        HWFout[:] = HWF.reshape((nyears,)+space)
        HWDout[:] = HWD.reshape((nyears,)+space)
        CHWNout[:] = CHWN.reshape((nyears,)+space)
        CHWFout[:] = CHWF.reshape((nyears,)+space)
        CHWDout[:] = CHWD.reshape((nyears,)+space)
        AHWNout[:] = AHWN.reshape((nyears,)+space)
        AHWFout[:] = AHWF.reshape((nyears,)+space)
        AHWDout[:] = AHWD.reshape((nyears,)+space)
        AHW1Nout[:] = AHW1N.reshape((nyears,)+space)
        AHW1Fout[:] = AHW1F.reshape((nyears,)+space)
        AHW1Dout[:] = AHW1D.reshape((nyears,)+space)
        AHW2Fout[:] = AHW2F.reshape((nyears,)+space)
        AHW2Dout[:] = AHW2D.reshape((nyears,)+space)
    yearlyout.close()


save_yearly(HWN_EHF,HWF_EHF,HWD_EHF,CHWN_EHF,CHWF_EHF,CHWD_EHF,AHWN_EHF,AHWF_EHF,AHWD_EHF,AHW1N_EHF,AHW1F_EHF,AHW1D_EHF,AHW2F_EHF,AHW2D_EHF,'EHF%spct'%(str(pcntl).rstrip('0').rstrip('.')),'t')
if options.t90pc:
    save_yearly(HWN_tx,HWF_tx,HWD_tx,CHWN_tx,CHWF_tx,CHWD_tx,AHWN_tx,AHWF_tx,AHWD_tx,AHW1N_tx,AHW1F_tx,AHW1D_tx,AHW2F_tx,AHW2D_tx,'tx%spct'%(str(pcntl).rstrip('0').rstrip('.')),'tx')
    save_yearly(HWN_tn,HWF_tn,HWD_tn,CHWN_tn,CHWF_tn,CHWD_tn,AHWN_tn,AHWF_tn,AHWD_tn,AHW1N_tn,AHW1F_tn,AHW1D_tn,AHW2F_tn,AHW2D_tn,'tn%spct'%(str(pcntl).rstrip('0').rstrip('.')),'tn')

###########DEFINTION-el3-bmax3-e2l3-dm6####################
# Inputs
bmax = 3 #maximum length of a break
e2length = 3 #required length of heat waves that can compound
elength = 3 #required length of first heat wave(>=e2length)
dmin = 6 #required minimum number of hot days for compound event

# Split by latitude
HWN_EHF, HWF_EHF, HWD_EHF, CHWN_EHF, CHWF_EHF, CHWD_EHF, AHWN_EHF, AHWF_EHF, AHWD_EHF, AHW1N_EHF, AHW1F_EHF, AHW1D_EHF, AHW2F_EHF, AHW2D_EHF = \
        split_hemispheres(EHF)
if options.t90pc:
    HWN_tx, HWF_tx, HWD_tx, CHWN_tx, CHWF_tx, CHWD_tx, AHWN_tx, AHWF_tx, AHWD_tx, AHW1N_tx, AHW1F_tx, AHW1D_tx, AHW2F_tx, AHW2D_tx = \
            split_hemispheres(txexceed)
    HWN_tn, HWF_tn, HWD_tn, CHWN_tn, CHWF_tn, CHWD_tn, AHWN_tn, AHWF_tn, AHWD_tn, AHW1N_tn, AHW1F_tn, AHW1D_tn, AHW2F_tn, AHW2D_tn = \
            split_hemispheres(tnexceed)

# Save to netCDF
try:
    experiment = str(elength).rstrip('0').rstrip('.')+str(bmax).rstrip('0').rstrip('.')+str(e2length).rstrip('0').rstrip('.')+str(dmin).rstrip('0').rstrip('.')
except AttributeError:
    experiment = ''
save_yearly(HWN_EHF,HWF_EHF,HWD_EHF,CHWN_EHF,CHWF_EHF,CHWD_EHF,AHWN_EHF,AHWF_EHF,AHWD_EHF,AHW1N_EHF,AHW1F_EHF,AHW1D_EHF,AHW2F_EHF,AHW2D_EHF,'EHF%spct'%(str(pcntl).rstrip('0').rstrip('.')),'t')
if options.t90pc:
    save_yearly(HWN_tx,HWF_tx,HWD_tx,CHWN_tx,CHWF_tx,CHWD_tx,AHWN_tx,AHWF_tx,AHWD_tx,AHW1N_tx,AHW1F_tx,AHW1D_tx,AHW2F_tx,AHW2D_tx,'tx%spct'%(str(pcntl).rstrip('0').rstrip('.')),'tx')
    save_yearly(HWN_tn,HWF_tn,HWD_tn,CHWN_tn,CHWF_tn,CHWD_tn,AHWN_tn,AHWF_tn,AHWD_tn,AHW1N_tn,AHW1F_tn,AHW1D_tn,AHW2F_tn,AHW2D_tn,'tn%spct'%(str(pcntl).rstrip('0').rstrip('.')),'tn')

###########DEFINTION-el6-bmax2-e2l1-dm7####################
# Inputs
bmax = 2 #maximum length of a break
e2length = 1 #required length of heat waves that can compound
elength = 6 #required length of first heat wave(>=e2length)
dmin = 7 #required minimum number of hot days for compound event

# Split by latitude
HWN_EHF, HWF_EHF, HWD_EHF, CHWN_EHF, CHWF_EHF, CHWD_EHF, AHWN_EHF, AHWF_EHF, AHWD_EHF, AHW1N_EHF, AHW1F_EHF, AHW1D_EHF, AHW2F_EHF, AHW2D_EHF = \
        split_hemispheres(EHF)
if options.t90pc:
    HWN_tx, HWF_tx, HWD_tx, CHWN_tx, CHWF_tx, CHWD_tx, AHWN_tx, AHWF_tx, AHWD_tx, AHW1N_tx, AHW1F_tx, AHW1D_tx, AHW2F_tx, AHW2D_tx = \
            split_hemispheres(txexceed)
    HWN_tn, HWF_tn, HWD_tn, CHWN_tn, CHWF_tn, CHWD_tn, AHWN_tn, AHWF_tn, AHWD_tn, AHW1N_tn, AHW1F_tn, AHW1D_tn, AHW2F_tn, AHW2D_tn = \
            split_hemispheres(tnexceed)

# Save to netCDF
try:
    experiment = str(elength).rstrip('0').rstrip('.')+str(bmax).rstrip('0').rstrip('.')+str(e2length).rstrip('0').rstrip('.')+str(dmin).rstrip('0').rstrip('.')
except AttributeError:
    experiment = ''
save_yearly(HWN_EHF,HWF_EHF,HWD_EHF,CHWN_EHF,CHWF_EHF,CHWD_EHF,AHWN_EHF,AHWF_EHF,AHWD_EHF,AHW1N_EHF,AHW1F_EHF,AHW1D_EHF,AHW2F_EHF,AHW2D_EHF,'EHF%spct'%(str(pcntl).rstrip('0').rstrip('.')),'t')
if options.t90pc:
    save_yearly(HWN_tx,HWF_tx,HWD_tx,CHWN_tx,CHWF_tx,CHWD_tx,AHWN_tx,AHWF_tx,AHWD_tx,AHW1N_tx,AHW1F_tx,AHW1D_tx,AHW2F_tx,AHW2D_tx,'tx%spct'%(str(pcntl).rstrip('0').rstrip('.')),'tx')
    save_yearly(HWN_tn,HWF_tn,HWD_tn,CHWN_tn,CHWF_tn,CHWD_tn,AHWN_tn,AHWF_tn,AHWD_tn,AHW1N_tn,AHW1F_tn,AHW1D_tn,AHW2F_tn,AHW2D_tn,'tn%spct'%(str(pcntl).rstrip('0').rstrip('.')),'tn')
