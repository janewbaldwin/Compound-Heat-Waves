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
parser.add_option('-d', '--daily', action="store_true", dest='daily', 
        help='output daily EHF values and heatwave indicators')
parser.add_option('--dailyonly', action="store_true", dest='dailyonly',
        help='output only daily values and suppress yearly output')
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
    bpstart = int(options.bp[1:4])
    bpend = int(options.bp[6:9])
# Percentile
pcntl = options.pcntl
# climpact/python/matlab
qtilemethod = options.qtilemethod
# season (winter/summer)
season = options.season
if (season!='summer')&(season!='winter'):
    print "Use either summer or winter."
    sys.exit(2)
# save daily EHF output
yearlyout = True
dailyout = options.daily
if options.dailyonly: 
    dailyout = True
    yearlyout = False

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
bpdayone = num2date(bptime[0], bptime.units, calendar=calendar)
bpdaylast = num2date(bptime[-1], bptime.units, calendar=calendar)
bpdates = pd.period_range(str(bpdayone), str(bpdaylast))
dates_base = bpdates[(bpstart<=bpdates.year)&(bpdates.year<=bpend)]
tmax = tmaxnc.variables[vname][(bpstart<=dates.year)&(dates.year<=bpend)]
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

# Remove leap days in gregorian calendars
if (calendar=='gregorian')|(calendar=='proleptic_gregorian')|\
            (calendar=='standard')|(calendar=='julian'):
    tave_base = tave_base[(dates_base.month!=2)|(dates_base.day!=29),...]
    tmax = tmax[(dates_base.month!=2)|(dates_base.day!=29),...]
    tmin = tmin[(dates_base.month!=2)|(dates_base.day!=29),...]
    del dates_base

# Remove incomplete starting year
first_year = dayone.year
if (dayone.month!=1)|(dayone.day!=1):
    first_year = dayone.year+1
    start = np.argmax(dates.year==first_year)
    tave = tave[start:,...]
    tmax = tmax[start:,...]
    tmin = tmin[start:,...]

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
    txpct[day,...] = percentile(tmax[window,...], pcntl, parameter)
    tnpct[day,...] = percentile(tmin[window,...], pcntl, parameter)
    window = np.roll(window,1)
del tave_base
del window

# Save to netCDF
try:
    experiment = tmaxnc.__getattribute__('experiment')
    model = tmaxnc.__getattribute__('model_id')
    realization = tmaxnc.__getattribute__('parent_experiment_rip')
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

def save_yearly(definition):
    yearlyout = Dataset('threshold_%s_%spcntl.nc'%(options.bp,str(pcntl).rstrip('0').rstrip('.')),'w')
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
    setattr(yearlyout, "definition", definition)
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
    otime = yearlyout.createVariable('time', 'f8', 'time', fill_value=-999.99)
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
    otxpct = yearlyout.createVariable('txpct', 'f8', 
	    ('day','lat','lon'), fill_value=-999.99)
    setattr(otxpct, 'long_name', 'tmax percentile')
    setattr(otxpct, 'units', 'degC')
    setattr(otxpct, 'description', 
            '90th percentile of %s-%s'%(str(bpstart),str(bpend)))
    otnpct = yearlyout.createVariable('tnpct', 'f8', 
	    ('day','lat','lon'), fill_value=-999.99)
    setattr(otnpct, 'long_name', 'tmin percentile')
    setattr(otnpct, 'units', 'degC')
    setattr(otnpct, 'description', 
            '90th percentile of %s-%s'%(str(bpstart),str(bpend)))
    otime[:] = range(first_year, daylast.year+1)
    olat[:] = tmaxnc.variables[latname][:]
    olon[:] = tmaxnc.variables[lonname][:]
    dummy_array = np.ones((daysinyear,)+original_shape[1:])*np.nan
    if options.maskfile:
        dummy_array[:,mask] = txpct
        dummy_array[np.isnan(dummy_array)] = -999.99
        otxpct[:] = dummy_array.copy()
        dummy_array[:,mask] = tnpct
        dummy_array[np.isnan(dummy_array)] = -999.99
        otnpct[:] = dummy_array.copy()
    else:
        otxpct[:] = txpct
	otnpct[:] = tnpct
    yearlyout.close()

if yearlyout:
    save_yearly("txn90pct")
