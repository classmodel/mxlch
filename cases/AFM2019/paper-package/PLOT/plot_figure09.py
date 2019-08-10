from pylab import *
import numpy
import datetime
import matplotlib.pylab as pl

# opening files MXLCH
isot      = loadtxt('../AFM19-PAPER/iso_delta',  skiprows=3)
co2var    = loadtxt('../AFM19-PAPER/output_ags', skiprows=3)
outiso    = loadtxt('../AFM19-PAPER/output_sca',  skiprows=3)
isob      = loadtxt('../AFM19-PAPER/iso_budget', skiprows=3)
frac      = loadtxt('../AFM19-PAPER/fraction', skiprows=3)

tm        = isot[:,0] - 5.0
deltac13  = isot[:,2]
deltac18  = isot[:,3]
deltaw18  = isot[:,4]
c13delta  = isot[:,5]
deltax18  = isot[:,8]
deltae18  = isot[:,9]
deltae18t = isot[:,10]
delta_lw  = frac[:,6]
rhsurf    = isot[:,11]
tcan      = isot[:,12]
wca       = co2var[:,2]
wcr       = co2var[:,3]
co2       = outiso[:,8]
wc18      = isob[:,7]
c18_p     = isob[:,8]
c18_r     = isob[:,9]
c18_c     = isob[:,10]
c18_d     = isob[:,11]
wc18s     = isob[:,17]
c18_ps    = isob[:,18]

nee     = wca + wcr

zeroline      = numpy.zeros(720)

rs18     = 0.0020052
need13   = deltac13*(co2/1000.)*40.9   # flux units permil micromol/m2s 
need18   = deltac18*(co2/1000.)*40.9   # flux units permil micromol/m2s 
wc18pm   = wc18*(co2/1000.)*40.9       # flux units permil micromol/m2s 
wc18pm_p = c18_p*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc18pm_r = c18_r*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc18pm_c = c18_c*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc18pm_d = c18_d*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc18pm_s = wc18s*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc18pm_ps= c18_ps*(co2/1000.)*40.9     # flux units permil micromol/m2s 

# From Wehr

wneed13 = (c13delta/1000.) * nee * (co2/1000.)*40.9

# conversion factor from mgC/m2S (MXLCH) to observations (micromol/m2s)
mwco2    = 44.
conv_fc  = (1000/mwco2)

neemol   = conv_fc * nee              # flux units CO2 microm/m2 s

def convert_datetime(value):
    time = datetime.datetime.strptime(value.decode(), '%d/%m/%Y_%H:%M:%S')
    return time.hour + time.minute/60. + time.second/3600. 

# Reading Harvard data .csv files
class readdata:
    def __init__(self2,path):
        file = numpy.loadtxt(path,delimiter=',', skiprows=2, converters={0:convert_datetime})
        self2.time_data = file[:,0]
        self2.fco2      = file[:,2]
        self2.fco2s     = file[:,3]
        self2.fd13      = file[:,4]
        self2.fd13s     = file[:,5]
        self2.fd18      = file[:,6]
        self2.fd18s     = file[:,7]

# Create object per MXL run
run1 = readdata('DATAHAR/ISOTOPES/hf2009-03-16092011_Fluxes_arizona.csv')
run2 = readdata('DATAHAR/ISOTOPES/hf2009-03-17092011_Fluxes_arizona.csv')
run3 = readdata('DATAHAR/ISOTOPES/hf2009-03-18092011_Fluxes_arizona.csv')
run4 = readdata('DATAHAR/ISOTOPES/hf2009-03-19092011_Fluxes_arizona.csv')

# Reading Harvard data .csv files for Relative Humidity surface
class readdata:
	  def __init__(self,path):
	    file = numpy.loadtxt(path,delimiter=',',skiprows=1)
	    self.time_sim  = file[:,3]
	    self.rh       = file[:,22]

# Open data file
arun1 = readdata('DATAHAR/hf004-01-160911-final01.csv')
arun2 = readdata('DATAHAR/hf004-01-170911-final01.csv')
arun3 = readdata('DATAHAR/hf004-01-180911-final01.csv')
arun4 = readdata('DATAHAR/hf004-01-190911-final01.csv')

def convert_datetime(value):
    time = datetime.datetime.strptime(value.decode(), '%Y-%m-%dT%H:%M')
    return time.hour + time.minute/60.

avrhsurf = np.mean(np.array([arun1.rh  ,arun2.rh  ,arun3.rh  ,arun4.rh  ]),axis=0)
strhsurf = np.std (np.array([arun1.rh  ,arun2.rh  ,arun3.rh  ,arun4.rh  ]),axis=0)

# data in hours (local time)
lttime = 0.
hour   = 24*(arun1.time_sim - arun1.time_sim[0]) + lttime

# Average array
# CO2 
avfco2 = np.mean(np.array([run1.fco2,run2.fco2,run3.fco2,run4.fco2]),axis=0)
stfco2 = np.std (np.array([run1.fco2,run2.fco2,run3.fco2,run4.fco2]),axis=0)

nee01  = run1.fco2 + run1.fco2s
nee02  = run2.fco2 + run2.fco2s
nee03  = run3.fco2 + run3.fco2s
nee04  = run4.fco2 + run4.fco2s
avnee = np.mean(np.array([nee01,nee02,nee03,nee04]),axis=0)
stnee = np.std (np.array([nee01,nee02,nee03,nee04]),axis=0)

avfco2s = np.mean(np.array([run1.fco2s,run2.fco2s,run3.fco2s,run4.fco2s]),axis=0)
stfco2s = np.std (np.array([run1.fco2s,run2.fco2s,run3.fco2s,run4.fco2s]),axis=0)
# COO18 
avfd18 = np.mean(np.array([run1.fd18,run2.fd18,run3.fd18,run4.fd18]),axis=0)
stfd18 = np.std (np.array([run1.fd18,run2.fd18,run3.fd18,run4.fd18]),axis=0)

d1801  = run1.fd18 + run1.fd18s
d1802  = run2.fd18 + run2.fd18s
d1803  = run3.fd18 + run3.fd18s
d1804  = run4.fd18 + run4.fd18s
avd18 = np.mean(np.array([d1801,d1802,d1803,d1804]),axis=0)
std18 = np.std (np.array([d1801,d1802,d1803,d1804]),axis=0)

avfd18s = np.mean(np.array([run1.fd18s,run2.fd18s,run3.fd18s,run4.fd18s]),axis=0)
stfd18s = np.std (np.array([run1.fd18s,run2.fd18s,run3.fd18s,run4.fd18s]),axis=0)

#plotting figures
fig1=figure(2, figsize=(10,10))
fig1.subplots_adjust(0.09,0.08,0.91,0.94,0.02,0.02)
#ax = fig1.add_subplot(311)
sub1=subplot(311)
sub1.xaxis.tick_top()
sub1.set_xlabel('LT [hours]')
sub1.xaxis.set_label_position('top')
sub1.grid()
sub1.plot (tm, 100.*rhsurf,'b',   label= 'RH$_{surf}$')
sub1.plot (hour, avrhsurf,'b^',   label= '')
sub1.fill_between(hour, avrhsurf-strhsurf, avrhsurf+strhsurf, color= 'b', label= '', alpha=0.1)
sub1.leg1 = legend(loc=1)
sub1.leg1.draw_frame(False)
sub1.set_ylabel('[%]')
#axis([6.,17.,0.,100.])
ax2 = sub1.twinx()
ax2.plot(tm, tcan, 'k',   label= 'T$_c$')
ax2.set_ylabel('T$_c$ [K]')
ax2.set_ylim(282,293.999)
ax2.leg1 = legend(loc=4)
ax2.leg1.draw_frame(False)
sub1.axis([7.,16.,00.,100.])
annotate('a',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

sub2 = subplot(312)
sub2.yaxis.tick_left()
#sub2.yaxis.tick_right()
sub2.yaxis.set_label_position('left')
sub2.grid()
sub2.plot (tm, deltae18,    'k',    label= '$\delta_{e}^{18}$ Eq. A25')
sub2.plot (tm, delta_lw,    'b--',  label= '$\delta_{lw}$: 1st r.h.s Eq. A25')
sub2.plot (tm, deltae18t,   'k--',  label= '2nd r.h.s Eq. A25')
sub2.leg1 = legend(loc=3, ncol=3)
sub2.leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('[permil]')
setp(sub2.get_xticklabels(), visible=False)
annotate('b',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')
axis([7.,16.,-20.,59.99])

sub3=subplot(313)
sub3.yaxis.tick_left()
sub3.yaxis.set_label_position('left')
sub3.grid()
plot (tm, deltax18,   'k-',  label= '$\delta_{eff}$ Eq. B.2')
leg1 = legend(loc=3)
leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('[permil]')
annotate('c',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')
axis([7.,16.,00.,59.99])

savefig('figure08.png',dpi=300)
#show()
