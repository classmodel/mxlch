from pylab import *
import numpy
import datetime
import matplotlib.pylab as pl

# opening files MXLCH
isot      = loadtxt('../AFM19-PAPER/iso_delta',  skiprows=3)
co2var    = loadtxt('../AFM19-PAPER/output_ags', skiprows=3)
outiso    = loadtxt('../AFM19-PAPER/output_sca',  skiprows=3)
isob      = loadtxt('../AFM19-PAPER/iso_budget', skiprows=3)

tm        = isot[:,0] - 5.0
deltac13  = isot[:,2]
deltac18  = isot[:,3]
deltaw18  = isot[:,4]
c13delta  = isot[:,5]
wca       = co2var[:,2]
wcr       = co2var[:,3]
co2       = outiso[:,8]
wc13      = isob[:,2]
c13_p     = isob[:,3]
c13_r     = isob[:,4]
c13_c     = isob[:,5]
c13_d     = isob[:,6]
wc18s     = isob[:,17]

nee       = wca + wcr

zeroline  = numpy.zeros(720)

need13   = deltac13*(co2/1000.)*40.9   # flux units permil micromol/m2s 
need18   = deltac18*(co2/1000.)*40.9   # flux units permil micromol/m2s 
wc13pm   = wc13*(co2/1000.)*40.9       # flux units permil micromol/m2s 
wc13pm   = wc13*(co2/1000.)*40.9       # flux units permil micromol/m2s 
wc13pm_p = c13_p*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc13pm_r = c13_r*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc13pm_c = c13_c*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc13pm_d = c13_d*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc18pm_s = wc18s*(co2/1000.)*40.9      # flux units permil micromol/m2s

# From Wehr

wneed13 = (c13delta/1000.) * nee * (co2/1000.)*40.9

# conversion factor from mgC/m2S (MXLCH) to observations (micromol/m2s)
mwco2    = 44.
conv_fc  = (1000/mwco2)

neemol   = conv_fc * nee              # flux units CO2 microm/m2 s

def convert_datetime(value):
#    time = datetime.datetime.strptime(value.decode(), '%Y-%m-%dT%H:%M:%S')
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

run1 = readdata('DATAHAR/ISOTOPES/hf2009-03-16092011_Fluxes_arizona.csv')
run2 = readdata('DATAHAR/ISOTOPES/hf2009-03-17092011_Fluxes_arizona.csv')
run3 = readdata('DATAHAR/ISOTOPES/hf2009-03-18092011_Fluxes_arizona.csv')
run4 = readdata('DATAHAR/ISOTOPES/hf2009-03-19092011_Fluxes_arizona.csv')

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
# COO13
avfd13 = np.mean(np.array([run1.fd13,run2.fd13,run3.fd13,run4.fd13]),axis=0)
stfd13 = np.std (np.array([run1.fd13,run2.fd13,run3.fd13,run4.fd13]),axis=0)

d1301  = run1.fd13 + run1.fd13s
d1302  = run2.fd13 + run2.fd13s
d1303  = run3.fd13 + run3.fd13s
d1304  = run4.fd13 + run4.fd13s
avd13 = np.mean(np.array([d1301,d1302,d1303,d1304]),axis=0)
std13 = np.std (np.array([d1301,d1302,d1303,d1304]),axis=0)

avfd13s = np.mean(np.array([run1.fd13s,run2.fd13s,run3.fd13s,run4.fd13s]),axis=0)
stfd13s = np.std (np.array([run1.fd13s,run2.fd13s,run3.fd13s,run4.fd13s]),axis=0)

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

# new data set
def convert_datetime(value):
    time = datetime.datetime.strptime(value.decode(), '%d/%m/%Y_%H:%M:%S')
    return time.hour + time.minute/60.

class readdata2:
          def __init__(self,path):
            file = numpy.loadtxt(path,delimiter=',', skiprows=1, converters={0:convert_datetime})
	    self.time_data = file[:,0]
#	    self.neehar    = file[:,11]
	    self.gpp       = -1.* file[:,5]
	    self.res       = file[:,6]

har1 = readdata2('DATAHAR/ISOTOPES/GPPRES/HF-PartitionedFluxes-180911.csv')
har2 = readdata2('DATAHAR/ISOTOPES/GPPRES/HF-PartitionedFluxes-190911.csv')

#plotting figures
fig1=figure(2, figsize=(9,10))
fig1.subplots_adjust(0.09,0.08,0.91,0.94,0.02,0.02)
#ax = fig1.add_subplot(311)
sub1=subplot(211)
#ax.get_xaxis().tick_bottom()
sub1.xaxis.tick_top()
sub1.set_xlabel('LT [hours]')
sub1.xaxis.set_label_position('top')
sub1.grid()
plot (run1.time_data, avfco2,             'go', label= 'EC')
fill_between(run1.time_data, avfco2-stfco2, avfco2+stfco2, color= 'g', label= '', alpha=0.2)
plot (run1.time_data, avfco2s,         'ro', label= 'ST')
fill_between(run1.time_data, avfco2s-stfco2s, avfco2s+stfco2s, color= 'r', label= '', alpha=0.1)
plot (run1.time_data, avnee,           'bo', label= 'NEE')
fill_between(run1.time_data, avnee-stnee, avnee+stnee, color= 'b', label= '', alpha=0.1)
#pl.plot (run1.time_data, run1.fco2,  'go', label= 'Fco2_e')
#pl.plot (run1.time_data, run1.fco2s, 'ro', label= 'Fco2_s')
#pl.plot (run1.time_data, run1.fco2 + run1.fco2s, 'ko', label= 'NEE')
plot (tm, neemol, 'b', label= 'NEE MXL', linewidth=1)
plot(tm, zeroline, color='k',ls='dashed')
leg1 = legend(loc=3)
leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('F CO$_2$ [$\mu$mol m$^{-2}$s$^{-1}$]')
axis([7.,16.,-20.,5.])
annotate('a',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

sub2=subplot(212)
sub2.grid()
sub2.plot(tm, conv_fc*wcr, 'r-', label='RES')
#plot (run2.ime_data + lttime , run2.res, 'ro')
plot(tm, conv_fc*wca, 'g-', label= 'GEP')
plot(tm,zeroline, color='k',ls='dashed')
plot (har1.time_data, har1.res,           'ro', label= 'RES-18/09')
plot (har1.time_data, har1.gpp,           'go', label= 'GEP-18/09')
plot (har2.time_data, har2.res,           'r^', label= 'RES-19/09')
plot (har2.time_data, har2.gpp,           'g^', label= 'GEP-19/09')
xlabel('LT [hours]')
ylabel('F CO$_2$  [$\mu$mol m$^{-2}$s$^{-1}$]')
leg1 = legend(loc=10, ncol=3)
leg1.draw_frame(False)
axis([7.,16.,-28.,8.])
annotate('b',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

savefig('figure05.png',dpi=300)
#show()
