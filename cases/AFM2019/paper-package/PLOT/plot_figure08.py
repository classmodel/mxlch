from pylab import *
import numpy
import datetime
import matplotlib.pylab as pl

# opening files MXLCH (original before change -22 to -28 file har240218)
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
sub1=subplot(211)
#ax.get_xaxis().tick_bottom()
sub1.xaxis.tick_top()
sub1.set_xlabel('LT [hours]')
sub1.xaxis.set_label_position('top')
sub1.grid()
sub1.plot (run1.time_data, avfd18,             'go', label= 'EC')
sub1.fill_between(run1.time_data, avfd18-stfd18, avfd18+stfd18, color= 'g', label= '', alpha=0.2)
sub1.plot (run1.time_data, avfd18s,         'ro', label= 'ST')
sub1.fill_between(run1.time_data, avfd18s-stfd18s, avfd18s+stfd18s, color= 'r', label= '', alpha=0.1)
sub1.plot (run1.time_data, avd18,           'bo', label= 'NEE')
sub1.fill_between(run1.time_data, avd18-std18, avd18+std18, color= 'b', label= '', alpha=0.1)
#plot (run1.time_data, run1.fd18,  'go', label= 'Fd18_e')
#plot (run1.time_data, run1.fd18s, 'ro', label= 'Fd18_s')
#plot (run1.time_data, run1.fd18 + run1.fd18s, 'ko', label= 'NEEd18')
sub1.plot (tm, need18,   'b',   linewidth=1, label= 'NEE PAR1')
sub1.plot (tm, wc18pm_s, 'b--', linewidth=1, label= 'NEE PAR2')
sub1.plot(tm, zeroline, color='k',linestyle='.')
#plot (tm, wc18pm, 'b^', label= '')
#plot (tm, wneed13, 'b^', label= '')
sub1.leg1 = legend(loc=4)
sub1.leg1.draw_frame(False)
pl.xlabel('LT [hours]')
pl.ylabel('F C$^{18}$OO  [permil $\mu$mol m$^{-2}$s$^{-1}$]')
pl.axis([7.,16.,-600.,400.])
annotate('a',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

sub2=subplot(212)
sub2.grid()
#plot (tm, wc18pm,   'k', label= 'NEEd18 (Eq. 22)')
#plot (tm, wc18pm_s, 'k--', label= 'NEEd18 (Eq. 34)')
plot (tm, wc18pm_p, 'g', label= 'PLANT PAR1')
plot (tm, wc18pm_ps, 'g--', label= 'PLANT PAR2')
plot (tm, wc18pm_r, 'r', label= 'SOIL Eq. A32')
plot(tm, zeroline, 'k-')
leg1 = legend(loc=2)
leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('F C$^{18}$OO  [permil $\mu$mol m$^{-2}$s$^{-1}$]')
#axis([6.,17.,-8000.,6000.])
axis([7.,16.,-200.,399.9999])
annotate('b',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

savefig('figure08.png',dpi=300)
#savefig('fig_flc18_paper.png',dpi=300)
#show()
