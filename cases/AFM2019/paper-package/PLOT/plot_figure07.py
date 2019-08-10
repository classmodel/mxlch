from pylab import *
import numpy
import datetime
import matplotlib.pylab as pl

# opening files MXLCH (original file 081217 before chnage -22 to -28)
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

nee     = wca + wcr

zeroline      = numpy.zeros(720)

need13   = deltac13*(co2/1000.)*40.9   # flux units permil micromol/m2s 
need18   = deltac18*(co2/1000.)*40.9   # flux units permil micromol/m2s 
wc13pm   = wc13*(co2/1000.)*40.9       # flux units permil micromol/m2s 
wc13pm   = wc13*(co2/1000.)*40.9       # flux units permil micromol/m2s 
wc13pm_p = c13_p*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc13pm_r = c13_r*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc13pm_c = c13_c*(co2/1000.)*40.9      # flux units permil micromol/m2s 
wc13pm_d = c13_d*(co2/1000.)*40.9      # flux units permil micromol/m2s 

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

#plotting figures
fig1=figure(2, figsize=(9,10))
fig1.subplots_adjust(0.09,0.08,0.91,0.94,0.02,0.02)
sub1=subplot(211)
sub1.xaxis.tick_top()
sub1.set_xlabel('LT [hours]')
sub1.xaxis.set_label_position('top')
sub1.grid()
sub1.plot (run1.time_data, avfd13,             'go', label= 'EC')
sub1.fill_between(run1.time_data, avfd13-stfd13, avfd13+stfd13, color= 'g', label= '', alpha=0.2)
sub1.plot (run1.time_data, avfd13s,         'ro', label= 'STO')
sub1.fill_between(run1.time_data, avfd13s-stfd13s, avfd13s+stfd13s, color= 'r', label= '', alpha=0.1)
sub1.plot (run1.time_data, avd13,           'bo', label= 'NEE')
sub1.fill_between(run1.time_data, avd13-std13, avd13+std13, color= 'b', label= '', alpha=0.1)
#pl.plot (run1.time_data, run1.fd13,  'go', label= 'Fd13_e')
#pl.plot (run1.time_data, run1.fd13s, 'ro', label= 'Fd13_s')
#pl.plot (run1.time_data, run1.fd13 + run1.fd13s, 'ko', label= 'NEEd13')
sub1.plot (tm, need13, 'b', linewidth=1, label= 'NEE MXL')
sub1.plot(tm, zeroline, color='k',ls='dashed')
#plot (tm, wneed13, 'b^', label= '')
sub1.leg1 = legend(loc=1)
sub1.leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('F $^{13}$CO$_2$ [permil $\mu$mol m$^{-2}$ s$^{-1}$]')
#pl.ylabel('FD13 [permil micromol /m2 s]')
axis([7.,16.,-50.,600.])
annotate('a',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

sub2=subplot(212)
sub2.grid()
#sub2.plot (tm, wc13pm,   'k', label= 'NEEd13')
sub2.plot (tm, wc13pm_p, 'g', label= 'PLA')
sub2.plot (tm, wc13pm_r, 'r', label= 'RES')
sub2.plot(tm, zeroline, color='k',ls='dashed')
sub2.leg1 = legend(loc=1)
sub2.leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('F $^{13}$CO$_2$ [permil $\mu$mol m$^{-2}$ s$^{-1}$]')
axis([7.,16.,-200.,600.])
annotate('b',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')
#axis([0.,24.,-500.,500.])

#savefig('fig_flc13_paper.jpg',dpi=300)
savefig('figure07.png',dpi=300)
#show()
