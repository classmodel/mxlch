from pylab import *
import numpy
import datetime
import matplotlib.pylab as pl
from matplotlib import rc 

# opening files MXLCH
conc12    = loadtxt('../AFM19-PAPER/output_sca', skiprows=3)
conc13    = loadtxt('../AFM19-PAPER/C13', skiprows=3)
conc18    = loadtxt('../AFM19-PAPER/COO18', skiprows=3)
isot      = loadtxt('../AFM19-PAPER/iso_delta', skiprows=3)

tm        = conc12[:,0] - 5.0
c12       = conc12[:,8]/1000.   # ppm
c13       = conc13[:,2]/1000.   # ppm
c18       = conc18[:,2]/1000.   # ppm
deltac13  = isot[:,5]
deltac18  = isot[:,6]
deltaw18  = isot[:,7]

def convert_datetime(value):
    time = datetime.datetime.strptime(value.decode(), '%Y-%m-%dT%H:%M:%S')
    return time.hour + time.minute/60. + time.second/3600. 

# Reading Harvard data .csv files
class readdata:
    def __init__(self2,path):
        file = numpy.loadtxt(path,delimiter=',', skiprows=1, converters={0:convert_datetime})
        self2.time_data = file[:,0]
        self2.co2       = file[:,2]
        self2.d13       = file[:,3]
        self2.d18       = file[:,4]
        self2.pres      = file[:,9]

# Create object per MXL run
run1 = readdata('DATAHAR/ISOTOPES/hf209-03-ambient_co2-160911.csv')
run2 = readdata('DATAHAR/ISOTOPES/hf209-03-ambient_co2-170911.csv')
run3 = readdata('DATAHAR/ISOTOPES/hf209-03-ambient_co2-180911.csv')
run4 = readdata('DATAHAR/ISOTOPES/hf209-03-ambient_co2-190911.csv')

# Calculating averages and standard deviations
# CO2
avco2 = np.mean(np.array([run1.co2,run2.co2,run3.co2,run4.co2]),axis=0)
stco2 = np.std (np.array([run1.co2,run2.co2,run3.co2,run4.co2]),axis=0)
# delta13 
avd13 = np.mean(np.array([run1.d13,run2.d13,run3.d13,run4.d13]),axis=0)
std13 = np.std (np.array([run1.d13,run2.d13,run3.d13,run4.d13]),axis=0)
# delta18 
avd18 = np.mean(np.array([run1.d18,run2.d18,run3.d18,run4.d18]),axis=0)
std18 = np.std (np.array([run1.d18,run2.d18,run3.d18,run4.d18]),axis=0)

#
Rsc13       = 0.0110570   # VPDB Adapted from Wehr's paper 0.0111797 
Rsc18       = 0.0020052   # VSMOW

# Calculating the stable isotope mixing ratio
obsc13run1 = run1.co2*Rsc13 * (run1.d13/1000. + 1)
obsc13run2 = run2.co2*Rsc13 * (run2.d13/1000. + 1)
obsc13run3 = run3.co2*Rsc13 * (run3.d13/1000. + 1)
obsc13run4 = run4.co2*Rsc13 * (run4.d13/1000. + 1)

avc13 = np.mean(np.array([obsc13run1,obsc13run2,obsc13run3,obsc13run4]),axis=0)
stc13 = np.std (np.array([obsc13run1,obsc13run2,obsc13run3,obsc13run4]),axis=0)

obsc18run1 = 2.*run1.co2*Rsc18 * (run1.d18/1000. + 1) 
obsc18run2 = 2.*run2.co2*Rsc18 * (run2.d18/1000. + 1) 
obsc18run3 = 2.*run3.co2*Rsc18 * (run3.d18/1000. + 1) 
obsc18run4 = 2.*run4.co2*Rsc18 * (run4.d18/1000. + 1) 

avc18 = np.mean(np.array([obsc18run1,obsc18run2,obsc18run3,obsc18run4]),axis=0)
stc18 = np.std (np.array([obsc18run1,obsc18run2,obsc18run3,obsc18run4]),axis=0)

# c13 initial
d1300  = -7.80
d1800  = 38.8 
c1200  = 388000.

c1300  = c1200*Rsc13 * (d1300/1000. + 1.)
c1800  = 2.*c1200*Rsc18 * (d1800/1000. + 1.)

print 'c1300 =', c1300
print 'c1800 =', c1800

# Calculating the stable isotope mixing ratio
mxld13 = ((1/Rsc13)*(c13/c12) -1) * 1000.

#plotting figures
fig1=figure(2, figsize=(10,10))
fig1.subplots_adjust(0.09,0.08,0.91,0.94,0.02,0.02)
#ax = fig1.add_subplot(311)
sub1=subplot(311)
#ax.get_xaxis().tick_bottom()
sub1.xaxis.tick_top()
sub1.set_xlabel('LT [hours]')
sub1.xaxis.set_label_position('top')
sub1.grid()
#sub1.axis["bottom"].set_visible(True)
#sub1.set_xticklabels(column_labels, minor=False)
#sub1.xaxis.tick_bottom('0')
plot (run1.time_data, avco2,             'ko', label= '')
fill_between(run1.time_data, avco2-stco2, avco2+stco2, color= 'k', label= '', alpha=0.2)
plot (tm, c12, 'k', label= '', linewidth=1)
#setp(sub1.get_xticklabels(), visible=False)
#pl.leg1 = legend(loc=1)
#pl.leg1.draw_frame(False)
#xlabel('LT [hours]')
annotate('a',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')
ylabel('CO$_2$ [ppm]')
axis([7.,16.,372.,392.])

sub2 = subplot(312)
sub2.yaxis.tick_left()
#sub2.yaxis.tick_right()
sub2.yaxis.set_label_position('left')
sub2.grid()
plot (run1.time_data, avd13,             'ko', label= '')
fill_between(run1.time_data, avd13-std13, avd13+std13, color= 'k', label= '', alpha=0.2)
#pl.plot (run1.time_data, run1.d13, 'go', label= '')
#pl.plot (run2.time_data, run2.d13, 'ro', label= '')
#pl.plot (run3.time_data, run3.d13, 'bo', label= '')
#pl.plot (run4.time_data, run4.d13, 'ko', label= '')
#plot (tm, mxld13, 'g', label= '')
plot (tm, deltac13  , 'k'  , linewidth=1)
#plot (tm, deltac13  , 'k'  , label= 'CTL', linewidth=1)
#plot (tm, deltac1301, 'k--', label= 'ENT', linewidth=1)
#pl.leg1 = legend(loc=4)
#pl.leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('$\delta^{13}_a$ [permil]')
axis([7.,16.,-8.4,-7.2])
annotate('b',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom', fontsize=18)


sub3=subplot(313)
sub3.yaxis.tick_left()
sub3.yaxis.set_label_position('left')
sub3.grid()
plot (run1.time_data, avd18,             'ko', label= '')
fill_between(run1.time_data, avd18-std18, avd18+std18, color= 'k', label= '', alpha=0.2)
plot (tm, deltac18, 'k'    ,  linewidth=1)
#plot (tm, deltac18, 'k'    , label= 'CTL', linewidth=1)
#plot (tm, deltac1801, 'k--', label= 'ENT', linewidth=1)
#plot (tm, deltac18par2, 'k-.', label= 'PAR2', linewidth=1)
#pl.leg1 = legend(loc=4)
#pl.leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('$\delta^{18}_a$  [permil]')
axis([7.,16.,37, 39.999])
setp(sub2.get_xticklabels(), visible=False)
annotate('c',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

savefig('figure04.png',dpi=300)
#savefig('fig_isoc13_paper.eps',dpi=300)
#show()
