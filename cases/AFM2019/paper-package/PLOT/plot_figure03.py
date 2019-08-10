 # Simple plotting model FORTRAN output

from pylab import *
import copy
import numpy
import datetime

# opening files MXLCH
outdyn        = loadtxt('../AFM19-PAPER/output_dyn', skiprows=3)
outsca        = loadtxt('../AFM19-PAPER/output_sca', skiprows=3)

# Reading Harvard data .csv files
class readdata:
	  def __init__(self,path):
            file = numpy.loadtxt(path,delimiter=',',skiprows=1)
            self.time_sim  = file[:,3]
            self.theta     = file[:,27]
            self.rh        = file[:,22]
            self.co2       = file[:,35]

# Open data file 
run1 = readdata('DATAHAR/hf004-01-160911-final01.csv')
run2 = readdata('DATAHAR/hf004-01-170911-final01.csv')
run3 = readdata('DATAHAR/hf004-01-180911-final01.csv')
run4 = readdata('DATAHAR/hf004-01-190911-final01.csv')

def convert_datetime(value):
    time2= datetime.datetime.strptime(value.decode(), '%Y-%m-%dT%H:%M:%S')
    return time2.hour + time2.minute/60. + time2.second/3600.

# Reading Harvard data .csv files
class readdata2:
    def __init__(self2,path):
        file = numpy.loadtxt(path,delimiter=',', skiprows=0, converters={0:convert_datetime})
        self2.time2_data = file[:,0]
        self2.co2w       = file[:,2]

# Create object per MXL run
run5 = readdata2('DATAHAR/ISOTOPES/hf209-03-19092011.csv')

# Defining constants

g     = 9.81    #[m/s2]
epsilon = 0.622   #[g water / g air]
tabs  = 273.    #[K]
cp    = 1004.67 #[J/kg K]
Rd    = 287.053 #[J/kg K]
Rv    = 461.50  #[J/kg K]
Lv    = 2.5e6   #[J/kg]
e0    = 0.611   #[kPa]
gd    = -(g/cp) #[K/m]
gamma_dry = g/cp
mwh20     = 18.
rhoair    = 1.225
rholv     = 3013.5
mwco2    = 44.

height28 = 27.9

timeabl  =  15.
ablheight= 2045.
# Calculating potential tempertaure
theta28m01 = run1.theta + (g/cp) * height28 + tabs
theta28m02 = run2.theta + (g/cp) * height28 + tabs
theta28m03 = run3.theta + (g/cp) * height28 + tabs
theta28m04 = run4.theta + (g/cp) * height28 + tabs

# Calculating average
avtheta = np.mean(np.array([theta28m01, theta28m02, theta28m03, theta28m04]),axis=0)
sttheta = np.std (np.array([theta28m01, theta28m02, theta28m03, theta28m04]),axis=0)

# Calculating specific humidity
pressure = 992.                                           #specific 19 Sept 2017
es01       = e0*exp((Lv/Rv)*((1/tabs)-(1/(run1.theta+tabs))))
qs01       = (epsilon*es01)/pressure
e01        = (run1.rh/100)*es01
qt28m01    = (epsilon*e01)/(pressure/10)

es02       = e0*exp((Lv/Rv)*((1/tabs)-(1/(run2.theta+tabs))))
qs02       = (epsilon*es02)/pressure
e02        = (run2.rh/100)*es02
qt28m02    = (epsilon*e02)/(pressure/10)

es03       = e0*exp((Lv/Rv)*((1/tabs)-(1/(run3.theta+tabs))))
qs03       = (epsilon*es03)/pressure
e03        = (run2.rh/100)*es03
qt28m03    = (epsilon*e03)/(pressure/10)

es04       = e0*exp((Lv/Rv)*((1/tabs)-(1/(run4.theta+tabs))))
qs04       = (epsilon*es04)/pressure
e04        = (run4.rh/100)*es03
qt28m04    = (epsilon*e04)/(pressure/10)

# Calculating average
avqt = np.mean(np.array([qt28m01, qt28m02, qt28m03, qt28m04]),axis=0)
stqt = np.std (np.array([qt28m01, qt28m02, qt28m03, qt28m04]),axis=0)


# defining arrays
tmfor       = outdyn[:,0]
zi          = outdyn[:,2]
theta       = outdyn[:,4]
q           = outsca[:,3]

# data in hours (local time)

lttime = 5.0
tmfor  = tmfor - lttime 
hour   = 24*(run1.time_sim - run1.time_sim[0]) 


# plotting chemical species
fig1=figure(2, figsize=(8,10))
fig1.subplots_adjust(0.09,0.08,0.91,0.94,0.02,0.02)
#ax = fig1.add_subplot(311)
sub1=subplot(311)
#ax.get_xaxis().tick_bottom()
sub1.xaxis.tick_top()
sub1.set_xlabel('LT [hours]')
sub1.xaxis.set_label_position('top')
sub1.grid()
plot(tmfor, theta, 'k')
#plot (hour, theta28m01, 'ro')
#plot (hour, theta28m02, 'bo')
#plot (hour, theta28m03, 'ko')
#plot (hour, theta28m04, 'go')
plot (hour, avtheta,             'ko')
fill_between(hour, avtheta-sttheta, avtheta+sttheta, color= 'k', label= '', alpha=0.2)
ylabel('$\Theta$ [K]')
#leg1 = legend(loc=2)
#leg1.draw_frame(False)
annotate('a',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')
axis([7.,16.,280.,292.])

sub2 = subplot(312)
sub2.yaxis.tick_left()
#sub2.yaxis.tick_right()
sub2.yaxis.set_label_position('left')
sub2.grid()
plot(tmfor, q, 'k-', label= '')
#plot (hour, qt28m01*1000., 'bo')
#plot (hour, qt28m02*1000., 'ro')
#plot (hour, qt28m03*1000., 'go')
#plot (hour, qt28m04*1000., 'ko')
plot (hour, avqt*1000.,             'ko', label= '')
fill_between(hour, 1000.*(avqt-stqt), 1000.*(avqt+stqt), color= 'k', label= '', alpha=0.2)
xlabel('Time (h UTC)')
ylabel('q [g$_w$ kg$_a$$^{-1}$]')
setp(sub2.get_xticklabels(), visible=False)
#leg1 = legend(loc=1)
#leg1.draw_frame(False)
annotate('b',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')
axis([7.,16.,4.,7.49999])

sub3=subplot(313)
sub3.yaxis.tick_left()
sub3.yaxis.set_label_position('left')
sub3.grid()
plot(tmfor, zi, 'k')
plot(timeabl, ablheight, 'ko')
errorbar(timeabl, ablheight, yerr=100., color='k')
#plot(tmfor, gf, 'r')
xlabel('LT [hours]')
ylabel('h [m]')
text(14.6,1750,'19 LT',fontsize=14)
annotate('c',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')
axis([7.,16.,00.,2200.])

savefig('figure03.png',dpi=300)
#savefig('fig_land.eps',dpi=300)
