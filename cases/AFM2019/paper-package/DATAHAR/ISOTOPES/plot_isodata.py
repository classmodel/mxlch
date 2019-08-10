from pylab import *
import numpy
import datetime

def convert_datetime(value):
    time = datetime.datetime.strptime(value.decode(), '%Y-%m-%dT%H:%M:%S')
    return time.hour + time.minute/60. + time.second/3600. 

# Reading Harvard data .csv files
class readdata:
    def __init__(self2,path):
        file = numpy.loadtxt(path,delimiter=',', skiprows=0, converters={0:convert_datetime})
        self2.time_data = file[:,0]
        self2.co2       = file[:,2]
        self2.d13       = file[:,3]
        self2.d18       = file[:,4]
        self2.pres      = file[:,9]

# Create object per MXL run
run1 = readdata('hf209-03-19092011.csv')

#
Rsc13       = 0.0112372   # VPDB
Rsc18       = 0.0020671   # VSMOW

# Calculating the stable isotope mixing ratio
c13 = run1.co2*Rsc13 * (run1.d13/1000. + 1) 
c18 = run1.co2*Rsc18 * (run1.d18/1000. + 1) 

#plotting figures
figure(1, figsize=(12,10))
subplot(221)
plot (run1.time_data, run1.co2, 'b^', label= 'CO2')
leg1 = legend(loc=2)
leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('CO2')

subplot(222)
plot (run1.time_data, run1.d13, 'g^', label= 'd13')
#plot (run1.time_data, run1.d18, 'r^', label= 'd18')
leg1 = legend(loc=2)
leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('d13')
#axis([0.,24.,-10.,40.])

subplot(223)
plot (run1.time_data, c13, 'b^', label= 'c13 (ppm)')
leg1 = legend(loc=2)
leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('d18o')

subplot(224)
plot (run1.time_data, c18, 'g^', label= 'c18 (ppm)')
leg1 = legend(loc=2)
leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('P')

#show()
