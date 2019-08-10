from pylab import *
import numpy
import datetime

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
        self2.fd13      = file[:,4]
        self2.fd18      = file[:,6]

# Create object per MXL run
run1 = readdata('hf2009-03-19092011_Fluxes_arizona.csv')

#plotting figures
figure(1, figsize=(12,10))
subplot(221)
plot (run1.time_data, run1.fco2, 'g^', label= 'Fco2')
leg1 = legend(loc=2)
leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('Fco2')

subplot(222)
plot (run1.time_data, run1.fd13, 'b^', label= 'Fd13')
leg1 = legend(loc=2)
leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('Fd13')

subplot(223)
plot (run1.time_data, run1.fd18, 'r^', label= 'Fd18')
leg1 = legend(loc=2)
leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('Fd18')

subplot(224)
plot (run1.time_data, run1.fd18, 'r^', label= 'Fd18')
leg1 = legend(loc=2)
leg1.draw_frame(False)
xlabel('LT [hours]')
ylabel('Fd18')

#show()
