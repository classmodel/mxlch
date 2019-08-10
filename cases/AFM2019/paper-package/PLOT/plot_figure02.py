 # Simple plotting model FORTRAN output

from pylab import *
import copy
import numpy
import datetime

# opening files MXLCH
outdyn        = loadtxt('../AFM19-PAPER/output_land', skiprows=3)
outsca        = loadtxt('../AFM19-PAPER/output_sca',  skiprows=3)
co2var        = loadtxt('../AFM19-PAPER/output_ags',  skiprows=3)

# Reading Harvard data .csv files
class readdata:
	  def __init__(self,path):
            file = numpy.loadtxt(path,delimiter=',',skiprows=1)
            self.time_sim  = file[:,3]
            self.qn       = file[:,19]
            self.nee      = file[:,13]
            self.sh       = file[:,15]
            self.lh       = file[:,16]
            self.par12    = file[:,20]
            self.par      = file[:,21]

# Open data file 
run1 = readdata('DATAHAR/hf004-01-160911-final01.csv')
run2 = readdata('DATAHAR/hf004-01-170911-final01.csv')
run3 = readdata('DATAHAR/hf004-01-180911-final01.csv')
run4 = readdata('DATAHAR/hf004-01-190911-final01.csv')

def convert_datetime(value):
    time = datetime.datetime.strptime(value.decode(), '%Y-%m-%dT%H:%M')
    return time.hour + time.minute/60.

# net radiation 
avqn  = np.mean(np.array([run1.qn ,run2.qn ,run3.qn ,run4.qn ]),axis=0)
stqn  = np.std (np.array([run1.qn ,run2.qn ,run3.qn ,run4.qn ]),axis=0)

# sensible heat flux 
avsh  = np.mean(np.array([run1.sh ,run2.sh ,run3.sh ,run4.sh ]),axis=0)
stsh  = np.std (np.array([run1.sh ,run2.sh ,run3.sh ,run4.sh ]),axis=0)

# latent heat flux 
avlh  = np.mean(np.array([run1.lh ,run3.lh ,run4.lh ]),axis=0)
stlh  = np.std (np.array([run1.lh ,run3.lh ,run4.lh ]),axis=0)

# removing the 0 values
index     = [10,11]
newavlh   = np.delete(avlh,index)
newstlh   = np.delete(stlh,index)
# Reading Harvard data .csv files
class readdata2:
          def __init__(self,path):
            file = numpy.loadtxt(path,delimiter=',', skiprows=1, converters={0:convert_datetime})
            self.time_data = file[:,0]
            self.nee       = file[:,11]
            self.res       = file[:,12]
            self.gee       = file[:,13]
            self.partot    = file[:,20]

# Create object per MXL run
run5 = readdata2('DATAHAR/hf004-02-160911.csv')
run6 = readdata2('DATAHAR/hf004-02-170911.csv')
run7 = readdata2('DATAHAR/hf004-02-180911.csv')
run8 = readdata2('DATAHAR/hf004-02-190911.csv')

def convert_datetime(value):
#    time = datetime.datetime.strptime(value.decode(), '%Y-%m-%dT%H:%M:%S')
    time3 = datetime.datetime.strptime(value.decode(), '%d/%m/%Y_%H:%M:%S')
    return time3.hour + time3.minute/60. + time3.second/3600.

# Reading Harvard data .csv files
class readdata3:
	  def __init__(self2,path):
	    file = numpy.loadtxt(path,delimiter=',', skiprows=2, converters={0:convert_datetime})
	    self2.time_data3= file[:,0]
	    self2.fco2      = file[:,2]
	    self2.fco2s     = file[:,3]
	    self2.fd13      = file[:,4]
	    self2.fd13s     = file[:,5]
	    self2.fd18      = file[:,6]
	    self2.fd18s     = file[:,7]

# Create object per MXL run
run9 = readdata3('DATAHAR/ISOTOPES/hf2009-03-19092011_Fluxes_arizona.csv')

# Conversion between observations and model results
# SH and LE in W/m2
# NEE micromol/m2 sec
# conversion factor from milimol/m2s to W/m2 (latent heat flux)

mwh20     = 18.
rhoair    = 1.225
rholv     = 3013.5

conv_le   = (mwh20*rholv)/(1000.*rhoair)

# conversion factor from mgC/m2S (MXLCH) to observations (micromol/m2s)

mwco2    = 44.

conv_fc  = (1000/mwco2)

# equivalnce PAR to w/2

par01     = run1.par/4.6 + 0.1*(run1.par/4.6)   # adding 10% diffuse contribution
par02     = run2.par/4.6 + 0.1*(run2.par/4.6)   # adding 10% diffuse contribution
par03     = run3.par/4.6 + 0.1*(run3.par/4.6)   # adding 10% diffuse contribution
par04     = run4.par/4.6 + 0.1*(run4.par/4.6)   # adding 10% diffuse contribution
# average PAR 
avpar  = np.mean(np.array([par01 ,par02 ,par03 ,par04 ]),axis=0)
stpar  = np.std (np.array([par01 ,par02 ,par03 ,par04 ]),axis=0)

partot01= run5.partot/4.6
partot02= run6.partot/4.6
partot03= run7.partot/4.6
partot04= run8.partot/4.6

# soil respiration 
avres  = np.mean(np.array([run5.res ,run6.res ,run7.res ,run8.res ]),axis=0)
stres  = np.std (np.array([run5.res ,run6.res ,run7.res ,run8.res ]),axis=0)

# soil gee 
avgee  = np.mean(np.array([run5.gee ,run6.gee ,run7.gee ,run8.gee ]),axis=0)
stgee  = np.std (np.array([run5.gee ,run6.gee ,run7.gee ,run8.gee ]),axis=0)

# data in hours (local time)
lttime = 0. 
hour   = 24*(run1.time_sim - run1.time_sim[0]) + lttime

index     = [10,11]
newhour   = np.delete(hour,index)
# defining arrays
tmfor       = outdyn[:,0] - 5.0
swi         = outdyn[:,2]
qn          = outdyn[:,6]
sh          = outdyn[:,7]
le          = outdyn[:,8]
gf          = outdyn[:,9]
wca         = co2var[:,2]
wcr         = co2var[:,3]

nee = wca + wcr

zeroline      = numpy.zeros(24)

# plotting chemical species
fig1=figure(2, figsize=(8,10))
fig1.subplots_adjust(0.09,0.08,0.91,0.94,0.02,0.02)
sub1=subplot(311)
sub1.xaxis.tick_top()
sub1.set_xlabel('LT [hours]')
sub1.xaxis.set_label_position('top')
sub1.grid()
#plot (run5.time_data + lttime, partot, 'ro', label= 'par12')
plot(tmfor, qn, 'k', label= 'Qn')
#plot (hour, run1.qn, 'bo', label= '')
#plot (hour, run2.qn, 'ko', label= '')
#plot (hour, run3.qn, 'go', label= '')
#plot (hour, run4.qn, 'ro', label= '')
plot (hour, avqn,             'ko', label= '')
fill_between(hour, avqn-stqn, avqn+stqn, color= 'k', label= '', alpha=0.2)
#plot(tmfor, gf, 'k', label= 'G')
plot(tmfor, 0.5*swi, 'k--', label= 'PAR')
#plot (hour, par01, 'ro', label= '')
#plot (hour, par02, 'r^', label= '')
#plot (hour, par03, 'g^', label= '')
#plot (hour, par04, 'k^', label= '')
plot (hour, avpar,             'ko', label= '')
fill_between(hour, avpar-stpar, avpar+stpar, color= 'k', label= '', alpha=0.2)
xlabel('LT [hours]')
ylabel('Radiative fluxes [Wm$^{-2}$]')
leg1 = legend(loc=1)
leg1.draw_frame(False)
annotate('a',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')
axis([7.,16.,00.,800.])

sub2 = subplot(312)
sub2.yaxis.tick_left()
#sub2.yaxis.tick_right()
sub2.yaxis.set_label_position('left')
sub2.grid()
plot(tmfor, sh, 'k-', label= '')
#plot (hour, run1.sh, 'r^')
#plot (hour, run2.sh, 'b^')
#plot (hour, run3.sh, 'k^')
#plot (hour, run4.sh, 'g^')
plot (hour, avsh,             'ko', label= '')
fill_between(hour, avsh-stsh, avsh+stsh, color= 'k', label= '', alpha=0.2)
xlabel('LT [hours]')
ylabel('SH [Wm$^{-2}$]')
#leg1 = legend(loc=1)
#leg1.draw_frame(False)
setp(sub2.get_xticklabels(), visible=False)
annotate('b',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')
axis([7.,16.,00.,299.99])


sub3=subplot(313)
sub3.yaxis.tick_left()
sub3.yaxis.set_label_position('left')
sub3.grid()
plot(tmfor, le, 'k-', label= '')
#plot (hour, conv_le*run1.lh, 'ro')
#plot (hour, conv_le*run2.lh, 'bo')
#plot (hour, conv_le*run3.lh, 'ko')
#plot (hour, conv_le*run4.lh, 'go')
#plot (hour, conv_le*avlh,             'ko', label= '')
plot (newhour, conv_le*newavlh,             'ko', label= '')
fill_between(newhour, conv_le*(newavlh-newstlh), conv_le*(newavlh+newstlh), color= 'k', label= '', alpha=0.2)
#fill_between(hour, conv_le*(avlh-stlh), conv_le*(avlh+stlh), color= 'k', label= '', alpha=0.2)
xlabel('LT [hours]')
ylabel('LE [Wm$^{-2}$]')
#leg1 = legend(loc=1)
#leg1.draw_frame(False)
annotate('c',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')
axis([7.,16.,00.,299.99])

#subplot(224)
#plot(tmfor, conv_fc*wcr, 'r-', label='Res')
#plot (run5.time_data + lttime , run5.res, 'ro')
#plot (run6.time_data + lttime , run6.res, 'g^')
#plot (run7.time_data + lttime , run7.res, 'k^')
#plot (run8.time_data + lttime , run8.res, 'b^')
#plot (hour, avres,             'ro', label= '')
#fill_between(hour, avres-stres, avres+stres, color= 'r', label= '', alpha=0.2)

#plot(tmfor, conv_fc*wca, 'g-', label= 'GEE')
#plot (run5.time_data + lttime , run5.gee, 'go')
#plot (run6.time_data + lttime , run6.gee, 'ro')
#plot (run7.time_data + lttime , run7.gee, 'ko')
#plot (run8.time_data + lttime , run8.gee, 'bo')
#plot (hour, avgee,             'go', label= '')
#fill_between(hour, avgee-stres, avgee+stgee, color= 'g', label= '', alpha=0.2)

#plot(tmfor, conv_fc*nee, 'k-', label= 'NEE')
#plot (run5.time_data + lttime , run5.nee, 'b^')
#plot (run9.time_data3, run9.fco2 + run9.fco2s, 'ko', label= '')
#plot(run5.time_data, zeroline, color='k',ls='dashed')
#xlabel('LT [hours]')
#ylabel('Carbon Fluxes [micromol/m2s]')
#leg1 = legend(loc=4)
#leg1.draw_frame(False)
#axis([7.,16.,-30.,8.])
#savefig('figure02.png',dpi=300)
savefig('figure02.png',dpi=300)
#savefig('fig_land.eps',dpi=300)
