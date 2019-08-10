 # Simple plotting model FORTRAN output

from pylab import *
import copy
import numpy
import datetime

# opening files MXLCH (original har261117)
co2rat    = loadtxt('../AFM19-PAPER/conc_co2', skiprows=3)

# constants and atmosphere variables
R = 8.314   # J/mol/K
T = 288     # K
P =100000 # Pa

# defining arrays
tmfor          = co2rat[:,0] - 5.0
co2s           = co2rat[:,2]
co2cs          = co2rat[:,3]
co2i           = co2rat[:,4]
co2a           = co2rat[:,5]
ratco2s        = co2rat[:,6]
ratco2cs       = co2rat[:,7]
ratco2i        = co2rat[:,8]
gm             = co2rat[:,9]

# plotting chemical species
fig1=figure(2, figsize=(9,10))
fig1.subplots_adjust(0.09,0.08,0.91,0.94,0.02,0.02)
sub1=subplot(211)
sub1.xaxis.tick_top()
sub1.set_xlabel('LT [hours]')
sub1.xaxis.set_label_position('top')
#sub1.grid()
sub1.plot(tmfor, co2a , 'k',   label= 'CO$_2$$_{atm}$')
sub1.plot(tmfor, co2i,  'k--', label= 'CO$_2$$_{leaf}$')
sub1.plot(tmfor, co2cs, 'k-.', label= 'CO$_2$$_{cs}$')
sub1.plot(tmfor, co2s,  'k:',  label= 'CO$_2$$_{soil}$')
xlabel('LT [hours]')
ylabel('CO$_2$ [ppm]')
leg1 = legend(loc=2, ncol=2)
leg1.draw_frame(False)
axis([7.,16.,100.,1000.])
annotate('n',xy=(0.1,2.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

sub2=subplot(212)
#sub2.grid()
sub2.plot(tmfor, ratco2i , 'k--', label= 'CO$_2$$_{leaf}$ /CO$_2$$_{atm}$')
sub2.plot(tmfor, ratco2cs, 'k-.', label= 'CO$_2$$_{cs}$   /CO$_2$$_{atm}$')
sub2.plot(tmfor, ratco2s,  'k:',  label= 'CO$_2$$_{soil}$ /CO$_2$$_{atm}$')
xlabel('LT [hours]')
ylabel(' [-]')
leg1 = legend(loc=2, ncol=2)
leg1.draw_frame(False)
axis([7.,16.,0.,2.999])
annotate('a',xy=(16.8,2.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

savefig('figure06.png',dpi=300)

fig1=figure(3, figsize=(9,10))
sub1=subplot(211)
xlabel('LT [hours]')
ylabel(' [m/s]')
plot(tmfor, gm , 'k', label= 'gm')
sub1=subplot(212)
plot(tmfor, (1000*gm*P)/(R*T) , 'k', label= 'gm')
xlabel('LT [hours]')
ylabel(' [mmol/m2s]')
#savefig('figure06.png',dpi=300)
#savefig('fig_land.eps',dpi=300)
