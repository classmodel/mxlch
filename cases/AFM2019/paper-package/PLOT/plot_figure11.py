from pylab import *
import numpy
import datetime
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)
import matplotlib.pylab as pl
matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)

# opening files MXLCH (original before change -22 to -28 file har240218)
sca01    = loadtxt('../AFM19-PAPER/iso_delta',  skiprows=3)
sca02    = loadtxt('../AFM19-EDD/iso_delta',  skiprows=3)
sca03    = loadtxt('../AFM19-ELC/iso_delta',  skiprows=3)
sca04    = loadtxt('../AFM19-ECM/iso_delta',  skiprows=3)
chm01    = loadtxt('../AFM19-PAPER/chem_conc',  skiprows=3)
chm02    = loadtxt('../AFM19-EDD/chem_conc',  skiprows=3)
chm03    = loadtxt('../AFM19-ELC/chem_conc',  skiprows=3)
chm04    = loadtxt('../AFM19-ECM/chem_conc',  skiprows=3)
lan01    = loadtxt('../AFM19-PAPER/output_land',  skiprows=3)
lan02    = loadtxt('../AFM19-EDD/output_land',  skiprows=3)
lan03    = loadtxt('../AFM19-ELC/output_land',  skiprows=3)

conv      = 3013.5 # conversion factor

deltac1801= sca02[:,6]
deltac1802= sca02[:,6]
deltac1803= sca03[:,6]
deltac1804= sca04[:,6]
deltaw1801= sca01[:,7]
deltaw1802= sca02[:,7]
deltaw1803= sca03[:,7]
deltaw1804= sca04[:,7]
wq_s01    = lan01[:,8]
wq_s02    = lan02[:,8]
wq_s03    = lan03[:,8]
q_01      = sca01[:,3]
q_02      = sca02[:,3]
q_03      = sca03[:,3]
tm        = sca01[:,0] - 5.0
delta_e01 = sca01[:,9]
delta_e02 = sca02[:,9]
delta_e03 = sca03[:,9]


wc1801= sca01[:,3]
wc1802= sca02[:,3]
wc1803= sca03[:,3]
wc1804 = sca04[:,3]
wh1801= sca01[:,4]
wh1802= sca02[:,4]
wh1803= sca03[:,4]
wh1804= sca04[:,4]

co201     = chm01[:,29]
co202     = chm02[:,29]
co203     = chm03[:,29]
co204     = chm04[:,29]
h2o01     = chm01[:,16]
h2o02     = chm02[:,16]
h2o03     = chm03[:,16]
h2o04     = chm04[:,16]

mfluxc1801 = wc1801*(co201/1000.)*40.9
mfluxc1802 = wc1802*(co202/1000.)*40.9
mfluxc1803 = wc1803*(co203/1000.)*40.9
mfluxc1804 = wc1804*(co204/1000.)*40.9

mfluxw1801 = wh1801*(h2o01/1000000.)*40.9
mfluxw1802 = wh1802*(h2o02/1000000.)*40.9
mfluxw1803 = wh1803*(h2o03/1000000.)*40.9
mfluxw1804 = wh1804*(h2o04/1000000.)*40.9

fig1=figure(3, figsize=(17,8))
fig1.subplots_adjust(0.09,0.08,0.91,0.94,0.02,0.02)
#fig1.subplots_adjust(bottom=0.15,wspace=0.05)
sub1=subplot(1,2,1)
sub1.plot(deltaw1801[0:423], deltac1801[0:423], 'ko',   linewidth=1, label= 'CTL')
sub1.plot(deltaw1802[0:423], deltac1802[0:423], 'ro',   linewidth=1, label= 'EDA')
sub1.plot(deltaw1803[0:423], deltac1803[0:423], 'go',   linewidth=1, label= 'ELC')
sub1.plot(deltaw1804[0:423], deltac1804[0:423], 'bo',   linewidth=1, label= 'ECM')
sub1.leg1 = legend(loc=3)
sub1.leg1.draw_frame(False)
sub1.grid()
sub1.set_xlabel('$\delta^{18}_{wa}$  [permil]', fontsize=16)
sub1.set_ylabel('$\delta^{18}_a$  [permil]', fontsize=16)
sub1.set_xlim(-23.5,-21.51111)
sub1.set_ylim(38.2,38.8)
#axis([-23.5, -20.49999, 38.2,39.0])
sub1.annotate('a',xy=(0.05,0.95),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom', fontsize=20)
sub1.annotate('Event', xy=(-21.7, 38.28), xytext=(-22.1, 38.28), fontsize=16,
          arrowprops=dict(facecolor='black', shrink=0.05),
			                )
# inset delta_lw
#ax = fig1.add_subplot(122)
ax1 = plt.axes([1,1,1,1]) 
ip  = InsetPosition(sub1, [0.65, 0.65, 0.3,0.3])
ax1.set_axes_locator(ip)
ax1.plot(tm,deltaw1801, 'k',   linewidth=1, label= 'CTL')
ax1.plot(tm,deltaw1802, 'r',   linewidth=1, label= 'EDA')
#ax1.plot(tm,deltaw1803[0:423], 'g',   linewidth=1, label= 'ELC')
ax1.set_xlabel('LT [hours]', fontsize=16)
ax1.set_ylabel('$\delta_{wa}^{18}$ [permil]', fontsize=16)
ax1.grid()
ax1.set_xlim(7,10.)
ax1.set_ylim(-22.9999999,-21.0)
ax1.set_xticks(np.arange(7,10.,1.))
ax1.set_yticks(np.arange(-23.,-21.,0.5))
#axis([8.0,8.4,44.,47.])
# end inset

sub2=subplot(1,2,2)
#ax = fig1.add_subplot(122)
sub2.yaxis.tick_right()
sub2.yaxis.set_label_position('right')
sub2.plot(mfluxw1801[0:423], mfluxc1801[0:423], 'ko',   linewidth=1, label= 'CTL')
sub2.plot(mfluxw1802[0:423], mfluxc1802[0:423], 'ro',   linewidth=1, label= 'EDA')
sub2.plot(mfluxw1803[0:423], mfluxc1803[0:423], 'go',   linewidth=1, label= 'ELC')
sub2.plot(mfluxw1804[0:423], mfluxc1804[0:423], 'bo',   linewidth=1, label= 'ECM')
sub2.leg1 = legend(loc=4)
sub2.leg1.draw_frame(False)
sub2.grid()
sub2.set_xlabel('F H$_2$$^{18}$O [permil mmol m$^{-2}$s$^{-1}$]]', fontsize=16)
sub2.set_ylabel('F C$^{18}$OO [permil $\mu$mol m$^{-2}$s$^{-1}$]', fontsize=16)
sub2.axis([10.000001, 90., -130.,-60.0])
sub2.annotate('b',xy=(0.95,0.95),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom', fontsize=20)
#sub2.annotate('Event', xy=(22,-118), xytext=(20.1, -126),
sub2.annotate('Event', xy=(22,-118), xytext=(18.6, -126), fontsize=16,
		          arrowprops=dict(facecolor='black', shrink=0.01),
			                                         )
# inset delta_lw
#ax = fig1.add_subplot(122)
ax1 = plt.axes([0,0,1,1]) 
ip  = InsetPosition(sub2, [0.15, 0.65, 0.3,0.3])
ax1.set_axes_locator(ip)
ax1.plot(tm,delta_e01, 'k',   linewidth=1, label= 'CTL')
ax1.plot(tm,delta_e02, 'r',   linewidth=1, label= 'EDA')
#ax1.plot(tm,delta_e03, 'g',   linewidth=1, label= 'ELC')
ax1.set_xlabel('LT [hours]', fontsize=16)
ax1.set_ylabel('$\delta_e^{18}$ [permil]', fontsize=16)
ax1.grid()
ax1.set_xlim(8,8.4)
ax1.set_ylim(44.000001,47)
ax1.set_xticks(np.arange(8,8.4,0.2))
ax1.set_yticks(np.arange(45.,47,1))
#axis([8.0,8.4,44.,47.])
# end inset
#fig1.tight_layout()
#savefig('fig11_scatter.pdf',dpi=300)
savefig('figure11.png',dpi=300)
