 # Simple plotting model FORTRAN output

from pylab import *
import copy
import numpy
import datetime

# Reading Harvard data 
outrad09   = loadtxt('HF_Sept_2009_rad.dat', skiprows=0)
outrad10   = loadtxt('HF_Sept_2010_rad.dat', skiprows=0)
outrad12   = loadtxt('HF_Sept_2013_rad.dat', skiprows=0)

hour       = outrad09[848:896,5]
minute     = outrad09[848:896,6]
swi        = outrad09[848:896,7]
swo        = outrad09[848:896,8]
lwi        = outrad09[848:896,9]
lwo        = outrad09[848:896,10]

hour10     = outrad10[374:422,5]
minute10   = outrad10[374:422,6]
swi10      = outrad10[374:422,7]
swo10      = outrad10[374:422,8]
lwi10      = outrad10[374:422,9]
lwo10      = outrad10[374:422,10]

hour12     = outrad12[683:731,5]
minute12   = outrad12[683:731,6]
swi12      = outrad12[683:731,7]
swo12      = outrad12[683:731,8]
lwi12      = outrad12[683:731,9]
lwo12      = outrad12[683:731,10]

time = hour + minute/60. 
qn   = swi - swo + lwi - lwo
time10 = hour10 + minute10/60. 
qn10   = swi10 - swo10 + lwi10 - lwo10
time12 = hour12 + minute12/60. 
qn12   = swi12 - swo12 + lwi12 - lwo12

# plotting chemical species
figure(2, figsize=(12,10))
subplot(221)
plot(time, swi, 'r')
plot(time, swo, 'k')
plot(time10, swi10, 'r^')
plot(time10, swo10, 'k^')
plot(time12, swi12, 'ro')
plot(time12, swo12, 'ko')
xlabel('Time (h UTC)')
ylabel('[SW] (Wm-2)')
annotate('a',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

subplot(222)
plot(time, swo/swi, 'b')
plot(time10, swo10/swi10, 'b^')
plot(time12, swo12/swi12, 'bo')
ylabel('[albedo] (-)')
xlabel('Time (h UTC)')
axis([10.,24.,0,0.4])
annotate('b',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

subplot(223)
plot(time, lwi, 'r')
plot(time, lwo, 'k')
plot(time10, lwi10, 'r^')
plot(time10, lwo10, 'k^')
plot(time12, lwi12, 'ro')
plot(time12, lwo12, 'ko')
xlabel('Time (h UTC)')
ylabel('[LW] (Wm-2)')
annotate('c',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

subplot(224)
plot(time, qn,  'r-')
plot(time10, qn10,'r^')
plot(time12, qn12,'ro')
xlabel('Time (h UTC)')
ylabel('[Qn] (micromol/m2s)')
savefig('fig_rad.png',dpi=300)
#savefig('fig_land.eps',dpi=300)
