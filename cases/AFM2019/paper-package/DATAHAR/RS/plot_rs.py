 # Simple plotting model FORTRAN output

from pylab import *
import copy
import numpy
import datetime

# Reading Harvard data 
rs12   = loadtxt('110919_12.txt', skiprows=1)
rs00   = loadtxt('110920_00.txt', skiprows=1)

pressure   = rs00[0:29,0]
temp       = rs00[0:29,1]
dtemp      = rs00[0:29,2]
pressure12 = rs12[0:29,0]
temp12     = rs12[0:29,1]
dtemp12    = rs12[0:29,2]

# create array
vlev    = 29
z       = zeros(vlev)
rho     = zeros(vlev)
dz      = zeros(vlev)
z12     = zeros(vlev)
rho12   = zeros(vlev)
dz12    = zeros(vlev)

# height sounding
z0      = 16.

# constants
dp      = 0.5
g       = 9.81    #[m/s2]
epsilon = 0.622   #[g water / g air]
tabs    = 273.15    #[K]
cp      = 1004.67 #[J/kg K]
Rd      = 287.053 #[J/kg K]
Rv      = 461.50  #[J/kg K]
Lv      = 2.5e6   #[J/kg]
e0      = 0.611   #[kPa]
gd      = -(g/cp) #[K/m]
gamma_dry = g/cp

tempK  = temp + tabs
dtempK = dtemp + tabs
tempK12  = temp12 + tabs
dtempK12 = dtemp12 + tabs

theta   = tempK  *(pressure[0]  /pressure)  **(Rd/cp)
theta12 = tempK12*(pressure12[0]/pressure12)**(Rd/cp)

# Parcel calculation
for i in range(0,vlev):
  if i == 0:
     dz[i]    = 0.
     z[i]     = z0
     dz12[i]    = 0.
     z12[i]     = z0
  else:
     rho[i]   = pressure[i]/(Rd*tempK[i])
     dp       = pressure[i] - pressure[i-1] 
     dz[i]    = -dp/(g*rho[i])
     z[i]     = z[i-1]+dz[i]
     rho12[i] = pressure12[i]/(Rd*tempK12[i])
     dp12     = pressure12[i] - pressure12[i-1] 
     dz12[i]  = -dp12/(g*rho12[i])
     z12[i]   = z12[i-1]+dz12[i]

# calculation RH and q
es    = e0*exp((Lv/Rv)*((1/tabs)-(1/tempK)))
e     = e0*exp((Lv/Rv)*((1/tabs)-(1/dtempK)))
rh    = e/es
rv    = epsilon * e /(pressure/10. -e)
q     = rv/(1 + rv) 
#
es12  = e0*exp((Lv/Rv)*((1/tabs)-(1/tempK12)))
e12   = e0*exp((Lv/Rv)*((1/tabs)-(1/dtempK12)))
rh12  = e12/es12
rv12  = epsilon * e12 /(pressure12/10. -e12)
q12   = rv12/(1 + rv12) 

# plotting chemical species
figure(2, figsize=(12,10))
subplot(121)
plot(tempK, z, 'r')
plot(theta, z, 'k', label='18 LT') 
plot(theta12, z12, 'k--', label='06 LT') 
xlabel('T (K)')
ylabel('height (m)')
leg1 = legend(loc=1)
leg1.draw_frame(False)
axis([270.,310.,0,3000])
annotate('a',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

subplot(122)
plot(1000. * q, z, 'b', label='18 LT')
plot(1000. * q12, z12, 'b--', label='06 LT')
xlabel('q (g/kg)')
ylabel('height (m)')
leg1 = legend(loc=1)
leg1.draw_frame(False)
axis([0.,10.,0,3000])
annotate('b',xy=(0.01,0.9),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

savefig('fig_rs.png',dpi=300)
#savefig('fig_land.eps',dpi=300)
