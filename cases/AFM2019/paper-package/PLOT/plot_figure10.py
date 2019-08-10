 # Simple plotting model FORTRAN output

from pylab import *
import copy
import numpy
from scipy.interpolate import spline
from scipy.interpolate import interpolate
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)

# opening files MXLCH (original file har181217 before change -22 to -28)
outdyn        = loadtxt('../AFM19-PAPER/output_dyn', skiprows=3)
outsca        = loadtxt('../AFM19-PAPER/output_sca', skiprows=3)
cabl          = loadtxt('../AFM19-PAPER/output_sca',  skiprows=3)
cftr13        = loadtxt('../AFM19-PAPER/C13',   skiprows=3)
cftr18        = loadtxt('../AFM19-PAPER/COO18',   skiprows=3)
isot          = loadtxt('../AFM19-PAPER/iso_delta',  skiprows=3)
co2var        = loadtxt('../AFM19-PAPER/output_ags', skiprows=3)
isob          = loadtxt('../AFM19-PAPER/iso_budget', skiprows=3)
fliso         = loadtxt('../AFM19-PAPER/flux_iso', skiprows=3)

# Defining constants

g         = 9.81    #[m/s2]
epsilon   = 0.622   #[g water / g air]
tabs      = 273.    #[K]
cp        = 1004.67 #[J/kg K]
Rd        = 287.053 #[J/kg K]
Rv        = 461.50  #[J/kg K]
Lv        = 2.5e6   #[J/kg]
e0        = 0.611   #[kPa]
gd        = -(g/cp) #[K/m]
gamma_dry = g/cp
rhoair    = 1.2
rholv     = 3013.5
MW_Air    = 28.97
MW_H2O    = 18
MW_CO2    = 44
height28  = 27.9
local_time= 5.
Rsc13     = 0.011057
Rscoo18   = 2.*0.0020052
timestep  = 60.
conv_hour = 3600. * 1000. # per hour and per mil

# defining arrays
tmfor       = outdyn[:,0] - local_time
h           = outdyn[:,2]
we          = outdyn[:,3]
theta       = outdyn[:,4]
q           = outsca[:,3]
wca         = co2var[:,2]
wcr         = co2var[:,3]
co2         = cabl[:,8]/1000.                               #ppm
co2ftr      = (cabl[:,9] + cabl[:,8])/1000.                 #ppm
c13abl      = cftr13[:,2]/1000.                             #ppm
c13ftr      = cftr13[:,3]/1000.                             #ppm
c18abl      = cftr18[:,2]/1000.                             #ppm
c18ftr      = cftr18[:,3]/1000.                             #ppm
deltac13    = isot[:,5]
deltac18    = isot[:,6]
wc13        = isob[:,2]
wc13_p      = isob[:,3]
wc13_r      = isob[:,4]
wc18        = isob[:,7]
wc18_p      = isob[:,8]
wc18_r      = isob[:,9]
wc18_par2   = isob[:,17]
wc18_p_par2 = isob[:,18]
neec13      = fliso[:,3]/1000.                               # in ppm m/s
neec18      = fliso[:,6]/1000.                               # in ppm m/s
neec18_p    = fliso[:,7]/1000.                               # in ppm m/s
neec18_r    = fliso[:,8]/1000.                               # in ppm m/s
nee         = (wca + wcr)*(MW_Air/MW_CO2)*(1.0/rhoair)       # in ppm m/s
wca_p       = (wca)*(MW_Air/MW_CO2)*(1.0/rhoair)             # in ppm m/s
wca_r       = (wcr)*(MW_Air/MW_CO2)*(1.0/rhoair)             # in ppm m/s
zeroline  = numpy.zeros(720)
# Additional arrays
#zeroline    = numpy.zeros(43200)
zeroline    = numpy.zeros(720)
tmform      = tmfor[:-1]
tmform3     = tmfor[:-3]
tmform5     = tmfor[:-5]
deltac13m   = deltac13[:-1]
deltac18m   = deltac18[:-1]
co2m        = co2[:-1]
hm          = h[:-1]
wem         = we[:-1]

# calculating budget terms

##(c13abl/co2)=(Rsc13*((deltac13/1000.)+1))

## surface #################################################
#delta-c13
sur1 = (1/(Rsc13*co2))*(1/h)*((c13abl/co2)*nee)
sur2 = (1/(Rsc13*co2))*(1/h)* (Rsc13*co2*(wc13/1000.))
sur  = sur1 + sur2 
surm = sur[:-1]

## surface
#delta
sur1   = (1/(Rsc13*co2))*(1/h)*((c13abl/co2)*nee)
#1 term
sur1_p = (1/(Rsc13*co2))*(1/h)*((c13abl/co2)*wca_p)
sur1_r = (1/(Rsc13*co2))*(1/h)*((c13abl/co2)*wca_r)
#
#2 term
sur2   = (1/(Rsc13*co2))*(1/h)* (Rsc13*co2*(wc13/1000.))
sur2_p = (1/(Rsc13*co2))*(1/h)* (Rsc13*co2*(wc13_p/1000.))
sur2_r = (1/(Rsc13*co2))*(1/h)* (Rsc13*co2*(wc13_r/1000.))

sur       = sur1 + sur2
surtotal  = sur1 + sur2 - (1/(Rsc13*co2))*(1/h)*(Rsc13*((deltac13/1000.)+1))*nee
#3 term
sur3_p  =  - (1/(Rsc13*co2))*(1/h)*(Rsc13*((deltac13/1000.)+1))*wca_p
sur3_r  =  - (1/(Rsc13*co2))*(1/h)*(Rsc13*((deltac13/1000.)+1))*wca_r
surm = sur[:-1]

#total term plant
sur_plantc13 = sur1_p + sur2_p + sur3_p
#total term soil
sur_soilc13 = sur1_r + sur2_r + sur3_r


#delta-c18
sur1c18   = (1/(Rscoo18*co2))*(1/h)*((c18abl/co2)*nee)
#1 term
sur1c18_p = (1/(Rscoo18*co2))*(1/h)*((c18abl/co2)*wca_p)
sur1c18_r = (1/(Rscoo18*co2))*(1/h)*((c18abl/co2)*wca_r)
#2 term
sur2c18        = (1/(Rscoo18*co2))*(1/h)* (Rscoo18*co2*(wc18/1000.))
sur2c18_par2   = (1/(Rscoo18*co2))*(1/h)* (Rscoo18*co2*(wc18_par2/1000.))
sur2c18_p      = (1/(Rscoo18*co2))*(1/h)* (Rscoo18*co2*(wc18_p/1000.))
sur2c18_p_par2 = (1/(Rscoo18*co2))*(1/h)* (Rscoo18*co2*(wc18_p_par2/1000.))
sur2c18_r      = (1/(Rscoo18*co2))*(1/h)* (Rscoo18*co2*(wc18_r/1000.))
#3 term
surco2_c18   = - (1/(Rscoo18*co2))*(1/h)*(Rscoo18*((deltac18/1000.)+1))*nee
surco2_c18_p = - (1/(Rscoo18*co2))*(1/h)*(Rscoo18*((deltac18/1000.)+1))*wca_p
surco2_c18_r = - (1/(Rscoo18*co2))*(1/h)*(Rscoo18*((deltac18/1000.)+1))*wca_r
#print surco2_c18 - surco2_c18_p - surco2_c18_r
#
surc18  = sur1c18 + sur2c18 
surc18  = sur1c18 + sur2c18 
surc18_par2  = sur1c18 + sur2c18_par2
surmc18 = surc18[:-1]

#total term plant
sur_plant        = sur1c18_p + sur2c18_p      + surco2_c18_p
sur_plant_par2   = sur1c18_p + sur2c18_p_par2 + surco2_c18_p
#total term soil 
sur_soil  = sur1c18_r + sur2c18_r + surco2_c18_r

#co2
co2sur = (1/h)*nee
co2surm= co2sur[:-1] 

#c13
c13sur = (1/h)*neec13
c13surm= c13sur[:-1] 

#c18
c18sur = (1/h)*neec18
c18sur_p = (1/h)*neec18_p
c18sur_r = (1/h)*neec18_r
c18surm= c18sur[:-1] 

## entrainment #################################################
#delta-c13
ent    = (1/(Rsc13*co2))*(1/h)  *we*(c13ftr-c13abl)
entm   = ent[:-1] 
#delta-c18
entc18 = (1/(Rscoo18*co2))*(1/h)*we*(c18ftr-c18abl)
entmc18= entc18[:-1] 
#co2
co2ent = (1/h)*we*(co2ftr-co2)
co2entm= co2ent[:-1] 
#c13 
c13ent = (1/h)*we*(c13ftr-c13abl)
c13entm= c13ent[:-1] 
#c13 
c18ent = (1/h)*we*(c18ftr-c18abl)
c18entm= c18ent[:-1] 

## surface co2 for c13 #############################################
surco2     = - (1/(Rsc13*co2))*(1/h)*(Rsc13*((deltac13/1000.)+1))*nee
## entrainment co2
entco2     = - (1/(Rsc13*co2))*(1/h)*(Rsc13*((deltac13/1000.)+1))*we*(co2ftr-co2) 
## surface co2 for c18 #############################################
surco2_c18 = - (1/(Rscoo18*co2))*(1/h)*(Rscoo18*((deltac18/1000.)+1))*nee
## entrainment co2
entco2_c18 = - (1/(Rscoo18*co2))*(1/h)*(Rscoo18*((deltac18/1000.)+1))*we*(co2ftr-co2) 

#####NEW ENTRAINMENT EXPRESSION
Rc13ft      = c13ftr/co2ftr
deltac13ft  = ((Rc13ft/Rsc13) - 1.) * 1000.
newentc13   = (1/h)* (we * (co2ftr/co2) * (deltac13ft - deltac13)) * (1/1000.)

#####NEW ENTRAINMENT EXPRESSION
Rc18ft      = c18ftr/co2ftr
deltac18ft  = ((Rc18ft/Rscoo18) - 1.) * 1000.
newentc18   = (1/h)* (we * (co2ftr/co2) * (deltac18ft - deltac18)) * (1/1000.)

## storage ################################################################## 
#delta13
sto     = -(1/(Rsc13*co2m))  *Rsc13*((deltac13m/1000.) + 1 )  *np.diff(co2)/timestep
#delta13
sto_c18 = -(1/(Rscoo18*co2m))*Rscoo18*((deltac18m/1000.) + 1 )*np.diff(co2)/timestep

## temporal delta ######################################################### 
#c13
tem     = np.diff((deltac13/1000))/timestep
#c18
tem_c18 = np.diff((deltac18/1000))/timestep

#co2
co2sto   = np.diff(co2)/timestep

#c13
c13sto   = np.diff(c13abl)/timestep
#c18
c18sto   = np.diff(c18abl)/timestep

## residual
#c13
res     = sur + ent + surco2 + entco2 
newres  = sur + surco2 + newentc13 
#c18
def moving_average(entc18, n=5) :
    ret = np.cumsum(entc18, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / float(n)
maentc18 = moving_average(entc18,n=6)
def moving_average(entco2_c18, n=5) :
    ret = np.cumsum(entco2_c18, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / float(n)
maentco2_c18 = moving_average(entco2_c18,n=6)

resc18    = surc18 + entc18     + surco2_c18 + entco2_c18 
newresc18 = surc18 + surco2_c18 + newentc18 
#maresc18  = surc18 + maentc18 + surco2_c18 + maentco2_c18 
#co2
co2resm = co2surm   + co2entm
co2res  = co2sur    + co2ent
#c13
c13res  = c13surm   + c13entm
#c18
c18resm = c18surm   + c18entm
c18res  = c18sur    + c18ent

# mpplying moving average to smooth curves
def moving_average(ent, n=5) :
    ret = np.cumsum(ent, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / float(n)
maent    = moving_average(ent,n=6)
def moving_average(entco2, n=5) :
    ret = np.cumsum(entco2, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / float(n)
maentco2 = moving_average(entco2,n=6)

def moving_average(res, n=5) :
    ret = np.cumsum(res, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / float(n)

mares= moving_average(res,n=6)

def moving_average(resc18, n=5) :
    ret = np.cumsum(resc18, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / float(n)

maresc18 = moving_average(resc18,n=6)

# plotting 
fig1=figure(2, figsize=(12,10))
fig1.subplots_adjust(0.09,0.08,0.91,0.94,0.02,0.02)
#ax = fig1.add_subplot(311)
sub1=subplot(211)
#ax.get_xaxis().tick_bottom()
sub1.xaxis.tick_top()
sub1.set_xlabel('LT [hours]')
sub1.xaxis.set_label_position('top')
sub1.grid()
sub1.plot(tmfor ,conv_hour*(sur + surco2)    , 'g',      label= 'sur')
sub1.plot(tmfor ,conv_hour*(sur_plantc13)       , 'g--',    label= 'sur$_{plant}$')
sub1.plot(tmfor ,conv_hour*(sur_soilc13)        , 'r--',    label= 'sur$_{soil}$')
#plot(tmfor ,conv_hour*(ent+entco2)      , 'b--',      label= 'ent')
sub1.plot(tmfor ,conv_hour*(newentc13)      , 'b',      label= 'ent')
#plot(tmfor ,conv_hour*(ent)      , 'b^',      label= 'ent')
#plot(tmfor ,conv_hour*(entco2)      , 'bo',      label= 'ent')
#plot(tmform5,conv_hour*(maent+maentco2)  , 'b',      label= 'ent')
#plot(tmfor, conv_hour*res               , 'k--',      label= 'tem')
#plot(tmform5, conv_hour*mares               , 'k',      label= 'tem')
sub1.plot(tmfor, conv_hour*newres              , 'k',     label= 'tem')
sub1.leg1 = legend(loc=4)
sub1.leg1.draw_frame(False)
sub1.set_xlabel('LT (hour)')
sub1.set_ylabel('($\partial \delta_a^{13}$/$\partial t$) [permil hour$^{-1}$]')
sub1.axis([7.,16.,-0.6,0.599])
sub1.annotate('a',xy=(0.95,0.95),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

# inset figure
ax1 = plt.axes([1,1,1,1])
ip  = InsetPosition(sub1, [0.40, 0.10, 0.3,0.3])
ax1.set_axes_locator(ip)
ax1.plot(tmfor,deltac13ft-deltac13, 'b',   linewidth=1)
ax1.plot(tmfor,zeroline, color='k',linewidth=1)
ax1.set_xlabel('LT [hours]', fontsize=14)
ax1.set_ylabel('$\Delta$ $\delta^{13}$ [permil]', fontsize=14)
ax1.grid()
ax1.set_xlim(7,9.)
ax1.set_ylim(-0.299999,0.3)
ax1.set_xticks(np.arange(7,9.,0.5))
ax1.set_yticks(np.arange(-.2,.200000001,0.2))

sub2=subplot(212)
sub2.grid()
plot(tmfor ,conv_hour*(surc18 + surco2_c18)    , 'g',  label= 'sur')
plot(tmfor ,conv_hour*(sur_plant)              , 'g--',label= 'sur$_{plant}$')
plot(tmfor ,conv_hour*(sur_soil)               , 'r--',label= 'sur$_{soil}$')
#plot(tmform5 ,conv_hour*(maentc18 + maentco2_c18)    , 'b',  label= 'ent')
plot(tmfor ,conv_hour*(newentc18)      , 'b',      label= 'ent')
plot(tmfor, conv_hour*newresc18          , 'k',     label= 'tem')
#plot(tmform5, conv_hour*maresc18   , 'k',  label= 'tem')
leg1 = legend(loc=1, ncol=1)
leg1.draw_frame(False)
xlabel('LT (hour)')
ylabel('($\partial \delta_a^{18}$/$\partial t$) [permil hour$^{-1}$]')
axis([7.,16.,-0.50,0.999])
annotate('b',xy=(0.05,0.95),xycoords='axes fraction',horizontalalignment='left', verticalalignment='bottom')

#savefig('fig_budcom.png',dpi=300)
savefig('figure10.png',dpi=300)
