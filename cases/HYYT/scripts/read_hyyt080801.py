#!/usr/bin/env python

######################    initializing    ########################
from matplotlib.pyplot import *
from numpy import *
from StringIO import StringIO
from matplotlib.ticker import OldScalarFormatter
from matplotlib.mathtext import *
import pickle 

######################    reading    ###############################
mdyn 	= loadtxt('/home/ruud/MXLCHgit/outhyyt/output_dyn', skiprows=3, comments='S')
msca 	= loadtxt('/home/ruud/MXLCHgit/outhyyt/output_sca', skiprows=3, comments='S')
mchem	= loadtxt('/home/ruud/MXLCHgit/outhyyt/chem_conc', skiprows=3)
mphoto  = loadtxt('/home/ruud/MXLCHgit/outhyyt/chem_photo', skiprows=3) 
mkeff   = loadtxt('/home/ruud/MXLCHgit/outhyyt/keff_cbl', skiprows=1) 
mftr    = loadtxt('/home/ruud/MXLCHgit/outhyyt/chem_ftr', skiprows=3) 
msoa    = loadtxt('/home/ruud/MXLCHgit/outhyyt/soa_part', skiprows=3) 

mhyyt 	   = loadtxt('/home/ruud/MXLCHgit/cases/HYYT/data/hyyt_meteo_08082001.txt', skiprows=1)
mhyyt_ams  = loadtxt('/home/ruud/MXLCHgit/cases/HYYT/data/Hyytiala2005_QAMS.txt', skiprows=1)
mhyyt_amsb = loadtxt('/home/ruud/MXLCHgit/cases/HYYT/data/Hyytiala2005_QAMS_binned.txt', skiprows=1)

# model 
at      = mdyn[:,0]+2 # convert UTC to LT
ah	= mdyn[:,2]
atheta  = mdyn[:,4]
awths	= mdyn[:,7]
awe     = mdyn[:,3]
adtheta = mdyn[:,5]
abeta   = mdyn[:,8]
aws     = mdyn[:,13]

at2	= msca[:,0]+2 # convert UTC to LT
ah2	= msca[:,2]
aq	= msca[:,3]
awqs	= msca[:,6]
adq     = msca[:,4]

at3	= mchem[:,0]+2 # convert UTC to LT
aO3	= mchem[:,3]
aO3_ft  = mftr[:,3]
aNO	= mchem[:,5]
aNO2	= mchem[:,6]
aOH	= mchem[:,13]
aISO	= mchem[:,11]
aINERT	= mchem[:,2]
aISORO2	= mchem[:,12]
aHO2	= mchem[:,14]
aMVK	= mchem[:,10]
aTERP   = mchem[:,24]
aCiT    = mchem[:,25]
aOAbg   = mchem[:,26]

aCoa     = msoa[:,3]
aNO2_ft  = mftr[:,6]
aTERP_ft = mftr[:,24]
aOAbg_ft = mftr[:,26]

aphi	= mphoto[:,5]
ajNO2	= mphoto[:,3]	#jno2*60(min)
ajO3	= mphoto[:,2]	#jo3*1000*60
azen	= mphoto[:,6]
asoz	= mphoto[:,8]

#aTref_cbl= mkeff[:,2]
#ak5	= mkeff[:,7]
#ak21	= mkeff[:,23]
#ak28    = mkeff[:,30]
#ak29    = mkeff[:,31]

# observations
ahour_h = mhyyt[:,3]
amin_h	= mhyyt[:,4]
asec_h	= mhyyt[:,5]
at_h	= ahour_h + amin_h/60 + asec_h/3600 
aT_h	= mhyyt[:,25]
atheta_h= mhyyt[:,38]
awths_h = mhyyt[:,9]
aq_h	= mhyyt[:,50]
awqs_h	= mhyyt[:,10]
aO3_h	= mhyyt[:,30]
aNO_h	= mhyyt[:,28]
aNOx_h	= mhyyt[:,29]
aNO2_h  = mhyyt[:,29]-mhyyt[:,28]
aH_h	= mhyyt[:,6]
aLE_h	= mhyyt[:,7]
aRn_h	= mhyyt[:,8]

ayear_ams   = mhyyt_ams[:,0]
amonth_ams  = mhyyt_ams[:,1]
aday_ams    = mhyyt_ams[:,2]
ahour_ams   = mhyyt_ams[:,3] + 2 # convert UTC to LT
amin_ams    = mhyyt_ams[:,4]
asec_ams    = mhyyt_ams[:,5]
at_ams      = ahour_ams + amin_ams/60 + asec_ams/3600 

aOOA1   = mhyyt_ams[:,6]
aOOA2   = mhyyt_ams[:,7]
aSO4    = mhyyt_ams[:,8]
aNH4    = mhyyt_ams[:,9]
aNO3    = mhyyt_ams[:,10]
aCF     = mhyyt_ams[:,11]

at_amsb      = mhyyt_amsb[:,0] 
aOOA2binmean = mhyyt_amsb[:,1]
aOOA2binmed  = mhyyt_amsb[:,4]

# parameters & constants
R         = 8.3145                      # gas constant (J mol-1 k-1)
T0        = 298.                        # temperature (K) at which the lab experiments and fitting are done, Tsimpidi'10
Na        = 6.02214e23                  # Avogadro's number (mol-1)
dHvap     = 30.e3   *1.0                # enthalpy of vaporization (J/mol), Lane'08

# balloon measurements
at_b    = [9.+30./60., 10.+50./60., 13., 15.]
aISO_b  = [0.0308,0.0121,0.0572,0.0136] # ug/m3
aapin_b = [0.00446,0.0914,0.0611,0.0896] # ug/m3
abpin_b = [0.0,0.0231,0.0130,0.0199] # ug/m3
atric_b = [nan,0.00412,0.00452,0.00583] # ug/m3
acamp_b = [0.0131,0.0133,0.0102,0.0172] # ug/m3
asabi_b = [nan,0.0101,0.00425,0.0198] # ug/m3
amyrc_b = [nan,nan,0.00268,0.00603] # ug/m3
acare_b = [0.00584,0.0333,0.0213,0.0364] # ug/m3
alimo_b = [0.0100,0.0319,0.0356,0.0715] # ug/m3
aTERP_b = [0.0334,0.207,0.153,0.266] # ug/m3 sum of all measured (but not all present) terpenes

######################   calculations    ##########################
## calculate phi from model & data
#ak5_d   = ak5[0] *ones((48))
#ak21_d  = ak21[0]*ones((48)) 
#ak21_d  = 3.0e-12*exp(-1500/aT_h) # convert units!
#ak5_d   = 1.67e-2*exp(-0.575/cos(1.1)) # convert units!

#aphi_m  = (aNO*aO3*ak21)    / (aNO2*ak5)
#aphi_h  = (aNO_h*aO3_h*ak21_d)/((aNOx_h-aNO_h)*ak5_d)

## convert kinematic to normal heat flux
rho_Cp  = 1.231e3 
aH      = awths*rho_Cp # W/m2
rho_Lv  = 3013.5
aLE     = awqs*rho_Lv   # W/m2

## convert modeled OH-concentration from ppb to molec./cm3
R	= 8.3145     # gas constant (J mol-1 k-1)
Na	= 6.02214e23 # Avogadro's number (mol-1)
p       = 99500	     # pressure (Pa)
T	= 290	     # temperature (K)	
V	= 1e-6       # volume (m3)
cf_OH   = 1e-9*((p*V*Na)/(R*T)) # factor to convert from ppb to molecules/cm3
aOHmr   = aOH # OH mixing ratio
aOH     = multiply(aOH,cf_OH) # molecules/cm3

# TERP-flux
mmTERP = 10.*12.01 + 16.*1.008 # molecular mass C10H16
fcTERP = mmTERP*(p/(R*T)*3600*1e-3) # flux conversion factor ppb m s-1 -> ug m-2 h-1
F_TERP = 15.e-3*fcTERP
print('F_TERP=',round(F_TERP),'ug m-2 h-1')
# ISO-flux
ISOmm   = 5.*12.01+8.*1.008 # molar mass C5H8 (g/mol)
fcISO = ISOmm*(p/(R*T)*3600*1e-3) # flux conversion factor ppb m s-1 -> ug m-2 h-1
F_ISO = 0.01*fcISO
print('F_ISO=',round(F_ISO),'ug m-2 h-1')
## convert balloon ISO-concentration from ug/m3 to ppb
cf_ISO  = (R*T)/(p*ISOmm) * 1e3  # factor to convert isoprene from ug/m3 to ppb
#aISO_b  = multiply(aISO_b,cf_ISO) # ppb
## convert model ISO-concentration from ppb to ug/m3 
mr2mc_ISO = (p*ISOmm)/(R*T) * 1e-3  # factor to convert terpenes from ppb to ug/m3 
aISOmc  = multiply(aISO,mr2mc_ISO) # ug/m3
## convert balloon terpene-concentration from ug/m3 to ppb
TERPmm   = 10.*12.01+16.*1.008 # molar mass C10H16 (g/mol)
mc2mr_TERP = (R*T)/(p*TERPmm) * 1e3  # factor to convert terpenes from ug/m3 to ppb
aTERP_bmr  = multiply(aTERP_b,mc2mr_TERP) # ppb
aapin_bmr  = multiply(aapin_b,mc2mr_TERP) # ppb
## convert model terpene-concentration from ppb to ug/m3 
mr2mc_TERP = (p*TERPmm)/(R*T) * 1e-3  # factor to convert terpenes from ppb to ug/m3 
aTERPmc  = multiply(aTERP,mr2mc_TERP) # ug/m3
#aapinmc  = multiply(aapin,mr2mc_TERP) # ug/m3
# calculate dthetav
adthetav = adtheta + 0.61*((aq/1e3)*adtheta + atheta*(adq/1e3) + adtheta*(adq/1e3))
adthetav_t1 = 0.61*((aq/1e3)*adtheta )
adthetav_t2 = 0.61*(atheta*(adq/1e3))
adthetav_t3 = 0.61*(adtheta*(adq/1e3))

# calculate wthvs
awthvs = awths+0.61*atheta*awqs*1.0e-3

# O3 deposition velocity
daylength = 43200/3600
O3_em = -0.25
awO3s = (O3_em * sin(pi * ((at-6.8)/daylength)))
aVdO3 = -awO3s/aO3

# filter out negative NO-measurements
aNO_h
for i in range(0,size(aNO_h,0)):
   if aNO_h[i] < 0.:
      aNO_h[i] = NaN
   #endif
#endfor
#std = sqrt(mean(abs(x - x.mean())**2))

## calculate Damkoehler numbers Da = t_t/t_c
#g          = 9.81           # gravitational constant (m s-2)
#aw_st      = ((g/atheta)*awthvs*ah)**(1./3.)
#at_t       = ah/aw_st
#akTERPO3   = ak28 * aO3
#akTERPOH   = ak29 * aOHmr
#aDa_TERPO3 = at_t*akTERPO3
#aDa_TERPOH = at_t*akTERPOH

###############################################################################################
#####################    plotting     ##############################

figure(1, figsize=(11.5,9)) # dynamics: wths, wqs, h, theta, q
#suptitle('hyyt080801 dynamics')
subplot(221) #
plot(at, aH, 'k-', at2, aLE, 'b--', at_h, aH_h, 'k.', at_h, aLE_h, 'b+')
xlabel('time LT (h)', fontsize=16)
ylabel(u'heat flux (W m\u207B\u00B2)', fontsize=16)
#xlim([0,24])
axis([6, 20, -50, 350])
text(6.5,305,'(a)', fontsize=16)
grid(True)
legend(('H', 'LE'), 'best')
subplot(222) # h
plot(at, ah ,'k-',(12.+20./60.),1000,'kx')
xlabel('time LT (h)', fontsize=16)
ylabel('h (m)', fontsize=16)
#xlim([0, 24])
axis([6, 20,0,1800])
text(6.5,1600,'(b)', fontsize=16)
grid(True)
subplot(223) # theta
plot(at,atheta,'k-', at_h, atheta_h, 'k.', lw=1)
xlabel('time LT (h)', fontsize=16)
ylabel(u'<\u03B8> (K)', fontsize=16)
#xlim([0, 24])
axis([6, 20,287,293])
text(6.5,292.2,'(c)', fontsize=16)
grid(True)
subplot(224) #
plot(at, aq ,'k-', at_h, aq_h, 'k.', lw=1)
xlabel('time LT (h)', fontsize=16)
ylabel(u'<q> (g kg\u207B\u00B9)', fontsize=16)
xlim([6,20])
axis([6, 20,5.5,8.5])
text(18.5,8.1,'(d)', fontsize=16)
grid(True)

figure(2, figsize=(11.5,9)) # chemistry : O3, NOx, OH
#suptitle('hyyt 080801 gas-phase chemistry')
subplot(2,2,1) #
plot(at, aO3,'k-',at_h, aO3_h,'k.')
xlabel('time LT (h)', fontsize=16)
ylabel(u'<O\u2083> (ppb)', fontsize=16)
#xlim([0, 24])
xlim([6, 20])
text(6.5,38.9,'(a)', fontsize=16)
grid(True)
subplot(2,2,2) #
plot(at3, aNO,'k-', at_h, aNO_h,'k.')
#errorbar(at_h, aNO_h, yerr=std(aNO_h), xerr=None, fmt='k.', ecolor='Black')
xlabel('time LT (h)', fontsize=16)
ylabel('<NO> (ppb)', fontsize=16)
#xlim([0, 24])
xlim([6, 20])
text(18.5,0.072,'(b)', fontsize=16)
grid(True)
subplot(2,2,3) #
plot(at3, aOH,'k-')
#errorbar(at_b, aOH_b, yerr=std(aOH_b), xerr=None, fmt='kx', ecolor='Black')
gca().yaxis.set_major_formatter(OldScalarFormatter())
ylabel(u'<OH> (molec. cm\u207B\u00B3)', fontsize=16)
xlabel('time LT (h)', fontsize=16)
#xlim([0, 24])
xlim([6, 20])
text(18.5,4.1e6,'(c)', fontsize=16)
grid(True)
subplot(2,2,4) #
plot(at3, aNO2,'k-', at_h,aNO_h*4,'k.') #at_h, aNO2_h,'k.',
#errorbar(at_h, aNO_h*5, yerr=std(aNO_h*5), xerr=None, fmt='k+', ecolor='Black')
xlabel('time LT (h)', fontsize=16)
ylabel(u'<NO\u2082> (ppb)', fontsize=16)
#xlim([0, 24])
xlim([6, 20])
#ylim([0,0.4])
text(18.5,0.27,'(d)', fontsize=16)
grid(True)

fig = figure(3, figsize=(6,5)) # TERP
##xlabel('time LT (h)', fontsize=16)
#ax1 = fig.add_subplot(121)     # mc
#ax2 = ax1.twinx()              # mr
#def aa(bb):
#    return (bb*mc2mr_TERP*1000)
#def update_ax2(ax1):
#   y1, y2 = ax1.get_ylim()
#   ax2.set_ylim(aa(y1), aa(y2))
#   ax2.figure.canvas.draw()
## automatically update ylim of ax2 when ylim of ax1 changes.
#ax1.callbacks.connect("ylim_changed", update_ax2)
#ax1.plot(at,aTERPmc,'k-')
#ax1.errorbar(at_b[0], aapin_b[0], yerr=multiply(aapin_b[0],1.0), xerr=None, fmt='kx', ecolor='Black')
#ax1.errorbar(at_b[1:4], aapin_b[1:4], yerr=multiply(aapin_b[1:4],0.5), xerr=None, fmt='kx', ecolor='Black')
#ax1.set_ylabel(u'<\u03B1-pinene> (\u00B5g m\u207B\u00B3)', fontsize=12)
#ax1.grid(True)
#xlim([6,20])
#ax2.set_ylabel(u'<\u03B1-pinene> (ppt)', fontsize=12)
#text(6.5,23,'(a)', fontsize=12)
xlabel('time LT (h)', fontsize=12)
ax1 = fig.add_subplot(111)     # mc
ax2 = ax1.twinx()              # mr
def aa(bb):
    return (bb*mc2mr_TERP*1000)
def update_ax2(ax1):
   y1, y2 = ax1.get_ylim()
   ax2.set_ylim(aa(y1), aa(y2))
   ax2.figure.canvas.draw()
# automatically update ylim of ax2 when ylim of ax1 changes.
ax1.callbacks.connect("ylim_changed", update_ax2)
ax1.plot(at,aTERPmc,'k-')
ax1.errorbar(at_b[0], aTERP_b[0], yerr=multiply(aTERP_b[0],0.5), xerr=None, fmt='kx', ecolor='Black')
ax1.errorbar(at_b[1:4], aTERP_b[1:4], yerr=multiply(aTERP_b[1:4],0.3), xerr=None, fmt='kx', ecolor='Black')
ax1.set_ylabel(u'<TERP> (\u00B5g m\u207B\u00B3)', fontsize=12)
ax1.grid(True)
xlim([6,20])
ax2.set_ylabel(u'<TERP> (ppt)', fontsize=12)
text(6.5,57,'(b)', fontsize=12)

figure(4) # OA
plot(at, aCoa,'r-', label='0.2')
grid(True)
leg1 = legend(loc=1)
leg1.draw_frame(False)
xlabel('time LT [h]')
ylabel(u'<Coa> (\u00B5g m\u207B\u00B3)')
axis([4.,20.,0.0,1.0])
a = axes([0.21,0.17,0.20,0.20], axisbg='w')
plot(at_amsb, aOOA2binmean,'r.-')#, at_amsb, aOOA2binmed,'b.-')
axis([0.0, 24.0, 0.2, 1.0])
leg1.draw_frame(False)
xlabel('time LT [h]')
ylabel(u'SV-OOA (\u00B5g m\u207B\u00B3)')
grid(True)
xticks(arange(0,30,6))
yticks(arange(0.2,1.1,0.2))
setp=(a)

#OAout = (aCoa,aX1,aX2,aX3,aX4,aB)  
#pickleFileName = 'OAout.pkl'
#pickleFile = open(pickleFileName, 'wb')
#pickle.dump(OAout, pickleFile, pickle.HIGHEST_PROTOCOL)
#pickleFile.close()

show()
