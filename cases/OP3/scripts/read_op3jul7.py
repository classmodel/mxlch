#!/usr/bin/env python

######################    initializing    ########################
from matplotlib.pyplot import *
from numpy import *
from StringIO import StringIO
from matplotlib.ticker import OldScalarFormatter
from matplotlib.mathtext import *
#import pickle 
from scipy.interpolate import interp1d
from scipy import linspace, polyval, polyfit, sqrt, stats
#from meanbin import meanbin

######################    reading model data   ###############################
mdyn 	= loadtxt('/home/ruud/MXLCH_SOA/outop3/output_dyn', skiprows=3, comments='S')
msca 	= loadtxt('/home/ruud/MXLCH_SOA/outop3/output_sca', skiprows=3, comments='S')
mchem	= loadtxt('/home/ruud/MXLCH_SOA/outop3/chem_conc', skiprows=3)
mphoto  = loadtxt('/home/ruud/MXLCH_SOA/outop3/chem_photo', skiprows=3) 
mkeff   = loadtxt('/home/ruud/MXLCH_SOA/outop3/keff_cbl', skiprows=1) 
mftr    = loadtxt('/home/ruud/MXLCH_SOA/outop3/chem_ftr', skiprows=3) 
msoa    = loadtxt('/home/ruud/MXLCH_SOA/outop3/soa_part', skiprows=3) 
mland   = loadtxt('/home/ruud/MXLCH_SOA/outop3/output_land', skiprows=3) 
mbvoc   = loadtxt('/home/ruud/MXLCH_SOA/outop3/voc_em', skiprows=3) 
#msoaft  = loadtxt('/home/ruud/MXLCH_SOA/outop3/soa_ftr', skiprows=3) 
mPLOH   = loadtxt('/home/ruud/MXLCH_SOA/outop3/PL/OH', skiprows=3) 
mPLTERP = loadtxt('/home/ruud/MXLCH_SOA/outop3/PL/TERP', skiprows=3)
mPLISO  = loadtxt('/home/ruud/MXLCH_SOA/outop3/PL/ISO', skiprows=3)

######################    reading measurement data   ###############################
mNOx	= loadtxt('/home/ruud/MXLCH_SOA/cases/OP3/data/ceh-nox-gradient_bukit-atur_20080625.na', skiprows=40) # 30, 45, 60, 75m
mO3     = loadtxt('/home/ruud/MXLCH_SOA/cases/OP3/data/ceh-o3-gradient_bukit-atur_20080625.na',skiprows=34) # 30, 45, 60, 75m
mmeteo  = loadtxt('/home/ruud/MXLCH_SOA/cases/OP3/data/man-sonic-2_bukit-atur_20080621_micromet.na', skiprows=42) # 45m
mPTRMS  = loadtxt('/home/ruud/MXLCH_SOA/cases/OP3/data/lanc-ptr-ms_bukit-atur_20080622.na', skiprows=35) # 75m
mAMS    = loadtxt('/home/ruud/MXLCH_SOA/cases/OP3/data/man-ams_bukit-atur_20080624_r3.na', skiprows=24) # 33m
mP      = loadtxt('/home/ruud/MXLCH_SOA/cases/OP3/data/cam-diracfmw_bukit-atur_20080616.na', skiprows=33)
mFAGE   = loadtxt('/home/ruud/MXLCH_SOA/cases/OP3/data/leeds-fage_bukit-atur_20080706_version3.na', skiprows=37) # 5m
mPMF    = loadtxt('/home/ruud/MXLCH_SOA/cases/OP3/data/OP3PMFfacs.txt', skiprows=1)
mO3_5   = loadtxt('/home/ruud/MXLCH_SOA/cases/OP3/data/york-o3_bukit-atur_20080624.na', skiprows=24) # 5m
mNOx_5  = loadtxt('/home/ruud/MXLCH_SOA/cases/OP3/data/york-noxy_bukit-atur_20080624.na', skiprows=41) # 5m
mBVOC   = loadtxt('/home/ruud/MXLCH_SOA/cases/OP3/data/BVOC_MIXING RATIOS_OP3_3_LANCASTER.txt', skiprows=1) # 5m

### LEAP YEAR!!! convert to Day Of Year and Local Time ###
at_AMS  = mAMS[:,0] + (8./24.) # UTC to UTC+8 (LT)
aOA     = mAMS[:,1]
aSO4    = mAMS[:,2]
aNO3    = mAMS[:,3]
aNH4    = mAMS[:,4]
aChl    = mAMS[:,5]

at_PMF  = mPMF[:,4]
a82fac  = mPMF[:,5]
aOOA1   = mPMF[:,6]
aOOA2   = mPMF[:,7]
a91fac  = mPMF[:,8]

at_meteo= mmeteo[:,0] + (8./24.)   # from UTC to LT
aH      = mmeteo[:,12]
aLE     = mmeteo[:,13]
aT      = mmeteo[:,7] + 273.15 # from C to K
aRH     = mmeteo[:,8]
auv     = mmeteo[:,1]
awdir   = mmeteo[:,2]

at_P    = mP[:,0]+(8./24.)   # from UTC+8 to LT
aP      = mP[:,1]*10.        # from kPa to hPa

at_O3   = mO3[:,0] - (39624-177) # convert from days since 1900-01-01 00:00:00 UTC+8 to days since 2008-01-01 00:00:00 LT 
aO3_30  = mO3[:,1]   # CEH 30m
aO3_45  = mO3[:,2]   # CEH 45m
aO3_60  = mO3[:,3]   # CEH 60m
aO3_75  = mO3[:,4]   # CEH 75m
aO3 = aO3_75

at_NOx  = mNOx[:,0] - (39624-177) # convert from days since 1900-01-01 00:00:00 UTC+8 to days since 2008-01-01 00:00:00 LT 
aNO_30  = mNOx[:,1]  # CEH 30m
aNO_45  = mNOx[:,3]  # CEH 45m
aNO_60  = mNOx[:,5]  # CEH 60m
aNO_75  = mNOx[:,7]  # CEH 75m
aNO = aNO_75

aNO2_30 = mNOx[:,2]  # CEH 30m
aNO2_45 = mNOx[:,4]  # CEH 45m
aNO2_60 = mNOx[:,6]  # CEH 60m
aNO2_75 = mNOx[:,8]  # CEH 75m
aNO2 = aNO2_75

at_O3_5 = mO3_5[:,0] + (8./24.)   # from UTC to LT
aO3_5   = mO3_5[:,1] # York 5m 

at_NOx_5 = mNOx_5[:,0] + (8./24.)   # from UTC to LT
aNO_5    = mNOx_5[:,1] # York 5m
aNO2_5   = mNOx_5[:,3] # York 5m

at_PTRMS= mPTRMS[:,0] - (39621-174) # convert from days since 1900-01-01 00:00:00 UTC+8 to days since 2008-01-01 00:00:00 LT 
aISO    = mPTRMS[:,1]
aFISO   = mPTRMS[:,2]
aTERP   = mPTRMS[:,5]
aFTERP  = mPTRMS[:,6]

at_BVOC = mBVOC[:,5] # DoY LT (UTC+8)
aMVK    = mBVOC[:,8] # MVK + MACR

at_FAGE = mFAGE[:,0] # LT
at_OH   = copy(at_FAGE)
at_HO2  = copy(at_FAGE)
aOH     = mFAGE[:,1]
aHO2    = mFAGE[:,2]

atm     = mdyn[:,0]+8 # convert UTC to LT
ahm	= mdyn[:,2]
athetam  = mdyn[:,4]
awthsm	= mdyn[:,7]
awem     = mdyn[:,3]
adthetam = mdyn[:,5]
abetam   = mdyn[:,8]
awsm     = mdyn[:,13]

at2m	= msca[:,0]+8 # convert UTC to LT
ah2m	= msca[:,2]
aqm	= msca[:,3]
awqsm	= msca[:,6]
adqm     = msca[:,4]

at3m	= mchem[:,0]+8 # convert UTC to LT
aO3m	= mchem[:,3]
aO3_ftm  = mftr[:,3]
aNOm	= mchem[:,5]
aNO2m	= mchem[:,6]
aOHm	= mchem[:,13]
aISOm	= mchem[:,11]
aINERTm	= mchem[:,2]
aISORO2m= mchem[:,12]
aHO2m	= mchem[:,14]
aMVKm	= mchem[:,10]
aTERPm    = mchem[:,24]

aOH_LR09 = -mPLOH[:,5] 
aOH_PR19 =  mPLOH[:,15]

aTERP_ftm = mftr[:,24]
atm_pl    = mPLTERP[:,0]+8 # convert UTC to LT
aTERPm_l    = -mPLTERP[:,5]
aTERPm_LR29 = -mPLTERP[:,3]
aTERPm_LR30 = -mPLTERP[:,4]

aISO_ftm = mftr[:,11]
aISOm_l    = -mPLISO[:,4]
aISOm_LR09 = -mPLISO[:,3]

aOAbgm    = msoa[:,2]
aCoam     = msoa[:,3]
aC1Tm     = msoa[:,4]
aC2Tm     = msoa[:,5]
aC3Tm     = msoa[:,6]
aC4Tm     = msoa[:,7]
aC1Im     = msoa[:,8]
aC2Im     = msoa[:,9]
aC3Im     = msoa[:,10] 
aX1m      = msoa[:,11]
aX2m      = msoa[:,12]
aX3m      = msoa[:,13]
aX4m      = msoa[:,14] 
amr2mc_PRODT = msoa[:,15]
amr2mc_PRODI = msoa[:,16]
abeta_terp = msoa[:,17]
abeta_iso  = msoa[:,18]

aSwin     = mland[:,2]

aFiso_meg     = mbvoc[:,2]
aFterp_meg    = mbvoc[:,3]
#aOAbg_ftm = msoaft[:,2] # ug m-3!
#aCoam_ftm = msoaft[:,3]
#aC1T_ftm = msoaft[:,4]
#aC2T_ftm = msoaft[:,5]
#aC3T_ftm = msoaft[:,6]
#aC4T_ftm = msoaft[:,7]
#aC1I_ftm = msoaft[:,8]
#aC2I_ftm = msoaft[:,9]
#aC3I_ftm = msoaft[:,10]

aphim	= mphoto[:,5]
ajNO2m	= mphoto[:,3]	#jno2*60(min)
ajO3m	= mphoto[:,2]	#jo3*1000*60
azenm	= mphoto[:,6]
asozm	= mphoto[:,7]

# parameters & constants
R         = 8.3145                      # gas constant (J mol-1 k-1)
T0        = 298.                        # temperature (K) at which the lab experiments and fitting are done, Tsimpidi'10
Na        = 6.02214e23                  # Avogadro's number (mol-1)
dHvap     = 30.e3   *1.0                # enthalpy of vaporization (J/mol), Lane'08

######################   calculations    ##########################

## convert kinematic to normal heat flux
rho_Cp  = 1.231e3 
aHm      = awthsm*rho_Cp # W/m2
rho_Lv  = 3013.5
aLEm     = awqsm*rho_Lv   # W/m2

# recycling OH in iso oxidation
ar_OH = aOH_PR19/-aOH_LR09

## convert modeled OH- and HO2-concentration from ppb to molec./cm3
R	= 8.3145     # gas constant (J mol-1 k-1)
Na	= 6.02214e23 # Avogadro's number (mol-1)
p       = 95000	     # pressure (Pa)
T	= 302	     # temperature (K)	
V	= 1e-6       # volume (m3)
cf_OH   = 1e-9*((p*V*Na)/(R*T)) # factor to convert from ppb to molecules/cm3
aOHmrm   = copy(aOHm) # OH mixing ratio
aOHm     = multiply(aOHm,cf_OH) # molecules/cm3
aHO2mrm  = copy(aHO2m) # OH mixing ratio
aHO2m    = multiply(aHO2m,cf_OH) # molecules/cm3

# TERP-flux
daylength = 42000./3600. # not the same as run lenght!!
mmTERP = 10*12.01 + 16*1.008 # molecular mass C10H16
fcTERP = mmTERP*(p/(R*T)*3600*1e-6) # flux conversion factor ppb m s-1 -> mg m-2 h-1
TERPmc_em  = 0.04*fcTERP
aFTERPm  = (TERPmc_em * sin(pi * ((atm-6.5)/daylength)))
## convert model ISO-concentration from ppb to ug/m3 
mmISO   = 5*12.01+8*1.008 # molar mass C5H8 (g/mol)
mr2mc_ISO = (p*mmISO)/(R*T) * 1e-3  # factor to convert terpenes from ppb to ug/m3 
aISOmcm  = multiply(aISOm,mr2mc_ISO) # ug/m3
fcISO = mmISO*(p/(R*T)*3600*1e-6) # flux conversion factor ppb m s-1 -> mg m-2 h-1
ISOmc_em  = 0.35*fcISO
aFISOm   = (ISOmc_em * sin(pi * ((atm-6.5)/daylength)))
## convert model SVOC concentration from ppb to ug/m3 
aC1Tmcm  = multiply(aC1Tm,amr2mc_PRODT) # ug/m3
aC2Tmcm  = multiply(aC2Tm,amr2mc_PRODT) # ug/m3
aC3Tmcm  = multiply(aC3Tm,amr2mc_PRODT) # ug/m3
aC4Tmcm  = multiply(aC4Tm,amr2mc_PRODT) # ug/m3
aC1Imcm  = multiply(aC1Im,amr2mc_PRODI) # ug/m3
aC2Imcm  = multiply(aC2Im,amr2mc_PRODI) # ug/m3
aC3Imcm  = multiply(aC3Im,amr2mc_PRODI) # ug/m3

aC1mcm   = sum((aC1Tmcm, aC1Imcm),0)  
aC2mcm   = sum((aC2Tmcm, aC2Imcm),0)
aC3mcm   = sum((aC3Tmcm, aC3Imcm),0)
aC4mcm   = aC4Tmcm 

# dTERP/dt = wTERPs/h + wTERPe/h + RTERP
awTERPs    = aFterp_meg/ahm
awTERPe    = (awem * (aTERP_ftm - aTERPm))/ahm
adTERP_dt  = aTERPm_l + awTERPs[1:-1] + awTERPe[1:-1]
adTERP_dt2 = (aTERPm[0:-1-1] - aTERPm[1:-1]) / (atm[0:-1-1] - atm[1:-1])

# dISO/dt = wISOs/h + wISOe/h + RISO
awISOs    = aFiso_meg/ahm
awISOe    = (awem * (aISO_ftm - aISOm))/ahm
adISO_dt  = aISOm_l + awISOs[1:-1] + awISOe[1:-1]
adISO_dt2 = (aISOm[0:-1-1] - aISOm[1:-1]) / (atm[0:-1-1] - atm[1:-1])

#####     calculations    #####################################################################################

# replace negative nrs. with NaN
aH_null  = copy(aH)
aLE_null = copy(aLE)
for i in range(0,size(aH_null,0)):
   if aH_null[i] == 9999:
      aH_null[i] = NaN
   if aLE_null[i] == 9999:
      aLE_null[i] = NaN   
   if aT[i] >= 9999:
      aT[i] = NaN      
   if aRH[i] >= 9999:
      aRH[i] = NaN 
      
aEF = aLE_null/(aH_null+aLE_null)

for i in range(0,size(aEF,0)):
   if aEF[i] <0. or aEF[i] >1:
      aEF[i] = NaN               
for i in range(0,size(aP,0)):
   if aP[i] >= 9999:
      aP[i] = NaN
for i in range(0,size(at_AMS,0)):
   if aOA[i] >= 9999:
      aOA[i] = NaN  
      at_AMS[i] = NaN  
   if aSO4[i] >= 9999:
      aSO4[i] = NaN 
   if aNO3[i] >= 9999:
      aNO3[i] = NaN 
   if aNH4[i] >= 9999:
      aNH4[i] = NaN 
   if aChl[i] >= 9999:
      aChl[i] = NaN   
for i in range(0,size(at_FAGE,0)):
   if aOH[i] >= 1.e10:
      aOH[i] = NaN   
      at_OH[i] = NaN   
   if aHO2[i] >= 1.e10:
      aHO2[i] = NaN  
      at_HO2[i] = NaN
            
# remove NAN's (every 2nd element)
at_AMS2 = at_AMS[~isnan(aOA)]
aOA2    = aOA[~isnan(aOA)]
aSO42   = aSO4[~isnan(aOA)]
aNO32   = aNO3[~isnan(aOA)]
aNH42   = aNH4[~isnan(aOA)]
aChl2   = aChl[~isnan(aOA)]

at_AMS  = copy(at_AMS2)
aOA     = copy(aOA2)
aSO4    = copy(aSO42)
aNO3    = copy(aNO32)
aNH4    = copy(aNH42)
aChl    = copy(aChl2)
aTOT    = sum([aOA,aSO4,aNO3,aNH4,aChl],0) 

at_OH2   = at_OH[~isnan(aOH)]
at_HO22  = at_HO2[~isnan(aHO2)]
aOH2     = aOH[~isnan(aOH)]
aHO22    = aHO2[~isnan(aHO2)]
at_OH    = copy(at_OH2)
aOH      = copy(aOH2)
at_HO2   = copy(at_HO22)
aHO2     = copy(aHO22)

for i in range(0,size(at_PTRMS,0)):
   if aISO[i] >= 9999:
      aISO[i] = NaN   
   if aFISO[i] >= 9999:
      aFISO[i] = NaN   
   if aTERP[i] >= 9999:
      aTERP[i] = NaN   
   if aFTERP[i] >= 9999:
      aFTERP[i] = NaN     
      
for i in range(0,size(at_O3,0)):
   if aO3[i] >= 9999:
      aO3[i] = NaN   
   if aNO[i] >= 9999:
      aNO[i] = NaN   
   if aNO2[i] >= 9999:
      aNO2[i] = NaN   
      
for i in range(0,size(at_NOx_5,0)):
   if aNO_5[i] >= 9999:
      aNO_5[i] = NaN   
   if aNO2_5[i] >= 9999:
      aNO2_5[i] = NaN         
                       
# interpolation pressure observations
at_Pip   = at_meteo
f_P      = interp1d(at_P, aP, kind='linear')
aPip     = f_P(at_meteo)
aP       = copy(aPip)

#theta & q & VPD from T, RH & p
atheta   = aT*(1000./aP)**0.286 
aq       = ((aRH/100.)*(0.622/(aP))*(6.1078*exp(17.2693882*(aT-273.15)/((aT-273.15)+237.3))))*1000.
aqsat    = aq/(aRH/100.)
aVPD     = aqsat-aq

at_oaho  = mod(at_PMF, 1.0)*24.
at_amho  = mod(at_AMS, 1.0)*24.
at_meho  = mod(at_meteo, 1.0)*24.
at_gcho  = mod(at_O3, 1.0)*24.  
at_voho  = mod(at_PTRMS, 1.0)*24.
at_ohho  = mod(at_OH, 1.0)*24.
at_hoho  = mod(at_HO2, 1.0)*24.
at_o5ho  = mod(at_O3_5, 1.0)*24.
at_n5ho  = mod(at_NOx_5, 1.0)*24
at_bvho  = mod(at_BVOC, 1.0)*24.

# select day 189: 7 July 2008
at200    = at_meho[721:768]
aH       = aH_null[721:768]
aLE      = aLE_null[721:768]
aT       = aT[721:768]
atheta   = atheta[721:768]
aq       = aq[721:768]
aRH      = aRH[721:768]

at_gcho  = at_gcho[1066:1161]
aO3      = aO3[1066:1161]
aO3_30   = aO3_30[1066:1161]
aO3_45   = aO3_45[1066:1161]
aO3_60   = aO3_60[1066:1161]
aO3_75   = aO3_75[1066:1161]
aNO      = aNO[1066:1161]
aNO_30   = aNO_30[1066:1161]
aNO_45   = aNO_45[1066:1161]
aNO_60   = aNO_60[1066:1161]
aNO_75   = aNO_75[1066:1161]
aNO2     = aNO2[1066:1161]
aNO2_30  = aNO2_30[1066:1161]
aNO2_45  = aNO2_45[1066:1161]
aNO2_60  = aNO2_60[1066:1161]
aNO2_75  = aNO2_75[1066:1161]
at_ohho  = at_ohho[172:334]   
aOH      = aOH[172:334]
at_hoho  = at_hoho[180:344]   
aHO2     = aHO2[180:344] 
at_oaho  = at_oaho[1713:1853]
aOOA2    = aOOA2[1713:1853]
aOOA1    = aOOA1[1713:1853]
a82fac   = a82fac[1713:1853]
a91fac   = a91fac[1713:1853]
at_voho  = at_voho[720:767]
aISO     = aISO[720:767] 
aTERP    = aTERP[720:767] 
aFISO    = aFISO[720:767]
aFTERP   = aFTERP[720:767]
at_amho  = at_amho[1730:1872]
aOA      = aOA[1730:1872]
aSO4     = aSO4[1730:1872]
aNO3     = aNO3[1730:1872]
aNH4     = aNH4[1730:1872]
aChl     = aChl[1730:1872]
aTOT     = aTOT[1730:1872]
at_o5ho  = at_o5ho[17529:18968]
aO3_5    = aO3_5[17529:18968]
at_n5ho  = at_n5ho[1728:1871]
aNO_5    = aNO_5[1728:1871]
aNO2_5   = aNO2_5[1728:1871]
at_bvho  = at_bvho[720:767]
aMVK     = aMVK[720:767]

###############################################################################################
#####################    plotting     ##############################

figure(1, figsize=(11.5,9)) # dynamics: wths, wqs, h, theta, q
subplot(221) #
plot(atm, aHm, 'k-', at2m, aLEm, 'b--', at200, aH, 'k.', at200, aLE, 'b+')
xlabel('time LT (h)', fontsize=16)
ylabel(u'heat flux (W m\u207B\u00B2)', fontsize=16)
xlim([0, 24])
text(2,1050,'(a)', fontsize=12)
grid(True)
legend(('H', 'LE'), 'best')
subplot(222) # h
plot(atm, ahm ,'k-')#,at_lidar,ahlidar,'k*')
xlabel('time LT (h)', fontsize=16)
ylabel('h (m)', fontsize=16)
xlim([0, 24])
text(2,920,'(b)', fontsize=12)
grid(True)
subplot(223) # theta
plot(atm,athetam,'k-', at200, atheta, 'k.', lw=1)
xlabel('time LT (h)', fontsize=16)
ylabel(u'<\u03B8> (K)', fontsize=16)
xlim([0, 24])
text(2,305.2,'(c)', fontsize=12)
grid(True)
subplot(224) #
plot(atm, aqm ,'k-', at200, aq, 'k.', lw=1)
xlabel('time LT (h)', fontsize=16)
ylabel(u'<q> (g kg\u207B\u00B9)', fontsize=16)
xlim([0, 24])
text(2,14.1,'(d)', fontsize=12)
grid(True)

fig = figure(2, figsize=(11.5,9)) # O3, NOx, HOx
subplot(221) #
plot(atm, aO3m,'k-',at_gcho,aO3,'k.')#,at200, aO3ip,'r.')
ylabel(u'<O\u2083> (ppb)', fontsize=16)
xlim([6,14])
text(6.3,19.2,'(a)', fontsize=12)
grid(True)

ax1 = fig.add_subplot(222)
ax1.plot(at3m, aOHm,'k-',at_ohho,aOH,'k.')
gca().yaxis.set_major_formatter(OldScalarFormatter())
ax1.set_ylabel(u'<OH> (molec. cm\u207B\u00B3)', fontsize=16)
ax1.grid(True)
for tl in ax1.get_yticklabels():
    tl.set_color('k')

ax2 = ax1.twinx()
ax2.plot(at3m,aHO2m,'b-',at_hoho,aHO2,'b+')
gca().yaxis.set_major_formatter(OldScalarFormatter())
ax2.set_ylabel(u'<HO\u2082> (molec. cm\u207B\u00B3)', fontsize=16, color='b')
for tl in ax2.get_yticklabels():
    tl.set_color('b')
xlim([6,14])
text(6.3,3.7e8,'(b)', fontsize=12)

subplot(223) #
plot(at3m, aNOm,'k-',at_gcho,aNO,'k.')#, at200, aNOip,'r.')
ylabel('<NO> (ppb)', fontsize=16)
xlim([6,14])
xlabel('time LT (h)', fontsize=16)
text(6.3,0.22,'(c)', fontsize=12)
grid(True)

subplot(224) #
plot(at3m, aNO2m,'k-',at_gcho,aNO2,'k.')#,at200,aHO2ip,'r.')
ylabel(u'<NO\u2082> (ppb)', fontsize=16)
xlabel('time LT (h)', fontsize=16)
xlim([6,14])
ylim([0,0.4])
text(6.3,0.36,'(d)', fontsize=12)
grid(True)


fig = figure(3, figsize=(11.5,9)) # VOCs + OA
ax1 = fig.add_subplot(221)
ax1.plot(atm,aFISOm,'k-',at_voho,aFISO,'k.')
ylim([0,4])
ax1.set_ylabel(u'F_ISO (mg m\u207B\u00B2 h\u207B\u00B1)', fontsize=16)
ax1.grid(True)
for tl in ax1.get_yticklabels():
    tl.set_color('k')
ax2 = ax1.twinx()
ax2.plot(atm,aFTERPm,'b-',at_voho,aFTERP,'b+')
ax2.set_ylabel(u'F_TERP (mg m\u207B\u00B2 h\u207B\u00B1)', fontsize=16, color='b')
for tl in ax2.get_yticklabels():
    tl.set_color('b')
xlim([6,14])
ylim([0,1])

ax1 = fig.add_subplot(222)
ax1.plot(atm,aISOm,'k-',at_voho,aISO,'k.')
ylim([0,4])
ax1.set_ylabel('<ISO> (ppb)', fontsize=16)
ax1.grid(True)
for tl in ax1.get_yticklabels():
    tl.set_color('k')
ax2 = ax1.twinx()
ax2.plot(atm,aTERPm,'b-',at_voho,aTERP,'b+')
ax2.set_ylabel(u'<TERP> (ppb)', fontsize=16, color='b')
for tl in ax2.get_yticklabels():
    tl.set_color('b')
xlim([6,14])
ylim([0,0.4])

subplot(223) #
plot(at3m,aMVKm,'k-',at_bvho,aMVK,'k.')
ylabel(u'<MVK> (ppb)', fontsize=16)
xlabel('time LT (h)', fontsize=16)
xlim([6,14])
grid(True)

subplot(224)
plot(atm,aCoam-0.58,'k-',at_oaho,aOOA2,'k.')
ylabel('<Coa> (ug m-3)', fontsize=16)
xlabel('time LT (h)', fontsize=16)
xlim([6,14])
legend(('TSOA + ISOA','OOA2'),'upper left')
grid(True)


figure(6, figsize=(11.5,9))
title('OP3 7 July 2008')
subplot(221)
plot(atm,aFiso_meg*fcISO,'k-',at_voho,aFISO,'k.')#,at_oaho,aOOA2,'k.')
ylabel('F_ISO (mg m-2 h-1)', fontsize=16)
xlabel('time LT (h)', fontsize=16)
xlim([6,14])
legend(('MEGAN','obs.'),'upper left')
grid(True)

subplot(222)
plot(atm,aFterp_meg*fcTERP,'k-',at_voho,aFTERP,'k.')#,at_oaho,aOOA2,'k.')
ylabel('F_TERP (mg m-2 h-1)', fontsize=16)
xlabel('time LT (h)', fontsize=16)
xlim([6,14])
legend(('MEGAN','obs.'),'upper left')
grid(True)


'''
fig = figure(4, figsize=(11.5,9)) # sensitivity subsidence & advection
subplot(221)
plot(atm,aCoam-0.58,'-',at_oaho,aOOA2,'k.')
ylabel(u'C_OA (\u00B5g m\u207B\u00B3)', fontsize=16)
xlim([6,14])
text(6.3,0.12,'(a)', fontsize=12)
grid(True)
subplot(222)
plot(atm,ahm,'-')
ylabel('h (m)', fontsize=16)
xlim([6,14])
text(6.3,1030,'(b)', fontsize=12)
grid(True)
subplot(223)
plot(atm,athetam,at200,atheta,'k.')
xlabel('time LT (h)', fontsize=16)
ylabel(u'<\u03B8> (K)', fontsize=16)
text(6.3,301.3,'(c)', fontsize=12)
xlim([6,14])
grid(True)
subplot(224)
plot(atm,aqm,at200,aq,'k.')
xlabel('time LT (h)', fontsize=16)
ylabel(u'<q> (g kg\u207B\u00B9)', fontsize=16)
text(6.3,13.2,'(d)', fontsize=12)
xlim([6,14])
ylim([11,14])
grid(True)
'''
'''
fig = figure(5, figsize=(11.5,9)) # OH-recycling
subplot(221)
plot(atm,aCoam-0.58,'-',at_oaho,aOOA2,'k.')
ylabel(u'C_OA (\u00B5g m\u207B\u00B3)', fontsize=16)
xlim([6,14])
text(6.3,0.22,'(a)', fontsize=12)
grid(True)
subplot(222)
plot(atm,aOHm,'-',at_ohho,aOH,'k.')
gca().yaxis.set_major_formatter(OldScalarFormatter())
ylabel('<OH> (molec cm-3)', fontsize=16)
xlim([6,14])
text(6.3,7e6,'(b)', fontsize=12)
grid(True)
subplot(223)
plot(atm,aISOm,'-',at_voho,aISO,'k.')
ylabel('<ISO> (ppb)', fontsize=16)
xlabel('time LT (h)', fontsize=16)
xlim([6,14])
text(6.3,3.7,'(c)', fontsize=12)
grid(True)
subplot(224)
plot(atm,aTERPm,'-',at_voho,aTERP,'k.')
ylabel('<TERP> (ppb)', fontsize=16)
xlabel('time LT (h)', fontsize=16)
xlim([6,14])
text(6.3,0.37,'(d)', fontsize=12)
grid(True)
'''
figure(7, figsize=(8.3,6.1)) # VOC-tendency
subplot(211)
plot(atm[1:-1],adTERP_dt2,'r-',atm, awTERPs*3600, 'b-', atm, awTERPe*3600, 'm-', atm_pl,aTERPm_l*3600,'g',atm[1:-1],adTERP_dt*3600,'k--',)
ylabel(u'dTERP/dt (ppb h-1)', fontsize=16)
legend(('total','emission','entrainment','chemistry','sum'), 'lower left')
xlim([6,14])
grid(True)

subplot(212)
plot(atm[1:-1],adISO_dt2,'r-',atm, awISOs*3600, 'b-', atm, awISOe*3600, 'm-', atm_pl,aISOm_l*3600,'g',atm[1:-1],adISO_dt*3600,'k--',)
ylabel(u'dISO/dt (ppb h-1)', fontsize=16)
legend(('total','emission','entrainment','chemistry','sum'), 'lower left')
xlim([6,14])
grid(True)

show()

#awTERPs    = aFterp_meg/ahm
#awTERPe    = (awem * (aTERP_ftm - aTERPm))/ahm
#adTERP_dt  = aTERPm_l + awTERPs[1:-1] + awTERPe[1:-1]
#adTERP_dt2 = (aTERPm[0:-1-1] - aTERPm[1:-1]) / (atm[0:-1-1] - atm[1:-1])
