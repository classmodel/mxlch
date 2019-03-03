program bulk_model
use modchem
! USE DFLIB
implicit none

! --------------------------------------------
! Defines platform (MS Windows vs Linux/OSX)
  logical              :: windows = .false.
! --------------------------------------------

! From Tennekes and Driedonks,
! Boundary Layer Meteorology (1981), 515-531
!
! variables of the system
! thetam : potential temperature mixing layer value
! dtheta: potential temperature jump
! zi: boundary layer height
! -- wind ------------------------------------------
! um : initial horizontal windspeed ( x - direction)
! vm : initial horizontal windspeed  ( y - direction)
! um0 : initial horizontal windspeed ( x - direction)
! vm0 : initial horizontal windspeed  ( y - direction)
! ug: geostrophic wind in x-direction
! vg: geostrophic wind in y-direction
! --Specific humidity-------------------------
! qm : specific humidity mixing layer value
! dq : specific humidity jump
! zi : boundary layer height
! ws : Subsidence velocity
! prescribed value
!  --Potential temperature-------------------------
! beta=-(wthetave/wthetavs) : entrainment to surface flux ratio
! wthetas : surface heat flux
! gamma: lapse rate
!  --Specific humidity-------------------------
! betaq=-(wqe/wqs) : entrainment to surface flux ratio
! wq : surface humidity flux
! gammaq: lapse rate specific moisture
! -- wind ------------------------------------------
! um : initial horizontal windspeed ( x - direction)
! vm : initial horizontal windspeed  ( y - direction)
!
! initial value
!  --Potential temperature-------------------------
! z0: initial boundary layer height
! thetam0: initial potential temperature mixing layer value
! dtheta0: initial potential temperature jump
!  --Specific humidity-------------------------
! qm0: initial specific humidity mixing layer value
! dq0: initial specific humidity jump
!  --Carbon dioxide----------------------------------
! cm0: initial carbon dioxide mixing layer value
! dc0: initial carbon dioxide jump
! -- wind ------------------------------------------
! um0 : initial horizontal windspeed ( x - direction)
! vm0 : initial horizontal windspeed  ( y - direction)
! ug: geostrophic wind in x-direction
! vg: geostrophic wind in y-direction

!
! constant for coding
! time = second
!
! solving using eulerien's way
!
! declaration
!	 dynamics
! beta=-(wthetae/wthetas) : entrainment to surface flux ratio
! wthetas : surface heat flux
! gamma: lapse rate
!  --Specific humidity-------------------------
! betaq=-(wqe/wqs) : entrainment to surface flux ratio
! wq : surface humidity flux
! gammaq: lapse rate specific moisture
! gammac
!
!
! solving using eulerien's way
!
! Reaction_1   O3 + NO -> NO2 + O2          4.75E-4
! Reaction_2   HO + CO -> HO2 + CO2         6.0E-3
! Reaction_3   HO + RH -> HO2 + products    1.8
! Reaction_4   HO2 + NO -> HO + NO2         2.10E-1
! Reaction_5   HO2 + O3 -> HO + 2 O2        5.00E-5
! Reaction_6   2HO2 -> H2O2 + O2            7.25E-2
! Reaction_7   HO + NO2 -> HNO3             2.75E-1
! Reaction_8   HO + O3 -> HO2 + O2          1.75E-3
! Reaction_9   HO + HO2 -> H2O + O2         2.75
!
! Reaction_10  O3 -> 2 HO + O2              2.70E-6
! Reaction_11   NO2 -> NO + O3              8.90E-3
! declaration
!   dynamics

  logical :: c_ustr=.true.,c_wth=.false.,c_fluxes=.false., lencroachment=.false., ladvecFT=.false., lgamma=.false., lenhancedentrainment=.false., lfixedlapserates=.false., lfixedtroposphere=.false. !c_fluxes replaces c_wth
  logical :: lradiation=.false.,lsurfacelayer=.false.,llandsurface=.false.,lrsAgs=.false., lCO2Ags=.false.
  double precision :: zi(2),zi0 = 200 ,thetam(2), dtheta(2),pressure = 1013.0, wthetae
  double precision temp_cbl, temp_ft
  integer :: runtime,t, time=24*3600.0,tt
  double precision :: beta = 0.2 ,wthetas=0.0,gamma = 0.006,thetam0 = 295,dtheta0 = 4,wthetav=0.0,dthetav
  real :: dtime = 1
  double precision :: z0 = 0.03, kappa, zp, alpha,z0m=0.03,z0h=0.03, hcrit = 10000, gamma2 = 0.006
!! ROUGHNESS LENGTH
!! Terrain Description                                     ZO  (m)
!! Open sea, fetch at least 5km                            0.0002
!! Open flat terrain; grass, few isolated obstacles        0.03
!! Low crops, occasional large obstacles; L/h > 20         0.10
!! High crops, scattered obstacles,  15 < L/h < 20         0.25
!! Parkland, bushes, numerous obstacles,  L/h < 10         0.50
!! Regular large obstacle coverage (suburb, forest)        0.50 - 1.0
!! http://www.webmet.com/met_monitoring/663.html
!!
  real isec,mins,ttt, sec
  integer :: day = 80
  real :: hour = 0
  integer ierr,iter
  double precision pi

  double precision :: wthetasmax = 0.0 , wqsmax = 0.0 , wcsmax = 0.0
 
  double precision :: uws=0,vws=0,uws0 = 0.0 , vws0 = 0.0
  double precision ustar
  double precision :: um(2), vm(2), um0 = 0.0 , vm0 = 0.0, ueff, wstar
  double precision :: ug = 0.0 , vg = 0.0
  double precision du(2), dv(2)
  double precision :: ws,wsls=0.0, wsubs = 0.0 !wsubs = subsiding motions (subsidence + convective mass flux compensation) affecting inversions through the lapse rates
  double precision uwe,vwe
  double precision :: f,gammau = 0.0 , gammav = 0.0
  double precision we
  double precision :: advq = 0.0 ,lsq              ! large scale advection moisture (units (g/Kg)/s)!
  double precision :: advtheta = 0.0 ,lstheta, advtheta_now, advq_now    ! large scale advection heat (units K/s)!
  double precision :: cc=0.0,Qtot=400.0,albedo=0.2 ! cloud cover and incoming energy,albedo, incoming radiation
  double precision :: DeltaF = 0.0, Rdistr = 1.0   ! Radiative gradient due to aerosols and clouds (neg. = absorption, pos. = emission): amount (W/m2) and distribution (1 = top, 0 = over entire mixed layer)
  double precision :: DeltaFsw=0.0, DeltaFlw=0.0   ! Shortwave component (e.g. aerosols) and longwave component (e.g. clouds)
  double precision :: wf                           ! Resulting boundary-layer development (following the lapse rate)
  double precision :: Ts !Skin temperature
  double precision :: thetasurf, qsurf, thetavsurf, thetav, Rib, L, L0, fx, Lstart, Lend, fxdif
  double precision :: T2m, q2m, u2m, v2m, esat2m , e2m, rh2m
  double precision :: Tr,Ta,costh !Needed for radiation calculation
  double precision :: Swin,Swout,Lwin,Lwout !Calculated radiations
  double precision :: Cs = -1, Constm
  double precision :: esatsurf,qsatsurf,cq,rs=1.e6,ra,zsl,esat,qsat,rssoil
  double precision :: desatdT, dqsatdT, efinal, f1, f2, f3, f4, C1, C2, wgeq
  double precision :: w2=0.42,wfc=0.491,wwilt=0.314,gD=0.0,rsmin=0, LAI = 1.0
  double precision :: wsat=0.6,CLa=0.083, CLb=11.4, CLc=12.0, C1sat=0.342, C2ref=0.3
  double precision :: wg=0.40,rssoilmin=0,cliq,Wlmx,Wl=0.0,Wmax=2.0e-4,cveg=1.0
  double precision :: Lambda = 5.9, Tsoil=285, T2=285, Tsoiltend, Wltend, CGsat=3.6e-6
  double precision :: LEveg, LEliq, LEsoil, LE, SH, GR, CG, wgtend
  double precision :: tau_ta !advection relaxation time term (units s for Nagori 2019)

  double precision qm(2), dq(2),wqe
  double precision :: betaq,wqs=0,gammaq = 0.0 ,qm0 = 0.0,dq0 =0.0

  double precision :: cm(2), dc(2), wce, CO2ags, CO2comp, CO2comp298=68.5, Q10CO2=1.5
  double precision :: gm298 = 7, Ammax298 = 2.2, Q10gm = 2, T1gm = 278, T2gm = 301 
  double precision :: Q10Am = 2, T1Am = 281, T2Am = 311, f0 = 0.89, ad = 0.07, cfrac, co2abs, Ammax, betaw, PAR
  double precision :: alpha0 = 0.017,Kx = 0.7, gmin = 2.5e-4, gm, fmin0, nuco2q=1.6, fmin, Ds, D0, ci, fstr, Am, Rdark
  double precision :: alphac, tempy, An, E1, AGSa1, Dstar, gcco2,rsAgs,rsCO2, Resp, fw, wco2
  double precision :: factorial
  double precision :: Cw = 1.6e-3,wsmax=0.55,wsmin=0.005,R10=0.23,Eact0=53.3e3
  double precision :: betac,wcs,gammac = 0.0,cm0 = 0.0,dc0 = 0.0

  logical          :: lsea=.false.  ! sea surface fluxes switch
  double precision :: sst=285.0     ! sea surface temperature


  ! variables for Shallow cumulus 
  logical          :: lscu=.false.,lrelaxdz=.false.
  double precision :: tau=7200,dz=200,ca=0.5        ! if(lrelaxdz): dz relaxes as ddz/dt = 1/tau * zlcl-h
  double precision :: q2=0,ac=0,wm=0,wqm=0
  double precision :: sig_sc2=0, wscm=0, dsc=0, wsce=0      
  double precision :: Ptop,Ttop,estop,etop,qstop,qqs
  double precision :: ev,tempd,templcl,zlcl,RHlcl
  integer          :: ii
  double precision, allocatable :: wscm2(:)

  ! IEPOX uptake resistor model as per Gaston 2014 and Hu 2016 (Nagori 2019)
  real, parameter :: alpha_gmaIPX = 0.1
  
  real :: Rgas_g, Rp, epsln, Hconc, aH, Cnuc, Cga, Dgas, Daq, Haq, kH, knuc, kga, Dorg, SurfaceArea
  double precision :: k_aq, omega, HorgDorg, IEPOX_SOA, Rc_gIpx
  double precision :: k_org, q_aq, q_org, q_org_new, num, denom, h, F_gIpx
  double precision :: inv_gas, inv_Torq_aq, inv_gamma, inv_alpha, gamma_uptake, k_het
 

  

 ! define variables for G/P-partitioning, added a-pinene, b-pinene and limonene per oxidant (Nagori 2019)
  logical          :: lvbs=.false. 
  integer m, low_high_NOx
  real :: alpha1_APOHlow, alpha2_APOHlow, alpha3_APOHlow,   &
     alpha4_APOHlow, alpha1_BPOHlow, alpha2_BPOHlow, alpha3_BPOHlow,   &
     alpha4_BPOHlow, alpha1_APO3low, alpha2_APO3low, alpha3_APO3low,   &
     alpha4_APO3low, alpha1_BPO3low, alpha2_BPO3low, alpha3_BPO3low,   &
     alpha4_BPO3low, alpha1_APNO3, alpha2_APNO3, alpha3_APNO3,   &
     alpha4_APNO3, alpha1_BPNO3, alpha2_BPNO3, alpha3_BPNO3,   &
     alpha4_BPNO3, alpha1_APOHhigh, alpha2_APOHhigh, alpha3_APOHhigh, &
     alpha4_APOHhigh, alpha1_BPOHhigh, alpha2_BPOHhigh,   &
     alpha3_BPOHhigh, alpha4_BPOHhigh, alpha1_APO3high,   &
     alpha2_APO3high, alpha3_APO3high, alpha4_APO3high,   &
     alpha1_BPO3high, alpha2_BPO3high, alpha3_BPO3high,   &
     alpha4_BPO3high,   &
     alpha1_LIOHlow, alpha2_LIOHlow, alpha3_LIOHlow, alpha4_LIOHlow,   &
     alpha2_ISO3low, alpha3_ISO3low, alpha4_ISO3low, alpha1_LIO3low,   &
     alpha2_LIO3low, alpha3_LIO3low, alpha4_LIO3low, alpha2_ISNO3,   &
     alpha3_ISNO3, alpha4_ISNO3, alpha1_LINO3, alpha2_LINO3,  &
     alpha3_LINO3, alpha4_LINO3,  &
     alpha1_LIOHhigh, alpha2_LIOHhigh,   &
     alpha3_LIOHhigh, alpha4_LIOHhigh, alpha2_ISO3high, alpha3_ISO3high,   &
     alpha4_ISO3high, alpha1_LIO3high, alpha2_LIO3high, alpha3_LIO3high,   &
     alpha4_LIO3high
  real :: APOHmc, APO3mc, APNO3mc, BPOHmc, BPO3mc, BPNO3mc, &
    LIOHmc, LIO3mc, LINO3mc, ISO3mc, ISNO3mc, beta_iso, beta_terp
  real :: epsi, C1Tmc, C2Tmc, C3Tmc, C4Tmc, C1_st, C2_st, C3_st, C4_st, X1, &
      X2, X3, X4, T_hh, Coait, Coait_new, Coa, OAbgb
  real :: C1AHlmc,C1AHhmc,C1BHlmc,C1BHhmc,C1LHlmc,C1LHhmc,C1AOlmc,C1AOhmc,C1BOlmc,C1BOhmc, &
          C1LOlmc,C1LOhmc,C1Anmc,C1BNmc,C1LNmc,C2AHlmc,C2AHhmc,C2BHlmc,C2BHhmc,C2AOlmc,C2AOhmc, &
          C2BOlmc,C2BOhmc,C2LHlmc,C2LHhmc,C2IOlmc,C2IOhmc,C2LOlmc,C2LOhmc,C2ANmc, &
          C2BNmc,C2INmc,C2LNmc,C3AHlmc,C3AHhmc,C3BHlmc,C3BHhmc,C3AOlmc,C3AOhmc,C3BOlmc,C3BOhmc, &
          C3LHlmc,C3LHhmc,C3IOlmc,C3IOhmc,C3LOlmc,C3LOhmc,C3ANmc,C3BNmc,C3INmc, &
          C3LNmc,C4AHlmc,C4AHhmc,C4BHlmc,C4BHhmc,C4AOlmc,C4AOhmc,C4BOlmc,C4BOhmc, &
          C4LHlmc,C4LHhmc,C4IOlmc,C4IOhmc,C4LOlmc,C4LOhmc,C4ANmc,C4BNmc,C4INmc,C4LNmc
  real :: C1AHlmc_ft,C1AHhmc_ft,C1BHlmc_ft,C1BHhmc_ft,C1LHlmc_ft,C1LHhmc_ft,C1AOlmc_ft,C1AOhmc_ft, &
          C1BOlmc_ft,C1BOhmc_ft,C1LOlmc_ft,C1LOhmc_ft,C1ANmc_ft,C1BNmc_ft,C1LNmc_ft,C2AHlmc_ft, &
          C2AHhmc_ft,C2BHlmc_ft,C2BHhmc_ft,C2AOlmc_ft,C2AOhmc_ft,C2BOlmc_ft,C2BOhmc_ft, &
          C2LHlmc_ft,C2LHhmc_ft,C2IOlmc_ft,C2IOhmc_ft,C2LOlmc_ft,C2LOhmc_ft,C2ANmc_ft, &
          C2BNmc_ft,C2INmc_ft,C2LNmc_ft,C3AHlmc_ft,C3AHhmc_ft,C3BHlmc_ft,C3BHhmc_ft,C3AOlmc_ft, &
          C3AOhmc_ft,C3BOlmc_ft,C3BOhmc_ft,C3LHlmc_ft,C3LHhmc_ft,C3IOlmc_ft, &
          C3IOhmc_ft,C3LOlmc_ft,C3LOhmc_ft,C3ANmc_ft,C3BNmc_ft,C3INmc_ft,C3LNmc_ft,C4AHlmc_ft, &
          C4AHhmc_ft,C4BHlmc_ft,C4BHhmc_ft,C4AOlmc_ft,C4AOhmc_ft,C4BOlmc_ft,C4BOhmc_ft, &
          C4LHlmc_ft,C4LHhmc_ft,C4IOlmc_ft,C4IOhmc_ft,C4LOlmc_ft,C4LOhmc_ft,C4ANmc_ft, &
          C4BNmc_ft,C4INmc_ft,C4LNmc_ft
  real :: C1mc, C2mc, C3mc, C4mc, C1mc_ft, C2mc_ft, C3mc_ft, C4mc_ft
  real :: APOHmc_ft, APO3mc_ft, APNO3mc_ft, BPOHmc_ft, &
          BPO3mc_ft, BPNO3mc_ft, LIOHmc_ft, LIO3mc_ft, LINO3mc_ft, &
          ISO3mc_ft, ISNO3mc_ft
  real :: SVOC_alp1, SVOC_alp2, SVOC_alp3, SVOC_alp4, SVOC_bep1, SVOC_bep2, SVOC_bep3, SVOC_bep4, SVOC_lim1, SVOC_lim2, SVOC_lim3, SVOC_lim4, SVOC_iso2, SVOC_iso3, SVOC_iso4, ISOOH, ISONO3
  real :: SVOC_alp1_ft, SVOC_alp2_ft, SVOC_alp3_ft, SVOC_alp4_ft, SVOC_bep1_ft, SVOC_bep2_ft, SVOC_bep3_ft, SVOC_bep4_ft, SVOC_lim1_ft, SVOC_lim2_ft, SVOC_lim3_ft, SVOC_lim4_ft, SVOC_iso2_ft, SVOC_iso3_ft, SVOC_iso4_ft
  real :: alp_OH1, alp_OH2, alp_OH3, alp_OH4, alp_O31, alp_O32, alp_O33, alp_O34, alp_NO31, alp_NO32, alp_NO33, alp_NO34, bep_OH1, bep_OH2, bep_OH3, bep_OH4
  real :: bep_O31, bep_O32, bep_O33, bep_O34, bep_NO31, bep_NO32, bep_NO34, lim_OH1, lim_OH2, lim_OH3, lim_OH4, lim_O31, lim_O32, lim_O33, lim_O34, lim_NO31, lim_NO32, lim_NO33, lim_NO34
  real :: isoo32, isoo33, isoo34, isono32, isono33, isono34
  real :: alp_OH1_ft, alp_OH2_ft, alp_OH3_ft, alp_OH4_ft, alp_O31_ft, alp_O32_ft, alp_O33_ft, alp_O34_ft, alp_NO31_ft, alp_NO32_ft, alp_NO33_ft, alp_NO34_ft, bep_OH1_ft, bep_OH2_ft, bep_OH3_ft, bep_OH4_ft
  real :: bep_O31_ft, bep_O32_ft, bep_O33_ft, bep_O34_ft, bep_NO31_ft, bep_NO32_ft, bep_NO34_ft, lim_OH1_ft, lim_OH2_ft, lim_OH3_ft, lim_OH4_ft, lim_O31_ft, lim_O32_ft, lim_O33_ft, lim_O34_ft, lim_NO31_ft, lim_NO32_ft, lim_NO33_ft, lim_NO34_ft
  real :: isoo32_ft, isoo33_ft, isoo34_ft, isono32_ft, isono33_ft, isono34_ft
 
  real :: epsi_ft, C1Tmc_ft, C2Tmc_ft, C3Tmc_ft, C4Tmc_ft, &
          C1_st_ft, C2_st_ft, C3_st_ft, C4_st_ft, X1_ft, X2_ft, X3_ft, &
          X4_ft, T_ft, OAbgb_ft, Coait_ft, Coait_new_ft, Coa_ft
  real :: mr2mc_OAbg, mr2mc_PRODT, mr2mc_PRODI, C10_st=0.1, &
          C20_st=1., C30_st=10., C40_st=100.
  real :: kISORO2NO, kISORO2HO2, kTERPRO2NO, kTERPRO2HO2
  real, parameter :: T0      = 298.    ! temperature (K) at which the lab experiments and fitting are done, Tsimpidi'10
  real, parameter :: dHvap   = 130.e3   ! enthalpy of vaporization (J/mol), Lane'08
  real, parameter :: PRODTmm  = 200.   ! molar mass terp-oxidation products (g/mol)
  real, parameter :: PRODImm  = 136.   ! molar mass iso-oxidation products (g/mol)
  real, parameter :: mmOAbg = 250.     ! molar mass background OA (g/mol)

  real conv_cbl,conv_ft, Rfact


  ! define variables for bVOC-emission
  logical          :: lBVOC =.false.
  double precision :: BaserateISO = 0.0, BaserateTER = 0.0
  double precision :: dwg = 0.06, gammaLAI, w1, gammasm
  double precision :: beta_MT, Tref, gammaT_MT, gammace_MT, gamma_MT, eps_MT, rho_MT, ER_MT, mm_MT, cf_MT, F_TERP 
  double precision :: Thr, Tdaily, Topt, C_T1, C_T2, x, Eopt, gammaT_ISO, aa, Pac, Pdaily, Ptoa, phi, gammap, gammace_ISO, gamma_ISO, eps_ISO, rho_ISO, ER_ISO, mm_ISO, cf_ISO, F_ISO
  
  !Define variables for the 'saturation level' program
  integer :: i,a,n,aver3,sat_lev
  double precision p0,epsilon,Tabs,Cp,Rd,Rv,Lv,e0,gammad,gammam
  double precision np,p,lcl,lcl0
  double precision z_s,T_s,qt_s,qs_s,es_s
  double precision z,Tair,qt,qs,es,e_saturation,Td
  double precision Tparcel_s,Tparcel
  double precision sat_time,sat_height

  integer :: atime = 60, aver1, itime
  double precision g, Cf,eta, CT, C_m

  !Define variables for vertical profiles
  integer ::  n_vert, atime_vert = 1800, aver2
  double precision :: h_max = 3000 , eps, inf
  
  !constants
!  double precision,parameter :: rvrd = 461.5/287.04  !  Rv/Rd
  real,parameter :: rvrd = 0.61     !Rd/Rv    
  real,parameter :: bolz =  5.67e-8 !Stefan-Boltzmann constant [-]
  real,parameter :: rhow =  1000.   !density of water [kg m-3]
  real,parameter :: rho  =    1.2   !density of air [kg m-3]
  real,parameter :: Rgas = 8.3145   ! gas constant (J mol-1 k-1)  
  real           :: S0   =  1368.   !Incoming shortwave radiation [W m-2]

  !chemistry
  integer k
  integer startdaytime,enddaytime,daylength,daystart,daytime_start,daytime_end,prevdaysec
  double precision photo
  real getth
  real psim,psih
  real factor, factor1, factor2, factor3, printhour
  real a0, a1, a2
  
  

  ! general
  character(len=25) inputchemfile
  character(len=80) formatstring
  character(len=80) formatstring2
  character(len=40) filename
  integer j
  character(len=2),allocatable :: dummy(:)
  character(len=25) :: outdir = 'RUN00'

  !Adapted surface fluxes
  !Times in s after start of run
  integer  starttime_wT , endtime_wT , starttime_wq, endtime_wq
  integer  starttime_chem, endtime_chem, starttime_adv_T, endtime_adv_T
  integer starttime_adv_q, endtime_adv_q
  !Offsets in standard flux units (K m s-1 & g kg-1 m s-1)
  double precision offset_wT, offset_wq
  !Functions of fluxes
  ! 0  =  no flux
  ! 1  =  constant flux
  ! 2  =  sinoid flux from starttime to endtime
  ! 3  =  constant flux froms tarttime to endtime
  integer function_wT, function_wq   

  character(len=1) dirsep
  character(len=5) kopie

  ! !!!!!!!!!!!!!!!!!!!!
  ! No longer used, see below declaration namelist
  ! BvS, nov2011
  ! !!!!!!!!!!!!!!!!!!!!
  !DEC$ IF DEFINED (LINUX)
  !  character(len=1) :: dirsep ='/'
  !  character(len=3) :: kopie = 'cp '
  !DEC$ ELSE !WINDOWS
  !  character(len=1) :: dirsep ='\'
  !  character(len=5) :: kopie = 'copy '
  !DEC$ ENDIF
  
  ! general options
  namelist/NAMRUN/ &
    outdir, &
    time, &
    dtime, &
    atime, &
    atime_vert, &
    h_max, &
    latt, &
    long, &
    day, &
    hour

  ! option for the dynamics
  namelist/NAMDYN/ &
    zi0, &
    beta, &
    lenhancedentrainment, &
    wsls, &
    lfixedlapserates, &
    lfixedtroposphere, &
    wthetasmax, &
    c_wth, &
    c_fluxes, &
    gamma, &
    lgamma, &
    hcrit,&
    gamma2,&
    thetam0, &
    dtheta0, &
    pressure, &
    wqsmax, &
    gammaq, &
    qm0, &
    dq0, &
    wcsmax, &
    gammac, &
    cm0, &
    dc0, &
    c_ustr, &
    z0, &
    uws0, &
    vws0, &
    gammau, &
    gammav, &
    um0, &
    vm0, &
    ug, &
    vg, &
    advq, &
    advtheta, &
    ladvecFT, &
    lencroachment, &
    lscu, &
    lrelaxdz, &
    tau

  namelist/NAMSURFLAYER/ &
    lsurfacelayer,&
    z0m,&
    z0h

  namelist/NAMRAD/ &
    lradiation,& !radiation scheme to determine Q and SW
    cc,&         !cloud cover
    S0,&         !Incoming radiation
    DeltaFsw,&   !Absorbed radiation by e.g. aerosols (neg. value)
    DeltaFlw,&   !Emitted radiation by e.g. clouds (pos. value)
    Rdistr,&     !Distribution of absorbing aerosols (see Barbaro et al., 2013)
    albedo       !Surface albedo

  namelist/NAMSURFACE/ &
    llandsurface,& !switch to use interactive landsurface
    Qtot,&         !Incoming energy
    lsea,&         !Using a sea surface instead of land
    sst,&          !Sea surface temperature
    Ts,&           !Initial surface temperature [K]
    wwilt,&        !wilting point
    w2,&           !Volumetric water content deeper soil layer
    wg,&           !Volumetric water content top soil layer
    wfc,&          !Volumetric water content field capacity
    wsat,&         !Saturated volumetric water content ECMWF config
    CLa,&          !Clapp and Hornberger retention curve parameter a
    CLb,&          !Clapp and Hornberger retention curve parameter b
    CLc,&          !Clapp and Hornberger retention curve parameter c
    C1sat,&        !Coefficient force term moisture
    C2ref,&        !Coefficient restore term moisture
    gD,&           !VPD correction factor for rs
    rsmin,&        !Minimum resistance of transpiration
    rssoilmin,&    !Minimum resistance of soiltranspiration
    LAI,&          !Leaf area index
    cveg,&         !Vegetation fraction
    Tsoil,&        !Temperature top soil layer
    T2,&           !Temperature deeper soil layer
    Wl,&           !Equivalent water layer depth for wet vegetation
    Lambda,&       !Thermal diffusivity skin layer
    CGsat,&        !Saturated soil conductivity for heat
    lrsAgs,&       !Switch to use A-gs model for surface resistances
    lCO2Ags,&      !Switch to use A-gs model for CO2 flux
    CO2comp298,&   !CO2 compensation concentration [mg m-3]
    Q10CO2,&       !function parameter to calculate CO2 compensation concentration [-]
    gm298   ,&     !mesophyill conductance at 298 K [mm s-1]
    Ammax298,&     !CO2 maximal primary productivity [mg m-2 s-1]
    Q10gm   ,&     !function parameter to calculate mesophyll conductance [-] 
    T1gm    ,&     !reference temperature to calculate mesophyll conductance gm [K]
    T2gm    ,&     !reference temperature to calculate mesophyll conductance gm [K]
    Q10Am   ,&     !function parameter to calculate maximal primary profuctivity Ammax
    T1Am    ,&     !reference temperature to calculate maximal primary profuctivity Ammax [K] 
    T2Am    ,&     !reference temperature to calculate maximal primary profuctivity Ammax [K]
    f0      ,&     !maximum value Cfrac [-]
    ad      ,&     !regression coefficient to calculate Cfrac [kPa-1]
    alpha0  ,&     !initial low light conditions [mg J-1]
    Kx      ,&     !extinction coefficient PAR [-]
    gmin    ,&     !cuticular (minimum) conductance [m s-1]
    Cw      ,&     !constant water stress correction (eq. 13 Jacobs et al. 2007) [-]
    wsmax   ,&     !upper reference value soil water [-]
    wsmin   ,&     !lower reference value soil water [-]
    R10     ,&     !respiration at 10 C [mg CO2 m-2 s-1]
    Eact0   ,&     !activation energy [53.3 kJ kmol-1]
    lBVOC   ,&     !Enable the calculation of BVOC (isoprene, terpene) emissions
    BaserateIso ,& !Base emission rate for isoprene emissions [microg m^4 h^-1]
    BaserateTer    !Base emission rate for terprene emissions [microg m^4 h^-1]

  ! option for the chemistry
  namelist/NAMCHEM/ &
    lchem, &
    lwritepl, &
    lcomplex, &
    ldiuvar,&
    h_ref  ,&
    lflux, &
    fluxstart, &
    fluxend, &
    pressure_ft ,&
    lchconst ,&
    t_ref_cbl ,&
    p_ref_cbl ,&
    q_ref_cbl ,&
    t_ref_ft  ,&
    p_ref_ft  ,&
    q_ref_ft

  ! options for changes in the surface fluxes
  namelist/NAMFLUX/ &
    starttime_wT  , &
    endtime_wT    , &
    offset_wT     , &
    starttime_wq  , &
    endtime_wq    , &
    offset_wq     , &
    starttime_chem, &
    endtime_chem  , &
    starttime_adv_T , &
    endtime_adv_T   , &
    starttime_adv_q , &
    endtime_adv_q   , &
    function_wT   , &
    function_wq   

  ! options for the gas/particle partitioning leading to SOA formation, using the Volatility Basis Set (VBS)
  ! terpene volatilty bins split into volatillity bins for alpha-pinene, beta-pinene, limonene and BVOC reactions with NO3 added
  ! (Nagori 2019)
  namelist/NAMSOA/ &
    lvbs          , &
    low_high_NOx  , &
    alpha1_APOHlow , &
    alpha2_APOHlow , &
    alpha3_APOHlow , &
    alpha4_APOHlow , &
    alpha1_BPOHlow , &
    alpha2_BPOHlow , &
    alpha3_BPOHlow , &
    alpha4_BPOHlow , &
    alpha1_APO3low , &
    alpha2_APO3low , &
    alpha3_APO3low , &
    alpha4_APO3low , &
    alpha1_BPO3low , &
    alpha2_BPO3low , &
    alpha3_BPO3low , &
    alpha4_BPO3low , &
    alpha1_APNO3 , &
    alpha2_APNO3 , &
    alpha3_APNO3 , &
    alpha4_APNO3 , &
    alpha1_BPNO3 , &
    alpha2_BPNO3 , &
    alpha3_BPNO3 , &
    alpha4_BPNO3 , &
    alpha1_APOHhigh , &
    alpha2_APOHhigh , &
    alpha3_APOHhigh , &
    alpha4_APOHhigh , &
    alpha1_BPOHhigh , &
    alpha2_BPOHhigh , &
    alpha3_BPOHhigh , &
    alpha4_BPOHhigh , &
    alpha1_APO3high , &
    alpha2_APO3high , &
    alpha3_APO3high , &
    alpha4_APO3high , &
    alpha1_BPO3high , &
    alpha2_BPO3high , &
    alpha3_BPO3high , &
    alpha4_BPO3high , &
    alpha1_LIOHlow , &
    alpha2_LIOHlow , &
    alpha3_LIOHlow , &
    alpha4_LIOHlow , &
    alpha2_ISO3low , &
    alpha3_ISO3low , &
    alpha4_ISO3low , &
    alpha1_LIO3low , &
    alpha2_LIO3low , &
    alpha3_LIO3low , &
    alpha4_LIO3low , &
    alpha2_ISNO3 , &
    alpha3_ISNO3 , &
    alpha4_ISNO3 , &
    alpha1_LINO3 , &
    alpha2_LINO3 , &
    alpha3_LINO3 , &
    alpha4_LINO3 , &
    alpha1_LIOHhigh , &
    alpha2_LIOHhigh , &
    alpha3_LIOHhigh , &
    alpha4_LIOHhigh , &
    alpha2_ISO3high , &
    alpha3_ISO3high , &
    alpha4_ISO3high , &
    alpha1_LIO3high , &
    alpha2_LIO3high , &
    alpha3_LIO3high , &
    alpha4_LIO3high

  if (windows) then
    dirsep = '\'
    kopie = 'copy '
  else
    dirsep ='/'
    kopie = 'cp '
  end if

  ! some variables
  inf     = 1.0e-9        ! Correction to avoid divide by 0 if z0=0
  eps     = 1.0e-4        ! Correction to avoid crashing with plots
                          !   in visual BASIC
  aver3 = 300             ! saturation level will be calculated every 300s

  pi=acos(-1.)

  Rfact= 8.314e-2 ! mbar*m3 /K*mol

  lcl = 3000.    !dummy value for first run
  a = 0
  n = 0
  g = 9.81         !m/s2
  Cf = 0.2         ! constants from D.Pino
  CT = 4.          ! J. At. Sci. Vol. 60
  eta = 2.         ! 1913 - 1926
  C_m = 0.7        !

  !    Define more constants for the 'saturation level'program

  p0 = 100.00      !kPa
  epsilon = 0.622  !g water / g air = Rd/Rv
  Tabs = 273.0     !K
  Cp = 1004.67     !J/ kg K
  Rd = 287.053     !J/ kg K
  Rv = 461.50      !J/ kg K
  Lv = 2.501*10**6 !J/kg
  e0 = 0.611       !kPa
  gammad = -(g/Cp) !K/m

  kappa = 0.4     ! Von Karman constant
  zp = 10.        ! Height for the calculation of the logarithmic equation for u* (m)

  ! Number of vertical points
  n_vert = 4

  ! Number of vertical levels of the saturation plot
  sat_lev= 50
  ! Fitting function isoprene flux Gaussian
    a0 = 0.653289
!   a0= 0.653289 a0 = 0.8775 (sim04) a0 = 0.42458 (sim05)
    a1 = 42705.1
    a2 = 7999.29 

  !with use DFLIB NARGS() visual fortran 6 on windows
  !DEC$ IF DEFINED (VISUAL_COMPILER)
  !  if ( NARGS ()<= 1 ) then
  !    inputchemfile = 'chem.inp'
  !  else
  !    call getarg(1,inputchemfile)
  !  endif
  !DEC$ ENDIF

  !New Intel compilers also take GET_COMMAND_ARGUMENT
  !DEC$ IF DEFINED (INTEL_COMPILER )
  !  if ( iargc ()< 1 ) then
  !    inputchemfile = 'chem.inp'
  !  else
  !    call getarg(1,inputchemfile,n)
  !  endif
  !DEC$ ENDIF

  !GNU fortran on Windows and Linux ?
  ! Works both with ifort and gnu fortran, don't know about visual
  ! BvS, nov2011
    if ( COMMAND_ARGUMENT_COUNT () > 0) then
      CALL GET_COMMAND_ARGUMENT(1,inputchemfile)
    else
      inputchemfile = 'chem.inp'
    endif

  open (1, file='namoptions')
  read (1,NAMRUN,iostat=ierr)
  close(1)
  open (1, file='namoptions')
  read (1,NAMDYN,iostat=ierr)
  if (ierr > 0) stop 'ERROR: Problem in namoptions'
  close(1)
  open (1, file='namoptions')
  read (1,NAMRAD,iostat=ierr)
  if (ierr > 0) stop 'ERROR: Problem in namoptions'
  close(1)

  pressure_ft = pressure
  Ts          = thetam0
  z0m         = z0
  z0h         = z0

  open (1, file='namoptions')
  read (1,NAMCHEM,iostat=ierr)
  if (ierr > 0) stop 'ERROR: Problem in namoptions'
  close(1)
  open (1, file='namoptions')
  read (1,NAMSURFACE,iostat=ierr)
  if (ierr > 0) stop 'ERROR: Problem in namoptions'
  close(1)
  open (1, file='namoptions')
  read (1,NAMSURFLAYER,iostat=ierr)
  if (ierr > 0) stop 'ERROR: Problem in namoptions'
  close(1)
  open (1, file='namoptions')
  read (1,NAMSOA,iostat=ierr)
  if (ierr > 0) stop 'ERROR: Problem in namoptions'
  close(1)

  if (lfixedtroposphere) lfixedlapserates = .true.
  if ( lCO2Ags ) then
    lrsAgs       = .true.
    llandsurface = .true.
  endif
  if ( .not. lchem ) lvbs = .false.
  if ( lrsAgs .and. (.not. llandsurface) ) then
    print *,""
    print *,""
    print *,"!!!!!!!!!!!!!!!!!!!!!"
    print *,""
    print *,"You enabled lrsAgs without enabling llandsurface!"
    print *,"Please check your namelist options."
    print *,""
    print *,"!!!!!!!!!!!!!!!!!!!!!"
    print *,""
    print *,""
  endif
  if ( lrsAgs .and. (.not. lsurfacelayer) ) then
    lsurfacelayer = .true.
    print *,""
    print *,""
    print *,"!!!!!!!!!!!!!!!!!!!!!"
    print *,""
    print *,"You enabled lrsAgs without enabling lsurfacelayer!"
    print *,"lsurfacelayer is automatically enabled"
    print *,""
    print *,"!!!!!!!!!!!!!!!!!!!!!"
    print *,""
    print *,""
  endif
  if ( lrsAgs .and. (.not. lradiation) ) then
    lradiation = .true.
    print *,""
    print *,""
    print *,"!!!!!!!!!!!!!!!!!!!!!"
    print *,""
    print *,"You enabled lrsAgs without enabling lradiation!"
    print *,"lradiation is automatically enabled"
    print *,""
    print *,"!!!!!!!!!!!!!!!!!!!!!"
    print *,""
    print *,""
  endif
  if (lchem .and. (dtime .gt. 2.)) then
    dtime = 2.
    print *,""
    print *,""
    print *,"!!!!!!!!!!!!!!!!!!!!!"
    print *,""
    print *,"You enabled lchem with dtime > 2 seconds!"
    print *,"This might result in unstable chemistry results"
    print *,"dtime set to 2 seconds"
    print *,""
    print *,"!!!!!!!!!!!!!!!!!!!!!"
    print *,""
    print *,""
  endif


  write(formatstring,'(a,a)')'mkdir ',trim(outdir)
  call system(formatstring)
  write(formatstring,'(a,a)')'mkdir ',trim(outdir)//dirsep//'PL'
  call system(formatstring)
  write(formatstring,'(a,a)')kopie//'namoptions '//trim(outdir)
  call system(formatstring)
  write(formatstring,'(a,a)')kopie//inputchemfile//trim(outdir)
  call system(formatstring)
  if (lcomplex) then
    write(formatstring,'(a,a)')kopie//'chemicals.txt '//trim(outdir)
    call system(formatstring)
  endif

  !estimate the daylength
  daystart=0
  do i=1,24*3600
    thour = i/3600.
      zenith = getth(1.0*day,latt,long,i/3600.)
    if (cos(zenith) > 0.0 .and. daystart == 0 ) then
      daystart=1
      startdaytime = i
    endif
    if (cos(zenith) <= 0.0 .and. daystart == 1 ) then
      enddaytime = i
      exit
    endif
  enddo
  daylength = enddaytime - startdaytime

  write (*,*)'LCHEM=',lchem,'LDIUVAR=',ldiuvar,'LCONST=',lchconst,'LFLUX=',lflux,'LVBS=',lvbs, 'LBVOC=', lBVOC
  write (*,*) ' long',long
  write (*,*) 'latt',latt
  write (*,*) ' day',day
  write (*,*) 'hour',hour
  write (*,*) 'daystart',startdaytime/3600.
  write (*,*) 'dayend  ',enddaytime/3600.
  write (*,*) 'daylength',daylength/3600.

  ! DEPRECATED: USE NAMELIST NAMFLUX to set individually
  ! if lfux then we use startflux and endflux times
  if (lflux .eqv. .true.) then
    startdaytime = fluxstart * 3600
    enddaytime  = fluxend * 3600
    daylength = enddaytime - startdaytime
    write (*,*)'**** LFLUX = TRUE  *****'
    write (*,*) 'daystart',startdaytime/3600.
    write (*,*) 'dayend  ',enddaytime/3600.
    write (*,*) 'daylength',daylength/3600.
  endif

  if ( hour*3600 > enddaytime )then   !we start in previous night
    prevdaysec    = nint((24 - hour) * 3600)
    daytime_start = nint((24 - hour) * 3600 + startdaytime) !seconds after starttime model
    daytime_end   = nint((24 - hour) * 3600 + enddaytime)
  else
    daytime_start = nint(startdaytime - hour * 3600) !seconds after starttime model
    daytime_end   = nint(enddaytime - hour *3600)
  endif

  write(*,*) 'Day time starts',daytime_start,'seconds after start model'
  write(*,*) 'Day time end',daytime_end,'seconds after start model'


  starttime_wT   = daytime_start 
  endtime_wT     = daytime_end
  starttime_wq   = daytime_start
  endtime_wq     = daytime_end
  starttime_chem = daytime_start
  endtime_chem   = daytime_end
  starttime_adv_T  = daytime_start
  endtime_adv_T    = daytime_end
  starttime_adv_q  = daytime_start
  endtime_adv_q    = daytime_end
  offset_wT      = 0.0
  offset_wq      = 0.0
  !!Function:    just like chemistry functions:
  !!             0  =  no flux
  !!             1  =  constant flux
  !!             2  =  sinoid flux from starttime to endtime
  !!             3  =  constant flux froms tarttime to endtime
  function_wT    = 2 
  function_wq    = 2

  open (1, file='namoptions',iostat=ierr)
  read (1,NAMFLUX,iostat=ierr)
  if (ierr > 0) stop 'ERROR: Problem in namoptions'
  close(1)

  ! initialisation dynamics

  if (lencroachment) then
    beta    = 0.0
    dtheta0 = 0.0
    dq0     = 0.0
  else
    dthetav = dtheta0 + 1.e-03*rvrd*(qm0*dtheta0+thetam0*dq0+dtheta0*dq0)
    if (dthetav .le. 0.0) then
      write(*,*) 'Initial dthetav <= 0, but encroachment switch not activated!'
      write(*,*) 'dq is set to 0.0 and dtheta to 0.01 K'
      dtheta0 = 0.01
      dq0     = 0.0
    endif
    if (beta .le. 0.0) then
      write(*,*) 'WARNING: beta is set to ',beta
      write(*,*) 'This could result in unrealistic situations (negative inversions)'
    endif
  endif

  f = 2.*pi/(24.*3600.)*2.*sin(pi*latt/180.) ! Coriolis parameter

  runtime=nint(time/dtime)

  zi(1)=zi0

  thetam(1)=thetam0
  dtheta(1)=dtheta0

  um(1)= um0
  vm(1)= vm0
  du(1)= ug-um0
  dv(1)= vg-vm0

  qm(1) = qm0
  dq(1) = dq0

  cm(1) = cm0
  dc(1) = dc0

  DeltaF = DeltaFsw + DeltaFlw

  aver1=atime/dtime
  aver2=atime_vert/dtime

  ! initialise chemistry
  !
  allocate(dummy(mrpcc))
  if (lchem) then
    if (lcomplex) then
      call inputchem_mozart(inputchemfile,outdir,dirsep)
    else
      call inputchem_simple(inputchemfile,outdir,dirsep)
    endif
  endif

  ! Shallow cumulus ventilation allocation (2017)
  allocate(wscm2(nchsp))
  wscm2 = 0.
 

  open (20, file=trim(outdir)//dirsep//'output_dyn')
    write (20,'(a4)') 'TIME'
    write (20,'(I4)') time/atime
    write (20,'(19a14)') 'UTC(hours)','RT(hours)','zi(m)','we(m/s)', &
            'thetam(K)','dtheta(K)','wte(Km/s)','wts(Km/s)', &
            'beta','um(ms-1)','du(ms-1)','vm(ms-1)','dv(ms-1)','ws(ms-1)','ac(-)','wm(ms-1)','dz(m)','wstar', 'dthetav'

  open (28, file=trim(outdir)//dirsep//'t_prof')
    write (28,'(a4)') 'VERT'
    write (28,'(I3)') time/atime_vert
    write (28,'(I2)') n_vert
    write (28,'(3a14)') 'z (m)','theta (K)','wtheta'

  open (29, file=trim(outdir)//dirsep//'q_prof')
    write (29,'(a4)') 'VERT'
    write (29,'(I3)') time/atime_vert
    write (29,'(I2)') n_vert
    write (29,'(3a14)') 'z (m)','q (g/kg)','wq'

  open (30, file=trim(outdir)//dirsep//'c_prof')
    write (30,'(a4)') 'VERT'
    write (30,'(I3)') time/atime_vert
    write (30,'(I2)') n_vert
    write (30,'(3a14)') 'z (m)','c (ppm)','wc'

  open (31, file=trim(outdir)//dirsep//'u_prof')
    write (31,'(a4)') 'VERT'
    write (31,'(I3)') time/atime_vert
    write (31,'(I2)') n_vert
    write (31,'(3a14)') 'z (m)','u (m/s)','wu'

  open (32, file=trim(outdir)//dirsep//'v_prof')
    write (32,'(a4)') 'VERT'
    write (32,'(I3)') time/atime_vert
    write (32,'(I2)') n_vert
    write (32,'(3a14)') 'z (m)','v (m/s)','wv'
    
 
  open (50, file=trim(outdir)//dirsep//'beta_shear')
     write (50,*) 'UTC(hours) ustar  uws  vws  uwe  vwe  du   dv   dVe '

  open (60, file=trim(outdir)//dirsep//'output_sca')
    write (60,'(a4)') 'TIMS'
    write (60,'(I4)') time/atime
    write (60,'(15a14)') &
            'UTC(hours)','RT(h-ours)','zi(m)','qm(g/kg)','dq(g/kg)','wqe', &
            'wqs','   betaq','  cm(ppm)','dc(ppm)', 'wce','wcs','  betac','q2','wqm'
   
 
  if(llandsurface)then
  open (61, file=trim(outdir)//dirsep//'output_land')
    write (61,'(a4)') 'TIME'
    write (61,'(I4)') time/atime
    write (61,'(12a14)') 'UTC(hours)','RT(hours)','SWin(W/m2)','SWout(w/m2)', &
            'Lwin(W/m2)','LWout(W/m2)','Qtot(W/m2)','SH(W/m2)','LE(W/m2)', &
            'GR(W/m2)','ra(s/m)','rs(s/m)' 
    if(lCO2Ags)then
    open (63, file=trim(outdir)//dirsep//'output_ags')
      write (63,'(a4)') 'TIME'
      write (63,'(I4)') time/atime
      write (63,'(12a14)') 'UTC(hours)','RT(hours)','An(mgC/m2s)','Resp(mgC/m2s)'
    endif
  endif

  if(lchem)then
    write(formatstring,'(A,i2,A)')'(',nchsp+2,'A15)'    !nchsp + 2 places for the 2 times
    open (40, file=trim(outdir)//dirsep//'chem_conc')
      write (40, '(a4)') 'CHEM'
      write (40,'(I4)') time/atime
      write (40,formatstring) 'UTC(hours)','RT(hours)',(PL_scheme(k)%name,k=1,nchsp)

    open (41, file=trim(outdir)//dirsep//'chem_entr')
      write (41, '(a4)') 'CHEM'
      write (41,'(I4)') time/atime
      write (41,formatstring) 'UTC(hours)','RT(hours)',(PL_scheme(k)%name,k=1,nchsp)

    open (42, file=trim(outdir)//dirsep//'chem_beta')
       write (42, '(a4)') 'CHEM'
       write (42,'(I4)') time/atime
       write (42,formatstring) 'UTC(hours)','RT(hours)',(PL_scheme(k)%name,k=1,nchsp)

    if(lvbs)then
      open (43, file=trim(outdir)//dirsep//'soa_part')
       write (43, '(a4)') 'SOA'
       write (43,'(I4)') time/atime
       write (43,formatstring)  'UTC(hours)','RT(hours)','OAbgb','Coa','C1mc','C2mc','C3mc','C4mc','X1','X2','X3','X4', &
          'beta_terp', 'beta_iso'

      open (44, file=trim(outdir)//dirsep//'soa_ftr')
       write (44, '(a4)') 'SOA'
       write (44,'(I4)') time/atime
       write (44,formatstring) 'UTC(hours)','RT(hours)','OAbgb','Coa','C1mc','C2mc','C3mc','C4Tmc','C4mc','X1','X2','X3','X4'
    endif
    
    if(lBVOC)then
      open (45, file=trim(outdir)//dirsep//'voc_em')
         write (45, '(a4)') 'VOC'
         write (45,'(I4)') time/atime
         write (45,formatstring) 'UTC(hours)','RT(hours)','F_ISO', 'F_TERP','gamma_ISO', 'gammaT_ISO', 'gammap', 'gammace_ISO', &
         'gamma_MT', 'gammaT_MT', 'gammace_MT', 'gammaLAI', 'gammasm', 'PAR', 'Ts' 
    endif

    open (46, file=trim(outdir)//dirsep//'initial_chem')
       write (46,*) '# name initial_c_CBL c_ft0    emission function'
       do k=1,nchsp
          write (46,'(2(A4),3(1X,F9.4),i5)') '#   ',PL_scheme(k)%name, c_cbl(k), c_ft(k), Q_cbl(k),Q_func(k)
       enddo

    open (47, file=trim(outdir)//dirsep//'chem_photo')
       write (47, '(a4)') 'CHEM'
       write (47,'(I4)') time/atime
       write (47,'(9a14)') 'UTC(hours)','RT(hours)','jo3*1000*60','jno2*60(min)','jch2o*100*60(min)',&
       'photo','zenith','RH_emiss', 'angle'

    open (48, file=trim(outdir)//dirsep//'chem_ftr')
       write (48, '(a4)') 'CHEM'
       write (48,'(I4)') time/atime
       write (48,formatstring) 'UTC(hours)','RT(hours)',(PL_scheme(k)%name,k=1,nchsp)

    open (49, file=trim(outdir)//dirsep//'SVOC')
      write (49, '(a4)') 'SVOC'
      write (49, '(I4)') time/atime
        write (49,formatstring) 'UTC(hours)','RT(hours)', 'alp_OH1', 'alp_OH2', 'alp_OH3', 'alp_OH4', 'alp_O31', 'alp_O32', 'alp_O33', 'alp_O34', 'alp_NO31', 'alp_NO32', 'alp_NO33', &
         'alp_NO34', 'bep_OH1', 'bep_OH2', 'bep_OH3', 'bep_OH4', 'bep_O31', 'bep_O32', 'bep_O33', 'bep_O34', 'bep_NO31', 'bep_NO32', 'bep_NO33', 'bep_NO34', 'lim_OH1', 'lim_OH2', 'lim_OH3', &
         'lim_OH4', 'lim_O31', 'lim_O32', 'lim_O33', 'lim_O34', 'lim_NO31', 'lim_NO32', 'lim_NO33', 'lim_NO34', 'isoo32', 'isoo33', 'isoo34', 'isono32', 'isono33', 'isono34'

    open (53, file=trim(outdir)//dirsep//'SVOC_ft')
      write (53, '(a4)') 'SVOC_ft'
      write (53, '(I4)') time/atime
        write (53,formatstring) 'UTC(hours)','RT(hours)', '1', '2', '1', '2', '1', '2', '1', '2'

    open (52, file=trim(outdir)//dirsep//'uptake_IEPOX')
      write (52,'(a4)') 'TIMS'
      write (52, '(I4)') time/atime
      write (52, '(16a14)') 'UTC(hours)','RT(hours)', 'uptakeIEPOX', 'q_aq', 'inv_gas', 'inv_Torq_aq', 'F_gIpx', 'k_org', 'omega', 'Rc_gIPX', 'T_hh', 'denom', 'q_org_new', 'HorgDorg', 'q_org', 'k_het'

    open ( 62, file =trim(outdir)//dirsep//'keff_cbl')
      write(formatstring,'(A,i3,A4)')'(',tnor+3,'A13)' !+3 id for 2 times and tem_cbl
      write (62,formatstring) 'UTC(hours)','RT(hours)','Temp_cbl',(RC(k)%rname,k=1,tnor)
 
  if(lscu)then
    write(formatstring2,'(A,i2,A)')'(',nchsp+2,'A15)'    !nchsp + 2 places for the 2 times
   open(70, file=trim(outdir)//dirsep//'shallow_cumulus')
    write (70, '(a4)') 'TIME'
    write (70, '(I4)') time/atime
    write (70, formatstring2) 'UTC(hours)','RT(h-ours)',(PL_scheme(k)%name,k=1,nchsp)
   endif

    if (lwritepl) then
      do i=1,nchsp
        if((.not. lcomplex).or.(PL_scheme(i)%prin)) then
          filename = PL_scheme(i)%name
          filename = trim(outdir)//dirsep//'PL'//dirsep//trim(PL_scheme(i)%name)
          do j=1,PL_scheme(i)%nr_PL
            if (PL_scheme(i)%PL(j)%PorL == PRODUCTION ) then
              dummy(j)='P_'
            else
              dummy(j)='L_'
            endif
          enddo
          open(100+i,FILE=trim(filename))
          write(formatstring,'(A,i3,A4)')'(',PL_scheme(i)%nr_PL+4,'A17)'
          write(100+i,formatstring) 'RT(hours)','['//trim(PL_scheme(i)%name)//']',(dummy(j)//RC(PL_scheme(i)%PL(j)%r_nr)%rname,j=1,PL_scheme(i)%nr_PL),'tot_loss','tot_prod'
        endif  
      enddo
    endif !lwritepl  
  endif !lchem

! initialisation

! checking init
  write (*,*) ' checking namrun t=0'
  write (*,*) 'time',time
  write (*,*) 'dtime',dtime
  write (*,*) 'the following only used for the complex'

  write (*,*) ' checking the initialisation dyn t=0'
  write (*,*) 'entrainment_ratio_(beta)= ',beta
  write (*,*) 'max surface_heat_flux_(wthetas)= ',wthetasmax
  write (*,*) 'exchange_temperature_coeff._in_FT_(gamma)= ',gamma
  write (*,*) 'initial_boundary_layer_height_(zi0)= ',zi0
  write (*,*) 'mixing-layer_pot._temp._(thetam0)= ',thetam0
  write (*,*) 'initial_pot._temp._jump_(dtheta0)= ',dtheta0
  write (*,*) 'timestep_(dtime)=',dtime
  write (*,*) 'uws=',uws
  write (*,*) 'vws=',vws
  write (*,*) 'ustar=',ustar
  write (*,*) 'initial_x.-dir._windspeed_(um0)= ',um0
  write (*,*) 'initial_y-dir._windspeed_(vm0)= ',vm0
  write (*,*) 'du=',du
  write (*,*) 'dv=',dv
  write (*,*) 'uwe=',uwe
  write (*,*) 'vwe=',vwe
  write (*,*) ' '
  write (*,*) ' checking the initialisation c_cbl c_ft0 emis at t=0'

! write (*,*) specname(i_R), c_cbl(i_R), c_ft0(i_R ),Q(i_R )

! dynamics

  tt=0
! run
  do t=1, runtime
    tt=tt+1
    sec = t * dtime    !number of seconds from start
    printhour=t * dtime/3600.
    thour= hour+t*dtime/3600.

    if ( thour > 24. ) then
      thour = thour - 24 * (ceiling(thour/24)-1)
    endif

    zsl = 0.1*zi(1) !Height of surface layer

!   Start with radiation calculation
    if (lradiation) then
      costh = max(0.0,cos(getth(1.0*day,latt,long,thour))) 
      Ta    = thetam(1)*((((100*pressure)-zsl*rho*g)/(100*pressure))**(Rd/Cp))!pressure*100 to compensate for SI, 0.1 to get T at top of the SL
      Tr    = (0.6 + 0.2 * costh) * (1 - 0.4 * cc)

      Swin  = S0 * Tr * costh + DeltaFsw
      Swout = albedo * Swin
      Lwin  = 0.8 * bolz * (Ta ** 4)
      Lwout = bolz * (Ts ** 4)
     ! WRITE(*,*) bolz, Ts
      Qtot  = Swin - Swout + Lwin - Lwout
    endif

    thetav  = thetam(1) * (1. + 0.61 * qm(1) * 1.e-3)
    wstar   = ( (g / thetav) * zi(1) * wthetav ) ** ( 1. / 3. )
    ueff    = sqrt(um(1)**2.+vm(1)**2.+wstar**2.)
!    ueff    = sqrt(um(1)**2.+vm(1)**2.)
    if (lsurfacelayer) then
      if(ueff .lt. 1.e-2) stop 'effective wind velocity below 1 cm/s' !changed from 1e-2
      if(t==1) then
        thetasurf = Ts 
      else
        thetasurf = thetam(1) + wthetas / (Cs * ueff)
      endif
      esatsurf    = 0.611e3 * exp(17.2694 * (thetasurf - 273.16) / (thetasurf - 35.86))
      qsatsurf    = 0.622 * esatsurf / (pressure*100)
      cq          = 0.0
      if(t/=1) cq = (1. + Cs * ueff * rs) ** (-1.)
      if(lsea) cq = 1.0
      qsurf       = (1. - cq) * qm(1) + cq * qsatsurf * 1.e3 !HGO factor for qsatsurf which is in kg/kg

      thetavsurf  = thetasurf * (1. + 0.61 * qsurf * 1.e-3)

      Rib         = min(0.2,g/thetav * zsl * (thetav-thetavsurf) / (ueff ** 2.))
      L           = sign(dble(0.01),Rib)
      L0          = sign(dble(0.1),Rib)

      iter        = 0
      do while(.true.)
        iter      = iter + 1
        L0        = L
        fx        = Rib - zsl / L * (log(zsl / z0h) - psih(zsl / L) + psih(z0h / L)) / (log(zsl / z0m) - psim(zsl / L) + psim(z0m / L)) ** 2.
        Lstart    = L - 0.001*L
        Lend      = L + 0.001*L
        fxdif     = ( (- zsl / Lstart * (log(zsl / z0h) - psih(zsl / Lstart) + psih(z0h / Lstart)) / (log(zsl / z0m) - psim(zsl / Lstart) + psim(z0m / Lstart)) ** 2.) &
                  - (-zsl / Lend * (log(zsl / z0h)- psih(zsl / Lend) + psih(z0h / Lend)) / (log(zsl / z0m) - psim(zsl / Lend) + psim(z0m / Lend)) ** 2.) ) / (Lstart - Lend)
        L         = L - fx / fxdif
        L         = sign(min(abs(L),1.e6),L)!capping L

        if(abs((L - L0)/L) < 1e-4) exit 
        if(abs((L - L0)) < 1e-3) exit 
      enddo

      Constm =  kappa ** 2. / (log(zsl / z0m) - psim(zsl / L) + psim(z0m / L)) ** 2.
      Cs     =  kappa ** 2. / (log(zsl / z0m) - psim(zsl / L) + psim(z0m / L)) / (log(zsl / z0h) - psih(zsl / L) + psih(z0h / L))

      ustar  = sqrt(Constm) * ueff
      if(ustar .le. 0) stop "ustar has to be greater than 0"
      uws    = - Constm * ueff * um(1)
      vws    = - Constm * ueff * vm(1)
      
      T2m    = thetasurf - wthetas / ustar / kappa * (log(2. / z0h) - psih(2. / L) + psih(z0h / L))
      q2m    = qsurf     - wqs     / ustar / kappa * (log(2. / z0h) - psih(2. / L) + psih(z0h / L))
      u2m    =           - uws     / ustar / kappa * (log(2. / z0m) - psim(2. / L) + psim(z0m / L))
      v2m    =           - vws     / ustar / kappa * (log(2. / z0m) - psim(2. / L) + psim(z0m / L))
      esat2m = 0.611e3 * exp(17.2694 * (T2m - 273.16) / (T2m - 35.86))
      e2m    = q2m * 1.e-3 * (100*pressure) / 0.622  !HGO factor for q2m which is in g/kg
      rh2m   = e2m / esat2m

    else !lsurfacelayer
      ! Two options are considered regarding the momentum surface fluxes.
      ! Introduce time-dependent uws ,vws and ustar (z0 constant) or constant momentum fluxes.
      
      ! --------------------------------------------------------------

      if (c_ustr) then
        if (t == 1) then
             write(*,*) 'Constant friction velocity'
        end if
        uws=uws0
        vws=vws0
        ustar = ((uws**2.)+(vws**2.))**0.25
      else
        ustar= kappa*sqrt(um(1)**2.+vm(1)**2.)/log(zp/z0)

        if (um(1).ne.0.) then
          alpha=atan(vm(1)/um(1))
        else
          if (vm(1).ne.0.) then
              alpha=pi/2.
          endif
        endif

        uws=-ustar**2.*cos(alpha)
        vws=-ustar**2.*sin(alpha)

      endif

    endif !lsurfacelayer

!   Introduce time-dependent momentum, humidity and heat surface fluxes.
!   We assume a sinusoidal behaviour of all three variables.

    if (c_wth .or. c_fluxes) then
      if (llandsurface .and. (t == 1)) then
        print *,""
        print *,""
        print *,"!!!!!!!!!!!!!!!!!!!!!"
        print *,""
        print *,"You enabled constant fluxes (c_wth/c_fluxes) AND the interactive"
        print *,"land surface scheme (llandsurface)!"
        print *,"The constant fluxes override the land surface model."
        print *,""
        print *,"!!!!!!!!!!!!!!!!!!!!!"
        print *,""
        print *,""
      endif
      if (t == 1) then
        write(*,*) 'Constant surface flux'
      end if
      wthetas = wthetasmax
      wqs     = wqsmax
      wcs     = wcsmax
      wthetav = wthetasmax+rvrd*thetam(1)*wqsmax*1.0e-03
      do i=1,nchsp
        Q_cbl(i) = Q_init(i)
       
      enddo
    else ! c_wth
      if (llandsurface) then
        if(ustar .le. 0) stop "ustar has to be greater than 0"
        ra       = ueff / (ustar ** 2.)

        if (lsea) then ! Sea surface
          esatsurf     = 0.611e3 * exp(17.2694 * (sst - 273.16) / (sst - 35.86))
          qsatsurf     = 0.622 * esatsurf / (pressure*100)

          LE           = rho * Lv / ra * (qsatsurf - qm(1) * 1.0e-03)
          SH           = rho * Cp / ra * (sst      - thetam(1))
        else           ! Land surface
          esat     = 0.611e3 * exp(17.2694 * (thetam(1) - 273.16) / (thetam(1) - 35.86))
          qsat     = 0.622 * esat / (pressure*100)
          desatdT  = esat * (17.2694 / (thetam(1) - 35.86) - 17.2694 * (thetam(1) - 273.16) / (thetam(1) - 35.86)**2.)
          dqsatdT  = 0.622 * desatdT / (pressure*100)
          efinal   = qm(1) * 1.0e-3 * (pressure*100) / 0.622 

          if (lradiation) then
                  f1 = 1.0 / ((0.004 * Swin + 0.05) / (0.81 * (0.004 * Swin + 1.)))
          else
                  f1 = 1.0
          endif

          f2     = 1.e8
          if (w2 .gt. wwilt) then
            f2   = (wfc - wwilt) / (w2 - wwilt)
          endif
          f2     = max(1.0,min(1.e8, f2))

          if (lsurfacelayer) then
            f3   = 1.0 / exp( - gD * (esat2m - e2m)/100.0 )
            f4   = 1.0 / (1.0 - 0.0016 * (298.0 - T2m) ** 2.)
          else !Then just use BL-averaged values instead of 2 m values
            f3   = 1.0 / exp( - gD * (esat - efinal)/100.0)
            f4   = 1.0 / (1.0 - 0.0016 * (298.0 - thetam(1)) ** 2.)
          endif

          rs     = (rsmin / LAI) * f1 * f2 * f3 * f4 

          if (lrsAgs) then
            if(lchem .and. (CO2%loc .gt. 0) ) then
              CO2ags = c_cbl(CO2%loc) / 1000.0 !CO2ags in ppm
            else
              CO2ags = cm(1) / 1000.0 !CO2ags in ppm
            endif

            ! calculate surface resistances using plant physiological (A-gs) model
            ! calculate CO2 compensation concentration
            if (thetasurf > 360.) then
                    write(*,*) thetasurf
                    thetasurf = 320.
            endif
            CO2comp  = CO2comp298 * Q10CO2 ** ( 0.1 * (thetasurf - 298.0) )
            CO2comp  = CO2comp * rho

            ! calculate mesophyll conductance
            gm       = gm298 * Q10gm ** (0.1 * (thetasurf - 298.0) ) / ( (1. + exp(0.3 * (T1gm - thetasurf))) * (1. + exp(0.3 * (thetasurf - T2gm))))
            gm       = gm / 1000   ! conversion from mm s-1 to m s-1

            ! calculate CO2 concentration inside the leaf (ci)
            fmin0    = gmin/nuco2q - (1.0/9.0) * gm
            fmin     = (-fmin0 + ( fmin0 ** 2.0 + 4 * gmin/nuco2q * gm ) ** (0.5)) / (2. * gm)

            esatsurf = 0.611e3 * exp(17.2694 * (Ts - 273.16) / (Ts - 35.86))
            Ds       = (esatsurf - efinal) / 1000.0 ! in kPa
            D0       = (f0 - fmin) / ad

            cfrac    = f0 * (1.0 - Ds/D0) + fmin * (Ds/D0)
            co2abs   = CO2ags*(MW_CO2/MW_Air)*rho
            ci       = cfrac * (co2abs - CO2comp) + CO2comp

            ! calculate maximal gross primary production in high light conditions (Ag)
            Ammax    = Ammax298 * Q10Am ** ( 0.1 * (thetasurf - 298.0) ) / ( (1. + &
                       exp(0.3 * (T1Am - thetasurf))) * (1. + exp(0.3 * (thetasurf - T2Am))))

            ! calculate effect of soil moisture stress on gross assimilation rate
            betaw    = max(1.0e-3,min(1.0,(w2 - wwilt)/(wfc - wwilt)))

            ! calculate stress function
            fstr     = betaw

            ! calculate gross assimilation rate (Am)
            Am       = Ammax * (1 - exp( -(gm * (ci - CO2comp) / Ammax) ) )

            Rdark    = (1.0/9) * Am

            PAR      = 0.50*max(0.1,Swin) 

            ! calculate  light use efficiency
            alphac   = alpha0 * (co2abs - CO2comp) / (co2abs + 2 * CO2comp)

            ! 1.- calculate upscaling from leaf to canopy: net flow CO2 into the plant (An) 
            tempy    = alphac * Kx * PAR / (Am + Rdark)
            An       = (Am + Rdark) * (1 - 1.0 / (Kx * LAI) * (E1( tempy * exp(-Kx * LAI)) - E1(tempy)))


            ! 2.- calculate upscaling from leaf to canopy: CO2 conductance at canopy level
            AGSa1    = 1.0 / (1 - f0)
            Dstar    = D0 / (AGSa1 * (f0-fmin))

            gcco2    = LAI * (gmin/nuco2q + AGSa1 * fstr * An / ((co2abs - CO2comp) * (1 + Ds / Dstar)))

            ! calculate surface resistance for moisture and carbon dioxide
            rsAgs    = 1.0 / (1.6 * gcco2)
            rsCO2    = 1.0 / gcco2

            rs       = rsAgs

            ! calculate net flux of CO2 into the plant (An)
            An       = -(co2abs - ci) / (ra + rsCO2)

            ! CO2 soil respiration surface flux
            fw       = Cw * wsmax / (w2 + wsmin)
!            Resp     = R10 * (1 - fw)*(1 + ustar) * exp( Eact0 / (283.15 * 8.314) * (1. - 283.15 / (Tsoil))) ! substitute thetasurf_>T soil
            Resp     = R10 * (1 - fw) * exp( Eact0 / (283.15 * 8.314) * (1. - 283.15 / (Tsoil))) ! substitute thetasurf_>T soil ustar is also removed

            wco2     = (An + Resp)*(MW_Air/MW_CO2)*(1.0/rho) !in ppm m/s

          endif

          ! recompute f2 using wg instead of w2
          f2     = 1.e8
          if (wg .gt. wwilt) then
            f2   = (wfc - wwilt) / (wg - wwilt)
          endif
          f2     = max(1.0,min(1.e8, f2))

          rssoil = rssoilmin * f2

          Wlmx   = LAI * Wmax
          cliq   = min(1.0, Wl / Wlmx)
          Ts     = (Qtot + rho * Cp / ra * thetam(1) + cveg * (1.0-cliq) * rho &
                 * Lv / (ra + rs) * (dqsatdT * thetam(1) - qsat + qm(1) * 1.0e-3) &
                 + (1.0 - cveg) * rho * Lv / (ra + rssoil) * (dqsatdT * thetam(1) - qsat + qm(1) * 1.0e-3) &
                 + cveg * cliq * rho * Lv / ra * (dqsatdT * thetam(1) - qsat + qm(1) * 1.0e-3) + Lambda * Tsoil) &
                 * (rho * Cp / ra + cveg * (1. - cliq) * rho * Lv / (ra + rs) * dqsatdT + (1. - cveg) &
                 * rho * Lv / (ra + rssoil) * dqsatdT + cveg * cliq * rho * Lv / ra * dqsatdT + Lambda) ** (-1.)
                 
          esatsurf = 0.611e3 * exp(17.2694 * (Ts - 273.16) / (Ts - 35.86))
          qsatsurf = 0.622 * esatsurf / (pressure*100)

          LEveg  = (1.0 - cliq) * cveg  * rho * Lv / (ra + rs)     * (dqsatdT * (Ts - thetam(1)) + qsat - qm(1) * 1.0e-3)
          LEliq  =        cliq  * cveg  * rho * Lv /  ra           * (dqsatdT * (Ts - thetam(1)) + qsat - qm(1) * 1.0e-3)
          LEsoil =         (1.0 - cveg) * rho * Lv / (ra + rssoil) * (dqsatdT * (Ts - thetam(1)) + qsat - qm(1) * 1.0e-3)

          Wltend = - LEliq / (rhow * Lv)
          Wl     = Wl + Wltend * dtime

          LE     = LEsoil + LEveg + LEliq
          SH     = rho * Cp / ra * (Ts - thetam(1))
          GR     = Lambda * (Ts - Tsoil)
          
          if (t>3200 .and. SH .le. 0) then
             print *,""
             print *,""
             print *,""
             print *,""
             print *,"Sensible heat flux is below 0"     
             print *,""
             print *,""
             print *,""
             stop    
         endif         

          Tsoiltend = CG * GR - 2.0 * pi / 86400.0 * (Tsoil - T2)
          Tsoil  = Tsoil + Tsoiltend * dtime

          C1     = C1sat * (wsat / wg) ** (CLb / 2.0 + 1.0)
          C2     = C2ref * (w2 / (wsat - w2) )
          wgeq   = w2 - wsat * CLa * ( (w2 / wsat) ** CLc * (1.0 - (w2 / wsat) ** (8.0 * CLc)) )
          wgtend = - C1 / (rhow * 0.1) * LEsoil / Lv - C2 / 86400 * (wg - wgeq)
          wg     = wg + wgtend * dtime
         ! WRITE(*,*)t, LE, SH, GR 
        endif ! lsea

        wthetas= SH / (rho * Cp)
        wqs    = LE / (rho * Lv * 1.0e-3)
        wcs    = wcsmax

        if (lBVOC) then  ! BVOC emissions
          if (llandsurface .eqv. .false. .or. lrsAgs .eqv. .false.) then
            print *,""
            print *,""
            print *,"!!!!!!!!!!!!!!!!!!!!!"
            print *,""
            print *,"You enabled interactive VOC emissions without enabling"
            print *,"the land surface scheme (llandsurface) and/or A-gs (lrsAgs)"
            print *,"VOC emissions can not be calculated due to missing input variables!"
            print *,""
            print *,"!!!!!!!!!!!!!!!!!!!!!"
            print *,""
            print *,""
            stop 'switch land surface & A-gs on'
          endif

           ! MT-emission following Guenther'06, using parameterized canopy environment emission activity (PCEEA) algorithm, and Sakulnayontvittaya'08      
           gammaLAI   = 0.49*LAI/((1.+0.2*LAI**2.)**0.5)
           
           w1         = wwilt + dwg  ! (eq. 20 Guenther06)
           if (wg .gt. w1) then  
           gammasm    = 1.  
           elseif (wg .lt. wwilt) then
           gammasm    = 0.
           else
           gammasm    = (wg-wwilt)/dwg
           endif

           beta_MT    = 0.1 ! 0.13 is from (K-1) (Sakulnayontvittaya 08)
           Tref       = 303  ! reference temperature for bVOC BaseEmissionRate (K) 
           gammaT_MT  = exp(beta_MT*(Ts-Tref))

           gammace_MT = gammaLAI * gammaT_MT ! emission activity factor canopy environment (-) (eq.10, Guenther96)
           gamma_MT   = gammace_MT * gammasm ! emission activity factor terpene (-) (eq.2, Guenther96) (no light and leaf age dependence)
             
           eps_MT     = BaserateTer ! emission factor representing net in-canopy emission rate at standard conditions (μg m-2 h-1 at 303 K)
           rho_MT     = 1. ! factor that accounts for chemical production and loss within plant canopies (-)
           ER_MT      = eps_MT * gamma_MT * rho_MT ! monoterpene emission rate (μg m-2 h-1)
            
           mm_MT      = 10.*12.01 + 16.*1.008    ! molecular mass C10H16
           cf_MT      = mm_MT*((pressure*100.)/(Rgas*Ts)*3600*1e-3) ! flux conversion factor ppb m s-1 -> ug m-2 h-1
           F_TERP     = ER_MT / cf_MT             ! isoprene emission rate (ppb m s-1)

           !ISO-emission following Guenther'06, using parameterized canopy environment emission activity (PCEEA) algorithm         
           Thr        = thetam(2)  ! changed to hourly average temperature (K) as per Alex Guenther's recommendation (Nagori 2019)
           Tdaily     = 298. ! daily average temperature (K)
           Topt       = 313 + (0.6*(Tdaily-297)) ! (adapted from eq. 8 Guenther06)
           C_T1       = 80.  ! empirical coefficient
           C_T2       = 200. ! empirical coefficient
           x          = ((1/Topt)-(1/Thr))/0.00831
           Eopt       = 1.75*exp(0.08*(Tdaily-297)) 
           gammaT_ISO = Eopt * C_T2 * exp(C_T1*x)/(C_T2-C_T1*(1-exp(C_T2*x)))

           aa         = -zenith*(180.0/pi) + 90. ! solar angle relative to surface = -zenith angle + 90 (deg.)
           if (t==1) then              
             aa = 26.57              
           endif
           !!!!!! Done since solar angle wrong in the first time step (2017 for Nagori 2019)


           Pac        = PAR * 4.766 !above canopy PPFD (umol m-2 s-1) 4.766 is the conversion from W/m2 to micromol/m2/s (Nagori 2019)
           Pdaily     = 500 !! daily average above canopy PPFD (umol m-2 s-1) as recommended by Alex Guenther (Nagori 2019)
           Ptoa       = 3000. + 99*cos((2*pi*((day-10.)/365.))*pi/180.0) ! 
           phi        = Pac/(sin(aa*(pi/180.0))*Ptoa)! above canopy PPFD transmission (-) ! 
           gammap     = sin(aa*(pi/180.0))*(2.46*(1.+0.0005*(Pdaily-400.))*phi*0.9*phi**2.) ! 

           gammace_ISO= gammaLAI * gammaT_ISO * (gammap)
           gamma_ISO  = (gammace_ISO) * gammasm
               
           eps_ISO    = BaserateIso  ! emission factor representing net in-canopy emission rate at standard conditions (μg m-2 h-1 at 303 K)
           rho_ISO    = 0.96     
           ER_ISO     = eps_ISO * gamma_ISO * rho_ISO ! isoprene emission rate (μg m-2 h-1)
           mm_ISO     = 5.*12.01 + 8.*1.008      ! molecular mass C5H8
           cf_ISO     = mm_ISO*((pressure*100.)/(Rgas*Ts)*3600*1e-3)! flux conversion factor ppb m s-1 -> ug m-2 h-1
           F_ISO      = ER_ISO / cf_ISO           ! monoterpene emission rate (ppb m s-1)
           
        endif ! bVOC

      else !llandsurface

        if (sec .le. daytime_start .or. sec > daytime_end) then
                wcs = 0.
        else
                wcs = wcsmax * sin(pi * (sec - daytime_start)/daylength)
        end if

        select case(function_wT)
          case(0)
            wthetas = 0.
          case(1)
            wthetas = wthetasmax
          case(2)
            if ((sec .le. starttime_wT) .or. (sec .ge. endtime_wT)) then
              wthetas = 0.
            else
              wthetas = wthetasmax * sin(pi * (sec - starttime_wT)/(endtime_wT - starttime_wT)) 
            end if
          case(3)
            if ((sec .le. starttime_wT) .or. (sec .ge. endtime_wT)) then
              wthetas = 0.
            else
              wthetas = wthetasmax 
            end if
          case(4)
            if ((sec .le. starttime_wT) .or. (sec .ge. endtime_wT)) then
              wthetas = 0.
            else
              wthetas = (wthetasmax/2) * (1 - cos(2*pi*(sec - starttime_wT)/(endtime_wT - starttime_wT)))
            end if
          case default
            if (t==1) print *,'Flag for the function of wT is invalid: ',function_wT
            stop 'change function'
        end select
        wthetas = wthetas + offset_wT

        select case(function_wq)
          case(0)
            wqs = 0.
          case(1)
            wqs = wqsmax
          case(2)
            if ((sec .le. starttime_wq) .or. (sec .ge. endtime_wq)) then
              wqs = 0.
            else
              wqs = wqsmax * sin(pi * (sec - starttime_wq)/(endtime_wq - starttime_wq)) 
            end if
          case(3)
            if ((sec .le. starttime_wq) .or. (sec .ge. endtime_wq)) then
              wqs = 0.
            else
              wqs = wqsmax 
            end if
          case(4)
            if ((sec .le. starttime_wq) .or. (sec .ge. endtime_wq)) then
              wqs = 0.
            else
              wqs = (wqsmax/2) * (1 - cos(2*pi*(sec - starttime_wq)/(endtime_wq - starttime_wq)))
            end if
          case default
            if (t==1) print *,'Flag for the function of wq is invalid: ',function_wq
            stop 'change function'
        end select
        wqs = wqs + offset_wq

      endif !llandsurface

      wthetav = wthetas + rvrd*thetam(1)*wqs*1.0e-03

      !G/p partitioning in the CBL moved up to separate the deposition for SVOC and aerosol (2017 for Nagori et al., 2019)

      if(lvbs) then
        ! G/P-partitioning in the BL introduced earlier to separate for deposition 
        T_hh    = thetam(2)+gammad*(zi(2)/2.)         ! T at h/2; dry adiabatic lapse rate T = (theta-g/Cp)*z
        T_ft = (thetam(2)+dtheta(2))+gammad*zi(2)  ! temperature at the inversion: <theta> + dtheta
    
        mr2mc_PRODT = ((pressure*100.)*PRODTmm)/(Rgas*T_hh) * 1e-3 ! factor to convert TERP products (SVOCs) from ppb to ug/m3
        mr2mc_PRODI = ((pressure*100.)*PRODImm)/(Rgas*T_hh) * 1e-3 ! factor to convert ISO products (SVOCs) from ppb to ug/m3
        mr2mc_OAbg  = ((pressure*100.)*mmOAbg)/(Rgas*T_hh) * 1e-3  ! factor to convert OAbg from ppb to ug/m3 
        !print *, mr2mc_OAbg
        ! Done per oxidant + BVOC (a-pinnene, b-pinene and limonene for easier tracking (Nagori et al., 2019)
        OAbgb = c_cbl(OAbg%loc)  * mr2mc_OAbg  ! convert mixing ratio (ppb) to mass concentration (ug m-3)
        APOHmc = c_cbl(APOH%loc) * mr2mc_PRODT  
        APO3mc = c_cbl(APO3%loc) * mr2mc_PRODT  
        APNO3mc = c_cbl(APNO3%loc) * mr2mc_PRODT  
        BPOHmc = c_cbl(BPOH%loc) * mr2mc_PRODT  
        BPO3mc = c_cbl(BPO3%loc) * mr2mc_PRODT  
        BPNO3mc = c_cbl(BPNO3%loc) * mr2mc_PRODT  
        LIOHmc = c_cbl(LIOH%loc) * mr2mc_PRODT  
        LIO3mc = c_cbl(LIO3%loc) * mr2mc_PRODT  
        LINO3mc = c_cbl(LINO3%loc) * mr2mc_PRODT
        ISO3mc = c_cbl(ISO3%loc) * mr2mc_PRODI
        ISNO3mc = c_cbl(ISNO3%loc) * mr2mc_PRODI
                
        C1AHlmc  = alpha1_APOHlow  * APOHmc ! low NOx yield per ox per terp A=apin H=OH
        C2AHlmc  = alpha2_APOHlow  * APOHmc
        C3AHlmc  = alpha3_APOHlow  * APOHmc
        C4AHlmc  = alpha4_APOHlow  * APOHmc
        C1BHlmc  = alpha1_BPOHlow  * BPOHmc ! low NOx yield per ox per terp B=bpin H=OH
        C2BHlmc  = alpha2_BPOHlow  * BPOHmc
        C3BHlmc  = alpha3_BPOHlow  * BPOHmc
        C4BHlmc  = alpha4_BPOHlow  * BPOHmc
        C1AOlmc  = alpha1_APO3low  * APO3mc ! low NOx yield per ox per terp A=apin O=O3
        C2AOlmc  = alpha2_APO3low  * APO3mc
        C3AOlmc  = alpha3_APO3low  * APO3mc
        C4AOlmc  = alpha4_APO3low  * APO3mc
        C1BOlmc  = alpha1_BPO3low  * BPO3mc ! low NOx yield per ox per terp B=bpin O=O3
        C2BOlmc  = alpha2_BPO3low  * BPO3mc
        C3BOlmc  = alpha3_BPO3low  * BPO3mc
        C4BOlmc  = alpha4_BPO3low  * BPO3mc
        C1ANmc  = alpha1_APNO3  * APNO3mc ! only one for NO3 ox; A=apin N=NO3
        C2ANmc  = alpha2_APNO3  * APNO3mc
        C3ANmc  = alpha3_APNO3  * APNO3mc
        C4ANmc  = alpha4_APNO3  * APNO3mc
        C1BNmc  = alpha1_BPNO3  * BPNO3mc ! only one for NO3 ox; B=bpin N=NO3
        C2BNmc  = alpha2_BPNO3  * BPNO3mc
        C3BNmc  = alpha3_BPNO3  * BPNO3mc
        C4BNmc  = alpha4_BPNO3  * BPNO3mc

        C1AHhmc  = alpha1_APOHhigh  * APOHmc ! high NOx yield per ox per terp A=apin H=OH
        C2AHhmc  = alpha2_APOHhigh  * APOHmc
        C3AHhmc  = alpha3_APOHhigh  * APOHmc
        C4AHhmc  = alpha4_APOHhigh  * APOHmc
        C1BHhmc  = alpha1_BPOHhigh  * BPOHmc ! high NOx yield per ox per terp B=bpin H=OH
        C2BHhmc  = alpha2_BPOHhigh  * BPOHmc
        C3BHhmc  = alpha3_BPOHhigh  * BPOHmc
        C4BHhmc  = alpha4_BPOHhigh  * BPOHmc
        C1AOhmc  = alpha1_APO3high  * APO3mc ! high NOx yield per ox per terp A=apin O=O3
        C2AOhmc  = alpha2_APO3high  * APO3mc
        C3AOhmc  = alpha3_APO3high  * APO3mc
        C4AOhmc  = alpha4_APO3high  * APO3mc
        C1BOhmc  = alpha1_BPO3high  * BPO3mc ! high NOx yield per ox per terp B=bpin O=O3
        C2BOhmc  = alpha2_BPO3high  * BPO3mc
        C3BOhmc  = alpha3_BPO3high  * BPO3mc
        C4BOhmc  = alpha4_BPO3high  * BPO3mc

        C1LHlmc  = alpha1_LIOHlow  * LIOHmc ! low NOx yield per ox per terp L=limo H=OH
        C2LHlmc  = alpha2_LIOHlow  * LIOHmc
        C3LHlmc  = alpha3_LIOHlow  * LIOHmc
        C4LHlmc  = alpha4_LIOHlow  * LIOHmc
        C2IOlmc  = alpha2_ISO3low  * ISO3mc ! low NOx yield per ox for isoprene; O=O3        
        C3IOlmc  = alpha3_ISO3low  * ISO3mc
        C4IOlmc  = alpha4_ISO3low  * ISO3mc
        C1LOlmc  = alpha1_LIO3low  * LIO3mc ! low NOx yield per ox per terp L=limo O=O3
        C2LOlmc  = alpha2_LIO3low  * LIO3mc
        C3LOlmc  = alpha3_LIO3low  * LIO3mc
        C4LOlmc  = alpha4_LIO3low  * LIO3mc
        C2INmc  = alpha2_ISNO3  * ISNO3mc ! N=NO3        
        C3INmc  = alpha3_ISNO3  * ISNO3mc
        C4INmc  = alpha4_ISNO3  * ISNO3mc
        C1LNmc  = alpha1_LINO3  * LINO3mc ! L=limo N=NO3
        C2LNmc  = alpha2_LINO3  * LINO3mc
        C3LNmc  = alpha3_LINO3  * LINO3mc
        C4LNmc  = alpha4_LINO3  * LINO3mc

        C1LHhmc  = alpha1_LIOHhigh  * LIOHmc ! high NOx yield per ox per terp L=limo H=OH
        C2LHhmc  = alpha2_LIOHhigh  * LIOHmc
        C3LHhmc  = alpha3_LIOHhigh  * LIOHmc
        C4LHhmc  = alpha4_LIOHhigh  * LIOHmc
        C2IOhmc  = alpha2_ISO3high  * ISO3mc ! high NOx yield per ox for isoprene; O=O3        
        C3IOhmc  = alpha3_ISO3high  * ISO3mc
        C4IOhmc  = alpha4_ISO3high  * ISO3mc
        C1LOhmc  = alpha1_LIO3high  * LIO3mc ! high NOx yield per ox per terp L=limo O=O3
        C2LOhmc  = alpha2_LIO3high  * LIO3mc
        C3LOhmc  = alpha3_LIO3high  * LIO3mc
        C4LOhmc  = alpha4_LIO3high  * LIO3mc

        
      endif


      do i=1,nchsp
        select case (Q_func(i))
          case (0)
            Q_cbl(i) = 0.
          case (1)
            Q_cbl(i) = Q_init(i)
          case (2)
            if ((sec .le. starttime_chem) .or. (sec .ge. endtime_chem)) then
              Q_cbl(i) = 0.
            else
              Q_cbl(i) = Q_init(i) * sin(pi * (sec - starttime_chem)/(endtime_chem - starttime_chem))
            end if
          case (3)
            if ((sec .le. starttime_chem) .or. (sec .ge. endtime_chem)) then
              Q_cbl(i) = 0.
            else
              Q_cbl(i) = Q_init(i)
            endif
          case (4)
            if ((sec .le. starttime_chem) .or. (sec .ge. endtime_chem)) then
              Q_cbl(i) = 0. 
            else
              Q_cbl(i) = (Q_init(i)/2) * (1 - cos(2*pi*(sec - starttime_chem)/(endtime_chem - starttime_chem)))
            endif
          case (5)
              Q_cbl(i) = - Q_init(i) * c_cbl(i)
          case (6)
           if ((sec .le. starttime_chem) .or. (sec .ge. endtime_chem)) then
              Q_cbl(i) = 0. 
           else
              Q_cbl(i) = -(Q_init(i) * (sin(pi*(sec - starttime_chem)/(endtime_chem - starttime_chem)))) * c_cbl(i)
           endif
          !case (7) 
          !  if (lBVOC) then
          !    if(lchem .and. (TERP%loc .gt. 0) ) then
          !    Q_cbl(TERP%loc) = F_TERP
          !    !write (*,*) Q_cbl(TERP%loc)
          !   endif
          !    if(lchem .and. (ISO%loc .gt. 0) ) then
          !    Q_cbl(ISO%loc) = F_ISO
          !    !endif
          !    !write (*,*)  Q_cbl(ISO%loc)
          !    endif
          !  else
          !    print *, "!!!!!!!!"
          !    print *,"MEGAN emissions used even though lBVOC is off"
          !    print *, "!!!!!!!!"
          !    stop 'Switch lBVOC on'
          !  endif
          case (8) !SVOC deposition in gas-phase only (done to not include deposition once in aerosol phase); done per SVOC (Nagori
                  ! et al, 2019)
            if(lvbs)then
              Q_cbl(APOH%loc) = -Q_init(APOH%loc)* ( (1-X1)*((1-beta_terp)*C1AHlmc + beta_terp*C1AHhmc)/mr2mc_PRODT + (1-X2)*((1-beta_terp)*C2AHlmc + beta_terp*C2AHhmc)/mr2mc_PRODT +  (1-X3)*((1-beta_terp)*C3AHlmc + beta_terp*C3AHhmc)/mr2mc_PRODT + (1-X4)*((1-beta_terp)*C4AHlmc + beta_terp*C4AHhmc)/mr2mc_PRODT)
              Q_cbl(APO3%loc) = -Q_init(APO3%loc)* ( (1-X1)*((1-beta_terp)*C1AOlmc + beta_terp*C1AOhmc )/mr2mc_PRODT + (1-X2)*((1-beta_terp)*C2AOlmc +  beta_terp*C2AOhmc )/mr2mc_PRODT  + (1-X3)*((1-beta_terp)*C3AOlmc +  &
                    beta_terp*C3AOhmc )/mr2mc_PRODT + (1-X4)*((1-beta_terp)*C4AOlmc +  beta_terp*C4AOhmc )/mr2mc_PRODT)
              Q_cbl(APNO3%loc) = -Q_init(APNO3%loc)* ( (1-X1)*C1ANmc/mr2mc_PRODT + (1-X2)*C2ANmc/mr2mc_PRODT + (1-X3)*C3ANmc/mr2mc_PRODT + (1-X4)*C4ANmc/mr2mc_PRODT)
              Q_cbl(BPOH%loc) = -Q_init(BPOH%loc)* ( (1-X1)*((1-beta_terp)*C1BHlmc + beta_terp*C1BHhmc )/mr2mc_PRODT + (1-X2)*((1-beta_terp)*C2BHlmc +  beta_terp*C2BHhmc )/mr2mc_PRODT  + (1-X3)*((1-beta_terp)*C3BHlmc +  &
                beta_terp*C3BHhmc )/mr2mc_PRODT + (1-X4)*((1-beta_terp)*C4BHlmc +  beta_terp*C4BHhmc )/mr2mc_PRODT)
              Q_cbl(BPO3%loc) = -Q_init(BPO3%loc)* ( (1-X1)*((1-beta_terp)*C1BOlmc + beta_terp*C1BOhmc )/mr2mc_PRODT + (1-X2)*((1-beta_terp)*C2BOlmc + beta_terp*C4BOhmc )/mr2mc_PRODT +   (1-X3)*((1-beta_terp)*C3BOlmc + beta_terp*C3BOhmc )/mr2mc_PRODT + (1-X4)*((1-beta_terp)*C4BOlmc + beta_terp*C4BOhmc )/mr2mc_PRODT)
              Q_cbl(BPNO3%loc) = -Q_init(BPNO3%loc)* ( (1-X1)*C1BNmc/mr2mc_PRODT + (1-X2)*C2BNmc/mr2mc_PRODT + (1-X3)*C3BNmc/mr2mc_PRODT + (1-X4)*C1BNmc/mr2mc_PRODT)
              Q_cbl(LIOH%loc) = -Q_init(LIOH%loc)* ( (1-X1)*( (1-beta_terp)*C1LHlmc + beta_terp*C1LHhmc)/mr2mc_PRODT + (1-X2)*( (1-beta_terp)*C2LHlmc + beta_terp*C2LHhmc)/mr2mc_PRODT +  (1-X3)*( (1-beta_terp)*C3LHlmc + beta_terp*C3LHhmc)/mr2mc_PRODT + (1-X4)*( (1-beta_terp)*C4LHlmc + beta_terp*C4LHhmc)/mr2mc_PRODT)
              Q_cbl(LIO3%loc) = -Q_init(LIO3%loc)* ( (1-X1)*((1-beta_terp)*C1LOlmc + beta_terp*C1LOhmc)/mr2mc_PRODT + (1-X2)*((1-beta_terp)*C2LOlmc  + beta_terp*C2LOhmc)/mr2mc_PRODT + (1-X3)*((1-beta_terp)*C3LOlmc &
                      + beta_terp*C3LOhmc)/mr2mc_PRODT + (1-X4)*((1-beta_terp)*C4LOlmc  + beta_terp*C4LOhmc)/mr2mc_PRODT)
              Q_cbl(LINO3%loc) = -Q_init(LINO3%loc)* ( (1-X1)*C1LNmc/mr2mc_PRODT + (1-X2)*C2LNmc/mr2mc_PRODT + (1-X3)*C3LNmc/mr2mc_PRODT + (1-X4)*C4LNmc/mr2mc_PRODT)
              Q_cbl(ISO3%loc) = -Q_init(ISO3%loc)* ( (1-X2)*((1-beta_iso)*C2IOlmc + beta_iso*C2IOhmc)/mr2mc_PRODI + (1-X3)*((1-beta_iso)*C3IOlmc + beta_iso*C3IOhmc)/mr2mc_PRODI + (1-X4 )*((1-beta_iso)*C4IOlmc + beta_iso*C4IOhmc)/ mr2mc_PRODI)
              Q_cbl(ISNO3%loc) = -Q_init(ISNO3%loc)* ( (1-X2)*C2INmc/mr2mc_PRODI + (1-X3)*C3INmc/mr2mc_PRODI + (1-X4)*C4INmc/mr2mc_PRODI)
            endif
         
          case default
            if (t==1) print *,'Emission function for species ', i, 'is undefined'
            Q_cbl(i) = 0.
        end select
      enddo

      if (lCO2Ags) then
        if(lchem .and. (CO2%loc .gt. 0) ) then
          Q_cbl(CO2%loc) = wco2 * 1000 !conversion from ppm m/s to ppb m/s
        else
          wcs            = wco2 * 1000 !conversion from ppm m/s to ppb m/s
        endif
      endif
      !HGO if lCO2Ags, then calculate wCO2 @ CO2%loc apart lchem and CO2%loc
      !.gt. 0
      if (lBVOC) then       
      !deposition if BVOC emissions are interactively calculated and TERP emissions are divided over a-pinene, b-pinene and limonene
      !(2017), based on observations at SOAS campaign (Nagori et al., 2019)
        !if(lchem .and. (TERP%loc .gt. 0) ) then
        if(lchem .and. (APIN%loc .gt. 0) ) then
          Q_cbl(APIN%loc) = F_TERP*0.45
        endif
        if(lchem .and. (BPIN%loc .gt. 0) ) then
          Q_cbl(BPIN%loc) = F_TERP*0.45
        endif
        if(lchem .and. (LIMO%loc .gt. 0) ) then
          Q_cbl(LIMO%loc) = F_TERP*0.1
        endif
        if(lchem .and. (ISO%loc .gt. 0) ) then
          Q_cbl(ISO%loc) = F_ISO
        endif
      endif
    endif !(c_wth)

    if (wthetav .lt. 0.0) wthetav = 0.0

    if (t==1) then ! write in initial_chem the used fluxes
      write (46,*)' '
      write (46,*)'  Q_cbl(INERT%loc)= Q_cbl * sin(pi * (sec-daytime_start)/daylenght)'
      write (46,*)'  Q_cbl(RH%loc)   = Q_cbl * sin(pi * (sec-daytime_start)/daylenght)'
    endif

!   subsidence velocity (large-scale advection)
    ws=-wsls*zi(1)
!
    lsq     = 0.0
    lstheta = 0.0
    if ((sec .ge. starttime_adv_T) .and. (sec .lt. endtime_adv_T)) then
      lstheta = advtheta
    else if (sec .ge. endtime_adv_T) then
      tau_ta = 600.                                                                                           
      advtheta_now = advtheta * exp((endtime_adv_T-sec)/ tau_ta)
      lstheta = advtheta_now
   endif
    
    if ((sec .ge. starttime_adv_q) .and. (sec .lt. endtime_adv_q)) then
      lsq = advq
    else if (sec .ge. endtime_adv_q) then
      tau_ta = 600.                                                                                           
      advq_now = advq * exp((endtime_adv_q-sec)/ tau_ta)
      lsq = advq_now
   endif    
!Add advection relaxation term (Nagori 2019)

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!   shallow cumulus
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(lscu .and. sec>7200) then    !start cloud formation 2 hours after model starts (Nagori 2017)
      ! 1. calculate saturation deficit at zi 
      Ptop     = (pressure*100.) / exp((g * zi(1))/(Rd * thetam(1)))
      Ttop     = (thetam(1) / (((p0*1000) / Ptop) ** (Rd / Cp)))  
      estop    = 611. * exp((Lv / Rv) * ((1. / 273.15)-(1. / Ttop)))
      qstop    = (0.622 * estop) / Ptop
      qqs      = (qm(1)/1000.) - qstop

      ! 2. iterative solution LCL
        !! solution differs from saturation routine MXLCH !!
        !! Here we use mixed-layer averaged values for theta
        !! saturation routine uses surface temp for start parcel
        !! and thus yields a higer lcl
      if(t.eq.1) then   ! Put random guess for lcl, RHlcl
        zlcl   = zi(1)
        RHlcl  = 0.5;
      else              ! Perturb RHlcl to start iteration
        RHlcl  = 0.9998;
      end if
      do ii=1,30        ! Puts limit on max amount of iterations, just to be sure...
        zlcl   = zlcl + (1.-RHlcl) * 1000.
        Ptop   = (pressure*100.) / exp((g * zlcl)/(Rd * thetam(1)))
        Ttop   = (thetam(1) / (((p0*1000) / Ptop) ** (Rd / Cp)))  
        estop  = 611. * exp((Lv / Rv) * ((1. / 273.15)-(1. / Ttop)))
        etop   = (qm(1)/1000.) * Ptop / 0.622 
        RHlcl  = etop / estop;
        if((RHlcl .ge. 0.9999) .and. (RHlcl .le. 1.0001)) exit
      end do

      ! 3. iterate wqM and ac
      wqm      = 0.
      q2       = 1e-10
      wm       = 0.
      do ii=1,30        ! to-do: smarter iteration
        if(dq(1) < 0 .and. wstar > 0) q2 = -((-we*(dq(1)/1000.))+wqm) * (dq(1)/1000.) * zi(1) / (dz * wstar)
        if(q2 < 0)    q2 = 1e-10
        ac     = max(0.,(0.5 + (0.36 * atan(1.55 * (((qm(1)/1000.) - qstop) / sqrt(q2))))))
        wm     = ac * wstar
        wqm    = ca * sqrt(q2) * wm
      end do

      ! 4. relaxe dz if(lrelaxdz)
      if(ac>0 .or. ((zlcl - zi(1)) .lt. 300.0)) then
       if(lrelaxdz) then
        dz    = dz + ((zlcl - zi(1)) - dz) * dtime/tau
       else
        dz    = min(max(zlcl - zi(1),100.0),500.0)
       end if
      end if
      
      ! 5. iterate wscm2 (ventilation due to clouds); wscm2 = ventilation flux, dsc = jump of scalar, wsce = entrainment flux (2017)
      sig_sc2   = 1e-5        !variance of scalar at top of ABL
      ac     = max(0.,(0.5 + (0.36 * atan(1.55 * (((qm(1)/1000.) - qstop) / sqrt(q2))))))
      wm     = ac * wstar
        if(lchem) then
          do k=1, nchsp
            wscm2(k)=.0                           
            dsc       = c_ft(k)-c_cbl(k)           
            wsce      = we * dsc                 
           
            if(wthetav > 0.)  then
              sig_sc2 = (wsce-wscm2(k)) * dsc * zi(1) / (dz * wstar)
            else
               sig_sc2 = 0.
            endif   

            if(sig_sc2 < 0.) then
              sig_sc2 = 1e-5
            endif
 
            if(dsc<0.)then
              wscm2(k)   =  sqrt(sig_sc2) * wm
            else
              wscm2(k) = 0.
            endif
           enddo
        endif
    end if

!   dynamics eulerien time-step
!   Potential temperature

!   switch for applying different lapse rate 
    if (lgamma .eqv. .true.) then
       if (zi(2)>hcrit) then
         gamma = gamma2
       endif
    endif      
    
    if (lencroachment .eqv. .false.) then
            
!   jump depending on the virtual temperature jump
      dthetav = dtheta(1) + &
      1.e-03*rvrd*(qm(1)*dtheta(1)+thetam(1)*dq(1)+dtheta(1)*dq(1))
      we = (beta*wthetav/dthetav)
      if (lenhancedentrainment) then
        we = we + 5. * (ustar**3.) * thetav / (g * zi(1) * dthetav)
      endif
      wf = Rdistr * DeltaF / ( rho * Cp * dtheta(1) )
!      we = we + wf  ! Aerosol absorption basically undoes part of we

      if (lfixedtroposphere) wsubs = ws - wm
      
      zi(2)  =zi(1)+(we+wf+ws-wm)*dtime            !eq (8)
      hcrit    = hcrit + ws * dtime
      dtheta(2)=dtheta(1)+ ((gamma*(we + wf + wsubs))- &                            !extra subsidence by wm only affects dtheta (through gamma) if fixed tropospheric profiles are considered
              (1/(zi(1)+inf))*(wthetas+(we+wf)*dtheta(1) - DeltaF/(rho*Cp) ))*dtime !eq. (3)
      if (.not. ladvecFT) dtheta(2) = dtheta(2) - lstheta*dtime
      thetam(2)=thetam(1)+ &
            ((1/(zi(1)+inf))*(wthetas+(we+wf)*dtheta(1) - DeltaF/(rho*Cp) ))*dtime + &        !eq  (1)
            lstheta*dtime
    else
!     Calculation boundary layer height.
!     Based on Garratt (6.18) page 155 and
!     on the article of Betts 1973 and 1974
      beta   = 0.
      we     = 0. 

      zi(2)   = zi(1) + ((1/(gamma*(zi(1)+inf)))*((1+beta)*wthetav) &
              + ws)*dtime
     
      hcrit   = hcrit + ws * dtime
      dtheta(2)=0.
      thetam(2)=thetam(1)+ &
            ((1/(zi(1)+inf))*wthetas)*dtime + lstheta*dtime         !eq  (1)
    endif

!    thetam(2)=thetam(1)+ &
!            ((1/(zi(1)+inf))*(wthetas+we*dtheta(1)))*dtime          !eq  (1)

!    we = (zi(2)-zi(1))/dtime - ws

! ---- Specific humidity
    if (beta /= 0) then
        wqe = -(we+wf)*dq(1)
        dq(2)=dq(1)+ ((gammaq*(we + wf + wsubs))-(1/(zi(1)+inf))*(wqs - wqe - (wqm*1000.)))*dtime   !eq. (3)
        if (.not. ladvecFT) dq(2) = dq(2) - lsq*dtime
    else
        wqe  = 0.
        dq(2)= 0.
    endif

    qm(2)=qm(1)+((1/(zi(1)+inf))*(wqs - wqe - (wqm*1000.)))*dtime+lsq*dtime          !eq  (1)
! -----------Carbon Dioxide

    if (beta /= 0) then
      wce = -(we+wf)*dc(1)
      dc(2)=dc(1)+ ((gammac*(we + wf + wsubs))-(1/(zi(1)+inf))*(wcs - wce))* dtime   !eq. (3)
    else
      wce    = 0.
      dc(2)  =0.
    endif

    cm(2)=cm(1)+((1/(zi(1)+inf))*(wcs - wce))*dtime        !eq  (1)

! ---- Windspeed

    if (beta /= 0) then
      uwe = -(we+wf)*du(1)
      vwe = -(we+wf)*dv(1)
      du(2)=du(1)+ ((gammau*(we + wf + wsubs))-(1/(zi(1)+inf))*(uws - uwe)+ f*dv(1))* dtime                 !eq. (3)
      dv(2)=dv(1)+ ((gammav*(we + wf + wsubs))-(1/(zi(1)+inf))*(vws - vwe) - f*du(1))*dtime                 !eq. (3)
    else
      uwe = 0.
      vwe = 0.
      du(2)=0.
      dv(2)=0.
    endif

    um(2)=um(1)+(-f*(dv(1))+(1/(zi(1)+inf))*(uws-uwe))*dtime
    vm(2)=vm(1)+(f*(du(1))+(1/(zi(1)+inf))*(vws-vwe))*dtime

!   closure assumption
    wthetae= - (we+wf) * dtheta(1)

    if (wqs.ne.0.) then
        betaq = -wqe/wqs
    else
        betaq=0.
    endif

    if (wcs.ne.0.) then
        betac = -wce/wcs
    else
        betac=0.
    endif

!   Calculate saturation level and vertical profiles of the air temperature,
!   the dry and moist adiabatic lapse rate of a lifted parcel and the total and
!   saturated specific humidity

!#############################################################################################################
!    if (mod(tt,aver3) == 0.) then    !calculate it only every 300 s
    if (mod(floor(sec),aver3) == 0 ) then    !this is also wrong if dtime < .5 sec

!   Define surface values

      T_s = thetam(2)              !Mixed-layer potential temp in K
      qt_s = qm(2)/1000            !Mixed-layer specific humidity in g/g

      z_s = 0
      es_s = e0 * exp((Lv/Rv)*((1/Tabs)-(1/T_s)))
      qs_s = (epsilon*es_s)/p0
      Tparcel_s = T_s


!     Create new file named saturation

      if (zi(2).lt.lcl) then

        open (14, file=trim(outdir)//dirsep//'saturation')

        write(14,'(a4)') 'SATU'
        write(14,'(I3)') sat_lev
        write(14,'(6a14)') 'z (m)','p (kPa)','T (K)','Tparcel (K)', &
        'qt (g/g)','qs (g/g)'

        write(14,'(6F14.4)') z_s,p0,T_s,Tparcel_s,qt_s,qs_s

!       Do loop to calculate variables

        np = (p0-75)/50.
        p = p0

        do i = 0, sat_lev-1

          p = p-np

!         If qs > qt then no saturation
          if ((qs_s-qt_s)> 0.0) then

            Tair = T_s*((p/p0)**(Rd/Cp))

            if (qm(2)>0) then
              qt = qm(2)/1000.        !Mixed-layer spec. humidity in g/Kg
            else
              qt = 0.000001
            end if

            z = (Tair - T_s)/gammad

!           Use of Clausius-Clapeyron

            e_saturation = (qt*p)/epsilon
            Td = ((1/Tabs)-((Rv/Lv)*log(e_saturation/e0)))**(-1)
            es = e0*exp((Lv/Rv)*((1/Tabs)-(1/Tair)))
            qs = (epsilon*es)/p

!           Dry adiabatic lapse rate, the last height z as a dry adiabat
!           represents saturation level

            Tparcel = T_s + gammad*z
!           print *, "Dry adiabat"

          else !((qs_s-qt_s)> 0.0)

            a = a + 1

            Tair = T_s*((p/p0)**(Rd/Cp))
            es = e0*exp((Lv/Rv)*((1/Tabs)-(1/Tair)))
            qs = (epsilon*es)/p

!           There is no information on the liquid water content (ql), qt and Tv are
!            therefore not conserved anymore in this program where the mixed-layer
!            model is used. The calculated height z is not correct either.

            qt = qs
            z = (Tair - T_s)/gammad

            if (a == 1) then
              lcl = z
            end if

!           Moist adiabatic lapse rate:

            gammam = gammad*((1+((qs*Lv)/(Rd*Tair)))/(1+(((Lv**2)*qs)/(Cp*Rv*(Tair**2)))))
            Tparcel = Tair + (z-lcl)*gammam

!           print *, "Moist adiabat"

          end if !((qs_s-qt_s)> 0.0)

          qt_s = qt
          qs_s = qs

          write(14,'(6F14.4)') z,p,Tair,Tparcel,qt,qs

        end do !i = 0, sat_lev-1

        a = 0

!       write(14,'(a25)')

        sat_time = printhour
        sat_height = lcl

        write(14,'(a6,F8.1, a15,F8.2,a3)') 'LCL at', &
                sat_height,'m reached after',sat_time,'hours running'
        close(14)

      else !(zi(2).lt.lcl)

        if (n == 0) then
          print *, 'saturation reached'
          n = 1
        end if

      end if !(zi(2).lt.lcl)

    end if !(mod(tt,aver3) == 0.)

!   chemistry

!   t is the number of sec since the beginning of the run (for dtime=1)
!   ttt= total number of seconds
    ttt=nint(hour*3600)+(t*dtime)
!   itme= number of full hours
    itime=nint(ttt/3600)
!   mins number of minutes of the current hour
    mins=nint((ttt-itime*3600)/60)
!   isec number of second of the current min
    isec=ttt-itime*3600 -mins*60
!   runing time in hours
    tday = 1.0*day

    if (lchem .and. (OAbg%loc .gt. 0))  then
      adv_chem_cbl(OAbg%loc) = adv_chem_cbl(OAbg%loc)* sin(pi * (sec - starttime_chem)/(endtime_chem - starttime_chem))
    endif


!    emission and entrainment each timestep + shallow cumulus (wscm2) added

    if(lchem) then
      do k=1, nchsp
        E(k)=-(we+wf)*(c_ft(k)-c_cbl(k))
        c_cbl(k)=c_cbl(k)+(1/zi(1))*(Q_cbl(k)-E(k)-wscm2(k))*dtime
        c_cbl(k)=c_cbl(k)+adv_chem_cbl(k)*dtime
        c_ft( k)=c_ft( k)+adv_chem_ft(k )*dtime
        c_cbl(k) = max(0.0, c_cbl(k) )
        c_ft( k) = max(0.0, c_ft( k) )

        if(k==45) then
          open  (85, file=trim(outdir)//dirsep//'Q_cbl')
          write (85, *) thour, printhour, (Q_cbl(k))
        endif

       
        if (Q_cbl(k) == 0.) then
          beta_ft(k)=0
        else
          beta_ft(k)=E(k)/Q_cbl(k)
        endif
      enddo
       
      if(lchconst .eqv. .true.)then
        temp_cbl = t_ref_cbl
        temp_ft  = t_ref_ft

        c_cbl(H2O%loc)= q_ref_cbl* MW_Air/MW_H2O * 1e6
        c_ft(H2O%loc) = q_ref_ft * MW_Air/MW_H2O * 1e6
      else
        temp_cbl = thetam(2)
        temp_ft  = thetam(2) + dtheta(2)

        c_cbl(H2O%loc)=qm(2)* MW_Air/MW_H2O * 1e6
        c_ft(H2O%loc)=(qm(2)+ dq(2))* MW_Air/MW_H2O * 1e6
      endif

      if(lcomplex) then
        call calc_K_mozart(pressure, real(temp_cbl + (zi(2)/2) * gammad), real(temp_ft + zi(2) * gammad) ) !Accounting for dry adiabatic lapse rate
      else
        call calc_K_simple(pressure, temp_cbl + (zi(2)/2) * gammad, temp_ft + zi(2) * gammad) !Accounting for dry adiabatic lapse rate
      endif

!     solve the chemistry gas phase
!     Boundary layer variables with beta closure
      c_current = c_cbl
      if(lcomplex) then
        call iter_mozart(CBL, c_cbl, c_current, dtime)
      else
        call iter_simple(CBL, c_cbl, c_current, dtime)
      endif

!     Free troposphere layer variables
      c_current = c_ft
      if(lcomplex) then
        call iter_mozart(FT, c_ft, c_current, dtime)
      else
        call iter_simple(FT, c_ft, c_current, dtime)
      endif

      if (RC(R_NO2%loc)%Keff_cbl == 0 .OR. c_cbl(NO2%loc)==0) then
          photo = 0
      else
          photo = (RC(R_NO%loc)%Keff_cbl * c_cbl(NO%loc)*c_cbl(O3%loc))/(RC(R_NO2%loc)%Keff_cbl*c_cbl(NO2%loc))
      endif
    endif !(lchem)

    ! Calculate Gamma IEPOX and heterogeneous reaction rate through resistor model as per Gaston (2014) and Hu et al (2016), used in
    ! Nagori et al., 2019
    if (lchem) then
      Rgas_g = 8.2057e-5                                  ! Ideal gas constant
      Rp = 140e-9                                         ! surface area weighted mean geometric radius of the particles in m 
      Rc_gIpx = 100e-9                                    ! radius of the aqueous core in m//this is the surface-area weighted particle radius of the inorganic core
      epsln = 0.1                                         ! factor from 0 to 1 used to decrease the Henry's law constant and diffusivity in the organic
      Hconc = 0.15                                        ! [H+] is the proton concentration (M)
      aH = 1.                                             ! proton activity
      Cnuc = 0.05                                         ! concentration of nucleophile (nitrate or sulfate), Gaston calculated it from E-AIM II; Unit is M
      cga = 0.59                                          ! concentration of HSO4 (this term is taken out for nitrates) Unit:  Gaston calculated it from E-AIM II, Unit is M
      Dgas = 1e-5                                         ! diffusivity in the gas phase in m2/s--from Antilla
      Daq = 1e-9                                          ! diffusivity in aqueous core in m2/s--value from Antilla,Griffiths, and Mentel
      Haq = 1.7e11                                        ! Henry’s law coefficient in the aqueous phase mol/m3/atm--value from Eddingsaas
      kH = 0.036                                          ! value from Eddingsaas (M-1s-1)
      knuc = 2e-4                                         ! value from Eddingsaas for sulfate and nitrate M-1 s-1
      kga = 7.3e-4                                        ! reaction rate due to HSO4, value from Eddingsaas for HSO4, unit M-1 s-1
      k_aq = (kH*Hconc)+(knuc*Cnuc*aH)+(kga*cga)          ! first order rate loss of IEPOX 
      SurfaceArea = 300e-6                                ! surface area
      
      omega = (8*8.3145*T_hh/(pi*0.118))**(1./2.)         ! Omega=Meanspeed*0.01; MeanSpeed=(8RT/pi*MolecularWeight) unit of molecular weight to be  kg/mol 
      inv_alpha = 1./alpha_gmaIPX
      HorgDorg = epsln*Haq*Daq                            ! Henry's law constant x diffusivity of the organic coating as a function of the aqueous phase in mol m2/s/L/atm;
      k_org = k_aq
      Dorg = 5e-10
      if(Rc_gIpx .EQ. 0) then
        q_aq = 0.
        q_org = Rp*((k_org/Dorg)**(1./2.))
        q_org_new = 0.
        num = 0.
        denom = 1.
        h = 0.
        F_gIpx = ((1./tanh (q_org)) + h)/(1 + ((1./tanh (q_org))*h))
        
      elseif(Rc_gIpx .GT. 0) then
        
        q_aq = Rc_gIpx*((k_aq/Daq)**(1./2.))
        q_org = Rp*((k_org/Dorg)**(1./2.))
        q_org_new = (Rc_gIpx/Rp)*q_org
        num = ((Haq*Daq)/(HorgDorg))*((q_aq*(1./tanh(q_aq)))-1)-((q_org_new*(1./tanh(q_org_new)))-1.)
        denom = ((Haq*Daq)/(HorgDorg))*((q_aq*(1./tanh (q_aq)))-1)-((q_org_new*(tanh (q_org_new)))-1.)
        h = (-tanh(q_org_new))*(num/denom)
        F_gIpx = ((1./tanh(q_org)) + h)/(1. + ((1./tanh(q_org))*h))
      endif
      ! //Resistor model accounting for gas phase diffusion, diffusion and reactivity in organic coating and aqueous core
      inv_gas = (Rp*omega)/(4*Dgas)
      inv_Torq_aq = (omega*Rp)/(4*Rgas_g*T_hh*HorgDorg*((q_org*F_gIpx)-1.))
      inv_gamma = inv_gas + inv_alpha + inv_Torq_aq
      gamma_uptake = 1./inv_gamma

      k_het = 1./4.*SurfaceArea*gamma_uptake*omega   ! heterogeneous reaction rate
      IEPOX_SOA = k_het * Q_cbl(IEPOX%loc)
      Q_cbl(IEPOXSOA%loc) = 0.11*IEPOX_SOA
      
      !print *,gamma_uptake
    endif
      
    ! Gas/Particle-partitioning using the Volatility Basis Set (Donahue'06, ES&T) with 4 bins
    if (lvbs) then
        kISORO2NO  = RC(R_ISORO2NO%loc)%Keff_cbl / convcbl  ! convert from cm3/molec*s to 1/ppb/s
        kISORO2HO2 = RC(R_ISORO2HO2%loc)%Keff_cbl / convcbl ! convert from cm3/molec*s to 1/ppb/s
        kTERPRO2NO = 2.2e-12 ! reaction rate for NO + BPRO2 at 298K (from MCM)
        kTERPRO2HO2= 2.1e-11 ! reaction rate for HO2 + BPRO2 at 298K (from MCM)
        
        select case(low_high_NOx)    ! low/high NOx branching of SVOC yields (Lane'08, ES&T)
          case(0) ! calculate branching ratio beta = rate_of_RO2+NO/(rate_of_RO2+NO + rate_of_RO2+HO2)
            beta_iso    = (kISORO2NO*c_cbl(NO%loc)) / (kISORO2NO*c_cbl(NO%loc) + kISORO2HO2*c_cbl(HO2%loc))
            beta_terp   = (kTERPRO2NO*c_cbl(NO%loc)) / (kTERPRO2NO*c_cbl(NO%loc) + kTERPRO2HO2*c_cbl(HO2%loc))
          case(1) ! only low NOx channel yields
            beta_terp = 0. 
            beta_iso  = 0.
          case(2) ! only high NOx channel yields
            beta_terp = 1. 
            beta_iso  = 1.
          case(3) ! high NOx channel yields for isoprene; both yields for monoterpene
            beta_terp = (kTERPRO2NO*c_cbl(NO%loc)) / (kTERPRO2NO*c_cbl(NO%loc) + kTERPRO2HO2*c_cbl(HO2%loc))
            beta_iso  = 1. 
          case default
            if (t==1) print *,'Flag for the function of high/low NOx branching is invalid: ',low_high_NOx
            stop 'change function'
        end select

 !SVOC per oxidant reaction (Easier to track per oxidant contribution) in the four volatilty basis set. Nagori et al., 2019

        SVOC_alp1 = (1-beta_terp)*C1AHlmc + beta_terp*C1AHhmc + (1-beta_terp)*C1AOlmc + beta_terp*C1AOhmc +  C1ANmc  
        SVOC_alp2 = (1-beta_terp)*C2AHlmc + beta_terp*C2AHhmc + (1-beta_terp)*C2AOlmc + beta_terp*C2AOhmc + C2ANmc 
        SVOC_alp3 = (1-beta_terp)*C3AHlmc + beta_terp*C3AHhmc + (1-beta_terp)*C3AOlmc + beta_terp*C3AOhmc + C3ANmc
        SVOC_alp4 = (1-beta_terp)*C4AHlmc + beta_terp*C4AHhmc + (1-beta_terp)*C4AOlmc + beta_terp*C4AOhmc + C4ANmc 
       
        SVOC_bep1 =  (1-beta_terp)*C1BHlmc + beta_terp*C1BHhmc + (1-beta_terp)*C1BOlmc + beta_terp*C1BOhmc + C1BNmc
        SVOC_bep2 = (1-beta_terp)*C2BHlmc + beta_terp*C2BHhmc + (1-beta_terp)*C2BOlmc + beta_terp*C2BOhmc + C2BNmc
        SVOC_bep3 =  (1-beta_terp)*C3BHlmc + beta_terp*C3BHhmc  + (1-beta_terp)*C3BOlmc + beta_terp*C3BOhmc + C3BNmc
        SVOC_bep4 = (1-beta_terp)*C4BHlmc + beta_terp*C4BHhmc + (1-beta_terp)*C4BOlmc + beta_terp*C4BOhmc + C4BNmc

        SVOC_lim1 = (1-beta_terp)*C1LHlmc + beta_terp*C1LHhmc + (1-beta_terp)*C1LOlmc   + beta_terp*C1LOhmc  + C1LNmc 
        SVOC_lim2 = (1-beta_terp)*C2LHlmc + beta_terp*C2LHhmc + (1-beta_terp)*C2LOlmc + beta_terp*C2LOhmc + C2LNmc 
        SVOC_lim3 = beta_terp*C3LHhmc  + (1-beta_terp)*C3LOlmc + beta_terp*C3LOhmc + C3LNmc
        SVOC_lim4 = beta_terp*C4LHhmc +  (1-beta_terp)*C4LOlmc + beta_terp*C4LOhmc  + C4LNmc
        
        SVOC_iso2 =  (1-beta_iso)*C2IOlmc + beta_iso*C2IOhmc +  C2INmc
        SVOC_iso3 = (1-beta_iso)*C3IOlmc + beta_iso*C3IOhmc + C3INmc
        SVOC_iso4 = (1-beta_iso)*C4IOlmc + beta_iso*C4IOhmc + C4INmc
        
        alp_OH1 = (1-beta_terp)*C1AHlmc + beta_terp*C1AHhmc
        alp_OH2 = (1-beta_terp)*C2AHlmc + beta_terp*C2AHhmc
        alp_OH3 = (1-beta_terp)*C3AHlmc + beta_terp*C3AHhmc
        alp_OH4 = (1-beta_terp)*C4AHlmc + beta_terp*C4AHhmc
        alp_O31 = (1-beta_terp)*C1AOlmc + beta_terp*C1AOhmc
        alp_O32 = (1-beta_terp)*C2AOlmc + beta_terp*C2AOhmc
        alp_O33 = (1-beta_terp)*C3AOlmc + beta_terp*C3AOhmc
        alp_O34 = (1-beta_terp)*C4AOlmc + beta_terp*C4AOhmc

        bep_OH1 = (1-beta_terp)*C1BHlmc + beta_terp*C1BHhmc
        bep_OH2 = (1-beta_terp)*C2BHlmc + beta_terp*C2BHhmc
        bep_OH3 = (1-beta_terp)*C3BHlmc + beta_terp*C3BHhmc
        bep_OH4 = (1-beta_terp)*C4BHlmc + beta_terp*C4BHhmc
        bep_O31 = (1-beta_terp)*C1BOlmc + beta_terp*C1BOhmc
        bep_O32 = (1-beta_terp)*C2BOlmc + beta_terp*C2BOhmc
        bep_O33 = (1-beta_terp)*C3BOlmc + beta_terp*C3BOhmc
        bep_O34 = (1-beta_terp)*C4BOlmc + beta_terp*C4BOhmc

        lim_OH1 = (1-beta_terp)*C1LHlmc + beta_terp*C1LHhmc
        lim_OH2 = (1-beta_terp)*C2LHlmc + beta_terp*C2LHhmc
        lim_OH3 = (1-beta_terp)*C3LHlmc + beta_terp*C3LHhmc
        lim_OH4 = (1-beta_terp)*C4LHlmc + beta_terp*C4LHhmc
        
        lim_O31 = (1-beta_terp)*C1LOlmc   + beta_terp*C1LOhmc
        lim_O32 = (1-beta_terp)*C2LOlmc   + beta_terp*C2LOhmc
        lim_O33 = (1-beta_terp)*C3LOlmc   + beta_terp*C3LOhmc
        lim_O34 = (1-beta_terp)*C4LOlmc   + beta_terp*C4LOhmc

        isoo32 = (1-beta_iso)*C2IOlmc + beta_iso*C2IOhmc
        isoo33 = (1-beta_iso)*C3IOlmc + beta_iso*C3IOhmc
        isoo34 = (1-beta_iso)*C4IOlmc + beta_iso*C4IOhmc

        C1mc   =  (1-beta_terp)*C1AHlmc + beta_terp*C1AHhmc + (1-beta_terp)*C1BHlmc +  &
        beta_terp*C1BHhmc + (1-beta_terp)*C1LHlmc + beta_terp*C1LHhmc + (1-beta_terp)*C1AOlmc & 
        + beta_terp*C1AOhmc + (1-beta_terp)*C1BOlmc + beta_terp*C1BOhmc + (1-beta_terp)*C1LOlmc &
        + beta_terp*C1LOhmc + C1ANmc + C1BNmc + C1LNmc
        C2mc   = (1-beta_terp)*C2AHlmc + beta_terp*C2AHhmc + (1-beta_terp)*C2BHlmc + &
        beta_terp*C2BHhmc + (1-beta_terp)*C2AOlmc + beta_terp*C2AOhmc + (1-beta_terp)*C2BOlmc &
        + beta_terp*C2BOhmc + (1-beta_terp)*C2LHlmc &
        + beta_terp*C2LHhmc + (1-beta_iso)*C2IOlmc + beta_iso*C2IOhmc + (1-beta_terp)*C2LOlmc &
        + beta_terp*C2LOhmc + C2ANmc + C2BNmc + C2INmc + C2LNmc
        C3mc   = (1-beta_terp)*C3AHlmc + beta_terp*C3AHhmc + (1-beta_terp)*C3BHlmc + beta_terp*C3BHhmc &
        + (1-beta_terp)*C3AOlmc + beta_terp*C3AOhmc + (1-beta_terp)*C3BOlmc + beta_terp*C3BOhmc &
        + beta_terp*C3LHhmc &
        + (1-beta_iso)*C3IOlmc + beta_iso*C3IOhmc + (1-beta_terp)*C3LOlmc + beta_terp*C3LOhmc + C3ANmc + C3BNmc + C3INmc + C3LNmc
        C4mc   = (1-beta_terp)*C4AHlmc + beta_terp*C4AHhmc + (1-beta_terp)*C4BHlmc + beta_terp*C4BHhmc &
        + (1-beta_terp)*C4AOlmc + beta_terp*C4AOhmc + (1-beta_terp)*C4BOlmc + beta_terp*C4BOhmc &
        + beta_terp*C4LHhmc &
        + (1-beta_iso)*C4IOlmc + beta_iso*C4Iohmc + (1-beta_terp)*C4LOlmc + beta_terp*C4LOhmc &
        + C4ANmc + C4BNmc + C4INmc + C4LNmc
        ! T-dependence following Clausius-Clapeyron (Sheehan and Bowman '01)
        C1_st    = C10_st * (T0/T_hh)*exp((dHvap/Rgas)*(1./T0-1./T_hh)) 
        C2_st    = C20_st * (T0/T_hh)*exp((dHvap/Rgas)*(1./T0-1./T_hh))
        C3_st    = C30_st * (T0/T_hh)*exp((dHvap/Rgas)*(1./T0-1./T_hh))
        C4_st    = C40_st * (T0/T_hh)*exp((dHvap/Rgas)*(1./T0-1./T_hh))
        m        = 0 
        Coait    = OAbgb  
        epsi = 1.
        do while (epsi>1e-10) ! set iteration tolerance here
             X1          = 1./((1.+(C1_st/Coait)))                                                                                         
             X2          = 1./((1.+(C2_st/Coait)))
             X3          = 1./((1.+(C3_st/Coait)))
             X4          = 1./((1.+(C4_st/Coait)))
             Coait_new   = X1*C1mc + X2*C2mc + X3*C3mc + X4*C4mc + OAbgb                  
             epsi         = abs(Coait-Coait_new) ! optimization parameter
             Coait       = Coait_new      
             m           = m+1  
             if (m>100) then
                write (*,*) 'Warning: no convergence (BL)!'  
             endif                   
        enddo ! while  
        Coa = Coait 


        


        ! G/P-partitioning FT
        OAbgb_ft = c_ft(OAbg%loc)* mr2mc_OAbg  ! convert to ug m-3
        APOHmc_ft = c_ft(APOH%loc) * mr2mc_PRODT  
        APO3mc_ft = c_ft(APO3%loc) * mr2mc_PRODT  
        APNO3mc_ft = c_ft(APNO3%loc) * mr2mc_PRODT  
        BPOHmc_ft = c_ft(BPOH%loc) * mr2mc_PRODT  
        BPO3mc_ft = c_ft(BPO3%loc) * mr2mc_PRODT  
        BPNO3mc_ft = c_ft(BPNO3%loc) * mr2mc_PRODT  
        LIOHmc_ft = c_ft(LIOH%loc) * mr2mc_PRODT  
        LIO3mc_ft = c_ft(LIO3%loc) * mr2mc_PRODT  
        LINO3mc_ft = c_ft(LINO3%loc) * mr2mc_PRODT
        ISO3mc_ft = c_ft(ISO3%loc) * mr2mc_PRODI
        ISNO3mc_ft = c_ft(ISNO3%loc) * mr2mc_PRODI
        C1AHlmc_ft   = alpha1_APOHlow  * APOHmc_ft  ! low NOx yield per ox per terp A=apin H=OH
        C2AHlmc_ft   = alpha2_APOHlow  * APOHmc_ft 
        C3AHlmc_ft   = alpha3_APOHlow  * APOHmc_ft 
        C4AHlmc_ft   = alpha4_APOHlow  * APOHmc_ft 
        C1BHlmc_ft   = alpha1_BPOHlow  * BPOHmc_ft  ! low NOx yield per ox per terp B=bpin H=OH
        C2BHlmc_ft   = alpha2_BPOHlow  * BPOHmc_ft 
        C3BHlmc_ft   = alpha3_BPOHlow  * BPOHmc_ft 
        C4BHlmc_ft   = alpha4_BPOHlow  * BPOHmc_ft 
        C1AOlmc_ft   = alpha1_APO3low  * APO3mc_ft  ! low NOx yield per ox per terp A=apin O=O3
        C2AOlmc_ft   = alpha2_APO3low  * APO3mc_ft 
        C3AOlmc_ft   = alpha3_APO3low  * APO3mc_ft 
        C4AOlmc_ft   = alpha4_APO3low  * APO3mc_ft 
        C1BOlmc_ft   = alpha1_BPO3low  * BPO3mc_ft  ! low NOx yield per ox per terp B=bpin O=O3
        C2BOlmc_ft   = alpha2_BPO3low  * BPO3mc_ft 
        C3BOlmc_ft   = alpha3_BPO3low  * BPO3mc_ft 
        C4BOlmc_ft   = alpha4_BPO3low  * BPO3mc_ft 
        C1ANmc_ft   = alpha1_APNO3  * APNO3mc_ft  ! only one for NO3 ox; A=apin N=NO3
        C2ANmc_ft   = alpha2_APNO3  * APNO3mc_ft 
        C3ANmc_ft   = alpha3_APNO3  * APNO3mc_ft 
        C4ANmc_ft   = alpha4_APNO3  * APNO3mc_ft 
        C1BNmc_ft   = alpha1_BPNO3  * BPNO3mc_ft  ! only one for NO3 ox; B=bpin N=NO3
        C2BNmc_ft   = alpha2_BPNO3  * BPNO3mc_ft 
        C3BNmc_ft   = alpha3_BPNO3  * BPNO3mc_ft 
        C4BNmc_ft   = alpha4_BPNO3  * BPNO3mc_ft 

        C1AHhmc_ft   = alpha1_APOHhigh  * APOHmc_ft  ! high NOx yield per ox per terp A=apin H=OH
        C2AHhmc_ft   = alpha2_APOHhigh  * APOHmc_ft 
        C3AHhmc_ft   = alpha3_APOHhigh  * APOHmc_ft 
        C4AHhmc_ft   = alpha4_APOHhigh  * APOHmc_ft 
        C1BHhmc_ft   = alpha1_BPOHhigh  * BPOHmc_ft  ! high NOx yield per ox per terp B=bpin H=OH
        C2BHhmc_ft   = alpha2_BPOHhigh  * BPOHmc_ft 
        C3BHhmc_ft   = alpha3_BPOHhigh  * BPOHmc_ft 
        C4BHhmc_ft   = alpha4_BPOHhigh  * BPOHmc_ft 
        C1AOhmc_ft   = alpha1_APO3high  * APO3mc_ft  ! high NOx yield per ox per terp A=apin O=O3
        C2AOhmc_ft   = alpha2_APO3high  * APO3mc_ft 
        C3AOhmc_ft   = alpha3_APO3high  * APO3mc_ft 
        C4AOhmc_ft   = alpha4_APO3high  * APO3mc_ft 
        C1BOhmc_ft   = alpha1_BPO3high  * BPO3mc_ft  ! high NOx yield per ox per terp B=bpin O=O3
        C2BOhmc_ft   = alpha2_BPO3high  * BPO3mc_ft 
        C3BOhmc_ft   = alpha3_BPO3high  * BPO3mc_ft 
        C4BOhmc_ft   = alpha4_BPO3high  * BPO3mc_ft 

        C1LHlmc_ft   = alpha1_LIOHlow  * LIOHmc_ft  ! low NOx yield per ox per terp L=limo H=OH
        C2LHlmc_ft   = alpha2_LIOHlow  * LIOHmc_ft 
        C3LHlmc_ft   = alpha3_LIOHlow  * LIOHmc_ft 
        C4LHlmc_ft   = alpha4_LIOHlow  * LIOHmc_ft 
        C2IOlmc_ft   = alpha2_ISO3low  * ISO3mc_ft  ! low NOx yield per ox for isoprene; O=O3        
        C3IOlmc_ft   = alpha3_ISO3low  * ISO3mc_ft 
        C4IOlmc_ft   = alpha4_ISO3low  * ISO3mc_ft 
        C1LOlmc_ft   = alpha1_LIO3low  * LIO3mc_ft  ! low NOx yield per ox per terp L=limo O=O3
        C2LOlmc_ft   = alpha2_LIO3low  * LIO3mc_ft 
        C3LOlmc_ft   = alpha3_LIO3low  * LIO3mc_ft 
        C4LOlmc_ft   = alpha4_LIO3low  * LIO3mc_ft 
        C2INmc_ft   = alpha2_ISNO3  * ISNO3mc_ft  ! N=NO3        
        C3INmc_ft   = alpha3_ISNO3  * ISNO3mc_ft 
        C4INmc_ft   = alpha4_ISNO3  * ISNO3mc_ft 
        C1LNmc_ft   = alpha1_LINO3  * LINO3mc_ft  ! L=limo N=NO3
        C2LNmc_ft   = alpha2_LINO3  * LINO3mc_ft 
        C3LNmc_ft   = alpha3_LINO3  * LINO3mc_ft 
        C4LNmc_ft   = alpha4_LINO3  * LINO3mc_ft 

        C1LHhmc_ft   = alpha1_LIOHhigh  * LIOHmc_ft  ! high NOx yield per ox per terp L=limo H=OH
        C2LHhmc_ft   = alpha2_LIOHhigh  * LIOHmc_ft 
        C3LHhmc_ft   = alpha3_LIOHhigh  * LIOHmc_ft 
        C4LHhmc_ft   = alpha4_LIOHhigh  * LIOHmc_ft 
        C2IOhmc_ft   = alpha2_ISO3high  * ISO3mc_ft  ! high NOx yield per ox for isoprene; O=O3        
        C3IOhmc_ft   = alpha3_ISO3high  * ISO3mc_ft 
        C4IOhmc_ft   = alpha4_ISO3high  * ISO3mc_ft 
        C1LOhmc_ft   = alpha1_LIO3high  * LIO3mc_ft  ! high NOx yield per ox per terp L=limo O=O3
        C2LOhmc_ft   = alpha2_LIO3high  * LIO3mc_ft 
        C3LOhmc_ft   = alpha3_LIO3high  * LIO3mc_ft 
        C4LOhmc_ft   = alpha4_LIO3high  * LIO3mc_ft

        SVOC_alp1_ft = (1-beta_terp)*C1AHlmc_ft + beta_terp*C1AHhmc_ft + (1-beta_terp)*C1AOlmc_ft + beta_terp*C1AOhmc_ft +  C1ANmc_ft  
        SVOC_alp2_ft = (1-beta_terp)*C2AHlmc_ft + beta_terp*C2AHhmc_ft + (1-beta_terp)*C2AOlmc_ft + beta_terp*C2AOhmc_ft + C2ANmc_ft 
        SVOC_alp3_ft = (1-beta_terp)*C3AHlmc_ft + beta_terp*C3AHhmc_ft + (1-beta_terp)*C3AOlmc_ft + beta_terp*C3AOhmc_ft + C3ANmc_ft
        SVOC_alp4_ft = (1-beta_terp)*C4AHlmc_ft + beta_terp*C4AHhmc_ft + (1-beta_terp)*C4AOlmc_ft + beta_terp*C4AOhmc_ft + C4ANmc_ft 
       
        SVOC_bep1_ft =  (1-beta_terp)*C1BHlmc_ft + beta_terp*C1BHhmc_ft + (1-beta_terp)*C1BOlmc_ft + beta_terp*C1BOhmc_ft + C1BNmc_ft
        SVOC_bep2_ft = (1-beta_terp)*C2BHlmc_ft + beta_terp*C2BHhmc_ft + (1-beta_terp)*C2BOlmc_ft + beta_terp*C2BOhmc_ft + C2BNmc_ft
        SVOC_bep3_ft =  (1-beta_terp)*C3BHlmc_ft + beta_terp*C3BHhmc_ft  + (1-beta_terp)*C3BOlmc_ft + beta_terp*C3BOhmc_ft + C3BNmc_ft
        SVOC_bep4_ft = (1-beta_terp)*C4BHlmc_ft + beta_terp*C4BHhmc_ft + (1-beta_terp)*C4BOlmc_ft + beta_terp*C4BOhmc_ft + C4BNmc_ft

        SVOC_lim1_ft = (1-beta_terp)*C1LHlmc_ft + beta_terp*C1LHhmc_ft + (1-beta_terp)*C1LOlmc_ft   + beta_terp*C1LOhmc_ft  + C1LNmc_ft 
        SVOC_lim2_ft = (1-beta_terp)*C2LHlmc_ft + beta_terp*C2LHhmc_ft + (1-beta_terp)*C2LOlmc_ft + beta_terp*C2LOhmc_ft + C2LNmc_ft 
        SVOC_lim3_ft = beta_terp*C3LHhmc_ft  + (1-beta_terp)*C3LOlmc_ft + beta_terp*C3LOhmc_ft + C3LNmc_ft
        SVOC_lim4_ft = beta_terp*C4LHhmc_ft +  (1-beta_terp)*C4LOlmc_ft + beta_terp*C4LOhmc_ft  + C4LNmc_ft
        
        SVOC_iso2_ft = (1-beta_iso)*C2IOlmc_ft + beta_iso*C2IOhmc_ft +  C2INmc_ft
        SVOC_iso3_ft = (1-beta_iso)*C3IOlmc_ft + beta_iso*C3IOhmc_ft + C3INmc_ft
        SVOC_iso4_ft = (1-beta_iso)*C4IOlmc_ft + beta_iso*C4IOhmc_ft + C4INmc_ft

        C1mc_ft    = (1-beta_terp)*C1AHlmc_ft  + beta_terp*C1AHhmc_ft  + (1-beta_terp)*C1BHlmc_ft &
        + beta_terp*C1BHhmc_ft  + (1-beta_terp)*C1LHlmc_ft  + beta_terp*C1LHhmc_ft  + (1-beta_terp)*C1AOlmc_ft &
        + beta_terp*C1AOhmc_ft  + (1-beta_terp)*C1BOlmc_ft  + beta_terp*C1BOhmc_ft  + (1-beta_terp)*C1LOlmc_ft &
        + beta_terp*C1LOhmc_ft  + C1ANmc_ft  + C1BNmc_ft  + C1LNmc_ft 
        C2mc_ft    = (1-beta_terp)*C2AHlmc_ft  + beta_terp*C2AHhmc_ft  + (1-beta_terp)*C2BHlmc_ft &
        + beta_terp*C2BHhmc_ft  + (1-beta_terp)*C2AOlmc_ft  + beta_terp*C2AOhmc_ft  + (1-beta_terp)*C2BOlmc_ft &
        + beta_terp*C2BOhmc_ft   &
        + beta_terp*C2LHhmc_ft  + (1-beta_iso)*C2IOlmc_ft  + beta_iso*C2IOhmc_ft  + (1-beta_terp)*C2LOlmc_ft &
        + beta_terp*C2LOhmc_ft  + C2ANmc_ft  + C2BNmc_ft  + C2INmc_ft  + C2LNmc_ft 
        C3mc_ft    = (1-beta_terp)*C3AHlmc_ft  + beta_terp*C3AHhmc_ft  + (1-beta_terp)*C3BHlmc_ft &
        + beta_terp*C3BHhmc_ft  + (1-beta_terp)*C3AOlmc_ft  + beta_terp*C3AOhmc_ft  + (1-beta_terp)*C3BOlmc_ft &
        + beta_terp*C3BOhmc_ft   &
        + beta_terp*C3LHhmc_ft  + (1-beta_iso)*C3IOlmc_ft  + beta_iso*C3IOhmc_ft  + (1-beta_terp)*C3LOlmc_ft &
        + beta_terp*C3LOhmc_ft  + C3ANmc_ft  + C3BNmc_ft  + C3INmc_ft  + C3LNmc_ft 
        C4mc_ft    = (1-beta_terp)*C4AHlmc_ft  + beta_terp*C4AHhmc_ft  + (1-beta_terp)*C4BHlmc_ft &
        + beta_terp*C4BHhmc_ft  + (1-beta_terp)*C4AOlmc_ft  + beta_terp*C4AOhmc_ft  + (1-beta_terp)*C4BOlmc_ft &
        + beta_terp*C4BOhmc_ft  &
        + beta_terp*C4LHhmc_ft  + (1-beta_iso)*C4IOlmc_ft  + beta_iso*C4IOhmc_ft  + (1-beta_terp)*C4LOlmc_ft &
        + beta_terp*C4LOhmc_ft  + C4ANmc_ft  + C4BNmc_ft  + C4INmc_ft  + C4LNmc_ft 
     
        ! T-dependence following Clausius-Clapeyron (Sheehan and Bowman '01)
        C1_st_ft    = C10_st * (T0/T_ft)*exp((dHvap/Rgas)*(1./T0-1./T_ft)) 
        C2_st_ft    = C20_st * (T0/T_ft)*exp((dHvap/Rgas)*(1./T0-1./T_ft))
        C3_st_ft    = C30_st * (T0/T_ft)*exp((dHvap/Rgas)*(1./T0-1./T_ft))
        C4_st_ft    = C40_st * (T0/T_ft)*exp((dHvap/Rgas)*(1./T0-1./T_ft))
        m        = 0 
        Coait_ft = OAbgb_ft  
        epsi_ft = 1.
        do while (epsi>1e-10) ! set iteration tolerance here
             X1_ft          = 1./((1.+(C1_st_ft/Coait_ft)))
             X2_ft          = 1./((1.+(C2_st_ft/Coait_ft)))
             X3_ft          = 1./((1.+(C3_st_ft/Coait_ft)))
             X4_ft          = 1./((1.+(C4_st_ft/Coait_ft)))
             Coait_new_ft   = X1_ft*C1mc_ft + X2_ft*C2mc_ft + X3_ft*C3mc_ft + X4_ft*C4mc_ft + OAbgb_ft                  
             epsi_ft         = abs(Coait_ft-Coait_new_ft) ! optimization parameter
             Coait_ft       = Coait_new_ft      
             m           = m+1  
             if (m>100) then
                write (*,*) 'Warning: no convergence (FT)!'  
             endif                   
        enddo ! while  
        Coa_ft = Coait_ft
    end if ! VBS    

!!!!!!!!   OUTPUT   !!!!!!!!

!   output

    if ((mod(tt,aver1)) == 0) then
!     if (tt == nint(60/dtime)) then
!     tt=0
      write (20,'(19F14.5)') thour,printhour,zi(2), &
        we,thetam(2),dtheta(2),wthetae,wthetas,beta, &
        um(2),du(2),vm(2),dv(2),ws,ac,wm,dz,wstar,dthetav

      write (50,'(2F14.4,13F14.4)') thour,printhour, &
        ustar,uws, vws, uwe, vwe, du(2), dv(2), &
        sqrt(du(2)**2+dv(2)**2)

      write (60,'(2F14.4,5F14.4,8E15.5)') thour,printhour, zi(2) &
        ,qm(2), dq(2), wqe, wqs, betaq &
        ,cm(2), dc(2), wce, wcs, betac, qs*1000., wqm*1000

      if(llandsurface)then
        write (61,'(12F14.5)') thour,printhour,Swin, &
          Swout,Lwin,Lwout,Qtot,SH,LE,GR,ra,rs
        if(lCO2Ags)then
          write (63,'(4F14.5)') thour,printhour,An, Resp
        endif
      endif

      if(lchem)then
        write(formatstring,'(A,i2,A)') '(2E15.5E3,',nchsp ,'E15.5E3)'

        write (40,formatstring) &
          thour,printhour,(c_cbl(k),k=1,nchsp)

        write (41,formatstring) &
          thour,printhour,(E(k),k=1,nchsp)

        write (42,formatstring) &
          thour,printhour,(beta_ft(k),k=1,nchsp)
        
        if(lvbs) then
              write (43,formatstring) &
                     thour,printhour, OAbgb, Coa, C1mc, C2mc, C3mc, C4mc, X1, &
                     X2, X3, X4, beta_terp, beta_iso  
              write (44,formatstring) &
                                   thour,printhour, OAbgb_ft, Coa_ft, C1mc_ft, &
                                   C2mc_ft, C3mc_ft, C4mc_ft, X1_ft, X2_ft, &
                                   X3_ft, X4_ft  
             write (49, formatstring) & 
                    thour, printhour,   alp_OH1, alp_OH2, alp_OH3, alp_OH4, alp_O31, alp_O32, alp_O33, alp_O34, C1ANmc, C2ANmc, C3ANmc, &
         C4ANmc, bep_OH1, bep_OH2, bep_OH3, bep_OH4, bep_O31, bep_O32, bep_O33, bep_O34, C1BNmc, C2BNmc, C3BNmc, C4BNmc, lim_OH1, lim_OH2, lim_OH3, &
         lim_OH4, lim_O31, lim_O32, lim_O33, lim_O34, C1LNmc, C2LNmc, C3LNmc, C4LNmc, isoo32, isoo33, isoo34, C2INmc, C3INmc, C4INmc
                    write (53, formatstring) & 
                    thour, printhour, C2LHlmc, C2LHhmc, C2LOlmc, C2LNmc, C4LHlmc, C4LHhmc, C4LOlmc, C4LNmc
        end if
        
        if(lBVOC)then
          write (45,'(15F14.5)') &
          thour,printhour, F_ISO, F_TERP, gamma_ISO, gammaT_ISO, gammap, gammace_ISO, gamma_MT, gammaT_MT, gammace_MT, gammaLAI, gammasm, PAR, Ts
        endif

        write (47,'(2F14.4,6(E14.5),f8.1)') &
        thour,printhour,6000. * RC(R_1%loc)%Keff_cbl, 60. * RC(R_NO2%loc)%Keff_cbl,&
               6000.*RC(R_CH2O%loc)%Keff_cbl,photo, zenith,Q_cbl(RH%loc), zenith*180/pi

        write (48,formatstring) &
          thour,printhour,(c_ft(k),k=1,nchsp)

        write (52, '(16F15.8)') &
          thour,printhour, gamma_uptake, q_aq, inv_gas, inv_Torq_aq, F_gIpx, k_org, omega, Rc_gIPX, T_hh, denom, q_org_new, HorgDorg, q_org, k_het

        write(formatstring,'(A,i3,A)') '(3F13.4,',tnor,'E13.4)'
        write (62,formatstring),thour,printhour,temp_cbl, (RC(k)%Keff_cbl,k=1,tnor)

        if(lwritepl) then
          do i=1,nchsp
            if((.not. lcomplex).or.(PL_scheme(i)%prin)) then
              write (formatstring,'(A,i3,A6)')'(f8.3,f8.3,',PL_scheme(i)%nr_PL +3,'E17.8)'
              write(100+i,formatstring),thour,printhour,c_cbl(i),(productionloss(i,k),k=1,PL_scheme(i)%nr_PL+2)
            endif  
          enddo
        endif
      endif !lchem

      
      
      if(lscu)then
         write (70,'(36e15.5)') thour,printhour,(wscm2(k),k=1,nchsp)
      endif
    endif !((mod(tt,aver1)) == 0)
        
    if (mod(tt,aver2) == 0.) then
      write (28,'(2F14.4,E14.5)') 0., thetam(2), wthetas
      write (28,'(2F14.4,E14.5)') zi(2), thetam(2),-beta*wthetas
      write (28,'(2F14.4,E14.5)') zi(2), thetam(2)+dtheta(2)+ eps, 0.+eps
      write (28,'(2F14.5,E14.5)') h_max, thetam(2)+dtheta(2)+ gamma*(h_max-zi(2)), 0.

      write (29,'(2F14.4,E14.5)') 0., qm(2), wqs
      write (29,'(2F14.4,E14.5)') zi(2), qm(2), wqe
      write (29,'(2F14.4,E14.5)') zi(2), qm(2)+dq(2)+ eps, 0.+eps
      write (29,'(2F14.4,E14.5)') h_max, qm(2)+dq(2)+ gammaq*(h_max-zi(2)), 0.

      write (30,'(2F14.4,E14.5)') 0., cm(2), wcs
      write (30,'(2F14.4,E14.5)') zi(2), cm(2), wce
      write (30,'(2F14.4,E14.5)') zi(2), cm(2)+dc(2)+ eps, 0.+eps
      write (30,'(2F14.4,E14.5)') h_max, cm(2)+dc(2)+ gammac*(h_max-zi(2)), 0.

      write (31,'(2F14.4,E14.5)') 0., um(2), uws
      write (31,'(2F14.4,E14.5)') zi(2), um(2), uwe
      write (31,'(2F14.4,E14.5)') zi(2), um(2)+du(2)+ eps, 0.+eps
      write (31,'(2F14.4,E14.5)') h_max, um(2)+du(2)+ gammau*(h_max-zi(2)), 0.

      write (32,'(2F14.4,E14.5)') 0., vm(2), vws
      write (32,'(2F14.4,E14.5)') zi(2), vm(2), vwe
      write (32,'(2F14.4,E14.5)') zi(2), vm(2)+dv(2)+ eps, 0.+eps
      write (32,'(2F14.4,E14.5)') h_max, vm(2)+dv(2)+ gammav*(h_max-zi(2)), 0.
    endif

    zi(1)=zi(2)

    dtheta(1)=dtheta(2)
    thetam(1)=thetam(2)

    um(1)=um(2)
    vm(1)=vm(2)
    du(1)=du(2)
    dv(1)=dv(2)

    dq(1)=dq(2)
    qm(1)=qm(2)

    dc(1)=dc(2)
    cm(1)=cm(2)

    if(.not. lfixedlapserates) then
      gamma  = gamma  * (1 + wsls * dtime)
      gamma2 = gamma2 * (1 + wsls * dtime)
      gammaq = gammaq * (1 + wsls * dtime)
      gammac = gammac * (1 + wsls * dtime)
      gammau = gammau * (1 + wsls * dtime)
      gammav = gammav * (1 + wsls * dtime)
    endif

  enddo !t=1, runtime

  write(28,'(a26,I5,a1)') 'Time int. between profiles', atime_vert,'s'
  write(29,'(a26,I5,a1)') 'Time int. between profiles', atime_vert,'s'
  write(30,'(a26,I5,a1)') 'Time int. between profiles', atime_vert,'s'
  write(31,'(a26,I5,a1)') 'Time int. between profiles', atime_vert,'s'
  write(32,'(a26,I5,a1)') 'Time int. between profiles', atime_vert,'s'

!  if (n == 0) then
!    write(20,'(a28)') 'Saturation level NOT reached'
!    write(60,'(a28)') 'Saturation level NOT reached'
!  else
!    write(20,'(a28)') 'Saturation level reached'
!    write(60,'(a28)') 'Saturation level reached'
!  end if

! print status
  close (20)
  close (28)
  close (29)
  close (30)
  close (31)
  close (32)
  close (40)
  close (41)
  close (42)
  if(lvbs) then
    close (43)
    close (44)
    close (49)
    close (53)
  endif
  if(lBVOC)then
    close (45)
  endif
  close (46)
  close (47)
  close (50)
  close (51)
  close (52)
  close (60)
  close (61)
  close (62)
  close (63)

end program

!c
!c ---- Function to calculate solar zenith angle for -65< lat <65 degrees north
!c------ and lon between -60 and + 60
real function getth(daynr,lat,lon,xhr)
implicit none
  real daynr
  real lat,lon
  real houra
  real obliq,deday,delta,lonr,latr
  real piby, xhr,pi
  pi   = acos(-1.)
  piby = pi/ 180.
  lonr = lon*piby
  latr = lat*piby
  obliq = 23.45 * piby
  deday = 4.88 + 2*pi/365* daynr
  delta = asin(sin(obliq)*sin(deday))
  houra = lonr - pi + xhr* (2.*pi/24.)
  getth = acos(sin(delta)*sin(latr) +  cos(delta)*cos(latr)*cos(houra))
  
  return
end

real function psim(zeta)
implicit none
  double precision zeta
  real x
  if(zeta <= 0) then
!   x    = (1. + 3.6 * abs(real(zeta)) ** (2./3.)) ** (-0.5)
!   psim = 3. * log( (1. + 1. / x) / 2.)
    x    = (1. - 16. * real(zeta)) ** (0.25)
    psim = 3.14159265 / 2. - 2. * atan(x) + log( (1.+x) ** 2. * (1. + x ** 2.) / 8.)
  else
    psim = -2./3. * (real(zeta) - 5./0.35) * exp(-0.35 * real(zeta)) - real(zeta) - (10./3.) / 0.35
  endif
  return
end

real function psih(zeta)
implicit none
  double precision zeta
  real x
  if(zeta <= 0) then
!   x    = (1. + 7.9 * abs(real(zeta)) ** (2./3.)) ** (-0.5)
!   psih = 3. * log( (1. + 1. / x) / 2.)
    x    = (1. - 16. * real(zeta)) ** (0.25)
    psih = 2. * log( (1. + x ** 2.) / 2. )
  else
    psih = -2./3. * (real(zeta) - 5./0.35) * exp(-0.35 * real(zeta)) - (1. + (2./3.) * real(zeta)) ** (1.5) - (10./3.) / 0.35 + 1.
  endif
  return
end

double precision function E1(x)
implicit none
  double precision x
  double precision E1sum, factorial
  integer k,t

  E1sum = 0.0
  do k=1,99
    E1sum = E1sum + (-1.0) ** (k + 0.0) * x ** (k + 0.0) / ( (k + 0.0) * factorial(k) )
  end do
  E1 = -0.57721566490153286060 - log(x) - E1sum
end

double precision function factorial(k)
implicit none
  integer k
  integer n
  factorial = 1.0
  do n = 2, k
    factorial = factorial * n
  enddo
end
