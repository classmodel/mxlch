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

  logical :: c_ustr=.true.,c_wth=.false.,c_fluxes=.false., c_jno2, lencroachment=.false., ladvecFT=.false. !c_fluxes replaces c_wth
  logical :: lradiation=.false.,lsurfacelayer=.false.,llandsurface=.false.,lrsAgs=.false., lCO2Ags=.false.
  double precision :: zi(2),zi0 = 200 ,thetam(2), dtheta(2),pressure = 1013.0, wthetae
  double precision temp_cbl, temp_ft
  integer :: runtime,t, time=24*3600.0,tt
  double precision :: beta = 0.2 ,wthetas=0.0,gamma = 0.006,thetam0 = 295,dtheta0 = 4,wthetav=0.0,dthetav
  real :: dtime = 1
  double precision :: z0 = 0.03, kappa, zp, alpha,z0m=0.03,z0h=0.03
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

  double precision :: wthetasmax = 0.0 , wqsmax = 0.0 , wcsmax = 0.0

  double precision :: uws=0,vws=0,uws0 = 0.0 , vws0 = 0.0
  double precision ustar
  double precision :: um(2), vm(2), um0 = 0.0 , vm0 = 0.0, ueff, wstar
  double precision :: ug = 0.0 , vg = 0.0
  double precision du(2), dv(2)
  double precision :: ws,wsls=0.0
  double precision uwe,vwe
  double precision :: f,gammau = 0.0 , gammav = 0.0
  double precision we
  double precision :: advq = 0.0 ,lsq   ! large scale advection moisture (units (g/Kg)/s)!
  double precision :: advtheta = 0.0 ,lstheta   ! large scale advection heat (units K/s)!
  double precision :: cc=0.0,Qtot=400.0,albedo=0.2 !cloud cover and incoming energy,albedo, incoming radiation
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

  double precision pi

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
  real,parameter :: rvrd = 0.61  !  Rd/Rv    
  real,parameter :: bolz =  5.67e-8 !Stefan-Boltzmann constant [-]
  real,parameter :: rhow =  1000.   !density of water [kg m-3]
  real,parameter :: rho  =    1.2   !density of air [kg m-3]
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
  character(len=40) filename
  integer j
  character(len=2),allocatable :: dummy(:)
  character(len=15) :: outdir = 'RUN00'

  !Adapted surface fluxes
  !Times in s after start of run
  integer  starttime_wT , endtime_wT , starttime_wq, endtime_wq
  integer  starttime_chem, endtime_chem, starttime_adv, endtime_adv
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
    wsls, &
    wthetasmax, &
    c_wth, &
    c_fluxes, &
    gamma, &
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
    lencroachment

  namelist/NAMSURFLAYER/ &
    lsurfacelayer,&
    z0m,&
    z0h

  namelist/NAMRAD/ &
    lradiation,& !radiation scheme to determine Q and SW
    cc,&   !cloud cover
    S0,&   !Incoming radiation
    albedo !Surface albedo

  namelist/NAMSURFACE/ &
    llandsurface,& !switch to use interactive landsurface
    Qtot,&         !Incoming energy
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
    Eact0          !activation energy [53.3 kJ kmol-1]

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
    starttime_adv , &
    endtime_adv   , &
    function_wT   , &
    function_wq   

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

  if ( lCO2Ags ) then
    lrsAgs       = .true.
    llandsurface = .true.
  endif
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

  write (*,*)'LCHEM=',lchem,'LDIUVAR=',ldiuvar,'LCONST=',lchconst,'LFLUX=',lflux
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
  starttime_adv  = daytime_start
  endtime_adv    = daytime_end
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

  open (20, file=trim(outdir)//dirsep//'output_dyn')
    write (20,'(a4)') 'TIME'
    write (20,'(I4)') time/atime
    write (20,'(14a14)') 'UTC(hours)','RT(hours)','zi(m)','we(m/s)', &
            'thetam(K)','dtheta(K)','wte(Km/s)','wts(Km/s)', &
            'beta','um(ms-1)','du(ms-1)','vm(ms-1)','dv(ms-1)','ws(ms-1)'

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
    write (60,'(13a14)') &
            'UTC(hours)','RT(h-ours)','zi(m)','qm(g/kg)','dq(g/kg)','wqe', &
            'wqs','   betaq','  cm(ppm)','dc(ppm)', 'wce','wcs','  betac'

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

    open ( 62, file =trim(outdir)//dirsep//'keff_cbl')
      write(formatstring,'(A,i3,A4)')'(',tnor+3,'A13)' !+3 id for 2 times and tem_cbl
      write (62,formatstring) 'UTC(hours)','RT(hours)','Temp_cbl',(RC(k)%rname,k=1,tnor)

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

      Swin  = S0 * Tr * costh
      Swout = albedo * Swin
      Lwin  = 0.8 * bolz * (Ta ** 4)
      Lwout = bolz * (Ts ** 4)

      Qtot  = Swin - Swout + Lwin - Lwout
    endif

    thetav  = thetam(1) * (1. + 0.61 * qm(1) * 1.e-3)
    wstar   = ( (g / thetav) * zi(1) * wthetav ) ** ( 1. / 3. )
    ueff    = sqrt(um(1)**2.+vm(1)**2.+wstar**2.)
!    ueff    = sqrt(um(1)**2.+vm(1)**2.)
    if (lsurfacelayer) then
      if(ueff .lt. 1.e-2) stop 'effective wind velocity below 1 cm/s'
      if(t==1) then
        thetasurf = Ts 
      else
        thetasurf = thetam(1) + wthetas / (Cs * ueff)
      endif
      esatsurf    = 0.611e3 * exp(17.2694 * (thetasurf - 273.16) / (thetasurf - 35.86))
      qsatsurf    = 0.622 * esatsurf / (pressure*100)
      cq          = 0.0
      if(t/=1) cq = (1. + Cs * ueff * rs) ** (-1.)
      qsurf       = (1. - cq) * qm(1) + cq * qsatsurf * 1.e3 !HGO factor for qsatsurf which is in kg/kg

      thetavsurf  = thetasurf * (1. + 0.61 * qsurf * 1.e-3)

      Rib         = min(0.2,g/thetav * zsl * (thetav-thetavsurf) / (ueff ** 2.))
      L           = sign(0.01,Rib)
      L0          = sign(0.1,Rib)

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
        f2     = min(1.e8, f2)

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
          betaw    = max(1.0e-3,min(1.0,(wg - wwilt)/(wfc - wwilt)))

          ! calculate stress function
          fstr     = betaw

          ! calculate gross assimilation rate (Am)
          Am       = Ammax * (1 - exp( -(gm * (ci - CO2comp) / Ammax) ) )

          Rdark    = (1.0/9) * Am

          PAR      = 0.40*max(0.1,Swin*cveg)

          ! calculate  light use efficiency
          alphac   = alpha0 * (co2abs - CO2comp) / (co2abs + 2 * CO2comp)

          ! 1.- calculate upscaling from leaf to canopy: net flow CO2 into the plant (An) 
          tempy    = alphac * Kx * PAR / (Am + Rdark)
          An       = (Am + Rdark) * (1 - 1.0 / (Kx * LAI) * (E1( tempy * exp(-Kx * LAI)) - E1(tempy)))


          ! 2.- calculate upscaling from leaf to canopy: CO2 conductance at canopy level
          AGSa1    = 1.0 / (1 - f0)
          Dstar    = D0 / (AGSa1 - 1)

          gcco2    = LAI * (gmin/nuco2q + AGSa1 * fstr * An / ((co2abs - CO2comp) * (1 + Ds / Dstar)))

          ! calculate surface resistance for moisture and carbon dioxide
          rsAgs    = 1.0 / (1.6 * gcco2)
          rsCO2    = 1.0 / gcco2

          rs       = rsAgs

          ! calculate net flux of CO2 into the plant (An)
          An       = -(co2abs - ci) / (ra + rsAgs)

          ! CO2 soil respiration surface flux
          fw       = Cw * wsmax / (wg + wsmin)
          Resp     = R10 * (1 - fw) * exp( Eact0 / (283.15 * 8.314) * (1. - 283.15 / (thetasurf)))

          wco2     = (An + Resp)*(MW_Air/MW_CO2)*(1.0/rho) !in ppm m/s

        endif

        ! recompute f2 using wg instead of w2
        f2     = 1.e8
        if (wg .gt. wwilt) then
          f2   = (wfc - wwilt) / (wg - wwilt)
        endif
        f2     = min(1.e8, f2)

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


        CG     = CGsat * (wsat / w2) ** (CLb / (2.0 * log(10.0)))

        Tsoiltend = CG * GR - 2.0 * pi / 86400.0 * (Tsoil - T2)
        Tsoil  = Tsoil + Tsoiltend * dtime

        C1     = C1sat * (wsat / wg) ** (CLb / 2.0 + 1.0)
        C2     = C2ref * (w2 / (wsat - w2) )
        wgeq   = w2 - wsat * CLa * ( (w2 / wsat) ** CLc * (1.0 - (w2 / wsat) ** (8.0 * CLc)) )

        wgtend = - C1 / (rhow * 0.1) * LEsoil / Lv - C2 / 86400 * (wg - wgeq)
        wg     = wg + wgtend * dtime

        wthetas= SH / (rho * Cp)
        wqs    = LE / (rho * Lv * 1.0e-3)
        wcs    = wcsmax

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
!   inroducing advection mpositure
    lsq     = 0.0
    lstheta = 0.0
    if ((sec .ge. starttime_adv) .and. (sec .lt. endtime_adv)) then
      lsq     = advq
      lstheta = advtheta
    endif

!   dynamics eulerien time-step
!   Potential temperature

    if (lencroachment .eqv. .false.) then
            
!   jump depending on the virtual temperature jump
      dthetav = dtheta(1) + &
      1.e-03*rvrd*(qm(1)*dtheta(1)+thetam(1)*dq(1)+dtheta(1)*dq(1))
      we = (beta*wthetav/dthetav) 
      
      zi(2)  =zi(1)+((beta*wthetav/dthetav)+ws)*dtime            !eq (8)
      dtheta(2)=dtheta(1)+ ((gamma*beta*wthetav/dthetav)- &
              (1/(zi(1)+inf))*(wthetas+we*dtheta(1)))*dtime           !eq. (3)
      if (.not. ladvecFT) dtheta(2) = dtheta(2) - lstheta*dtime
      thetam(2)=thetam(1)+ &
            ((1/(zi(1)+inf))*(wthetas+we*dtheta(1)))*dtime + &        !eq  (1)
            lstheta*dtime
    else
!     Calculation boundary layer height.
!     Based on Garratt (6.18) page 155 and
!     on the article of Betts 1973 and 1974
      beta   = 0.
      we     = 0. 

      zi(2)   = zi(1) + ((1/(gamma*(zi(1)+inf)))*((1+beta)*wthetav) &
              + ws)*dtime
      dtheta(2)=0.
      thetam(2)=thetam(1)+ &
            ((1/(zi(1)+inf))*wthetas)*dtime + lstheta*dtime         !eq  (1)
    endif

!    thetam(2)=thetam(1)+ &
!            ((1/(zi(1)+inf))*(wthetas+we*dtheta(1)))*dtime          !eq  (1)

!    we = (zi(2)-zi(1))/dtime - ws

! ---- Specific humidity
    if (beta /= 0) then
        wqe = -we*dq(1)
        dq(2)=dq(1)+ ((gammaq*we)-(1/(zi(1)+inf))*(wqs - wqe))*dtime   !eq. (3)
        if (.not. ladvecFT) dq(2) = dq(2) - lsq*dtime
    else
        wqe  = 0.
        dq(2)= 0.
    endif

    qm(2)=qm(1)+((1/(zi(1)+inf))*(wqs - wqe))*dtime+lsq*dtime          !eq  (1)
! -----------Carbon Dioxide

    if (beta /= 0) then
      wce = -we*dc(1)
      dc(2)=dc(1)+ ((gammac*we)-(1/(zi(1)+inf))*(wcs - wce))* dtime   !eq. (3)
    else
      wce    = 0.
      dc(2)  =0.
    endif

    cm(2)=cm(1)+((1/(zi(1)+inf))*(wcs - wce))*dtime        !eq  (1)

! ---- Windspeed

    if (beta /= 0) then
      uwe = -we*du(1)
      vwe = -we*dv(1)
      du(2)=du(1)+ ((gammau*we)-(1/(zi(1)+inf))*(uws - uwe)+ f*dv(1))* dtime                 !eq. (3)
      dv(2)=dv(1)+ ((gammav*we)-(1/(zi(1)+inf))*(vws - vwe) - f*du(1))*dtime                 !eq. (3)
    else
      uwe = 0.
      vwe = 0.
      du(2)=0.
      dv(2)=0.
    endif

    um(2)=um(1)+(-f*(dv(1))+(1/(zi(1)+inf))*(uws-uwe))*dtime
    vm(2)=vm(1)+(f*(du(1))+(1/(zi(1)+inf))*(vws-vwe))*dtime

!   closure assumption
    wthetae= - beta *wthetav

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

!    emission and entrainment each timestep

    if(lchem) then
      do k=1, nchsp
        E(k)=-we*(c_ft(k)-c_cbl(k))
        c_cbl(k)=c_cbl(k)+(1/zi(1))*(Q_cbl(k)-E(k))*dtime
        c_cbl(k)=c_cbl(k)+adv_chem_cbl(k)*dtime
        c_ft( k)=c_ft( k)+adv_chem_ft(k )*dtime
        c_cbl(k) = max(0.0, c_cbl(k) )
        c_ft( k) = max(0.0, c_ft( k) )

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

!   output

    if ((mod(tt,aver1)) == 0) then
!     if (tt == nint(60/dtime)) then
!     tt=0
      write (20,'(14F14.5)') thour,printhour,zi(2), &
        we,thetam(2),dtheta(2),wthetae,wthetas,beta, &
        um(2),du(2),vm(2),dv(2),ws

      write (50,'(2F14.4,13F14.4)') thour,printhour, &
        ustar,uws, vws, uwe, vwe, du(2), dv(2), &
        sqrt(du(2)**2+dv(2)**2)

      write (60,'(2F14.4,5F14.4,6E15.5)') thour,printhour, zi(2) &
        ,qm(2), dq(2), wqe, wqs, betaq &
        ,cm(2), dc(2), wce, wcs, betac

      if(lchem)then
        write(formatstring,'(A,i2,A)') '(2E15.5E3,',nchsp ,'E15.5E3)'

        write (40,formatstring) &
          thour,printhour,(c_cbl(k),k=1,nchsp)

        write (41,formatstring) &
          thour,printhour,(E(k),k=1,nchsp)

        write (42,formatstring) &
          thour,printhour,(beta_ft(k),k=1,nchsp)

        write (47,'(2F14.4,6(E14.5),f8.1)') &
        thour,printhour,6000. * RC(R_1%loc)%Keff_cbl, 60. * RC(R_NO2%loc)%Keff_cbl,&
               6000.*RC(R_CH2O%loc)%Keff_cbl,photo, zenith,Q_cbl(RH%loc), zenith*180/pi

        write (48,formatstring) &
          thour,printhour,(c_ft(k),k=1,nchsp)

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

  enddo !t=1, runtime

  write(28,'(a26,I5,a1)') 'Time int. between profiles', atime_vert,'s'
  write(29,'(a26,I5,a1)') 'Time int. between profiles', atime_vert,'s'
  write(30,'(a26,I5,a1)') 'Time int. between profiles', atime_vert,'s'
  write(31,'(a26,I5,a1)') 'Time int. between profiles', atime_vert,'s'
  write(32,'(a26,I5,a1)') 'Time int. between profiles', atime_vert,'s'

  if (n == 0) then
    write(20,'(a28)') 'Saturation level NOT reached'
    write(60,'(a28)') 'Saturation level NOT reached'
  else
    write(20,'(a28)') 'Saturation level reached'
    write(60,'(a28)') 'Saturation level reached'
  end if

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
  close (46)
  close (47)
  close (50)
  close (51)
  close (60)
  close (62)

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
