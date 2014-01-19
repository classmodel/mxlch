! version 1.0 3 chemicals on input and 4 on output but only integer coefficients.
! version 1.1 4 chemicals on input and 4 on output integers on input but real coefficients on output side allowed
!
! Adjust reaction speeds
!
subroutine calc_K_simple(pressure,temp_cbl, temp_ft)
use modchem
implicit none
double precision temp_cbl, temp_ft
double precision pressure
integer i
real coszen, factor
real getth
real K,K1,K2,K3,k4,k5
real conv_cbl,conv_ft
real Rfact


Rfact= 8.314e-2 ! mbar*m3 /K*mol

if(lchconst) then
   conv_cbl=6.023e8 * p_ref_cbl / (Rfact*temp_cbl)  ! =1 /6.77e-16/1000/60 = (6.023e23/e15)*P/RT
   conv_ft =6.023e8 * p_ref_ft / (Rfact*temp_ft)  ! =1 /6.77e-16/1000/60 = (6.023e23/e15)*P/RT
else
   conv_cbl=6.023e8 * pressure    / (Rfact*temp_cbl)  ! =1 /6.77e-16/1000/60 = (6.023e23/e15)*P/RT
   conv_ft =6.023e8 * pressure_ft / (Rfact*temp_ft)  ! =1 /6.77e-16/1000/60 = (6.023e23/e15)*P/RT
endif

convcbl = conv_cbl
convft  = conv_ft

!c
!c   Calculation solar zenith angle according to LES
!c

zenith    = getth(tday,latt,long,thour)
coszen = max(0.0, cos(zenith))

if (coszen > 0.00) then !it is day
  lday = .true.
  if (dayswitch .eqv. .false.) then
    dayswitch = .true.
    write(*,*)'The SUN is UP at time=',thour
  endif
  if (ldiuvar .eqv. .false.) then ! we have to fix the sza to the fixed hour h_ref
    zenith = getth(tday,latt,long,h_ref)
    coszen = max(0.0,cos(zenith))
  endif
else
  lday = .false.
  if(dayswitch .eqv. .true.) then
    dayswitch = .false.
    write(*,*)'The SUN is DOWN at time=',thour
  endif
endif

! adjust the Kreact depending on the func1 code

 do i=1, tnor
  if (RC(i)%RadDep == 1 ) then
    if (lday .eqv. .false. ) then
      RC(i)%Keff_cbl = 0.0
      RC(i)%Keff_ft  = 0.0
    else
      select case ( RC(i)%func1 )
      case (0) ! constant independent of sza
        RC(i)%Keff_cbl = RC(i)%KReact
        RC(i)%Keff_ft  = RC(i)%KReact
      case (1) ! constant independent of sza
        RC(i)%Keff_cbl = RC(i)%A
        RC(i)%Keff_ft  = RC(i)%A
      case (2)! exponential function
        RC(i)%Keff_cbl = RC(i)%A * exp(RC(i)%B /coszen )
        RC(i)%Keff_ft  = RC(i)%A * exp(RC(i)%B /coszen )
      case (3) ! powerfunction
        RC(i)%Keff_cbl = RC(i)%A * coszen ** RC(i)%B
        RC(i)%Keff_ft  = RC(i)%A * coszen ** RC(i)%B
      case (4) ! powerfunction but special for JO3
      !need [H20] in kg/kg so c_cbl[H2O]*1e-9
         K = RC(i)%A * coszen ** RC(i)%B
         RC(i)%Keff_cbl = K * RC(i)%D *  c_cbl(H2O%loc)*1.e-9 / &
              (RC(i)%D * c_cbl(H2O%loc)*1.e-9  + RC(i)%E * (1.- c_cbl(H2O%loc)*1.e-9))
         RC(i)%Keff_ft = K * RC(i)%D *  c_ft(H2O%loc)* 1.e-9 / &
              (RC(i)%D * c_cbl(H2O%loc)* 1.e-9 + RC(i)%E * (1.- c_cbl(H2O%loc)*1.e-9))
       case default !if someone put by mistake a number
         RC(i)%Keff_cbl = 1
         RC(i)%Keff_ft  = 1
      end select
     endif
  else   ! normal (no photolysis) reactions func1 can has a different meaning
    select case ( RC(i)%func1 )
    case(0)
      !do nothing K is in PPB and constant
    case(1) ! K in cm3/molecules*sec
      RC(i)%Keff_cbl = RC(i)%A * conv_cbl
      RC(i)%Keff_ft  = RC(i)%A * conv_ft
    case(2) !temperature depence of K for both cbl and ft
      RC(i)%Keff_cbl = RC(i)%A * exp(RC(i)%B / temp_cbl)* conv_cbl
      RC(i)%Keff_ft  = RC(i)%A * exp(RC(i)%B / temp_ft )* conv_ft
    case (3) !more complex temperature dependence
      RC(i)%Keff_cbl = RC(i)%A * (temp_cbl/RC(i)%B)**RC(i)%C * exp(RC(i)%D / temp_cbl)* conv_cbl
      RC(i)%Keff_ft  = RC(i)%A * (temp_ft /RC(i)%B)**RC(i)%C * exp(RC(i)%D / temp_ft )* conv_ft
    case(4:5) !use Lindemann-Hinshelwood equation  4 = second order 5 = first order so no conv_XXX factor
      !first CBL
      k1=RC(i)%A * (temp_cbl/300)**RC(i)%B * exp(RC(i)%C/temp_cbl) * conv_cbl * 1e9
      k2=RC(i)%D * (temp_cbl/300)**RC(i)%E * exp(RC(i)%F/temp_cbl)
      K = k1*k2/(k1+k2) * RC(i)%G
      if (RC(i)%func1 == 4) then
          RC(i)%Keff_cbl = K * conv_cbl
      else
          RC(i)%Keff_cbl = K
      endif
      !for FT
      k1=RC(i)%A * (temp_ft/300)**RC(i)%B * exp(RC(i)%C/temp_ft) * conv_ft * 1e9
      k2=RC(i)%D * (temp_ft/300)**RC(i)%E * exp(RC(i)%F/temp_ft)
      K = k1*k2/(k1+k2) * RC(i)%G
      if (RC(i)%func1 == 4) then
          RC(i)%Keff_ft = K * conv_ft
      else
          RC(i)%Keff_ft = K
      endif
    case(6)!special function of reaction 2H02->H202
       !first CBL
       k1 = RC(i)%A * exp(RC(i)%B / temp_cbl)* conv_cbl
       k2 = RC(i)%C * exp(RC(i)%D / temp_cbl)* conv_cbl**2 *1e9
       k3 = RC(i)%E * exp(RC(i)%F / temp_cbl)* conv_cbl * c_cbl(H2O%loc)
       RC(i)%Keff_cbl = (k1+k2)*(1+k3)
       !for FT
       k1 = RC(i)%A * exp(RC(i)%B / temp_ft)* conv_ft
       k2 = RC(i)%C * exp(RC(i)%D / temp_ft)* conv_ft**2 *1e9
       k3 = RC(i)%E * exp(RC(i)%F / temp_ft)* conv_ft * c_ft(H2O%loc)
       RC(i)%Keff_ft = (k1+k2)*(1+k3)
     case(7) ! same as 3 but third order so conv_XXX**2
       RC(i)%Keff_cbl = RC(i)%A * (temp_cbl/RC(i)%B)**RC(i)%C * exp(RC(i)%D / temp_cbl)* (conv_cbl**2)
       RC(i)%Keff_ft  = RC(i)%A * (temp_ft /RC(i)%B)**RC(i)%C * exp(RC(i)%D / temp_ft )* (conv_ft**2)
     case default !if someone put by mistake a different number
       write (*,*) 'FUNCTIONS GREATER THEN 7 NOT IMPLEMENTED'
       STOP
    end select
  endif
end do !tnor

end subroutine

subroutine iter_simple(switch,ynew,current,dt)
use modchem
implicit none

! -----------------------------------------------------
! This is the CHEMISTRY module of the MXL model
! It calulates the chemistry according to
! Krol et al. (2003) and Vinuesa et al. (2006)
! It also calculates the diurnal variation
! of photolysis rate (now in calc_K )

! Definition variables

  integer switch
  integer t !t is the number of sec since the beginning of the run (for dtime=1)
  double precision, dimension (0:nchsp) :: ynew, current
  real dt

  integer n,i,j,k
  integer iiter,niter
  double precision,target :: YP, YL
  double precision, POINTER :: YPL_pointer
  double precision kreact
  double precision YPL

! start Gauss Seidel iterations
! Gauss-Seidel iterations

  niter=4
  do iiter=1,niter
    do n=1,nchsp
     if (PL_scheme(n)%active .EQV. .TRUE.)then
!       if (PL_scheme(n)%name == CO%name) cycle  ! don't do calculations for CO
       if (PL_scheme(n)%name == H2O%name) cycle ! don't do calculations for H2O
       if (PL_scheme(n)%name == PRODUC%name) cycle
   
       YL = 0.
       YP = 0.
       do j = 1, PL_scheme(n)%nr_PL
   
         if (PL_scheme(n)%PL(j)%PorL == PRODUCTION ) then  ! this the production
           YPL_pointer => YP
         else
           YPL_pointer => YL
         endif

         if ( switch == CBL ) then
            kreact = PL_scheme(n)%PL(j)%coef * RC(PL_scheme(n)%PL(j)%r_nr)%Keff_cbl
         else
            kreact = PL_scheme(n)%PL(j)%coef * RC(PL_scheme(n)%PL(j)%r_nr)%Keff_ft
         endif
   
         select case (PL_scheme(n)%PL(j)%formula)
           case (0)
             YPL = kreact
           case (1)
             YPL = kreact * ynew(PL_scheme(n)%PL(j)%comp1)
           case (2)
             YPL = kreact * ynew(PL_scheme(n)%PL(j)%comp1) * ynew(PL_scheme(n)%PL(j)%comp2)
           case (3)
             YPL = kreact * (ynew(PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1)
           case (4)
             YPL = kreact * (ynew(PL_scheme(n)%PL(j)%comp1) ** PL_scheme(n)%PL(j)%exp1) &
                          * (ynew(PL_scheme(n)%PL(j)%comp2) ** PL_scheme(n)%PL(j)%exp2)
           case (5)
             YPL = kreact * ynew(PL_scheme(n)%PL(j)%comp1) * ynew(PL_scheme(n)%PL(j)%comp2) *ynew(PL_scheme(n)%PL(j)%comp3)
           case (6)
             YPL = kreact * (ynew(PL_scheme(n)%PL(j)%comp1)** PL_scheme(n)%PL(j)%exp1) &
                          * (ynew(PL_scheme(n)%PL(j)%comp2)** PL_scheme(n)%PL(j)%exp2) &
                          * (ynew(PL_scheme(n)%PL(j)%comp3)** PL_scheme(n)%PL(j)%exp3) 
           case (7)
             YPL = kreact * (ynew(PL_scheme(n)%PL(j)%comp1)** PL_scheme(n)%PL(j)%exp1) &
                          * (ynew(PL_scheme(n)%PL(j)%comp2)** PL_scheme(n)%PL(j)%exp2) &
                          * (ynew(PL_scheme(n)%PL(j)%comp3)** PL_scheme(n)%PL(j)%exp3) &
                          * (ynew(PL_scheme(n)%PL(j)%comp4)** PL_scheme(n)%PL(j)%exp4) 
         end select
   
         YPL_Pointer = YPL_Pointer + YPL
         if(lwritepl) then
          if( switch == CBL )then
            if (iiter == 1 .and. n == 1) then
              productionloss = 0.
            endif
            if (iiter == niter) then
              if (PL_scheme(n)%PL(j)%PorL == PRODUCTION ) then
                productionloss(n,j) = YPL
              else
                productionloss(n,j) = YPL * ynew(n)
              endif
            endif
          endif
         endif
       enddo !j=1, PL_scheme(n)%nr_PL
   
       ynew(n) = max(0.0, (current(n)+dt*YP)/(1.0+dt*YL))
       
       if (lwritepl) then
         if (switch == CBL)then
           productionloss(n, PL_scheme(n)%nr_PL + 1)= YL*ynew(n)
           productionloss(n, PL_scheme(n)%nr_PL + 2)= YP
         endif
       endif
   
     endif  !active == true
   enddo !i=1,nchsp
  enddo !iiter
  return
end
