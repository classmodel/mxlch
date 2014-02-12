Module modchem
! version 1.0 3 chemicals on input and output but only integer coefficients.
! version 1.1 4 chemicals on input and 4 on output integers on input but real coefficients on output side allowed
! version 1.2 4 chem on input and NCCAA(zie modchem_mxl) on output

implicit none
save

!c***********************************************************************
!c
!c  parameters for a set of 1st/2nd order chemical reactions
!c
!c***********************************************************************

!cACP - BEGIN: Constants related to chemistry !min(nlimit)=1!
integer    nchsp,tnor  !nchsp = NumberCHemiscalSPecies   tnor=TotalNumberOfReactions
real    :: latt = 0, long = 0, tday,thour,zenith  !used for calculation solar zenith angle

integer mrpcc	!max number of reactions in which a particular chemical is used size is not important
				!but should be large enough
parameter (mrpcc = 99)
integer,parameter :: NCCBA = 4   !Number Chemical Components Before Arrow !!!!!!!  4 is MAXIMUM !!!!!
integer,parameter :: NCCAA = 10   !Number Chemical Components After Arrow
integer,parameter :: NNSPEC = 2*(NCCBA + NCCAA) - 1

integer CBL,FT
parameter( CBL = 1, FT = 0)
integer PRODUCTION, LOSS
parameter (PRODUCTION=1, LOSS=2)
logical :: lday = .false.
logical :: dayswitch = .false.
real, parameter :: MW_Air = 28.97
real, parameter :: MW_H2O = 18
real, parameter :: MW_CO2 = 44


logical :: lchem = .false.
logical :: lwritepl = .true.
logical :: lcomplex = .false.

logical :: ldiuvar =.true.
real    :: h_ref   = 12.

logical :: lflux   = .false.
real    :: fluxstart = 0.
real    :: fluxend   = 0.

logical :: lchconst  = .false.
real :: t_ref_cbl    = 298.
real :: p_ref_cbl    = 1013.5
real :: q_ref_cbl    = 10.
real :: t_ref_ft     = 298.
real :: p_ref_ft     = 1013.5
real :: q_ref_ft     = 10.
real :: pressure_ft
real :: convcbl, convft

type RCdef
	character*6 rname
	integer*1 raddep	       !1 if reaction = radiation dependend
	double precision Kreact	   !reaction konstant from input file
	double precision Keff_cbl  !to use with special circumstances cq with radiation and/or temperature depend reactions
	double precision Keff_ft   !can be different from Keff_cbl different Temp or moisture cq cloud
	integer order;
	integer func1
	real A
	real B
	real C
	real D
	real E
	real F
	real G
end type RCdef

type (RCdef), allocatable :: RC(:) !tnor

type Form
	integer*1 formula	!number of the formula to use
	integer r_nr		!reaction number, index to RC
	integer PorL		!0=> not active   1=>Production  2=>Loss
	real coef		    !coefficient in formula
	integer comp1		!index to c_cbl
	integer*1 exp1
	integer comp2		!index to c_cbl
	integer*1 exp2
	integer comp3		!index to c_cbl
	integer*1 exp3
	integer comp4		!index to c_cbl
	integer*1 exp4
end type Form

type Name_Number
	character (len=8) name	!name of chemical
	logical active			!active=1 else 0
	logical prin            !active=1 else 0
	integer chem_number		!number (not really used)
	integer nr_PL			!total number of reactions in which this chemical is used
	type (Form) PL(mrpcc)	!stucture holding the reaction components, reaction number etc
end type Name_Number

type (Name_Number), allocatable ::PL_scheme(:)   !(nchsp)

type Chem
	real coeff
	character (len=8) name
	integer chem_nr
	integer index_sv0
end type Chem

type Reaction
	character*6 name
	real kr			!kn2rd
	integer RadDep   !reaction is radiation dependend
	integer Order	!order of reaction
	integer nr_chem	!nr of chemicals in reaction (including non active species
	integer nr_chem_inp !nr of chem on input
	integer nr_chem_outp !nr of chem on output
	type (Chem) inp(NCCBA)
	type (Chem) outp(NCCAA)
end type Reaction

type (Reaction),allocatable ::reactions(:)

type (Name_Number), allocatable :: PL_temp(:)
integer, allocatable ::plot(:)
double precision, allocatable :: c_cbl(:),c_ft(:)  ! 0:nchsp to create dummy space
double precision, allocatable :: adv_chem_cbl(:),adv_chem_ft(:)  ! 0:nchsp to create dummy space
double precision, allocatable :: beta_ft(:),beta_ft_p(:)
double precision, allocatable :: Q_cbl(:), Q_init(:),E(:)
double precision, allocatable :: c_current(:)
integer, allocatable :: Q_func(:)

real ,allocatable ::productionloss(:,:)

type location
	character (len = 6) name
	integer ::  loc = 0
end type location

type (location) :: INERT, PRODUC , O3, O1D, NO2, NO, NO3, N2O5, HNO3, R, ISO, RO2, H2O2, HO2, HO, CO, CO2, H2O, NH3, H2SO4, CH2O, CH3O2, MVK, TERP, OAbg, CiT, CiI
type (location) :: R_O3,  R_NO,  R_NO2,  R_RH,  R_CO,  R_HNO3,  R_23,  R_25,  R_26a,  R_28,  R_43,  R_45,  R_54A,  R_57A,  R_58A,  R_61A,  R_62A,  R_63Aa,  R_63Ab, R_CH2O, R_1 , R_ISORO2NO, R_ISORO2HO2   
end module
