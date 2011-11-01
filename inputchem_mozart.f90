subroutine inputchem_mozart(inputchemfile,outdir,dirsep)
! version 1.0 3 chemicals on input and 4 on output but only integer coefficients.
! version 1.1 4 chemicals on input and NCAA on output integers on input but real coefficients on output side allowed
!  
! special version for MOZART
!

use modchem
implicit none

  character*25 inputchemfile
  character*15 outdir
  character*1  dirsep

  integer i,j,k,l,react,max_reactions
  integer*2 number,nr_raddep_react
  integer react_nr
  real coefficient
  real reactconst
  real fact(7),num
  integer raddep,func1,nr_chemcomp,nr_active_chemicals
  integer count_raddep
  character*12 spec(25)
  character*6 rname
  character*255 line,line1
  character*255 scalarline
  character*8, allocatable::chem_name(:)
  logical prod,found
  character (len=8) tempname
  character (len=8) name
  character (len=3) react_str
  character (len=4) coef_str
!  character (len=3) comp1_str
!  character (len=3) comp2_str
  integer icoeff(4)

  react = 0
  number = 0
  count_raddep = 0

  open (unit=10,file='chemicals.txt',err=100,status='old',form='formatted')

  do while(.true.)
    read(10,'(a)',err=100) line
    if (line(1:1)=='#')then
      print *, line
    elseif (line(1:1)=='%')then
      read(line(2:10),*)nchsp        !read number of chemicals
      call allocate_arrays1()         ! allocate space
      allocate(chem_name(nchsp))
    elseif (len(trim(line))== 0)then
       !empty line do nothing
    elseif (line(1:1) == '$') then
      !end of chemicals
      exit  !do while
    else
      number = number + 1
      read(line,*)chem_name(number),c_cbl(number),c_ft(number),Q_init(number),Q_func(number),plot(number)
    endif  
  end do   
  close(10)
  Q_cbl= Q_init
  call read_chem_mozart(chem_name)
  
  
  write(*,*)'Number of chemicals = ',number
  nchsp=number !the real number of chemicals 
  open(unit=11,file='mozart.log',status='unknown')
  open (unit=10,file=inputchemfile,err=100,status='old',form='formatted')
 
  do while(.true.)
    read(10,'(a)',err=100) line
    if (line(1:1)=='#')then
      print *, line
    elseif (line(1:1)=='%')then
      read(line(2:10),*)tnor        !read number of reactions
      call allocate_arrays2()       ! allocate space for reactions
    elseif (len(trim(line))== 0)then
       !empty line do nothing
    elseif (line(1:1)=='$')then
       !end of chemical rections
       exit   !ends do while
    else
      react = react + 1
      if( react > tnor ) then
          write(6,*) 'Number of reactions is greater then specified',tnor
          STOP
      endif

      !    print *,'analysing'
      write(*,'(i03,2x,a)'),react,line
      write(11,'(i03,2x,a)'),react,line
       do i=1,25
        spec(i)='            '
      enddo
      read(line,*,end=200)(spec(j),j=1, 25)
  200 j=j-1
      nr_chemcomp = (j+1)/2
          ! the above 2 statements only work correctly with intel fortran
      i=1
      do while (len_trim(spec(i))>0 .and. i<25)
        ! print *,i, spec(i),len_trim(spec(i))
        i=i+1
      end do

      !print *,'number chemical components=', i,(i+1)/2
      nr_chemcomp = (i + 1)/2
      prod = .false.
      L=0
      
      fact = 0.
       
      read(10,'(a)',err=100) line1
      read(line1,*,end=300)rname,raddep,func1,(fact(j),j=1,7)
  300 i=0
      do i=j,7
       fact(i)=0.
      enddo
      
      reactions(react)%kr     = 0
      reactions(react)%name   = rname
      reactions(react)%RadDep = raddep

      RC(react)%Kreact   = 0
      RC(react)%Keff_cbl = 0
      RC(react)%Keff_ft  = 0
      RC(react)%rname   = rname
      RC(react)%RadDep = raddep
      RC(react)%func1 = func1
      RC(react)%order = 0;
      RC(react)%A = fact(1)
      RC(react)%B = fact(2)
      RC(react)%C = fact(3)
      RC(react)%D = fact(4)
      RC(react)%E = fact(5)
      RC(react)%F = fact(6)
      RC(react)%G = fact(7)

      if (raddep == 1) then
        count_raddep = count_raddep +1
      end if

      reactions(react)%nr_chem = nr_chemcomp

 
  !analyze reaction scheme,determine chem species and location in scalar.inp
  !and store in reactions(react)

      do j = 1, 2*nr_chemcomp - 1
        select case (spec(j))
          case ('!          ')
             ! stop searching for more comment starts
             exit  ! jump out
          case ('+          ')
             !print *,'found +'
          case ('->         ')
            prod = .true.
            reactions(react)%nr_chem_inp = L
            RC(react)%order = L - 1 !if 2 species we have to multiply conv if 3 species conv^2 with 4 species conv^3
                                    !only fro calc_K functions 1, 2, 3
            L=0
              !print *,'found ->'
          case default
            !print *, j,spec(j)
            L=L+1
            if ( spec(j)(1:1) == '(' ) then
              !non active species forget it
              L=L-1
            else
              do i=1,len(spec(j))
                if( spec(j)(i:i) .GT. '@' ) then
                  !starting name chemical component
                  tempname = spec(j)(i:len(spec(j)))
                  if( i == 1 ) then !nothing before chem comp
                    coefficient = 1.
                  else !we have numbers before
                    read( spec(j)(1:i-1),*)coefficient
                  endif
                  if ((prod .eqv. .false.) .and. ((coefficient +.0005)< 1.))then
                    write(*,*) 'Sorry, coefficient on input should be a multiply of 1'
                    STOP
                  endif
                  if ((prod .eqv. .false.) .and. (coefficient/int(coefficient) > 1.005) ) then
                    write(*,*) 'Sorry, coefficient on input should be a multiply of 1 found:',spec(j)
                    STOP
                  endif
                  exit
                else
                  if (i >= len(spec(j)) )then
                    write(*,*)'Probably space between coefficient and chemical component'
                    write(*,*) 'look between',spec(j),spec(j+1)
                    STOP
                  endif  
                endif
              enddo

              !find index in sv0
              i=1
              do while(tempname /= chem_name(i) )
                i= i+1
                if (i > nchsp) then
                  print *,'Name ',tempname, 'NOT FOUND in speciesline after @'
                  stop
                end if
              end do

              if (prod .EQV. .false.) then
                reactions(react)%inp(L)%name = tempname
                reactions(react)%inp(L)%coeff = coefficient
                reactions(react)%inp(L)%index_sv0 = i
              else
                reactions(react)%outp(L)%name = tempname
                reactions(react)%outp(L)%coeff = coefficient
                reactions(react)%outp(L)%index_sv0 = i
              endif
            endif
        end select
      enddo ! 1, 2*nr_chemcomp - 1
      reactions(react)%nr_chem_outp = L
    endif


  enddo ! end while(1)

  close(10)
  close(11)
  print *, 'Total number of reactions is',react,'of which', count_raddep,'is/are radiation dependent'

  !we now make tnor equal to real number of reactions it could have been to large
  tnor = react

  !make a list of chemical species and in which reaction number it is formed and/or losst
  k=0
  do i=1,react
    do j=1,reactions(i)%nr_chem_inp ! look only on input side of reaction
      name = reactions(i)%inp(j)%name
      found = .false.
      do L=1,k
        if(name == PL_scheme(L)%name) then
          found = .true.
          reactions(i)%inp(j)%chem_nr = L    !put chem component number in reaction
          PL_scheme(L)%nr_PL = PL_scheme(L)%nr_PL +1 !count number of reactions
          if ( PL_scheme(L)%nr_PL > mrpcc ) then
            print *, 'mrpcc to low, increase mrpcc in modchem'
            stop
          end if
          PL_scheme(L)%PL(PL_scheme(L)%nr_PL)%r_nr = i   !store reaction number index to RC
          PL_scheme(L)%PL(PL_scheme(L)%nr_PL)%PorL = LOSS   !this is a loss reaction for this component
          exit
        end if
      enddo
      if (found .EQV. .false.) then
        k=k+1
        PL_scheme(k)%name=name
        PL_scheme(k)%chem_number = k
        reactions(i)%inp(j)%chem_nr = k
        PL_scheme(L)%nr_PL = PL_scheme(L)%nr_PL +1
        if ( PL_scheme(L)%nr_PL > mrpcc ) then
          print *, 'mrpcc to low, increase mrpcc in modchem'
          stop
        end if
        PL_scheme(L)%PL(PL_scheme(L)%nr_PL)%r_nr = i   !store reaction number
        PL_scheme(L)%PL(PL_scheme(L)%nr_PL)%PorL = LOSS   !this is a loss reaction for this component
      endif
    enddo
  enddo

  do i=1,react
    do j=1,reactions(i)%nr_chem_outp  !this is for after the ->
      name = reactions(i)%outp(j)%name
      found = .false.
      do L=1,k
        if(name == PL_scheme(L)%name) then
          found = .true.
          reactions(i)%outp(j)%chem_nr = L
          PL_scheme(L)%nr_PL = PL_scheme(L)%nr_PL +1
          if ( PL_scheme(L)%nr_PL > mrpcc ) then
            print *, 'mrpcc to low, increase mrpcc in modchem'
            stop
          end if
          PL_scheme(L)%PL(PL_scheme(L)%nr_PL)%r_nr = i   !store reaction number
          PL_scheme(L)%PL(PL_scheme(L)%nr_PL)%PorL = PRODUCTION   !this is a production reaction for this component
          exit
        end if
      enddo
      if (found .EQV. .false.) then
        k=k+1
        PL_scheme(k)%name=name
        PL_scheme(k)%chem_number = k
        reactions(i)%outp(j)%chem_nr = k
        PL_scheme(L)%nr_PL = PL_scheme(L)%nr_PL +1
        if ( PL_scheme(L)%nr_PL > mrpcc ) then
          print *, 'mrpcc to low, increase mrpcc in modchem'
          stop
        end if
        PL_scheme(L)%PL(PL_scheme(L)%nr_PL)%r_nr = i   !store reaction number
        PL_scheme(L)%PL(PL_scheme(L)%nr_PL)%PorL = PRODUCTION   !this is a production reaction for this component
      endif
    enddo !j=1,reactions(i)%nr_chem_outp
  enddo !i=1,react

  nr_active_chemicals = k

  if (nr_active_chemicals < nchsp ) then
    print *, 'WARNING: More active chemicals specified in @ line then actually used.', nr_active_chemicals,' <', nchsp
  endif

  call reaction_location_mozart()  !only for ease of programming

  !Determine from the reactions which formula to use for all the producing and loss reactions
  !first do all reactions on the production side

  do i=1,nr_active_chemicals        !********************* misschien nchsp
    do j=1,PL_scheme(i)%nr_PL
      react_nr = PL_scheme(i)%PL(j)%r_nr
      if( PL_scheme(i)%PL(j)%PorL == PRODUCTION) then     !this is a PRODUCTION
        do k=1,reactions(react_nr)%nr_chem_outp
          if( reactions(react_nr)%outp(k)%name == PL_scheme(i)%name) then
            select case(reactions(react_nr)%nr_chem_inp)  !left of arrow 1 reactant
              case (1)
              icoeff(1) = int(reactions(react_nr)%inp(1)%coeff +0.05)
                 select case (icoeff(1))
                   case (1)
                     PL_scheme(i)%PL(j)%formula = 1
                     PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                   case (2)
                     PL_scheme(i)%PL(j)%formula = 2
                     PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                     PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(1)%index_sv0
                   case default
                     PL_scheme(i)%PL(j)%formula = 3
                     PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                     PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                  end select
                  PL_scheme(i)%PL(j)%coef =  reactions(react_nr)%outp(k)%coeff
              case(2)  !there are 2 reacting species
              icoeff(1) = int(reactions(react_nr)%inp(1)%coeff +0.05)
              icoeff(2) = int(reactions(react_nr)%inp(2)%coeff +0.05)
                if ((icoeff(1) == 1) .AND. (icoeff(2) == 1) ) then
                  PL_scheme(i)%PL(j)%formula = 2
                  PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                  PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                else
                  PL_scheme(i)%PL(j)%formula = 4
                  PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                  PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                  PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                  PL_scheme(i)%PL(j)%exp2  = icoeff(2)
                endif
                PL_scheme(i)%PL(j)%coef =  reactions(react_nr)%outp(k)%coeff
              case (3) !! there are 3 reacting species
               icoeff(1) = int(reactions(react_nr)%inp(1)%coeff +0.05)
               icoeff(2) = int(reactions(react_nr)%inp(2)%coeff +0.05)
               icoeff(3) = int(reactions(react_nr)%inp(3)%coeff +0.05)
               if( (icoeff(1)==1) .and. (icoeff(2) == 1) .and. (icoeff(3) == 1)) then
                  ! we don't need exponents keep it simple
                  PL_scheme(i)%PL(j)%formula = 5
                  PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                  PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                  PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(3)%index_sv0
                else
                  PL_scheme(i)%PL(j)%formula = 6
                  PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                  PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                  PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                  PL_scheme(i)%PL(j)%exp2  = icoeff(2)
                  PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(3)%index_sv0
                  PL_scheme(i)%PL(j)%exp3  = icoeff(3)
                endif
                PL_scheme(i)%PL(j)%coef =  reactions(react_nr)%outp(k)%coeff
              case (4)
               icoeff(1) = int(reactions(react_nr)%inp(1)%coeff +0.05)
               icoeff(2) = int(reactions(react_nr)%inp(2)%coeff +0.05)
               icoeff(3) = int(reactions(react_nr)%inp(3)%coeff +0.05)
               icoeff(4) = int(reactions(react_nr)%inp(4)%coeff +0.05)
               PL_scheme(i)%PL(j)%formula = 7
               PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
               PL_scheme(i)%PL(j)%exp1  = icoeff(1)
               PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
               PL_scheme(i)%PL(j)%exp2  = icoeff(2)
               PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(3)%index_sv0
               PL_scheme(i)%PL(j)%exp3  = icoeff(3)
               PL_scheme(i)%PL(j)%comp4 = reactions(react_nr)%inp(4)%index_sv0
               PL_scheme(i)%PL(j)%exp4  = icoeff(4)
               PL_scheme(i)%PL(j)%coef =  reactions(react_nr)%outp(k)%coeff
            end select
          endif
        enddo ! k=1,reactions(react_nr)%nr_chem_outp
      endif !( PL_scheme(i)%PL(j)%PorL == 1)
    enddo !j=1,PL_scheme(i)%nr_PL
  enddo !1,nr_active_chemicals

  !do all reactions on the loss side
  do i=1,nr_active_chemicals
    do j=1,PL_scheme(i)%nr_PL
      react_nr = PL_scheme(i)%PL(j)%r_nr
      if( PL_scheme(i)%PL(j)%PorL == LOSS) then     !This is LOSS
      do k=1,reactions(react_nr)%nr_chem_inp
        icoeff(1) = int(reactions(react_nr)%inp(1)%coeff +0.05)
        icoeff(2) = int(reactions(react_nr)%inp(2)%coeff +0.05)
        icoeff(3) = int(reactions(react_nr)%inp(3)%coeff +0.05)
        icoeff(4) = int(reactions(react_nr)%inp(4)%coeff +0.05)
        select case(reactions(react_nr)%nr_chem_inp)
		case (1) !the loss comp is the only reactant
          select case (icoeff(1))
            case (1)
              PL_scheme(i)%PL(j)%formula = 0
            case (2)
              PL_scheme(i)%PL(j)%formula = 1
              PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
            case (3)
              PL_scheme(i)%PL(j)%formula = 3
              PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
              PL_scheme(i)%PL(j)%exp1  = icoeff(1) - 1
          end select
          PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(1)%coeff
        case (2) ! we have 2 components in which one is current species
          if( reactions(react_nr)%inp(k)%name == PL_scheme(i)%name) then ! current selected
            select case (k)
            case(1)
              if( (icoeff(1) ==1) .and. (icoeff(2) == 1) ) then
                PL_scheme(i)%PL(j)%formula = 1
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(2)%index_sv0
              else if((icoeff(1) ==1) .and. (icoeff(2) > 1) ) then
                PL_scheme(i)%PL(j)%formula = 3
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(2)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(2)
              else if((icoeff(1) ==2) .and. (icoeff(2) == 1) ) then
                PL_scheme(i)%PL(j)%formula = 2
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
              else
                PL_scheme(i)%PL(j)%formula = 4
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(1) - 1
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                PL_scheme(i)%PL(j)%exp2  = icoeff(2)
              endif
              PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(1)%coeff
            case(2)
              if( (icoeff(1) ==1) .and. (icoeff(2) == 1) ) then
                PL_scheme(i)%PL(j)%formula = 1
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
              else if((icoeff(2) ==1) .and. (icoeff(1) > 1)) then
                PL_scheme(i)%PL(j)%formula = 3
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(1)
              else if((icoeff(2) ==2) .and. (icoeff(1) == 1) ) then
                PL_scheme(i)%PL(j)%formula = 2
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
              else
                PL_scheme(i)%PL(j)%formula = 4
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                PL_scheme(i)%PL(j)%exp2  = icoeff(2) - 1
              endif
              PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(2)%coeff
            end select
          endif
        case (3) !we have 3 components on input
          if( reactions(react_nr)%inp(k)%name == PL_scheme(i)%name) then ! current selected
            select case (k)
              case (1)
                PL_scheme(i)%PL(j)%formula = 6
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(2)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(2)
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(3)%index_sv0
                PL_scheme(i)%PL(j)%exp2  = icoeff(3)
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(1)%coeff
                if ( icoeff(1) == 1 ) then
                  PL_scheme(i)%PL(j)%formula = 4
                else
                  PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(1)%index_sv0
                  PL_scheme(i)%PL(j)%exp3  = icoeff(1) - 1
                endif
              case (2)
                PL_scheme(i)%PL(j)%formula = 6
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(3)%index_sv0
                PL_scheme(i)%PL(j)%exp2  = icoeff(3)
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(2)%coeff
                if ( icoeff(2) == 1 ) then
                  PL_scheme(i)%PL(j)%formula = 4
                else
                  PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(2)%index_sv0
                  PL_scheme(i)%PL(j)%exp3  = icoeff(2) -1
                endif
              case (3)
                PL_scheme(i)%PL(j)%formula = 6
                PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
                PL_scheme(i)%PL(j)%exp1  = icoeff(1)
                PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
                PL_scheme(i)%PL(j)%exp2  = icoeff(2)
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(3)%coeff
                if ( icoeff(3) == 1 ) then
                  PL_scheme(i)%PL(j)%formula = 4
                else
                  PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(3)%index_sv0
                  PL_scheme(i)%PL(j)%exp3  = icoeff(3) -1
                endif
            end select
		   endif
        case (4) !we have 4 components on input
          if( reactions(react_nr)%inp(k)%name == PL_scheme(i)%name) then ! current selected
            PL_scheme(i)%PL(j)%formula = 7
            PL_scheme(i)%PL(j)%comp1 = reactions(react_nr)%inp(1)%index_sv0
            PL_scheme(i)%PL(j)%exp1  = icoeff(1)
            PL_scheme(i)%PL(j)%comp2 = reactions(react_nr)%inp(2)%index_sv0
            PL_scheme(i)%PL(j)%exp2  = icoeff(2)
            PL_scheme(i)%PL(j)%comp3 = reactions(react_nr)%inp(3)%index_sv0
            PL_scheme(i)%PL(j)%exp3  = icoeff(3)
            PL_scheme(i)%PL(j)%comp4 = reactions(react_nr)%inp(4)%index_sv0
            PL_scheme(i)%PL(j)%exp4  = icoeff(4)
            select case(k)
              case (1)
                PL_scheme(i)%PL(j)%exp1  = icoeff(1) - 1
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(1)%coeff
              case (2)
                PL_scheme(i)%PL(j)%exp2  = icoeff(2) - 1
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(2)%coeff
              case (3)
                PL_scheme(i)%PL(j)%exp3  = icoeff(3) - 1
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(3)%coeff
              case (4)
                PL_scheme(i)%PL(j)%exp4  = icoeff(4) - 1
                PL_scheme(i)%PL(j)%coef = reactions(react_nr)%inp(4)%coeff
            end select
          endif
 	    end select
      enddo
      endif
    enddo
  enddo

!  open (15,file='PL_scheme.bin',form='binary')
 ! write(15)PL_scheme
 ! close(15)
 ! l=0
 ! open (15,file='PL_scheme.bin',form='binary')
 ! read(15)PL_temp
!  close(15)

  nr_raddep_react = l

  !refill PL_scheme %chem_number
  do i=1,nchsp
    do j=1,nchsp
      if(PL_scheme(i)%name == chem_name(j)) then
        !we don't need the orignal chem_number, so now use it as index to C_cbl
        PL_scheme(i)%chem_number = j
        exit
      end if
    end do
  end do

  ! the order of chemicals is in the order in which they appear in the chem reactions but we want them
  ! in the order of the @ line which is stored in chem_name()
  do i=1,nchsp
    found = .false.
    do j=1,nchsp
      if (chem_name(i) == PL_scheme(j)%name ) then
        PL_temp(i)= PL_scheme(j)
        found = .true.
        exit
      endif
    enddo
    if (found .EQV. .false.) then
      PL_temp(i)%name = chem_name(i)
    endif
  enddo

  PL_scheme = PL_temp

  deallocate (PL_temp)
  ! check which chemicals are really used
  do i=1, nchsp
    if (PL_scheme(i)%nr_PL == 0 ) then
      PL_scheme(i)%active = .FALSE.
    else
      PL_scheme(i)%active = .TRUE.
    endif
  enddo


  do i=1,nchsp
    PL_scheme(i)%prin = plot(i)
  enddo 

  !Print out the reaction schemes
  open (15,file=trim(outdir)//dirsep//'reaction_scheme',RECL=132)


  do i=1,nchsp
    write(15,'(a8,a4,i02)')PL_scheme(i)%name,'    ',PL_scheme(i)%nr_PL
  end do
  
  write(15,*)' '  
  do i=1,nchsp
  write(15,*)'---------------------------------------'
  write(15,*)' '
 if (PL_scheme(i)%active .EQV. .FALSE. ) then
    write(15,'(a6,a2,i2,a,a)')PL_scheme(i)%name,'(',i,')','NO REACTION'
    write(15,*)' '
  else
    write(15,'(a6,a2,i2,a)')PL_scheme(i)%name,'(',i,')'
    write(15,'(a5)') 'YP = '
    do j=1,PL_scheme(i)%nr_PL
      write(coef_str,'(f4.2)')PL_scheme(i)%PL(j)%coef
      if (PL_scheme(i)%PL(j)%PorL == PRODUCTION ) then
        select case (PL_scheme(i)%PL(j)%formula)
        case (0)
          write(15,*) '    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),')  ( F=',PL_scheme(i)%PL(j)%formula,')'
        case (1)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']','( F=',PL_scheme(i)%PL(j)%formula,')'
        case (2)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']', &
                    '* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ( F=',PL_scheme(i)%PL(j)%formula,')'
        case (3)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),'] ** (',PL_scheme(i)%PL(j)%exp1,')','( F=',PL_scheme(i)%PL(j)%formula,')'
        case(4)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),'] ** (',PL_scheme(i)%PL(j)%exp1,')', &
                    ' * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),']** (',PL_scheme(i)%PL(j)%exp2,') ( F=',PL_scheme(i)%PL(j)%formula,')'
        case(5)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']', &
                    '* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ( F=',PL_scheme(i)%PL(j)%formula,')'
        case(6)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']** (',PL_scheme(i)%PL(j)%exp1, &
                    ')* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ** (',PL_scheme(i)%PL(j)%exp2,')*  Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ** (',PL_scheme(i)%PL(j)%exp3,')( F=',PL_scheme(i)%PL(j)%formula,')'

        case(7)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']** (',PL_scheme(i)%PL(j)%exp1, &
                    ')* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ** (',PL_scheme(i)%PL(j)%exp2,')*  Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ** (',PL_scheme(i)%PL(j)%exp3, &
                    ')* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp4)),'] ** (',PL_scheme(i)%PL(j)%exp4,')( F=',PL_scheme(i)%PL(j)%formula,')'
        end select
      endif
    enddo

    write(15,'(a5)') 'YL = '
    do j=1,PL_scheme(i)%nr_PL
      if (PL_scheme(i)%PL(j)%PorL == LOSS ) then
        write(coef_str,'(f4.2)')PL_scheme(i)%PL(j)%coef
        select case (PL_scheme(i)%PL(j)%formula)
        case (0)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),')  ( F=',PL_scheme(i)%PL(j)%formula,')'
        case (1)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']','( F=',PL_scheme(i)%PL(j)%formula,')'
        case (2)
          write(15,*) '    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']', &
                    '* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ( F=',PL_scheme(i)%PL(j)%formula,')'
        case (3)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),'] ** (',PL_scheme(i)%PL(j)%exp1,') ( F=',PL_scheme(i)%PL(j)%formula,')'
        case(4)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),'] ** (',PL_scheme(i)%PL(j)%exp1,')', &
                    ' * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),']** (',PL_scheme(i)%PL(j)%exp2,') ( F=',PL_scheme(i)%PL(j)%formula,')'
        case(5)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']', &
                    '* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ( F=',PL_scheme(i)%PL(j)%formula,')'
        case(6)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']** (',PL_scheme(i)%PL(j)%exp1, &
                    ')* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ** (',PL_scheme(i)%PL(j)%exp2,')*  Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ** (',PL_scheme(i)%PL(j)%exp3,')( F=',PL_scheme(i)%PL(j)%formula,')'
        case(7)
          write(15,*)'    +',coef_str,' * K(',trim(RC(PL_scheme(i)%PL(j)%r_nr)%rname),') * Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp1)),']** (',PL_scheme(i)%PL(j)%exp1, &
                    ')* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp2)),'] ** (',PL_scheme(i)%PL(j)%exp2,')*  Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp3)),'] ** (',PL_scheme(i)%PL(j)%exp3, &
                    ')* Y[',trim(chem_name(PL_scheme(i)%PL(j)%comp4)),'] ** (',PL_scheme(i)%PL(j)%exp4,')( F=',PL_scheme(i)%PL(j)%formula,')'
        end select
      endif
    enddo
    write(15,*) ' '
  endif
  enddo

!  do i=1,nchsp
!  write(15,*)'---------------------------------------'
!  write(15,*)' '
!  if (PL_scheme(i)%active .EQV. .FALSE. ) then
!    write(15,*)PL_scheme(i)%name,'(',i,')'
!      write(15,*)' '
!  else
!    write(15,*)PL_scheme(i)%name,'(',i,')'
!    do j=1,PL_scheme(i)%nr_PL
!      if (PL_scheme(i)%PL(j)%PorL == 1 ) then
!        select case (PL_scheme(i)%PL(j)%formula)
!        case (0)
!          write(15,*)'F=',PL_scheme(i)%PL(j)%formula, '   YPtemp = YPtemp +',PL_scheme(i)%PL(j)%coef,' * K(',PL_scheme(i)%PL(j)%r_nr,')'
!        case (1)
!          write(15,*)'F=',PL_scheme(i)%PL(j)%formula, '   YPtemp = YPtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),']'
!        case (2)
!          write(15,*)'F=',PL_scheme(i)%PL(j)%formula, '   YPtemp = YPtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),']', &
!                    '* Y[',chem_name(PL_scheme(i)%PL(j)%comp2),']'
!        case (3)
!          write(15,*)'F=',PL_scheme(i)%PL(j)%formula, '   YPtemp = YPtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),'] ** (',PL_scheme(i)%PL(j)%exp1,')'
!        case(4)
!          write(15,*)'F=',PL_scheme(i)%PL(j)%formula, '   YPtemp = YPtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),'] ** (',PL_scheme(i)%PL(j)%exp1,')', &
!                    ' * Y[',chem_name(PL_scheme(i)%PL(j)%comp2),']** (',PL_scheme(i)%PL(j)%exp2,')'
!        end select
!      endif
!    enddo
!    write(15,*)'YP(',PL_scheme(i)%name,') = YPtemp'
!
!    do j=1,PL_scheme(i)%nr_PL
!      if (PL_scheme(i)%PL(j)%PorL == 2 ) then
!        select case (PL_scheme(i)%PL(j)%formula)
!        case (0)
!          write(15,*)'F=',PL_scheme(i)%PL(j)%formula, '   YLtemp = YLtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,')'
!        case (1)
!          write(15,*)'F=',PL_scheme(i)%PL(j)%formula, '   YLtemp = YLtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),']'
!        case (2)
!          write(15,*) 'F=',PL_scheme(i)%PL(j)%formula,'   YLtemp = YLtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),']', &
!                    '* Y[',chem_name(PL_scheme(i)%PL(j)%comp2),']'
!        case (3)
!          write(15,*)'F=',PL_scheme(i)%PL(j)%formula, '   YLtemp = YLtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),'] ** (',PL_scheme(i)%PL(j)%exp1,')'
!        case(4)
!          write(15,*)'F=',PL_scheme(i)%PL(j)%formula, '   YLtemp = YLtemp +',PL_scheme(i)%PL(j)%coef,' * K(',RC(PL_scheme(i)%PL(j)%r_nr)%rname,') * Y[',chem_name(PL_scheme(i)%PL(j)%comp1),'] ** (',PL_scheme(i)%PL(j)%exp1,')', &
!                    ' * Y[',chem_name(PL_scheme(i)%PL(j)%comp2),']** (',PL_scheme(i)%PL(j)%exp2,')'
!        end select
!      endif
!    enddo
!    write(15,*)'YL(',PL_scheme(i)%name,') = YLtemp'
!    write(15,*) ' '
!  endif
!  enddo
  close (15)

  goto 120
  100  print *,'error in inputchem filename= ',inputchemfile
    stop
  120  continue
  end SUBROUTINE

  subroutine read_chem_mozart(chem_name)
  use modchem
  implicit none

  character*8 ,dimension(nchsp)::chem_name
  character*255 scalarline
  integer i,j,status

  INERT%name  = 'INERT'
  PRODUC%name = 'PRODUC'
  O3%name     = 'O3'
  NO%name     = 'NO'
  NO2%name    = 'NO2'
  NO3%name    = 'NO3'
  N2O5%name   = 'N2O5'
  HNO3%name   = 'HNO3'
  HO2%name    = 'HO2'
  HO%name     = 'HO'
  H2O2%name   = 'H2O2'
  H2O%name    = 'H2O'
  CO%name     = 'CO'
  CO2%name    = 'CO2'
  RH%name     = 'ISOP'
  R%name      = 'R'
  NH3%name    = 'NH3'
  H2SO4%name  = 'H2SO4'
  ISO%name    = 'ISO'

  !set all 0 elements to 1. incase we do calculations with unknown componets
  c_cbl(0)=1.
  c_ft(0)=1.
  beta_ft(0)=1.
  Q_cbl(0)=1.
  E(0)=1.
  c_current(0)=1

  do i=1, nchsp
    if (O3%name    == chem_name(i)) then ; O3%loc   = i;  cycle; endif
    if (NO%name    == chem_name(i)) then ; NO%loc   = i;  cycle; endif
    if (NO2%name   == chem_name(i)) then ; NO2%loc  = i;  cycle; endif
    if (NO3%name   == chem_name(i)) then ; NO3%loc  = i;  cycle; endif
    if (N2O5%name  == chem_name(i)) then ; N2O5%loc = i;  cycle; endif
    if (HNO3%name  == chem_name(i)) then ; HNO3%loc = i;  cycle; endif
    if (HO2%name   == chem_name(i)) then ; HO2%loc  = i;  cycle; endif
    if (HO%name    == chem_name(i)) then ; HO%loc   = i;  cycle; endif
    if (H2O2%name  == chem_name(i)) then ; H2O2%loc = i;  cycle; endif
    if (H2O%name   == chem_name(i)) then ; H2O%loc  = i;  cycle; endif
    if (CO%name    == chem_name(i)) then ; CO%loc   = i;  cycle; endif
    if (CO2%name   == chem_name(i)) then ; CO2%loc  = i;  cycle; endif
    if (RH%name    == chem_name(i)) then ; RH%loc   = i;  cycle; endif
    if (ISO%name   == chem_name(i)) then ; ISO%loc  = i;  cycle; endif
    if (R%name     == chem_name(i)) then ; R%loc    = i;  cycle; endif
    if (NH3%name   == chem_name(i)) then ; NH3%loc  = i;  cycle; endif
    if (H2SO4%name == chem_name(i)) then ; H2SO4%loc  = i;  cycle; endif
    if (INERT%name == chem_name(i)) then ; INERT%loc  = i;  cycle; endif
    if (PRODUC%name == chem_name(i)) then ; PRODUC%loc  = i;  cycle; endif
  enddo

!  read(10,'(a)',err=400)scalarline
 ! read(scalarline,*)(c_cbl(j),j=1,nchsp)

!  read(10,'(a)',err=400)scalarline
!  read(scalarline,*)(c_ft(j),j=1,nchsp)

!  read(10,'(a)',err=400)scalarline
!  read(scalarline,*)(Q_init(j),j=1,nchsp)
!  Q_cbl = Q_init

!  read(10,'(a)',err=400)scalarline
!  read(scalarline,*)(Q_func(j),j=1,nchsp)

  do i=1,nchsp
    write(6,'(i3,x,a5,x,3E13.5)') i,chem_name(i),c_cbl(i),Q_cbl(i),c_ft(i)
  enddo

  goto 500
  400  print *, 'error in reading chem species in inputchem'
  500  print *,''

  end


  subroutine allocate_arrays1()
  use modchem
  implicit none
  integer i

  ! the arrays are declared from 0 to have a extra space beacuse all unknown components and reactions point with the %loc to 0
  ! so if we use an unidentified reactions and/or component in a reaction we don't have problems. The value of concentrations will be set to 1
  allocate (plot(0:nchsp))
  allocate (c_cbl(0:nchsp),c_ft(0:nchsp))
  allocate (beta_ft(0:nchsp))
  allocate (Q_cbl(0:nchsp),Q_init(0:nchsp),E(0:nchsp))
  allocate (Q_func(0:nchsp))
  allocate (c_current(0:nchsp))
  allocate (PL_scheme(nchsp), PL_temp(nchsp))
  allocate (productionloss(nchsp,mrpcc+2))!+2 for total production and loss terms

  PL_scheme(1)%name = '     '
  PL_scheme(1)%active = .false.
  PL_scheme(1)%prin = .false.
  PL_scheme(1)%chem_number = 0
  PL_scheme(1)%nr_PL = 0

  do i=1,mrpcc
    PL_scheme(1)%PL(i)%formula = 0
    PL_scheme(1)%PL(i)%r_nr = 0
    PL_scheme(1)%PL(i)%PorL = 0
    PL_scheme(1)%PL(i)%coef = 0
    PL_scheme(1)%PL(i)%comp1 = 0
    PL_scheme(1)%PL(i)%exp1 = 0
    PL_scheme(1)%PL(i)%comp2 = 0
    PL_scheme(1)%PL(i)%exp2 = 0
    PL_scheme(1)%PL(i)%comp3 = 0
    PL_scheme(1)%PL(i)%exp3 = 0
    PL_scheme(1)%PL(i)%comp4 = 0
    PL_scheme(1)%PL(i)%exp4 = 0
  enddo

  do i=2, nchsp
      pl_scheme(i)=pl_scheme(1)
  enddo
  pl_temp=pl_scheme
  end subroutine

  subroutine allocate_arrays2()
  use modchem
  implicit none
  integer i
    allocate (reactions(tnor))
    allocate (RC(0:tnor))
  end subroutine


  subroutine reaction_location_mozart()
  use modchem
  implicit none

  integer i

  RC(0)%Kreact=1.
  RC(0)%Keff_cbl=1.
  RC(0)%Keff_ft=1.

  R_O3%name   = 'R_O3'
  R_NO%name   = 'R_NO'
  R_NO2%name  = 'R_NO2'
  R_RH%name   = 'R_RH'
  R_CO%name   = 'R_CO'
  R_HNO3%name = 'R_HNO3'
  R_1%name    = 'R_1'
  R_23%name   = 'R_23'
  R_25%name   = 'R_25'
  R_26a%name  = 'R_26'
  R_28%name   = 'R_28'
  R_43%name   = 'R_43'
  R_45%name   = 'R_45'
  R_54A%name  = 'R_54A'
  R_57A%name  = 'R_57A'
  R_58A%name  = 'R_58A'
  R_61A%name  = 'R_61A'
  R_62A%name  = 'R_62A'
  R_63Aa%name = 'R_63A'
  R_63Ab%name = 'R_63Ab'
  R_CH2O%name = 'R_CH2O'

  do i=1, tnor
    if( reactions(i)%name == R_O3%name )   then ; R_O3%loc   = i;  cycle; endif
    if( reactions(i)%name == R_NO%name )   then ; R_NO%loc   = i;  cycle; endif
    if( reactions(i)%name == R_NO2%name )  then ; R_NO2%loc  = i;  cycle; endif
    if( reactions(i)%name == R_RH%name )   then ; R_RH%loc   = i;  cycle; endif
    if( reactions(i)%name == R_CO%name )   then ; R_CO%loc   = i;  cycle; endif
    if( reactions(i)%name == R_HNO3%name ) then ; R_HNO3%loc = i;  cycle; endif
    if( reactions(i)%name == R_1%name )    then ; R_1%loc    = i;  cycle; endif
    if( reactions(i)%name == R_23%name )   then ; R_23%loc   = i;  cycle; endif
    if( reactions(i)%name == R_25%name )   then ; R_25%loc   = i;  cycle; endif
    if( reactions(i)%name == R_26a%name )  then ; R_26a%loc  = i;  cycle; endif
    if( reactions(i)%name == R_28%name )   then ; R_28%loc   = i;  cycle; endif
    if( reactions(i)%name == R_43%name )   then ; R_43%loc   = i;  cycle; endif
    if( reactions(i)%name == R_45%name )   then ; R_45%loc   = i;  cycle; endif
    if( reactions(i)%name == R_54A%name )  then ; R_54A%loc  = i;  cycle; endif
    if( reactions(i)%name == R_57A%name )  then ; R_57A%loc  = i;  cycle; endif
    if( reactions(i)%name == R_58A%name )  then ; R_58A%loc  = i;  cycle; endif
    if( reactions(i)%name == R_61A%name )  then ; R_61A%loc  = i;  cycle; endif
    if( reactions(i)%name == R_62A%name )  then ; R_62A%loc  = i;  cycle; endif
    if( reactions(i)%name == R_63Aa%name ) then ; R_63Aa%loc = i;  cycle; endif
    if( reactions(i)%name == R_63Ab%name ) then ; R_63Ab%loc = i;  cycle; endif
    if( reactions(i)%name == R_CH2O%name ) then ; R_CH2O%loc = i;  cycle; endif
  enddo

  end subroutine

