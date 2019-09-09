 ! *********************************************************************************************
 ! This module defines the trial function, the drift force and the local energy of the 
 ! washboard potential with few wells
 ! ******************************************************************************************* *

MODULE functions

    Use mpi
    Implicit none

    INTEGER, PARAMETER    ::  N1 = 100100, N2 = 100, ncorr0=2000  
    REAL(8), PUBLIC       ::  pi=3.14159265358979323846
    INTEGER ::  Lx 
    INTEGER, PARAMETER    ::        Ly = 1
    INTEGER, PUBLIC  ::     nwalk, ncorr 

    REAL(8), DIMENSION(N2)     :: ist
    REAL(8), PUBLIC, DIMENSION(N1)     :: mag, mag_new, En, Es
    REAL(8), PUBLIC, DIMENSION(:), ALLOCATABLE     :: Eo, Eo_r, Eo_l, En_l, En_r    ! Modified by metro hidden
    REAL(8), PUBLIC, DIMENSION(:,:), ALLOCATABLE  :: spin, spin_new
    REAL(8), PUBLIC, DIMENSION(:,:), ALLOCATABLE  :: lspin, rspin

    REAL(8)         :: e, p, q, smin, smax, ds
    REAL(8)         :: dt_rel, dt_in, lweight, rweight, wi, et, mu

    INTEGER         :: igg, iwalk, idelta,  nwalk1, imult, idx1, idx2,ispins,i2,i, j, k, usedsize, tau, Nspins
    INTEGER, PUBLIC :: nstep1, nstep2, igen, count, Maxsons, nsample, checkb, checka


    INTEGER,PUBLIC :: m, mes, wal
    REAL(8),PUBLIC :: a_mag, a_eng, hlong, hfield!,der_locE_r, der_locE_s, der_locE_rs 
    REAL(8),DIMENSION(4),PUBLIC :: a, wa
    REAL(8),PUBLIC :: dt !, beta_r, beta_s, Jrs
    LOGICAL,PUBLIC :: key 

    REAL(8),  PUBLIC ::  Jo  
    PUBLIC :: epot, ediff, metropolis_real, metropolis_hidden
	 

    PUBLIC :: initialize ,vmc, is_weight, sgd, rand

    INTEGER, ALLOCATABLE, DIMENSION(:)  :: SEED

    integer :: seed1                                                                                                   
    integer(kind=4) errcode
    !REAL(8), EXTERNAL :: rand  

    REAL(kind=4) :: aran = 0.e0,bran=1.e0
    INTEGER(KIND=4) :: isone4, Lxpo4,Lx4
    INTEGER :: isone

    INTEGER, PARAMETER :: nnearest = 2
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: isnear

    REAL(8), DIMENSION(:,:), ALLOCATABLE ::sconftime
    INTEGER :: itaumax = 100,igroup=100
    REAL(8), DIMENSION(:), ALLOCATABLE :: Ctau
    INTEGER :: imeasured, idum_seed, cnt
    REAL(8), PUBLIC :: avg_clas,ls_avg_clas,rs_avg_clas
    REAL(8), PUBLIC :: ls_avg_mcla,rs_avg_mcla,ls_avg_eloc,rs_avg_eloc

    real(8), public, dimension(3) :: lambda, der_lambda

    real(8), public :: beta_r,beta_s,Jrs 
    real(8), public, dimension(12)  :: params


    ! MPI stuff
    integer :: ierr
    integer :: numtasks, rank
    integer :: loc_size, low_bound, up_bound, mid_bound, rest, offset

    CONTAINS

    ! ********************************************************************************************
    SUBROUTINE initialize()
        INTEGER :: iaux, iaux2
        Real(8) :: lEo, rEo  
        
        cnt = 0
        ! Boundaries of the histogram
        smin = -1.d0
        smax =  1.d0

        ! Width of a bin in the histogram
        ds = (smax-smin)/DBLE(N2)

        ! Parameters initialization

        if (rank == 0) then 
            OPEN(UNIT=17,FILE='parameters.txt', STATUS='unknown')
            READ(17,*)
            !READ(17,*) Jo, hfield, hlong, nstep1, nstep2, nwalk, mu, Lx, seed1 , beta_r, beta_s, Jrs
            read(17,*) params(:)
            CLOSE(17)
            !Debug:  ! write(*,*) params
        end if

        call MPI_BCAST( params, 12, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr )

        !if (ierr .ne. 0 ) then
        !    call MPI_ABORT(MPI_COMM_WORLD, "Error in params BCAST: ", ierr)
        !end if

        Jo = params(1)
        hfield = params(2)
        hlong = params(3)
        nstep1 = int(params(4))
        nstep2 = int(params(5))
        nwalk = int(params(6))
        mu = params(7)
        Lx = int(params(8))
        seed1 = int(params(9))
        beta_r = params(10)
        beta_s = params(11)
        Jrs = params(12)

        !if( rank ==0 ) then
        !    print *, ierr
        !end if
        !stop

        isone = 1
        isone4 = 1
        Lxpo4 = Lx+1
        Lx4 = Lx

        ncorr = nstep2/4

        IF (ALLOCATED(spin))     DEALLOCATE(spin)
        IF (ALLOCATED(lspin))     DEALLOCATE(lspin)
        IF (ALLOCATED(rspin))     DEALLOCATE(rspin)
        IF (ALLOCATED(isnear))     DEALLOCATE(isnear)

        ALLOCATE(spin(Lx,nwalk))
        ALLOCATE(lspin(Lx,nwalk))
        ALLOCATE(rspin(Lx,nwalk))
        ALLOCATE(isnear(nnearest,Lx))
        ALLOCATE(Eo(nwalk))
        ALLOCATE(Eo_l(nwalk))
        ALLOCATE(Eo_r(nwalk))




        DO iaux = 1, Lx
            iaux2 = iaux+1
            IF (iaux2 .GT. Lx) iaux2 = 1
            isnear(2, iaux) = iaux2
            iaux2 = iaux-1
            IF (iaux2 .LT. 1) iaux2 = Lx
            isnear(1, iaux) = iaux2
        END DO


        imeasured = 0

        if (rank == 0) then 
            write(*,*)""
            PRINT*, "The strength of nn spins interaction is:", Jo
            PRINT*, "The strength of the transverse field is:", hfield
            PRINT*, "The strength of the longitudinal field is:", hlong
            PRINT*, "Number of steps for MCs:", nstep1
            PRINT*, "Number of SGD steps:", nstep2
            PRINT*, "Target number of walkers:", nwalk
            PRINT*, "Length of the spin chain:", Lx
            write(*,*)""
        end if

        Nspins = Lx*Ly

        Maxsons = nwalk
        nwalk1 = nwalk
	  
        lambda(1) = beta_r
        lambda(2) = beta_s
        lambda(3) = Jrs

        idum_seed = seed1

	    !  Magnetization and energy initialization
        mag = 0.d0
        e = 0.d0
        En = 0.d0
        Eo = 0.d0
        lEo = 0.d0
        rEo = 0.d0

	    ! Initialize the matrix containing some important physical values 
        wa = 0.d0

        ! Initialize the arrays for histogram 
        ist = 0.d0

        count=0

    END SUBROUTINE initialize
!********************************************************************************************

!!CC{{{  Functions
    FUNCTION rand()
    INTEGER IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV, seed_rnd
    DOUBLE PRECISION rand,AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,&
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
    NTAB=32,NDIV=1+IMM1/NTAB,EPS=3.d-16,RNMX=1.d0-EPS)
    INTEGER idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/

    cnt = cnt + 1
    seed_rnd = idum_seed

    if ( idum_seed .le. 0 ) then
        idum_seed=max(-idum_seed,1)
        idum2=idum_seed

        do 11 j = NTAB+8, 1, -1

            k=idum_seed/IQ1
            idum_seed=IA1*(idum_seed-k*IQ1)-k*IR1
            if (idum_seed.lt.0) idum_seed=idum_seed+IM1
            if (j.le.NTAB) iv(j)=idum_seed

11      continue

        iy=iv(1)

    endif

    k=idum_seed/IQ1
    idum_seed=IA1*(idum_seed-k*IQ1)-k*IR1

    if (idum_seed.lt.0) idum_seed=idum_seed+IM1

    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2

    if (idum2.lt.0) idum2=idum2+IM2

    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum_seed

    if(iy.lt.1)iy=iy+IMM1

    rand=min(AM*iy,RNMX)
    !write(*,*)"seed : ", seed_rnd, "Random No: ", rand

    return
    end function
!C  (C) Copr. 1986-92 Numerical Recipes Software (9`3j32150.
!CC}}}

!********************************************************************************************
! This subroutine performs VMC
!********************************************************************************************

    subroutine vmc( beta_r, beta_s, Jrs, energy, energy_err, derivative )
        REAL(8), INTENT(IN) :: beta_r, beta_s, Jrs
        real(8), intent(out) :: energy, energy_err
        real(8), intent(out), dimension(3):: derivative 
        real(8) :: eebavg, eebavg2, var_E
        real(8) :: rn
        real(8), dimension(3) :: der, der_var_E
        integer :: iibavg, it

        eebavg=0.0
        eebavg2=0.0
        der = 0.0
        iibavg = 0
       
        DO it = 1, nstep1
            rn = rand()
        
            IF ( rn .lt. 0.5) THEN
            ! Move only shadow spins
                call metropolis_hidden( beta_r, beta_s, Jrs, var_E, der_var_E)
            ELSE 
            ! Move only physical spin
                CALL metropolis_real( beta_r, Jrs, var_E, der_var_E)
            END IF
        

            IF ( ( MOD( it , ncorr0 ) == 0 ) .AND. ( it .GT. ncorr0 ) )THEN
                eebavg  = eebavg  + var_E 
                eebavg2 = eebavg2 +(var_E)**2.d0
                der = der + der_var_E 
                iibavg  = iibavg  + 1
            END IF
      
            WRITE(9,'(i10,2F22.8)')  it, var_E
        END DO  
       
       ! Results 
        energy = eebavg/dble(iibavG)
        energy_err = sqrt((eebavg2/dble(iibavg) - (eebavg/dble(iibavg))**2)/dble(iibavg-1)) 

        derivative = der/dble(iibavg)

        print*, 'energy', energy, '+/-',energy_err 
    end subroutine vmc 


    ! ********************************************************************************************
    ! This subroutine performs Metropolis Algorithm on hidden spins
    ! ********************************************************************************************

    subroutine metropolis_hidden(beta_r, beta_s, Jrs, var_E, der_var_E)
        REAL(8), INTENT(IN) :: beta_s, beta_r, Jrs
        real(8), intent(out) :: var_E
        real(8), intent(out), dimension(3) :: der_var_E
        integer :: iwalk
        REAL(8) :: deltaE_r, deltaE_l,prob_l,prob_r
        REAL(8) :: ds_l, ds_r
        INTEGER :: stobemoved_l, stobemoved_r
        REAL(8) :: ls_eloc, rs_eloc ,avg_clas, ls_avg_clas, rs_avg_clas
        REAL(8) :: ls_avg_mcla,rs_avg_mcla,ls_avg_eloc,rs_avg_eloc
        Real(8) :: ecum1,ecum2,ecum3,ecum4,scum1,scum2,rscum1,rscum2
        real(8) ::  lweight,rweight,lEcl_rs,rEcl_rs
        real(8) :: rn
        Real(8) :: der_locE_r , der_locE_rs, der_locE_s

        avg_clas   =0.d0
        ls_avg_clas=0.d0
        rs_avg_clas=0.d0    
        ls_avg_mcla=0.d0
        rs_avg_mcla=0.d0   
        ls_avg_eloc=0.d0
        rs_avg_eloc=0.d0


        ecum1 = 0.d0 
        ecum2 = 0.d0
        ecum3 = 0.d0
        ecum4 = 0.d0
        scum1 = 0.d0
        scum2 = 0.d0

        rscum1 = 0.d0
        rscum2 = 0.d0

        !  print*, 'b', ecum1

        ! Move a walker
        DO iwalk = 1, nwalk

            stobemoved_l = INT((rand() * Nspins) + 1.d0) 

            lspin(stobemoved_l,iwalk) = -lspin(stobemoved_l,iwalk)

            deltaE_l = ediff(lspin(1:Nspins,iwalk),stobemoved_l)

            ds_l     = 2.d0* Jrs * spin(stobemoved_l,iwalk) * lspin(stobemoved_l,iwalk)

            prob_l   = exp(-beta_s * deltaE_l + ds_l)

            rn = rand()

            IF( prob_l .GE. 1.D0 )THEN
                Eo_l(iwalk) = Eo_l(iwalk) + deltaE_l
            ELSE IF(DBLE(rn) .LT. prob_l)THEN
                Eo_l(iwalk) = Eo_l(iwalk) + deltaE_l
            ELSE IF(DBLE(rn) .GE. prob_l)THEN
                lspin(stobemoved_l,iwalk) = -lspin(stobemoved_l,iwalk)
            END IF

            stobemoved_r = INT((rand() * Nspins) + 1.d0)

            rspin(stobemoved_r,iwalk) = -rspin(stobemoved_r,iwalk)

            deltaE_r = ediff(rspin(1:Nspins,iwalk),stobemoved_r)
            ds_r     = 2.d0*Jrs*spin(stobemoved_r,iwalk)*rspin(stobemoved_r,iwalk)
            prob_r   = exp(-beta_s*deltaE_r+ds_r)

            rn = rand()
            IF(prob_r .GE. 1.D0 )THEN
                Eo_r(iwalk) = Eo_r(iwalk) + deltaE_r
            ELSE IF(DBLE(rn) .LT. prob_r)THEN
                Eo_r(iwalk) = Eo_r(iwalk) + deltaE_r
            ELSE IF(DBLE(rn) .GE. prob_r)THEN
                rspin(stobemoved_r,iwalk) = -rspin(stobemoved_r,iwalk)
            END IF
      
            !***** CUMULATE DATA ****************************************
            CALL is_weight( iwalk, lweight, rweight, lEcl_rs, rEcl_rs, beta_r, Jrs)

            ls_eloc = -hfield*lweight   + Eo(iwalk)
            rs_eloc = -hfield*rweight   + Eo(iwalk)


            avg_clas  =    avg_clas  + Eo(iwalk)
            ls_avg_clas  = ls_avg_clas  + Eo_l(iwalk)
            rs_avg_clas  = rs_avg_clas  + Eo_r(iwalk)
         
            ls_avg_mcla  = ls_avg_mcla  + lEcl_rs
            rs_avg_mcla  = rs_avg_mcla  + rEcl_rs

            ls_avg_eloc = ls_avg_eloc  + ls_eloc
            rs_avg_eloc = rs_avg_eloc  + rs_eloc

      

            ecum1 = ecum1 + 0.5d0*(ls_eloc+rs_eloc)            
                    
            ecum2 = ecum2 + (Eo(iwalk))*rs_eloc
            ecum3 = ecum3 + (Eo(iwalk))*ls_eloc

            scum1 = scum1  + (Eo_l(iwalk))*rs_eloc
            scum2 = scum2  + (Eo_r(iwalk))*ls_eloc

            rscum1 = rscum1 + lEcl_rs*rs_eloc
            rscum2 = rscum2 + rEcl_rs*ls_eloc

        END DO


    der_locE_r  = (avg_clas*rs_avg_eloc)/dble((nwalk)**2) - (ecum2)/dble(nwalk) + &
    (avg_clas*ls_avg_eloc)/dble((nwalk)**2) - (ecum3)/dble(nwalk)


    der_locE_s  = (ls_avg_clas*rs_avg_eloc)/dble((nwalk)**2) - (scum1)/dble(nwalk) + &
    (rs_avg_clas*ls_avg_eloc)/dble((nwalk)**2) - (scum2)/dble(nwalk)

    der_locE_rs = (ls_avg_mcla*rs_avg_eloc)/dble((nwalk)**2) - (rscum1)/dble(nwalk) + &
    (rs_avg_mcla*ls_avg_eloc)/dble((nwalk)**2) - (rscum2)/dble(nwalk)


    var_E = ecum1/dble(nwalk)

    der_var_E(1) = der_locE_r
    der_var_E(2) = der_locE_s
    der_var_E(3) = der_locE_rs 

end subroutine metropolis_hidden

! ********************************************************************************************
! This subroutine performs Metropolis Algorithm on real spins
! ********************************************************************************************

SUBROUTINE metropolis_real(beta_r,Jrs,var_E, der_var_E)
REAL(8), INTENT(IN) :: beta_r,Jrs
real(8), intent(out) :: var_E
real(8), intent(out), dimension(3) :: der_var_E
INTEGER :: iwalk
REAL(8) :: enew,prob,deltaE,dshadow,lEcl_rs,rEcl_rs
INTEGER :: imoveact
REAL(8) :: prn
REAL(8) :: ls_eloc,rs_eloc ,avg_clas,ls_avg_clas,rs_avg_clas
REAL(8) :: ls_avg_mcla,rs_avg_mcla,ls_avg_eloc,rs_avg_eloc
Real(8) :: ecum1,ecum2,ecum3,ecum4,scum1,scum2,rscum1,rscum2 

Real(8) :: der_locE_r , der_locE_rs, der_locE_s
Real(8) :: magn1  

    avg_clas=0.d0
    ls_avg_clas=0.d0
    rs_avg_clas=0.d0    
    ls_avg_mcla=0.d0
    rs_avg_mcla=0.d0   
    ls_avg_eloc=0.d0
    rs_avg_eloc=0.d0


    ecum1 = 0.d0 
    ecum2 = 0.d0
    ecum3 = 0.d0
    ecum4 = 0.d0
    scum1 = 0.d0
    scum2 = 0.d0

    rscum1 = 0.d0
    rscum2 = 0.d0

    enew  = 0.d0
    magn1 = 0.d0
    count = 0
       

    ! Move a walker
    do iwalk = 1, nwalk
               
        imoveact = INT((rand() * Nspins) + 1.d0)

        spin(imoveact,iwalk) = -spin(imoveact,iwalk)

        deltaE  = ediff(spin(1:Nspins,iwalk),imoveact)    
        dshadow = 2.d0*Jrs*spin(imoveact,iwalk)*(lspin(imoveact,iwalk)+rspin(imoveact,iwalk))
 
        prob   = exp(-2.d0*beta_r*deltaE+dshadow) 
     
        prn = rand() 

        IF(prob .GE. 1.D0 )THEN
            ! mag(iwalk)  = mag(iwalk) + 2.d0*spin(imoveact,iwalk) 
            Eo(iwalk) = Eo(iwalk) + deltaE
            count = count + 1      
        ELSE IF(DBLE(prn) .LT. prob)THEN 
            ! mag(iwalk)  = mag(iwalk) + 2.d0*spin(imoveact,iwalk)
            Eo(iwalk) = Eo(iwalk) + deltaE
            count = count + 1
        ELSE IF(DBLE(prn) .GE. prob)THEN
            spin(imoveact,iwalk) = -spin(imoveact,iwalk)
        END IF  
      
        !***** CUMULATE DATA ****************************************
        CALL is_weight(iwalk,lweight,rweight,lEcl_rs,rEcl_rs,beta_r,Jrs)

        ls_eloc = -hfield*lweight   + Eo(iwalk)
        rs_eloc = -hfield*rweight   + Eo(iwalk)


        avg_clas  =    avg_clas  + Eo(iwalk)
        ls_avg_clas  = ls_avg_clas  + Eo_l(iwalk)
        rs_avg_clas  = rs_avg_clas  + Eo_r(iwalk)
     
        ls_avg_mcla  = ls_avg_mcla  + lEcl_rs
        rs_avg_mcla  = rs_avg_mcla  + rEcl_rs

        ls_avg_eloc = ls_avg_eloc  + ls_eloc
        rs_avg_eloc = rs_avg_eloc  + rs_eloc

      

        ecum1 = ecum1 + 0.5d0*(ls_eloc+rs_eloc)            
                
        ecum2 = ecum2 + (Eo(iwalk))*rs_eloc
        ecum3 = ecum3 + (Eo(iwalk))*ls_eloc

        scum1 = scum1  + (Eo_l(iwalk))*rs_eloc
        scum2 = scum2  + (Eo_r(iwalk))*ls_eloc

        rscum1 = rscum1 + lEcl_rs*rs_eloc
        rscum2 = rscum2 + rEcl_rs*ls_eloc
      
        end do

   
    der_locE_r  = (avg_clas*rs_avg_eloc)/dble((nwalk)**2) - (ecum2)/dble(nwalk) + &
    (avg_clas*ls_avg_eloc)/dble((nwalk)**2) - (ecum3)/dble(nwalk)


    der_locE_s  = (ls_avg_clas*rs_avg_eloc)/dble((nwalk)**2) - (scum1)/dble(nwalk) + &
    (rs_avg_clas*ls_avg_eloc)/dble((nwalk)**2) - (scum2)/dble(nwalk)


    der_locE_rs = (ls_avg_mcla*rs_avg_eloc)/dble((nwalk)**2) - (rscum1)/dble(nwalk) + &
    (rs_avg_mcla*ls_avg_eloc)/dble((nwalk)**2) - (rscum2)/dble(nwalk)


    var_E = ecum1/dble(nwalk)

    der_var_E(1)=der_locE_r
    der_var_E(2)=der_locE_s
    der_var_E(3)=der_locE_rs 

END SUBROUTINE metropolis_real



! ******************************************************************************************
! This subroutine performs stochastic gradient descend
! ****************************************************************************************

  subroutine sgd(beta_r,beta_s,Jrs,energy,energy_err,derivative, mu_t)
   REAL(8), INTENT(inout) :: beta_s, beta_r, Jrs
   real(8), intent(inout) :: energy, energy_err
   real(8), intent(inout), dimension(3) :: derivative 
   !real(8) :: mu_t
   real(8), intent(in) :: mu_t
 
   call vmc(beta_r,beta_s,Jrs,energy,energy_err,derivative)
   !mu_t    = (1.d0 - dble(count)/dble(nstep2)) !**(0.7)     
   beta_r  = beta_r - mu_t*derivative(1)
   beta_s  = beta_s - mu_t*derivative(2)
   Jrs     = Jrs    - mu_t*derivative(3)
   
  end subroutine sgd
! *********************************************************************************************


! *********************************************************************************************
! This subroutine computes the importance sampling weight
! *********************************************************************************************
SUBROUTINE is_weight(iw,lwt,rwt,lEcl_rs,rEcl_rs,beta_r,Jrs)
REAL(8), INTENT(OUT) :: lwt,rwt,lEcl_rs,rEcl_rs
REAL(8)              :: dE,lds,rds,beta_r,Jrs
INTEGER, INTENT(IN)  :: iw
INTEGER              :: is

    lwt   = 0.d0
    rwt   = 0.d0 
    lEcl_rs = 0.d0
    rEcl_rs = 0.d0

    do is = 1, Nspins

        spin(is,iw) = -spin(is,iw)
        dE = ediff( spin(1:Nspins,iw), is) 
        lds = spin(is,iw) * lspin(is,iw)
        rds = spin(is,iw) * rspin(is,iw) 
        lEcl_rs = lEcl_rs + lds  
        rEcl_rs = rEcl_rs + rds
        lwt  = lwt + exp(-beta_r*dE+2.d0*Jrs*lds)
        rwt  = rwt + exp(-beta_r*dE+2.d0*Jrs*rds)
        spin(is,iw) = -spin(is,iw)

    end do
 
  
  END SUBROUTINE is_weight
! ********************************************************************************************



!********************************************************************************************
!This function computes the potential energy for one copy of the system (walker)
!********************************************************************************************
REAL(8)  FUNCTION epot(spin,iwalk) RESULT(Y)
REAL(8), DIMENSION(:,:) , INTENT(IN) :: spin
INTEGER, INTENT(IN) :: iwalk 
INTEGER :: i 
REAL(8) ::  E

    E = 0.d0 

    !No Periodic Boundary Conditions
    DO i = 1, Lx-1
        E = E - Jo * (spin(i,iwalk)) * (spin(i+1,iwalk)) !- hlong*spin(i,iwalk)
    END DO

    Y = E !- Jo*(spin(Lx,iwalk))*(spin(1,iwalk))  !- hlong*spin(Lx,iwalk)

END FUNCTION epot

!********************************************************************************************
!This function computes the diiference in potential energy for one copy of the system
!(walker)
!********************************************************************************************
REAL(8)  FUNCTION ediff(spinsave,imoveact) RESULT(Y)

REAL(8), DIMENSION(Nspins) , INTENT(IN) :: spinsave
INTEGER :: j,iact,jact
REAL(8) ::  E
INTEGER, INTENT(IN) :: imoveact

    E = 0.d0

    iact = imoveact

    if ( iact == 1 )then
        E = E - Jo * (spinsave(iact)) * (spinsave(iact+1))
    else if (iact == Nspins) then
        E = E - Jo * (spinsave(iact)) * (spinsave(iact-1))
    else
        do j = 1, nnearest
            jact = isnear(j,iact)
            E = E - Jo*(spinsave(iact))*(spinsave(jact)) !- hlong*spinsave(iact)
        end do
    end if
        
    Y = 2.d0 * E
END FUNCTION ediff

! ********************************************************************************************
END MODULE functions



