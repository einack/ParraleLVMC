 ! *********************************************************************************************
 ! This module defines the trial function, the drift force and the local energy of the 
 ! washboard potential with few wells
 ! ******************************************************************************************* *

MODULE functions
    Implicit none

    INTEGER, PARAMETER      :: N2 = 100
    INTEGER, PARAMETER    ::  N1 = 100100, ncorr0=2000
    REAL(8), PUBLIC       ::  pi=3.14159265358979323846
    INTEGER, save,PUBLIC ::        Lx=30 
    INTEGER, PARAMETER    ::        Ly = 1
    INTEGER, PUBLIC  ::     nwalk, ncorr 

    REAL(8), DIMENSION(N2)     :: ist
    REAL(8), PUBLIC, DIMENSION(N1)     :: mag, mag_new, En, Eo, Es
    REAL(8), PUBLIC, DIMENSION(N1)     :: Eo_r, Eo_l, En_l, En_r
    REAL(8), PUBLIC, DIMENSION(:,:), ALLOCATABLE  :: spin, spin_new, spin_old
    REAL(8), PUBLIC, DIMENSION(:,:), ALLOCATABLE  :: lspin, rspin

    REAL(8)         :: e, p, q, smin, smax, ds
    REAL(8)         :: dt_rel, dt_in, lweight, rweight, wi, et, mu
    !REAL(8), PUBLIC :: ecum1, ecum2, magn1, ecum3, ecum4,  magn, magn2, magn4, Pa, hfi

    INTEGER         :: igg, iwalk, idelta,  nwalk1, imult, idx1, idx2,ispins,i2,i, j, k, usedsize, tau, Nspins
    INTEGER, PUBLIC :: nstep1, nstep2, igen, count, Maxsons, nsample, checkb, checka


    INTEGER,PUBLIC :: m, mes, wal
    REAL(8),PUBLIC :: a_mag, a_eng, hlong, hfield!,der_locE_r, der_locE_s, der_locE_rs 
    REAL(8),DIMENSION(4),PUBLIC :: a, wa
    REAL(8),PUBLIC :: dt !, beta_r, beta_s, Jrs
    LOGICAL,PUBLIC :: key 

    REAL(8),  PUBLIC ::  Jo  
    PUBLIC :: epot,  prob_accept,ediff, metropolis_real, metropolis_hidden
	 

    PUBLIC :: initialize ,vmc, histogram, output,is_weight,Jperp,sgd

    INTEGER, ALLOCATABLE, DIMENSION(:)  :: SEED

    integer :: seed1                                                                                                  
    integer(kind=4) errcode
    REAL(8), EXTERNAL :: ran2  

    REAL(kind=4) :: aran = 0.e0,bran=1.e0
    REAL(kind=4), DIMENSION(:),ALLOCATABLE :: ranv
    REAL(kind=4), DIMENSION(:),ALLOCATABLE :: pran,rnd1, rnd2
    REAL(kind=4), DIMENSION(:),ALLOCATABLE :: ranvtot
    INTEGER(kind=4), DIMENSION(:), ALLOCATABLE :: nbinv,indmove
    INTEGER(KIND=4) :: isone4, Lxpo4,Lx4
    INTEGER :: isone

    INTEGER, PARAMETER :: nnearest = 2
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: isnear

    REAL(8), DIMENSION(:,:), ALLOCATABLE ::sconftime
    INTEGER :: itaumax = 100,igroup=100
    REAL(8), DIMENSION(:), ALLOCATABLE :: Ctau
    INTEGER :: imeasured, idum
    REAL(8), PUBLIC :: avg_clas,ls_avg_clas,rs_avg_clas
    REAL(8), PUBLIC :: ls_avg_mcla,rs_avg_mcla,ls_avg_eloc,rs_avg_eloc

    real(8), public, dimension(3) :: lambda, der_lambda

    real(8),public :: beta_r,beta_s,Jrs 


    CONTAINS

    ! ********************************************************************************************
    SUBROUTINE initialize()
        !REAL(8), DIMENSION(4), INTENT(INOUT) :: a
        INTEGER :: iaux, iaux2
        Real(8) :: lEo, rEo  

        ! Boundaries of the histogram
        smin = -1.d0
        smax =  1.d0

        ! Width of a bin in the histogram
        ds = (smax-smin)/DBLE(N2)

        ! Parameters initialization

        OPEN(UNIT=17,FILE='parameters.txt', STATUS='unknown')
        READ(17,*)
        READ(17,*) Jo,hfield,hlong,nstep1,nstep2,nwalk,mu,Lx,seed1 ,beta_r,beta_s,Jrs
        CLOSE(17)

        isone = 1
        isone4 = 1
        Lxpo4 = Lx+1
        Lx4 = Lx

        ncorr = nstep2/4

        IF (ALLOCATED(ranv))     DEALLOCATE(ranv)
        IF (ALLOCATED(pran))     DEALLOCATE(pran)
        IF (ALLOCATED(rnd1))     DEALLOCATE(rnd1)
        IF (ALLOCATED(rnd2))     DEALLOCATE(rnd2)
        IF (ALLOCATED(ranvtot))  DEALLOCATE(ranvtot)
        IF (ALLOCATED(spin))     DEALLOCATE(spin)
        IF (ALLOCATED(indmove))  DEALLOCATE(indmove)
        IF (ALLOCATED(spin_old)) DEALLOCATE(spin_old)
        IF (ALLOCATED(nbinv))    DEALLOCATE(nbinv)
        IF (ALLOCATED(lspin))     DEALLOCATE(lspin)
        IF (ALLOCATED(rspin))     DEALLOCATE(rspin)

        ALLOCATE(ranv(2*N1))
        ALLOCATE(pran(2*N1))
        ALLOCATE(spin(Lx,N1))
        ALLOCATE(indmove(N1))
        ALLOCATE(ranvtot(3*Lx*N1))
        ALLOCATE(spin_old(Lx,N1))
        ALLOCATE(nbinv(2*N1))
        ALLOCATE(lspin(Lx,N1))
        ALLOCATE(rspin(Lx,N1))
        ALLOCATE(isnear(nnearest,Lx))

        ALLOCATE(rnd1(nstep1))
        ALLOCATE(rnd2(nstep2))



        DO iaux = 1, Lx
            iaux2 = iaux+1
            IF (iaux2 .GT. Lx) iaux2 = 1
            isnear(2, iaux) = iaux2
            iaux2 = iaux-1
            IF (iaux2 .LT. 1) iaux2 = Lx
            isnear(1, iaux) = iaux2
        END DO

        imeasured = 0
        PRINT*, "The strength of nn spins interaction is:", Jo
        PRINT*, "The strength of the transverse field is:", hfield
        PRINT*, "The strength of the longitudinal field is:", hlong
        PRINT*, "Number of steps for MCs:", nstep1
        PRINT*, "Number of SGD steps:", nstep2
        PRINT*, "Target number of walkers:", nwalk
        PRINT*, "Length of the spin chain:", Lx

        Nspins = Lx*Ly

        Maxsons = nwalk
        nwalk1 = nwalk
	  
        lambda(1) = beta_r
        lambda(2) = beta_s
        lambda(3) = Jrs

        idum = seed1

	    !  Magnetization and energy initialization
        mag = 0.d0
        e = 0.d0
        En = 0.d0
        Eo = 0.d0
        lEo = 0.d0
        rEo = 0.d0

	    ! Initialize the matrix containing some important physical values 
        !a = 0.d0
        wa = 0.d0

        ! Initialize the arrays for histogram 
        ist = 0.d0

        count=0

    END SUBROUTINE initialize

    subroutine vmc(beta_r,beta_s,Jrs,energy,energy_err,derivative)
        REAL(8), INTENT(IN) :: beta_r,beta_s,Jrs
        real(8), intent(out) :: energy,energy_err
        real(8), intent(out), dimension(3):: derivative 
        real(8) :: eebavg,eebavg2,var_E
        real(8) :: rn
        real(8), dimension(3) :: der,der_var_E
        integer :: iibavg, it

        eebavg=0.0
        eebavg2=0.0
        der = 0.0
        iibavg = 0
       
        DO it = 1, nstep1
            rn = ran2(idum)
          
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
        INTEGER :: stobemoved_l,stobemoved_r
        REAL(8) :: ls_eloc,rs_eloc ,avg_clas,ls_avg_clas,rs_avg_clas
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
            stobemoved_l = INT((ran2(idum) * Nspins) + 1.d0) 
            lspin(stobemoved_l,iwalk) = -lspin(stobemoved_l,iwalk)
            deltaE_l = ediff(lspin(1:Nspins,iwalk),stobemoved_l)
            ds_l     = 2.d0*Jrs*spin(stobemoved_l,iwalk)*lspin(stobemoved_l,iwalk)
            prob_l   = exp(-beta_s*deltaE_l+ds_l)

            rn = ran2(idum)

            IF( prob_l .GE. 1.D0 )THEN
                Eo_l(iwalk) = Eo_l(iwalk) + deltaE_l
            ELSE IF(DBLE(rn) .LT. prob_l)THEN
                Eo_l(iwalk) = Eo_l(iwalk) + deltaE_l
            ELSE IF(DBLE(rn) .GE. prob_l)THEN
                lspin(stobemoved_l,iwalk) = -lspin(stobemoved_l,iwalk)
            END IF

            stobemoved_r = INT((ran2(idum) * Nspins) + 1.d0)
            rspin(stobemoved_r,iwalk) = -rspin(stobemoved_r,iwalk)
            deltaE_r = ediff(rspin(1:Nspins,iwalk),stobemoved_r)
            ds_r     = 2.d0*Jrs*spin(stobemoved_r,iwalk)*rspin(stobemoved_r,iwalk)
            prob_r   = exp(-beta_s*deltaE_r+ds_r)

            rn = ran2(idum)
            IF(prob_r .GE. 1.D0 )THEN
                Eo_r(iwalk) = Eo_r(iwalk) + deltaE_r
            ELSE IF(DBLE(rn) .LT. prob_r)THEN
                Eo_r(iwalk) = Eo_r(iwalk) + deltaE_r
            ELSE IF(DBLE(rn) .GE. prob_r)THEN
                rspin(stobemoved_r,iwalk) = -rspin(stobemoved_r,iwalk)
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

    END DO


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
    DO iwalk = 1, nwalk
               
      imoveact = INT((ran2(idum) * Nspins) + 1.d0)

      spin(imoveact,iwalk) = -spin(imoveact,iwalk)

     deltaE  = ediff(spin(1:Nspins,iwalk),imoveact)    
     dshadow = 2.d0*Jrs*spin(imoveact,iwalk)*(lspin(imoveact,iwalk)+rspin(imoveact,iwalk))
 
      prob   = exp(-2.d0*beta_r*deltaE+dshadow) 
     
      prn = ran2(idum) 
      IF(prob .GE. 1.D0 )THEN
!        mag(iwalk)  = mag(iwalk) + 2.d0*spin(imoveact,iwalk) 
        Eo(iwalk) = Eo(iwalk) + deltaE
        count = count + 1      
      ELSE IF(DBLE(prn) .LT. prob)THEN 
!        mag(iwalk)  = mag(iwalk) + 2.d0*spin(imoveact,iwalk)
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
      
   End do

   
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


!********************************************************************************************
! This subroutine performs VMC
!********************************************************************************************


! ******************************************************************************************
! This subroutine performs stochastic gradient descend
! ****************************************************************************************

  subroutine sgd(beta_r,beta_s,Jrs,energy,energy_err,derivative,mu_t)
   REAL(8), INTENT(inout) :: beta_s,beta_r,Jrs
   real(8), intent(inout) :: energy,energy_err
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

  DO is = 1, Nspins
    spin(is,iw) = -spin(is,iw)
             dE = ediff(spin(1:Nspins,iw),is) 
            lds = spin(is,iw)*lspin(is,iw)
            rds = spin(is,iw)*rspin(is,iw) 
        lEcl_rs = lEcl_rs + lds  
        rEcl_rs = rEcl_rs + rds
           lwt  = lwt + exp(-beta_r*dE+2.d0*Jrs*lds)
           rwt  = rwt + exp(-beta_r*dE+2.d0*Jrs*rds)
    spin(is,iw) = -spin(is,iw)

  END DO
 
  
  END SUBROUTINE is_weight

  ! ********************************************************************************************
! This subroutines generates an histogram for the ground state distribution of
! diffusers
! ********************************************************************************************

    SUBROUTINE histogram(mes,key)

        INTEGER, INTENT(IN)    :: mes
        LOGICAL, INTENT(IN)    :: key
        REAL(8)                :: mps    ! magnetization per spin
        Integer :: idx

        ist = 0.d0
        igg = igen
        wi  = 1.d0/DBLE(igg)

        DO iwalk = 1, igg
            mps = mag(iwalk)/DBLE(Nspins)
            IF( mps>smin .AND. mps<smax ) THEN
                idx = INT((mps-smin)/ds)+1
                ist(idx) = ist(idx)+ wi/ds
            END IF
        END DO


        IF(key .EQV. .TRUE.) THEN

        IF(MOD(mes,nstep2)==0) THEN

        OPEN(UNIT=17,FILE='profile.txt', STATUS='unknown', &
        ACTION='write')

        DO j = 1, N2
            WRITE(17,'(2F12.7)') smin+ds*DBLE(j-1)+0.5d0*ds,ist(j)/(DBLE(nstep2)/DBLE(ncorr))
        END DO

        WRITE(17,*)
        WRITE(17,*)
        CLOSE(17)

        END IF

        END IF

    END SUBROUTINE histogram

! ********************************************************************************************
! This subroutine computes the average magnetization and energy 
! ********************************************************************************************

 !   SUBROUTINE data(iwalk,avg_clas,ls_avg_clas,rs_avg_clas,ls_avg_mcla, &
 !     rs_avg_mcla,ls_avg_eloc,rs_avg_eloc,ecum1,ecum2,ecum3,scum1,scum2,rscum1,rscum2)

  !    integer, intent(in) :: iwalk
      !real(8) ::  lweight,rweight,lEcl_rs,rEcl_rs, ls_eloc, rs_eloc
  !    real(8), intent(inout) :: avg_clas, ls_avg_clas, rs_avg_clas
  !    real(8), intent(inout) :: ls_avg_mcla, rs_avg_mcla, ls_avg_eloc, rs_avg_eloc
  !    real(8), intent(inout) :: ecum1, ecum2, ecum3, scum1, scum2, rscum1, rscum2

    !***** CUMULATE DATA ****************************************
!      CALL is_weight(iwalk,lweight,rweight,lEcl_rs,rEcl_rs)

!      ls_eloc = -hfield*lweight   + Eo(iwalk)
 !     rs_eloc = -hfield*rweight   + Eo(iwalk)


!         avg_clas  =    avg_clas  + Eo(iwalk)
!      ls_avg_clas  = ls_avg_clas  + Eo_l(iwalk)
!      rs_avg_clas  = rs_avg_clas  + Eo_r(iwalk)
     
!      ls_avg_mcla  = ls_avg_mcla  + lEcl_rs
!      rs_avg_mcla  = rs_avg_mcla  + rEcl_rs

!      ls_avg_eloc = ls_avg_eloc  + ls_eloc
!      rs_avg_eloc = rs_avg_eloc  + rs_eloc

      

!      ecum1 = ecum1 + 0.5d0*(ls_eloc+rs_eloc)            
                
 !     ecum2 = ecum2 + (Eo(iwalk))*rs_eloc
  !    ecum3 = ecum3 + (Eo(iwalk))*ls_eloc

!      scum1 = scum1  + (Eo_l(iwalk))*rs_eloc
 !     scum2 = scum2  + (Eo_r(iwalk))*ls_eloc

!     rscum1 = rscum1 + lEcl_rs*rs_eloc
 !    rscum2 = rscum2 + rEcl_rs*ls_eloc


!    END SUBROUTINE data
! ********************************************************************************************
! This write the physical quantities in a file
! ********************************************************************************************

   SUBROUTINE output(hfield)

    REAL(8), INTENT(IN) :: hfield
    REAL(8):: Eave, E2ave, Mave, M2ave, sigmaE, sigmaM, avg

    avg = (DBLE(nstep2)/DBLE(ncorr))

    OPEN(10,file='avge_out.dat',STATUS='unknown', &
       ACTION='write')
    OPEN(13,file='fluctuations_nw.dat',STATUS='unknown', &
      ACTION='write')
!    OPEN(12,file='mag_out.dat',STATUS='unknown', &
!      ACTION='write')

    Eave     = a(1)/avg                                ! Average energy per site      
    E2ave    = a(2)/avg
    Mave     = a(3)/avg                                ! Average magnetization per site
    M2ave    = a(4)/avg

!    b2       = b(1)/avg
!    b4       = b(3)/avg
!    be2      = b(2)/avg
!    be4      = b(4)/avg
!    binderc  = 1.d0 -(b4/b2**2)/3.d0                   ! Computation of the Binder cumulant

!    bb2      = SQRT((be2-(b2**2))/(avg-1.d0))           ! Error in computing b2
!    bb4      = SQRT((be4-(b4**2))/(avg-1.d0))           ! Error in computing b4

    sigmaE   = SQRT((E2ave-(Eave**2))/(avg-1.d0))      ! Energy error bar
    sigmaM   = SQRT((M2ave-(Mave**2))/(avg-1.d0))      ! Magnetization error bar 
!    sigmaB   = SQRT( bb4**2/(9.d0*b2**4) + (4.d0/9.d0)*(b4**2/b2**6)*bb2**2 ) ! Binder cumulant error (propagation)

    !sigmaB   = bb4/ABS(b4) + (2.d0/3.d0)*(bb2/ABS(b2)) ! Binder cumulant error (overestimation)

    ! Energy computed via weighted average
!    Etot     = wa(1)/wa(2)                                ! Average energy per site      
!    varE     = wa(3)/(dt*(wa(2))**2)
!    sigmaEt  = SQRT(varE)             
!    Wavg     = wa(2)/wa(4)

    WRITE(10,'(F12.3,3F15.8)') hfield/Jo, Eave, sigmaE !,Removed beta
    WRITE(13,'(F7.2,4F15.8)') hfield/Jo,  Mave, sigmaM ! Removed, beta 
!    WRITE(12,'(F7.2,4F15.8)') hfield/Jo, Mave, sigmaM, binderc,sigmaB

    CLOSE(10)
    CLOSE(13)
!    CLOSE(12)

    END SUBROUTINE output
!  *******************************************************************************************

! ********************************************************************************************
! This is the probability to accept a spin flip 
! ********************************************************************************************

      REAL(8)  FUNCTION prob_accept(hf) RESULT(Y)
      REAL(8), INTENT(IN) :: hf

      Y = (SINH(hf))/(EXP(hf))

      END FUNCTION prob_accept
! ********************************************************************************************
! This is the interaction between physical spins and shadow ones
! ********************************************************************************************

      REAL(8)  FUNCTION Jperp(hf) RESULT(Y)
      REAL(8), INTENT(IN) :: hf

      Y = -0.5d0*LOG(TANH(0.5d0*hf))

      END FUNCTION Jperp
! ********************************************************************************************
! This function computes the potential energy for one copy of the system (walker)
! ********************************************************************************************
      REAL(8)  FUNCTION epot(spin,iwalk) RESULT(Y)
      
      REAL(8), DIMENSION(Lx,N1) , INTENT(IN) :: spin
      INTEGER, INTENT(IN) :: iwalk 
      INTEGER :: i 
      REAL(8) ::  E
       E = 0.d0 

 !No Periodic Boundary Conditions
         DO i = 1, Lx-1
         E = E - Jo*(spin(i,iwalk))*(spin(i+1,iwalk)) !- hlong*spin(i,iwalk)
         END DO

         Y = E !- Jo*(spin(Lx,iwalk))*(spin(1,iwalk))  !- hlong*spin(Lx,iwalk)

      END FUNCTION epot

! ********************************************************************************************
! This function computes the diiference in potential energy for one copy of the system
! (walker)
! ********************************************************************************************
      REAL(8)  FUNCTION ediff(spinsave,imoveact) RESULT(Y)
      
      REAL(8), DIMENSION(Nspins) , INTENT(IN) :: spinsave
      INTEGER :: j,iact,jact
      REAL(8) ::  E
      INTEGER, INTENT(IN) :: imoveact

      E = 0.d0

      iact = imoveact

      IF(iact ==1 )then
        E = E - Jo*(spinsave(iact))*(spinsave(iact+1))
      ELSE IF (iact == Nspins) then
        E = E - Jo*(spinsave(iact))*(spinsave(iact-1))
      ELSE
        DO j = 1, nnearest
          jact = isnear(j,iact)
          E = E - Jo*(spinsave(iact))*(spinsave(jact)) !- hlong*spinsave(iact)
        END DO
      END if
          
      Y = 2.d0*E
      END FUNCTION ediff

 ! ********************************************************************************************
      END MODULE functions



