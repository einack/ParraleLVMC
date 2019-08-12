! ********************************************************************************************
! Main program for the 1D Quantum Ising Model
! ********************************************************************************************

PROGRAM  isingmodel
    Use mpi            
    ! Make use of a module
    USE functions

    Implicit none
    integer :: cnt, cn1, nbin, ipar
    REAL,dimension(3)    :: rn
    real    :: rn1, rn2
    REAL(8) :: eecum, ecum1, eecum2, ebavg,eebavg,eebavg2 
    INTEGER :: ecount, ecount2, cntmu, ibavg,iibavg,isgd
    real(8) :: sigma,learning_rate
    real(8) :: energy,energy_err
    real(8), dimension(3):: der, alpha

    integer :: ierr
    integer :: numtasks, rank
    integer :: loc_frame, rest, offset


    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

    !#ifdef __DEBUG
    write(*,*)'My rank is', rank
    write(*,*)
    !#endif

    if ( rank == 0 ) then 
        write(*,*)'My rank is', rank
        write(*,*)
        OPEN(unit=8,File='energy.dat',Status='unknown')
        OPEN(unit=9,File='energy_opt.dat',Status='unknown')
        OPEN(unit=10,File='optimized.dat',Status='unknown')
        OPEN(unit=11,File='poptimized.dat',Status='unknown')
    endif

    ! Initialize quantities by root rank 

    CALL initialize(a)

    write(*,*) "I am: ", rank ," ::::: ", a
    write(*,*)

    sigma = 0.5

    !**************** PARAMETERS INITIALIZATION ********************
    pi = 4.d0*atan(1.0d0)



    do ipar = 1, 3

        rn1 = ran2(idum)
        rn2 = ran2(idum)
        rn(ipar)  = SQRT(-2.d0*(sigma**2)*LOG(1.d0-rn1))*sin(2*pi*rn2)

    end do

    alpha(1) = dble(rn(1))
    alpha(2) = dble(rn(2))
    alpha(3) = dble(rn(3))

    ibavg  = 0
    iibavg = 0
    ebavg  = 0.d0
    eebavg = 0.d0
    eebavg2= 0.d0

    spin = 0.d0
    lspin = 0.d0
    rspin = 0.d0

    print*, 'beta_r', alpha(1)
    print*, 'beta_s',alpha(2)
    print*, 'Jrs',alpha(3)
 
    !**************   WALKERS INITIALIZATION ************************
    DO iwalk = 1, nwalk

        DO i = 1, Lx
            IF (ran2(idum) .LT. 0.5) THEN
               spin(i,iwalk) = 1.d0
            ELSE
               spin(i,iwalk) = -1.d0
            END IF

            IF (ran2(idum) .LT. 0.5) THEN
               lspin(i,iwalk) = 1.d0
            ELSE
               lspin(i,iwalk) = -1.d0
            END IF

            IF (ran2(idum) .LT. 0.5) THEN
               rspin(i,iwalk) = 1.d0
            ELSE
               rspin(i,iwalk) = -1.d0
            END IF

        END DO
     
        Eo(iwalk) = epot(spin,iwalk)

        Eo_l(iwalk) = epot(lspin,iwalk)
        Eo_r(iwalk) = epot(rspin,iwalk)

    END DO

    PRINT*, 'Begin optimization' 
    learning_rate = mu
 
    ecount  = 1
    ecount2 = 0
    eecum   =  ecum1/(DBLE(Nspins*nwalk))
    eecum2  = (ecum1/(DBLE(Nspins*nwalk)))**2.d0
 
 
    call vmc(alpha(1),alpha(2),alpha(3),energy,energy_err,der)

    DO isgd = 1, nstep2
        cnt = isgd
        call sgd(alpha(1),alpha(2),alpha(3),energy,energy_err,der,cnt,learning_rate)
        IF(MOD(isgd,10)==0 .OR. isgd==1) THEN
        ecount2 = ecount2 + 1

        print*, "i:",isgd, 'E:',energy, 'learning rate',learning_rate
        WRITE(8,'(i10,2F22.8)')  isgd, energy, energy_err
        WRITE(10,'(i10,4E22.8)') isgd, der(1), der(2),der(3)
        WRITE(11,'(i10,4E22.8)') isgd, alpha(1),alpha(2),alpha(3)

        END IF
    END DO
  
    print*, "End optimization"

    CLOSE(8)  
    CLOSE(10)
    CLOSE(11)

    call MPI_FINALIZE(ierr)

  END PROGRAM isingmodel

! ********************************************************************************************
! ********************************************************************************************

