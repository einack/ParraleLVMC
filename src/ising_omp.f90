! ********************************************************************************************
! Main program for the 1D Quantum Ising Model
! ********************************************************************************************

PROGRAM  isingmodel

    ! Make use of a module
    Use omp_lib
    USE functions_omp
    
    Implicit none

    integer :: ipar
    real(8)    :: rn1, rn2
    REAL(8) :: ebavg,eebavg,eebavg2 
    INTEGER :: ecount, ecount2, ibavg,iibavg,isgd
    real(8) :: sigma,learning_rate
    real(8) :: energy, energy_err
    real(8), dimension(3):: der, alpha
    real(8) :: timeinit, timef, time1, time3, time4 


    call cpu_time(timeinit)

    OPEN(unit=7,File='results_omp/spins.dat',Status='unknown')
    OPEN(unit=8,File='results_omp/energy.dat',Status='unknown')
    OPEN(unit=9,File='results_omp/energy_opt.dat',Status='unknown')
    OPEN(unit=10,File='results_omp/optimized.dat',Status='unknown')
    OPEN(unit=11,File='results_omp/poptimized.dat',Status='unknown')


    ! Initialize quantities 
    CALL initialize()

    sigma = 0.5
	!**************** PARAMETERS INITIALIZATION ********************
    pi = 4.d0*atan(1.0d0)

    do ipar = 1, 3
        rn1 = rand()
        !print*, ""
	 	!print *, "Rand1  = ", rn1, " Seed: ", idum 
        rn2 = rand()
	 	!print *, "Rand2  = ", rn2, " seed: ", idum 
        !print*, ""

        !stop(": in 1st do loop")
        alpha(ipar) = dble( SQRT(-2.d0*(sigma**2)*LOG(1.d0-rn1))*sin(2*pi*rn2) )
    end do

    
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
    ! spin, lspin and rspin buffers are populated

    buff_size = Lx * nwalk * 3

    allocate(randnumbers_1d( buff_size ))
    randnumbers_1d = 0.0

    do ipar = 1, buff_size 
        randnumbers_1d(ipar) = rand() 
    end do

    loc_size_spins = 3

    !$omp parallel private(low_bound, mid_bound, up_bound)
    nthreads = OMP_GET_NUM_THREADS() 
    Write(*,*) "Total Num of threads : ", nthreads
    Write(*,*) "Total Num of Rand calls: ", cnt

    do iwalk = 1, nwalk
    !$omp do 
        do i = 1, Lx

            low_bound = ( ( iwalk - 1 ) * loc_size_spins  * Lx ) + ( i * loc_size_spins) - (loc_size_spins - 1) 
            mid_bound = low_bound + 1
            up_bound = low_bound + loc_size_spins - 1 

            if ( randnumbers_1d(low_bound) .LT. 0.5) then
                spin(i,iwalk) = 1.d0
            else
                spin(i,iwalk) = -1.d0
            end if

            if ( randnumbers_1d(mid_bound) .LT. 0.5) then
                lspin(i,iwalk) = 1.d0
            else
                lspin(i,iwalk) = -1.d0
            end if

            if ( randnumbers_1d(up_bound) .LT. 0.5) then
                rspin(i,iwalk) = 1.d0
            else
                rspin(i,iwalk) = -1.d0
            end if

        end do
    !$omp end do 
    
        !$omp barrier    

        !$omp single
        ! Compute Potential energies for each walker
        Eo(iwalk) = epot(spin(1:Nspins,iwalk))
        Eo_l(iwalk) = epot(lspin(1:Nspins,iwalk))
        Eo_r(iwalk) = epot(rspin(1:Nspins,iwalk))
        !$omp end single

    end do
    !$omp end parallel

    do i=1, 10
        write(7,'(5(f4.1,1X))') spin(i,:)    
    end do


    PRINT*, 'Begin optimization' 
    learning_rate = mu
    ecount  = 1
    ecount2 = 0

    call cpu_time(time1)


    call vmc(alpha(1),alpha(2),alpha(3), energy, energy_err, der)

    call cpu_time(time3)


    !write(*,*) "here here 0"
    !stop ("Stoped after vmc")  ! Code stops here

    do isgd = 1, nstep2
        call sgd(alpha(1),alpha(2),alpha(3), energy, energy_err, der, learning_rate)

        if(MOD(isgd,10)==0 .OR. isgd==1) THEN
            ecount2 = ecount2 + 1

            print*, "i:",isgd, 'E:',energy, 'learning rate',learning_rate
            WRITE(8,'(i10,2F22.8)')  isgd, energy, energy_err
            WRITE(10,'(i10,4E22.8)') isgd, der(1), der(2),der(3)
            WRITE(11,'(i10,4E22.8)') isgd, alpha(1),alpha(2),alpha(3)

        end if
    end do

    call cpu_time(time4)
    
  
    print*, "End optimization"

    write(7,*)"Num of times rand is called: ", cnt

    CLOSE(7)  
    CLOSE(8)  
    CLOSE(10)
    CLOSE(11)
    
    !write(*,*)

    call cpu_time(timef)

    deallocate(spin)
    deallocate(lspin)
    deallocate(rspin)
    deallocate(isnear)

    write(*,fmt=771)

    write(*,fmt=777)'1', 'Before vmc', time1-timeinit 
    write(*,fmt=777)'2', 'Duration of vmc', time3-time1 
    write(*,fmt=777)'3', 'Duration of sgd:', time4-time3 
    write(*,fmt=777)'4', 'Total Runtime: ', timef-timeinit 

!770 format (5(f4.1,1X))
771 format('# No',1x,'Function',1x,'Duration')    
777 format(1a,1x,a15,1x,f12.6, 'secs')
END PROGRAM isingmodel

! ********************************************************************************************
! ********************************************************************************************

