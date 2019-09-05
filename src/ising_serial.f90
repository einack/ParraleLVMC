! ********************************************************************************************
! Main program for the 1D Quantum Ising Model
! ********************************************************************************************

PROGRAM  isingmodel

    ! Make use of a module
    USE functions
    
    Implicit none

    integer :: ipar
    real(8)    :: rn1, rn2
    REAL(8) :: eecum, eecum2, ecum1, ebavg,eebavg,eebavg2 
    INTEGER :: ecount, ecount2, ibavg,iibavg,isgd
    real(8) :: sigma,learning_rate
    real(8) :: energy,energy_err
    real(8), dimension(3):: der, alpha
    real(8) :: timeinit, timef, time1, time2, time3, time4 


    call cpu_time(timeinit)

    OPEN(unit=8,File='energy.dat',Status='unknown')
    OPEN(unit=9,File='energy_opt.dat',Status='unknown')
    OPEN(unit=10,File='optimized.dat',Status='unknown')
    OPEN(unit=11,File='poptimized.dat',Status='unknown')


    ! Initialize quantities 
    CALL initialize()

    sigma = 0.5
	!**************** PARAMETERS INITIALIZATION ********************
    pi = 4.d0*atan(1.0d0)
    do ipar = 1, 3
        rn1 = ran2(idum)
        rn2 = ran2(idum)
        alpha(ipar) = SQRT(-2.d0*(sigma**2)*LOG(1.d0-rn1))*sin(2*pi*rn2)
    end do


    !alpha(1) = dble(rn(1))
    !alpha(2) = dble(rn(2))
    !alpha(3) = dble(rn(3))

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

    do iwalk = 1, nwalk
        do i = 1, Lx
            if (ran2(idum) .LT. 0.5) then
                spin(i,iwalk) = 1.d0
            else
                spin(i,iwalk) = -1.d0
            end if

            if (ran2(idum) .LT. 0.5) then
                lspin(i,iwalk) = 1.d0
            else
                lspin(i,iwalk) = -1.d0
            end if

            if (ran2(idum) .LT. 0.5) then
                rspin(i,iwalk) = 1.d0
            else
                rspin(i,iwalk) = -1.d0
            end if

        end do

        Eo(iwalk) = epot(spin,iwalk)
        Eo_l(iwalk) = epot(lspin,iwalk)
        Eo_r(iwalk) = epot(rspin,iwalk)

    end do

    PRINT*, 'Begin optimization' 
    learning_rate = mu
    ecount  = 1
    ecount2 = 0
    eecum   =  ecum1/(DBLE(Nspins*nwalk))
    eecum2  = (ecum1/(DBLE(Nspins*nwalk)))**2.d0

    call cpu_time(time1)

    call cpu_time(time2)

    call vmc(alpha(1),alpha(2),alpha(3),energy,energy_err,der)

    call cpu_time(time3)


    !write(*,*) "here here 0"
    !write(*,*) "Debug, alph", alpha
    !stop 0  ! Code stops here

    do isgd = 1, nstep2
        call sgd(alpha(1),alpha(2),alpha(3),energy,energy_err,der,learning_rate)
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

    CLOSE(8)  
    CLOSE(10)
    CLOSE(11)
    
    call cpu_time(timef)
    write(*,fmt=771)

    write(*,fmt=777)'1', 'Before vmc', time1-timeinit 
    write(*,fmt=777)'2', 'Duration of vmc', time3-time2 
    write(*,fmt=777)'3', 'Duration of sgd:', time4-time3 
    write(*,fmt=777)'4', 'Total Runtime: ', timef-timeinit 

771 format('# No',1x,'Function',1x,'Duration')    
777 format(1a,1x,a15,1x,f12.6, 'secs')
END PROGRAM isingmodel

! ********************************************************************************************
! ********************************************************************************************

