
!    subroutine rand()
!    INTEGER :: IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
!    DOUBLE PRECISION :: AM,EPS,RNMX, rand_num
!    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,&
!    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
!    NTAB=32,NDIV=1+IMM1/NTAB,EPS=3.d-16,RNMX=1.d0-EPS )
!    INTEGER :: idum2,j,k,iv(NTAB),iy
!    SAVE :: iv,iy,idum2
!    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
!
!
!    if ( idum_seed .le. 0 ) then
!        idum_seed=max(-idum_seed,1)
!        idum2=idum_seed
!
!        do 11 j = NTAB+8, 1, -1
!
!            k=idum_seed/IQ1
!            idum_seed=IA1*(idum_seed-k*IQ1)-k*IR1
!            if (idum_seed.lt.0) idum_seed=idum_seed+IM1
!            if (j.le.NTAB) iv(j)=idum_seed
!
!11      continue
!            
!        iy=iv(1)
!
!    endif
!
!    k=idum_seed/IQ1
!    idum_seed=IA1*(idum_seed-k*IQ1)-k*IR1
!
!    if (idum_seed.lt.0) idum_seed=idum_seed+IM1
!
!    k=idum2/IQ2
!    idum2=IA2*(idum2-k*IQ2)-k*IR2
!
!    if (idum2.lt.0) idum2=idum2+IM2
!
!    j=1+iy/NDIV
!    iy=iv(j)-idum2
!    iv(j)=idum_seed
!
!    if(iy.lt.1)iy=iy+IMM1
!
!    rand_num=min(AM*iy,RNMX)
!
!    end subroutine rand
!
!
!
!!!CC{{{  Functions 
    FUNCTION rand()
    INTEGER IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    DOUBLE PRECISION rand,AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1.d0/IM1,IMM1=IM1-1,&
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
    NTAB=32,NDIV=1+IMM1/NTAB,EPS=3.d-16,RNMX=1.d0-EPS)
    INTEGER idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/

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

    return
    end function
!C  (C) Copr. 1986-92 Numerical Recipes Software (9`3j32150.
!CC}}}

