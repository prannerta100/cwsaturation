!      program main
    subroutine sample_z(ndim, zmat, zdiag, izmat, jzmat, kzmat, w)
!f2py intent(in) ndim
!f2py intent(in) zmat
!f2py intent(in) zdiag
!f2py intent(in) izmat
!f2py intent(in) jzmat
!f2py intent(in) kzmat
!f2py intent(out) w
!*********************************************************************72

!c SAMPLE_Z illustrates the use of ZGEXPV and ZHEXPV.

!  Discussion:

!    Hermitian problem (Example 6.2 in the Expokit report).
   
    implicit none
    include 'limits.inc'
    external scmvm !zgcoov, zgcrsv, zgccsv

    double precision :: tic, tac, clock

!  matrix data.
!  BEWARE: these values must match those in zgmatv.f

!    integer :: n, nz, nmax, nzmax
!    parameter( nmax = 5500, nzmax = 50000 )
    integer :: ndim, izmat(MXEL), jzmat(MXDIM+1), kzmat(MXDIM+1)
!    double precision :: x,y
    complex*16 :: x(MXDIM),y(MXDIM)
    double precision :: zmat(MXEL), zdiag(2,MXDIM)

!    integer :: ia(nzmax), ja(nzmax)
!    complex*16 a(nzmax)
!    common /CMAT/ a, ia, ja, nz, n

!  arguments variables.

    integer :: m, mmax, lwsp, liwsp, ndim
    parameter( mmax = 50 )
    parameter( lwsp = nmax*(mmax+2)+5*(mmax+2)**2+7, liwsp = mmax+2 )
    integer :: iwsp(liwsp)
    double precision :: t, tol, anorm, s1, s2
    complex*16 v(nmax), w(nmax), wsp(lwsp)

    integer :: i, j, nnz, itrace, iflag, iseed(4)
    complex*16 ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )

    double precision :: DLARAN
    intrinsic ABS, CMPLX, CONJG, DBLE

!  Executable statements.

!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'SAMPLE_Z'

!  Load a symmetric pattern.

!    n = nmax
!    nz = nzmax/2
!    call getpat ( 'bcspwr10.txt', n, nz, ia, ja )

!  For the purpose of the experiments, expand to COOrdinates.

!    nnz = nz
!    do j = n,1,-1
!        do i = 1,ja(j+1)-ja(j)
!            ja(nnz) = j
!            nnz = nnz-1
!        end do
!    end do

!  Fill-in an Hermitian matrix -- the conjugate part is included.

!    iseed(1) = 0
!    iseed(2) = 0
!    iseed(3) = 0
!    iseed(4) = 5
!    nnz = nz
!    do i = 1,nz
!        if ( ia(i) /= ja(i) ) then
!            s1 = 10.0d0*DLARAN( iseed ) - 5.0d0
!            s2 = 10.0d0*DLARAN( iseed ) - 5.0d0
!            a(i) = CMPLX( s1,s2 )
!            nnz = nnz + 1
!            a(nnz) = CONJG( a(i) )
!            ia(nnz) = ja(i)
!            ja(nnz) = ia(i)
!        else
!            s1 = 10.0d0*DLARAN( iseed ) - 5.0d0
!            a(i) = CMPLX( s1,0.0d0 )
!        endif
!    end do
!    nz = nnz

!  Compute the infinity norm of A.

!    do i = 1,n
!        wsp(i) = ZERO
!    end do
!    do i = 1,nz
!        wsp(ia(i)) = wsp(ia(i)) + ABS( a(i) )
!    end do
!    anorm = wsp(1)
!    do i = 2,n
!        if ( anorm < DBLE(wsp(i)) ) then
!            anorm =  wsp(i)
!        end if
!    end do

!  Convert from COO to CRS.

!      call zgcnvr( 'coo','crs','n', n,n, nz, ia, ja, a, iwsp )

!  Compute the infinity norm of A.

!      anorm = 0.0d0
!      do i = 1,n
!         s1 = 0.0d0
!         do j = ia(i),ia(i+1)-1
!            s1 = s1 + ABS( a(j) )
!         end do
!         if ( anorm.lt.tmp ) anorm = s1
!      end do

    write(UNIT=*,FMT='(A,E8.2)') '||A||_inf= ',anorm

!  The operand vector v is set to e_1 + e_n.

    do i = 1,n
        v(i) = ZERO
    end do
    v(1) = ONE
    v(n) = ONE

!  Set other input arguments.

    t = 1.0d0
    tol = 1.0d-5
    m = 30
    itrace = 0

!  Compute w = exp(t*A)v with ZGEXPV.

    tic = clock()
    call ZGEXPV( n, m, t,v,w, tol, anorm, &
    wsp,lwsp, iwsp,liwsp, zgcoov, itrace, iflag)
!     &             a, ia, ja, nz, n  )

    tac = clock()

    print 9001,'----------------------------------------------------'
    print 9001,'ZGEXPV has completed:'
    print 9001,'----------------------------------------------------'
    print 9001,'w(1:10) ='
    do i = 1,10
        print*,w(i)
    end do

!  Display some statistics if desired.

    print 9001,'final report----------------------------------------'
    print 9002,'runtime   = ',tac-tic
    print 9002,'||A||_inf = ',anorm
    print 9003,'nz        =',nz
    print 9003,'n         =',n
    print 9003,'m         =',m
    print 9003,'itrace    =',itrace
    print 9003,'iflag     =',iflag
    print 9003,'ibrkflag  =',iwsp(6)
    print 9003,'mbrkdwn   =',iwsp(7)
    print 9003,'nstep     =',iwsp(4)
    print 9003,'nreject   =',iwsp(5)
    print 9003,'nmult     =',iwsp(1)
    print 9003,'nexph     =',iwsp(2)
    print 9003,'nscale    =',iwsp(3)

    print 9002,'tol       = ',tol
    print 9002,'t         = ',t
    print 9002,'tbrkdwn   = ',DBLE( wsp(7) )
    print 9002,'step_min  = ',DBLE( wsp(1) )
    print 9002,'step_max  = ',DBLE( wsp(2) )
    print 9002,'max_round = ',DBLE( wsp(3) )
    print 9002,'sum_round = ',DBLE( wsp(4) )
    print 9002,'max_error = ',DBLE( wsp(5) )
    print 9002,'sum_error = ',DBLE( wsp(6) )
    print 9002,'hump      = ',DBLE( wsp(9) )
    print 9002,'scale-norm= ',DBLE( wsp(10) )
    write(*,*) 'nina=',nina
          
    9001 format(A)
    9002 format(A,E8.2)
    9003 format(A,I9)
    end subroutine sample_z
