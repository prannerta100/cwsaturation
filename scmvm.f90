!         zdiag    : diagonal elements of Z
!         zmat     : real and imaginary elements of Z upper diagonal
!         jzmat(i) : index of 1st imaginary element of row i in zmat
!         kzmat(i) : index of 1st real element of row i in zmat
!         izmat(i) : Z column index of ith element of zmat

!         arguments:

!         x        : input vector for matrix-vector multiplication
!         y        : resultant vector

!       scalars :

!             ndim  : number of rows in matrix

!       Local Variables:
!       ---------------

!             accr  : real part of dot product of a row of the
!                     matrix with the input vector
!             acci  : imaginary part of dot product of a row of the
!                     matrix with the input vector

!**********************************************************************

    subroutine scmvm(x,y,ndim, zmat, zdiag, izmat, jzmat, kzmat)

    include 'limits.inc'
!    include 'rndoff.inc'
!    include 'eprmat.inc'

    integer :: ndim, izmat(MXEL), jzmat(MXDIM+1), kzmat(MXDIM+1)
!    double precision :: x,y
!    complex*16 :: x(MXDIM),y(MXDIM)
    complex*16 :: x(*),y(*)
    double precision :: zmat(MXEL), zdiag(2,MXDIM)
!    dimension x(2,mxdim),y(2,mxdim)

    integer :: j,k,m,n,n1
    double precision :: accr,acci

!######################################################################

!----------------------------------------------------------------------
!     do diagonal elements first
!----------------------------------------------------------------------

    do 100 n=1,ndim
!        y(1,n)=zdiag(1,n)*x(1,n)-zdiag(2,n)*x(2,n), y(2,n)=zdiag(1,n)*x(2,n)+zdiag(2,n)*x(1,n)
        y(n)=cmplx(zdiag(1,n)*real(x(n))-zdiag(2,n)*aimag(x(n)), zdiag(1,n)*aimag(x(n))+zdiag(2,n)*real(x(n)))
    100 END DO

!----------------------------------------------------------------------
!     loop over rows (columns) of matrix for off-diagonal elements
!----------------------------------------------------------------------

    do 200 n=1,ndim
        n1=n+1
    
        accr=0.0D0
        acci=0.0D0
    
    !       imaginary matrix elements
    
        if (jzmat(n) /= jzmat(n1) ) then
            do 210 j=jzmat(n),jzmat(n1)-1
                m=izmat(j)
!                acci=acci+zmat(j)*x(1,m)
                acci=acci+zmat(j)*real(x(m))
!                y(2,m)=y(2,m)+zmat(j)*x(1,n)
!                accr=accr-zmat(j)*x(2,m)
                accr=accr-zmat(j)*aimag(x(m))
!                y(1,m)=y(1,m)-zmat(j)*x(2,n)
                y(m)=y(m)+cmplx(-zmat(j)*aimag(x(n)),zmat(j)*real(x(n)))
            210 END DO
        endif
    
    !       real matrix elements
    
        if (kzmat(n) /= kzmat(n1)) then
            do 220 k=kzmat(n),kzmat(n1)-1
                j = MXEL-k+1
                m=izmat(j)
                accr=accr+zmat(j)*real(x(m))
!                y(1,m)=y(1,m)+zmat(j)*x(1,n)
                acci=acci+zmat(j)*aimag(x(m))
!                y(2,m)=y(2,m)+zmat(j)*x(2,n)
                y(m)=y(m)+cmplx(zmat(j)*real(x(n)),zmat(j)*aimag(x(n)))
            220 END DO
        endif
    
!        y(1,n)=y(1,n)+accr
    ! djs        if (abs(y(1,n)).lt.rndoff) y(1,n)=0.0D0
!        y(2,n)=y(2,n)+acci
         y(n)=y(n)+cmplx(accr,acci)       
    ! djs        if (abs(y(2,n)).lt.rndoff) y(2,n)=0.0D0
    
    200 END DO

    return
    end subroutine scmvm
