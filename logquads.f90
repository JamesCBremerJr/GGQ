!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This program constructs a collection of quadrature rules which allow for the 
!  discretization of integral equations of the form
!
!    T[\sigma](x) =  \int    K(x,y) \sigma(y) dy,
!                       \Gamma
!
!  where \Gamma is a curve in the plane and K(x,y) behaves as 
!
!    K(x,y) ~ log (|x-y| ) S_1(x,y) + S_2(x,y)
!
!  as x --> y, with S_1 and S_2 smooth.
!
!  The resulting quadrature rules are written to a file on the disk called logsing.f90.
!
!  It is advisable to execute this code using extended precision arithmetic.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module     logquads_funs

use utils
use legendre
use legepw
use chebyshev_quad
use gaussian_quad
use makequad
use iso_c_binding


double precision, allocatable                :: x0

contains


subroutine logquad_sing(norder,x,nquad,xs,whts)
implicit double precision (a-h,o-z)
double precision                            :: xs(:), whts(:)
!
!  Construct a quadrature rule for evaluating integrals of the form
!
!         1 
!    \int    log|x-y| P(y) + Q(y)   dy, where
!        -1
!        
!  P and Q are polynonomials of degree less than or equal to norder and x
!  is a specified point in [-1,1].
!

double precision, allocatable :: xsquad(:),whtsquad(:)
type(c_ptr)                   :: userptr

x0     = x
npols  = norder+1
ifsing = 0

a      = -1.0d0-x0
b      =  1.0d0-x0
nfuns  = npols*2

call ggquad(ifsing,a,b,nfuns,funlog1,userptr,nquad,xsquad,whtsquad)
xsquad        = xsquad + x0
xs(1:nquad)   = xsquad
whts(1:nquad) = whtsquad

end subroutine


subroutine funlog1(nfuns,n,xs,vals,userptr)
implicit double precision (a-h,o-z)
double precision           :: xs(n)
double precision           :: vals(:,:)
type(c_ptr)                :: userptr

double precision, allocatable :: pols(:)

npols = nfuns/2
allocate(pols(npols+10))

do i=1,n
x = xs(i)
call leges(npols,x+x0,pols)

vals(i,1:npols)          = pols(1:npols) * log(abs(x))
vals(i,npols+1:2*npols)  = pols(1:npols)
end do

end subroutine


subroutine logquad_near(norder,nquad,xs,whts)
implicit double precision (a-h,o-z)
double precision, allocatable          :: xs(:), whts(:)
!
!  Construct a quadrature rule for evaluating integrals of the form
!
!         1 
!    \int    log|x-y| P(y) + Q(y)   dy, where
!        -1
!        
!  P and Q are polynonomials of degree less than or equal to norder and x
!  is a point in one of the interval [-2,-1] or [1,2].
!

double precision, allocatable         :: ab(:,:)

double precision, allocatable         :: xsquad(:),whtsquad(:)
double precision, allocatable         :: xslege(:),whtslege(:)
double precision, pointer             :: xs0(:,:)
type(c_ptr)                           :: userptr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nlege = norder+1
! call legendre_quad(nlege,xslege,whtslege)

! nxs = 6*nlege
! allocate(xs0(2,nxs))

! idx = 0
! a0  = -5.0d0
! b0  = -1.0d0
! do i=1,nlege
! idx        = idx+1
! xs0(1,idx) = (b0-a0)/2*xslege(i) + (b0+a0)/2
! xs0(2,idx) = 0

! end do

! a0  = -3.0d0
! b0  = -1.0d0
! do i=1,nlege
! idx        = idx+1
! xs0(1,idx) = (b0-a0)/2*xslege(i) + (b0+a0)/2
! xs0(2,idx) = 0
! end do

! a0  =-2.0d0
! b0  =-1.0d0
! do i=1,nlege
! idx        = idx+1
! xs0(1,idx) = (b0-a0)/2*xslege(i) + (b0+a0)/2
! xs0(2,idx) = 0
! end do

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! a0  = 1.0d0
! b0  = 2.0d0
! do i=1,nlege
! idx        = idx+1
! xs0(1,idx) = (b0-a0)/2*xslege(i) + (b0+a0)/2
! xs0(2,idx) = 0

! end do

! a0  = 1
! b0  = 3
! do i=1,nlege
! idx        = idx+1
! xs0(1,idx) = (b0-a0)/2*xslege(i) + (b0+a0)/2
! xs0(2,idx) = 0
! end do

! a0  = 1
! b0  = 5
! do i=1,nlege
! idx        = idx+1
! xs0(1,idx) = (b0-a0)/2*xslege(i) + (b0+a0)/2
! xs0(2,idx) = 0
! end do

! call prin2("xs0 = ",xs0)

!
!  Build the points in the target quadrature
!
! allocate(xs0(2,nlege*2))

! a0 = 1.0d0+eps
! b0 = 2.0d0
! do i=1,nlege
! xs0(1,i) = xslege(i)*(b0-a0)/2 + (b0+a0)/2
! xs0(2,i) = 0
! end do

! a0 = -2.0d0
! b0 = -1.0d0-eps
! do i=1,nlege
! xs0(1,i+nlege) = xslege(i)*(b0-a0)/2 + (b0+a0)/2
! xs0(2,i+nlege) = 0
! end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! nxs = 4*nlege
! allocate(xs0(2,4*nlege))

! idx = 0
! a0  = -eps
! b0  = eps

! do i=1,nlege
! idx        = idx+1
! xs0(1,idx) = 1+eps
! xs0(2,idx) = xslege(i)*(b0-a0)/2 + (a0+b0)/2
! end do

! do i=1,nlege
! idx        = idx+1
! xs0(1,idx) = 1-eps
! xs0(2,idx) = xslege(i)*(b0-a0)/2 + (a0+b0)/2
! end do

! do i=1,nlege
! idx        = idx+1
! xs0(1,idx) = xslege(i)*(b0-a0)/2 + (a0+b0)/2
! xs0(2,idx) = 1-eps
! end do

! do i=1,nlege
! idx        = idx+1
! xs0(1,idx) = xslege(i)*(b0-a0)/2 + (a0+b0)/2
! xs0(2,idx) = 1+eps
! end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nlege = 16
call legendre_quad(nlege,xslege,whtslege)

nints = 10
allocate(ab(2,nints))
do i=1,nints
ab(1,i) = 1 + 2.0d0**(-nints+i-1)
ab(2,i) = 1 + 2.0d0**(-nints+i)
end do

nxs = nlege*nints*2
idx = 0

allocate(xs0(2,nxs))
do int=1,nints
a = ab(1,int)
b = ab(2,int)
do i=1,nlege
x = (b-a)/2 * xslege(i) + (b+a)/2
idx = idx + 1
xs0(1,idx)  = x
xs0(2,idx)  = 0
idx = idx + 1
xs0(1,idx)  = -x
xs0(2,idx)  = 0
end do
end do

call prin2("in logquad_near, xs0 = ",xs0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nlege = norder+1
! call legendre_quad(nlege,xslege,whtslege)
! allocate(xs0(2,nlege))
! nxs = nlege
! a0 = -2
! b0 = -1
! do i=1,nlege
! xs0(1,i) = (b0-a0)/2*xslege(i) + (b0+a0)/2
! xs0(2,i) = 0
! end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

userptr = c_loc(xs0)
npols  = norder+1
ifsing = 0
a      = -1.0d0
b      =  1.0d0
nfuns  = nxs+1

call ggquad_prod(ifsing,a,b,nfuns,npols,funlog2,funpolys,userptr,userptr, &
  nquad,xs,whts)

call prin2("xs = ",xs)
call prin2("whts = ",whts)

end subroutine



subroutine funlog2(nfuns,n,xs,vals,userptr)
implicit double precision (a-h,o-z)
double precision           :: xs(n)
double precision           :: vals(:,:)
type(c_ptr)                :: userptr
double precision, pointer  :: xs0(:,:)

call c_f_pointer(userptr,xs0,[2,nfuns-1])

do i=1,n
x               = xs(i)

do j=1,nfuns-1
x0              = xs0(1,j)
y0              = xs0(2,j)
dd              = (x-x0)**2 + y0**2
vals(i,j)       = log(sqrt(dd))
end do

vals(i,nfuns)   = 1.0d0

end do


end subroutine


subroutine funpolys(nfuns,n,xs,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: xs(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
!
!  Evaluate the Legendre polynomials
!

do i=1,n
x  = xs(i)
call leges(nfuns,x,vals(i,:))
end do

end subroutine

end module logquads_funs


program logquads

use utils
use legendre
use legepw
use chebyshev_quad
use gaussian_quad
use makequad
use logquads_funs

implicit double precision (a-h,o-z)

integer, allocatable                   :: nquadssing(:)
double precision, allocatable          :: xssing(:,:), whtssing(:,:)
double precision, allocatable          :: xslege(:), whtslege(:)
double precision, allocatable          :: xs(:),whts(:)
double precision, allocatable          :: xsnear(:),whtsnear(:)
!character(len=14)                      :: name
!

!  Construct the quadrature rules
!



maxquad = 100
norder  = 19
nlege   = norder+1
call legendre_quad(nlege,xslege,whtslege)

allocate(nquadssing(nlege))
allocate(xssing(maxquad,nlege))
allocate(whtssing(maxquad,nlege))

do i=1,nlege
x = xslege(i)
call logquad_sing(norder,x,nquadssing(i),xssing(:,i),whtssing(:,i))
end do

call logquad_near(norder,nquadnear,xsnear,whtsnear)


!
!  Write them out to disk
!

iw = 20
open(iw,FILE="logsing.f90")

write(iw,"(A,I2.2,A)") "subroutine logquads",norder,"(nquads,xs,whts,nquadnear,xsnear,whtsnear)"
write(iw,"(A)") "implicit double precision (a-h,o-z)"
write(iw,"(A)") "double precision :: xs(:,:),whts(:,:)"
write(iw,"(A)") "double precision :: xsnear(:),whtsnear(:)"
write(iw,"(A)") "integer nquads(:)"
write(iw,"(A,I3)") "do i=1,",nlege

do i=1,nlege
write(iw,"(A,I3,A,I3)") "nquads(",i,")=",nquadssing(i)
write(iw,"(A,I3)") "do j=1,",nquadssing(i)
do j=1,nquadssing(i)
write(iw,"(A,I3,A,I3,A,D44.36)") "xs  (",j,",",i,") = ",xssing(j,i)
write(iw,"(A,I3,A,I3,A,D44.36)") "whts(",j,",",i,") = ",whtssing(j,i)

end do
write(iw,"(A)") "end do"
end do

write(iw,"(A,I3,A,I3)") "nquadnear = ",nquadnear
do j=1,nquadnear
write(iw,"(A,I3,A,D44.36)") "xsnear   (",j,") = ",xsnear(j)
write(iw,"(A,I3,A,D44.36)") "whtsnear (",j,") = ",whtsnear(j)
end do

write(iw,"(A)") "end do"
write(iw,"(A)") "end subroutine"




end program
