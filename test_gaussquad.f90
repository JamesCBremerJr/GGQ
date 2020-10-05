module test_gaussquad_funs
use iso_c_binding
use legendre

type     funalphas_data
integer                       :: nalphas
integer                       :: npols
double precision, allocatable :: alphas(:)
end type funalphas_data

double precision              :: dsub
contains

subroutine funalphas(nfuns,n,us,whts,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: us(n), whts(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
type(funalphas_data), pointer :: fundata


double precision              :: xs(n), dxs(n), pols(500)

xs   = us**dsub
dxs  = dsub * us**(dsub-1)

call c_f_pointer(userptr,fundata)

do j=1,nfuns
alpha     = fundata%alphas(j)
vals(:,j) = abs(xs)**alpha * dxs * sqrt(whts)
end do


end subroutine


subroutine disc_alphas(nalphas,alpha1,alpha2,fundata)
implicit double precision (a-h,o-z)
type(funalphas_data), pointer :: fundata

allocate(fundata)
allocate(fundata%alphas(nalphas))
fundata%nalphas = nalphas
do i=1,nalphas
dd                = (i-1.0d0)/(nalphas-1.0d0)
fundata%alphas(i) = alpha1 + (alpha2-alpha1) * dd
end do
end subroutine



end module


program test_gaussquad

use utils
use legendre
use linalg
use legepw
use chebyshev_quad
use gaussian_quad
use iso_c_binding
use test_gaussquad_funs


implicit double precision (a-h,o-z)

type(legepw_disc)             :: disc
type(c_ptr)                   :: userptr
type(funalphas_data), pointer :: fundata
double precision, allocatable :: xs0(:), whts0(:)
double precision, allocatable :: xs(:), whts(:)
double precision, allocatable :: us(:), uwhts(:)

integer, allocatable          :: ipivs(:)
double precision, allocatable :: vals(:,:), valsq(:,:), r(:,:), rnorms(:)
double precision, allocatable :: alphas(:), rinv(:,:)
double precision, pointer     :: alpha0

eps0    = epsilon(0.0d0)

epsdisc = 1.0d-13      ! precision for discretization
epsqr   = 1.0d-13      ! precision for the QR decomposition
epsnewt = 1.0d-7      ! precision for the Newton iterations

nlege   = 30
a       =-1.0d0
b       = 1.0d0
dsub    = 3

!
!  Adaptively construct a discretization scheme for the input functions |x|^alpha
!


nalphas = 100
nfuns   = nalphas
alpha1  =-0.50d0
alpha2  = 0.50d0
call disc_alphas(nalphas,alpha1,alpha2,fundata)


call prini("nfuns = ",nfuns)
call legepw_init(disc,nlege,a,b)
userptr = c_loc(fundata)

call elapsed(t1)
call legepw_adap(epsdisc,nfuns,funalphas,disc,userptr)
call elapsed(t2)

call legepw_quad(disc,nquad0,xs0,whts0)
call prin2("legepw_adap time = ",t2-t1)
call prini("nints = ",disc%nints)
call prini("nquad0 = ",nquad0)

!
!  Evaluate the input functions at the discretization nodes
!
allocate( vals(nquad0,nfuns) )
call elapsed(t1)
call funalphas(nfuns,nquad0,xs0,whts0,vals,userptr)
call elapsed(t2)
call prin2("funalphas time = ",t2-t1)


!
!  Orthonormalize the functions ...
!

call elapsed(t1)
call qrdecomp(epsqr,nquad0,nfuns,vals,krank,ipivs,valsq,r)
allocate(rnorms(krank))
do i=1,krank
rnorms(i) = abs(r(i,i))
end do
call elapsed(t2)
call prin2("qrdecomp time = ",t2-t1)
call prini("after qrdecomp, krank = ",krank)
call prin2("after qrdecomp, rnorms = ",rnorms)


!
!  Construct the Chebyshev quadrature  
!

call elapsed(t1)
call chebquad(disc,krank,valsq,nquad,us,uwhts)
call elapsed(t2)
call prin2("chebquad time  = ",t2-t1)
call prin2("chebquad us    = ",us)
call prin2("chebquad uwhts = ",uwhts)
call prina("")

!
!  Use Newton iterations to reduce the quadrature rule
!


call elapsed(t1)
iorder=1
call gaussquad(epsnewt,iorder,disc,krank,valsq,nquad,us,uwhts)
call elapsed(t2)
call prin2("gaussquad time = ",t2-t1)
call prini("guassquad nquad = ",nquad)
call prin2("gaussquad us = ",us)
call prin2("gaussquad uwhts = ",uwhts)


!
!  Account for the substituion
!

allocate(xs(nquad), whts(nquad) )
xs    = us**dsub
whts  = dsub * us**(dsub-1) * uwhts

call prina("")
call prini("krank          = ",krank)
call prini("final nquad    = ",nquad)
call prin2("final xs       = ",xs)
call prin2("final whts     = ",whts)






dmax = 0
do l=1,nalphas
alpha = fundata%alphas(l)

sum2 = 2/(1.0d0+alpha)
sum1 = 0

do i=1,nquad
x = xs(i)
wht = whts(i)
sum1 = sum1 +  abs(x)**alpha * wht
end do

dmax = max(abs(sum1-sum2),dmax)
end do

call prin2("maximum absolute error = ",dmax)

stop

!
!  Check the accuracy of the formula for nn random values of alpha
!  in the specified range by comparison with adaptive quadrature ---
!  this is slow but greatly reduces the possibility that the
!  initial discretization went awry
!

allocate(alpha0)
userptr = c_loc(alpha0)

nn = 1000
allocate(alphas(nn))
do i=1,nn
call random_number(dd)
alphas(i) = alpha1 + (alpha2-alpha1)*dd
end do
call quicksort(nn,alphas)

dmax = 0
do j=1,nn

alpha0 = alphas(j)
dsum1  = 0
dsum2  = 0

dsum1 = b**(1+alpha0)/(1+alpha0)
do i=1,nquad
x     = xs(i)
wht   = whts(i)
dsum2 = dsum2 + x**(alpha0)*wht 
end do
derr = abs(dsum1-dsum2)/abs(dsum1)
dmax = max(derr,dmax)
end do

call prin2("max int error = ",dmax)


end program
