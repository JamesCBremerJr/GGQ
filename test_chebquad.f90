module test_chebquad_funs
use iso_c_binding

type     funalphas_data
integer                       :: nalphas
double precision, allocatable :: alphas(:)
end type funalphas_data

contains

subroutine funalphas(nfuns,n,xs,whts,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: xs(n), whts(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
type(funalphas_data), pointer :: fundata
call c_f_pointer(userptr,fundata)
do j=1,nfuns
vals(:,j) = xs**fundata%alphas(j) * sqrt(whts)
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


program test_chebquad

use utils
use legendre
use linalg
use legepw
use chebyshev_quad
use iso_c_binding
use test_chebquad_funs


implicit double precision (a-h,o-z)

type(legepw_disc)             :: disc
type(c_ptr)                   :: userptr
type(funalphas_data), pointer :: fundata
double precision, allocatable :: xs(:), whts(:)
double precision, allocatable :: xscheb(:), whtscheb(:)
integer, allocatable          :: ipivs(:)
double precision, allocatable :: vals(:,:), valsq(:,:), r(:,:), rnorms(:)
double precision, allocatable :: alphas(:)
double precision, pointer     :: alpha0

eps0    = epsilon(0.0d0)
epsadap = eps0*2
epsdisc = 1.0d-13
epsqr   = 1.0d-13
nlege   = 24

a0      = 0.0d0
a       = eps0**2
b       = 1.0d0

!
!  Adaptively construct a discretization scheme for the input functions |x|^alpha
!

nalphas = 100
alpha1  =-0.5d0
alpha2  = 0.5d0


call disc_alphas(nalphas,alpha1,alpha2,fundata)

call prini("nalphas = ",nalphas)
!call prin2("alphas  = ",fundata%alphas)

call legepw_init(disc,nlege,a,b)
userptr = c_loc(fundata)

call elapsed(t1)
call legepw_adap(epsdisc,nalphas,funalphas,disc,userptr)
call elapsed(t2)

call legepw_quad(disc,nquad,xs,whts)

call prin2("legepw_adap time = ",t2-t1)
call prini("nlege = ",nlege)
call prini("nints = ",disc%nints)
call prini("nquad = ",nquad)

!
!  Evaluate the input functions at the discretization nodes
!
allocate( vals(nquad,nalphas) )

call elapsed(t1)
call funalphas(nalphas,nquad,xs,whts,vals,userptr)
call elapsed(t2)


call prin2("funalphas time = ",t2-t1)

!
!  Orthonormalize the functions ...
!

call elapsed(t1)
call qrdecomp(epsqr,nquad,nalphas,vals,krank,ipivs,valsq,r)
allocate(rnorms(krank))
do i=1,krank
rnorms(i) = abs(r(i,i))
end do
call elapsed(t2)
call prin2("qrdecomp time = ",t2-t1)
call prini("after qrdecomp, krank = ",krank)
call prin2("after qrdecomp, rnorms = ",rnorms)

!
! ... or perform interpolative decomposition
!

! call elapsed(t1)
! call intdecomp(epsqr,nquad,nalphas,vals,krank,ipivs,r)
! allocate(valsq(nquad,krank), dersq(nquad,krank) )
! valsq = vals(:,ipivs)
! dersq = ders(:,ipivs)
! call elapsed(t2)
! call prin2("intdecomp time = ",t2-t1)
! call prini("after intdecomp, rank = ",krank)

!
!  Construct the Chebyshev quadrature  
!

call elapsed(t1)
call chebquad(disc,krank,valsq,nquadcheb,xscheb,whtscheb)
call elapsed(t2)
call prin2("chebquad time  = ",t2-t1)
call prin2("chebquad xs    = ",xscheb)
call prin2("chebquad whts  = ",whtscheb)

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


do i=1,nquadcheb
x     = xscheb(i)
wht   = whtscheb(i)
dsum2 = dsum2 + abs(x)**alpha0*wht 
end do

derr = abs(dsum1-dsum2)/abs(dsum1)
dmax = max(derr,dmax)
!print *,alpha0,dsum1,dsum2,dsum1-dsum2

end do

call prin2("max int error = ",dmax)

end program
