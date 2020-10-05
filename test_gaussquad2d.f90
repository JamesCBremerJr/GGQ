module test_gaussquad2d_funs

use utils
use legendre
use bilege
use bilegepw
use gaussian_quad2d
use iso_c_binding

type           funuser_data
integer                          :: n1, n2

integer                          :: nrs, nthetas
double precision, allocatable    :: rs(:)
double precision, allocatable    :: thetas(:)

end type       funuser_data

contains



subroutine funuser(nfuns,n,xs,ys,whts,vals,userptr)
implicit double precision (a-h,o-z)
type (c_ptr)                :: userptr
double precision            :: xs(:),ys(:),whts(:)
double precision            :: vals(:,:)
type(funuser_data), pointer :: fundata

call c_f_pointer(userptr,fundata)

n1      = fundata%n1
n2      = fundata%n2
nthetas = fundata%nthetas
nrs     = fundata%nrs

do i=1,n
u   = xs(i)
v   = ys(i)
wht = whts(i)

idx   = 0
do i1=1,nthetas
do i2=1,nrs
r0     = fundata%rs(i1)
theta0 = fundata%thetas(i2)

theta  = theta0*u
dd     = r0*sin(theta0) / (sin(theta) - r0*sin(theta-theta0) )
r      = v*dd

do j1=1,n1
do j2=1,n2
idx         = idx+1
vals(i,idx) = dd**(j1+2)*v**(j1+1) * u**(j2) * sqrt(wht)
end do
end do
end do
end do
end do

end subroutine



subroutine funregion(x,y,wht,isvalid,userptr)
implicit double precision (a-h,o-z)
type (c_ptr)        :: userptr
isvalid = 1
if (x .gt. 1 .OR. x .lt. 0) isvalid = 0
if (y .gt. 1 .OR. y .lt. 0) isvalid = 0
end subroutine


end module


program test_gaussquad2d

use utils
use plot
use legendre
use bilege
use bilegepw
use gaussian_quad2d
use test_gaussquad2d_funs
use iso_c_binding

implicit double precision (a-h,o-z)

type(bilegepw_disc)           :: disc
type(funuser_data), pointer   :: fundata
type(c_ptr)                   :: userptr
double precision, allocatable :: xs0(:),ys0(:),whts0(:), rects(:,:)
double precision, allocatable :: vals(:,:), valsq(:,:), rnorms(:), r(:,:)
integer, allocatable          :: ipivs(:)
double precision, allocatable :: xs(:),ys(:),whts(:)

! double precision, allocatable :: vals(:,:), valsq(:,:), r(:,:), rnorms(:)
!,vals0(:,:),rints(:),rints0(:)
! double precision, allocatable :: coefs(:,:)
! integer, allocatable          :: ipivs(:)

pi         = acos(-1.0d0)

a          =  0.0d0
b          =  1.0d0
c          =  0.0d0
d          =  1.0d0

epsdisc    = 1.0d-13
epsqr      = 1.0d-12
epsnewt    = 1.0d-8

n1         = 8
n2         = 8
nrs        = 20
nthetas    = 20

r0         = 0.9d0
r1         = 1.0d0

theta0     = pi/2+0.1d0
theta1     = pi/2-0.1d0

nfuns      = n1*n2*nrs*nthetas

allocate(fundata)
userptr        = c_loc(fundata)
fundata%n1      = n1
fundata%n2      = n2
fundata%nrs     = nrs
fundata%nthetas = nthetas

allocate(fundata%rs(nrs))
allocate(fundata%thetas(nthetas))

do i=1,nrs
fundata%rs(i) = r0 + (r1-r0)*(i-1.0d0)/(nrs-1.0d0)
end do

do i=1,nthetas
fundata%thetas(i) = theta0 + (theta1-theta0)*(i-1.0d0)/(nrs-1.0d0)
end do


call bilegepw_init(20,disc,a,b,c,d)
call bilegepw_adap(epsdisc,disc,nfuns,funuser,userptr)

call bilegepw_rects(disc,nrects,rects)
call plot_rectangles("*", "disc.pdf", nrects, rects )
call bilegepw_quad(disc,nquad0,xs0,ys0,whts0)

call prini("nleaves = ",disc%nleaves)
call prini("nquad0 = ",nquad0)



allocate(vals(nquad0,nfuns) )
call funuser(nfuns,nquad0,xs0,ys0,whts0,vals,userptr)
call qrdecomp(epsqr,nquad0,nfuns,vals,krank,ipivs,valsq,r)

allocate(rnorms(krank))
do i=1,krank
rnorms(i) = r(i,i)
end do

call prini("krank  = ",krank)
call prin2("rnorms = ",rnorms(1:krank))


call chebquad2d_pw(disc,krank,valsq,nquad,xs,ys,whts)
call prini("after chebquad, nquad   = ",nquad)
call prin2("after chebquad, xs   = ",xs)
call prin2("after chebquad, ys   = ",ys)
call prin2("after chebquad, whts = ",whts)

call gaussquad2_pw(epsnewt,disc,krank,valsq(:,1:krank),nquad,xs,ys,whts)
call prini("after gaussquad, nfuns = ",nfuns)
call prini("after gaussquad, nquad = ",nquad)
call prin2("after gaussquad, xs   = ",xs)
call prin2("after gaussquad, ys   = ",ys)
call prin2("after gaussquad, whts = ",whts)

! dmax = 0

! do nn=0,norder
! do i=0,nn
! j = nn-i

! dsum1 = sum(xs**i * ys**j * whts)
! dsum2 = ( 1+(-1)**i )  * ( 1+(-1)**j ) / (1.0d0+i) * 1/(1.0d0+j)
! derr  = abs(dsum1-dsum2)
! dmax  = max(derr,dmax)
! end do
! end do

! call prin2("max integration error = ",dmax)


! norder = 13
! !call guassian_rect(norder,nquad,xs,ys,whts)

! call polys_rect(norder,nquad,xs,ys,whts)

end program
