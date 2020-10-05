module test_legepw_funs

use utils
use legendre
use iso_c_binding

integer                       :: nalphas
double precision, allocatable :: alphas(:)

contains

subroutine funalphas(nfuns,n,xs,whts,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: xs(n), whts(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
do j=1,nfuns
vals(:,j) = abs(xs)**alphas(j) * sqrt(whts)
end do

end subroutine


subroutine funtest(nquad,xs,whts,vals,ders)
implicit double precision (a-h,o-z)
double precision              :: xs(nquad), whts(nquad), vals(nquad), ders(nquad)

alpha1 = 0.25d0
alpha2 =-0.30d0

vals =         abs(xs)**(alpha1) + &
       3.0d0 * abs(xs)**(alpha2)
ders =        (alpha1) * abs(xs)**(alpha1-1)*sign(1.0d0,xs)  + &
       3.0d0 *(alpha2) * abs(xs)**(alpha2-1)*sign(1.0d0,xs)

vals = vals*sqrt(whts)
ders = ders*sqrt(whts)



end subroutine

end module

program test_legepw

use legendre
use legepw
use test_legepw_funs

implicit double precision (a-h,o-z)
type(c_ptr)                   :: userptr
type(legepw_disc)             :: disc
double precision, allocatable :: ab0(:,:)
double precision, allocatable :: xs(:), whts(:)
double precision, allocatable :: xs0(:), whts0(:)
double precision, allocatable :: vals(:), coefs(:), ders(:)
double precision, allocatable :: valsout(:), dersout(:)
double precision, allocatable :: valsout0(:), dersout0(:)


epsadap = 1.0d-13
epsqr   = 1.0d-13
nlege   = 30

nints0  = 2
allocate(ab0(2,nints0))
ab0(1,1) = -1.0d0
ab0(2,1) = -1.0d-15
ab0(1,2) =  1.0d-15
ab0(2,2) =  1.0d0

a       =-1.0d0
b       = 1.0d0

nalphas = 100
alpha1  =-0.50d0
alpha2  = 0.50d0

allocate(alphas(nalphas))
do i=1,nalphas
alphas(i) = alpha1 + (alpha2-alpha1) * (i-1.0d0)/(nalphas-1.0d0)
end do
call prin2("alphas = ",alphas)

call legepw_init(disc,nlege,nints0,ab0)



!
!  Discretize the functions
!

call elapsed(t1)
call legepw_adap(epsadap,nalphas,funalphas,disc,userptr)
call elapsed(t2)
call legepw_quad(disc,nquad,xs,whts)


call prin2("legepw time = ",t2-t1)
call prini("after legepw, nquad = ",nquad)
call prin2("after legepw, ab = ",disc%ab)

!
!  Construct a target quadrature for testing
!
nn = 1000
call legendre_quad(nn,xs0,whts0)
xs0   = (b-a)/2*xs0 + (b+a)/2
whts0 = (b-a)/2*whts0


allocate(valsout(nn), valsout0(nn) )
allocate(dersout(nn), dersout0(nn) )

!
!  Test the evaluation routines
!


allocate(vals(nquad), ders(nquad), coefs(nquad) )

vals = sqrt(whts)

call funtest(nquad,xs,whts,vals,ders)
call legepw_coefs(disc,vals,coefs)
call funtest(nn,xs0,whts0,valsout0,dersout0)


call elapsed(t1)
do i = 1,nn
x   = xs0(i)
wht = whts0(i)
call legepw_evalder(disc,coefs,x,val,der)
valsout(i) = val*sqrt(wht)
dersout(i) = der*sqrt(wht)
end do
call elapsed(t2)

derr1 = norm2(valsout-valsout0) / norm2(abs(valsout0))
derr2 = norm2(dersout-dersout0) / norm2(abs(dersout0))
tavg1 = (t2-t1)/(nn+0.0d0)

call prin2("L^2 eval error for f(x)  = ",derr1)
call prin2("L^2 eval error for f'(x) = ",derr2)
call prin2("avg eval time            = ",tavg1)
call prina("")

!
!  Test the interpolation routines
!

call legepw_diff(disc,vals,ders)

call elapsed(t1)
do i = 1,nn
x   = xs0(i)
wht = whts0(i)
call legepw_interp(disc,vals,x,val)
call legepw_interp(disc,ders,x,der)

valsout(i) = val*sqrt(wht)
dersout(i) = der*sqrt(wht)

end do
call elapsed(t2)
tavg2 = (t2-t1)/(nn+0.0d0)

derr3 = norm2(valsout-valsout0) / norm2(abs(valsout0))
derr4 = norm2(dersout-dersout0) / norm2(abs(dersout0))

call prin2("L^2 interp error for f(x)  = ",derr3)
call prin2("L^2 interp error for f'(x) = ",derr4)
call prin2("avg interp time            = ",tavg2)



end program

