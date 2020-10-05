module test_makequad_funs

use iso_c_binding
use adapquad
use makequad

! a structure containing data which is passed to the subroutines that
! supply the values of the input functions to the quadrature code
type     fundata
integer                       :: nalphas
double precision              :: a 
double precision              :: b 
double precision, allocatable :: alphas(:)
end type fundata

contains


subroutine test_quad1()
implicit double precision (a-h,o-z)
!
!  Test the ggquad routine by building a quadrature for integrals of the form
!
!         1
!     \int  |x|^alpha    P_j(x)
!        -1
!
!  for a range of alpha's.
!

type(fundata), pointer         :: data
type(c_ptr)                    :: userptr
double precision, allocatable  :: xs(:), whts(:)
double precision, target       :: alpha
double precision, allocatable  :: pols(:), vals(:), vals0(:)

eps0    = epsilon(0.0d0)

ifsing  = 1
a       = 0.00d0
b       = 1.00d0
alpha   = -0.5d0
norder  = 19
nfuns   = norder+1

userptr = c_loc(alpha)


call elapsed(t1)
call ggquad(ifsing,a,b,nfuns,funpolys2,userptr,nquad,xs,whts)
call elapsed(t2)

call prin2("ggquad time = ",t2-t1)
! call prini("after ggquad, nquad = ",nquad)
! call prin2("after ggquad, xs = ",xs)
! call prin2("after ggquad, whts = ",whts)


!
!  Check the accuracy of the formula
!


allocate(vals(0:norder),vals0(0:norder),pols(0:norder))
dmax = 0

do i=0,norder
if (mod(i,2) == 1) then
vals0(i) = 0
m = (i-1)/2

vals0(i) = (-1)**m *  gamma(m+0.5d0-0.5d0*alpha) * gamma(1.0d0+0.5d0*alpha)   / &
                  (2* gamma(0.5d0-0.5d0*alpha)*gamma(m+2d0+alpha/2)) 

else
m = i/2
vals0(i) = (-1)**m *  gamma(m-0.5d0*alpha) * gamma(0.5d0+0.5d0*alpha)   / &
                  (2* gamma(-0.5d0*alpha)*gamma(m+1.5d0+alpha/2)) 

endif
vals0(i) = vals0(i)*sqrt(i+0.5d0)
end do


vals=0
do l=1,nquad
x     = xs(l)
wht   = whts(l)
call leges(norder+1,x,pols)
vals = vals + pols*abs(x)**alpha*wht 
end do

derr = maxval(abs(vals-vals0))
call prin2("maximum absolute error in integrals = ",derr)

dtime = t2-t1

deallocate(vals,vals0,pols)

close(iw)


end subroutine


subroutine test_quad2()
implicit double precision (a-h,o-z)
!
!  Test the ggquad_prod routine by building a quadrature for integrals of the form
!
!         1
!     \int      (|x|^\alpha P(x) ) dx
!        -1
!
!  for alpha in a specified interval and P(X) a polynomial of degree less than
!  or equal to m. 
!

type(fundata), pointer         :: data
type(c_ptr)                    :: userptr
double precision, allocatable  :: xs(:), whts(:)
double precision, allocatable  :: alphas(:)
double precision, allocatable  :: pols(:), vals(:), vals0(:)

eps0    = epsilon(0.0d0)
epsadap = eps0*2

ifsing  = 1
a       =-1.0d0
b       = 1.0d0
norder  = 29

!
!  Discretize the interval [alpha1, alpha2] ... oversample it thoroughly
!


call elapsed(time1)
nalphas =  100
npolys  =  norder+1
alpha1  = -0.50d0 
alpha2  =  0.50d0

allocate(data)
allocate(data%alphas(nalphas))

data%nalphas = nalphas
data%a       = a
data%b       = b

!call legendre_quad(nalphas-1,data%alphas,whts)
call legendre_quad(nalphas-1,xs,whts)

data%alphas(1:nalphas-1) = xs*(alpha2-alpha1)/2 + (alpha1+alpha2)/2
data%alphas(nalphas)     = 0.0d0
userptr = c_loc(data)

call elapsed(t1)
call ggquad_prod(ifsing,a,b,nalphas,npolys,funalphas,funpolys,userptr,userptr,nquad,xs,whts)
call elapsed(t2)

call prin2("ggquad_prod time = ",t2-t1)
call prini("after ggquad_prod, nquad = ",nquad)
call prin2("after ggquad_prod, xs = ",xs)
call prin2("after ggquad_prod, whts = ",whts)

!
!  Check the accuracy of the formula for nn random values of alpha
!  in the specified range by comparison with a known formula for
!  the integral 
!
!               1
!           \int P_j(x) x^alpha  dx
!               0  
!


nn = 100000
allocate(alphas(nn))
do i=1,nn
call random_number(dd)
!dd = (i-1.0d0)/(nn-1.0d0)
alphas(i) = alpha1 + (alpha2-alpha1)*dd
end do


allocate(vals(0:norder),vals0(0:norder),pols(0:norder))
dmax = 0

do j=1,nn
alpha = alphas(j)

do i=0,npolys-1

if (mod(i,2) == 1) then
vals0(i) = 0
m = (i-1)/2
! vals0(i) = (-1)**m *  gamma(m+0.5d0-0.5d0*alpha) * gamma(1.0d0+0.5d0*alpha)   / &
!                   (2* gamma(0.5d0-0.5d0*alpha)*gamma(m+2d0+alpha/2)) 

vals0(i) = 0

else
m = i/2
vals0(i) = (-1)**m *  gamma(m-0.5d0*alpha) * gamma(0.5d0+0.5d0*alpha)   / &
                  (2* gamma(-0.5d0*alpha)*gamma(m+1.5d0+alpha/2)) 

vals0(i) = vals0(i)*2

endif

vals0(i) = vals0(i)*sqrt(i+0.5d0)

end do


vals=0
do l=1,nquad
x     = xs(l)
wht   = whts(l)
call leges(norder+1,x,pols)
vals = vals + pols*abs(x)**alpha*wht 
end do

derr = maxval(abs(vals-vals0))
dmax = max(derr,dmax)
end do

call prin2("maximum absolute error in integrals = ",dmax)
call prina("")
call elapsed(time2)

dd1 = nquad
dd2 = ceiling(krank/2.0d0)
dd  = dd1/dd2




end subroutine



subroutine test_quad3()
implicit double precision (a-h,o-z)

!
!  Test the code for constructing "interpolation quadratures"
!

type(fundata), pointer         :: data
type(c_ptr)                    :: userptr
double precision, allocatable  :: xs(:), whts(:), vals(:)
double precision, allocatable  :: alphas(:)
integer, target                :: j

eps0    = epsilon(0.0d0)

ifsing  = 1
a       =-1.0d0
b       = 1.0d0

!
!  Discretize the interval [alpha1, alpha2] ... oversample it thoroughly
!

nfuns   = 30
userptr = c_loc(data)

call elapsed(t1)
call ggquad_interp(a,b,nfuns,funcosines,userptr,nquad,xs,whts,rcond)
call elapsed(t2)

allocate(vals(nfuns))

userptr = c_loc(j)

dmax = 0 
do j=1,nfuns-1
!dsum1 = sum(sqrt(abs(xs))*cos((j-1)*xs)*whts)
dsum1 = sum(cos((j-1)*xs)*whts)

call  adapint(ier,eps0,a,b,funtest,userptr,dsum2,ntotal)
derr = abs(dsum1-dsum2)
dmax = max(derr,dmax)
end do
call prin2("max int error = ",dmax)

end subroutine



subroutine test_quad4()
implicit double precision (a-h,o-z)

!
!  Test the code for constructing "interpolation quadratures"
!

type(fundata), pointer         :: data
type(c_ptr)                    :: userptr
double precision, allocatable  :: xs(:), whts(:)
double precision, allocatable  :: alphas(:), vals(:)
double precision, target       :: alpha

eps0    = epsilon(0.0d0)

ifsing  = 1
a       =-1.0d0
b       = 1.0d0

nalphas = 100
alpha1  =-0.50d0
alpha2  = 0.50d0

allocate(data)
allocate(data%alphas(nalphas))
data%nalphas = nalphas
call legendre_quad(nalphas,data%alphas,whts)
data%alphas = data%alphas*(alpha2-alpha1)/2 + (alpha1+alpha2)/2

!
!  Discretize the interval [alpha1, alpha2] ... oversample it thoroughly
!

nfuns   = nalphas+1
userptr = c_loc(data)

call elapsed(t1)
call ggquad_interp_sing(nfuns,funalphas2,userptr,nquad,xs,whts,rcond)
call elapsed(t2)


nn = 100
allocate(alphas(nn))
do i=1,nn
call random_number(dd)
!dd = (i-1.0d0)/(nn-1.0d0)
alphas(i) = alpha1 + (alpha2-alpha1)*dd
end do


allocate(vals(nquad))
dmax = 0

do j=1,nn
!alpha   = data%alphas(j)
alpha   = alphas(j)
userptr =  c_loc(alpha)

call funtest2(nquad,xs,vals,userptr)
vals  = vals*whts
dsum1 = sum(vals)
call  adapint(ier,eps0,a,b,funtest2,userptr,dsum2,ntotal)

derr  = abs(dsum1-dsum2)
dmax = max(derr,dmax)

end do
call prin2("integration error max = ",dmax)


end subroutine



subroutine test_quad5()
implicit double precision (a-h,o-z)
!
!  Test teh ggquad_sing code.
! 

type(fundata), pointer         :: data
type(c_ptr)                    :: userptr
double precision, allocatable  :: xs(:), whts(:)
double precision, allocatable  :: alphas(:)
double precision, allocatable  :: pols(:), vals(:), vals0(:)

eps0    = epsilon(0.0d0)
a       =-1.0d0
b       = 1.0d0

!
!  Discretize the interval [alpha1, alpha2] ... oversample it thoroughly
!

iw = 20
open(iw)


alpha1=0

do ii=1,9
nalphas =  100
alpha1  = alpha1-0.1d0
alpha2  = 0.0d0

allocate(data)
allocate(data%alphas(nalphas))
data%nalphas = nalphas
data%a       = a
data%b       = b

call legendre_quad(nalphas,data%alphas,whts)
data%alphas = data%alphas*(alpha2-alpha1)/2 + (alpha1+alpha2)/2
userptr = c_loc(data)

call elapsed(t1)
call ggquad_sing(nalphas,funalphas,userptr,nquad,xs,whts)
call elapsed(t2)

! call prin2("ggquad_sing time = ",t2-t1)
! call prini("after ggquad_sing, nquad = ",nquad)
! call prin2("after ggquad_sing, xs = ",xs)
! call prin2("after ggquad_sing, whts = ",whts)

write (iw,"(A,I1,A)") "XX",ii,"={"

do i=1,nquad
if (i .eq. nquad) then
write (iw,"(A,F36.30,A,F36.30,A)") "{ ",log(abs(xs(i)))," , ", alpha1,"}};"
else
write (iw,"(A,F36.30,A,F36.30,A)") "{ ",log(abs(xs(i)))," , ", alpha1,"},"
endif
end do

deallocate(data)
end do


nalphas =  100
alpha1  =  -0.95d0
alpha2  =  0.00d0

allocate(data)
allocate(data%alphas(nalphas))
data%nalphas = nalphas
data%a       = a
data%b       = b

call legendre_quad(nalphas,data%alphas,whts)
data%alphas = data%alphas*(alpha2-alpha1)/2 + (alpha1+alpha2)/2
userptr = c_loc(data)

call elapsed(t1)
call ggquad_sing(nalphas,funalphas,userptr,nquad,xs,whts)
call elapsed(t2)


write (iw,"(A)") "XX10={"

do i=1,nquad
if (i .eq. nquad) then
write (iw,"(A,F36.30,A,F36.30,A)") "{ ",log(abs(xs(i)))," , ", alpha1,"}};"
else
write (iw,"(A,F36.30,A,F36.30,A)") "{ ",log(abs(xs(i)))," , ", alpha1,"},"
endif
end do

deallocate(data)


nalphas =  100
alpha1  =  -0.90d0
alpha2  =  -0.00d0

allocate(data)
allocate(data%alphas(nalphas))
data%nalphas = nalphas
data%a       = a
data%b       = b

call legendre_quad(nalphas,data%alphas,whts)
data%alphas = data%alphas*(alpha2-alpha1)/2 + (alpha1+alpha2)/2
userptr = c_loc(data)

call elapsed(t1)
call ggquad_sing(nalphas,funalphas,userptr,nquad,xs,whts)
call elapsed(t2)


write (iw,"(A)") "XX11={"

do i=1,nquad
if (i .eq. nquad) then
write (iw,"(A,F36.30,A,F36.30,A)") "{ ",log(abs(xs(i)))," , ", alpha1,"}};"
else
write (iw,"(A,F36.30,A,F36.30,A)") "{ ",log(abs(xs(i)))," , ", alpha1,"},"
endif
end do

deallocate(data)

close(iw)


stop


!
!  Check the accuracy of the formula for nn random values of alpha
!

nn = 1000
allocate(alphas(nn))
do i=1,nn
!call random_number(dd)
dd = (i-1.0d0)/(nn-1.0d0)
alphas(i) = alpha1 + (alpha2-alpha1)*dd
end do

dmax = 0
do j=1,nalphas
alpha = alphas(j)

sum1 = 0
do i=1,nquad
x = xs(i)
wht = whts(i)
sum1 = sum1 + abs(x)**alpha * wht
end do

sum2 = 2.0d0/(1.0d0+alpha)
derr = abs(sum1-sum2)
dmax = max(derr,dmax)
end do

call prin2("maximum error in integral = ",dmax)

end subroutine


subroutine test_quad6()
implicit double precision (a-h,o-z)
!
!  Test the ggquad_prod routine by building a quadrature for integrals of the form
!
!         b
!     \int      |x|^\alpha dx
!         a
!
!  with alpha close to -1
!

type(fundata), pointer         :: data
type(c_ptr)                    :: userptr
double precision, allocatable  :: xs(:), whts(:)
double precision, allocatable  :: alphas(:)
double precision, allocatable  :: pols(:), vals(:), vals0(:)

eps0    = epsilon(0.0d0)
a       =-1.0d0
b       = 1.0d0

!
!  Discretize the interval [alpha1, alpha2] ... oversample it thoroughly
!
nalphas =  100
norder  =  19
npolys  =  norder+1
alpha1  = -0.50d0 
alpha2  =  0.50d0

allocate(data)
allocate(data%alphas(nalphas))
data%nalphas = nalphas
data%a       = a
data%b       = b

call legendre_quad(nalphas,data%alphas,whts)
data%alphas = data%alphas*(alpha2-alpha1)/2 + (alpha1+alpha2)/2
userptr = c_loc(data)

call elapsed(t1)
call ggquad_singprod(nalphas,npolys,funalphas,funpolys,userptr,userptr,nquad,xs,whts)
call elapsed(t2)

call prin2("ggquad_singprod time = ",t2-t1)
! call prini("after ggquad_singprod, nquad = ",nquad)
! call prin2("after ggquad_singprod, xs = ",xs)
! call prin2("after ggquad_singprod, whts = ",whts)

!
!  Check the accuracy of the formula for nn random values of alpha
!  in the specified range by comparison with a known formula for
!  the integral 
!
!               1
!           \int P_j(x) |x|^alpha  dx
!              -1
!


nn = 1000
allocate(alphas(nn))
do i=1,nn
call random_number(dd)
!dd = (i-1.0d0)/(nn-1.0d0)
alphas(i) = alpha1 + (alpha2-alpha1)*dd
end do

allocate(vals(0:norder),vals0(0:norder),pols(0:norder))
dmax = 0

do j=1,nn
alpha = alphas(j)
do i=0,npolys-1

if (mod(i,2) == 1) then
vals0(i) = 0
! m = (i-1)/2
! vals0(i) = 2*(-1)**m *  gamma(m+0.5d0-0.5d0*alpha) * gamma(1.0d0+0.5d0*alpha)   / &
!                    (2* gamma(0.5d0-0.5d0*alpha)*gamma(m+2d0+alpha/2)) 

else
m = i/2
vals0(i) = 2*(-1)**m *  gamma(m-0.5d0*alpha) * gamma(0.5d0+0.5d0*alpha)   / &
                  (2* gamma(-0.5d0*alpha)*gamma(m+1.5d0+alpha/2)) 

endif

vals0(i) = vals0(i)*sqrt(i+0.5d0)

end do


vals=0
do l=1,nquad
x     = xs(l)
wht   = whts(l)
call leges(norder+1,x,pols)
vals = vals + pols*abs(x)**alpha*wht 
end do

derr = maxval(abs(vals-vals0))
dmax = max(derr,dmax)
end do

call prin2("maximum absolute error in integrals = ",dmax)
call prina("")


end subroutine



subroutine funalphas(nfuns,n,xs,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: xs(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
type(fundata), pointer        :: data

call c_f_pointer(userptr,data)
pi = acos(-1.0d0)
do j=1,nfuns
alpha     = data%alphas(j)
vals(:,j) = abs(xs)**alpha
end do

end subroutine


subroutine funpolys(nfuns,n,xs,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: xs(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
type(fundata), pointer        :: data
!
!  Evaluate the Legendre polynomials
!
call c_f_pointer(userptr,data)
a      = data%a
b      = data%b

do i=1,n
x  = xs(i)
xx = (2*x - (b+a) ) /(b-a)
call leges(nfuns,xx,vals(i,:))


! unnormalized
! do j=1,nfuns
! vals(i,j) = vals(i,j)/sqrt(j-0.5d0)
! end do
vals(i,:) = vals(i,:) * sqrt(2/(b-a))

end do

end subroutine




subroutine funalphas2(nfuns,n,xs,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: xs(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
type(fundata), pointer        :: data
type(c_ptr)                   :: userptr2
double precision, target      :: alpha

call c_f_pointer(userptr,data)
pi = acos(-1.0d0)

do j=1,nfuns-1
alpha     = data%alphas(j)
userptr2  = c_loc(alpha)
call funtest2(n,xs,vals(:,j),userptr2)
end do

vals(:,nfuns) = 1.0d0

end subroutine

subroutine funtest2(n,xs,vals,userptr)
implicit double precision (a-h,o-z)
double precision :: xs(n)
double precision :: vals(n)
type(c_ptr)      :: userptr
double precision, pointer :: alpha
call c_f_pointer(userptr,alpha)
vals = abs(xs)**alpha * (cos(xs)-sin(xs))
end subroutine




subroutine funcosines(nfuns,n,xs,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: xs(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
double precision, pointer     :: alpha

pi = acos(-1.0d0)

do j=1,nfuns-1
vals(:,j) = cos( (j-1)*xs )
! * sqrt(abs(xs))
end do

! the collection better include 1 !!!!
vals(:,nfuns) = 1 

end subroutine


subroutine funtest(n,xs,vals,userptr)
implicit double precision (a-h,o-z)
double precision :: xs(n)
double precision :: vals(n)
type(c_ptr)      :: userptr
integer, pointer :: j

call c_f_pointer(userptr,j)
vals = cos((j-1)*xs)*sqrt(abs(xs))

end subroutine




subroutine funpolys2(nfuns,n,xs,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: xs(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
double precision, pointer     :: alpha

call c_f_pointer(userptr,alpha)


do i=1,n
x  = xs(i)
xx = 2*x-1
call leges(nfuns,xx,vals(i,:))
vals(i,:) = vals(i,:) * abs(x)**alpha

! call leges(nfuns,x,vals(i,:))
! vals(i,:) = vals(i,:) * abs(x)**alpha

! call leges(nfuns,x,vals(i,:))
! vals(i,:) = vals(i,:) * log(abs(x))

end do

end subroutine



end module

program test_makequad
use makequad
use test_makequad_funs
implicit double precision (a-h,o-z)

!call test_quad1()
!call prina("")

!call test_quad2()
!call prina("")

!call test_quad3()
!call prina("")

call test_quad4()
call prina("")

!call test_quad5()
!call prina("")

! call test_quad6()

end program
