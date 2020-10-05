!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains routines for constructing generalized Gaussian quadrature rules.
!  They are wrappers around the codes in chebyshev_quad.f90 and gaussian_quad.f90.
!
!  The following subroutines are publicly callable:
!
!    ggquad - construct a generalized Gaussian quadrature rule for a collection of
!     input functions given on an interval (a,b) and specified via an external
!     subroutine
!
!    ggquad_prod - construct a generalized Gaussian quadrature ruke for a collection of
!     functions of the form
!
!        f_i(x) g_j(x)     i=1,...,n1,  j=1,...,m1,
!
!     where the f_i and g_j are specified by the user via external subroutines
!
!    ggquad_interp - use Kirill Serkh's procedure to construct an n-point quadrature
!      which integrates a collection of n smooth (or, very mildly singular)  functions 
!      and whose nodes are stable interpolation points for the collection
!
!      IMPORTANT NOTE: THE FUNCTION f(x) = 1 MUST BE AMONG THE INPUT FUNCTIONS OR
!      THE RESULTING QUADRATURE MIGHT NOT ACCURATELY INTEGRATE THE INPUT FUNCTIONS
!
!    ggquad_interp_sing - use  Kirill Serkh's procedure to construct an n-point 
!      "interpolation quadrature" for a collection of n functions given on the interval
!      (-1,1) and which are singular at 0 
!
!      IMPORTANT NOTE: THE FUNCTION f(x) = 1 MUST BE AMONG THE INPUT FUNCTIONS OR
!      THE RESULTING QUADRATURE MIGHT NOT ACCURATELY INTEGRATE THE INPUT FUNCTIONS

!    ggquad_sing - construct a generalized Gaussian quadrature rule for a collection
!     of functions given on the interval (-1,1) that are strongly singular at 0 
!     (e.g., |x|^alpha for alpha close to -1)
!
!    ggquad_singprod - construct a generalized Gaussian quadrature rule for a collection
!     of functions of the form
!
!        f_i(x) g_j(x)     i=1,...,n1, j=1,...,m1
!
!     where the f_i are given on the interval (-1,1) and are strongly singular at 0,
!     and the g_j(x) are nonsingular.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module makequad
use utils
use legendre
use legepw
use linalg
use chebyshev_quad
use gaussian_quad
use iso_c_binding



! set this to one to suppress output from these routines
integer, private, parameter                 :: ifsuppress = 0

interface

subroutine ggfun(nfuns,n,xs,vals,userptr)
!
!  Return the (n,nfuns) matrix vals whose jth column gives the
!  values of the jth input function at a collection of points
!  specified by the xs array.  The parameter userptr is a "void *"
!  pointer which is supplied by the user to the quadrature routines
!  and passed on to this routine.
!
import c_ptr
double precision           :: xs(n)
double precision           :: vals(:,:)
type(c_ptr)                :: userptr
end subroutine

end interface

type     ggwrapper_data
type(c_ptr)                       :: userptr
procedure(ggfun), pointer, nopass :: userfun
integer                           :: ifscale
integer                           :: ifodd
integer                           :: ifeven
double precision                  :: dsub
end type ggwrapper_data

contains


subroutine ggquad(ifsing,a,b,nfuns,fun,userptr,nquad,xs,whts)
implicit double precision (a-h,o-z)
procedure(ggfun)                           :: fun
double precision, allocatable, intent(out) :: xs(:), whts(:)
type(c_ptr)                                :: userptr
!
!   Construct a generalized Gaussian quadrature for a collection of functions
!   of the form
!
!        f_i(x),      i=1,...,n,
!
!   where the f_i are supplied by the user via an external subroutine.
!
!  Input parameters:
!    ifsing - an integer parameter which, when set to 1, indicates that the
!      input functions are given on an interval of the form (a,b) with
!      a <= 0 <=  b, and that they have mild singularities at 0
!
!    (a,b) - the interval on which the input functions are given
!    nfuns - the number of input functions
!    fun - an external subroutine conforming to the ggfun interface which
!      supplies the values of the input functions
!    userptr - a user-supplied "void *" pointer which is passed to the subroutine
!      fun
!
!  Output parameters:
!    (nquad,xs,whts) - the resulting quadrature rule
!

type(legepw_disc)                   :: disc
type(c_ptr)                         :: wrapptr
type(ggwrapper_data), pointer       :: wrapdata
double precision, allocatable       :: ab0(:,:)
double precision, allocatable       :: vals0(:,:),  vals(:,:), r(:,:)
double precision, allocatable       :: xs0(:), whts0(:)
double precision, allocatable       :: us(:), uwhts(:), rnorms(:)
integer, allocatable                :: ipivs(:)

eps0   = epsilon(0.0d0)

epsdisc  = 1.0d-13
epsqr    = 1.0d-13
epsnewt  = 1.0d-7

if (eps0 .lt. 1.0d-17) then
epsdisc  = 1.0d-15
epsqr    = 1.0d-15
epsnewt  = 1.0d-9
endif

if (eps0 .lt. 1.0d-30) then
epsdisc  = 1.0d-30
epsqr    = 1.0d-30
epsnewt  = 1.0d-15
endif

nlege = 30
dsub  = 1
a0    = a
b0    = b

if (ifsing .eq. 1) then
dsub   = 3
a0     =-abs(a)**(1.0d0/dsub)
b0     = b**(1.0d0/dsub)

if (a .lt. 0 .AND. b .gt. 0) then
nints0 = 2
allocate(ab0(2,nints0))
ab0(1,1) = a0
ab0(2,1) = 0.0d0

ab0(1,2) = 0.0d0
ab0(2,2) = b0

elseif (b .gt. 0) then
nints0 = 1
allocate(ab0(2,nints0))
ab0(1,1) = 0.0d0
ab0(2,1) = b0
else
nints0 = 1
allocate(ab0(2,nints0))
ab0(1,1) = a0
ab0(2,1) = 0.0d0
endif
else
nints0 = 1
allocate(ab0(2,nints0))
ab0(1,1) = a0
ab0(2,1) = b0
endif

!
!  Discretize the user-supplied collection of functions
!

allocate(wrapdata)
wrapdata%userptr  = userptr
wrapdata%userfun => fun
wrapdata%ifscale  = 1
wrapdata%dsub     = dsub
wrapptr           = c_loc(wrapdata)


call legepw_init(disc,nlege,nints0,ab0)
call legepw_adap(epsdisc,nfuns,ggwrapper,disc,wrapptr)
call legepw_quad(disc,nquad0,xs0,whts0)

if (ifsuppress .ne. 1) then
call prini("in ggquad, nints = ",disc%nints)
call prini("in ggquad, nquad0 = ",nquad0)
endif

!call prin2("in ggquad, ab = ",disc%ab)
!
!  Orthonormalize the input functions
!

allocate(vals0(nquad0,nfuns))
call ggwrapper(nfuns,nquad0,xs0,whts0,vals0,wrapptr)
call qrdecomp(epsqr,nquad0,nfuns,vals0,krank,ipivs,vals,r)
allocate(rnorms(krank))
do i=1,krank
rnorms(i) = abs(r(i,i))
end do

if (ifsuppress .ne. 1) then

call prini("in ggquad, nfuns = ",nfuns)
call prini("in ggquad, krank = ",krank)
call prin2("in ggquad, rnorms = ",rnorms)
endif


!
! Construct the oversampled "Chebyshev" quadrature
!
call chebquad(disc,krank,vals,nquad,us,uwhts)
if (ifsuppress .ne. 1) then
call prin2("in ggquad, chebxs = ",us)
call prin2("in ggquad, chebwhts = ",uwhts)
endif

!
!  Perform Newton iterations
!
iorder=1
call gaussquad(epsnewt,iorder,disc,krank,vals,nquad,us,uwhts)

!
!  Apply the substitution 
!

allocate(xs(nquad), whts(nquad))
xs   = us**dsub
whts = uwhts * dsub*us**(dsub-1)

if (ifsuppress .ne. 1) then
call prini("in ggquad, nfuns = ",nfuns)
call prini("in ggquad, krank = ",krank)
call prini("in ggquad, final nquad = ",nquad)
call prin2("in ggquad, xs = ",xs)
call prin2("in ggquad, whts = ",whts)
endif

end subroutine


subroutine ggquad_prod(ifsing,a,b,nfuns1,nfuns2,fun1,fun2,userptr1,userptr2, &
  nquad,xs,whts)
implicit double precision (a-h,o-z)
procedure(ggfun)                           :: fun1, fun2
double precision, allocatable, intent(out) :: xs(:), whts(:)
type(c_ptr)                                :: userptr1, userptr2
!
!   Construct a generalized Gaussian quadrature for a collection of functions
!   of the form
!
!        f_i(x) g_j(x)      i=1,...,n,  j=1,..,m
!
!   where the f_i and g_j are supplied by the user via an external subroutines.
!
!  Input parameters:
!    ifsing - an integer parameter which, when set to 1, indicates that the input 
!      functions are given on an interval of the form (a,b) with a <= 0 <= b
!      and that the f_i(x) might have a mild sinuglarities at 0.
!
!    (a,b) - the interval on which the input functions are given
!    nfuns1 - the number of f_i's
!    nfuns2 - the number of g_j's
!    fun1 - an external subroutine conforming to the ggfun interface which
!      supplies the values of the f_i's 
!    fun2 - an external subroutien conforming to the ggfun interface which
!      supplies the values of the g_j's
!    userptr1 - a user-supplied "void *" pointer which is passed to the subroutines
!      fun1
!    userptr2 - a user-supplied "void *" pointer which is passed to the subroutines
!      fun2
!
!  Output parameters:
!    (nquad,xs,whts) - the resulting quadrature rule
!

type(legepw_disc)                   :: disc
type(c_ptr)                         :: wrapptr
type(ggwrapper_data), pointer       :: wrapdata

double precision, allocatable       :: vals10(:,:), vals1(:,:), rnorms1(:)
double precision, allocatable       :: vals20(:,:), vals2(:,:), rnorms2(:)
double precision, allocatable       :: prods0(:,:), prods(:,:), rnorms(:), r(:,:)

double precision, allocatable       :: ab0(:,:)
double precision, allocatable       :: xs0(:), whts0(:)
double precision, allocatable       :: us0(:), uwhts0(:)

double precision, allocatable       :: us(:), uwhts(:)
integer, allocatable                :: ipivs(:)

eps0     = epsilon(0.0d0)

epsdisc  = 1.0d-13
epsqr    = 1.0d-13
epsnewt  = 1.0d-7

if (eps0 .lt. 1.0d-17) then
epsdisc  = 1.0d-15
epsqr    = 1.0d-15
epsnewt  = 1.0d-7
endif

if (eps0 .lt. 1.0d-30) then
epsdisc  = 1.0d-30
epsqr    = 1.0d-30
epsnewt  = 1.0d-15
endif


nlege = 30
dsub  = 1
a0    = a
b0    = b


if (ifsing .eq. 1) then
dsub   = 5
a0     = sign(1.0d0,a)*abs(a)**(1.0d0/dsub)
b0     = b**(1.0d0/dsub)

if (a .lt. 0 .AND. b .gt. 0) then
nints0 = 2
allocate(ab0(2,nints0))
ab0(1,1) = a0
ab0(2,1) = 0.0d0
ab0(1,2) = 0.0d0
ab0(2,2) = b0
elseif (b .gt. 0) then
nints0 = 1
allocate(ab0(2,nints0))
ab0(1,1) = 0.0d0
ab0(2,1) = b0
else
nints0 = 1
allocate(ab0(2,nints0))
ab0(1,1) = a0
ab0(2,1) = 0.0d0
endif
else
nints0 = 1
allocate(ab0(2,nints0))
ab0(1,1) = a0
ab0(2,1) = b0
endif

!
!  Discretize the two user-supplied collections of functions
!

call legepw_init(disc,nlege,nints0,ab0)

allocate(wrapdata)

wrapdata%dsub     = dsub
wrapptr           = c_loc(wrapdata)

wrapdata%ifscale  = 0
wrapdata%userfun => fun2
wrapdata%userptr  = userptr2
call legepw_adap(epsdisc,nfuns2,ggwrapper,disc,wrapptr)
call prini("in ggquad_prod, nints = ",disc%nints)

wrapdata%ifscale  = 1
wrapdata%userfun => fun1
wrapdata%userptr  = userptr1
call legepw_adap(epsdisc,nfuns1,ggwrapper,disc,wrapptr)
call prini("in ggquad_prod, nints = ",disc%nints)

call legepw_order(disc,nlege*2)
call legepw_quad(disc,nquad0,xs0,whts0)
call prini("in ggquad_prod, nquad0 = ",nquad0)
!call prin2("in ggquad_prod, ab = ",disc%ab)

!
!  Orthonormalize each collection of input functions separately
!

allocate(vals10(nquad0,nfuns1), vals20(nquad0,nfuns2) )

wrapdata%userfun => fun1
wrapdata%userptr  = userptr1
wrapdata%ifscale  = 1
call ggwrapper(nfuns1,nquad0,xs0,whts0,vals10,wrapptr)

wrapdata%userfun => fun2
wrapdata%userptr  = userptr2
wrapdata%ifscale  = 0
call ggwrapper(nfuns2,nquad0,xs0,whts0,vals20,wrapptr)

call qrdecomp(epsqr,nquad0,nfuns1,vals10,krank1,ipivs,vals1,r)
allocate(rnorms1(krank1))
do i=1,krank1
rnorms1(i) = abs(r(i,i))
end do

if (ifsuppress .ne. 1) then
call prini("in ggquad_prod, nfuns1 = ",nfuns1)
call prini("in ggquad_prod, krank1 = ",krank1)
call prin2("in ggquad_prod, rnorms1 = ",rnorms1)
endif

call qrdecomp(epsqr,nquad0,nfuns2,vals20,krank2,ipivs,vals2,r)
allocate(rnorms2(krank2))
do i=1,krank2
rnorms2(i) = abs(r(i,i))
end do

if (ifsuppress .ne. 1) then
call prini("in ggquad_prod, nfuns2 = ",nfuns2)
call prini("in ggquad_prod, krank2 = ",krank2)
call prin2("in ggquad_prod, rnorms2 = ",rnorms2)
endif

!
!  Scale each set of functions by their "sinuglar values"
!
do j=1,krank1
vals1(:,j) = vals1(:,j) * rnorms1(j)
end do

do j=1,krank2
vals2(:,j) = vals2(:,j) * rnorms2(j)
end do

!
!  Construct the products and orthonormalize 'em
!

nprods = krank1*krank2
allocate(prods0(nquad0,nprods) )
do i=1,krank1
do j=1,krank2
idx = (i-1)*krank2
prods0(:,idx+j) = vals1(:,i)*vals2(:,j)*1/sqrt(whts0)
end do
end do

call qrdecomp(epsqr,nquad0,nprods,prods0,krank,ipivs,prods,r)
allocate(rnorms(krank))
do i=1,krank
rnorms(i) = abs(r(i,i))
end do

if (ifsuppress .ne. 1) then
call prini("in ggquad_prod, krank = ",krank)
call prin2("in ggquad_prod, rnorms = ",rnorms)
endif


!
! Construct the oversampled "Chebyshev" quadrature
!

call chebquad(disc,krank,prods,nquad,us,uwhts)

!
!  Perform Newton iterations
!
iorder=1
call gaussquad(epsnewt,iorder,disc,krank,prods,nquad,us,uwhts)

!
!  Apply the substitution (if any)
!



allocate(xs(nquad), whts(nquad))
xs   = us**dsub
whts = uwhts * dsub*us**(dsub-1)

if (ifsuppress .ne. 1) then
call prini("in ggquad_prod, nprods = ",nprods)
call prini("in ggquad_prod, krank  = ",krank)
call prini("in ggquad_prod, nquad  = ",nquad)
call prin2("in ggquad_prod, xs    = ",xs)
call prin2("in ggquad_prod, whts  = ",whts)
call prina("")
endif

end subroutine

subroutine ggquad_interp(a,b,nfuns,fun,userptr,nquad,xs,whts,rcond)
implicit double precision (a-h,o-z)
procedure(ggfun)                           :: fun
double precision, allocatable, intent(out) :: xs(:), whts(:)
type(c_ptr)                                :: userptr
!
!  Use Kirill Serkh's procedure to construct am interpolation quadrature for a 
!  collection of user-supplied smooth (or very mildly singular) functions 
!
!  Input parameters:
!    (a,b) - the domain of definition of the functions
!    nfuns - the number of input functions
!    fun - an external subroutine conforming to the ggfun interface which
!      supplies the values of the input functions
!    userptr - a user-supplied "void *" pointer which is passed to the subroutines
!      fun
!
!  Output parameters:
!    (nquad,xs,whts) - the resulting quadrature rule
!    rcond - an upper bound on the 
!

type(legepw_disc)                   :: disc
type(c_ptr)                         :: wrapptr
type(ggwrapper_data), pointer       :: wrapdata

double precision, allocatable       :: ab0(:,:)
double precision, allocatable       :: xs0(:), whts0(:)
double precision, allocatable       :: vals0(:,:), vals(:,:), rnorms(:), r(:,:)
double precision, allocatable       :: prods0(:,:), prods(:,:), rnormsprods(:)
integer, allocatable                :: ipivs(:), idxs(:)

double precision, allocatable       :: amatr(:,:)


eps0     = epsilon(0.0d0)

epsdisc  = 1.0d-14
epsqr    = 1.0d-14
epsnewt  = 1.0d-7

if (eps0 .lt. 1.0d-17) then
epsdisc  = 1.0d-15
epsqr    = 1.0d-15
epsnewt  = 1.0d-9
endif

if (eps0 .lt. 1.0d-30) then
epsdisc  = 1.0d-30
epsqr    = 1.0d-30
epsnewt  = 1.0d-15
endif

nlege = 30
dsub  = 1
call legepw_init(disc,nlege,a,b)

allocate(wrapdata)

wrapdata%dsub     = dsub
wrapptr           = c_loc(wrapdata)
wrapdata%ifscale  = 0
wrapdata%userfun => fun
wrapdata%userptr  = userptr

call legepw_adap(epsdisc,nfuns,ggwrapper,disc,wrapptr)
call prini("in ggquad_interp, nints = ",disc%nints)

!call legepw_order(disc,nlege*2)
call legepw_quad(disc,nquad0,xs0,whts0)
call prini("in ggquad_interp, nquad0 = ",nquad0)
!call prin2("in ggquad_interp, ab = ",disc%ab)

allocate(vals0(nquad0,nfuns))
call ggwrapper(nfuns,nquad0,xs0,whts0,vals0,wrapptr)
call qrdecomp(epsqr,nquad0,nfuns,vals0,krank,ipivs,vals,r)

allocate(rnorms(krank))
do i=1,krank
rnorms(i) = abs(r(i,i))
end do

do j=1,krank
dd = rnorms(j)

if (dd .lt. epsnewt)  then
krank = j-1
exit
endif

vals(:,j) = vals(:,j)  * sqrt(dd)
end do


call prini("in ggquad_interp, nfuns  = ",nfuns)
call prini("in ggquad_interp, krank  = ",krank)
call prin2("in ggquad_interp, rnorms = ",rnorms(1:krank))

!
!  Construct the products and orthonormalize them
!

nprods = krank*krank
allocate(prods0(nquad0,nprods) )
do i=1,krank
do j=1,krank
idx             = (i-1)*krank
prods0(:,idx+j) = vals(:,i)*vals(:,j)/sqrt(whts0)
end do
end do

eps = epsnewt/10
call qrdecomp(eps,nquad0,nprods,prods0,krankprods,ipivs,prods,r)

allocate(rnormsprods(krankprods))
do i=1,krankprods
dd             = abs(r(i,i))
rnormsprods(i) = dd
if (dd .lt. epsnewt) then
krankprods = i-1
exit
endif
end do

!krankprods = min(2*krank,krankprods)
call prini("in ggquad_interp, krankprods = ",krankprods)
call prin2("in ggquad_interp, rnormsprods = ",rnormsprods(1:krankprods))


!
! Construct the oversampled "Chebyshev" quadrature
!

call chebquad(disc,krankprods,prods(:,1:krankprods),nquad,xs,whts)
call prin2("in ggquad_interp, xs = ",xs)
call prin2("in ggquad_interp, whts = ",whts)


!
!  Perform Newton iterations to reduce the quadrature rule
!

iorder=1
call gaussquad(epsnewt,iorder,disc,krankprods,prods(:,1:krankprods),nquad,xs,whts)

!
! Apply the substitution 
!

call prini("in ggquad_interp, krank = ",krank)
call prini("in ggquad_interp, krankprods = ",krankprods)
call prini("in ggquad_interp, nquad  = ",nquad)
call prin2("in ggquad_interp, xs = ",xs)
call prin2("in ggquad_interp, whts = ",whts)


!
!  Bound the condition number of interpolation
!

do j=1,krank
vals(:,j) = vals(:,j)/sqrt(rnorms(j))
end do

allocate(amatr(nquad,krank))

do i=1,nquad
x   = xs(i)
wht = whts(i)

call legepw_interp(disc,vals,x,amatr(i,:))
amatr(i,:) = amatr(i,:) *sqrt(wht)
end do

dd    = norm2(matmul(transpose(amatr),amatr)-eye(krank))
rcond = sqrt( (1+dd**2) /(1-dd**2) )

call prind("in ggquad_interp, ||A*A - I||_2 = ",dd)
call prind("in ggquad_interp, rcond = ",rcond)
call prina("")

deallocate(amatr)

end subroutine


subroutine ggquad_interp_sing(nfuns,fun,userptr,nquad,xs,whts,rcond)
implicit double precision (a-h,o-z)
procedure(ggfun)                           :: fun
double precision, allocatable, intent(out) :: xs(:), whts(:)
type(c_ptr)                                :: userptr
!
!  Use Kirill Serkh's procedure to construct am interpolation quadrature for a 
!  collection of user-supplied functions given on the interval (-1,1) and which
!  can have singularities at 0
!
!  Input parameters:
!    nfuns - the number of input functions
!    fun - an external subroutine conforming to the ggfun interface which
!      supplies the values of the input functions
!    userptr - a user-supplied "void *" pointer which is passed to the subroutines
!      fun
!
!  Output parameters:
!    (nquad,xs,whts) - the resulting quadrature rule
!    rcond - an upper bound on the "condition number of the basis"
!

type(legepw_disc)                   :: disc
type(c_ptr)                         :: wrapptr
type(ggwrapper_data), pointer       :: wrapdata

double precision, allocatable       :: ab0(:,:)
double precision, allocatable       :: xs0(:), whts0(:)
double precision, allocatable       :: vals0(:,:), vals(:,:), rnorms(:), r(:,:)
double precision, allocatable       :: prods0(:,:), prods(:,:), rnormsprods(:)
double precision, allocatable       :: us(:), uwhts(:)
integer, allocatable                :: ipivs(:), idxs(:)

double precision, allocatable       :: amatr(:,:), rints(:), rints0(:)


eps0     = epsilon(0.0d0)

epsdisc  = 1.0d-14
epsqr    = 1.0d-14
epsnewt  = 1.0d-10

if (eps0 .lt. 1.0d-17) then
epsdisc  = 1.0d-17
epsqr    = 1.0d-17
epsnewt  = 1.0d-10
endif

if (eps0 .lt. 1.0d-30) then
epsdisc  = 1.0d-30
epsqr    = 1.0d-30
epsnewt  = 1.0d-20
endif

nlege = 30

!
!  First construct a quadratue rule for the 
!

nints0 = 4
allocate(ab0(2,nints0))

ab0(1,1) = -70.0d0
ab0(2,1) =  -1.0d0

ab0(1,2) =  -1.0d0
ab0(2,2) =   0.0d0

ab0(1,3) =   0.0d0
ab0(2,3) =   1.0d0

ab0(1,4) =   1.0d0
ab0(2,4) =  70.0d0


call legepw_init(disc,nlege,nints0,ab0)

allocate(wrapdata)

wrapdata%dsub     = dsub
wrapptr           = c_loc(wrapdata)
wrapdata%userfun => fun
wrapdata%userptr  = userptr

call legepw_adap(epsdisc,nfuns,ggwrapper3,disc,wrapptr)
call prini("in ggquad_interp_sing, nints = ",disc%nints)

! promoting the order is unnecessary
!call legepw_order(disc,nlege*2)
call legepw_quad(disc,nquad0,xs0,whts0)
call prini("in ggquad_interp_sing, nquad0 = ",nquad0)
call prin2("in ggquad_interp_sing, ab = ",disc%ab)

allocate(vals0(nquad0,nfuns))
call ggwrapper3(nfuns,nquad0,xs0,whts0,vals0,wrapptr)
call qrdecomp(epsqr,nquad0,nfuns,vals0,krank,ipivs,vals,r)

allocate(rnorms(krank))
do i=1,krank
rnorms(i) = abs(r(i,i))
end do

do j=1,krank
dd = rnorms(j)

if (dd .lt. epsnewt)  then
krank = j-1
exit
endif

vals(:,j) = vals(:,j)  * sqrt(dd)
end do

call prini("in ggquad_interp_sing, nfuns  = ",nfuns)
call prini("in ggquad_interp_sing, krank  = ",krank)
call prin2("in ggquad_interp_sing, rnorms = ",rnorms(1:krank))

!
!  Construct the products and orthonormalize them
!

nprods = krank*krank
allocate(prods0(nquad0,nprods) )
do i=1,krank
do j=1,krank
idx             = (i-1)*krank
prods0(:,idx+j) = vals(:,i)*vals(:,j)/sqrt(whts0)
end do
end do

eps = epsnewt/100
call qrdecomp(eps,nquad0,nprods,prods0,krankprods,ipivs,prods,r)

allocate(rnormsprods(krankprods))
do i=1,krankprods
dd             = abs(r(i,i))
rnormsprods(i) = dd

if (dd .lt. epnewt) then
krankprods = i-1
exit
endif

end do

call prini("in ggquad_interp_sing, krankprods = ",krankprods)
call prin2("in ggquad_interp_sing, rnormsprods = ",rnormsprods(1:krankprods))

!
! Construct the oversampled "Chebyshev" quadrature
!

call chebquad(disc,krankprods,prods(:,1:krankprods),nquad,us,uwhts)

!
!  Perform Newton iterations to reduce the quadrature rule
!

iorder=1
call gaussquad(epsnewt,iorder,disc,krankprods,prods(:,1:krankprods),nquad,us,uwhts)

!
! Apply the substitution and sort the nodes
!

allocate(xs(nquad), whts(nquad), idxs(nquad))
xs   = exp(-abs(us))* sign(1.0d0,us)
whts = exp(-abs(us))*uwhts

call insort2(nquad,xs,idxs)
whts = whts(idxs)

call prini("in ggquad_interp_sing, krank = ",krank)
call prini("in ggquad_interp_sing, krankprods = ",krankprods)
call prini("in ggquad_interp_sing, nquad  = ",nquad)
call prin2("in ggquad_interp_sing, xs = ",xs)
call prin2("in ggquad_interp_sing, whts = ",whts)

!
!  Bound the condition number of interpolation
!

do j=1,krank
vals(:,j) = vals(:,j)/sqrt(rnorms(j))
end do

allocate(amatr(nquad,krank))

do i=1,nquad
x   = us(i)
wht = uwhts(i)

call legepw_interp(disc,vals,x,amatr(i,:))
amatr(i,:) = amatr(i,:) *sqrt(wht)
end do

dd    = norm2(matmul(transpose(amatr),amatr)-eye(krank))
rcond = sqrt(1+dd**2)/sqrt(1-dd**2)
deallocate(amatr)

call prind("in ggquad_interp_sing, rcond = ",rcond)

end subroutine


subroutine ggquad_sing(nfuns,fun,userptr,nquad,xs,whts)
implicit double precision (a-h,o-z)
procedure(ggfun)                           :: fun
double precision, allocatable, intent(out) :: xs(:), whts(:)
type(c_ptr)                                :: userptr
!
!   Construct a generalized Gaussian quadrature on the interval (-1,1) for 
!   a collection of the functions which are strongly sinuglar at 0.
!
!  Input parameters:
!    nfuns - the number of input functions
!    fun - an external subroutine conforming to the ggfun interface which
!      supplies the values of the input functions
!    userptr - a user-supplied "void *" pointer which is passed to the subroutine
!      fun
!
!  Output parameters:
!    (nquad,xs,whts) - the resulting quadrature rule
!

type(legepw_disc)                   :: disc
type(c_ptr)                         :: wrapptr
type(ggwrapper_data), pointer       :: wrapdata
double precision, allocatable       :: vals0(:,:),  vals(:,:), r(:,:), ab0(:,:)
double precision, allocatable       :: xs0(:), whts0(:)
double precision, allocatable       :: us(:), uwhts(:), rnorms(:)
integer, allocatable                :: ipivs(:), idxs(:)

eps0   = epsilon(0.0d0)

epsdisc  = 1.0d-13
epsqr    = 1.0d-13
epsnewt  = 1.0d-7

if (eps0 .lt. 1.0d-17) then
epsdisc  = 1.0d-17
epsqr    = 1.0d-17
epsnewt  = 1.0d-12
endif

if (eps0 .lt. 1.0d-30) then
epsdisc  = 1.0d-30
epsqr    = 1.0d-30
epsnewt  = 1.0d-15
endif

nlege  = 30
nints0 = 2
allocate(ab0(2,nints0))
ab0(1,1)   = -500.0d0
ab0(2,1)   = -1.0d0

ab0(1,2)   = -1.0d0
ab0(2,2)   =  0.0d0

ab0(1,3)   =  0.0d0
ab0(2,3)   =  1.0d0

ab0(1,4)   =  1.0d0
ab0(2,4)   =900.0d0


!  Discretize the user-supplied collection of functions
!

allocate(wrapdata)
wrapdata%ifscale  = 1
wrapdata%userptr  = userptr
wrapdata%userfun => fun
wrapptr           = c_loc(wrapdata)

call legepw_init(disc,nlege,nints0,ab0)
call legepw_adap(epsdisc,nfuns,ggwrapper2,disc,wrapptr)
call legepw_quad(disc,nquad0,xs0,whts0)

call prini("in ggquad_sing, nints = ",disc%nints)
call prini("in ggquad_sing, nquad0 = ",nquad0)

!
!  Orthonormalize the input functions
!

allocate(vals0(nquad0,nfuns))
call ggwrapper2(nfuns,nquad0,xs0,whts0,vals0,wrapptr)
call qrdecomp(epsqr,nquad0,nfuns,vals0,krank,ipivs,vals,r)
allocate(rnorms(krank))
do i=1,krank
rnorms(i) = abs(r(i,i))
end do

call prini("in ggquad_sing, krank = ",krank)
call prin2("in ggquad_sing, rnorms = ",rnorms)

!

! Construct the oversampled "Chebyshev" quadrature
!
call chebquad(disc,krank,vals,nquad,us,uwhts)


!
!  Perform Newton iterations
!
iorder=1
call gaussquad(epsnewt,iorder,disc,krank,vals,nquad,us,uwhts)

!
!  Apply the substitution 
!

allocate(xs(nquad), whts(nquad), idxs(nquad))
xs   = exp(-abs(us))* sign(1.0d0,us)
whts = exp(-abs(us))*uwhts


call insort2(nquad,xs,idxs)
whts = whts(idxs)

call prini("in ggquad_sing, nfuns = ",nfuns)
call prini("in ggquad_sing, krank = ",krank)
call prini("in ggquad_sing, nquad = ",nquad)
call prin2("in ggquad_sing, xs    = ",xs)
call prin2("in ggquad_sing, whts  = ",whts)
call prina("")

end subroutine



subroutine ggquad_singprod(nfuns1,nfuns2,fun1,fun2,userptr1,userptr2,nquad,xs,whts)
implicit double precision (a-h,o-z)
procedure(ggfun)                           :: fun1, fun2
double precision, allocatable, intent(out) :: xs(:), whts(:)
type(c_ptr)                                :: userptr1, userptr2
!
!   Construct a generalized Gaussian quadrature for evaluating integrals of the form
!
!        1
!    \int    f_i(x) g_j(x)  dx,       i=1,...,n,  j=1,..,m,
!       -1
! 
!   where the f_i(x) are singular at 0 and the g_j(x) are smooth.
!
!  Input parameters:
!    nfuns1 - the number of f_i's
!    nfuns2 - the number of g_j's
!    fun1 - an external subroutine conforming to the ggfun interface which
!      supplies the values of the f_i's 
!    fun2 - an external subroutien conforming to the ggfun interface which
!      supplies the values of the g_j's
!    userptr1 - a user-supplied "void *" pointer which is passed to the subroutines
!      fun1
!    userptr2 - a user-supplied "void *" pointer which is passed to the subroutines
!      fun2
!
!  Output parameters:
!    (nquad,xs,whts) - the resulting quadrature rule
!

type(legepw_disc)                   :: disc
type(c_ptr)                         :: wrapptr
type(ggwrapper_data), pointer       :: wrapdata

double precision, allocatable       :: vals10(:,:), vals1(:,:), rnorms1(:)
double precision, allocatable       :: vals20(:,:), vals2(:,:), rnorms2(:)
double precision, allocatable       :: prods0(:,:), prods(:,:), rnorms(:), r(:,:)

double precision, allocatable       :: xs0(:), whts0(:), ab0(:,:)
double precision, allocatable       :: us0(:), uwhts0(:)

double precision, allocatable       :: us(:), uwhts(:), dxs(:)
integer, allocatable                :: ipivs(:), idxs(:)

eps0     = epsilon(0.0d0)

epsdisc  = 1.0d-14
epsqr    = 1.0d-14
epsnewt  = 1.0d-10

if (eps0 .lt. 1.0d-17) then
epsdisc  = 1.0d-15
epsqr    = 1.0d-15
epsnewt  = 1.0d-12
endif

if (eps0 .lt. 1.0d-30) then
epsdisc  = 1.0d-30
epsqr    = 1.0d-30
epsnewt  = 1.0d-20
endif

nlege = 30

nints0 = 4
allocate(ab0(2,nints0))
ab0(1,1) = -100.0d0
ab0(2,1) = -1.0d0

ab0(1,2) = -1.0d0
ab0(2,2) = 0.0d0

ab0(1,3) = 0
ab0(2,3) = 1.0d0

ab0(1,4) = 1
ab0(2,4) = 100.0d0

!
!  Discretize the two user-supplied collections of functions
!

call legepw_init(disc,nlege,nints0,ab0)


allocate(wrapdata)
wrapptr           = c_loc(wrapdata)

wrapdata%ifscale  = 0
wrapdata%userfun => fun2
wrapdata%userptr  = userptr2
call legepw_adap(epsdisc,nfuns2,ggwrapper2,disc,wrapptr)
call prini("in ggquad_singprod, nints = ",disc%nints)


wrapdata%ifscale  = 1
wrapdata%userfun => fun1
wrapdata%userptr  = userptr1
call legepw_adap(epsdisc,nfuns1,ggwrapper2,disc,wrapptr)
call prini("in ggquad_singprod, nints = ",disc%nints)

call legepw_order(disc,nlege*2)
call legepw_quad(disc,nquad0,xs0,whts0)
call prini("in ggquad_singprod, nquad0 = ",nquad0)


!
!  Orthonormalize each collection of input functions separately
!

allocate(vals10(nquad0,nfuns1), vals20(nquad0,nfuns2) )


wrapdata%userfun => fun2
wrapdata%userptr  = userptr2
wrapdata%ifscale  = 0
call ggwrapper2(nfuns2,nquad0,xs0,whts0,vals20,wrapptr)

wrapdata%userfun => fun1
wrapdata%userptr  = userptr1
wrapdata%ifscale  = 1
call ggwrapper2(nfuns1,nquad0,xs0,whts0,vals10,wrapptr)

call qrdecomp(epsqr,nquad0,nfuns1,vals10,krank1,ipivs,vals1,r)
allocate(rnorms1(krank1))
do i=1,krank1
rnorms1(i) = abs(r(i,i))
end do
call prini("in ggquad_singprod, nfuns1 = ",nfuns1)
call prini("in ggquad_singprod, krank1 = ",krank1)
call prin2("in ggquad_singprod, rnorms1 = ",rnorms1)

call qrdecomp(epsqr,nquad0,nfuns2,vals20,krank2,ipivs,vals2,r)
allocate(rnorms2(krank2))
do i=1,krank2
rnorms2(i) = abs(r(i,i))
end do
call prini("in ggquad_singprod, nfuns2 = ",nfuns2)
call prini("in ggquad_singprod, krank2 = ",krank2)
call prin2("in ggquad_singprod, rnorms2 = ",rnorms2)

!
!  Construct the products and orthonormalize 'em
!



nprods = krank1*krank2
allocate(prods0(nquad0,nprods) )
do i=1,krank1
idx = (i-1)*krank2
do j=1,krank2
prods0(:,idx+j) = vals1(:,i)*vals2(:,j)*1/sqrt(whts0)*rnorms1(i)
end do
end do

call qrdecomp(epsqr,nquad0,nprods,prods0,krank,ipivs,prods,r)
allocate(rnorms(krank))
do i=1,krank
rnorms(i) = abs(r(i,i))
end do
call prini("in ggprod_singprod, krank = ",krank)
call prin2("in ggprod_singprod, rnorms = ",rnorms)


!
! Construct the oversampled "Chebyshev" quadrature
!
call chebquad(disc,krank,prods,nquad,us,uwhts)

!
!  Perform Newton iterations
!

iorder=1
call gaussquad(epsnewt,iorder,disc,krank,prods,nquad,us,uwhts)

!
!  Apply the substitution 
!

allocate(xs(nquad), whts(nquad), idxs(nquad))
do i=1,nquad
xs(nquad-i+1)   = exp(-abs(us(i)))*sign(1.0d0,us(i))
whts(nquad-i+1) = exp(-abs(us(i)))*uwhts(i)
end do

call insort2(nquad,xs,idxs)
whts = whts(idxs)

call prini("in ggprod_singprod, nprods = ",nprods)
call prini("in ggprod_singprod, krank  = ",krank)
call prini("in ggprod_singprod, nquad  = ",nquad)
call prin2("in ggprod_singprod, xs     = ",xs)
call prin2("in ggprod_singprod, whts   = ",whts)


end subroutine



subroutine ggwrapper(nfuns,n,us,uwhts,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: us(n)
double precision              :: uwhts(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
type(ggwrapper_data), pointer :: wrapdata
double precision, allocatable :: xs(:),dxs(:)

call c_f_pointer(userptr,wrapdata)

allocate(xs(n),dxs(n))

dsub              = wrapdata%dsub
ifscale           = wrapdata%ifscale
xs                = us**dsub
dxs               = dsub*us**(dsub-1)
call wrapdata%userfun(nfuns,n,xs,vals,wrapdata%userptr)

if (ifscale .eq. 1) then
do j=1,nfuns
vals(:,j) = vals(:,j) * dxs * sqrt(uwhts)
end do
else
do j=1,nfuns
vals(:,j) = vals(:,j) * sqrt(uwhts)
end do
endif


end subroutine


subroutine ggwrapper2(nfuns,n,us,uwhts,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: us(n)
double precision              :: uwhts(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
type(ggwrapper_data), pointer :: wrapdata

double precision, allocatable :: xs(:),dxs(:)

call c_f_pointer(userptr,wrapdata)

allocate(xs(n),dxs(n))


ifscale           = wrapdata%ifscale
xs                = exp(-abs(us))*sign(1.0d0,us)
dxs               = exp(-abs(us))


call wrapdata%userfun(nfuns,n,xs,vals,wrapdata%userptr)

if (ifscale .eq. 1) then

do j=1,nfuns
vals(:,j) = vals(:,j) * dxs * sqrt(uwhts)
end do

else

do j=1,nfuns
vals(:,j) = vals(:,j) * sqrt(uwhts)
end do

endif

end subroutine



subroutine ggwrapper3(nfuns,n,us,uwhts,vals,userptr)
implicit double precision (a-h,o-z)
double precision              :: us(n)
double precision              :: uwhts(n)
double precision              :: vals(:,:)
type(c_ptr)                   :: userptr
type(ggwrapper_data), pointer :: wrapdata
double precision, allocatable :: xs(:),whts(:),dxs(:)

call c_f_pointer(userptr,wrapdata)
allocate(xs(n),dxs(n),whts(n))


dsub              = wrapdata%dsub
ifscale           = wrapdata%ifscale

xs                = exp(-abs(us))*sign(1.0d0,us)
whts              = exp(-abs(us))*uwhts

call wrapdata%userfun(nfuns,n,xs,vals,wrapdata%userptr)

do j=1,nfuns
vals(:,j) = vals(:,j) * sqrt(whts) 
end do

end subroutine


end module
