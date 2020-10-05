!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for constructing "generalized Gaussian" and ``generalized
!  Chebyshev'' quadrature rules for collections of bivariate functions.
!
!  A generalized Chebyshev quadrature rule for a linearly independent collection of 
!  3n bivariate functions f_1,...,f_3n given on a region Omega is an quadrature rule 
!  of the form
!
!                                          3n
!    \int  \int   f(x, y) dx dy  \approx  \sum  f(x_j,y_j) w_j                               (1)
!        Omega                             j=1
!
!  which is exact for the functions f_1, ..., f_3n (or, at least, the formula (1) 
!  achieves specified accuracy when one of the f_j is substituted for f). 

!  A generalized Gaussian quadrature rule for a linearly independent 
!  collection of 3n bivariate functions f_1,...,f_3n is an quadrature rule of the form
!
!                                           n
!    \int  \int   f(x, y) dx dy  \approx  \sum  f(x_j,y_j) w_j                               (2)
!       Omega                              j=1
!
!  which is exact for the functions f_1, ..., f_3n (or, at least, the formula (2) 
!  achieves specified accuracy when one of the f_j is substituted for f). 
!
!  The following subroutines should be regarded as publicly callable:
!
!    chebquad2d - construct a generalized Chebyshev quadrature for a collection of
!      user-supplied input functions which are specified using an external subroutine
!
!    gaussquad2d - by applying Newton's method the obvious nonlinear system of 
!      equations satisfied by the rule (2), downsample an existing quadrature point-by-
!      point in the hopes of creating a generalized Gaussian quadrature rule for a 
!      collection of user-supplied functions which are specified via an external 
!      subroutine
!
!    chebquad2d_pw - construct a generalized Chebyshev quadrature rule for a collection
!      of functions specified as piecewise bivariate Legendre expansions 
!      (ala bilege_pw.f90)
!   
!    gaussquad2_pw - construct a generalized Gaussian quadrature rule for a collection of
!      functions specified as piecewise bivariate Legendre expansions 
!      (ala bilege_pw.f90) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gaussian_quad2d

use utils
use linalg
use bilege
use bilegepw
use iso_c_binding

interface

subroutine gaussfun2d(nfuns,n,xs,ys,vals,dersx,dersy,userptr)
import c_ptr
implicit double precision (a-h,o-z)
type(c_ptr)           :: userptr
double precision      :: xs(:),ys(:)
double precision      :: vals(:,:), dersx(:,:), dersy(:,:)
end subroutine

subroutine gaussfun2d_region(x,y,wht,isvalid,userptr)
import c_ptr
implicit double precision (a-h,o-z)
type(c_ptr)           :: userptr
double precision      :: x,y,wht
integer               :: isvalid
end subroutine

end interface

type      gq2dwrapper_data
type(bilegepw_disc), pointer  :: disc
double precision, allocatable :: coefs(:,:)
end type  gq2dwrapper_data

contains


subroutine chebquad2d(nfuns,funuser,userptr,nquad,xs,ys,whts)
implicit double precision (a-h,o-z)
procedure(gaussfun2d)                        :: funuser
type(c_ptr)                                  :: userptr
double precision, allocatable, intent(inout) :: xs(:), ys(:),whts(:)
!
!  Construct a generalized Chebyshev quadrature rule for a collection of user-supplied
!  functions represented via piecewise Legendre expansions.  This routine takes as
!  input the *scaled* values of the input functions at the nodes of the 
!  quadrature associated with a piecewise Legendre discretization scheme (see legepw.f90).
!
!  The size of the resulting quadrature rule will be equal to the numerical dimension
!  of the space they space.
!
!  Input parameters:
!    disc - a data structure describing the piecewise Legendre scheme used to represent
!     the input functions
!    nfuns - the number of input functions
!    funuser - an external subroutine supplying the values and first order derivatives
!      of the input functions
!
!  Output parameters: 
!    nquad - the number of nodes in the resulting quadrature, which will be 
!      the numerical dimension of the space spanned by the input functions
!    xs - an array specifying the nodes of the resulting quadrature rule
!    whts - an array specifying the weights of the resulting quadrature rule
!

double precision, allocatable :: rints(:), vals(:,:), dersx(:,:), dersy(:,:)
double precision, allocatable :: xsnew(:), ysnew(:), whtsnew(:)
double precision, allocatable :: a(:,:), b(:,:), r(:,:), rnorms(:), q(:,:)
integer, allocatable          :: ipivs(:), idxs(:)

eps0  = epsilon(0.0d0)

allocate(xsnew(nfuns),ysnew(nfuns),whtsnew(nfuns),rints(nfuns))

!
!  Evaluate the user-supplied function at the quadrature nodes
!
allocate(vals(nquad,nfuns),dersx(nquad,nfuns),dersy(nquad,nfuns))
call funuser(nfuns,nquad,xs,ys,vals,dersx,dersy,userptr)

do i=1,nquad
vals(i,:) = vals(i,:)*whts(i)
end do
rints = sum(vals,1)
do i=1,nquad
vals(i,:) = vals(i,:)/sqrt(whts(i))
end do

!
!  Form the weighted transpose
!

allocate(a(nfuns,nquad))
a = transpose(vals)
call intdecomp(eps0,nfuns,nquad,a,krank,ipivs,r)

!
!  Subsample the matrix and solve for the quadrature weights
!

nquadnew = krank
allocate(b(nfuns,krank))
b = a(:,ipivs)

call leastsq(nfuns,krank,b,rints,whtsnew)
xsnew   = xs(ipivs)
ysnew   = ys(ipivs)
whtsnew = whtsnew * sqrt(whts(ipivs))

!
!  Sort the quadrature nodes
!

allocate(idxs(krank))
call insort2(krank,xsnew,idxs)

deallocate(xs,ys,whts)
nquad = krank

allocate(xs(nquad),ys(nquad),whts(nquad))
xs   = xsnew
ys   = ysnew(idxs)
whts = whtsnew(idxs)


end subroutine


subroutine gaussquad2d(eps,nfuns,funuser,funregion,userptr,nquad,xs,ys,whts)
implicit double precision (a-h,o-z)
procedure(gaussfun2d)                                :: funuser
procedure(gaussfun2d_region)                         :: funregion
double precision, allocatable, intent(inout)         :: xs(:), ys(:), whts(:)
type(c_ptr)                                          :: userptr
!
!  Starting with an existing generalized Chebyshev quadrature rule for a collection 
!  of n user-supplied functions, attempt to reduce the formula point-by-point using 
!  Newton iterations in the hopes of constructing a generalized Gaussian quadrature 
!  rule for the collection.
!
!  Input parameters:
!    eps - the desired precision for the desired quadrature rule
!    nfuns - the number of input functions
!    funuser - a user-supplied external subroutine which conforms to the 
!      gaussfun2d interface and evaluates the input functions and their
!      derivatives  
!    funregion - an external subroutine which indicates which points are within
!      the integration domain
!    userptr - a "void *" pointer which is passed along to funuser and funregion
!    (nquad, xs, ys, whts) - 
!
!  Output parameters: 
!    (nquad, xs, ys, whts) - the resulting quadrature rule
!

double precision, allocatable   :: rints(:), signifs(:), rints0(:)
integer, allocatable            :: idxs(:)
double precision, allocatable   :: xsnew(:), ysnew(:), whtsnew(:)


eps0 = epsilon(0.0d0)

!
!  Compute the integrals of the input functions
!
allocate(rints(nfuns), rints0(nfuns) )
call gaussquad2d_rints(nfuns,funuser,userptr,nquad,xs,ys,whts,rints)
call prin2("in gaussquad2d, rints = ",rints)

!
!  Attempt to remove 2*nquad/3 points
!


do nremove=1,2*nquad/3
call elapsed(time2)
write (*,"(A,I3,A,I3,A)") "[--- gaussquad2d: nfuns = ",nfuns,", nquad = ",nquad," ---------------------------]"

!
!  Reorder the nodes of the existing quadrature rule as specified by the user
!  and attempt to remove each point in that order
!

allocate(signifs(nquad), idxs(nquad) )
allocate(xsnew(nquad), ysnew(nquad), whtsnew(nquad))
call gaussquad2d_signifs(nfuns,funuser,userptr,rints,nquad,xs,ys,whts,signifs)
!call prin2("in gaussquad2d, signifs = ",signifs)

call insort2(nquad,signifs,idxs)
xsnew   = xs(idxs)
ysnew   = ys(idxs)
whtsnew = whts(idxs)
xs      = xsnew
ys      = ysnew
whts    = whtsnew


call gaussquad2d_rints(nfuns,funuser,userptr,nquad,xs,ys,whts,rints0)
errl2 = norm2(rints-rints0)
deallocate(signifs, idxs, xsnew, ysnew, whtsnew)


nquadnew = nquad-1
allocate(xsnew(nquadnew), ysnew(nquadnew), whtsnew(nquadnew))

do ipt=1,nquad


write (*,"(A,I3,A)") "attempting to remove ipt = ",ipt
call elapsed(t1)

!
!  Remove ipt from the quadrature rule, and solve the obvious nonlinear system of equations via 
!  Newton's method to attempt to improve the accuracy of the resulting rule
!

idx = 0  
do i=1,nquad
if (i .eq. ipt) cycle
idx           = idx+1
xsnew(idx)    = xs(i)
ysnew(idx)    = ys(i)
whtsnew(idx)  = whts(i)
end do

call gaussquad2d_newton(ier,eps,nfuns,funuser,funregion,userptr,rints,nquadnew,xsnew,ysnew,whtsnew,errl2)

if (ier .eq. 0) exit
end do


call elapsed(t2)

! if we have removed a point, update the quadrature 
if (ier .eq. 0 ) then

deallocate(xs,ys,whts)
allocate(xs(nquadnew), ys(nquadnew), whts(nquadnew))
nquad = nquadnew
xs    = xsnew
ys    = ysnew
whts  = whtsnew
deallocate(xsnew,ysnew,whtsnew)

call prin2("xs = ", xs)
call prin2("ys = ", ys)
call prin2("whts = ", whts)

write (*,"(A,F9.2,A)") "successful, time = ",t2-t1, " seconds "
write (*,*) ""


else
write (*,"(A,F9.2,A)") "failed,     time = ",t2-t1, " seconds "
write (*,*) ""

! otherwise, we have the best formula we could find already
deallocate(xsnew,ysnew,whtsnew)
exit
endif

end do

 
allocate(ysnew(nquad),whtsnew(nquad),idxs(nquad))
ysnew   = ys
whtsnew = whts

call insort2(nquad,xs,idxs)
ys   = ysnew(idxs)
whts = whtsnew(idxs)

end subroutine



subroutine chebquad2d_pw(disc,nfuns,vals,nquad,xs,ys,whts)
implicit double precision (a-h,o-z)
type(bilegepw_disc), target                  :: disc
double precision                             :: vals(:,:)
double precision, allocatable, intent(out)   :: xs(:), ys(:),whts(:)
!
!  Construct a generalized Chebyshev quadrature rule for a collection of user-supplied
!  functions represented via piecewise Legendre expansions.  This routine takes as
!  input the *scaled* values of the input functions at the nodes of the 
!  quadrature associated with a piecewise Legendre discretization scheme (see legepw.f90).
!
!  The size of the resulting quadrature rule will be equal to the numerical dimension
!  of the space they space.
!
!  Input parameters:
!    disc - a data structure describing the discretizations cheme used to represent
!     the input functions
!    nfuns - the number of input functions
!    vals - the matrix whose jth column gives the scaled values of the jth input
!      function at the nodes of the discretization quadrature     
!
!  Output parameters: 
!    (nquad,xs,ys,whts) - the newly constructed quadrature rule
!

double precision, allocatable :: xs0(:), ys0(:), whts0(:), rints(:)
double precision, allocatable :: a(:,:), b(:,:), r(:,:)
integer, allocatable          :: ipivs(:), idxs(:)

eps0  = epsilon(0.0d0)


allocate(rints(nfuns),xs(nfuns),ys(nfuns),whts(nfuns))

!
!  Fetch the oversampled quadrature rule
!
call bilegepw_quad(disc,nquad0,xs0,ys0,whts0)

!
!  Compute the integrals of the input functions
!

do i=1,nquad0
vals(i,:) = vals(i,:)*sqrt(whts0(i))
end do
!rints = sum(vals(i1:i2,:),1)
rints = sum(vals,1)
do i=1,nquad0
vals(i,:) = vals(i,:)/sqrt(whts0(i))
end do



!
!  Form the weighted transpose
!

allocate(a(nfuns,nquad0))
a = transpose(vals)
call intdecomp(eps0,nfuns,nquad0,a,krank,ipivs,r)

!
!  Subsample the matrix and solve for the quadrature weights
!

nquad = krank
allocate(b(nfuns,krank))
b = a(:,ipivs)

call leastsq(nfuns,krank,b,rints,whts)
xs   = xs0(ipivs)
ys   = ys0(ipivs)
whts = whts * sqrt(whts0(ipivs))

!
!  Sort the quadrature nodes
!

allocate(idxs(krank))
call insort2(krank,xs,idxs)
ys   = ys(idxs)
whts = whts(idxs)

end subroutine




subroutine gaussquad2_pw(eps,disc,nfuns,vals,nquad,xs,ys,whts)
implicit double precision (a-h,o-z)
type(bilegepw_disc), target                          :: disc
double precision                                     :: vals(:,:)
double precision, allocatable, intent(inout)         :: xs(:), ys(:), whts(:)
!
!  Starting with an existing quadrature rule for a collection of n user-supplied
!  functions, attempt to reduce the formula point-by-point using Newton iterations
!  in the hopes of constructing a generalized Gaussian quadrature rule for the
!  collection.  The input functions are specified via their values 
!
!  Input parameters:
!    eps - the desired precision for the desired quadrature rule
!    nfuns - the number of input functions
!    funuser - a user-supplied external subroutine which conforms to the 
!      gaussfun2d interface and evaluates the input functions and their
!      derivatives  
!    userptr - a "void *" pointer which is passed along to funuser
!    (nquad, xs, ys, whts) - 
!
!  Output parameters: 
!    (nquad, xs, ys, whts) - the resulting quadrature rule
!
type(gq2dwrapper_data), pointer :: wrapdata
type(c_ptr)                     :: wrapptr
double precision, allocatable   :: xs0(:), ys0(:), whts0(:)

! double precision, allocatable :: coefs(:,:)
! double precision, allocatable :: xsnew(:), whtsnew(:)
! double precision, allocatable :: rints(:), rints0(:)
! double precision, allocatable :: signifs(:),amatr(:,:)
! integer, allocatable          :: idxs(:), idxs2(:)


eps0 = epsilon(0.0d0)
call bilegepw_quad(disc,nquad0,xs0,ys0,whts0)

allocate(wrapdata)
wrapptr          = c_loc(wrapdata)
wrapdata%disc   => disc
allocate(wrapdata%coefs(nquad0,nfuns))
call bilegepw_coefs(disc,vals,wrapdata%coefs)
call gaussquad2d(eps,nfuns,gaussquad2d_wrapper,gaussquad2d_funregion,wrapptr,nquad,xs,ys,whts)

end subroutine


subroutine gaussquad2d_wrapper(nfuns,n,xs,ys,vals,dersx,dersy,userptr)
implicit double precision (a-h,o-z)
type(c_ptr)                     :: userptr
double precision                :: xs(:),ys(:)
double precision                :: vals(:,:), dersx(:,:), dersy(:,:)
type(gq2dwrapper_data), pointer :: wrapdata

call c_f_pointer(userptr,wrapdata)

do i=1,n
x   = xs(i)
y   = ys(i)
call bilegepw_evalder(wrapdata%disc,wrapdata%coefs,x,y,vals(i,:),dersx(i,:),dersy(i,:))
end do

end subroutine


subroutine gaussquad2d_funregion(x,y,wht,isvalid,userptr)
implicit double precision (a-h,o-z)
type(c_ptr)                                          :: userptr
procedure(gaussfun2d)                                :: funuser
procedure(gaussfun2d_region)                         :: funregion
type(gq2dwrapper_data), pointer                      :: wrapdata
call c_f_pointer(userptr,wrapdata)


isvalid = 0
ntopboxes = wrapdata%disc%ntopboxes

do i=1,ntopboxes
x1 = wrapdata%disc%boxes(i)%x1
y1 = wrapdata%disc%boxes(i)%y1
x2 = wrapdata%disc%boxes(i)%x2
y2 = wrapdata%disc%boxes(i)%y2

ifcontains = 1
if (x .lt. x1) ifcontains = 0
if (x .gt. x2) ifcontains = 0
if (y .lt. y1) ifcontains = 0
if (y .gt. y2) ifcontains = 0

if (ifcontains .eq. 1) then
isvalid = 1
return
endif

end do


end subroutine


subroutine gaussquad2d_newton(ier,eps,nfuns,funuser,funregion,userptr,rints,nquad,xs,ys,whts,errl2)
implicit double precision (a-h,o-z)
procedure(gaussfun2d)                                :: funuser
procedure(gaussfun2d_region)                         :: funregion
type(c_ptr)                                          :: userptr
double precision                                     :: rints(:)
double precision, allocatable, intent(inout)         :: xs(:), ys(:),whts(:)
!
!  Apply Newton iterations to the nonlinear system of equations (2) in the
!  hopes of improving the accuracy of a quadrature rule for a collection of
!  input functions.
!
!  Input parameters:
!    eps - precision for the computations; more explicitly, the procedure is considered
!      to have converged when this precision is reached
!    disc - a data structure describing the piecewise Legendre scheme used to represent
!     the input functions
!    nfuns - the number of input functions
!    funuser - a user-supplied external subroutine which conforms to the 
!      gaussfun2d interface and evaluates the input functions and their
!      derivatives  
!    funregion - an external subroutine which indicates which points are within
!      the integration domain
!    userptr - a "void *" pointer which is passed along to funuser and funregion
!    rints - the integrals of the input functions
!    (nquad,xs,whts) -  the quadrature rule to regine
!
!  Output parameters:
!     ier = 0   indicates that the specified quadrature rule was successfully refined
!     ier = 4   means that the maximum number of iteration was exceeded during the
!               step-length control procedure AND the step would have moved the
!               quadrature outside of the interval
!     ier = 8   means the maximum number of iterations was exceeded during the
!               step-length control procedure
!     ier = 16  means the maximum number of Newton iterations was exceeded without
!               finding a sufficiently accuracy quadrature rule
!    
!     (nquad,xs,whts)  - assuming ier = 0, this will be the newly refined
!       quadrature rule
!
!     errl2 - the error in the newly found quadrature rule
!

double precision, allocatable :: amatr(:,:), rints0(:)
double precision, allocatable :: xx(:), yy(:)
double precision, allocatable :: xs0(:), ys0(:), whts0(:)

ier      = 0 

eps0     = epsilon(0.0d0)
maxiters = 30
nextra   = 2

allocate(amatr(nfuns,3*nquad), xx(3*nquad), yy(nfuns) )
allocate(xs0(nquad), ys0(nquad), whts0(nquad), rints0(nfuns))

!
!  Perform Newton iterations
!
do iter=1,maxiters

call elapsed(t1)
call gaussquad2d_newtstep(ier,nfuns,funuser,funregion,userptr,rints,nquad,xs,ys,whts,errl2,nsteps)
call elapsed(t2)

write (*, "(A,I2,A,I2,A,D10.2,A,F10.2)") "iter = ",iter,",  nsteps = ", nsteps, &
  ",  errl2 = ",errl2,",  time = ",t2-t1

if (errl2 .lt. eps) exit
if (ier .ne. 0) return

end do

!
!  Exit with an error if we have exceeded the maximum number of Newton
!  iterations
!
if (iter .gt. maxiters) then
ier = 16
return
endif

!
!  Take a few more Newton steps to refine the accuracy of the rule
!

do i=1,nextra
call gaussquad2d_newtstep(ier,nfuns,funuser,funregion,userptr,rints,nquad,xs,ys,whts,errl2,nsteps)
write (*, "(A,I2,A,I2,A,D10.2,A,F10.2)") "iter = ",iter+i,",  nsteps = ", nsteps, &
  ",  errl2 = ",errl2,",  time = ",t2-t1

end do


end subroutine



subroutine gaussquad2d_newtstep(ier,nfuns,funuser,funregion,userptr,rints,nquad,xs,ys,whts,errl2,nsteps)
implicit double precision (a-h,o-z)
procedure(gaussfun2d)                                :: funuser
procedure(gaussfun2d_region)                         :: funregion
type(c_ptr)                                          :: userptr
double precision                                     :: rints(:)
double precision, allocatable, intent(inout)         :: xs(:), ys(:), whts(:)
!
!  Perform a single Newton iteration.
!
!  Input parameters:
!    eps - precision for the computations; more explicitly, the procedure is considered
!      to have converged when this precision is reached
!    funuser - a user-supplied external subroutine which conforms to the 
!      gaussfun2d interface and evaluates the input functions and their
!      derivatives  
!    funregion - an external subroutine which indicates which points are within
!      the integration domain
!    userptr - a "void *" pointer which is passed along to funuser and funregion
!    rints - the integrals of the input functions
!    (nquad,xs,whts) -  the quadrature rule to regine
!
!  Output parameters:
!     ier = 0   indicates that the error in the quadrature rule was improved
!     ier = 4   means that an updated quadrature whose nodes lie inside the
!               interval could not be found
!     ier = 8   means the maximum number of steps in the step-length control
!               procedure was exceeded before any improvement was ssen

!    
!     (nquad,xs,ys,whts)  - assuming ier = 0, this will be the newly refined
!       quadrature rule
!
!     errl2 - the error in the newly found quadrature rule
!     nsteps - the number of steps taken by the step-length control procedure
!

double precision, allocatable :: amatr(:,:), rints0(:)
double precision, allocatable :: xx(:), yy(:)
double precision, allocatable :: xs0(:), ys0(:), whts0(:)

ier      = 0 
eps0     = epsilon(0.0d0)

maxsteps = 40
minsteps = 1
dscale   = 2d0

allocate(amatr(nfuns,3*nquad), xx(3*nquad), yy(nfuns) )
allocate(xs0(nquad), ys0(nquad), whts0(nquad), rints0(nfuns))

! form the coefficient matrix for the linear system
call gaussquad2d_linearized(nfuns,funuser,userptr,nquad,xs,ys,whts,amatr,rints0)
yy     = rints-rints0
errl2  = norm2(yy)

! solve the linearized system in a least squares sense
call leastsq(nfuns,3*nquad,amatr,yy,xx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! perform step-length control 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

errbest   = errl2
alphabest = 0.0d0
alpha     = dscale

do istep=1,maxsteps
nsteps = istep
alpha = alpha/dscale
xs0   = xs   + alpha*xx(1:nquad)
ys0   = ys   + alpha*xx(nquad+1:2*nquad)
whts0 = whts + alpha*xx(2*nquad+1:3*nquad)

! REJECT STEPS WHICH ARE OUTSIDE THE REGION
ifaccept = 1
do i=1,nquad
call funregion(xs0(i),ys0(i),whts0(i),ifaccept,userptr)
if (ifaccept .eq. 0) exit
end do
if (ifaccept .eq. 0) cycle

call gaussquad2d_rints(nfuns,funuser,userptr,nquad,xs0,ys0,whts0,rints0)
yy     = rints-rints0
errnew = norm2(yy)

if (errnew .lt. errbest) then
errbest   = errnew
alphabest = alpha
endif

if (errbest .lt. errl2 .AND. istep .ge. minsteps) exit

end do

!
!  Update the quadrature rule IF we have improved the accuracy, otherwise
!  we are stuck and we should exit with an error
!
if (errbest .lt. errl2) then
xs    = xs   + alphabest*xx(1:nquad)
ys    = ys   + alphabest*xx(nquad+1:2*nquad)
whts  = whts + alphabest*xx(2*nquad+1:3*nquad)
errl2 = errbest
return
endif

ier = 8

end subroutine




subroutine gaussquad2d_signifs(nfuns,funuser,userptr,rints,nquad,xs,ys,whts,signifs)
implicit double precision (a-h,o-z)
procedure(gaussfun2d)                                :: funuser
type(c_ptr)                                          :: userptr
double precision                                     :: rints(:)
double precision                                     :: xs(:), ys(:), whts(:)
double precision                                     :: signifs(:)
!
!  Compute the "significance" of each node of a quadrature rule.  The significance
!  of a node is the norm of the solution of the linerized problem when
!
!  The naive approach to computing the significance has O(n^4) complexity.  This
!  routine uses an approached based on the normal equations and the Sherman-Morrison-
!  Woodbury formula to perform the calculations in O(N^3) times.
!
!  Input parameters:
!    rints - an array specifying the values of the integrals of the input functions
!    (nquad,xs,whts) - the input quadrature rule
!
!  Output parameters:
!    signifs - an array of length nquad giving the significance of each node
!

double precision, allocatable :: a(:,:), b(:,:)
double precision, allocatable :: a0(:,:), b0(:,:)
double precision, allocatable :: xx(:), yy(:), zz(:), rints0(:)
!double precision, allocatable :: xsnew(:), whtsnew(:), rints0(:)


allocate( a(nfuns,3*nquad),  b(nfuns, nfuns)     )
allocate( a0(nfuns,3*nquad), b0(nfuns,nfuns)     )
allocate( rints0(nfuns), yy(nfuns), xx(3*nquad), zz(nfuns)  )


!
!  Compute the Jacobian A of the system, form the matrix AA^t, and invert 
!  AA^t.
!



call gaussquad2d_linearized(nfuns,funuser,userptr,nquad,xs,ys,whts,a,rints0)

b = matmul(a,transpose(a))
call invert(nfuns,b)

do ipt=1,nquad

a0 = a
b0 = b

call smw_update(nfuns,b0,-a0(:,ipt),a0(:,ipt))
call smw_update(nfuns,b0,-a0(:,nquad+ipt),a0(:,nquad+ipt))
call smw_update(nfuns,b0,-a0(:,2*nquad+ipt),a0(:,2*nquad+ipt))

a0(:,ipt)          = 0
a0(:,nquad+ipt)    = 0
a0(:,2*nquad+ipt)  = 0

do j=1,nquad
a0(:,j+2*nquad)  = a0(:,j+2*nquad)*whts(j)
end do

rints0 = sum(a0(:,2*nquad+1:3*nquad),2)

do j=1,nquad
a0(:,j+2*nquad) = a0(:,j+2*nquad)/whts(j)
end do


yy = rints-rints0
zz = matmul(b0,yy)
xx = matmul(transpose(a0),zz)
signifs(ipt) = norm2(xx)

end do


end subroutine


subroutine gaussquad2d_linearized(nfuns,funuser,userptr,nquad,xs,ys,whts,amatr,rints)
implicit double precision (a-h,o-z)
procedure(gaussfun2d)                                :: funuser
type(c_ptr)                                          :: userptr
double precision                                     :: xs(:), ys(:), whts(:)
double precision                                     :: amatr(:,:), rints(:)
!
!  Form the coefficient matrix for the linearization of the system (2) with
!
!    x_1,...,x_m, y1_, ..., w_m, w_1,...,w_m
!
!  a user-supplied quadrature rule.  Also, use the quadrature rule to approximate
!  the integrals of the input functions.
!
!  Input parameters:
!    (nquad,xs,whts) -  the quadrature rule to regine
!
!  Output parameters:
!    amatr - the (nfuns, 3*nquad) coefficient matricx
!    rints - an array giving the approximations of the integrals of the input
!      functions
!


double precision, allocatable :: vals(:,:), dersx(:,:), dersy(:,:)

!
!  Form the coefficient matrix for the linearized system
!


allocate(vals(nquad,nfuns), dersx(nquad,nfuns), dersy(nquad,nfuns) )
call funuser(nfuns,nquad,xs,ys,vals,dersx,dersy,userptr)

!
!  Scale the derivatives by the weights 
!

do j=1,nquad
x   = xs(j)
y   = ys(j)
wht = whts(j)
amatr(:,j)         = dersx(j,:)*wht
amatr(:,j+nquad)   = dersy(j,:)*wht
amatr(:,j+2*nquad) = vals(j,:)
end do


do i=1,nquad
vals(i,:) = vals(i,:) * whts(i)
end do
rints = sum(vals,1)

end subroutine



subroutine gaussquad2d_rints(nfuns,funuser,userptr,nquad,xs,ys,whts,rints)
implicit double precision (a-h,o-z)
procedure(gaussfun2d)                                :: funuser
type(c_ptr)                                          :: userptr
double precision                                     :: xs(:), ys(:), whts(:)
double precision                                     :: rints(:)
!
!  Use a quadrature to approximate the integrals of the input functions.
!
!  Input parameters:
!    nfuns - the number of input functions
!    (nquad,xs,whts) -  the quadrature rule
!
!  Output parameters:
!    rints - an array giving the approximations of the integrals of the input
!      functions
!
double precision, allocatable :: vals(:,:), dersx(:,:), dersy(:,:)


allocate(vals(nquad,nfuns), dersx(nquad,nfuns), dersy(nquad,nfuns) )
call funuser(nfuns,nquad,xs,ys,vals,dersx,dersy,userptr)
do i=1,nquad
vals(i,:) = vals(i,:) * whts(i)
end do

rints = sum(vals,1)

end subroutine

end module


