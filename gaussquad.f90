!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for constructing "generalized Gaussian" quadrature rules.
!  A generalized Gaussian quadrature rule for a linearly independent collection of 2n
!  functions f_1,...,f_2n is an quadrature rule of the form
!
!        b                    n
!    \int f(x) dx  \approx  \sum  f(x_j) w_j                                                 (1)
!        a                   j=1
!
!  which is exact for the functions f_1, ..., f_2n (or, at least, the formula (1) 
!  achieves specified accuracy when one of the f_j is substituted for f). 
!
!  This code takes as input an existing quadrature rule of length > n which integrates
!  the functions f_1, ..., f_2n.  It attempts to downsample this rule by removing
!  quadrature nodes one at a time, and solving the obvious system of n nonlinear 
!  equations in 2m variales
!
!      m                       b  
!    \sum  f_i(x_j) w_j  = \int f_i(x) dx,   i=1,...,2n                                      (2)
!     j=1                      a
!
!       
!  to (hopefully) restore the accuracy of the reduced quadrature rule.
!
!  The following subroutines should be regarded as publicly callable:
!
!    gaussquad - by applying Newton's method the obvious nonlinear system of equations
!      satisfied by the rule (1), downsample an existing quadrature point-by-point
!      in the hopes of creating a generalized Gaussian quadrature rule for a collection
!      of user-supplied functions which are represented via piecewise Legendre 
!      expansions ala legepw.f90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gaussian_quad

use utils
use legendre
use linalg
use legepw

  ! set this to 1 to suppress output from these routines
integer, private, parameter              ::  ifsuppress = 0

contains


subroutine gaussquad(eps,iorder,disc,nfuns,vals,nquad,xs,whts)
implicit double precision (a-h,o-z)
type(legepw_disc)                                    :: disc
double precision                                     :: vals(:,:)
double precision, allocatable, intent(inout)         :: xs(:), whts(:)
!
!  Starting with an existing quadrature rule for a collection of n user-supplied
!  functions, attempt to reduce the formula point-by-point using Newton iterations
!  in the hopes of constructing a generalized Gaussian quadrature rule for the
!  collection.  The input functions are represented via piecewise Legendre
!  expansions.
!
!  Input parameters:
!    eps - the desired precision for the desired quadrature rule
!    iorder - an integer parameter specifying the order in which it attempts to
!     remove quadrature nodes
!
!     iorder   = 0    means no reordering of nodes
!     iorder   = 1    means reorder the nodes based on "significance"
!     iorder   = 2    means order the nodes according to the magnitude of their weights
!     iorder   = 3    means use a random ordering of the nodes
!     iorder   = 4    means reorder the nodes based on "signifcance" with nodes
!                     associated with negative weights given priority
!
!    disc - a data structure describing the piecewise Legendre scheme used to represent
!     the input functions
!    nfuns - the number of input functions
!    vals - an array the jth column of which specifies the *scaled* Values of the jth
!     input function at the nodes of the piecewise Gauss-Legendre quadrature
!     described by the disc structure
!
!  Output parameters: 
!    nquad - the number of nodes in the resulting quadrature
!    xs - an array specifying the nodes of the resulting quadrature rule
!    whts - an array specifying the weights of the resulting quadrature rule
!

double precision, allocatable :: coefs(:,:)
double precision, allocatable :: xs0(:), whts0(:)
double precision, allocatable :: xsnew(:), whtsnew(:)
double precision, allocatable :: rints(:), rints0(:)
double precision, allocatable :: signifs(:),amatr(:,:)
integer, allocatable          :: idxs(:), idxs2(:)

double precision, allocatable :: coefsders(:,:)

eps0 = epsilon(0.0d0)
call elapsed(time1)

if (ifsuppress .ne. 1) then
write (*,*)
endif

!
!  Fetch the discretization quadrature and construct the matrix of 
!  expansion coefficients
!
call legepw_quad(disc,nquad0,xs0,whts0)
allocate(coefs(nquad0,nfuns))
call legepw_coefs(disc,vals,coefs)


!
!  Compute the integrals of the input functions
!

allocate(rints(nfuns),rints0(nfuns))
do i=1,nquad0
vals(i,:) = vals(i,:) * sqrt(whts0(i))
end do
rints = sum(vals,1)
do i=1,nquad0
vals(i,:) = vals(i,:) / sqrt(whts0(i))
end do


!
!  Attempt to remove nquad/2 points
!

do nremove=1,nquad/2
call elapsed(time2)
if (ifsuppress .ne. 1) then
write (*,"(A,I3,A,I3,A)") "[--- gaussquad: nfuns = ",nfuns,", nquad = ",nquad," ---------------------------]"
endif


!
!  Reorder the nodes of the existing quadrature rule as specified by the user
!  and attempt to remove each point in that order
!

if (iorder .gt. 0) then

allocate(signifs(nquad), idxs(nquad), idxs2(nquad) )
allocate(xsnew(nquad), whtsnew(nquad))

if (iorder .eq. 1 .OR. iorder .eq. 4) then
call gaussquad_signifs(disc,nfuns,coefs,rints,nquad,xs,whts,signifs)

elseif (iorder .eq. 2) then
signifs = whts
elseif (iorder .eq. 3) then
do j=1,nquad
idxs(j) = j
end do

call permute_randomly(nquad,idxs,idxs2)

do j=1,nquad
signifs(j) = idxs2(j)
end do

endif


call insort2(nquad,signifs,idxs)
xsnew   = xs(idxs)
whtsnew = whts(idxs)
xs      = xsnew
whts    = whtsnew


call gaussquad_rints(disc,nfuns,coefs,nquad,xs,whts,rints0)
errl2 = norm2(rints-rints0)
deallocate(signifs, idxs, xsnew, whtsnew, idxs2)

endif


nquadnew = nquad-1
allocate(xsnew(nquadnew), whtsnew(nquadnew))

do ipt=1,nquad

if (ifsuppress .ne. 1) then
write (*,"(A,I3,A)") "attempting to remove ipt = ",ipt
endif

call elapsed(t1)

!
!  Remove ipt from the quadrature rule, and solve the obvious nonlinear system of equations via 
!  Newton's method to attempt to improve the accuracy of the resulting rule
!


idx = 0  
do i=1,nquad
if (i .eq. ipt) then
x0 = xs(i)
wht0 = whts(i)
cycle
endif

idx           = idx+1
xsnew(idx)    = xs(i)
whtsnew(idx)  = whts(i)
end do
call gaussquad_newton(ier,eps,disc,nfuns,coefs,rints,nquadnew,xsnew,whtsnew,errl2)

if (ier .eq. 0) exit
end do


call elapsed(t2)

! if we have removed a point, update the quadrature 
if (ier .eq. 0 ) then

deallocate(xs,whts)
allocate(xs(nquadnew), whts(nquadnew))
nquad = nquadnew
xs    = xsnew
whts  = whtsnew
deallocate(xsnew,whtsnew)


if (ifsuppress .ne. 1) then
write (*,"(A,F9.2,A)") "successful, time = ",t2-t1, " seconds "
write (*,*) ""
endif

else
if (ifsuppress .ne. 1) then
write (*,"(A,F9.2,A)") "failed,     time = ",t2-t1, " seconds "
write (*,*) ""
endif

! otherwise, we have the best formula we could find already
deallocate(xsnew,whtsnew)
exit
endif

end do

 
allocate(xsnew(nquad),idxs(nquad))

call insort2(nquad,xs,idxs)
whtsnew = whts(idxs)
whts    = whtsnew



end subroutine



subroutine gaussquad_newton(ier,eps,disc,nfuns,coefs,rints,nquad,xs,whts,errl2)
implicit double precision (a-h,o-z)
type(legepw_disc)                                    :: disc
double precision                                     :: coefs(:,:)
double precision                                     :: rints(:)
double precision, allocatable, intent(inout)         :: xs(:), whts(:)
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
!    coefs - an array the jth column of which specifies the piecewise Gauss-Legendre
!     expansion coefficients of the jth input function 
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
double precision, allocatable :: xs0(:), whts0(:)

ier      = 0 

eps0     = epsilon(0.0d0)
maxiters = 20
nextra   = 2

nints    = disc%nints
a        = disc%ab(1,1)
b        = disc%ab(2,nints)

allocate(amatr(nfuns,2*nquad), xx(2*nquad), yy(nfuns) )
allocate(xs0(nquad), whts0(nquad), rints0(nfuns))

!
!  Perform Newton iterations
!
do iter=1,maxiters

call elapsed(t1)
call gaussquad_newtstep(ier,disc,nfuns,coefs,rints,nquad,xs,whts,errl2,nsteps)
if (ier .ne. 0) return
call elapsed(t2)


if (ifsuppress .ne. 1) then
write (*, "(A,I2,A,I2,A,D10.2,A,F10.2)") "iter = ",iter,",  nsteps = ", nsteps, &
  ",  errl2 = ",errl2,",  time = ",t2-t1
endif

if (errl2 .lt. eps) exit

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

call gaussquad_newtstep(ier,disc,nfuns,coefs,rints,nquad,xs,whts,errl2,nsteps)
if (ifsuppress .ne. 1) then
write (*, "(A,I2,A,I2,A,D10.2,A,F10.2)") "iter = ",iter+i,",  nsteps = ", nsteps, &
  ",  errl2 = ",errl2,",  time = ",t2-t1
endif

end do

end subroutine



subroutine gaussquad_newtstep(ier,disc,nfuns,coefs,rints,nquad,xs,whts,errl2,nsteps)
implicit double precision (a-h,o-z)
type(legepw_disc)                                    :: disc
double precision                                     :: coefs(:,:)
double precision                                     :: rints(:)
double precision, allocatable, intent(inout)         :: xs(:), whts(:)
!
!  Perform a single Newton iteration.
!
!  Input parameters:
!    eps - precision for the computations; more explicitly, the procedure is considered
!      to have converged when this precision is reached
!    disc - a data structure describing the piecewise Legendre scheme used to represent
!     the input functions
!    nfuns - the number of input functions
!    coefs - an array the jth column of which specifies the piecewise Gauss-Legendre
!     expansion coefficients of the jth input function 
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
!     (nquad,xs,whts)  - assuming ier = 0, this will be the newly refined
!       quadrature rule
!
!     errl2 - the error in the newly found quadrature rule
!     nsteps - the number of steps taken by the step-length control procedure
!

double precision, allocatable :: amatr(:,:), rints0(:)
double precision, allocatable :: xx(:), yy(:)
double precision, allocatable :: xs0(:), whts0(:)

ier      = 0 
eps0     = epsilon(0.0d0)

maxsteps = 20
minsteps = 1
dscale   = 4

nints    = disc%nints
a        = disc%ab(1,1)
b        = disc%ab(2,nints)

allocate(amatr(nfuns,2*nquad), xx(2*nquad), yy(nfuns) )
allocate(xs0(nquad), whts0(nquad), rints0(nfuns))

! form the coefficient matrix for the linear system
call gaussquad_linearized(disc,nfuns,coefs,nquad,xs,whts,amatr,rints0)

yy     = rints-rints0
errl2  = norm2(yy)

! solve the linearized system in a least squares sense
call leastsq(nfuns,2*nquad,amatr,yy,xx)

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
whts0 = whts + alpha*xx(nquad+1:2*nquad)


! REJECT STEPS WHICH ARE OUTSIDE THE INTERVAL
ifaccept = 1
do i=1,nquad
if (xs0(i) .gt. b .OR. xs0(i) .lt. a) ifaccept=0
end do
if (ifaccept .eq. 0) cycle

! MAKE SURE WE DON't GO OUT OF BOUNDS AS A LAST RESORT
! xs0 = max(xs0,a)
! xs0 = min(xs0,b)

call gaussquad_rints(disc,nfuns,coefs,nquad,xs0,whts0,rints0)
yy     = rints-rints0
errnew = norm2(yy)

if (errnew .lt. errbest) then
errbest   = errnew
alphabest = alpha
endif

if (errbest .lt. errl2 .AND. istep .ge. minsteps) exit

end do

if (ifaccept .eq. 0) then
ier = 4
return
endif



!
!  Update the quadrature rule IF we have improved the accuracy, otherwise
!  we are stuck and we should exit with an error
!
if (errbest .lt. errl2) then
xs    = xs   + alphabest*xx(1:nquad)
whts  = whts + alphabest*xx(nquad+1:2*nquad)
xs    = max(xs,a)
xs    = min(xs,b)
errl2 = errbest
return
endif

ier = 8

end subroutine




subroutine gaussquad_signifs(disc,nfuns,coefs,rints,nquad,xs,whts,signifs)
implicit double precision (a-h,o-z)
type(legepw_disc)                                    :: disc
double precision                                     :: coefs(:,:), rints(:)
double precision                                     :: xs(:), whts(:)
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
!    disc - a data structure describing the piecewise Legendre scheme used to represent
!     the input functions
!    nfuns - the number of input functions
!    coefs - an array the jth column of which specifies the piecewise Gauss-Legendre
!     expansion coefficients of the jth input function 
!    rints - an array specifying the values of the integrals of the input functions
!    (nquad,xs,whts) - the input quadrature rule
!
!  Output parameters:
!    signifs - an array of length nquad giving the significance of each node
!

double precision, allocatable :: a(:,:), b(:,:)
double precision, allocatable :: a0(:,:), b0(:,:)
double precision, allocatable :: xx(:), yy(:), zz(:)
double precision, allocatable :: xsnew(:), whtsnew(:), rints0(:)

!
!  Use the O(N^3) algorithm which might have poor constants
!  

call elapsed(t1)

allocate( a(nfuns,2*nquad),  b(nfuns, nfuns)     )
allocate( a0(nfuns,2*nquad), b0(nfuns,nfuns)     )
allocate( rints0(nfuns), yy(nfuns), xx(2*nquad), zz(nfuns)  )
!
!  Compute the Jacobian A of the system, form the matrix AA^t, and invert 
!  AA^t.
!

call gaussquad_linearized(disc,nfuns,coefs,nquad,xs,whts,a,rints0)
b = matmul(a,transpose(a))
call invert(nfuns,b)

do ipt=1,nquad


a0 = a
b0 = b

call smw_update(nfuns,b0,-a0(:,ipt),a0(:,ipt))
call smw_update(nfuns,b0,-a0(:,nquad+ipt),a0(:,nquad+ipt))
a0(:,ipt)        = 0
a0(:,nquad+ipt)  = 0


do j=1,nquad
a0(:,j+nquad)  = a0(:,j+nquad)*whts(j)
end do
rints0 = sum(a0(:,nquad+1:2*nquad),2)
do j=1,nquad
a0(:,j+nquad) = a0(:,j+nquad)/whts(j)
end do

yy = rints-rints0
zz = matmul(b0,yy)
xx = matmul(transpose(a0),zz)

signifs(ipt) = norm2(xx)

end do


end subroutine


subroutine gaussquad_linearized(disc,nfuns,coefs,nquad,xs,whts,amatr,rints0)
implicit double precision (a-h,o-z)
type(legepw_disc)                                    :: disc
double precision                                     :: coefs(:,:)
double precision                                     :: xs(:), whts(:)
double precision                                     :: amatr(:,:), rints0(:)
!
!  Form the coefficient matrix for the linearization of the system (2) with
!
!    x_1,...,x_m, w_1,...,w_m
!
!  a user-supplied quadrature rule.  Also, use the quadrature rule to approximate
!  the integrals of the input functions.
!
!  Input parameters:
!    disc - a data structure describing the piecewise Legendre scheme used to represent
!     the input functions
!    nfuns - the number of input functions
!    coefs - an array the jth column of which specifies the piecewise Gauss-Legendre
!     expansion coefficients of the jth input function 
!    (nquad,xs,whts) -  the quadrature rule to regine
!
!  Output parameters:
!    amatr - the (nfuns, 2*nquad) coefficient matricx
!    rints0 - an array giving the approximations of the integrals of the input
!      functions
!


!
!  Form the coefficient matrix for the linearized system
!

do j=1,nquad
x   = xs(j)
wht = whts(j)
call legepw_evalder(disc,coefs,x,amatr(:,j+nquad),amatr(:,j))

amatr(:,j)        = amatr(:,j)*wht
amatr(:,j+nquad)  = amatr(:,j+nquad)*wht
end do

rints0 = sum(amatr(:,nquad+1:2*nquad),2)

do j=1,nquad
amatr(:,j+nquad) = amatr(:,j+nquad)/whts(j)
end do

end subroutine



subroutine gaussquad_rints(disc,nfuns,coefs,nquad,xs,whts,rints0)
implicit double precision (a-h,o-z)
type(legepw_disc)                                    :: disc
double precision                                     :: coefs(:,:)
double precision                                     :: xs(:), whts(:)
double precision                                     :: rints0(:)
!
!  Use a quadrature to approximate the integrals of the input functions/
!
!  Input parameters:
!    disc - a data structure describing the piecewise Legendre scheme used to represent
!     the input functions
!    nfuns - the number of input functions
!    coefs - an array the jth column of which specifies the piecewise Gauss-Legendre
!     expansion coefficients of the jth input function 
!    (nquad,xs,whts) -  the quadrature rule to regine
!
!  Output parameters:
!    amatr - the (nfuns, 2*nquad) coefficient matricx
!    rints0 - an array giving the approximations of the integrals of the input
!      functions
!
double precision, allocatable :: b(:,:)

!
!  Form the coefficient matrix for the linearized system
!

allocate(b(nquad,nfuns))

do i=1,nquad
x   = xs(i)
wht = whts(i)
call legepw_eval(disc,coefs,x,b(i,:))
b(i,:) = b(i,:)*wht
end do

rints0 = sum(b,1)

end subroutine

end module


