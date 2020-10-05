!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for constructing and manipulating piecewise Legendre 
!  expansions.  By a piecewise Legendre expansion of order n on the subintervals
!
!   (a_1, b_1),  (a_2,b_2), ... (a_m,b_m),
!
!  we mean a sum of the form
!
!             m                         n-1                       
!    f(x) = \sum \chi     (x)          \sum   c     \tilde{P}     (x),                      (2)
!            j=1      [a_j, b_j)        i=0    i,j           i,j
!
!
!  where \chi_[a_j,b_j) is the characteristic function of the half-open interval 
!  [a_j, b_j) and \tilde{P}_i,j(x) denotes the L^2 normalized Legendre polynomial of 
!  degree i on the interval (a_j,b_j).  That is, 
!
!
!                           (  2n+1 )           (    2        a_j+b_j )
!   \tilde{P}_i,j(x)  = sqrt( ----- )       P_i ( ------- x + ------- ) ,
!                           (  b-a  )           ( b_j-a_j     a_j-b_j )
!
!
!  with P_i(x) the Legendre polynomial of degree i.  This module provides routines 
!  for handling both real-valued and complex-values expansions of this type.
!
!  Expansions are represented via either the vector
!
!     ( a_0,0         )
!     ( a_1,0         )
!     ( ...           )
!     ( a_{n-1},0     )
!     ( a_0,1         )
!     ( a_1,1         )
!     ( ...           )                                                                       (3)
!     ( a_{n-1},1     ) 
!     ( ...           )
!     ( a_{0},{m-1}   ) 
!     ( ...           )
!     ( a_{n-1},{m-1} )
!
!  of expansion coefficients or via the vector
!
!     (  f(x_1)  \sqrt{w_1}  )
!     (  f(x_n)  \sqrt{w_2}  )
!     (         ...          )                                                               (4)
!     (  f(x_mn) \sqrt{w_mn} )
!
!  of the *scaled* values of the expansion at the nodes of the mn-point quadrature
!  rule constructing by amalgamating the n-point Gauss-Legendre rules on each
!  of the intervals (a_j,b_j].  We refer to this mechanism for representing
!  functions as a ``piecewise Legendre discretization scheme'' or simply a
!  ``discretization scheme'' and the amalgamated quadrature as a piecewise Legendre 
!  quadrature or the discretization quadrature.
!
!  The following subroutines should be regarded as publicly callable:
!
!    legepw_init - initialize the data structure describing a piecewise
!      Legendre discretization scheme; the user must supply an  initial
!      interval or list of intervals for the scheme
!
!    legepw_uniform - refine a discretization scheme until the length of every
!      subinterval is at most equal to a specified quantity
!  
!    legepw_adap - refine a discretization scheme until it suffices to represent a
!      user-specified collection of functions to a specified precision
!
!    legepw_quad - return the piecewise Legendre quadrature rule associated with
!      a discretization scheme
!
!    legepw_coefs - given the *scaled* values of one or more piecewise Legendre
!      expansions at the nodes of the discretization quadrature, compute their 
!      piecewise Legendre coefficient expansions
!
!    legepw_diff - given the *scaled* values of one or more piecewise Legendre
!      expansions at the nodes of the discretization quadrature, compute the
!      *scaled* values of their first derivatives at the nodes of the
!      discretization quadrature 
!
!      **WARNING**  THE CONDITION NUMBER OF SPECTRAL DIFFERENTIATION 
!                   DETERIORATES WITH INCREASING ORDER
!
!    legepw_eval - given the vector of coefficients for one or more expansions
!      of the form (2), compute their values at a specified point
!
!    legepw_evalder - given the vector of coefficients for one or more expansions
!      of the form (2), compute their values and their first derivative
!      at a specified point
!
!    legepw_interp - given the *scaled* values of one or more piecewise Legendre
!      expansions at the nodes of the discretization quadrature, evaluate them at a
!      specified point 
!
!    legepw_order - increase the order of the Legendre expansions used to represent
!      functions on each subinterval; this routine is principally useful for
!      representing products of the form
!
!              f_i(x) g_j(x)    i = 1,....,n,  j = 1,...,m
!
!      by discretizing the collections { f_i } and { g_j } separately
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module legepw

use utils
use legendre
use iso_c_binding

type       legepw_disc

integer                       :: nlege           ! the length of the Gauss-Legendre quadrature
double precision, allocatable :: xslege(:)       ! the Gauss-Legendre quadrature nodes
double precision, allocatable :: whtslege(:)     ! the Gauss-Legendre quadrature weights
double precision, allocatable :: u(:,:)          ! the weighted values to coefs matrix
double precision, allocatable :: adiff(:,:)      ! the spectral differentiation matrix


! second quadrature for adaptive integration
integer                       :: kextra
integer                       :: nlege2           ! the length of the Gauss-Legendre quadrature
double precision, allocatable :: xslege2(:)       ! the Gauss-Legendre quadrature nodes
double precision, allocatable :: whtslege2(:)     ! the Gauss-Legendre quadrature weights
double precision, allocatable :: u2(:,:)          ! the weighted values to coefs matrix


integer                       :: maxints         ! the maximum number of subintervals
integer                       :: nints           ! the number of subintervals
double precision, allocatable :: ab(:,:)         ! the (2,nints) list of subintervals

end type   legepw_disc

interface

subroutine legepw_adapfun(nfuns,n,xs,whts,vals,userptr)
import c_ptr
double precision           :: xs(n)
double precision           :: whts(n)
double precision           :: vals(:,:)
type(c_ptr)                :: userptr
end subroutine

subroutine legepw_cadapfun(nfuns,n,xs,whts,vals,userptr)
import c_ptr
double precision           :: xs(n)
double precision           :: whts(n)
double complex             :: vals(:,:)
type(c_ptr)                :: userptr
end subroutine

end interface


interface     legepw_init
module procedure legepw_init1
module procedure legepw_init2
end interface legepw_init

interface     legepw_coefs
module procedure legepw_coefs1
module procedure legepw_coefs2
module procedure legepw_coefs3
module procedure legepw_coefs4
end interface legepw_coefs

interface     legepw_diff
module procedure legepw_diff1
module procedure legepw_diff2
module procedure legepw_diff3
module procedure legepw_diff4
end interface legepw_diff

interface     legepw_eval
module procedure legepw_eval1
module procedure legepw_eval2
module procedure legepw_eval3
module procedure legepw_eval4
end interface legepw_eval

interface     legepw_evalder
module procedure legepw_evalder1
module procedure legepw_evalder2
module procedure legepw_evalder3
module procedure legepw_evalder4
end interface legepw_evalder

interface     legepw_interp
module procedure legepw_interp1
module procedure legepw_interp2
module procedure legepw_interp3
module procedure legepw_interp4
end interface legepw_interp

contains


subroutine legepw_init1(disc,nlege,a,b)
implicit double precision (a-h,o-z)
type(legepw_disc)         :: disc
!
!  Initalize the structure which stores a piecewise Legendre discretization
!  scheme.  The initial scheme will consist of the single interval (a,b).
!
!  Input parameters:
!   nlege - the number of Legendre nodes to use on each subinterval
!   (a,b) - the initial interval in the discretization scheme
!
!  Output parameters:
!   disc - the data structure which specifies a discretization scheme on
!     the interval (a,b) will be initialized
!

maxints = 10000
kextra  = 8

if ( allocated(disc%ab) )        deallocate(disc%ab)

if ( allocated(disc%xslege) )    deallocate(disc%xslege)
if ( allocated(disc%whtslege) )  deallocate(disc%whtslege)
if ( allocated(disc%u) )         deallocate(disc%u)
if ( allocated(disc%adiff) )     deallocate(disc%adiff)

if ( allocated(disc%xslege2) )   deallocate(disc%xslege2)
if ( allocated(disc%whtslege2) ) deallocate(disc%whtslege2)
if ( allocated(disc%u2) )        deallocate(disc%u2)

call legendre_quad(nlege,disc%xslege,disc%whtslege)
call legendre_coefsmatrix(nlege,disc%xslege,disc%whtslege,disc%u)
call legendre_diffmatrix(nlege,disc%xslege,disc%whtslege,disc%adiff)

nlege2 = nlege+kextra
call legendre_quad(nlege2,disc%xslege2,disc%whtslege2)
call legendre_coefsmatrix(nlege2,disc%xslege2,disc%whtslege2,disc%u2)


disc%nints    = 1
disc%nlege    = nlege
disc%maxints  = maxints
disc%kextra   = kextra
disc%nlege2   = nlege2

allocate( disc%ab(2,1) )
disc%ab(1,1)  = a
disc%ab(2,1)  = b

end subroutine


subroutine legepw_init2(disc,nlege,nints,ab)
implicit double precision (a-h,o-z)
type(legepw_disc)         :: disc
double precision          :: ab(:,:)
!
!  Initialize the data structure describing a piecewise Legendre expansion.
!  The initial scheme will consist of a collection of intervals specified
!  by the user.
!
!  Input parameters:
!   nlege - the number of Legendre nodes to use on each subinterval
!   nints - the number of intervals in the initial 
!   ab - a (2,nints) array whose (1,j) and (2,j) entries give the left-
!     hand and right-hand endpoints of the jth interval
!
!     NOTE: THE INTERVALS MUST BE ORDERED FROM LEFT TO RIGHT
!
!  Output parameters:
!   disc - the data structure which specifies a discretization scheme on
!     the interval (a,b) will be initialized
!

maxints = 10000
kextra  = 8

if ( allocated(disc%ab) )       deallocate(disc%ab)

if ( allocated(disc%xslege) )   deallocate(disc%xslege)
if ( allocated(disc%whtslege) ) deallocate(disc%whtslege)
if ( allocated(disc%u) )        deallocate(disc%u)
if ( allocated(disc%adiff) )    deallocate(disc%adiff)

if ( allocated(disc%xslege2) )   deallocate(disc%xslege2)
if ( allocated(disc%whtslege2) ) deallocate(disc%whtslege2)
if ( allocated(disc%u2) )        deallocate(disc%u2)

call legendre_quad(nlege,disc%xslege,disc%whtslege)
call legendre_coefsmatrix(nlege,disc%xslege,disc%whtslege,disc%u)
call legendre_diffmatrix(nlege,disc%xslege,disc%whtslege,disc%adiff)

nlege2 = nlege+kextra
call legendre_quad(nlege2,disc%xslege2,disc%whtslege2)
call legendre_coefsmatrix(nlege2,disc%xslege2,disc%whtslege2,disc%u2)

disc%nints    = nints
disc%nlege    = nlege
disc%maxints  = maxints
disc%kextra   = kextra
disc%nlege2   = nlege2

allocate( disc%ab(2,nints) )
disc%ab       = ab(:,1:nints)

end subroutine




subroutine legepw_uniform(disc,dlen)
implicit double precision (a-h,o-z)
type(legepw_disc)                  :: disc
!
!  Refine a discretization scheme until each subinterval is of length at most
!  dlen.
!
!  Input parameters:
!    disc - the data structure describing an existing discretization scheme
!    dlen - the maximum length of an interval in the scheme.
!
!  Ouput parameters:
!    disc - the discretization described by this data structure is refined
!      as specified
!

double precision, allocatable      :: stack(:,:), ab0(:,:)

nstack   = 0
nints0   = 0
maxints  = disc%maxints
maxstack = maxints
nlege    = disc%nlege
nints    = disc%nints

allocate(stack(2,maxstack), ab0(2,maxints) )

nints0           = 0 
nstack           = nints
stack(:,1:nints) = disc%ab

do while (nstack .gt. 0) 

a      = stack(1,nstack)
b      = stack(2,nstack)
nstack = nstack-1

ifsplit = 0
if (b-a .gt. dlen) ifsplit = 1

if (ifsplit .eq. 1) then
if (nstack .ge. maxstack) then
call prina("in legepw_uniform, stack overflowed")
stop
endif

nstack               = nstack+1
stack(1,nstack)      = a
stack(2,nstack)      = (a+b)/2

nstack               = nstack+1
stack(1,nstack)      = (a+b)/2
stack(2,nstack)      = b

else

if (nints0 .eq. maxints) then
call prina("in legepw_uniform, maximum number of intervals exceeded")
stop
endif

nints0               = nints0+1
ab0(1,nints0)        = a
ab0(2,nints0)        = b

endif
end do

call quicksort(nints0*2,ab0)
deallocate(disc%ab)
allocate(disc%ab(2,nints0))
disc%nints = nints0
disc%ab    = ab0(:,1:nints0)

end subroutine


subroutine legepw_gradedmesh(disc)
implicit double precision (a-h,o-z)
type(legepw_disc)                  :: disc
!
!  Refine a discretization scheme so that every interval is at a distance
!  from 0 which is less than or equal to its length.
!
!  Input parameters:
!    disc - the data structure describing an existing discretization scheme
!
!  Ouput parameters:
!    disc - the discretization described by this data structure is refined
!      as specified
!

double precision, allocatable      :: stack(:,:), ab0(:,:)

nstack   = 0
nints0   = 0
maxints  = disc%maxints
maxstack = maxints
nlege    = disc%nlege
nints    = disc%nints

allocate(stack(2,maxstack), ab0(2,maxints) )

nints0           = 0 
nstack           = nints
stack(:,1:nints) = disc%ab

do while (nstack .gt. 0) 

a      = stack(1,nstack)
b      = stack(2,nstack)
nstack = nstack-1

dlen   = b-a

ifsplit = 0
if (dlen .gt. abs(a)) ifsplit=1

if (ifsplit .eq. 1) then
if (nstack .ge. maxstack) then
call prina("in legepw_uniform, stack overflowed")
stop
endif

nstack               = nstack+1
stack(1,nstack)      = a
stack(2,nstack)      = (a+b)/2

nstack               = nstack+1
stack(1,nstack)      = (a+b)/2
stack(2,nstack)      = b

else

if (nints0 .eq. maxints) then
call prina("in legepw_uniform, maximum number of intervals exceeded")
stop
endif

nints0               = nints0+1
ab0(1,nints0)        = a
ab0(2,nints0)        = b

endif

end do

call quicksort(nints0*2,ab0)
deallocate(disc%ab)
allocate(disc%ab(2,nints0))
disc%nints = nints0
disc%ab    = ab0(:,1:nints0)

end subroutine


subroutine legepw_adap(eps,nfuns,fun,disc,userptr)
implicit double precision (a-h,o-z)
type(legepw_disc)                  :: disc
procedure(legepw_adapfun)          :: fun
type(c_ptr)                        :: userptr
!
!  Adaptively refine a discretization scheme until a user-specified collection
!  of square integrable functions can be represented to a specified precision
!  using it.
!
!  Input parameters:
!    eps - the desired precision for the discretization
!    nfuns - the number of input functions
!    fun - an external subroutine conforming to the legepw_adapfun interface
!      which returns the matrix giving the values of the collection of input 
!      functions at a specified collection of input points
!    userptr - a "void *" pointer which is passed to the user-supplied function
!    disc - the data structure describing an existing discretization scheme
!
!  Input/Output parameters:
!    disc - the discretization described by this data structure is refined
!      until is suffices to represent the input functions
!

double precision, allocatable      :: vals(:,:), sums1(:), sums2(:)
double precision, allocatable      :: stack(:,:), ab0(:,:), coefs0(:)
double precision, allocatable      :: coefs(:,:), ratios(:)
double precision, allocatable      :: xs0(:), whts0(:)

nstack   = 0
nints0   = 0
maxints  = disc%maxints
maxstack = maxints
nlege    = disc%nlege
nints    = disc%nints
kk       = 8
epssq    = eps**2

!allocate(stack(2,maxstack), vals(nlege2,nfuns), coefs0(kextra) )
! allocate(xs0(nlege2), whts0(nlege2), ab0(2,maxints), sums1(nfuns), sums2(nfuns) )

allocate(stack(2,maxstack), vals(nlege,nfuns), coefs0(nlege), coefs(nlege,nfuns) )
allocate(xs0(nlege), whts0(nlege), ab0(2,maxints), sums1(nfuns), sums2(nfuns), ratios(nfuns) )


nstack           = nints
do i=1,nstack
stack(:,nints-i+1) = disc%ab(:,i)
end do

do while(nstack .gt. 0)
a       = stack(1,nstack)
b       = stack(2,nstack)
nstack  = nstack-1

xs0     = (b-a)/2 * disc%xslege + (a+b)/2
whts0   = (b-a)/2 * disc%whtslege

call fun(nfuns,nlege,xs0,whts0,vals,userptr)
coefs = matmul(disc%u, vals)
coefs = coefs**2

sums1   = sum(coefs(1:nlege-kk,:),2)
sums2   = sum(coefs(nlege-kk+1:nlege,:),2)
sums1   = sums1+sums2+1.0d0
ratios  = sums2 / sums1


ifsplit = 0
if (maxval(ratios) .gt. epssq) ifsplit = 1


if (ifsplit .eq. 1) then

if (nstack .ge. maxstack) then
call prina("in legepw_adap, stack overflowed")
stop
endif

nstack               = nstack+1
stack(1,nstack)      = (a+b)/2
stack(2,nstack)      = b

nstack               = nstack+1
stack(1,nstack)      = a
stack(2,nstack)      = (a+b)/2


else

if (nints0 .eq. maxints) then
call prina("in legepw_adap, maximum number of intervals exceeded")
stop
endif

nints0               = nints0+1
ab0(1,nints0)        = a
ab0(2,nints0)        = b

endif

end do



! call quicksort(nints0*2,ab0)
deallocate(disc%ab)
allocate(disc%ab(2,nints0))
disc%nints = nints0
disc%ab    = ab0(:,1:nints0)


end subroutine


subroutine legepw_cadap(eps,nfuns,fun,disc,userptr)
implicit double precision (a-h,o-z)
type(legepw_disc)                  :: disc
procedure(legepw_cadapfun)         :: fun
type(c_ptr)                        :: userptr
!
!  Adaptively refine a discretization scheme until a user-specified collection
!  of square integrable functions can be represented to a specified precision
!  using it.
!
!  Input parameters:
!    eps - the desired precision for the discretization
!    nfuns - the number of input functions
!    fun - an external subroutine conforming to the legepw_adapfun interface
!      which returns the matrix giving the values of the collection of input 
!      functions at a specified collection of input points
!    userptr - a "void *" pointer which is passed to the user-supplied function
!    disc - the data structure describing an existing discretization scheme
!
!  Input/Output parameters:
!    disc - the discretization described by this data structure is refined
!      until is suffices to represent the input functions
!

double complex, allocatable        :: vals(:,:), coefs(:,:), sums1(:), sums2(:)
double precision, allocatable      :: stack(:,:), ab0(:,:)
double precision, allocatable      :: xs0(:), whts0(:)

nstack   = 0
nints0   = 0
maxints  = disc%maxints
maxstack = maxints
nlege    = disc%nlege
nints    = disc%nints

allocate(stack(2,maxstack), vals(nlege,nfuns), coefs(nlege,nfuns) )
allocate(xs0(nlege), whts0(nlege), ab0(2,maxints), sums1(nfuns), sums2(nfuns) )

nstack           = nints
stack(:,1:nints) = disc%ab

do while(nstack .gt. 0)
a       = stack(1,nstack)
b       = stack(2,nstack)
nstack  = nstack-1
xs0     = (b-a)/2 * disc%xslege + (a+b)/2
whts0   = (b-a)/2 * disc%whtslege

call fun(nfuns,nlege,xs0,whts0,vals,userptr)

ifsplit = 0

coefs = matmul(disc%u,vals)
coefs = coefs**2
sums1 = sum(coefs(1:nlege/2,:),1)
sums2 = sum(coefs(nlege/2+1:nlege,:),1)
sums1 = sqrt(sums2/(sums1+sums2))

dd1 = maxval(abs(sums1))
if (dd1 .gt. eps) ifsplit = 1

if (ifsplit .eq. 1) then
if (nstack .ge. maxstack) then
call prina("in legepw_cadap, stack overflowed")
stop
endif

nstack               = nstack+1
stack(1,nstack)      = a
stack(2,nstack)      = (a+b)/2

nstack               = nstack+1
stack(1,nstack)      = (a+b)/2
stack(2,nstack)      = b

else

if (nints0 .eq. maxints) then
call prina("in legepw_cadap, maximum number of intervals exceeded")
stop
endif

nints0               = nints0+1
ab0(1,nints0)        = a
ab0(2,nints0)        = b

endif


end do

call quicksort(nints0*2,ab0)
deallocate(disc%ab)
allocate(disc%ab(2,nints0))
disc%nints = nints0
disc%ab    = ab0(:,1:nints0)

end subroutine



subroutine legepw_quad(disc,nquad,xs,whts)
implicit double precision (a-h,o-z)
type(legepw_disc)                              :: disc
double precision, allocatable, intent(out)     :: xs(:), whts(:)
!
!  Return the discretization quadrature rule.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!
!  Output parameters:
!    nquad - the number of nodes in the quadrature rule
!    xs - the nodes in the quadruature rule
!    whts - the weights in the quadrature rule
!

nints = disc%nints
nlege = disc%nlege
nquad = nints*nlege

allocate(xs(nquad), whts(nquad) )

i2 = 0
do int = 1,nints
a  = disc%ab(1,int)
b  = disc%ab(2,int)
i1 = i2+1
i2 = i1+nlege-1

xs(i1:i2)   = disc%xslege*(b-a)/2 + (b+a)/2
whts(i1:i2) = disc%whtslege*(b-a)/2

end do

end subroutine


subroutine legepw_coefs1(disc,vals,coefs)
implicit double precision (a-h,o-z)
type(legepw_disc)              :: disc
double precision               :: vals(:), coefs(:)
!
!  Given the vector (4) of *scaled* values of a real-valued expansion of the form
!  (2) at the nodes of the discretization quadrature, compute the vector (3) of 
!  expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    vals - the array specifying the scaled values of the expansion
!
!  Output parameters:
!    coefs - the array specifying the vecotr of coefficients
!
!

nlege = disc%nlege
nints = disc%nints

i2 = 0
do int=1,nints
a  = disc%ab(1,int)
b  = disc%ab(2,int)
i1 = i2+1
i2 = i1+nlege-1
coefs(i1:i2)  = matmul(disc%u,vals(i1:i2))
end do

end subroutine

subroutine legepw_coefs2(disc,vals,coefs)
implicit double precision (a-h,o-z)
type(legepw_disc)              :: disc
double complex                 :: vals(:), coefs(:)
!
!  Given the vector (4) of *scaled* values of a complex-valued expansion of the form
!  (2) at the nodes of the discretization quadrature, compute the vector (3) of 
!  expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    vals - the array specifying the scaled values of the expansion
!
!  Output parameters:
!    coefs - the array specifying the vecotr of coefficients
!

nlege = disc%nlege
nints = disc%nints

i2 = 0
do int=1,nints
a  = disc%ab(1,int)
b  = disc%ab(2,int)
i1 = i2+1
i2 = i1+nlege-1
coefs(i1:i2)  = matmul(disc%u,vals(i1:i2))
end do

end subroutine


subroutine legepw_coefs3(disc,vals,coefs)
implicit double precision (a-h,o-z)
type(legepw_disc)              :: disc
double complex                 :: vals(:,:), coefs(:,:)
!
!  Given the vector (4) of *scaled* values of a collection of complex-valued expansions
!  of the form (2) at the nodes of the discretization quadrature, their compute the vectors 
!  (3) of  expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    vals - the matrix whose jth column gives the scaled values of the jth
!      input expansion
!
!  Output parameters:
!    coefs - the matrix whose jth column gives the coefficinets of the jth
!      input expansion
!

nlege = disc%nlege
nints = disc%nints

i2 = 0
do int=1,nints
a  = disc%ab(1,int)
b  = disc%ab(2,int)
i1 = i2+1
i2 = i1+nlege-1
coefs(i1:i2,:)  = matmul(disc%u,vals(i1:i2,:))
end do

end subroutine


subroutine legepw_coefs4(disc,vals,coefs)
implicit double precision (a-h,o-z)
type(legepw_disc)              :: disc
double precision               :: vals(:,:), coefs(:,:)
!
!  Given the vector (4) of *scaled* values of a collection of real-valued expansions
!  of the form (2) at the nodes of the discretization quadrature, their compute the vectors 
!  (3) of  expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    vals - the matrix whose jth column gives the scaled values of the jth
!      input expansion
!
!  Output parameters:
!    coefs - the matrix whose jth column gives the coefficinets of the jth
!      input expansion
!

nlege = disc%nlege
nints = disc%nints

i2 = 0
do int=1,nints
a  = disc%ab(1,int)
b  = disc%ab(2,int)
i1 = i2+1
i2 = i1+nlege-1
coefs(i1:i2,:)  = matmul(disc%u,vals(i1:i2,:))
end do

end subroutine


subroutine legepw_diff1(disc,vals,ders)
implicit double precision (a-h,o-z)
type(legepw_disc)              :: disc
double precision               :: vals(:), ders(:)
!
!  Given the vector (4) of *scaled* values of a real-valued expansion of the form
!  (2) at the nodes of the discretization quadrature, compute the vector (3) of the scaled
!  values of its derivative at the same.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    vals - the array specifying the scaled values of the expansion
!
!  Output parameters:
!    ders - the array specifying the scaled values of the expansion's derviative
!
!

nlege = disc%nlege
nints = disc%nints


i2 = 0
do int=1,nints
a  = disc%ab(1,int)
b  = disc%ab(2,int)
i1 = i2+1
i2 = i1+nlege-1
ders(i1:i2)  = matmul(disc%adiff,vals(i1:i2))
ders(i1:i2)  = ders(i1:i2)*2/(b-a)
end do

end subroutine

subroutine legepw_diff2(disc,vals,ders)
implicit double precision (a-h,o-z)
type(legepw_disc)              :: disc
double complex                 :: vals(:), ders(:)
!
!  Given the vector (4) of *scaled* values of a real-valued expansion of the form
!  (2) at the nodes of the discretization quadrature, compute the vector (3) of the scaled
!  values of its derivative at the same.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    vals - the array specifying the scaled values of the expansion
!
!  Output parameters:
!    ders - the array specifying the scaled values of the expansion's derviative
!
!

nlege = disc%nlege
nints = disc%nints

i2 = 0
do int=1,nints
a  = disc%ab(1,int)
b  = disc%ab(2,int)
i1 = i2+1
i2 = i1+nlege-1
ders(i1:i2)  = matmul(disc%adiff,vals(i1:i2))
ders(i1:i2)  = ders(i1:i2)*2/(b-a)

end do

end subroutine


subroutine legepw_diff3(disc,vals,ders)
implicit double precision (a-h,o-z)
type(legepw_disc)              :: disc
double precision               :: vals(:,:), ders(:,:)
!
!  Given the vector (4) of *scaled* values of a collection of real-valued expansions of the form
!  (2) at the nodes of the discretization quadrature, compute the (3) scaled values
!  of their derivative at the same.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    vals - the array specifying the scaled values of the expansion
!
!  Output parameters:
!    ders - the array specifying the scaled values of the expansion's derviative
!
!

nlege = disc%nlege
nints = disc%nints

i2 = 0
do int=1,nints
a  = disc%ab(1,int)
b  = disc%ab(2,int)
i1 = i2+1
i2 = i1+nlege-1
ders(i1:i2,:)  = matmul(disc%adiff,vals(i1:i2,:))
ders(i1:i2,:)  = ders(i1:i2,:)*2/(b-a)

end do

end subroutine

subroutine legepw_diff4(disc,vals,ders)
implicit double precision (a-h,o-z)
type(legepw_disc)              :: disc
double complex                 :: vals(:,:), ders(:,:)
!
!  Given the vector (4) of *scaled* values of a collection of complex-valued expansions of the form
!  (2) at the nodes of the discretization quadrature, compute the (3) scaled values
!  of their derivative at the same.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    vals - the array specifying the scaled values of the expansion
!
!  Output parameters:
!    ders - the array specifying the scaled values of the expansion's derviative
!
!

nlege = disc%nlege
nints = disc%nints

i2 = 0
do int=1,nints
a  = disc%ab(1,int)
b  = disc%ab(2,int)
i1 = i2+1
i2 = i1+nlege-1
ders(i1:i2,:)  = matmul(disc%adiff,vals(i1:i2,:))
ders(i1:i2,:)  = ders(i1:i2,:)*2/(b-a)
end do

end subroutine


subroutine legepw_eval1(disc,coefs,x,valout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double precision             :: coefs(:)
double precision             :: valout
!
!  Evaluate a real-valued expansion of the form (2) at a specified point
!  given its vector (3) of expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    coefs - the vector (3) of expansion coefficients
!    x - the point at which to evaluate the input expansion
!
!  Output parameters:
!    valout - the value of the expansion
!

nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)

i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_eval(nlege,coefs(i1:i2),xx,valout)
valout = valout * sqrt(2/(b-a))

end subroutine


subroutine legepw_eval2(disc,coefs,x,valout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double complex               :: coefs(:)
double complex               :: valout
!
!  Evaluate a complex-valued expansion of the form (2) at a specified point
!  given its vector (3) of expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    coefs - the vector (3) of expansion coefficients
!    x - the point at which to evaluate the input expansion
!
!  Output parameters:
!    valout - the value of the expansion
!

nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!

do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)

i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_eval(nlege,coefs(i1:i2),xx,valout)
valout = valout * sqrt(2/(b-a))

end subroutine


subroutine legepw_eval3(disc,coefs,x,valsout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double complex               :: coefs(:,:)
double complex               :: valsout(:)
!
!  Evaluate a collection of complex-valued expansions of the form (2) 
!  at a specified point given the vectors (3) of their expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    coefs - a matrix whose jth column contains the expansion coefficients
!      for the jth input expansion
!    x - the point at which to evaluate the input expansion
!
!  Output parameters:
!    valsout - the jth entry of this array will contain the value of the
!      jth expansion at x
!

nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)

i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_eval(nlege,coefs(i1:i2,:),xx,valsout)
valsout = valsout * sqrt(2/(b-a))

end subroutine


subroutine legepw_eval4(disc,coefs,x,valsout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double precision             :: coefs(:,:)
double precision             :: valsout(:)
!
!  Evaluate a collection of real-valued expansions of the form (2) 
!  at a specified point given the vectors (3) of their expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    coefs - a matrix whose jth column contains the expansion coefficients
!      for the jth input expansion
!    x - the point at which to evaluate the input expansion
!
!  Output parameters:
!    valsout - the jth entry of this array will contain the value of the
!      jth expansion at x
!

nlege = disc%nlege
nints = disc%nints


!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)

i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_eval(nlege,coefs(i1:i2,:),xx,valsout)
valsout = valsout * sqrt(2/(b-a))

end subroutine


subroutine legepw_evalder1(disc,coefs,x,valout,derout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double precision             :: coefs(:)
double precision             :: valout
double precision             :: derout
!
!  Evaluate a real-valued expansion of the form (2) at a specified point
!  given its vector (3) of expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    coefs - the vector (3) of expansion coefficients
!    x - the point at which to evaluate the input expansion
!
!  Output parameters:
!    valout - the value of the expansion at the point x
!    derout - the derivative of the expansion at the point x
!

nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)

i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_evalder(nlege,coefs(i1:i2),xx,valout,derout)
valout = valout * sqrt(2/(b-a))
derout = derout * 2/(b-a)*sqrt(2.0d0/(b-a))


end subroutine


subroutine legepw_evalder2(disc,coefs,x,valout,derout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double complex               :: coefs(:)
double complex               :: valout
double complex               :: derout
!
!  Evaluate a complex-valued expansion of the form (2) at a specified point
!  given its vector (3) of expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    coefs - the vector (3) of expansion coefficients
!    x - the point at which to evaluate the input expansion
!
!  Output parameters:
!    valout - the value of the expansion
!

nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)

i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_evalder(nlege,coefs(i1:i2),xx,valout,derout)
valout = valout * sqrt(2/(b-a))
derout = derout * 2/(b-a)*sqrt(2.0d0/(b-a))

end subroutine


subroutine legepw_evalder3(disc,coefs,x,valsout,dersout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double complex               :: coefs(:,:)
double complex               :: valsout(:)
double complex               :: dersout(:)

!
!  Evaluate a collection of complex-valued expansions of the form (2) 
!  at a specified point given the vectors (3) of their expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    coefs - a matrix whose jth column contains the expansion coefficients
!      for the jth input expansion
!    x - the point at which to evaluate the input expansion
!
!  Output parameters:
!    valsout - the jth entry of this array will contain the value of the
!      jth expansion at x
!

nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)

i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_eval(nlege,coefs(i1:i2,:),xx,valsout)
valsout = valsout * sqrt(2/(b-a))
dersout = dersout * 2/(b-a)*sqrt(2.0d0/(b-a))

end subroutine


subroutine legepw_evalder4(disc,coefs,x,valsout,dersout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double precision             :: coefs(:,:)
double precision             :: valsout(:)
double precision             :: dersout(:)
!
!  Evaluate a collection of real-valued expansions of the form (2) 
!  at a specified point given the vectors (3) of their expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the discretization scheme
!    coefs - a matrix whose jth column contains the expansion coefficients
!      for the jth input expansion
!    x - the point at which to evaluate the input expansion
!
!  Output parameters:
!    valsout - the jth entry of this array will contain the value of the
!      jth expansion at x
!


nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .le. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)

i1 = 1 + (int-1)*nlege
i2 = nlege*int
call legendre_evalder(nlege,coefs(i1:i2,:),xx,valsout,dersout)
valsout = valsout * sqrt(2/(b-a))
dersout = dersout * 2/(b-a)*sqrt(2.0d0/(b-a))

end subroutine

subroutine legepw_interp1(disc,vals,x,valout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double precision             :: vals(:)
double precision             :: valout
!
!
!  Input parameters:
!
!  Output parameters:
!

nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)
i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_interp(nlege,disc%xslege,disc%whtslege,vals(i1:i2),xx,valout)
valout = valout * sqrt(2/(b-a))

end subroutine


subroutine legepw_interp2(disc,vals,x,valout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double complex               :: vals(:)
double complex               :: valout
!
!
!  Input parameters:
!
!  Output parameters:
!

nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)

i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_interp(nlege,disc%xslege,disc%whtslege,vals(i1:i2),xx,valout)
valout = valout * sqrt(2/(b-a))
end subroutine


subroutine legepw_interp3(disc,vals,x,valsout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double complex               :: vals(:,:)
double complex               :: valsout(:)
!
!
!  Input parameters:
!
!  Output parameters:
!

nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)


i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_interp(nlege,disc%xslege,disc%whtslege,vals(i1:i2,:),xx,valsout)
valsout = valsout * sqrt(2/(b-a))
end subroutine


subroutine legepw_interp4(disc,vals,x,valsout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double precision             :: vals(:,:)
double precision             :: valsout(:)
!
!
!  Input parameters:
!
!  Output parameters:
!

nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)


i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_interp(nlege,disc%xslege,disc%whtslege,vals(i1:i2,:),xx,valsout)
valsout = valsout * sqrt(2/(b-a))
end subroutine


subroutine legepw_evalder(disc,coefs,x,valout,derout)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
double precision             :: coefs(:)
double precision             :: valout, derout
!
!
!  Input parameters:
!
!  Output parameters:
!

nlege = disc%nlege
nints = disc%nints

!
!  Find the interval containing the point
!
do int=1,nints-1
if (x .lt. disc%ab(2,int)) exit
end do

a  = disc%ab(1,int)
b  = disc%ab(2,int)
xx = (2*x - (b+a) ) /(b-a)

i1 = 1 + (int-1)*nlege
i2 = nlege*int

call legendre_evalder(nlege,coefs(i1:i2),xx,valout,derout)
valout = valout * sqrt(2.0d0/(b-a))
derout = derout * 2/(b-a)*sqrt(2.0d0/(b-a))

end subroutine


subroutine legepw_order(disc,nlege)
implicit double precision (a-h,o-z)
type(legepw_disc)            :: disc
!
!  Change the order of the Legendre expansions used to represent functions
!  on each interval.
!
!  Input parameters:
!    disc - the data structure describing an existing discretization scheme
!    nlege - the number of terms in the expansion which will be used to 
!     represent functions on each discretization subinterval 
!
!  Output parameters:
!    disc - the data structure describing the discretization scheme is
!     updated
!

if ( allocated(disc%u) )        deallocate(disc%u)
if ( allocated(disc%adiff) )    deallocate(disc%adiff)
if ( allocated(disc%xslege) )   deallocate(disc%xslege)
if ( allocated(disc%whtslege) ) deallocate(disc%whtslege)

disc%nlege = nlege
allocate(disc%xslege(nlege), disc%whtslege(nlege) )

call legendre_quad(nlege, disc%xslege, disc%whtslege)
call legendre_coefsmatrix(nlege,disc%xslege,disc%whtslege,disc%u)
call legendre_diffmatrix(nlege,disc%xslege,disc%whtslege,disc%adiff)

end subroutine



end module
