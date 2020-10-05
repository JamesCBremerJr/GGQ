!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for constructing "generalized Chebyshev" quadrature rules.
!  A generalized Chebyshev quadrature rule for a linearly independent collection of n
!  functions f_1,...,f_n is an quadrature rule of the form
!
!        b                    n
!    \int f(x) dx  \approx  \sum  f(x_j) w_j                                                 (1)
!        a                   j=1
!
!  which is exact for  the functions f_1, ..., f_n (or, at least, the formula (1) achieves
!  near machine precision accuracy when one of the f_j is substituted for f). 
!
!  The following subroutines should be regarded as publicly callable:
!
!    chebquad - construct a "generalized Chebyshev quadrature" for a user-supplied
!      collection of n functions represented by piecewise Legendre expansions ala
!      legepw.f90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module chebyshev_quad

use legendre
use linalg
use legepw

contains

subroutine chebquad(disc,nfuns,vals,nquad,xs,whts)
implicit double precision (a-h,o-z)
type(legepw_disc)                          :: disc
double precision                           :: vals(:,:)
double precision, allocatable, intent(out) :: xs(:), whts(:)

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
!    vals - an array the jth column of which specifies the *scaled* Values of the jth
!     input function at the nodes of the piecewise Gauss-Legendre quadrature
!     described by the disc structure
!
!  Output parameters: 
!    nquad - the number of nodes in the resulting quadrature, which will be 
!      the numerical dimension of the space spanned by the input functions
!    xs - an array specifying the nodes of the resulting quadrature rule
!    whts - an array specifying the weights of the resulting quadrature rule
!
!

double precision, allocatable :: xs0(:),whts0(:), rints(:)
double precision, allocatable :: a(:,:), b(:,:), r(:,:)
integer, allocatable          :: ipivs(:), idxs(:)

eps0  = epsilon(0.0d0)

!
!  Construct the oversampled quadrature rule from the user-supplied
!  discretization
!
call legepw_quad(disc,nquad0,xs0,whts0)

!
!  Form the transpose of the input matrix and scale it by the square roots
!  of the oversampled quadrature
!

allocate(a(nfuns,nquad0), rints(nfuns) )
a = transpose(vals)

do j=1,nquad0
a(:,j) = a(:,j)*sqrt(whts0(j))
end do

rints = sum(a,2)

do j=1,nquad0
a(:,j) = a(:,j)/sqrt(whts0(j))
end do

!
!  Subsample the quadrature nodes
!
call intdecomp(eps0,nfuns,nquad0,a,krank,ipivs,r)

!
!  Subsample the matrix and solve for the quadrature weights
!

nquad = krank
allocate(b(nfuns,krank), xs(krank) ,whts(krank) )
b = a(:,ipivs)

call leastsq(nfuns,krank,b,rints,whts)
xs   = xs0(ipivs)
whts = whts * sqrt(whts0(ipivs))
!whts = whts * whts0(ipivs)

!
!  Sort the quadrature nodes
!
allocate(idxs(nquad))
call insort2(nquad,xs,idxs)
whts = whts(idxs)


end subroutine



end module
