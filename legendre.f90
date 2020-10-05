!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for constructing Guass-Legendre quadrature rules, and for
!  forming and manipulating univariate Legendre expansions.  A univariate Legendre 
!  expansion is a sum of the form 
!
!              n-1
!      f(x) = \sum   a_i \tilde{P}_i(x),                                                   (1)
!              i=0
!
!  where \tilde{P}_i denotes the L^2 normalized Legendre polynomial of the first
!  kind of degree i.  This code provides routines for both real-valued and complex-valued
!  expansions.
!
!  Expansions of the form (1) are represented either via the vector of coefficients
!   
!    ( a_0     )
!    ( a_1     )
!    ( ...     )                                                                            (2)
!    ( a_{n-1} )
!
!  or via the vector 
!
!    ( f(x_1) \sqrt{w_1} )
!    ( f(x_2) \sqrt{w_2} )
!    (      ...          )                                                                 (3)
!    ( f(x_n) \sqrt{w_n} )
!
!  of their *scaled* values at the nodes of the n-point Gauss-Legendre quadrature 
!  rule.  Note that the matrix taking (2) to (3) is orthogonal.
!
!  The following routines should be regarded as public:
!    
!    legendre_quad - return the nodes and weights of the n-point Gauss-Legendre 
!      quadrature rule; note that this routine uses an O(n^2) algorithm and is
!      intended for constructing quadrature rules of relatively small sizes (n < 200)
!
!    legendre_coefs - compute the vector (2) of the coefficients of one of more 
!       expansions of the form (1) given the vectors of their scaled values (3)
!
!    legendre_interp - use barycentric interpolation to evaluate one of more expansions 
!     of the form (1) given their *scaled* values at nodes of the n-point
!     Gauss-Legendre quadrature rule
!
!    legendre_eval - evaluate one or more expansions of the form (1) at a specified
!      point given their coefficient vectors (2)
!
!    legendre_evalder - evaluate an expansion of the form (1) and its
!      derivative at a specified point given the vector (2) of coefficients
!
!    legendre_interpmatrix - return an (m,n) matrix which takes the vector (3) of the
!      *scaled* values of an expansion at the nodes of the n-point Gauss-Legendre rule
!      to *scaled* values at a user-specified collection of points
!
!    legendre_coefsmatrix - return the (n,n) matrix which takes the vector (3) of 
!      the *scaled* values of an expansion of the form (1) to the vector (2) of its
!      coefficients -- note that this matrix is orthogonal and its transpose
!      takes the vector (2) to the vector (3)
!
!    legendre_diffmatrix - return the (n,n) ``spectral differentiation matrix'' which
!      takes the vector (3) of scaled values of an expansion of the form (1) to
!      the vector of *scaled* values of the derivative of the expansion
!
!      **WARNING** SPECTRAL DIFFERENTIATION BECOME INCREASINGLY ILL-CONDITIONED AS THE
!                  ORDER OF THE EXPANSIONS GROW
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module legendre

interface     legendre_interp
module procedure legendre_interp1
module procedure legendre_interp2
module procedure legendre_interp3
module procedure legendre_interp4
end interface legendre_interp

interface     legendre_coefs
module procedure legendre_coefs1
module procedure legendre_coefs2
module procedure legendre_coefs3
module procedure legendre_coefs4
end interface legendre_coefs

interface     legendre_eval
module procedure legendre_eval1
module procedure legendre_eval2
module procedure legendre_eval3
module procedure legendre_eval4
end interface legendre_eval

interface     legendre_evalder
module procedure legendre_evalder1
module procedure legendre_evalder2
module procedure legendre_evalder3
module procedure legendre_evalder4
end interface legendre_evalder

contains

subroutine legendre_quad(n,xslege,whtslege)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out)            :: xslege(:),whtslege(:)
!
!  Form the n-point Gauss-Legendre quadrature rule.  Note that this routine uses a
!  O(n^2) algorithm and is inteded for use in the case in which n is relatively
!  small (n < 100 or so).
!
!  Input parameters:
!    n - the length of the desired quadrature rule
!
!  Output parameters:
!    xslege - the nodes of the desired quadrature rule
!    whtlege - the weights of the desired quadrature rule
!

double precision :: pols(n+1)

data pi / 3.14159265358979323846264338327950288d0 /
allocate( xslege(n), whtslege(n) )

maxiters = 12
eps0     = epsilon(0.0d0)

if (n .gt. 100) then
call legequad2(n,xslege,whtslege)
return
endif

!
!   Use Newton's method and the recurrence relation to find half of the
!   roots --- the others are obtained via symmetry.
!
!   Note that we also store the value of the derivative at each of the obtained
!   roots for use in computing the weights below.
!
ifodd = 0
nn = (n+1)/2
if (nn /= n/2) then
ifodd=1
nn=nn-1
endif

!
!  Use Tricomi's formula to  estimate the roots of P_{n+1}
!

do i =nn+1,n
   dk = i-nn-1
   dn = n
   theta = (4*(ceiling(dn/2)-dk)-1) / (4*dn+2) * pi
   x0 = 1.0d0 - (dn-1)/(8*dn**3)-1.0d0/(384.0d0*dn**4)* &
        (39.0d0-28.0d0/sin(theta)**2)
   xslege(i)=x0*cos(theta)
enddo

!
!  Perform Newton iterations in order to refine the estimates.
!
do iter = 1,maxiters

!
!  Evaluate the Legendre polynomial of degree n at each point; save
!  the values of the derivative in the whts array.
!
do i=nn+1,n
call legeder(n,xslege(i),pols(i),whtslege(i))
end do

!
!  Perform one Newton iteration
!
pols(nn+1:n)     = pols(nn+1:n)/whtslege(nn+1:n)
xslege(nn+1:n)   = xslege(nn+1:n) - pols(nn+1:n)

if(norm2(pols(nn+1:n)) < eps0) then
exit
endif

end do

if (iter == maxiters)  then
print *,"legemdre_quad: newton iterations failed!"
stop
end if

!
! Compute the weights using the derivatives we stored above.
!
do j=nn+1,n
x           = xslege(j)
dd          = 2.0d0/(1.0d0-x**2)
whtslege(j) = dd/(whtslege(j)**2) * (n+0.5d0)
end do

!
! Reflect the quadrature nodes.
!
do j=1,nn
xslege(j)   = -xslege(n-j+1)
whtslege(j) = whtslege(n-j+1)
end do

!
! Handle the root at 0 if n is odd.
!

if (ifodd .eq. 1) then
x0          = 0
call legeder(n,x0,pol,der)
xslege(nn+1)   = x0
whtslege(nn+1) = 2.0d0/(der**2) * (n+0.5d0)
endif

end subroutine





subroutine legendre_interp1(n,xslege,whtslege,vals,x,valout)
implicit double precision (a-h,o-z)
double precision           :: xslege(n), whtslege(n)
double complex             :: vals(:), valout
!
!  Use barycentric Lagrange interpolation to evaluate a complex-valued expansion of
!  the form (1) at a specified points.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
!      quadrature rule
!    vals - the *scaled* values of (1) at the n Gauss-Legendre nodes
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valout - the value of (1) at the point x
!
double precision :: whts(n)
double complex   :: dsum1

eps0  = epsilon(0.0d0)
dsign = 1.0d0
do i=1,n
dd      = x - xslege(i)
if ( abs(dd) .lt. eps0) then
valout = vals(i)/sqrt(whtslege(i))
return
endif
whts(i) = dsign/dd * sqrt(1.0d0-xslege(i)**2)
dsign   = -dsign
end do


dsum1 = 0
dsum2 = 0

do i=1,n
dsum1 = dsum1 + whts(i) * vals(i)
dsum2 = dsum2 + whts(i)*sqrt(whtslege(i))
end do

valout = dsum1/dsum2

end subroutine


subroutine legendre_interp2(n,xslege,whtslege,vals,x,valout)
implicit double precision (a-h,o-z)
double precision           :: xslege(n), whtslege(n)
double precision           :: vals(:), valout
!
!  Use barycentric Lagrange interpolation to evaluate a real-valued expansion of
!  the form (1) at a specified points.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
!      quadrature rule
!    vals - the *scaled* values of (1) at the n Gauss-Legendre nodes
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valout - the value of (1) at the point x
!
double precision :: whts(n)
double precision :: dsum1

eps0  = epsilon(0.0d0)
dsign = 1.0d0
do i=1,n
dd      = x - xslege(i)
if ( abs(dd) .lt. eps0) then
valout = vals(i)/sqrt(whtslege(i))
return
endif
whts(i) = dsign/dd * sqrt(1.0d0-xslege(i)**2)
dsign   = -dsign
end do


dsum1 = 0
dsum2 = 0

do i=1,n
dsum1 = dsum1 + whts(i) * vals(i)
dsum2 = dsum2 + whts(i)*sqrt(whtslege(i))
end do

valout = dsum1/dsum2

end subroutine

subroutine legendre_interp3(n,xslege,whtslege,vals,x,valsout)
implicit double precision (a-h,o-z)
double precision           :: xslege(n), whtslege(n)
double precision           :: vals(:,:), valsout(:)
!
!  Use barycentric Lagrange interpolation to evaluate a collection of real-valued
!  expansions of the form (1) at a specified point.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
!      quadrature rule
!    vals - a matrix whose jth column gives the *scaled* values (3) of the jth
!      input expanion
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valsout - an array whose jth entry gives the value of the jth input expansion
!      at the point x
!
double precision              :: whts(n)
double precision, allocatable :: sums(:)

l = size(vals,2)
allocate(sums(l))

! compute the barcyentric weights and sum them
eps0  = epsilon(0.0d0)
dsign = 1.0d0
do i=1,n
dd      = x - xslege(i)
if ( abs(dd) .lt. eps0) then
valsout = vals(i,:)/sqrt(whtslege(i))
return
endif
whts(i) = dsign/dd * sqrt(1.0d0-xslege(i)**2)
dsign   = -dsign
end do

! compute the sums
do i=1,l
vals(1:n,i) = vals(1:n,i) * whts 
sums(i)     = sum (vals(:,i))
vals(1:n,i) = vals(1:n,i) / whts 
end do

! evaluate the formula
whts    = whts * sqrt(whtslege)
valsout = sums/sum(whts)

end subroutine


subroutine legendre_interp4(n,xslege,whtslege,vals,x,valsout)
implicit double precision (a-h,o-z)
double precision           :: xslege(n), whtslege(n)
double complex             :: vals(:,:), valsout(:)
!
!  Use barycentric Lagrange interpolation to evaluate a collection of complex-valued
!  expansions of the form (1) at a specified point.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
!      quadrature rule
!    vals - a matrix whose jth column gives the *scaled* values (3) of the jth
!      input expanion
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    valsout - an array whose jth entry gives the value of the jth input expansion
!      at the point x
double precision              :: whts(n)
double complex, allocatable   :: sums(:)

l = size(vals,2)
allocate(sums(l))

! compute the barcyentric weights and sum them
eps0  = epsilon(0.0d0)
dsign = 1.0d0
do i=1,n
dd      = x - xslege(i)
if ( abs(dd) .lt. eps0) then
valsout = vals(i,:)/sqrt(whtslege(i))
return
endif
whts(i) = dsign/dd * sqrt(1.0d0-xslege(i)**2)
dsign   = -dsign
end do

! compute the sums
do i=1,l
vals(1:n,i) = vals(1:n,i) * whts 
sums(i)     = sum (vals(:,i))
vals(1:n,i) = vals(1:n,i) / whts 
end do

whts    = whts * sqrt(whtslege)
valsout = sums/sum(whts)

end subroutine



subroutine legendre_coefs1(n,xslege,whtslege,vals,coefs)
implicit double precision (a-h,o-z)
double precision            :: xslege(n), whtslege(n)
double complex              :: vals(:),coefs(:)
!
!  Given the vector (3) of scaled values of a complex-valued expansion of the form
!  (1), compute the vector (2) of expansion coefficients.
!  
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    xslege - the nodes of the n-point Gauss-Legendre quadrature
!    whtslege - the weights of the n-point Gauss-Legendre quadrature
!    vals - the vector (3) of scaled values
!
!  Output parameters:
!    coefs - the vector (2) of coefficients
!
double precision :: pols(n)

coefs = 0

do j=1,n
x   = xslege(j)
wht = whtslege(j)
call leges(n,x,pols)
pols  = pols * sqrt(wht)
coefs = coefs + pols*vals(j)
end do

end subroutine


subroutine legendre_coefs2(n,xslege,whtslege,vals,coefs)
implicit double precision (a-h,o-z)
double precision            :: xslege(n), whtslege(n)
double precision            :: vals(:),coefs(:)
!
!  Given the vector (3) of scaled values of a real-valued expansion of the form
!  (1), compute the vector (2) of expansion coefficients.
!  
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    xslege - the nodes of the n-point Gauss-Legendre quadrature
!    whtslege - the weights of the n-point Gauss-Legendre quadrature
!    vals - the vector (3) of scaled values
!
!  Output parameters:
!    coefs - the vector (2) of coefficients
!

double precision :: pols(n)

coefs = 0

do j=1,n
x   = xslege(j)
wht = whtslege(j)
call leges(n,x,pols)
pols  = pols * sqrt(wht)
coefs = coefs + pols*vals(j)
end do

end subroutine


subroutine legendre_coefs3(n,xslege,whtslege,vals,coefs)
implicit double precision (a-h,o-z)
double precision            :: xslege(n), whtslege(n)
double precision            :: vals(:,:),coefs(:,:)
!
!  Given the vectors (3) of scaled values for a collection of real-valued expansions
!  of the form (1), compute the vectors (2) of their expansion coefficients.
!  
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    xslege - the nodes of the n-point Gauss-Legendre quadrature
!    whtslege - the weights of the n-point Gauss-Legendre quadrature
!    vals - a matrix the jth column of which gives the vector (3) of scaled
!      values of the jth input expansion
!
!  Output parameters:
!    coefs - a matrix the jth column of which gives the vector (2) of 
!      expansion coefficients for the jth input expansion
!

double precision :: pols(n)

l     = size(vals,2)
coefs = 0

do j=1,n
x   = xslege(j)
wht = whtslege(j)
call leges(n,x,pols)
pols  = pols * sqrt(wht)
do i=1,l
coefs(:,i) = coefs(:,i) + pols*vals(j,i)
end do
end do


end subroutine


subroutine legendre_coefs4(n,xslege,whtslege,vals,coefs)
implicit double precision (a-h,o-z)
double precision            :: xslege(n), whtslege(n)
double complex              :: vals(:,:),coefs(:,:)
!
!  Given the vectors (3) of scaled values for a collection of complex-valued
!  expansions of the form (!), compute the vectors (2) of their expansion coefficients.
!  
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    xslege - the nodes of the n-point Gauss-Legendre quadrature
!    whtslege - the weights of the n-point Gauss-Legendre quadrature
!    vals - a matrix the jth column of which gives the vector (3) of scaled
!      values of the jth input expansion
!
!  Output parameters:
!    coefs - a matrix the jth column of which gives the vector (2) of 
!      expansion coefficients for the jth input expansion
!
double precision :: pols(n)

l     = size(vals,2)
coefs = 0

do j=1,n
x   = xslege(j)
wht = whtslege(j)
call leges(n,x,pols)
pols  = pols * sqrt(wht)
do i=1,l
coefs(:,i) = coefs(:,i) + pols*vals(j,i)
end do
end do


end subroutine


subroutine legendre_eval1(n,coefs,x,val)
implicit double precision (a-h,o-z)
double complex                 :: coefs(:), val
!
!  Evaluate a complex-valued expansion of the form (1) given the vector
!  (2) of its  expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion at the point x
! 

double precision :: pols(n)

call leges(n,x,pols)

val = 0
do i=1,n
val = val + pols(i)*coefs(i)
end do

!
!  Use the Clenshaw algorithm ... it is probably better to use the
!  alternative algorithm above as there is then only one point of
!  failure when evaluating the Legendrep polynomials
!
! p0 = 1
! p1 = x

! !d0 = 0
! !d1 = 1

! dd1 = sqrt(0.5d0)
! dd2 = sqrt(1.5d0)

! val = coefs(1)*dd1+coefs(2)*dd2*x
! !der = coefs(2)

! do k=1,n-2
! p2 = (2*k+1.0d0)*x*p1-k*p0
! p2 = p2/(k+1.0d0)

! !d2 = (2*k+1.0d0)*x*d1+(2*k+1)*p1-k*d0
! !d2 = d2/(k+1.0d0)

! dd  = sqrt(k+1.5d0)
! val = val+p2*coefs(k+2)*dd
! !der = der+d2*coefs(k+2)

! p0 = p1
! p1 = p2

! !d0 = d1
! !d1 = d2
! end do

end subroutine


subroutine legendre_eval2(n,coefs,x,val)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:), val
!
!  Evaluate a real-valued expansion of the form (1) given the vector
!  (2) of its  expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion at the point x
! 

double precision :: pols(n)

call leges(n,x,pols)

val = 0
do i=1,n
val = val + pols(i)*coefs(i)
end do

!  Clenshaw algorithm
! p0 = 1
! p1 = x

! !d0 = 0
! !d1 = 1

! dd1 = sqrt(0.5d0)
! dd2 = sqrt(1.5d0)

! val = coefs(1)*dd1+coefs(2)*dd2*x
! !der = coefs(2)

! do k=1,n-2
! p2 = (2*k+1.0d0)*x*p1-k*p0
! p2 = p2/(k+1.0d0)

! !d2 = (2*k+1.0d0)*x*d1+(2*k+1)*p1-k*d0
! !d2 = d2/(k+1.0d0)

! dd  = sqrt(k+1.5d0)
! val = val+p2*coefs(k+2)*dd
! !der = der+d2*coefs(k+2)

! p0 = p1
! p1 = p2

! !d0 = d1
! !d1 = d2
! end do

end subroutine


subroutine legendre_eval3(n,coefs,x,vals)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:,:), vals(:)
!
!  Evaluate a collection of real-valued expansions of the form (1) given
!  their expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - a matrix whose jth column gives the coefficient vector (2)
!      for the jth input expansion
!    x - the point at which to evaluate the expansions
!
!  Output parameters:
!    vals - the jth entry of this array will give the value of the jth
!       input expansion
! 

double precision :: pols(n+1)

call leges(n,x,pols)

vals = 0
do i=1,n
vals = vals + pols(i)*coefs(i,:)
end do


! !
! !  Use the Clenshaw algorithm ...
! !

! p0 = 1
! p1 = x

! !d0 = 0
! !d1 = 1

! dd1 = sqrt(0.5d0)
! dd2 = sqrt(1.5d0)

! vals = coefs(1,:)*dd1+coefs(2,:)*dd2*x

! do k=1,n-2
! p2 = (2*k+1.0d0)*x*p1-k*p0
! p2 = p2/(k+1.0d0)

! !d2 = (2*k+1.0d0)*x*d1+(2*k+1)*p1-k*d0
! !d2 = d2/(k+1.0d0)

! dd   = sqrt(k+1.5d0)
! vals = vals+p2*coefs(k+2,:)*dd
! !der = der+d2*coefs(k+2)

! p0 = p1
! p1 = p2

! !d0 = d1
! !d1 = d2
! end do

end subroutine


subroutine legendre_eval4(n,coefs,x,vals)
implicit double precision (a-h,o-z)
double complex                  :: coefs(:,:), vals(:)
!
!  Evaluate a collection of complex-valued expansions of the form (1) given
!  their expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - a matrix whose jth column gives the coefficient vector (2)
!      for the jth input expansion
!    x - the point at which to evaluate the expansions
!
!  Output parameters:
!    vals - the jth entry of this array will give the value of the jth
!       input expansion
! 

double precision :: pols(n)

call leges(n,x,pols)

vals = 0
do i=1,n
vals = vals + pols(i)*coefs(i,:)
end do

! !
! !  Use the Clenshaw algorithm ... 
! !

! p0 = 1
! p1 = x

! !d0 = 0
! !d1 = 1

! dd1 = sqrt(0.5d0)
! dd2 = sqrt(1.5d0)

! vals = coefs(1,:)*dd1+coefs(2,:)*dd2*x

! do k=1,n-2
! p2 = (2*k+1.0d0)*x*p1-k*p0
! p2 = p2/(k+1.0d0)

! !d2 = (2*k+1.0d0)*x*d1+(2*k+1)*p1-k*d0
! !d2 = d2/(k+1.0d0)

! dd   = sqrt(k+1.5d0)
! vals = vals+p2*coefs(k+2,:)*dd
! !der = der+d2*coefs(k+2)

! p0 = p1
! p1 = p2

! !d0 = d1
! !d1 = d2
! end do

end subroutine


subroutine legendre_evalder1(n,coefs,x,val,der)
implicit double precision (a-h,o-z)
double complex                 :: coefs(:), val, der
!
!  Evaluate an expansion of the form (1) and its derivative given the expansion
!  coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion
! 

double precision :: pols(n), ders(n)

call legeders(n,x,pols,ders)

val = 0
der = 0

do i=1,n
val = val + pols(i)*coefs(i)
der = der + ders(i)*coefs(i)
end do

end subroutine


subroutine legendre_evalder2(n,coefs,x,val,der)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:)
!
!  Evaluate an expansion of the form (1) and its derivative given the expansion
!  coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - an array specifying the coefficients
!    x - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion
! 

double precision :: pols(n), ders(n)

call legeders(n,x,pols,ders)

val = 0
der = 0

do i=1,n
val = val + pols(i)*coefs(i)
der = der + ders(i)*coefs(i)
end do

end subroutine


subroutine legendre_evalder3(n,coefs,x,vals,ders)
implicit double precision (a-h,o-z)
double precision                 :: coefs(:,:), vals(:), ders(:)
!
!  Evaluate a collection of expansions of the form (1) and their derivatives at
!  a specified point x given the vector (3) of coefficients for each expansion.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - a matrix whose jth column gives the expansion coefficients for the
!      jth input expansion
!    x - the point at which to evaluate the expansions and their derivatives
!
!  Output parameters:
!    vals - an array whose jth entry will contain the value of the jth input
!      expansion at x
!    ders - an array whose jth entry will contain the value of the derivative
!      of the jth expansion at x
! 

double precision :: pols(n), ders0(n)

call legeders(n,x,pols,ders0)

vals = 0
ders = 0

do i=1,n
vals = vals + pols(i)*coefs(i,:)
ders = ders + ders0(i)*coefs(i,:)
end do

end subroutine


subroutine legendre_evalder4(n,coefs,x,vals,ders)
implicit double precision (a-h,o-z)
double complex                 :: coefs(:,:), vals(:), ders(:)
!
!  Evaluate a collection of expansions of the form (1) and their derivatives at
!  a specified point x given the vector (3) of coefficients for each expansion.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    coefs - a matrix whose jth column gives the expansion coefficients for the
!      jth input expansion
!    x - the point at which to evaluate the expansions and their derivatives
!
!  Output parameters:
!    vals - an array whose jth entry will contain the value of the jth input
!      expansion at x
!    ders - an array whose jth entry will contain the value of the derivative
!      of the jth expansion at x
! 

double precision :: pols(n), ders0(n)

call legeders(n,x,pols,ders0)

vals = 0
ders = 0

do i=1,n
vals = vals + pols(i)*coefs(i,:)
ders = ders + ders0(i)*coefs(i,:)
end do

end subroutine


subroutine legendre_coefsmatrix(n,xslege,whtslege,umatr)
implicit double precision (a-h,o-z)
double precision                           :: xslege(n), whtslege(n)
double precision, allocatable, intent(out) :: umatr(:,:)
!
!  Return the (n,n) matrix which takes the *scaled* values of an expansion of the
!  form (1) at the nodes of the n-point Gauss-Legendre quadrature to the
!  expansion coefficients.
!
!  Input parameters:
!    n - the number of terms in the expansion
!    xslege - the nodes of the n-point Gauss-Legendre quadrature rule
!    whtslege - the weights of the n-point Gauss-Legendre rule
!
!  Output parameters:
!    umatr - the (n,n) matrix which takes values to coefficients
! 

allocate(umatr(n,n))

do i=1,n
call leges(n,xslege(i),umatr(:,i))
end do

do j=1,n
umatr(j,:) = umatr(j,:) * sqrt(whtslege)
end do

end subroutine


subroutine legendre_interpmatrix(n,xslege,whtslege,m,xsout,whtsout,ainterp)
implicit double precision (a-h,o-z)
dimension xslege(n), whtslege(n), xsout(m), whtsout(m)
double precision, allocatable, intent(out)  :: ainterp(:,:)
!
!  Construct the matrix which takes the *scaled* values of (1) at the Gauss-Legendre
!  nodes to its scaled values at the nodes of a user-specified quadrature rule.
!
!  Input parameters:
!    n - the number of terms in the expansion (1)
!    (xslege,whtslege) - the nodes and weights of the n-point Guass-Legendre
!      quadrature rule
!    m - the number of points in the output target quadrature rule
!    (xsout,whtsout) - the nodes and weights of the output quadrature rule
!
!  Output parameters:
!    ainterp - the (m,n) interpolation matrix
!

dimension whts(n)

eps0  = epsilon(0.0d0)

allocate(ainterp(m,n))
ainterp = 0

do j=1,m
x     = xsout(j)
wht   = whtsout(j)
idx   = 0

dsign = 1.0d0
do i=1,n
dd      = x - xslege(i)
if ( abs(dd) .lt. eps0) then
idx = i 
exit
endif

whts(i) = dsign/dd * sqrt( ( 1.0d0-xslege(i)**2 ) )
! * whtslege(i)
dsign   = -dsign
end do

if (idx .ne. 0) then
ainterp(j,idx) = 1 / sqrt(whtslege(idx)) * sqrt(wht)
cycle
endif

dsum = sum(whts*sqrt(whtslege))

do i=1,n
ainterp(j,i) = whts(i)/dsum * sqrt(wht)
end do

end do

end subroutine


subroutine legendre_diffmatrix(n,xslege,whtslege,adiff)
implicit double precision (a-h,o-z)
double precision                           :: xslege(:), whtslege(:)
double precision, allocatable, intent(out) :: adiff(:,:)
!
!  Return the spectral differentiation matrix which takes the scaled
!  values of an expansion of the form (1) at the Gauss-Legendre nodes
!  to the scaled values of its derivatves at the same nodes.
!
!  Input parameters:
!    n - the number of terms in the Legendre expansion
!    (xs,whts) - the n-point Gauss-Legendre quadrature rule
!
!  Output parameters:
!    adiff - the (n,n) spectral differentiation matrix
!
double precision         :: pols(0:n),ders(0:n),u(n,n),vals(n,n)

allocate(adiff(n,n))

do i=1,n
x   = xslege(i)
wht = whtslege(i)
call legeders(n-1,x,pols,ders)
do j=1,n
u(j,i)    = pols(j-1)*sqrt(wht)
vals(i,j) = ders(j-1)*sqrt(wht)
end do
end do

adiff = matmul(vals,u)

end subroutine


subroutine legeder(n,x,pol,der)
implicit double precision (a-h,o-z)
integer n
double precision x,pol
!
!  Evaluate the L^2 normalized Legendre polynomial of degree n and its derivative
!  at the point x.
!
!  Input parameters:
!    n - the degree of the polynomial to evaluate
!    x - the point at which the polynomial is to be evaluated
!
!  Output parameters:
!    pol - the value of P_n(x)
!    der - the value of P_n'(x)
!

! if (abs(x-1.0d0) .lt. 1.0d-3) then
! dn = n
! call lege_p1(dn,x,pol,der)
! return
! endif

if ( x == 1.0d0) then
pol = 1
der = (n)*(n+1.0d0)/2
goto 1000
endif

if ( x == -1.0d0) then
pol = (-1.0d0)**n
der = -pol*(n)*(n+1.0d0)/2
goto 1000
endif


if (n == 0) then
pol = 1
der = 0
else if (n == 1) then
pol = x
der = 1.0d0
else
p1 = 1
p2 = x

do j=2,n
   p  = ((2*j-1)*x*p2-(j-1)*p1)/j
   p1 = p2
   p2 = p
end do
!
pol = p

!
! Compute the derivative using another well-known formula.
!
der=n*(x*p2-p1)/(x**2-1.0d0)
endif

!
!  Normalize the polynomial
!

1000 continue

dd  = sqrt((n + 0.5d0))
pol = pol * dd
der = der * dd

end subroutine



subroutine leges(n,x,pols)
implicit double precision (a-h,o-z)
double precision           :: pols(:)
!
!  Return the values of the L^2 normalized Legendre polynomials of the
!  first kind of degrees 0 through n-1 at a specified point.
!
!  Input parameters:
!     n - the number of polynomials to evaluate
!     x - the point at which to evaluate them
!
!  Output parameters:
!     pols - the array of length n containing the values of the first n 
!      Legendre polynomials
!

if (x == 1.0d0) then
do i=1,n
pols(i) = 1.0d0
end do
goto 1000
endif

if (x == -1.0d0) then
dsign = 1.0d0
do i=1,n
pols(i) = dsign
dsign   = -dsign
end do
goto 1000
endif


pols(1) = 1.0d0
if (n == 1) goto 1000
pols(2) = x
if (n == 2) goto 1000

do j=2,n-1
pols(j+1) = ((2*j-1)*x*pols(j)-(j-1)*pols(j-1))/j
end do

!
!  Normalize the polynomials
!
1000 continue

do j=1,n
dd      = sqrt(j - 0.5d0)
pols(j) = pols(j) * dd
end do



end subroutine


subroutine legeders(n,x,pols,ders)
implicit double precision (a-h,o-z)
integer                       :: n
double precision              :: x
double precision, intent(out) :: pols(:),ders(:)
!
!  Evaluate the n Legendre polynomials of degree 0 through n-1 at the
!  point x using the standard 3-term recurrence relation.  Return the values
!  of their derivative at the point x as well.
!
!  Input parameters:
!    n - an integer specifying the number of polynomials to be evaluated
!    x - the point at which to evaluate the polynomials and their derivatives 
!
!  Output parameters:
!    pols - the ith entry of this user-supplied array of length n
!      will contain the value of normalized Legendre polynomial of degree
!      i-1 at x
!    ders - the ith entry of this user-supplied array will contain the
!      value of the derivative of the normalized Legendre polynomial at x
!
!
double precision :: pols2(n)


if ( x == 1.0d0) then
do i=1,n
pols(i) = 1.0d0
ders(i) = (i-1)*i/2
end do
goto 1000
endif

if ( x == -1.0d0) then
dsign = 1.0d0
do i=1,n
pols(i) = dsign
ders(i) = -dsign*(i-1)*i/2
dsign   = -dsign
end do
goto 1000
endif

pols(1)  = 1
ders(1)  = 0

if (n .gt. 1) then
pols(2)  = x
ders(2)  = 1
end if

!
!  Calculate the values of the unnormalized polynomials
!
do j=2,n-1
pols(j+1)=((2.0d0*j-1.0d0)*x*pols(j)-(j-1.0d0)*pols(j-1))/j
end do

!
!  Compute the derivatives of the unnormalized polynomials
!
d=x**2.0d0-1.0d0
do j=3,n
ders(j)=(j-1.0d0)*(x*pols(j)-pols(j-1))/d
end do

!
!  Normalize the polynomials and scale the derivatives
!
1000 continue

do i=1,n
dd      = sqrt(i-0.5d0)
pols(i) = pols(i) * dd
ders(i) = ders(i) * dd
end do

end subroutine



subroutine legequad2(n,xs,whts)
implicit double precision (a-h,o-z)
dimension xs(n),whts(n)
data pi /3.14159265358979323846264338327950288d0/
!
!  Use Newton's method and local Taylor expansions (i.e., the Glaser-Liu Rokhlin method)
!  to compute the n roots of the Legendre polynomial of degree n.  This routine 
!  is O(n) in the number of nodes n.

maxiters = 12
eps0     = epsilon(0.0d0)*10
k        = 30
!
!  Find the roots in [0,1]
!
        ifodd = 0
        nn = (n+1)/2
        if (nn .ne. n/2) then
        ifodd=1
        nn=nn-1
        endif
!
!  Use a negative value of x0 to indicate that the procedure has
!  not been initialized.
!
!
        x0 = -1
!
        do 2000 i=nn+1,n
!
!       Use Chebyshev node as initial guess for roots ...
!
!        x1 = cos(-pi+(2*i-1)*pi/(2*n))
!
!       ... or use this somewhat better approximation.
!
        dk = i-nn-1
        dn = n
        theta = (4*(ceiling(dn/2)-dk)-1) / (4*dn+2) * pi

        x1 = 1.0d0 - (dn-1)/(8*dn**3)-1.0d0/(384.0d0*dn**4)*   &
            (39.0d0-28.0d0/sin(theta)**2)
!        
        x1=x1*cos(theta)
!
!       Conduct Newton iterations.
!
        do 2100 iter=1,maxiters
!
!       Evaluate the nth polynomial at x1 using a Taylor expansion
!       and the recurrence relation.  The first evaluation must be
!       handled specially.
!
        if (x0 .lt. 0) then
        call legeder(n,x1,pol,der)
        else
        dx = x1-x0       
        call legetayl(n,k,x0,dx,pol,der)
        endif
!
!       Newton iteration.
!
        x0 = x1
        dd = -pol/der
        x1 = x0+dd
!
        if (abs(dd) .lt. eps0) then
        xs(i)=x1
        whts(i)=der
        goto 2000
        endif
!
 2100 continue
!
!       This doesn't happen.
!
        print *,"legeroots bombed!"
        stop
 2000 continue
!
!       Compute the weights using the derivatives we stored above.
!        
        do 3000 j=nn+1,n
        x       = xs(j)
        dd      = 2.0d0/(1-x**2)
        whts(j) = dd/whts(j)**2
 3000 continue
!
!       Reflect the quadrature on [-1,0] to [0,1].
!
        do 4000 j=1,nn
        xs(j)   = -xs(n-j+1)
        whts(j) = whts(n-j+1)
 4000 continue
!
!       Handle the root at 0 if n is odd.
!
        if (ifodd .eq. 1) then
        x0          = 0
        call legeder(n,x0,pol,der)
!
        xs(nn+1)   = x0
        whts(nn+1) = 2.0d0/der**2
        endif
!
end subroutine


        subroutine legetayl(n,k,x,dx,pol,der)
        implicit double precision (a-h,o-z)
!
!       Evaluate the Legendre polynomial of order n using a Taylor
!       expansion and the recurrence relation for the power series
!       coefficients.
!
        k0 = min(k,n)
        p1     = pol
        p2     = der*dx
        sum    = p1+p2
        sumder = p2/dx
        do 1000 j=0,k0-2
        d = 2*x*(j+1.0d0)**2/dx*p2-(n*(n+1.0d0)-j*(j+1.0d0))*p1
        d = d/(j+1.0d0)/(j+2.0d0)*dx**2/(1.0d0-x**2)
        sum    = sum+d
        sumder = sumder+d*(j+2)/dx
        p1 = p2
        p2 = d
 1000 continue
        pol = sum
        der = sumder
        end subroutine


end module
