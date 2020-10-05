!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This file contains code for constructing and manipulating bivariate Legendre
!  expansions --- that is, expressions of the form
!
!               
!    f(x,y) = \sum          a_ij \tilde{P}_i(x) \tilde{P}_j(y),                            (1)
!            0 <= i+j <= n
!
!  where \tilde{P}_i denotes the L^2 normalized Legendre polynomial of degree i.
!  The expansions (1) can be either real or complex-valued.
!
!  Expansions of the form (1) are represented in two different ways: via the
!  vector
!
!     ( a_00         )
!     ( a_10         )
!     ( a_01         )
!     ( a_20         )
!     ( a_11         )                                                                     (2)
!     ( a_02         )
!     ( a_30         )
!        ...
!     ( a_n0         )
!     ( a_0n         )
!
!  of coefficients (note the ordering of the coefficients), and via the vector
!
!     ( f(x_1,y_1) \sqrt{w1} )
!     ( f(x_2,y_1) \sqrt{w1} )
!     ( f(x_3,y_1) \sqrt{w1} )
!            ...
!     ( f(x_n,y_1) \sqrt{w1} )                                                             (3)
!     ( f(x_2,y_1) \sqrt{w1} )                                                             
!     ( f(x_2,y_2) \sqrt{w1} )                                                            
!            ...
!     ( f(x_n,y_n) \sqrt{w1} ),
!
!  where (x_1,y_1),...(x_n,y_n), w_1,...,w_n are the nodes and weights of a
!  generalized Gaussian quadrature rule which integrates polynomials of degree
!  somewhat higher than n-1.
!
!  The quadrature rules and the matrices which takes the vector (3) to (2) are
!  precomputed (the precomputed rules and matrices are stored in the file 
!  sqquads.f90).  Consequently, only certain possible values of the order n, 
!  which controls the size of the expansion, are permitted (see the subroutine
!  bilege_quad).
!
!  The following subroutines are public:
!
!    bilege_quad - return the generaized Gaussian quadrature rule used to
!     represent expansions of the form (1)
!
!    bilege_interpmatrix - return the matrix taking the vector (3) of *scaled*
!      values of the expansion (1) to the *scaled* values of (1) at the nodes
!      of a user-specified quadrature rule
!
!    bilege_coefsmatrix - return the matrix taking the vector (3) of *scaled*
!      values of the expansion (1) to the vector (2) of coefficients -- note that
!      this matrix is orthogonal so its transpose is the matrix which takes (2) to 
!      (3)
!
!    bilege_coefs - given the vector (3) of scaled values of an expansion of
!      the form (1), compute the vector (2) of coefficients in the series (1)
!
!    bilege_eval - evaluate one or more expansions of the form (1) given the 
!      their expansion coefficients 
!
!    bilege_evalder - evaluate one or more expansions of the form (1) and their
!      derivatives given the expansions' coefficients
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bilege

use utils
use legendre
use sqquads

interface      bilege_coefs
module procedure bilege_coefs1
module procedure bilege_coefs2
module procedure bilege_coefs3
module procedure bilege_coefs4
end interface  bilege_coefs

interface      bilege_eval
module procedure bilege_eval1
module procedure bilege_eval2
module procedure bilege_eval3
module procedure bilege_eval4
end interface  bilege_eval

interface      bilege_evalder
module procedure bilege_evalder1
module procedure bilege_evalder2
module procedure bilege_evalder3
module procedure bilege_evalder4
end interface  bilege_evalder


contains

subroutine bilege_quad(n,nquad,xs,ys,whts)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out) ::  xs(:),ys(:),whts(:)
!
!  Return a quadrature rule used to discretize expansions of a given
!  order.
!
!  Input parameters:
!    n - the desired order of the expansion (1), which must be 4, 6, 8, 10, 12,
!      14, 16, 18, 20, 22, 24 or 30
!
!  Output parameters:
!    nquad - the number of points in the discretization quadrature rule
!    xs - an array specifying the x-coordinates of the quadrature nodes
!    ys - an array specifying the y-coordinates of the quadrature nodes
!    whts - an array specifying the quadrature weights
!
!double precision, allocatable :: xslege(:), whtslege(:)

if (n == 4) then
nquad  = 15
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs4
ys    = sqys4
whts  = sqwhts4
return
endif

if (n == 6) then
nquad  = 28
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs6
ys    = sqys6
whts  = sqwhts6
return
endif

if (n == 8) then
nquad  = 45
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs8
ys    = sqys8
whts  = sqwhts8
return
endif

if (n == 10) then
nquad  = 66
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs10
ys    = sqys10
whts  = sqwhts10
return
endif

if (n == 12) then
nquad  = 91
!m      = 20
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs12
ys    = sqys12
whts  = sqwhts12
return
endif

if (n == 14) then
nquad  = 120
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs14
ys    = sqys14
whts  = sqwhts14
return
endif

if (n == 16) then
nquad  = 153
!m      = 27
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs16
ys    = sqys16
whts  = sqwhts16
return
endif

if (n == 18) then
nquad  = 190
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs18
ys    = sqys18
whts  = sqwhts18
return
endif

if (n == 20) then
nquad  = 231
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs20
ys    = sqys20
whts  = sqwhts20
return
endif


if (n == 22) then
nquad  = 276
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs22
ys    = sqys22
whts  = sqwhts22
return
endif

if (n == 24) then
nquad  = 325
!m      = 34
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs24
ys    = sqys24
whts  = sqwhts24
return
endif

if (n == 30) then
nquad  = 496
allocate( xs(nquad), ys(nquad), whts(nquad) )
xs    = sqxs30
ys    = sqys30
whts  = sqwhts30
return
endif


call prini("bilege_quad: invalid order n = ",n)
stop

end subroutine


subroutine bilege_interpmatrix(n,nout,xsout,ysout,whtsout,amatr)
implicit double precision (a-h,o-z)
double precision                            :: xsout(nout),ysout(nout),whtsout(nout)
double precision, allocatable, intent(out)  :: amatr(:,:)
!
!  Construct the (nout,nquad) matrix which takes the vector (3) of scaled values
!  of an expansion of the form (1) at the quadrature nodes to the scaled
!  values of (1) at the nodes of a user-specified quadrature rule.
!
!  Input parameters:
!    n - the order of the bivarite expansion
!    nout - the number of nodes in the target quadrature rule
!    xsout - the x-coords of the target quadrature nodes
!    ysout - the y-coords of the target quadrature nodes
!    whtsout - the weights of the target quadrature rule
!
!  Output parameters:
!    amatr - the (nout,nquad) interpolation matrix
!

double precision              :: polsx(0:n),polsy(0:n)
double precision, allocatable :: bmatr(:,:)

nquad = 0

if (n==4)   nquad=15
if (n==6)   nquad=28
if (n==8)   nquad=45
if (n==10)  nquad=66
if (n==12)  nquad=91
if (n==14)  nquad=120
if (n==16)  nquad=153
if (n==18)  nquad=190
if (n==20)  nquad=231
if (n==22)  nquad=276
if (n==24)  nquad=325
if (n==30)  nquad=496

if (nquad==0) then
call prini("bilege_intmatrix: invalid order n = ",n)
stop
endif

allocate(bmatr(nout,nquad))

!
!  Build the matrix consisting of the scaled values of the basis functions 
!  at the nodes of the output quadrature
!

do i=1,nout
x   = xsout(i)
y   = ysout(i)
wht = whtsout(i)

call leges(n+1,x,polsx)
call leges(n+1,y,polsy)

j=0
do nn=0,n
do jj=0,nn
ii=nn-jj
j=j+1
bmatr(i,j) = polsx(ii)*polsy(jj)*sqrt(wht)
end do
end do
end do

if (n==4)   amatr = matmul(bmatr,squmatr4)
if (n==6)   amatr = matmul(bmatr,squmatr6)
if (n==8)   amatr = matmul(bmatr,squmatr8)
if (n==10)  amatr = matmul(bmatr,squmatr10)
if (n==12)  amatr = matmul(bmatr,squmatr12)
if (n==14)  amatr = matmul(bmatr,squmatr14)
if (n==16)  amatr = matmul(bmatr,squmatr16)
if (n==18)  amatr = matmul(bmatr,squmatr18)
if (n==20)  amatr = matmul(bmatr,squmatr20)
if (n==22)  amatr = matmul(bmatr,squmatr22)
if (n==24)  amatr = matmul(bmatr,squmatr24)
if (n==30)  amatr = matmul(bmatr,squmatr30)

end subroutine



subroutine bilege_coefsmatrix(n,umatr)
implicit double precision (a-h,o-z)
double precision, allocatable, intent(out) :: umatr(:,:)
!
!  Return the (nquad,nquad) matrix which takes the vector (3) of scaled
!  values of the expansion (1) to the vector (2) of expansion coefficients.
!
!  Input parameters:
!    n - the order of the expansion (1)
!
!  Output parameters:
!    umatr - the (nquad,nquad) matrix which takes (3) to (2)
!



if (n==4) then
allocate(umatr(15,15))
umatr = squmatr4
return
endif

if (n==6) then
allocate(umatr(28,28))
umatr = squmatr6
return
endif

if (n==8) then
allocate(umatr(45,45))
umatr = squmatr8
return
endif

if (n==10) then
allocate(umatr(66,66))
umatr = squmatr10
return
endif

if (n==12) then
allocate(umatr(91,91))
umatr = squmatr12
return
endif

if (n==14) then
allocate(umatr(120,120))
umatr = squmatr14
return
endif


if (n==16) then
allocate(umatr(153,153))
umatr = squmatr16
return
endif

if (n==18) then
allocate(umatr(190,190))
umatr = squmatr18
return
endif

if (n==20) then
allocate(umatr(231,231))
umatr = squmatr20
return
endif

if (n==22) then
allocate(umatr(276,276))
umatr = squmatr22
return
endif

if (n==24) then
allocate(umatr(325,325))
umatr = squmatr24
return
endif

if (n==30) then
allocate(umatr(496,496))
umatr = squmatr30
return
endif

call prini("bilege_coefsmatrix: invalid order n = ",n)
stop

end subroutine


subroutine bilege_coefs1(n,vals,coefs)
implicit double precision (a-h,o-z)
double complex                             :: vals(:),coefs(:)
!
!  Compute the coefficients in a complex-valued expansion of the form (1) given the 
!  vector (3) of its *scaled* values at the nodes of the tensor product
!  quadrature rule.
!
!  Input parameters:
!    n - the order of the expansion
!    vals - the vector (3) of scaled values of the expansion
!
!  Output parameters:
!    coefs - the vector (2) of expansion coefficients
!

double precision         :: polsx(n), polsy(n)

if (n == 4) then
coefs = matmul(squmatr4,vals) 
return
endif

if (n == 6) then
coefs = matmul(squmatr6,vals) 
return
endif

if (n == 8) then
coefs = matmul(squmatr8,vals) 
return
endif

if (n == 10) then
coefs = matmul(squmatr10,vals) 
return
endif

if (n == 12) then
coefs = matmul(squmatr12,vals) 
return
endif

if (n == 14) then
coefs = matmul(squmatr14,vals) 
return
endif

if (n == 16) then
coefs = matmul(squmatr16,vals)
return
endif

if (n == 18) then
coefs = matmul(squmatr18,vals)
return
endif

if (n == 20) then
coefs = matmul(squmatr20,vals) 
return
endif

if (n == 22) then
coefs = matmul(squmatr22,vals) 
return
endif

if (n == 24) then
coefs = matmul(squmatr24,vals) 
return
endif

if (n == 30) then
coefs = matmul(squmatr30,vals) 
return
endif

call prini("bilege_coefs: invalid order n = ",n)
stop

end subroutine


subroutine bilege_coefs2(n,vals,coefs)
implicit double precision (a-h,o-z)
double complex                             :: vals(:,:),coefs(:,:)
!
!  Compute the coefficients in a collection of expansions (1) given the matrix
!  whose columns are the vectors  (3) of *scaled* values of the expansions at the 
!  nodes of the tensor product quadrature rule.
!
!  Input parameters:
!    n - the order of the expansion
!    vals - the matrix whose jth column is the vector (3) of scaled values of the
!      jth input expansion
!
!  Output parameters:
!    coefs - the matrix whose jth column is the vector (2) of expansion coefficients
!      of the jth input expansion
!

double precision         :: polsx(n), polsy(n)

if (n == 4) then
coefs = matmul(squmatr4,vals) 
return
endif

if (n == 6) then
coefs = matmul(squmatr6,vals) 
return
endif

if (n == 8) then
coefs = matmul(squmatr8,vals) 
return
endif

if (n == 10) then
coefs = matmul(squmatr10,vals) 
return
endif

if (n == 12) then
coefs = matmul(squmatr12,vals) 
return
endif

if (n == 14) then
coefs = matmul(squmatr14,vals) 
return
endif

if (n == 16) then
coefs = matmul(squmatr16,vals)
return
endif

if (n == 18) then
coefs = matmul(squmatr18,vals)
return
endif

if (n == 20) then
coefs = matmul(squmatr20,vals) 
return
endif

if (n == 22) then
coefs = matmul(squmatr22,vals) 
return
endif

if (n == 24) then
coefs = matmul(squmatr24,vals) 
return
endif

if (n == 30) then
coefs = matmul(squmatr30,vals) 
return
endif

call prini("bilege_coefs: invalid order n = ",n)
stop

end subroutine


subroutine bilege_coefs3(n,vals,coefs)
implicit double precision (a-h,o-z)
double precision                             :: vals(:),coefs(:)
!
!  Compute the coefficients in a real-valued expansion of the form (1) given the 
!  vector (3) of its *scaled* values at the nodes of the tensor product
!  quadrature rule.
!
!  Input parameters:
!    n - the order of the expansion
!    vals - the vector (3) of scaled values of the expansion
!
!  Output parameters:
!    coefs - the vector (2) of expansion coefficients
!

double precision         :: polsx(n), polsy(n)

if (n == 4) then
coefs = matmul(squmatr4,vals) 
return
endif

if (n == 6) then
coefs = matmul(squmatr6,vals) 
return
endif

if (n == 8) then
coefs = matmul(squmatr8,vals) 
return
endif

if (n == 10) then
coefs = matmul(squmatr10,vals) 
return
endif

if (n == 12) then
coefs = matmul(squmatr12,vals) 
return
endif

if (n == 14) then
coefs = matmul(squmatr14,vals) 
return
endif

if (n == 16) then
coefs = matmul(squmatr16,vals)
return
endif

if (n == 18) then
coefs = matmul(squmatr18,vals)
return
endif

if (n == 20) then
coefs = matmul(squmatr20,vals) 
return
endif

if (n == 22) then
coefs = matmul(squmatr22,vals) 
return
endif

if (n == 24) then
coefs = matmul(squmatr24,vals) 
return
endif

if (n == 30) then
coefs = matmul(squmatr30,vals) 
return
endif

call prini("bilege_coefs: invalid order n = ",n)
stop

end subroutine


subroutine bilege_coefs4(n,vals,coefs)
implicit double precision (a-h,o-z)
double precision                             :: vals(:,:),coefs(:,:)
!
!  Compute the coefficients in a collection of expansions (1) given the matrix
!  whose columns are the vectors  (3) of *scaled* values of the expansions at the 
!  nodes of the tensor product quadrature rule.
!
!  Input parameters:
!    n - the order of the expansion
!    vals - the matrix whose jth column is the vector (3) of scaled values of the
!      jth input expansion
!
!  Output parameters:
!    coefs - the matrix whose jth column is the vector (2) of expansion coefficients
!      of the jth input expansion
!

double precision         :: polsx(n), polsy(n)

if (n == 4) then
coefs = matmul(squmatr4,vals) 
return
endif

if (n == 6) then
coefs = matmul(squmatr6,vals) 
return
endif

if (n == 8) then
coefs = matmul(squmatr8,vals) 
return
endif

if (n == 10) then
coefs = matmul(squmatr10,vals) 
return
endif

if (n == 12) then
coefs = matmul(squmatr12,vals) 
return
endif

if (n == 14) then
coefs = matmul(squmatr14,vals) 
return
endif

if (n == 16) then
coefs = matmul(squmatr16,vals)
return
endif

if (n == 18) then
coefs = matmul(squmatr18,vals)
return
endif

if (n == 20) then
coefs = matmul(squmatr20,vals) 
return
endif

if (n == 22) then
coefs = matmul(squmatr22,vals) 
return
endif

if (n == 24) then
coefs = matmul(squmatr24,vals) 
return
endif

if (n == 30) then
coefs = matmul(squmatr30,vals) 
return
endif

call prini("bilege_coefs: invalid order n = ",n)
stop

end subroutine



subroutine bilege_eval1(n,coefs,x,y,val)
implicit double precision (a-h,o-z)
double precision                :: xslege(n), whtslege(n)
double complex                  :: coefs(:),val
!
!  Given the vector (2) of coefficients of a complex-valued expansion of the
!  form (1), evaluate the expansion at a specified point.
!
!  Input parameters:
!    n - the order of the bivariate expansion
!    coefs - the vector (2) of expansion coefficients
!    (x,y) - the point at which to evaluate the expansion (1)
!
!  Output parameters:
!   val - the value of the expansion
!

double precision          :: polsx(0:n), polsy(0:n)

call leges(n+1,x,polsx)
call leges(n+1,y,polsy)

val = 0
idx = 1
do nn=0,n
do jj=0,nn
ii = nn-jj
val = val + coefs(idx)*polsx(ii)*polsy(jj)
idx = idx+1
end do
end do

end subroutine

subroutine bilege_eval2(n,coefs,x,y,val)
implicit double precision (a-h,o-z)
double precision                :: xslege(n), whtslege(n)
double precision                :: coefs(:),val
!
!  Given the vector (2) of coefficients of a real-valued expansion of the
!  form (1), evaluate the expansion at a specified point.
!
!  Input parameters:
!    n - the order of the bivariate expansion
!    coefs - the vector (2) of expansion coefficients
!    (x,y) - the point at which to evaluate the expansion (1)
!
!  Output parameters:
!   val - the value of the expansion
!

double precision          :: polsx(0:n), polsy(0:n)

call leges(n+1,x,polsx)
call leges(n+1,y,polsy)

val = 0
idx = 1
do nn=0,n
do jj=0,nn
ii = nn-jj
val = val + coefs(idx)*polsx(ii)*polsy(jj)
idx = idx+1
end do
end do

end subroutine


subroutine bilege_eval3(n,coefs,x,y,vals)
implicit double precision (a-h,o-z)
double precision                :: xslege(n), whtslege(n)
double precision                :: coefs(:,:),vals(:)
!
!  Given the expansion coefficients of a collection of real-valued expansions of the
!  form (1), evaluate each  at a specified point.
!
!  Input parameters:
!    n - the order of the bivariate expansion
!    coefs - the vector (2) of expansion coefficients
!    (x,y) - the point at which to evaluate the expansion (1)
!
!  Output parameters:
!   val - the value of the expansion
!

double precision          :: polsx(0:n), polsy(0:n)

call leges(n+1,x,polsx)
call leges(n+1,y,polsy)

vals = 0
idx = 1
do nn=0,n
do jj=0,nn
ii   = nn-jj
vals = vals + coefs(idx,:)*polsx(ii)*polsy(jj)
idx  = idx+1
end do
end do

end subroutine

subroutine bilege_eval4(n,coefs,x,y,vals)
implicit double precision (a-h,o-z)
double precision                :: xslege(n), whtslege(n)
double complex                  :: coefs(:,:),vals(:)
!
!  Given the expansion coefficients of a collection of real-valued expansions of the
!  form (1), evaluate each  at a specified point.
!
!  Input parameters:
!    n - the order of the bivariate expansion
!    coefs - the vector (2) of expansion coefficients
!    (x,y) - the point at which to evaluate the expansion (1)
!
!  Output parameters:
!   val - the value of the expansion
!

double precision          :: polsx(0:n), polsy(0:n)

call leges(n+1,x,polsx)
call leges(n+1,y,polsy)

vals = 0
idx = 1
do nn=0,n
do jj=0,nn
ii   = nn-jj
vals = vals + coefs(idx,:)*polsx(ii)*polsy(jj)
idx  = idx+1
end do
end do

end subroutine


subroutine bilege_evalder1(n,coefs,x,y,val,derx,dery)
implicit double precision (a-h,o-z)
double precision                :: xslege(n), whtslege(n)
double complex                  :: coefs(:),val,derx,dery
!
!  Given the vector (2) of coefficients of a complex-valued expansion of the
!  form (1), evaluate the expansion at a specified point.
!
!  Input parameters:
!    n - the order of the bivariate expansion
!    coefs - the vector (2) of expansion coefficients
!    (x,y) - the point at which to evaluate the expansion (1)
!
!  Output parameters:
!   val - the value of the expansion
!

double precision          :: polsx(0:n), polsy(0:n)
double precision          :: dersx(0:n), dersy(0:n)

call legeders(n+1,x,polsx,dersx)
call legeders(n+1,y,polsy,dersy)

val  = 0
derx = 0
dery = 0
idx  = 1
do nn=0,n
do jj=0,nn
ii = nn-jj
val  = val  + coefs(idx)*polsx(ii)*polsy(jj)
derx = derx + coefs(idx)*dersx(ii)*polsy(jj)
dery = dery + coefs(idx)*polsx(ii)*dersy(jj)
idx = idx+1
end do
end do

end subroutine



subroutine bilege_evalder2(n,coefs,x,y,val,derx,dery)
implicit double precision (a-h,o-z)
double precision                :: xslege(n), whtslege(n)
double precision                :: coefs(:),val,derx,dery
!
!  Given the vector (2) of coefficients of a real-valued expansion of the
!  form (1), evaluate the expansion at a specified point.
!
!  Input parameters:
!    n - the order of the bivariate expansion
!    coefs - the vector (2) of expansion coefficients
!    (x,y) - the point at which to evaluate the expansion (1)
!
!  Output parameters:
!   val - the value of the expansion
!

double precision          :: polsx(0:n), polsy(0:n)
double precision          :: dersx(0:n), dersy(0:n)

call legeders(n+1,x,polsx,dersx)
call legeders(n+1,y,polsy,dersy)

val  = 0
derx = 0
dery = 0
idx  = 1
do nn=0,n
do jj=0,nn
ii = nn-jj
val  = val  + coefs(idx)*polsx(ii)*polsy(jj)
derx = derx + coefs(idx)*dersx(ii)*polsy(jj)
dery = dery + coefs(idx)*polsx(ii)*dersy(jj)
idx = idx+1
end do
end do

end subroutine


subroutine bilege_evalder3(n,coefs,x,y,vals,dersx,dersy)
implicit double precision (a-h,o-z)
double precision                :: xslege(n), whtslege(n)
double precision                :: coefs(:,:),vals(:),dersx(:),dersy(:)
!
!  Given the coefficients of a collection of real-valued expansions of the
!  form (1), evaluate them and their derivatives  at a specified point.
!
!  Input parameters:
!    n - the order of the bivariate expansion
!    coefs - the vector (2) of expansion coefficients
!    (x,y) - the point at which to evaluate the expansion (1)
!
!  Output parameters:
!   val - the value of the expansion
!

double precision          :: polsx(0:n), polsy(0:n)
double precision          :: poldersx(0:n), poldersy(0:n)

call legeders(n+1,x,polsx,poldersx)
call legeders(n+1,y,polsy,poldersy)

vals  = 0
dersx = 0
dersy = 0
idx  = 1
do nn=0,n
do jj=0,nn
ii = nn-jj
vals  = vals  + coefs(idx,:)*polsx(ii)*polsy(jj)
dersx = dersx + coefs(idx,:)*poldersx(ii)*polsy(jj)
dersy = dersy + coefs(idx,:)*polsx(ii)*poldersy(jj)
idx   = idx+1
end do
end do

end subroutine


subroutine bilege_evalder4(n,coefs,x,y,vals,dersx,dersy)
implicit double precision (a-h,o-z)
double complex                  :: xslege(n), whtslege(n)
double complex                  :: coefs(:,:),vals(:),dersx(:),dersy(:)
!
!  Given the coefficients of a collection of complex-valued expansions of the
!  form (1), evaluate them and their derivatives  at a specified point.
!
!  Input parameters:
!    n - the order of the bivariate expansion
!    coefs - the vector (2) of expansion coefficients
!    (x,y) - the point at which to evaluate the expansion (1)
!
!  Output parameters:
!   val - the value of the expansion
!

double precision          :: polsx(0:n), polsy(0:n)
double precision          :: poldersx(0:n), poldersy(0:n)

call legeders(n+1,x,polsx,poldersx)
call legeders(n+1,y,polsy,poldersy)

vals  = 0
dersx = 0
dersy = 0
idx  = 1
do nn=0,n
do jj=0,nn
ii = nn-jj
vals  = vals  + coefs(idx,:)*polsx(ii)*polsy(jj)
dersx = dersx + coefs(idx,:)*poldersx(ii)*polsy(jj)
dersy = dersy + coefs(idx,:)*polsx(ii)*poldersy(jj)
idx   = idx+1
end do
end do

end subroutine

end module
