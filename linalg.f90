!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for performing certain linear algrebraic operations.  Most
!  of the subroutines are wrappers around LAPACK and the ID library of Tygert, et. al.
!
!  The following subroutines are publicly callable:
!
!    linsolve - solve the system of equations Ax = y with A is a square matrix
!
!    leastsq - solve the system Ax = y, where A is rectangular, in a least
!      square sense
!
!    invert - construct the inverse of a square matrix
!
!    smw_update - use the Sherman-Morrison-Woodbury formula to compute the inverse of
!       (A + UV^T) given the inverse of A
!
!    intdecomp - construct a rank-revealing column interpolative decomposition of 
!      an input matrix
!
!    qrdecomp - construct a rank-revealing QR decomposition of an input matrix
!
!    svd - construct a rank-revealing SVD decomposition of an input matrix
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module     linalg
use utils

contains

subroutine linsolve(n,a,y,x)
implicit double precision (a-h,o-z)
double precision       :: a(:,:), x(:), y(:)
!
!  Solve the equation Ax = y, where A is a square matrix.
!
!  **THE INPUT MATRIX IS DESTROYED BY THIS ROUTINE**
!
!  Input parameters:
!    n - the dimensions of the input matrix
!    a - the (n,n) input matrix, WHICH IS DESTROYED BY THIS ROUTINE
!    y - the right-hand side of the equation
!
!  Output parameters:
!    x - the solution of the equation
!
double precision, allocatable :: q(:,:), r(:,:)
integer, allocatable :: ipivs(:)

eps0 = epsilon(0.0d0)
if (eps0 .lt. 1.0d-17) then
x = y
call linalg_qrsolve(a,n,x)
return
endif

!
!  Use LAPACK
!

x = y
allocate(ipivs(n))
call dgesv(n,1,a,n,ipivs,x,n,info)

end subroutine


subroutine leastsq(n,m,a,y,x)
implicit double precision (a-h,o-z)
double precision           :: a(:,:), x(:), y(:)

!
!  Solve the equation Ax = y, where A is a rectangular matrix, in a least 
!  squares sense.  The input matrix can be overdetermined or underdetermined.
!  In the overdetermined case, the solution with the smallest L^2 norm is
!  returned.
!
!  **THE INPUT MATRIX IS DESTROYED BY THIS ROUTINE**
!
!  Input parameters:
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    y - the right-hand side of the system
!
!  Output parameters:
!    x - the least squares solution
!
!

double precision, allocatable :: work(:)
integer, allocatable          :: ipivs(:)

double precision, allocatable :: q1(:,:), r1(:,:), b(:,:)
double precision, allocatable :: q2(:,:), r2(:,:), z(:), r(:,:)
double precision, allocatable :: u(:)

eps0 = epsilon(0.0d0)

! factor the input matrix as a = U R2 V^* with U, V orthogonal
! using the QR decomposition and use the factorization to solve
! the problem
if (eps0 .lt. 1.0d-17) then
call qrdecomp0(eps0,n,m,a,krank,ipivs,q1,r1)
allocate(b(m,krank))
b = transpose(r1)
call qrdecomp0_np(m,krank,b,q2,r2)
deallocate(b)
allocate(z(krank),r(krank,krank),u(m))
z  = matmul(transpose(q1),y)
r  = transpose(r2)
call  linalg_qrsolve(r,krank,z)
u = matmul(q2,z)
x(ipivs) = u
return
endif

k = max(n,m)
allocate(ipivs(m),z(k))
lwork  = -1 
z(1:n) = y
call dgelsy(n,m,1,a,n,z,k,ipivs,rcond,krank,ww,lwork,info) 
lwork = ww
allocate(work(lwork))
call dgelsy(n,m,1,a,n,z,k,ipivs,rcond,krank,work,lwork,info) 
x(1:m) = z

end subroutine


subroutine invert(n,a)
implicit double precision (a-h,o-z)
integer                                        :: n
double precision                               :: a(n,n)
!
!  Invert a user-supplied matrix.
!
!  Input parameters:
!    n - the dimensions of the matrix to invert
!    a - the (n,n) matrix to invert which will replaced with its inverse
!
!  Output parameters:
!    a - the inverse of A
!
!

double precision, allocatable                  :: work(:), u(:,:), sigma(:), vt(:,:), b(:,:)
integer, allocatable                           :: ipiv(:)

eps0 = epsilon(0.0d0)
if (eps0 .lt. 1.0d-17) then
call linalg_orthom(a,n)
return
endif


allocate(ipiv(n))
call dgetrf(n,n,a,n,ipiv,info)
lwork = -1
call dgetri(n,a,n,ipiv,ww,lwork,info)
lwork = ww
allocate(work(lwork))
call dgetri(n,a,n,ipiv,work,lwork,info)

end subroutine


subroutine smw_update(n,ainv,u,v)
implicit double precision (a-h,o-z)
double precision                           :: ainv(:,:)
double precision                           :: u(:)
double precision                           :: v(:)
!
!  Use the Sherman-Morrison-Woodbury formula to compute
!
!          (A + uv^t)^(-1),
!
!  where u and v are vectors, given the inverse of A.
!
!  Input parameters:
!    n - the dimension of the matrix a
!    ainv - the (n,n) inverse of a
!    u, v - a pair of vectors
!
!  Output parameters:
!    ainv - the inverse of the updated matrix
!

double precision, allocatable :: z(:), w(:)


allocate(z(n),w(n))

eps0  = epsilon(0.0d0)
z     = matmul(ainv,u)
denom = 1 + dot_product(v,z)
w     = matmul(v,ainv)/denom

if (eps0 .lt. 1.0d-17) then
do i=1,n
do j=1,n
ainv(i,j) = ainv(i,j) -z(i)*w(j)
end do
end do
else
call dger(n,n,-1.0d0,z,1,w,1,ainv,n)
endif

end subroutine


subroutine intdecomp(eps,n,m,a,krank,ipivs,r,krankest)
implicit double precision (a-h,o-z)
integer                                    :: n,m
double precision                           :: a(:,:)
double precision, allocatable, intent(out) :: r(:,:)
integer, allocatable, intent(out)          :: ipivs(:)
integer, optional                          :: krankest
!
!  Construct a rank-revaling column interpolation decomposition of an input matrix.
!  That is, factor the input matrix A as
!
!    A = A(:,ipivs) * R
!  
!  where ipivs is an integer array of length krank and R is a (krank,m)
!  matrix.  The rank of the decomposition (krank) is chosen so that 
!  the decomposition roughly achieves a level of relative accuracy in the L^2 norm
!  specified by the user.
!
!  Input parameters:
!    eps - the desired precision for the factorization
!    (n,m) - the dimensions of the input matrix a
!    a - the input matrix, which is not destroyed by this routine
!    krankest - an optional paramter specifying an estimate of the rank
!      of the matrix --- providing a *GOOD* estimate can greatly
!      accelerate the algorithm
!
!  Output parameters:
!    krank - the rank of the factorization 
!    ipivs - an array of length krank specifying the columns of A which
!      form a basis for the column space of A
!    r - the (n,krank-m) matrix appearing in the factorization
!

double precision, allocatable :: b(:,:)
double precision, allocatable :: wrand(:), wrand2(:)


if ( PRESENT(krankest) ) then
krankmax = krankest + 20
else
krankmax = 100
endif


! use the nonrandomized version of the routine if krankmax is too large
if (krankmax .gt. n-20) then
call intdecomp0(eps,n,m,a,krank,ipivs,r)
return
endif

!
!  Apply a random transform to a so as to form b
!

allocate( wrand(101*n+1000) )
call idd_frmi(n,n2,wrand)
allocate(b(n2,m))


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(j,wrand2)
allocate(wrand2(101*n+1000))
wrand2 = wrand
!$OMP DO
do j=1,m
call idd_frm(n,n2,wrand2,a(:,j),b(:,j))
end do
!$OMP END DO
!$OMP END PARALLEL



!
!  Compute a QR decomposition of the reduced matrix with
!  the number of rows increasing ...
!

ifdone = 0

do while (ifdone .eq. 0)

! resort to the nonrandomized routine if the rank is too large
if (krankmax .gt. n2-20 .OR. krankmax .gt. n-20) then
call intdecomp0(eps,n,m,a,krank,ipivs,r)
return
endif

call intdecomp0(eps,krankmax,m,b(1:krankmax,:),krank,ipivs,r)
if (krank .gt. krankmax-20) then
krankmax = 2*krankmax
else
ifdone = 1
endif
end do

end subroutine


subroutine intdecomp0(eps,n,m,a,krank,ipivs,r)
implicit double precision (a-h,o-z)
double precision                           :: a(:,:)
double precision, allocatable, intent(out) :: r(:,:)
integer, allocatable, intent(out)          :: ipivs(:)
!
!  Constructthe interpolative decomposition of a user-supplied matrix
!  WITHOUT randomization.
!
double precision, allocatable              :: b(:,:), rnorms(:)
double precision, allocatable              :: tau(:), work(:), rr(:,:)
integer, allocatable                       :: jpivs(:)

eps0 = epsilon(0.0d0)


allocate(b(n,m))
b = a

!
!  Use the ID library when the code is run in extended precision arithmetic.
!
if (eps0 .lt. 1.0d-17) then
allocate(jpivs(m+n),rnorms(m+n))
call iddp_id(eps,n,m,b,krank,jpivs,rnorms)
allocate(rr(krank,m-krank))
call linalg_dcopy(krank*(m-krank),b,1,rr,1)

else

!
!  Otherwise, use LAPACK.
!

k     = min(n,m)
allocate(tau(k), jpivs(m))

jpivs = 0
lwork = -1
call dgeqp3(n,m,b,n,jpivs,tau,ww,lwork,info)

lwork = ww
allocate(work(lwork))
call dgeqp3(n,m,b,n,jpivs,tau,work,lwork,info)

dd    = eps * (eps+abs(b(1,1)))
do krank=0,k-1
if (abs(b(krank+1,krank+1)) .lt. dd)  exit
end do

allocate(rr(krank,m-krank))
rr = b(1:krank,krank+1:m)
call dtrtrs('U','N','N',krank,m-krank,b,n,rr,krank,info)

endif

allocate(r(krank,m), ipivs(krank))

do i=1,krank
r(:,jpivs(i)) = 0
r(i,jpivs(i)) = 1
end do

r(:,jpivs(krank+1:m)) = rr
ipivs = jpivs(1:krank)

end subroutine


subroutine qrdecomp(eps,n,m,a,krank,ipivs,q,r,krankest)
implicit double precision (a-h,o-z)
integer                                    :: n,m
double precision                           :: a(:,:)
double precision, allocatable, intent(out) :: q(:,:), r(:,:)
integer, allocatable, intent(out)          :: ipivs(:)
integer, optional                          :: krankest
!
!  Construct a rank-revealing QR decomposition of a real-valued input matrix A.  That is,
!  factor A as 
!
!    AP = QR = Q [R1 R2] 
!
!  where P is an (m,m)  permutation matrix, Q is an (m,krank) matrix whose columns are 
!  orthonormal, and R1 and an (krank,krank) upper triangular matrix.  
!
!  The rank of the decomposition (krank) is chosen so that the decomposition
!  roughly achieves a specified relative accuracy in the L^2 norm.
!
!  Input parameters:
!    eps - desired precision for the factorization
!    (n,m) - the dimensions of the input matrix a
!    a - the input matrix, which will NOT be destroyed by this routine
!    krankest - an OPTIONAL argument which gives an estimate for the rank
!      of the matrix
!
!  Output parameters:
!    krank - the rank of the decomposition constructed
!    ipivs - an integer array specifying the columns of the 
!    q - the matrix q in the decomposition
!    r - the matrix r in the decomposition
!

double precision, allocatable :: b(:,:), r1(:,:), r3(:,:), r4(:,:)
double precision, allocatable :: a2(:,:), r2(:,:)

!
!  Compute the interpolative decomposition of the input matrix
!


call intdecomp(eps,n,m,a,krank,ipivs,r1,krankest)
allocate(r3(krank,m))
r3 = r1(:,ipivs)

!
!  Construct the QR decomposition of the submatrix of A
!
allocate(a2(n,krank))
a2 = a(:,ipivs)

call qrdecomp0_np(n,krank,a2,q,r2)
allocate(r(krank,m))
r = matmul(r2,r3)

end subroutine


subroutine qrdecomp0_np(n,m,a,q,r)
implicit double precision (a-h,o-z)
double precision                           :: a(:,:)
double precision, allocatable, intent(out) :: q(:,:)
double precision, allocatable, intent(out) :: r(:,:)
!
!  Use a nonrandomized algorithm to construct the QR decomposition of a user-supplied
!  matrix via nonpivoted algorithm IN THE CASE m <=n.
!
double precision, allocatable :: b(:,:), tau(:)
double precision, allocatable :: work(:)
integer, allocatable          :: ipivs(:)

eps0 = epsilon(0.0d0)

! use GS with reorthogonalization if the code is being run in extended precision
if (eps0 .lt. 1.0d-17) then
allocate(q(n,m), r(m,m) )
q = a
call linalg_gsnopiv(n,m,q)
r = matmul(transpose(q),a)
return
endif

allocate(b(n,m), ipivs(m), tau(m) )
b    = a

do i=1,m
ipivs(i)=i
end do

lwork = -1
call dgeqp3(n,m,b,n,ipivs,tau,ww,lwork,info)

lwork = ww
allocate(work(lwork))
call dgeqp3(n,m,b,n,ipivs,tau,work,lwork,info)

allocate(q(n,m),r(m,m))

r = 0
do j=1,m
r(1:j,j) = b(1:j,j)
end do

call dorgqr(n, m, m, b, n, tau, work, lwork, info)
q = b

end subroutine

subroutine qrdecomp0(eps,n,m,a,krank,ipivs,q,r)
implicit double precision (a-h,o-z)
double precision                           :: a(:,:)
integer, allocatable, intent(out)          :: ipivs(:)
double precision, allocatable, intent(out) :: q(:,:)
double precision, allocatable, intent(out) :: r(:,:)
!
!  Use a nonrandomized algorithm to construct the QR decomposition of
!  a user-supplied matrix.  This version of the routine uses a
!  pivoted algorithm.
!

double precision, allocatable :: b(:,:), rnorms(:), tau(:)
double precision, allocatable :: work(:), r2(:,:)


allocate(b(n,m))
b = a

eps0 = epsilon(0.0d0)


if (eps0 .lt. 1.0d-17) then
allocate(ipivs(m),rnorms(m))
call linalg_gspiv(b,n,m,eps,rnorms,ipivs,krank)
allocate(q(n,krank), r(krank,m), r2(krank,m))
q  = b(:,1:krank)
r2 = matmul(transpose(q),a)
r  = r2(:,ipivs)
do i=1,krank
end do
return
endif

k = min(n,m)

allocate(ipivs(m), tau(k) )

ipivs = 0 
lwork = -1
call dgeqp3(n,m,b,n,ipivs,tau,ww,lwork,info)

lwork = ww
allocate(work(lwork))
call dgeqp3(n,m,b,n,ipivs,tau,work,lwork,info)

!
!  Find the rank and extract Q and R
!

dd = abs(b(1,1))
krank = 0
do krank=0,k-1
if (abs(b(krank+1,krank+1)) .lt. eps*dd) exit
end do

allocate(q(n,krank),r(krank,m))

r = 0
do j=1,krank
r(1:j,j) = b(1:j,j)
end do

do j=krank+1,m
r(1:krank,j) = b(1:krank,j)
end do

call dorgqr(n, krank, krank, b, n, tau, work, lwork, info)
q = b(1:n,1:krank)

end subroutine


subroutine svd(eps,n,m,a,krank,u,sigma,vt,krankest)
implicit double precision (a-h,o-z)
double precision                           :: a(:,:)
double precision, allocatable, intent(out) :: u(:,:)
double precision, allocatable, intent(out) :: vt(:,:)
double precision, allocatable, intent(out) :: sigma(:)
integer, optional                          :: krankest
!
!  Compute a rank-revealing singular value decomposition.  That is,
!  factor a as
!
!    a(n,m) = u(n,krank) sigma(krank,krank) vt(krank,m)
!
!  where the columns of u and the rows of v are orthonormal,
!  and sigma is the diagonal matrix whose nonzero entries are the
!  krank largest singular values of a, in descending order.
!
!  The rank of the decomposition, krank, is chosen so that the
!  factorization is roughly accurate to a specified precision
!  in the relative L^2 norm.
!
!
!  Input parameters:
!    eps - the precision for the factorization
!    (n,m) - the dimensions of the input matrix
!    a - the input matrix
!    krankest - an optional paramter specifying an estimate of the rank
!      of the matrix --- providing a *GOOD* estimate can greatly
!      accelerate the algorithm
!
!  Output parameters:
!    krank - the apparent rank of a to the specified precision
!    u - the (n,krank) matrix appearing in the factorization
!    sigma - an array of length krank giving the largest krank
!     singular values of a --- in other words, the array of
!     diagonal entries of the matrix sigma
!    vt - the (krank,m) matrix appearing the factorization
!

double precision, allocatable :: b(:,:), r(:,:)

integer, allocatable          :: jpivs(:)
double precision, allocatable :: q1(:,:), r1(:,:)
double precision, allocatable :: q2(:,:), r2(:,:)
double precision, allocatable :: u0(:,:), vt0(:,:)

!
!  Compute the interpolative decomposition of the input matrix
!

call intdecomp(eps,n,m,a,krank,jpivs,r,krankest)

!
!  Factor the neutered matrix as q1 * r1
!
allocate(b(n,krank))
b = a(:,jpivs)
call qrdecomp0_np(n,krank,b,q1,r1)
deallocate(b)

!
!  Factor the transpose of R as q2 * r2
!
allocate(b(m,krank))
b = transpose(r)
call qrdecomp0_np(m,krank,b,q2,r2)
deallocate(b)

!
!  Compute the SVD of R1 * R2^T
!
allocate(b(krank,krank))
b = matmul(r1,transpose(r2))
call svd00(krank,krank,b,u0,sigma,vt0)
deallocate(b)


!
!  Form u and vt 
!

allocate(u(n,krank), vt(krank,m) )

u  = matmul(q1,u0)
vt = matmul(vt0,transpose(q2))

end subroutine

subroutine svd0(eps,n,m,a,krank,u,sigma,vt)
implicit double precision (a-h,o-z)
double precision                           :: a(:,:)
double precision, allocatable, intent(out) :: u(:,:)
double precision, allocatable, intent(out) :: vt(:,:)
double precision, allocatable, intent(out) :: sigma(:)
!
!  Compute a rank-revealing singular value decomposition using a nonrandomized
!  algorithm.
!
double precision, allocatable :: sings(:), u0(:,:), vt0(:,:)
double precision, allocatable :: work(:), b(:,:)

eps0 = epsilon(0.0d0)

if (eps0 .lt. 1.0d-17) then
call prina("svd0: extended precision version of the SVD not yet implemented")
stop
endif

allocate(b(n,m))
b = a

k = min(n,m)
allocate(sings(k),u0(n,k),vt0(k,m))

lwork = -1
call dgesvd('S','S',n,m,b,n,sings,u0,n,vt0,m,ww,lwork,info)
lwork = ww
allocate(work(lwork))
call dgesvd('S','S',n,m,b,n,sings,u0,n,vt0,k,work,lwork,info)

dd = eps*sings(1)
do krank=0,k-1
if (sings(krank+1) .lt. dd) exit
end do

allocate(u(n,krank), vt(krank,m), sigma(krank) )
u     = u0(:,1:krank)
vt    = vt0(1:krank,:)
sigma = sings(1:krank)

end subroutine



subroutine svd00(n,m,a,u,sigma,vt)
implicit double precision (a-h,o-z)
double precision                           :: a(:,:)
double precision, allocatable, intent(out) :: u(:,:)
double precision, allocatable, intent(out) :: vt(:,:)
double precision, allocatable, intent(out) :: sigma(:)
!
!  Compute a singular value decomposition using a nonrandomized
!  algorithm.
!
double precision, allocatable :: sings(:), u0(:,:), vt0(:,:)
double precision, allocatable :: work(:), b(:,:)

eps0 = epsilon(0.0d0)

if (eps0 .lt. 1.0d-17) then
call prina("svd0: extended precision version of the SVD not yet implemented")
stop
endif

allocate(b(n,m))
b = a

k = min(n,m)
allocate(sigma(k),u(n,k),vt(k,m))

lwork = -1
call dgesvd('S','S',n,m,b,n,sigma,u,n,vt,k,ww,lwork,info)
lwork = ww
allocate(work(lwork))
call dgesvd('S','S',n,m,b,n,sigma,u,n,vt,k,work,lwork,info)

end subroutine



subroutine singvals00(eps,n,m,a,sigma)
implicit double precision (a-h,o-z)
double precision                           :: a(:,:)
double precision, allocatable, intent(out) :: sigma(:)
!
!  Compute the singular values of a matrix.
!
double precision, allocatable :: sings(:)
double precision, allocatable :: work(:), b(:,:)

eps0 = epsilon(0.0d0)

if (eps0 .lt. 1.0d-17) then
call prina("svd0: extended precision version of the SVD not yet implemented")
stop
endif

allocate(b(n,m))
b = a

k = min(n,m)
allocate(sigma(k))

lwork = -1
call dgesvd('N','N',n,m,b,n,sigma,u0,n,vt0,k,ww,lwork,info)
lwork = ww
allocate(work(lwork))
call dgesvd('N','N',n,m,b,n,sigma,u0,n,vt0,k,work,lwork,info)

end subroutine

subroutine linalg_gspiv(b,n,m,eps,rnorms,ipivots,ncols)
implicit double precision (a-h,o-z)
dimension rnorms(m),ipivots(m)
double precision b(n,m),cd
!
!  Perform Gram-Schmidt with pivoting and reorthogonalization on the input
!  matrix.
!
double precision, allocatable :: v(:)
allocate(v(n))

!      
! Set rank to zero just in case.
!     
ncols=0

!
! Initialize the array of pivots
! 
do i=1,m
ipivots(i)=i
end do

! 
! Prepare the array of values of norms of columns
! 
dtot=0
do i=1,m
d         = dot_product(b(:,i),b(:,i))
dtot      = dtot+d
rnorms(i) = sqrt(d)
end do

if (dtot .lt. eps**2) then
ncols = 0
return
endif

thresh=dtot*eps**2
thresh=sqrt(thresh)

! 
!  Conduct GS
! 
do i=1,min(m,n)
ipivot=i
rn=rnorms(i)
do j=i+1,m
if(rnorms(j) .le. rn) cycle
rn=rnorms(j)
ipivot=j
end do


! 
! Put the column number ipivot in the i-th place
!
v            = b(:,i)
b(:,i)       = b(:,ipivot)
b(:,ipivot)  = v
 
ii              = ipivots(i)
ipivots(i)      = ipivots(ipivot)
ipivots(ipivot) = ii
 
d=rnorms(ipivot)
rnorms(ipivot)=rnorms(i)
rnorms(i)=d
! 
! orthogonalize the i-th column to all preceding ones
! 
do j=1,i-1
cd = dot_product(b(:,i), b(:,j))
b(:,i) = b(:,i) - cd*b(:,j)
end do
! 
! normalize the i-th column
!
cd = dot_product(b(:,i),b(:,i))
d  = sqrt(cd)
! 
if(d .lt. thresh ) return
 
ncols=i
b(:,i) = b(:,i) / d
 
if(i .eq. m) return

! 
! orthogonalize everything else to it
! 
do j=i+1,m
if(rnorms(j) .lt. thresh) cycle
cd = dot_product(b(:,i), b(:,j))
b(:,j) = b(:,j) - cd*b(:,i)
dd     = norm2(b(:,j))
rnorms(j)=sqrt(dd)
end do
end do


return
end subroutine


subroutine linalg_gsnopiv(n,m,b)
implicit double precision (a-h,o-z)
double precision            :: b(:,:)
!
!  Orthonormalize the columns of b via the Gram-Schmidt algorithm with
!  reorthogonalization.  NOTE THAT THIS IS A NONPIVOTED VERSION OF
!  THIS ROUTINE.
!

l = min(n,m)

do i=1,l
       
do j=1,i-1
cd = dot_product(b(:,i),b(:,j))
b(:,i) = b(:,i) - b(:,j)*cd
end do
        
d = dot_product(b(:,i),b(:,i))
b(:,i) = b(:,i) /sqrt(d)

do j=i+1,m
cd = dot_product(b(:,i),b(:,j))
b(:,j)    = b(:,j) - b(:,i)*cd
end do
end do
end subroutine



subroutine linalg_gspiv_partial(b,n,m,i1,i2,eps,rnorms,ipivots,ncols)
implicit double precision (a-h,o-z)
dimension rnorms(m),ipivots(m)
double precision b(n,m),cd
!
!  Orthonormalize a subblock b(i1:i2,:) of the input matrix, while performing
!  the same operations on the entire column.
!
double precision, allocatable :: v(:)
allocate(v(n))

!      
! Set rank to zero just in case.
!     
ncols=0

!
! Initialize the array of pivots
! 
do i=1,m
ipivots(i)=i
end do

! 
! Prepare the array of values of norms of columns
! 
dtot=0
do i=1,m
d         = dot_product(b(i1:i2,i),b(i1:i2,i))
dtot      = dtot+d
rnorms(i) = sqrt(d)
end do

if (dtot .lt. eps**2) then
ncols = 0
return
endif

thresh=dtot*eps**2
thresh=sqrt(thresh)

! 
!  Conduct GS
! 
do i=1,min(m,n)
ipivot=i
rn=rnorms(i)
do j=i+1,m
if(rnorms(j) .le. rn) cycle
rn=rnorms(j)
ipivot=j
end do


! 
! Put the column number ipivot in the i-th place
!
v            = b(:,i)
b(:,i)       = b(:,ipivot)
b(:,ipivot)  = v
 
ii              = ipivots(i)
ipivots(i)      = ipivots(ipivot)
ipivots(ipivot) = ii
 
d=rnorms(ipivot)
rnorms(ipivot)=rnorms(i)
rnorms(i)=d
! 
! orthogonalize the i-th column to all preceding ones
! 
do j=1,i-1
cd = dot_product(b(i1:i2,i), b(i1:i2,j))
b(:,i) = b(:,i) - cd*b(:,j)
end do
! 
! normalize the i-th column
!
cd = dot_product(b(i1:i2,i),b(i1:i2,i))
d  = sqrt(cd)
! 
if(d .lt. thresh ) return
 
ncols=i
b(:,i) = b(:,i) / d
 
if(i .eq. m) return

! 
! orthogonalize everything else to it
! 
do j=i+1,m
if(rnorms(j) .lt. thresh) cycle
cd = dot_product(b(i1:i2,i), b(i1:i2,j))
b(:,j) = b(:,j) - cd*b(:,i)
dd     = norm2(b(i1:i2,j))
rnorms(j)=sqrt(dd)
end do
end do

return
end subroutine

subroutine linalg_dcopy(n,dx,incx,dy,incy)
implicit double precision (a-h,o-z)
integer          :: n,incb,incc
double precision :: dx(1),dy(1)

ix = 1
iy = 1
do i=1,n
dy(iy) = dx(ix)
ix = ix + incx
iy = iy + incy
end do

end subroutine


subroutine linalg_qrsolve(a,n,rhs)
implicit double precision (a-h,o-z)
integer            :: n
double precision   :: a(n,n),rhs(n)
!
!  This subroutine uses a QR-decomposition to solve the equation
!  A x = b.  Both the input matrix a and the right-hand side are destroyed
!  by the routine.
!
!  Input parameters:
!    a - the (n,n) matrix of coefficients
!    n - an integer specifying the size of the system of 
!    rhs - a vector of length n speciying the rhs of the system
!
!  Output parameters:
!   rhs - upon return, the solution of the linear system


double precision :: aa(2),u(2,2)

! 
! transpose the input matrix a 
!

size22=0
do i=1,n
do j=1,i
d=a(j,i)
a(j,i)=a(i,j)
a(i,j)=d
size22=size22+a(j,i)**2
size22=size22+a(i,j)**2
end do
end do

!
!  Reduce to upper triangular 
!
do i=1,n-1
do j=n,i+1,-1
aa(1)=a(i,j-1)
aa(2)=a(i,j)

u22=-aa(1)
u12=aa(2)
d=u22**2+u12**2
if(d .lt. size22*1.0d-66) then
u(2,2)=1
u(1,2)=0
u(1,1)=1
u(2,1)=0
else
d=sqrt(d)
u(2,2)=u22/d
u(1,2)=u12/d
u(1,1)=-u(2,2)
u(2,1)=u(1,2)
endif

do ii=i,n
d1=u(1,1)*a(ii,j-1)+u(1,2)*a(ii,j)
d2=u(2,1)*a(ii,j-1)+u(2,2)*a(ii,j) 
a(ii,j-1)=d1
a(ii,j)=d2
end do

d1=u(1,1)*rhs(j-1)+u(1,2)*rhs(j)
d2=u(2,1)*rhs(j-1)+u(2,2)*rhs(j)
rhs(j-1)=d1
rhs(j)=d2
end do
end do

!
!  Apply the inverse of the triangular matrix
! 

rhs(n)=rhs(n)/a(n,n)
do i=n-1,1,-1
d=0
do j=n,i+1,-1
d=d+a(j,i)*rhs(j)
end do
rhs(i)=(rhs(i)-d)/a(i,i)
end do

return
end subroutine


subroutine linalg_orthom(a,n)
implicit double precision (a-h,o-z)
double precision a(n,n)
double precision, allocatable :: b(:,:), work(:)
allocate(b(n,n), work(n))
done=1
zero=0

b = 0.0d0
do i=1,n
b(i,i) = 1.0d0
end do
        
! GS the input matrix
do i=1,n
do j=1,i-1
cd = dot_product(a(i,:),a(j,:))
a(i,:) = a(i,:) - a(j,:)*cd
b(i,:) = b(i,:) - b(j,:)*cd
end do

d = 1/norm2(a(i,:))
a(i,:) = a(i,:)*d
b(i,:) = b(i,:)*d
do j=i+1,n
cd = dot_product(a(i,:),a(j,:))
a(j,:) = a(j,:)-cd*a(i,:)
b(j,:) = b(j,:)-cd*b(i,:)
end do
end do
!
!  Multiply the adjoint of the resulting orthogonal matrix by the
!  triangular one to inverse a  
!
do i=1,n
do j=1,n
cd = dot_product(a(:,i),b(:,j))
work(j)=cd
end do
a(:,i) = work
end do
!
!  Transpose A
!
b = transpose(a)
a = b

return
end subroutine


function eye(n) result(a)
implicit double precision (a-h,o-z)
double precision, allocatable :: a(:,:)
integer n

!
!  Return an identity matrix which is dimensioned (n,n).
!

allocate(a(n,n))
a = 0
do i=1,n
a(i,i)=1.0d0
end do

end function


end module
