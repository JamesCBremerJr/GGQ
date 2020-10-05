module     test_linalg_funs
use utils
contains
end module test_linalg_funs

program test_linalg

use utils
use linalg

implicit double precision (a-h,o-z)

double precision, allocatable       :: a(:,:),q(:,:),r(:,:),xx(:,:),yy(:,:),b(:,:)
double precision, allocatable       :: u(:,:),vt(:,:),sigma(:), rnorms(:)
double precision, allocatable       :: x(:), y(:), x0(:), uu(:), vv(:)

integer, allocatable                :: ipivs(:)
integer                             :: iseed(4)

eps = 1.0d-14

!
!  Test the factorization routines
!

n   = 20000
m   = 1000
k   = 77

allocate(a(n,m) )
allocate(xx(n,k), yy(k,m))

call id_srand(n*k,xx)
call id_srand(m*k,yy)
xx = xx / norm2(xx)
yy = yy / norm2(yy)

do j=1,k
xx(:,j) = xx(:,j) * 1.5d0**(-j)
end do
a = matmul(xx,yy)
deallocate(xx,yy)

call prini("n = ",n)
call prini("m = ",m)

call elapsed(t1)
call intdecomp(eps,n,m,a,krank,ipivs,r)
call elapsed(t2)
call prini("krank = ",krank)
call prin2("intdecomp time = ",t2-t1)
derr = norm2(matmul(a(:,ipivs),r)-a) / norm2(a)
call prin2("error in intdecomp = ",derr)
call prina("")

call elapsed(t1)
call qrdecomp(eps,n,m,a,krank,ipivs,q,r)
call elapsed(t2)
allocate(rnorms(krank))
do i=1,krank
rnorms(i) = abs(r(i,i))
end do



call prini("krank = ",krank)
call prin2("qrdecomp time = ",t2-t1)
derr = norm2(matmul(q,r)-a(:,ipivs)) / norm2(a)
call prin2("error in qrdecomp = ",derr)
call prina("")

call elapsed(t1)
call svd(eps,n,m,a,krank,u,sigma,vt)
call elapsed(t2)
call prini("krank = ",krank)
call prin2("svd time = ",t2-t1)
do i=1,krank
u(:,i) = u(:,i) * sigma(i)
end do
derr = norm2(matmul(u,vt)-a) / norm2(a)
call prin2("error in svd = ",derr)
call prina("")


! call prin2("sings = ",sigma)
! call prin2("rnorms = ",rnorms)

!
!  Test the linear solve routines
!
deallocate(a)

n  = 1000
allocate(a(n,n),x(n),y(n),x0(n),b(n,n))
a  = 0
x0 = 0

do j=1,n
do i=1,j
a(i,j) = 1
end do
end do
b = a

call id_srand(n,x0)
y = matmul(a,x0)

call elapsed(t1)
call linsolve(n,a,y,x)
call elapsed(t2)
a = b

derr  = norm2(x-x0)/norm2(x0)
call prini("n = ",n)
call prin2("linsolve time = ",t2-t1)
call prin2("linsolve error = ",derr)
call prina("")


!
!  Test inversion and the SMW formula
!
call elapsed(t1)
call invert(n,a)
call elapsed(t2)
call prini("n = ",n)

x = matmul(a,y)
derr  = norm2(x-x0)/norm2(x0)

call prin2("invert time = ",t2-t1)
call prin2("invert error = ",derr)


allocate(uu(n), vv(n))
call id_srand(n,uu)
call id_srand(n,vv)

uu = uu / norm2(uu)
vv = vv / norm2(vv)

do i=1,n
do j=1,n
b(i,j) = b(i,j) + uu(i)*vv(j)
end do
end do

call invert(n,b)
call elapsed(t1)
call smw_update(n,a,uu,vv)
call elapsed(t2)
derr = norm2(a-b)
call prin2("SMW time = ",t2-t1)
call prin2("SMW error = ",derr)
call prina("")

!
!  Least squares
!
n = 2000
m = 1000

deallocate(a,x,y,x0)

allocate(a(n,m),x(m),y(n),x0(m))

call id_srand(m,x0)
call id_srand(n*m,a)

y = matmul(a,x0)

call elapsed(t1)
call leastsq(n,m,a,y,x)
call elapsed(t2)

derr  = norm2(x-x0)/norm2(x0)
call prini("n = ",n)
call prini("m = ",m)
call prin2("leastsq time = ",t2-t1)
call prin2("leastsq error = ",derr)

end program
