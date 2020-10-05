program test_bilege
use utils
use legendre
use bilege

implicit double precision (a-h,o-z)

double precision, allocatable    :: xs(:),ys(:),whts(:)

double precision, allocatable    :: xslege(:),whtslege(:)
double precision, allocatable    :: xsout(:),ysout(:),whtsout(:)


double precision, allocatable    :: errsint(:,:), umatr(:,:)
double precision, allocatable    :: amatr(:,:)
double complex, allocatable      :: vals(:),coefs(:)
double complex, allocatable      :: valsout(:),valsout0(:)
double complex                   :: val, val0, ima, derx, dery, derx0, dery0

double precision, allocatable    :: polsx(:),polsy(:),sings(:),work(:)


ima = (0.0d0,1.0d0)

n = 20
call bilege_quad(n,nquad,xs,ys,whts)
call prini("n     = ",n)
call prini("nquad = ",nquad)

if (n==4)  m = 7
if (n==6)  m = 11
if (n==8)  m = 14
if (n==10) m = 18
if (n==12) m = 20
if (n==14) m = 23      
if (n==16) m = 27
if (n==18) m = 30
if (n==20) m = 34
if (n==24) m = 41
if (n==30) m = 50      

!
!  Check that the quadrature integrates the appropriate set of
!  polynomials
!

allocate(errsint(0:m,0:m))

do nn=0,m
do jj=0,nn
ii = nn-jj
val          = sum(xs**ii*ys**jj*whts)
val0         = (1.0d0 + (-1.0d0)**ii) / (ii+1.0d0) *  (1.0d0 + (-1.0d0)**jj) / (jj+1.0d0)
errsint(ii,jj) = abs(val-val0)
end do
end do
call prin2("max poly integration error = ",maxval(errsint))


!
!  Test the coefficient expansions
!

allocate(vals(nquad),coefs(nquad))
vals = xs*ys*sqrt(whts)

call bilege_coefs(n,vals,coefs)

nn = 10
do i=1,nn
do j=1,nn
x = -1.0d0 + 2*(i-1.0d0)/(nn-1.0d0)
y = -1.0d0 + 2*(j-1.0d0)/(nn-1.0d0)
call bilege_eval(n,coefs,x,y,val)
val0 = x*y
derr = abs(val-val0)
dmax = max(derr,dmax)

end do
end do

call prin2("maximum evaluation error = ",dmax)

!
!  Test the values to coefficient matrix
!

call bilege_coefsmatrix(n,umatr)
coefs = matmul(umatr,vals)

nn = 10
do i=1,nn
do j=1,nn
x = -1.0d0 + 2*(i-1.0d0)/(nn-1.0d0)
y = -1.0d0 + 2*(j-1.0d0)/(nn-1.0d0)

call bilege_evalder(n,coefs,x,y,val,derx,dery)
val0  = x*y
derx0 = y
dery0 = x

derr = abs(val-val0) + abs(derx0-derx) + abs(dery0-dery)
dmax = max(derr,dmax)

end do
end do

call prin2("maximum evaluation error = ",dmax)

!
!  Test the construction of interpolation matrices
!


nlege = 30
nout  = nlege*nlege

allocate(xsout(nout),ysout(nout),whtsout(nout))
allocate(valsout(nout),valsout0(nout))
call legendre_quad(nlege,xslege,whtslege)

idx=0
do i=1,nlege
do j=1,nlege
xsout(idx)   = xslege(i)
ysout(idx)   = xslege(j)
whtsout(idx) = whtslege(i)*whtslege(j)
end do
end do

call bilege_interpmatrix(n,nout,xsout,ysout,whtsout,amatr)
vals     = xs**2*sqrt(whts)
valsout0 = xsout**2*sqrt(whtsout)
valsout  = matmul(amatr,vals)
derr = maxval(abs(valsout-valsout0))
call prin2("interp matrix error = ",derr)


end program
