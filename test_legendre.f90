program test_legendre

use utils
use legendre

implicit double precision (a-h,o-z)

double precision, allocatable  :: xslege(:), whtslege(:), xs(:)
double precision, allocatable  :: errsint(:), ainterp(:,:), xsout(:)
double precision, allocatable  :: whtsout(:), umatr(:,:), adiff(:,:)
!double precision, allocatable  :: u(:,:), vt(:,:), sigma(:)

double complex, allocatable    :: vals(:), valsout(:), valsout0(:), coefs(:)
double complex, allocatable    :: dersout(:), dersout0(:), ders(:), ders0(:)
double complex                 :: val0, val, ima, der, der0
double complex, allocatable    :: vals2(:,:), coefs2(:,:), valsout2(:,:)


ima = (0.0d0,1.0d0)
pi  = acos(-1.0d0)

!
!  Construct the Gauss-Legendre quadrature
!

n = 24
call elapsed(t1)
call legendre_quad(n,xslege,whtslege)
call elapsed(t2)
call prin2("legendre_quad time = ",t2-t1)


!
!  Check that the appropriate polynomials are integrated correctly
!

allocate(errsint(2*n))
do i=0,2*n-1
val  = 1.0d0/(i+1.0d0) - (-1.0d0)**(i+1)/(i+1.0d0)
val0 = sum(xslege**i*whtslege)
errsint(i+1) = abs(val-val0)
end do
call prin2("maximum polynomial integration error = ",maxval(errsint))
deallocate(errsint)

!
!  Test barcyentric interpolation
!

allocate(vals(n))
vals = exp(3*ima*xslege) * sqrt(whtslege)
x    = 0.1112128912d0
call legendre_interp(n,xslege,whtslege,vals,x,val)
val0 = exp(3*ima*x)
derr = abs(val-val0)
call prin2("interpolation error = ",derr)

x    = xslege(3)
call legendre_interp(n,xslege,whtslege,vals,x,val)
val0 = exp(3*ima*x)
derr = abs(val-val0)
call prin2("interpolation error = ",derr)
deallocate(vals)

!
!  Test the construction of interpolation matrices
!

m = 66
call legendre_quad(m,xsout,whtsout)
allocate(vals(n), valsout(m), valsout0(m)  )
vals     = exp(3*ima*xslege) * sqrt(whtslege)

call legendre_interpmatrix(n,xslege,whtslege,m,xsout,whtsout,ainterp)
valsout  = matmul(ainterp,vals)
valsout0 = exp(3*ima*xsout)*sqrt(whtsout)
call prin2("max interp matrix errors = ",maxval(abs(valsout-valsout0)))
deallocate(valsout, valsout0)

!
!  Test the coefficients routine
!
allocate(coefs(n))
vals  = exp(3*ima*xslege)/(2+xslege**2)*sqrt(whtslege)
call legendre_coefs(n,xslege,whtslege,vals,coefs)
x = 0.7327382732d0
val0 = exp(3*ima*x)/(2+x**2)
der0 = ((0.0d0,1.0d0)*exp(3*ima*x)*(6 + x*((0.0d0,2.0d0) + 3*x)))/(2 + x**2)**2
call legendre_eval(n,coefs,x,val)
derr = abs(val-val0)
call prin2("eval error = ",derr)
call legendre_evalder(n,coefs,x,val,der)
derr = abs(val-val0) + abs(der-der0)
call prin2("evalder evaluation error = ",derr)


!
!  Test coefficients matrix
!
call legendre_coefsmatrix(n,xslege,whtslege,umatr)
vals  = exp(3*ima*xslege)/(2+xslege**2)*sqrt(whtslege)
coefs = matmul(umatr,vals)

x = 0.7327382732d0
val0 = exp(3*ima*x)/(2+x**2)
call legendre_eval(n,coefs,x,val)
derr = abs(val-val0)
call prin2("coefs matrix evaluation error = ",derr)

!
!  Test the vectorized versions of the evaluation routine
!

nfuns = 10000
allocate(vals2(n,nfuns), coefs2(n,nfuns), valsout(nfuns), valsout0(nfuns) )
allocate(dersout(nfuns), dersout0(nfuns) )

x0 = 0.5d0

do j=1,nfuns
vals2(:,j)  = cos(j/(nfuns+0.0d0)*xslege) * sqrt(whtslege)
valsout0(j) = cos(j/(nfuns+0.0d0)*x0)
dersout0(j) = -j/(nfuns+0.0d0)*sin(j/(nfuns+0.0d0)*x0)
end do

call legendre_coefs(n,xslege,whtslege,vals2,coefs2)
!call legendre_eval(n,coefs2,x0,valsout)
call legendre_evalder(n,coefs2,x0,valsout,dersout)
derr = maxval(abs(valsout-valsout0))
call prin2("error in vectorized eval = ",derr)
derr = maxval(abs(dersout-dersout0))
call prin2("error in vectorized evalder = ",derr)


print *,""

call elapsed(t1)
coefs2 = matmul(umatr,vals2)
call elapsed(t2)
print *,"matmul time       = ",t2-t1

call elapsed(t1)
call legendre_coefs(n,xslege,whtslege,vals2,coefs2)
call elapsed(t2)
print *,"vectorized time   = ",t2-t1

call elapsed(t1)
do j=1,nfuns
call legendre_coefs(n,xslege,whtslege,vals2(:,j),coefs2(:,j))
end do
call elapsed(t2)
print *,"unvectorized time = ",t2-t1
print *,""

m          = 1
xsout(1)   = x0
whtsout(1) = 1.0d0
call legendre_interpmatrix(n,xslege,whtslege,m,xsout,whtsout,ainterp)

call elapsed(t1)
valsout2 = matmul(ainterp,vals2)
call elapsed(t2)
print *,"matvec time       = ",t2-t1

call elapsed(t1)
do j=1,nfuns
call legendre_eval(n,coefs2(:,j),x,valsout(j))
end do
call elapsed(t2)
print *,"unvec eval time   = ",t2-t1

call elapsed(t1)
call legendre_eval(n,coefs2,x0,valsout)
call elapsed(t2)
print *,"vec eval time     = ",t2-t1
call prina("")

!
!  Test the spectral differentiation matrix
!


eps = 1.0d-15
call legendre_diffmatrix(n,xslege,whtslege,adiff)
! call svd(eps,n,n,adiff,krank,u,sigma,vt)
! rcond = sigma(1)/sigma(krank)
! call prin2("sigma = ",sigma)
! call prin2("rcond = ",rcond)

allocate(ders(n),ders0(n))
vals  = cos(xslege)  * sqrt(whtslege)
ders0 = -sin(xslege) * sqrt(whtslege)

ders  = matmul(adiff,vals)
derr = norm2(abs(ders-ders0))
call prin2("spectral diff err =",derr)

end program
