module gausssq_funs

use utils
use legendre
use bilege
use bilegepw
use gaussian_quad2d
use iso_c_binding

contains


subroutine polys_rect(norder,nquad,xs,ys,whts)
implicit double precision (a-h,o-z)
double precision, allocatable :: xs(:),ys(:),whts(:)

!
!  Construct a generalized Gaussian quadrature for bivariate polynomials of
!  a specified order given on [-1,1] x [-1,1].
!

double precision, allocatable :: xslege(:), whtslege(:)
integer, target               :: norder
type(c_ptr)                   :: userptr


pi     = acos(-1.0d0)
eps    = 1.0d-13
nfuns  = (norder+1)*(norder+2)/2

call prini("norder = ",norder)
call prini("nfuns  = ",nfuns)

!
!  Start with a product Legendre quadrature rule
!

nlege = norder+1
call legendre_quad(nlege,xslege,whtslege)
nquad = nlege*nlege
allocate(xs(nquad),ys(nquad),whts(nquad))

idx = 0
do i=1,nlege
do j=1,nlege
idx       = idx+1
xs(idx)   = xslege(i)
ys(idx)   = xslege(j)
whts(idx) = whtslege(i)*whtslege(j)
end do
end do

userptr = c_loc(norder)

call chebquad2d(nfuns,funpols,userptr,nquad,xs,ys,whts)
call prini("after chebquad2d, nquad = ",nquad)
call prin2("after chebquad2d, xs = ",xs)
call prin2("after chebquad2d, ys = ",ys)
call prin2("after chebquad2d, whts = ",whts)
call prina("")

call gaussquad2d(eps,nfuns,funpols,funregion,userptr,nquad,xs,ys,whts)

call prini("optimal nquad  = ",ceiling(nfuns/3.0d0))
call prini("nquad = ",nquad)
call prin2("xs    = ",xs)
call prin2("ys    = ",ys)
call prin2("whts  = ",whts)


deallocate(xs,ys,whts)
!end do 

end subroutine


subroutine funregion(x,y,wht,isvalid,userptr)
implicit double precision (a-h,o-z)
type (c_ptr)        :: userptr

isvalid = 1
if (x .gt. 1 .OR. x .lt. -1) isvalid = 0
if (y .gt. 1 .OR. y .lt. -1) isvalid = 0

end subroutine


subroutine funpols(nfuns,n,xs,ys,vals,dersx,dersy,userptr)
implicit double precision (a-h,o-z)
type (c_ptr)        :: userptr
double precision    :: xs(:),ys(:)
double precision    :: vals(:,:), dersx(:,:),dersy(:,:)

integer, pointer              :: norder
double precision, allocatable :: polsx(:), polsy(:)
double precision, allocatable :: poldersx(:), poldersy(:)

call c_f_pointer(userptr,norder)

npols = norder+1
allocate(polsx(0:norder), poldersx(0:norder))
allocate(polsy(0:norder), poldersy(0:norder))

do i=1,n
x   = xs(i)
y   = ys(i)

call legeders(npols,x,polsx,poldersx)
call legeders(npols,y,polsy,poldersy)

idx = 0
do nn=0,norder
do j1=0,nn
j2 = nn-j1
idx          = idx + 1
vals(i,idx)  = polsx(j1)*polsy(j2)
dersx(i,idx) = poldersx(j1)*polsy(j2) 
dersy(i,idx) = polsx(j1)*poldersy(j2)
end do
end do

end do
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!subroutine bilegepw_adapfun(nfuns,n,xs,ys,whts,vals,userptr)

subroutine funpols1(nfuns,n,xs,ys,whts,vals,userptr)
implicit double precision (a-h,o-z)
type (c_ptr)        :: userptr
double precision    :: xs(:),ys(:),whts(:)
double precision    :: vals(:,:)

integer, pointer              :: norder
double precision, allocatable :: polsx1(:),    polsy1(:)
double precision, allocatable :: polsx2(:),    polsy2(:)

call c_f_pointer(userptr,norder)
npols = norder+1


allocate(polsx1(0:norder), polsx2(0:norder))
allocate(polsy1(0:norder), polsy2(0:norder))


do i=1,n
x   = xs(i)
y   = ys(i)
wht = whts(i)

call leges(npols,x,polsx1)
call leges(npols,-x,polsx2)

call leges(npols,y,polsy1)
call leges(npols,-y,polsy2)

idx = 0
do nn=0,norder
do j1=0,nn
j2 = nn-j1
idx          = idx + 1
! vals(i,idx)  = ( polsx1(j1)*polsy1(j2) + polsx2(j1)*polsy1(j2) +  &
!                  polsx1(j1)*polsy2(j2) + polsx2(j1)*polsy2(j2) ) * sqrt(wht)/4

vals(i,idx)  = ( x**(j1)*y**(j2) + (-x)**j1*y**j2 +  &
                 x**(j1)*(-y)**j2 + (-x)**j1*(-y)**j2 ) * sqrt(wht)/4

end do
end do

end do
end subroutine


subroutine funregion1(x,y,wht,isvalid,userptr)
implicit double precision (a-h,o-z)
type (c_ptr)        :: userptr

isvalid = 1
if (x .gt. 1 .OR. x .lt.0) isvalid = 0
if (y .gt. 1 .OR. y .lt.0) isvalid = 0

end subroutine


! subroutine funpols2(nfuns,n,xs,ys,whts,vals,userptr)
! implicit double precision (a-h,o-z)
! type (c_ptr)        :: userptr
! double precision    :: xs(:),ys(:),whts(:)
! double precision    :: vals(:,:)
! integer, pointer              :: norder
! double precision, allocatable :: polsx(:), polsy(:)
! double precision, allocatable :: poldersx(:), poldersy(:)

! call c_f_pointer(userptr,norder)

! npols = norder+1
! allocate(polsx(0:norder), poldersx(0:norder))
! allocate(polsy(0:norder), poldersy(0:norder))

! do i=1,n
! x   = xs(i)
! y   = ys(i)
! wht = whts(i)
! call legeders(npols,x,polsx,poldersx)
! call legeders(npols,y,polsy,poldersy)

! idx = 0
! do nn=0,norder
! do j1=0,nn
! j2 = nn-j1
! idx          = idx + 1
! vals(i,idx)  = polsx(j1)*polsy(j2)*sqrt(wht)
! ! dersx(i,idx) = poldersx(j1)*polsy(j2) 
! ! dersy(i,idx) = polsx(j1)*poldersy(j2)
! end do
! end do

! end do
! end subroutine


subroutine funpols2(nfuns,n,xs,ys,vals,dersx,dersy,userptr)
implicit double precision (a-h,o-z)
type (c_ptr)        :: userptr
double precision    :: xs(:),ys(:)
double precision    :: vals(:,:), dersx(:,:),dersy(:,:)

integer, pointer              :: norder
double precision, allocatable :: polsx(:), polsy(:)
double precision, allocatable :: poldersx(:), poldersy(:)

call c_f_pointer(userptr,norder)

npols = norder+1
allocate(polsx(0:norder), poldersx(0:norder))
allocate(polsy(0:norder), poldersy(0:norder))

do i=1,n
x   = xs(i)
y   = ys(i)

call legeders(npols,x,polsx,poldersx)
call legeders(npols,y,polsy,poldersy)

idx = 0
do nn=0,norder
do j1=0,nn
j2 = nn-j1
idx          = idx + 1
vals(i,idx)  = polsx(j1)*polsy(j2)
dersx(i,idx) = poldersx(j1)*polsy(j2) 
dersy(i,idx) = polsx(j1)*poldersy(j2)
end do
end do

end do
end subroutine


end module



program gausssq

use utils
use plot
use legendre
use bilege
use bilegepw
use gaussian_quad2d
use gausssq_funs
use iso_c_binding

implicit double precision (a-h,o-z)
double precision, allocatable :: xs(:),ys(:),whts(:)

norder = 12
call polys_rect(norder,nquad,xs,ys,whts)

end program
