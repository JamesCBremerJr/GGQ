module test_bilegepw_funs

use iso_c_binding

contains


subroutine funtest(nfuns,n,xs,ys,whts,vals,userptr)
implicit double precision (a-h,o-z)
double precision           :: xs(:)
double precision           :: ys(:)
double precision           :: whts(:)
double precision           :: vals(:,:)
type(c_ptr)                :: userptr

vals(:,1) = (1.0d0+exp(-44*(xs**2+ys**2)))*sqrt(whts)

end subroutine

end module

program test_bilegepw

use utils
use plot
use bilege
use bilegepw
use test_bilegepw_funs
use iso_c_binding

implicit double precision (a-h,o-z)

type(bilegepw_disc)             :: disc
type(c_ptr)                     :: userptr
double precision, allocatable   :: rects(:,:)
double precision, allocatable   :: xs(:), ys(:), whts(:)
double precision, allocatable   :: vals(:), coefs(:)
integer, allocatable            :: idxs(:),idxs2(:)

eps    =  1.0d-13

a      = -4.0d0
b      =  4.0d0
c      = -4.0d0
d      =  4.0d0

nfuns  = 1
ngrid  = 4

!
!  Adaptively discretize the test functions
!
norder = 20
call bilegepw_init(norder,disc,ngrid,a,b,c,d)
call bilegepw_adap(eps,disc,nfuns,funtest,userptr)

call bilegepw_rects(disc,nrects,rects)
call plot_rectangles("*", "disc.pdf", nrects, rects )
call bilegepw_quad(disc,nquad,xs,ys,whts)


call prini("nboxes = ",disc%nboxes)
call prini("nquad = ",nquad)


!
!  Check the routine which returns indices
!

ibox = 1
x1   = disc%boxes(ibox)%x1
y1   = disc%boxes(ibox)%y1
x2   = disc%boxes(ibox)%x2
y2   = disc%boxes(ibox)%y2

call bilegepw_indices(disc,ibox,nidxs,idxs)

do i=1,nidxs
idx = idxs(i)
x   = xs(idx)
y   = ys(idx)

if (x .gt. x2 .OR. x .lt. x1) then
print *,"out of box!"
print *,idx,x1,x2,y1,y2,x,y
stop
endif

if (y .gt. y2 .OR. y .lt. y1) then
print *,"out of box!"
print *,idx,x1,x2,y1,y2,x,y
stop
endif

end do


!
!  Test the accuracy of the discretization
!

allocate(vals(nquad),coefs(nquad))
vals = (xs**2+ys**2+xs*ys)*sqrt(whts)
!vals = (xs**2-ys**2)*sqrt(whts)
call bilegepw_coefs(disc,vals,coefs)

x     = 3.0d0
y     = 3.0d0
val0  = x**2+y**2 + x*y
derx0 = 2*x + y
dery0 = 2*y + x

call bilegepw_evalder(disc,coefs,x,y,val,derx,dery)

print *,val-val0
print *,derx-derx0
print *,dery-dery0


end program
