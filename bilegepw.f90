!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for constructing and manipulating piecewise bivariate 
!  Legendre expansions.  By a piecewise bivariate Legendre expansion of order n on the
!  rectangles
!
!   R_j = [a_j,b_j] x [c_j, d_j],  j=1,...,m,
!
!  we mean a sum of the form
!
!             m                 n-1                       
!    f(x) = \sum \chi   (x)    \sum     c        \tilde{P}     (x,y),                        (1)
!            l=1     R_j     0<=i+j<=n    i,j,l           i,j,l
!
!
!  where \chi_R_j is the characteristic function on the rectangle R_j and 
!  \tilde{P}_i,j,l(x) is defined via
!
!
!                               (  2n+1   )      (    2        a_j+b_j )
!   \tilde{P}_i,j,l(x,y)  = sqrt( ------- )  P_i ( ------- x + ------- )   *
!                               ( b_j-a_j )      ( b_j-a_j     a_j-b_j )
!
!
!                               (  2n+1   )      (    2        d_j+c_j )
!                           sqrt( ------- )  P_j ( ------- y + ------- ) 
!                               ( d_j-c_j )      ( d_j-c_j     d_j-c_j )
!
!
!  with P_i(x) the Legendre polynomial of degree i.  This module provides routines 
!  for handling both real-valued and complex-values expansions of this type.
!
!  The rectangles R_j are the leaf nodes in a collection of quadtree, the root nodes
!  of which are rectangles in the plane.  Expansions are represented either via the 
!  vector
!
!     ( a_0,0         )
!     ( a_1,0         )
!     ( ...           )
!     ( a_{n-1},0     )
!     ( a_0,1         )
!     ( a_1,1         )
!     ( ...           )                                                                      (3)
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
!  of the *scaled* values of the expansion at the nodes of the quadrature rule
!  constructed by amalgamating a collection of quadrature rules, one for
!  each R_j.  A list of all the leaf nodes in each quadtree is maintained
!  the ordering of the coefficients and nodes of the quadrature rule corresponds
!  to the order of leaf nodes in this list.
!
!  The quadrature rule used on each rectangle is a scaled version
!  of the generalized Gaussian rule used in bilege.f90 to represent bivariate
!  Legendre expansions on [-1,1] x [-1,1].
!
!  We refer to this mechanism for representing functions as a ``piecewise 
!  bivariate Legendre discretization scheme'' or simply a ``discretization scheme'' 
!  and the amalgamated quadrature as the discretization quadrature.
!
!  The following subroutines are publicly callable:
!
!    bilegepw_init - initialize the data structure describing a piecewise bivariate
!     Legendre discretization scheme
!
!    bilegepw_uniform - refine a discretization scheme until every rectangle in
!     the scheme has sidelength less than a specified quantity
!
!    bilegepw_adap - refine a discretization scheme until it is sufficient to
!     represent a collection of user-supplied functions to a specified precision
!
!    bilegepw_quad - return the `discretization quadrature rule' associated
!     with a discretization scheme
! 
!    bilegepw_reorder_box - reorder the list of leaf nodes so that the
!     n discretization nodes in a specified *top level* box appear as nodes
!     1 to n in the ordering of all discretzation nodes
!
!     WARNING: THIS CHANGES THE ORDERING OF THE DISCRETIZATION NODES AND HENCE
!     THE MECHANISM USED TO REPRESENT FUNCTIONS
!
!    bilegepw_rects - return the coordinates of all of the leaf rectangles in the
!     quadtree
!
!    bilegepw_coefs - compute the coefficients in one or more expansion of the
!     form (1) given their *scaled* values at the nodes of the discretization
!     quadrature rule
!
!    bilegepw_eval - evaluate one or more expansions of the form (1) at a 
!     specified point
!
!    bilegepw_evalder - evaluate one or more expansions of the form (1) and their 
!     derivatives at a specified point
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    bilegepw_reorder - reorder the leaf nodes in each quadtree so that they
!     correspond to a depth-first search with the children of each box
!     ordered as: upper left, upper right, lower right, lower left.
!
!     WARNING: THIS CHANGES THE ORDERING OF THE DISCRETIZATION NODES AND HENCE
!     THE MECHANISM USED TO REPRESENT FUNCTIONS
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module bilegepw

use bilege
use iso_c_binding

!
!  Store the 
!

type     bilegepw_disc
integer                         :: norder
integer                         :: nquad
double precision, allocatable   :: xs(:)
double precision, allocatable   :: ys(:)
double precision, allocatable   :: whts(:)
double precision, allocatable   :: u(:,:)

integer                         :: nleaves
integer                         :: maxboxes
integer                         :: nboxes
integer                         :: ntopboxes
type(bilegepw_box), allocatable :: boxes(:)
integer, allocatable            :: ileaves(:)
end type bilegepw_disc



type     bilegepw_box
double precision              :: x1, y1           ! coords of the lower left corner of 
                                                  ! the box
double precision              :: x2, y2           ! coords of the upper right corner of
                                                  ! the box

integer                       :: ilevel           ! level of the quadtree
integer                       :: ileaf            ! index in the list of leaves or 0 if not a leaf

integer                       :: ilr              ! index of the lower right child box or 0
integer                       :: ill              ! index of the lower left child box
integer                       :: iur
integer                       :: iul
end type bilegepw_box



interface

subroutine bilegepw_adapfun(nfuns,n,xs,ys,whts,vals,userptr)
import c_ptr
double precision           :: xs(:)
double precision           :: ys(:)
double precision           :: whts(:)
double precision           :: vals(:,:)
type(c_ptr)                :: userptr
end subroutine

end interface

interface     bilegepw_init
module procedure bilegepw_init1
module procedure bilegepw_init2
module procedure bilegepw_init3
end interface bilegepw_init


interface     bilegepw_coefs
module procedure bilegepw_coefs1
module procedure bilegepw_coefs2
module procedure bilegepw_coefs3
module procedure bilegepw_coefs4
end interface bilegepw_coefs

interface     bilegepw_eval
module procedure bilegepw_eval1
module procedure bilegepw_eval2
module procedure bilegepw_eval3
module procedure bilegepw_eval4
end interface bilegepw_eval

interface     bilegepw_evalder
module procedure bilegepw_evalder1
module procedure bilegepw_evalder2
module procedure bilegepw_evalder3
module procedure bilegepw_evalder4
end interface bilegepw_evalder

contains


subroutine bilegepw_init1(norder,disc,a,b,c,d)
implicit double precision (a-h,o-z)
type(bilegepw_disc), intent(out)  :: disc
!
!  Initialize the data structure describing a piecewise Legendre discretization
!  scheme and construct a scheme which consists of a single top-level box
!  [a,b] x [c,d].
!
!  Input parameters:
!    norder - the order for the bivariate Legendre expansions used on
!      each subrectangle; see bilege.f90 for a list of possible orders
!    [a,b] x [c,d] - the rectangle on which the discretization scheme is
!      defined
!
!  Output parameters:
!    disc - the data structure describing the scheme
!

maxboxes    = 10000
disc%norder = norder

call bilege_quad(norder,disc%nquad,disc%xs,disc%ys,disc%whts)
call bilege_coefsmatrix(norder,disc%u)

disc%ntopboxes = 1
disc%maxboxes  = maxboxes
disc%nboxes    = 1
disc%nleaves   = 1

allocate(disc%boxes(maxboxes))
allocate(disc%ileaves(maxboxes))

disc%ileaves(1)      = 1
disc%boxes(1)%x1     = a
disc%boxes(1)%y1     = c
disc%boxes(1)%x2     = b
disc%boxes(1)%y2     = d
disc%boxes(1)%iur    = 0
disc%boxes(1)%iul    = 0
disc%boxes(1)%ilr    = 0
disc%boxes(1)%ill    = 0
disc%boxes(1)%ileaf  = 1

end subroutine

subroutine bilegepw_init2(norder,disc,nrects,rects)
implicit double precision (a-h,o-z)
type(bilegepw_disc), intent(out)  :: disc
double precision                  :: rects(:,:)
!
!  Initialize the data structure describing a piecewise Legendre discretization
!  scheme and construct a scheme consisting of a collection of user-specified
!  rectangles with disjoint interiors.
!
!  Input parameters:
!    norder - the order for the bivariate Legendre expansions used on
!      each subrectangle; see bilege.f90 for a list of possible orders
!   nrects - the number of top level boxes
!   rects - a (4,nrects) array each column of which gives the coordinates
!     of one top level box [a,b] x [c,d]; the ordering of the coordiates
!     is (a,c,b,d)
!
!  Output parameters:
!    disc - the data structure describing the scheme
!

maxboxes    = 10000
disc%norder = norder

call bilege_quad(norder,disc%nquad,disc%xs,disc%ys,disc%whts)
call bilege_coefsmatrix(norder,disc%u)

disc%ntopboxes = nrects
disc%maxboxes  = maxboxes
disc%nboxes    = nrects
disc%nleaves   = nrects

allocate(disc%boxes(maxboxes))
allocate(disc%ileaves(maxboxes))

do i=1,nrects

disc%ileaves(i)      = i
disc%boxes(i)%x1     = rects(1,i)
disc%boxes(i)%y1     = rects(2,i)
disc%boxes(i)%x2     = rects(3,i)
disc%boxes(i)%y2     = rects(4,i)
disc%boxes(i)%iur    = 0
disc%boxes(i)%iul    = 0
disc%boxes(i)%ilr    = 0
disc%boxes(i)%ill    = 0
disc%boxes(i)%ileaf  = i

end do

end subroutine


subroutine bilegepw_init3(norder,disc,ngrid,a,b,c,d)
implicit double precision (a-h,o-z)
type(bilegepw_disc), intent(out)  :: disc
!
!  Initialize the data structure describing a piecewise Legendre discretization
!  scheme and construct a scheme 
!
!  Input parameters:
!    norder - the order for the bivariate Legendre expansions used on
!      each subrectangle; see bilege.f90 for a list of possible orders
!    ngrid - the size of the discretization grid
!    (a,b) x (c,d) - the region to decompose 
!
!  Output parameters:
!    disc - the data structure describing the scheme
!

double precision, allocatable :: rects(:,:)

nrects = ngrid*ngrid
allocate(rects(4,ngrid*ngrid))


idx = 0
do i=1,ngrid
do j=1,ngrid

a1 = a + (b-a)/(ngrid-0.0d0) * (i-1)
b1 = a + (b-a)/(ngrid-0.0d0) * i

c1 = c + (d-c)/(ngrid-0.0d0) * (j-1)
d1 = c + (d-c)/(ngrid-0.0d0) * j

idx          = idx+1
rects(1,idx) = a1
rects(2,idx) = c1
rects(3,idx) = b1
rects(4,idx) = d1

end do
end do

call bilegepw_init2(norder,disc,nrects,rects)

end subroutine


subroutine bilegepw_splitleaf(disc,ileaf,iul,iur,ilr,ill)
implicit double precision (a-h,o-z)
type(bilegepw_disc)            :: disc
!
!  Split a specified leaf box.  The box is reference by its index in the list of
!  leaves, and NOT its index into boxes array.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    ileaf - the index of the leaf box to split
!
!  Output parameters:
!    iul - the index in the list of leaves of the upper left child of the box
!    iur - the index in the list of leaves of the upper right child of the box
!    ilr - the index in the list of leaves of the lower right child of the box
!    ill - the index in the list of leaves of the lower left child of the box
!

nboxes   = disc%nboxes
maxboxes = disc%maxboxes
nleaves  = disc%nleaves

ibox     = disc%ileaves(ileaf)
x1       = disc%boxes(ibox)%x1
y1       = disc%boxes(ibox)%y1
x2       = disc%boxes(ibox)%x2
y2       = disc%boxes(ibox)%y2
ilevel   = disc%boxes(ibox)%ilevel
ileaf    = disc%boxes(ibox)%ileaf
! iul      = disc%boxes(ibox)%iul
! iur      = disc%boxes(ibox)%iur
! ilr      = disc%boxes(ibox)%ilr
! ill      = disc%boxes(ibox)%ill

if (nboxes+4 .gt. maxboxes) then
call prina("in bilegepw_split, maximum number of boxes exceeded")
stop
endif

disc%boxes(ibox)%ileaf  = 0
disc%boxes(ibox)%iul    = nboxes+1
disc%boxes(ibox)%iur    = nboxes+2
disc%boxes(ibox)%ilr    = nboxes+3
disc%boxes(ibox)%ill    = nboxes+4

disc%ileaves(ileaf)     = nboxes+1
disc%ileaves(nleaves+1) = nboxes+2
disc%ileaves(nleaves+2) = nboxes+3
disc%ileaves(nleaves+3) = nboxes+4

iul                     = ileaf
iur                     = nleaves+1
ilr                     = nleaves+2
ill                     = nleaves+3


x0 = (x1+x2)/2
y0 = (y1+y2)/2

! upper left box
nboxes = nboxes+1
disc%boxes(nboxes)%x1     = x1
disc%boxes(nboxes)%y1     = y0
disc%boxes(nboxes)%x2     = x0
disc%boxes(nboxes)%y2     = y2
disc%boxes(nboxes)%ileaf  = ileaf
disc%boxes(nboxes)%ilevel = ilevel+1
disc%boxes(nboxes)%iul    = 0
disc%boxes(nboxes)%iur    = 0
disc%boxes(nboxes)%ilr    = 0
disc%boxes(nboxes)%ill    = 0

! upper right box
nboxes = nboxes+1
disc%boxes(nboxes)%x1     = x0
disc%boxes(nboxes)%y1     = y0
disc%boxes(nboxes)%x2     = x2
disc%boxes(nboxes)%y2     = y2
disc%boxes(nboxes)%ileaf  = nleaves+1
disc%boxes(nboxes)%ilevel = ilevel+1
disc%boxes(nboxes)%iul    = 0
disc%boxes(nboxes)%iur    = 0
disc%boxes(nboxes)%ilr    = 0
disc%boxes(nboxes)%ill    = 0


! lower right box
nboxes = nboxes+1
disc%boxes(nboxes)%x1     = x0
disc%boxes(nboxes)%y1     = y1
disc%boxes(nboxes)%x2     = x2
disc%boxes(nboxes)%y2     = y0
disc%boxes(nboxes)%ileaf  = nleaves+2
disc%boxes(nboxes)%ilevel = ilevel+1
disc%boxes(nboxes)%iul    = 0
disc%boxes(nboxes)%iur    = 0
disc%boxes(nboxes)%ilr    = 0
disc%boxes(nboxes)%ill    = 0

! lower left box
nboxes = nboxes+1
disc%boxes(nboxes)%x1     = x1
disc%boxes(nboxes)%y1     = y1
disc%boxes(nboxes)%x2     = x0
disc%boxes(nboxes)%y2     = y0
disc%boxes(nboxes)%ileaf  = nleaves+3
disc%boxes(nboxes)%ilevel = ilevel+1
disc%boxes(nboxes)%iul    = 0
disc%boxes(nboxes)%iur    = 0
disc%boxes(nboxes)%ilr    = 0
disc%boxes(nboxes)%ill    = 0

nleaves                 = nleaves+3
disc%nboxes             = nboxes
disc%nleaves            = nleaves

end subroutine


subroutine bilegepw_indices(disc,ibox,nidxs,idxs)
implicit double precision (a-h,o-z)
type(bilegepw_disc)               :: disc
integer, allocatable, intent(out) :: idxs(:)
!
!  Return the indices of every discretization node in a specified box in the 
!  tree.
!  
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate Legendre 
!      discretization scheme
!    ibox - the index of the box 
!
!  Output parameters:
!   nidxs - the number  of such indices
!   idxs - a list of the indices
!
integer, allocatable      :: istack(:), ileaves(:)

maxboxes = disc%maxboxes
nquad0   = disc%nquad

!
!  Make a list of all leaves in the box
!

allocate(istack(maxboxes),ileaves(maxboxes))

nleaves        = 0
nstack         = 1
istack(nstack) = ibox

do while (nstack > 0 )
ibox   = istack(nstack)
nstack = nstack-1

ileaf  = disc%boxes(ibox)%ileaf
iul    = disc%boxes(ibox)%iul
iur    = disc%boxes(ibox)%iur
ilr    = disc%boxes(ibox)%ilr
ill    = disc%boxes(ibox)%ill
x1     = disc%boxes(ibox)%x1
y1     = disc%boxes(ibox)%y1
x2     = disc%boxes(ibox)%x2
y2     = disc%boxes(ibox)%y2



if (ileaf .eq. 0) then

nstack = nstack+1
istack(nstack) = iul

nstack = nstack+1
istack(nstack) = iur

nstack = nstack+1
istack(nstack) = ilr

nstack = nstack+1
istack(nstack) = ill

else
nleaves          = nleaves+1
ileaves(nleaves) = ileaf
endif

end do

!
!  Now make the list of indices
!

nidxs = nleaves*nquad0
allocate(idxs(nidxs))

idx   = 0
do i=1,nleaves
ileaf = ileaves(i)
do j=1,nquad0
idx                  = j + (ileaf-1)*nquad0
idxs(j+(i-1)*nquad0) = idx
end do
end do


end subroutine



subroutine bilegepw_reorder_box(disc,itopbox,npts)
implicit double precision (a-h,o-z)
type(bilegepw_disc)             :: disc
!
!  Reorder the list of leaf nodes so that the n discretization nodes in a specified top 
!  level box appear as nodes 1 to n in the ordering of all discretzation nodes.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    itopbox - the index of one of the top
!
!  Ouput parameters:
!    disc - the updated data structure
!    npts - the number of discretization nodes in the top-level box
!
!
integer, allocatable :: istack(:,:), ilist1(:), ilist2(:)

maxboxes  = disc%maxboxes
nleaves   = disc%nleaves
ntopboxes = disc%ntopboxes
nquad0    = disc%nquad

allocate( istack(2,maxboxes), ilist1(maxboxes), ilist2(maxboxes) )

!
!  Make lists of all leaves in the top-level box of interest and a list
!  of leaves in all other top boxes
!

nlist1         = 0
nlist2         = 0

nstack = ntopboxes
do i=1,nstack
istack(1,i) =  i
if (i .eq. itopbox) then
istack(2,i) =  1
else
istack(2,i) =  0
endif
end do


do while (nstack > 0 )
ibox   = istack(1,nstack)
intop  = istack(2,nstack)

nstack = nstack-1

ileaf  = disc%boxes(ibox)%ileaf
iul    = disc%boxes(ibox)%iul
iur    = disc%boxes(ibox)%iur
ilr    = disc%boxes(ibox)%ilr
ill    = disc%boxes(ibox)%ill
x1     = disc%boxes(ibox)%x1
y1     = disc%boxes(ibox)%y1
x2     = disc%boxes(ibox)%x2
y2     = disc%boxes(ibox)%y2

if (ileaf .eq. 0) then

nstack = nstack+1
istack(1,nstack) = iul
istack(2,nstack) = intop

nstack = nstack+1
istack(1,nstack) = iur
istack(2,nstack) = intop

nstack = nstack+1
istack(1,nstack) = ilr
istack(2,nstack) = intop

nstack = nstack+1
istack(1,nstack) = ill
istack(2,nstack) = intop

else
if (intop .eq. 1) then
nlist1           = nlist1+1
ilist1(nlist1)   = ibox
else
nlist2           = nlist2+1
ilist2(nlist2)   = ibox
endif

endif

end do


!
!  Reorder the list of leaves
!

do i=1,nlist1
ibox                   = ilist1(i)
disc%ileaves(i)        = ibox
disc%boxes(ibox)%ileaf = i
end do

do i=1,nlist2
ibox                   = ilist2(i)
disc%ileaves(nlist1+i) = ibox
disc%boxes(ibox)%ileaf = nlist1+i
end do

npts = nlist1*nquad0

end subroutine


subroutine bilegepw_uniform(disc,dlen)
implicit double precision (a-h,o-z)
type(bilegepw_disc)             :: disc
!
!  Recursively subdivide rectangles in the tree until the sidelengths of each
!  box is less than or equal to a specified value.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    dlen - the maximum side length for a box
!
!  Output parameters:
!    disc - the updated version of the data structure
!

integer, allocatable :: istack(:)

nleaves  = disc%nleaves
maxboxes = disc%maxboxes

allocate(istack(maxboxes))
nstack = nleaves
do i=1,nstack
istack(i) = i
end do

do while(nstack > 0) 
ileaf  = istack(nstack)
nstack = nstack-1
ibox   = disc%ileaves(ileaf)


x1     = disc%boxes(ibox)%x1
y1     = disc%boxes(ibox)%y1
x2     = disc%boxes(ibox)%x2
y2     = disc%boxes(ibox)%y2

! ilevel = disc%boxes(ibox)%ilevel
! ileaf  = disc%boxes(ibox)%ileaf
! iul    = disc%boxes(ibox)%iul
! iur    = disc%boxes(ibox)%iur
! ilr    = disc%boxes(ibox)%ilr
! ill    = disc%boxes(ibox)%ill


dd     = min(x2-x1,y2-y1)

if (dd .gt. dlen) then
call bilegepw_splitleaf(disc,ileaf,iul,iur,ilr,ill)


nstack = nstack+1
istack(nstack) = iul

nstack = nstack+1
istack(nstack) = iur

nstack = nstack+1
istack(nstack) = ilr

nstack = nstack+1
istack(nstack) = ill

endif

end do


end subroutine


subroutine bilegepw_adap(eps,disc,nfuns,fun,userptr)
implicit double precision (a-h,o-z)
type(bilegepw_disc)             :: disc
procedure(bilegepw_adapfun)     :: fun
type(c_ptr)                     :: userptr
!
!  Refine a discretization scheme until a collection of user-supplied functions is
!  represented to a specified precision.
!
!  Input parameters:
!    eps - precision for the discretization
!    disc - the data structure describing the piecewise bivariate Legendre discretization 
!      scheme
!    nfuns - the number of input functions
!    fun - an external subroutine supplying the input functions
!    userptr - a "void *" pointer which is passed to fun
!
!  Output parameters:
!    disc - the updated version of the data structure
!
!
!

integer, allocatable          :: istack(:)
double precision, allocatable :: xs0(:), ys0(:), whts0(:)
double precision, allocatable :: vals(:,:), coefs(:)

nleaves  = disc%nleaves
maxboxes = disc%maxboxes
norder   = disc%norder
nquad    = disc%nquad
kk       = norder
epssq    = eps**2


allocate(xs0(nquad),ys0(nquad),whts0(nquad))
allocate(vals(nquad,nfuns),coefs(nquad))

allocate(istack(maxboxes))
nstack = nleaves
do i=1,nstack
istack(i) = i
end do


do while(nstack > 0) 
ileaf  = istack(nstack)
nstack = nstack-1
ibox   = disc%ileaves(ileaf)

x1     = disc%boxes(ibox)%x1
y1     = disc%boxes(ibox)%y1
x2     = disc%boxes(ibox)%x2
y2     = disc%boxes(ibox)%y2

xs0    = (x2-x1)/2 * disc%xs + (x2+x1)/2
ys0    = (y2-y1)/2 * disc%ys + (y2+y1)/2
whts0  = (x2-x1)/2 * (y2-y1)/2 * disc%whts

ifsplit = 0
call fun(nfuns,nquad,xs0,ys0,whts0,vals,userptr)

do j=1,nfuns

coefs = matmul(disc%u,vals(:,j))

coefs = coefs**2
dd1   = sum(coefs(1:nquad-kk))
dd2   = sum(coefs(nquad-kk+1:nquad))
dd    = dd2/(dd1+dd2+1)

if (dd .gt. epssq) then
ifsplit = 1
exit
endif

end do


if (ifsplit .eq. 1) then
call bilegepw_splitleaf(disc,ileaf,iul,iur,ilr,ill)

nstack = nstack+1
istack(nstack) = iul

nstack = nstack+1
istack(nstack) = iur

nstack = nstack+1
istack(nstack) = ilr

nstack = nstack+1
istack(nstack) = ill

endif


end do


end subroutine



subroutine bilegepw_quad(disc,nquad,xs,ys,whts)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double precision, allocatable, intent(out) :: xs(:), ys(:), whts(:)
!
!  Return the discretization quadrature.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!
!  Output parameters:
!    (nquad,xs,ys,whts) - the discretization quadrature
!

nquad0  = disc%nquad
nleaves = disc%nleaves
nquad   = nquad0*nleaves

allocate(xs(nquad),ys(nquad),whts(nquad))

do i=1,nleaves
ibox = disc%ileaves(i)
x1   = disc%boxes(ibox)%x1
y1   = disc%boxes(ibox)%y1
x2   = disc%boxes(ibox)%x2
y2   = disc%boxes(ibox)%y2

i1   = 1 + (i-1)*nquad0
i2   = i*nquad0

xs(i1:i2)   = (x2-x1)/2 * disc%xs + (x2+x1)/2
ys(i1:i2)   = (y2-y1)/2 * disc%ys + (y2+y1)/2
whts(i1:i2) = (x2-x1)/2 * (y2-y1)/2 * disc%whts

end do

end subroutine


subroutine bilegepw_rects(disc,nrects,rects)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double precision, allocatable, intent(out) :: rects(:,:)
!
!  Return an array giving the coordinates of every leaf rectangle in the
!  each quadtree.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!
!  Output parameters:
!    nrects - the number of leaf rectangles
!    rects - a (4,nrects) array each column of which gives the coordinates
!       of one leaf rectangle
!


nrects = disc%nleaves


allocate(rects(4,nrects))
do i=1,nrects
ibox = disc%ileaves(i)
x1   = disc%boxes(ibox)%x1
y1   = disc%boxes(ibox)%y1
x2   = disc%boxes(ibox)%x2
y2   = disc%boxes(ibox)%y2


rects(1,i) = x1
rects(2,i) = y1
rects(3,i) = x2
rects(4,i) = y2

end do

end subroutine




subroutine bilegepw_coefs1(disc,vals,coefs)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double precision                           :: vals(:)
double precision                           :: coefs(:)
!
!  Given the vector (3) of *scaled* values of a real-valued expansion of the form (1),
!  compute the vector (2) of its expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    vals - the vector of *scaled* valused
!
!  Output parameters:
!    coefs- the vector of coefficients
!


norder  = disc%norder
nquad0  = disc%nquad
nleaves = disc%nleaves
nquad   = nquad0*nleaves


do i=1,nleaves
ibox = disc%ileaves(i)
x1   = disc%boxes(ibox)%x1
y1   = disc%boxes(ibox)%y1
x2   = disc%boxes(ibox)%x2
y2   = disc%boxes(ibox)%y2

i1   = 1 + (i-1)*nquad0
i2   = i*nquad0

coefs(i1:i2) = matmul(disc%u,vals(i1:i2))

end do

end subroutine


subroutine bilegepw_coefs2(disc,vals,coefs)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double complex                             :: vals(:), coefs(:)
!
!  Given the vector (3) of *scaled* values of a complex-valued expansion of the form (1),
!  compute the vector (2) of its expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    vals - the vector of *scaled* valused
!
!  Output parameters:
!    coefs- the vector of coefficients
!


norder  = disc%norder
nquad0  = disc%nquad
nleaves = disc%nleaves
nquad   = nquad0*nleaves


do i=1,nleaves
ibox = disc%ileaves(i)
x1   = disc%boxes(ibox)%x1
y1   = disc%boxes(ibox)%y1
x2   = disc%boxes(ibox)%x2
y2   = disc%boxes(ibox)%y2

i1   = 1 + (i-1)*nquad0
i2   = i*nquad0

coefs(i1:i2) = matmul(disc%u,vals(i1:i2))

end do

end subroutine

subroutine bilegepw_coefs3(disc,vals,coefs)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double complex                             :: vals(:,:), coefs(:,:)
!
!  Given the *scaled* values of a collection of complex-valued expansions of the form 
!  (1), compute their expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    vals - a matrix whose jth column gives the scaled values of the jth
!      input expansion
!
!  Output parameters:
!    coefs - a matrix whose jth column gives the expansion coefficients of the
!      jth input expansion
!
!

norder  = disc%norder
nquad0  = disc%nquad
nleaves = disc%nleaves
nquad   = nquad0*nleaves


do i=1,nleaves
ibox = disc%ileaves(i)
x1   = disc%boxes(ibox)%x1
y1   = disc%boxes(ibox)%y1
x2   = disc%boxes(ibox)%x2
y2   = disc%boxes(ibox)%y2

i1   = 1 + (i-1)*nquad0
i2   = i*nquad0

coefs(i1:i2,:) = matmul(disc%u,vals(i1:i2,:))

end do

end subroutine


subroutine bilegepw_coefs4(disc,vals,coefs)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double precision                           :: vals(:,:), coefs(:,:)
!
!  Given the *scaled* values of a collection of real-valued expansions of the form 
!  (1), compute their expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    vals - a matrix whose jth column gives the scaled values of the jth
!      input expansion
!
!  Output parameters:
!    coefs - a matrix whose jth column gives the expansion coefficients of the
!      jth input expansion
!

norder  = disc%norder
nquad0  = disc%nquad
nleaves = disc%nleaves
nquad   = nquad0*nleaves


do i=1,nleaves
ibox = disc%ileaves(i)
x1   = disc%boxes(ibox)%x1
y1   = disc%boxes(ibox)%y1
x2   = disc%boxes(ibox)%x2
y2   = disc%boxes(ibox)%y2

i1   = 1 + (i-1)*nquad0
i2   = i*nquad0

coefs(i1:i2,:) = matmul(disc%u,vals(i1:i2,:))

end do

end subroutine



subroutine bilegepw_eval1(disc,coefs,x,y,val)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double precision                           :: coefs(:)
double precision                           :: val
!
!  Evaluate a real-valued expansion of the form (1) at a specified point given the
!  vector of its expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    coefs - the vector of expansion coefficients
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion
!


norder    = disc%norder
ntopboxes = disc%ntopboxes
nquad0    = disc%nquad
nleaves   = disc%nleaves
nquad     = nquad0*nleaves

ibox = 0

do i=1,ntopboxes
x1 = disc%boxes(i)%x1
y1 = disc%boxes(i)%y1
x2 = disc%boxes(i)%x2
y2 = disc%boxes(i)%y2


ifcontains = 1
if ( x .gt. x2 .OR. x .lt. x1) ifcontains = 0
if ( y .gt. y2 .OR. y .lt. y1) ifcontains = 0

if (ifcontains .eq. 1) then
ibox = i
exit
endif

end do

if (ibox .eq. 0) then
call prina("in bilegepw_eval, point out of bounds")
call prin2("x = ",x)
call prin2("y = ",y)
stop
endif


do while (1 .eq. 1)

x1    = disc%boxes(ibox)%x1
y1    = disc%boxes(ibox)%y1
x2    = disc%boxes(ibox)%x2
y2    = disc%boxes(ibox)%y2
iul   = disc%boxes(ibox)%iul
iur   = disc%boxes(ibox)%iur
ilr   = disc%boxes(ibox)%ilr
ill   = disc%boxes(ibox)%ill
ileaf = disc%boxes(ibox)%ileaf

x0    = (x1+x2)/2
y0    = (y1+y2)/2

if (ileaf .gt. 0) exit

if (x .gt. x0 ) then

  if (y .gt. y0 ) then
  ibox = iur
  else
  ibox = ilr
  endif

else

  if (y .gt. y0 ) then
  ibox = iul
  else
  ibox = ill
  endif

endif


end do


xx = (2*x - (x2+x1) ) /(x2-x1)
yy = (2*y - (y2+y1) ) /(y2-y1)

i1 = 1+(ileaf-1)*nquad0
i2 = ileaf*nquad0

call bilege_eval(norder,coefs(i1:i2),xx,yy,val)

dd  = sqrt(2/(x2-x1))*sqrt(2/(y2-y1))
val = val * (dd)

end subroutine


subroutine bilegepw_eval2(disc,coefs,x,y,val)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double complex                             :: coefs(:), val
!
!  Evaluate a complex-valued expansion of the form (1) at a specified point given the
!  vector of its expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    coefs - the vector of expansion coefficients
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the expansion
!


norder    = disc%norder
ntopboxes = disc%ntopboxes
nquad0    = disc%nquad
nleaves   = disc%nleaves
nquad     = nquad0*nleaves

ibox = 0
do i=1,ntopboxes
x1 = disc%boxes(i)%x1
y1 = disc%boxes(i)%y1
x2 = disc%boxes(i)%x2
y2 = disc%boxes(i)%y2


ifcontains = 1
if ( x .gt. x2 .OR. x .lt. x1) ifcontains = 0
if ( y .gt. y2 .OR. y .lt. y1) ifcontains = 0

if (ifcontains .eq. 1) then
ibox = i
exit
endif

end do

if (ibox .eq. 0) then
ibox = 1
! call prina("in bilegepw_eval, point out of bounds")
! call prin2("x = ",x)
! call prin2("y = ",y)
! stop
endif


do while (1 .eq. 1)

x1    = disc%boxes(ibox)%x1
y1    = disc%boxes(ibox)%y1
x2    = disc%boxes(ibox)%x2
y2    = disc%boxes(ibox)%y2
iul   = disc%boxes(ibox)%iul
iur   = disc%boxes(ibox)%iur
ilr   = disc%boxes(ibox)%ilr
ill   = disc%boxes(ibox)%ill
ileaf = disc%boxes(ibox)%ileaf

x0    = (x1+x2)/2
y0    = (y1+y2)/2

if (ileaf .gt. 0) exit

if (x .gt. x0 ) then

  if (y .gt. y0 ) then
  ibox = iur
  else
  ibox = ilr
  endif

else

  if (y .gt. y0 ) then
  ibox = iul
  else
  ibox = ill
  endif

endif


end do


xx = (2*x - (x2+x1) ) /(x2-x1)
yy = (2*y - (y2+y1) ) /(y2-y1)

i1 = 1+(ileaf-1)*nquad0
i2 = ileaf*nquad0

call bilege_eval(norder,coefs(i1:i2),xx,yy,val)


dd  = sqrt(2/(x2-x1))*sqrt(2/(y2-y1))
val = val * (dd)

end subroutine


subroutine bilegepw_eval3(disc,coefs,x,y,vals)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double complex                             :: coefs(:,:), vals(:)
!
!  Evaluate a collection of complex-valued expansions of the form (1) at a 
!  specified point given their expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    coefs - the matrix whose jth column gives the expansion coefficients of the
!      jth input function
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    vals - the jth entry gives the value of the jth input expansion
!


norder    = disc%norder
ntopboxes = disc%ntopboxes
nquad0    = disc%nquad
nleaves   = disc%nleaves
nquad     = nquad0*nleaves

ibox = 0
do i=1,ntopboxes
x1 = disc%boxes(i)%x1
y1 = disc%boxes(i)%y1
x2 = disc%boxes(i)%x2
y2 = disc%boxes(i)%y2


ifcontains = 1
if ( x .gt. x2 .OR. x .lt. x1) ifcontains = 0
if ( y .gt. y2 .OR. y .lt. y1) ifcontains = 0

if (ifcontains .eq. 1) then
ibox = i
exit
endif

end do

if (ibox .eq. 0) then
ibox = 1
! call prina("in bilegepw_eval, point out of bounds")
! call prin2("x = ",x)
! call prin2("y = ",y)
! stop
endif


do while (1 .eq. 1)

x1    = disc%boxes(ibox)%x1
y1    = disc%boxes(ibox)%y1
x2    = disc%boxes(ibox)%x2
y2    = disc%boxes(ibox)%y2
iul   = disc%boxes(ibox)%iul
iur   = disc%boxes(ibox)%iur
ilr   = disc%boxes(ibox)%ilr
ill   = disc%boxes(ibox)%ill
ileaf = disc%boxes(ibox)%ileaf

x0    = (x1+x2)/2
y0    = (y1+y2)/2

if (ileaf .gt. 0) exit

if (x .gt. x0 ) then

  if (y .gt. y0 ) then
  ibox = iur
  else
  ibox = ilr
  endif

else

  if (y .gt. y0 ) then
  ibox = iul
  else
  ibox = ill
  endif

endif


end do


xx = (2*x - (x2+x1) ) /(x2-x1)
yy = (2*y - (y2+y1) ) /(y2-y1)

i1 = 1+(ileaf-1)*nquad0
i2 = ileaf*nquad0

call bilege_eval(norder,coefs(i1:i2,:),xx,yy,vals)
dd   = sqrt(2/(x2-x1))*sqrt(2/(y2-y1))
vals = vals * (dd)

end subroutine


subroutine bilegepw_eval4(disc,coefs,x,y,vals)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double precision                           :: coefs(:,:), vals(:)
!
!  Evaluate a collection of real-valued expansions of the form (1) at a 
!  specified point given their expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    coefs - the matrix whose jth column gives the expansion coefficients of the
!      jth input function
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    vals - the jth entry gives the value of the jth input expansion
!


norder    = disc%norder
ntopboxes = disc%ntopboxes
nquad0    = disc%nquad
nleaves   = disc%nleaves
nquad     = nquad0*nleaves

ibox = 0
do i=1,ntopboxes
x1 = disc%boxes(i)%x1
y1 = disc%boxes(i)%y1
x2 = disc%boxes(i)%x2
y2 = disc%boxes(i)%y2


ifcontains = 1
if ( x .gt. x2 .OR. x .lt. x1) ifcontains = 0
if ( y .gt. y2 .OR. y .lt. y1) ifcontains = 0

if (ifcontains .eq. 1) then
ibox = i
exit
endif

end do

if (ibox .eq. 0) then
ibox = 1
! call prina("in bilegepw_eval, point out of bounds")
! call prin2("x = ",x)
! call prin2("y = ",y)
! stop
endif


do while (1 .eq. 1)

x1    = disc%boxes(ibox)%x1
y1    = disc%boxes(ibox)%y1
x2    = disc%boxes(ibox)%x2
y2    = disc%boxes(ibox)%y2
iul   = disc%boxes(ibox)%iul
iur   = disc%boxes(ibox)%iur
ilr   = disc%boxes(ibox)%ilr
ill   = disc%boxes(ibox)%ill
ileaf = disc%boxes(ibox)%ileaf

x0    = (x1+x2)/2
y0    = (y1+y2)/2

if (ileaf .gt. 0) exit

if (x .gt. x0 ) then

  if (y .gt. y0 ) then
  ibox = iur
  else
  ibox = ilr
  endif

else

  if (y .gt. y0 ) then
  ibox = iul
  else
  ibox = ill
  endif

endif

end do


xx = (2*x - (x2+x1) ) /(x2-x1)
yy = (2*y - (y2+y1) ) /(y2-y1)

i1 = 1+(ileaf-1)*nquad0
i2 = ileaf*nquad0

call bilege_eval(norder,coefs(i1:i2,:),xx,yy,vals)
dd   = sqrt(2/(x2-x1))*sqrt(2/(y2-y1))
vals = vals * (dd)

end subroutine


subroutine bilegepw_evalder1(disc,coefs,x,y,val,derx,dery)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double precision                           :: coefs(:),val,derx,dery
!
!  Evaluate a real-valued expansions of the form (1) and its
!  first order derivatives at a specified point given its
!  expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    coefs - the vector of expansion coefficients
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the input expansion
!    derx - the value of the derivative w.r.t of input expansion
!    dery - the value of the derivative w.r.t of  input expansion
!


norder    = disc%norder
ntopboxes = disc%ntopboxes
nquad0    = disc%nquad
nleaves   = disc%nleaves
nquad     = nquad0*nleaves


ibox = 0
do i=1,ntopboxes
x1 = disc%boxes(i)%x1
y1 = disc%boxes(i)%y1
x2 = disc%boxes(i)%x2
y2 = disc%boxes(i)%y2


ifcontains = 1
if ( x .gt. x2 .OR. x .lt. x1) ifcontains = 0
if ( y .gt. y2 .OR. y .lt. y1) ifcontains = 0

if (ifcontains .eq. 1) then
ibox = i
exit
endif

end do

if (ibox .eq. 0) then
! call prina("in bilegepw_eval, point out of bounds")
! call prin2("x = ",x)
! call prin2("y = ",y)
! stop
ibox = 1
endif


do while (1 .eq. 1)

x1    = disc%boxes(ibox)%x1
y1    = disc%boxes(ibox)%y1
x2    = disc%boxes(ibox)%x2
y2    = disc%boxes(ibox)%y2
iul   = disc%boxes(ibox)%iul
iur   = disc%boxes(ibox)%iur
ilr   = disc%boxes(ibox)%ilr
ill   = disc%boxes(ibox)%ill
ileaf = disc%boxes(ibox)%ileaf

x0    = (x1+x2)/2
y0    = (y1+y2)/2

if (ileaf .gt. 0) exit

if (x .gt. x0 ) then

  if (y .gt. y0 ) then
  ibox = iur
  else
  ibox = ilr
  endif

else

  if (y .gt. y0 ) then
  ibox = iul
  else
  ibox = ill
  endif

endif


end do


xx = (2*x - (x2+x1) ) /(x2-x1)
yy = (2*y - (y2+y1) ) /(y2-y1)

i1 = 1+(ileaf-1)*nquad0
i2 = ileaf*nquad0

call bilege_evalder(norder,coefs(i1:i2),xx,yy,val,derx,dery)

dd   = sqrt(2/(x2-x1))*sqrt(2/(y2-y1))
val  = val * dd
derx = derx * dd * 2/(x2-x1)
dery = dery * dd * 2/(y2-y1)

end subroutine



subroutine bilegepw_evalder2(disc,coefs,x,y,val,derx,dery)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double complex                             :: coefs(:),val,derx,dery
!
!  Evaluate a complex-valued expansions of the form (1) and its
!  first order derivatives at a specified point given its
!  expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    coefs - the vector of expansion coefficients
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    val - the value of the input expansion
!    derx - the value of the derivative w.r.t of input expansion
!    dery - the value of the derivative w.r.t of  input expansion
!



norder    = disc%norder
ntopboxes = disc%ntopboxes
nquad0    = disc%nquad
nleaves   = disc%nleaves
nquad     = nquad0*nleaves

ibox = 0
do i=1,ntopboxes
x1 = disc%boxes(i)%x1
y1 = disc%boxes(i)%y1
x2 = disc%boxes(i)%x2
y2 = disc%boxes(i)%y2


ifcontains = 1
if ( x .gt. x2 .OR. x .lt. x1) ifcontains = 0
if ( y .gt. y2 .OR. y .lt. y1) ifcontains = 0

if (ifcontains .eq. 1) then
ibox = i
exit
endif

end do

if (ibox .eq. 0) then
! call prina("in bilegepw_eval, point out of bounds")
! call prin2("x = ",x)
! call prin2("y = ",y)
! stop
ibox = 1
endif

do while (1 .eq. 1)

x1    = disc%boxes(ibox)%x1
y1    = disc%boxes(ibox)%y1
x2    = disc%boxes(ibox)%x2
y2    = disc%boxes(ibox)%y2
iul   = disc%boxes(ibox)%iul
iur   = disc%boxes(ibox)%iur
ilr   = disc%boxes(ibox)%ilr
ill   = disc%boxes(ibox)%ill
ileaf = disc%boxes(ibox)%ileaf

x0    = (x1+x2)/2
y0    = (y1+y2)/2

if (ileaf .gt. 0) exit

if (x .gt. x0 ) then

  if (y .gt. y0 ) then
  ibox = iur
  else
  ibox = ilr
  endif

else

  if (y .gt. y0 ) then
  ibox = iul
  else
  ibox = ill
  endif

endif


end do


xx = (2*x - (x2+x1) ) /(x2-x1)
yy = (2*y - (y2+y1) ) /(y2-y1)

i1 = 1+(ileaf-1)*nquad0
i2 = ileaf*nquad0

call bilege_evalder(norder,coefs(i1:i2),xx,yy,val,derx,dery)

dd   = sqrt(2/(x2-x1))*sqrt(2/(y2-y1))
val  = val * dd
derx = derx * dd * 2/(x2-x1)
dery = dery * dd * 2/(y2-y1)

end subroutine



subroutine bilegepw_evalder3(disc,coefs,x,y,vals,dersx,dersy)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double precision                           :: coefs(:,:),vals(:),dersx(:),dersy(:)
!
!  Evaluate a collection of real-valued expansions of the form (1) and their
!  first order derivatives at a  specified point given their expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    coefs - the matrix whose jth column gives the expansion coefficients of the
!      jth input function
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    vals - the jth entry gives the value of the jth input expansion
!    dersx - the jth entry gives the derivative w.r.t. x of the jth expansion
!    dersy - the jth entry gives the derivative w.r.t. y of the jth expansion
!


norder    = disc%norder
ntopboxes = disc%ntopboxes
nquad0    = disc%nquad
nleaves   = disc%nleaves
nquad     = nquad0*nleaves

ibox = 0
do i=1,ntopboxes
x1 = disc%boxes(i)%x1
y1 = disc%boxes(i)%y1
x2 = disc%boxes(i)%x2
y2 = disc%boxes(i)%y2


ifcontains = 1
if ( x .gt. x2 .OR. x .lt. x1) ifcontains = 0
if ( y .gt. y2 .OR. y .lt. y1) ifcontains = 0

if (ifcontains .eq. 1) then
ibox = i
exit
endif

end do

if (ibox .eq. 0) then
! call prina("in bilegepw_eval, point out of bounds")
! call prin2("x = ",x)
! call prin2("y = ",y)
! stop
ibox = 1
endif

do while (1 .eq. 1)

x1    = disc%boxes(ibox)%x1
y1    = disc%boxes(ibox)%y1
x2    = disc%boxes(ibox)%x2
y2    = disc%boxes(ibox)%y2
iul   = disc%boxes(ibox)%iul
iur   = disc%boxes(ibox)%iur
ilr   = disc%boxes(ibox)%ilr
ill   = disc%boxes(ibox)%ill
ileaf = disc%boxes(ibox)%ileaf

x0    = (x1+x2)/2
y0    = (y1+y2)/2

if (ileaf .gt. 0) exit

if (x .gt. x0 ) then

  if (y .gt. y0 ) then
  ibox = iur
  else
  ibox = ilr
  endif

else

  if (y .gt. y0 ) then
  ibox = iul
  else
  ibox = ill
  endif

endif


end do


xx = (2*x - (x2+x1) ) /(x2-x1)
yy = (2*y - (y2+y1) ) /(y2-y1)

i1 = 1+(ileaf-1)*nquad0
i2 = ileaf*nquad0

call bilege_evalder(norder,coefs(i1:i2,:),xx,yy,vals,dersx,dersy)

dd    = sqrt(2/(x2-x1))*sqrt(2/(y2-y1))
vals  = vals * dd
dersx = dersx * dd * 2/(x2-x1)
dersy = dersy * dd * 2/(y2-y1)

end subroutine


subroutine bilegepw_evalder4(disc,coefs,x,y,vals,dersx,dersy)
implicit double precision (a-h,o-z)
type(bilegepw_disc)                        :: disc
double complex                             :: coefs(:,:),vals(:),dersx(:),dersy(:)
!
!  Evaluate a collection of complex-valued expansions of the form (1) and their
!  first order derivatives at a  specified point given their expansion coefficients.
!
!  Input parameters:
!    disc - the data structure describing the piecewise bivariate
!      Legendre discretization scheme
!    coefs - the matrix whose jth column gives the expansion coefficients of the
!      jth input function
!    (x,y) - the point at which to evaluate (1)
!
!  Output parameters:
!    vals - the jth entry gives the value of the jth input expansion
!    dersx - the jth entry gives the derivative w.r.t. x of the jth expansion
!    dersy - the jth entry gives the derivative w.r.t. y of the jth expansion
!


norder    = disc%norder
ntopboxes = disc%ntopboxes
nquad0    = disc%nquad
nleaves   = disc%nleaves
nquad     = nquad0*nleaves

ibox = 0
do i=1,ntopboxes
x1 = disc%boxes(i)%x1
y1 = disc%boxes(i)%y1
x2 = disc%boxes(i)%x2
y2 = disc%boxes(i)%y2


ifcontains = 1
if ( x .gt. x2 .OR. x .lt. x1) ifcontains = 0
if ( y .gt. y2 .OR. y .lt. y1) ifcontains = 0

if (ifcontains .eq. 1) then
ibox = i
exit
endif

end do

if (ibox .eq. 0) then
! call prina("in bilegepw_eval, point out of bounds")
! call prin2("x = ",x)
! call prin2("y = ",y)
! stop
ibox = 1
endif


do while (1 .eq. 1)

x1    = disc%boxes(ibox)%x1
y1    = disc%boxes(ibox)%y1
x2    = disc%boxes(ibox)%x2
y2    = disc%boxes(ibox)%y2
iul   = disc%boxes(ibox)%iul
iur   = disc%boxes(ibox)%iur
ilr   = disc%boxes(ibox)%ilr
ill   = disc%boxes(ibox)%ill
ileaf = disc%boxes(ibox)%ileaf

x0    = (x1+x2)/2
y0    = (y1+y2)/2

if (ileaf .gt. 0) exit

if (x .gt. x0 ) then

  if (y .gt. y0 ) then
  ibox = iur
  else
  ibox = ilr
  endif

else

  if (y .gt. y0 ) then
  ibox = iul
  else
  ibox = ill
  endif

endif


end do


xx = (2*x - (x2+x1) ) /(x2-x1)
yy = (2*y - (y2+y1) ) /(y2-y1)

i1 = 1+(ileaf-1)*nquad0
i2 = ileaf*nquad0

call bilege_evalder(norder,coefs(i1:i2,:),xx,yy,vals,dersx,dersy)

dd    = sqrt(2/(x2-x1))*sqrt(2/(y2-y1))
vals  = vals * dd
dersx = dersx * dd * 2/(x2-x1)
dersy = dersy * dd * 2/(y2-y1)

end subroutine


end module
