module test_adapquad_funs

contains

subroutine testfun(n,xs,ys,vals,userptr)
use    iso_c_binding
implicit double precision (a-h,o-z)
double precision         :: xs(n),ys(n)
type(c_ptr)              :: userptr
double precision         :: vals(n)
vals = log(xs**2+ys**2)
end subroutine



subroutine testfun2(n,xs,vals,userptr)
use    iso_c_binding
implicit double precision (a-h,o-z)
double precision         :: xs(n)
type(c_ptr)              :: userptr
double precision         :: vals(n)
vals = log(xs**2)
end subroutine


end module


program test_adapquad
use utils
use adapquad
use iso_c_binding
use test_adapquad_funs

implicit double precision (a-h,o-z)
type(c_ptr)              :: userptr

x1 = 0.0d0
y1 = 0.0d0

x2 = 1.0d0
y2 = 1.0d0


eps = epsilon(0.0d0)
eps = eps*100


call adaprect(ier,eps,x1,y1,x2,y2,testfun,userptr,val,ntotal)
val0 = -0.7360564926451580713514461869020719898259151660d0


call prind("integral value = ",val)
call prin2("relative error in value = ",abs((val-val0)/val0))
call prini("nodes in quadrature = ",ntotal)
call prina("")

x1 = 0.0d0
y1 = 0.0d0
call adapint(ier,eps,x1,x2,testfun2,userptr,val,ntotal)
val0 = -2.0d0

call prind("integral value = ",val)
call prin2("relative error in value = ",abs((val-val0)/val0))
call prini("nodes in quadrature = ",ntotal)
call prina("")

end program
