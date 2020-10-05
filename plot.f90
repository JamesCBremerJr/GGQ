!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains code for creating various plots using python and matplotlib.
!  Its routine operate by constructing files which contain python source code for
!  producing PDF files.
!
!  For each routine, the name of the python source file can be specified in case the
!  user wished to make modifications to it.  If "*" is specified for this name, then
!  a random name is generated for the python source file and this file is deleted after
!  the plot has been produced.
!
!  Similarly, each routine takes as input the name of the PDF file to produce.
!  If "*" is specified by the user, then no PDF file is produced and matplotlib
!  is instructed to show the plot in a window.
!
!  The following subroutines are publicly callable:
!
!    plot_rectangles - show a collection of rectangles in the plane
!
!    plot1d_begin - begin the process of createing a plot of one-dimensional
!      functions
!
!    plot1d_add - add a single function to the p
!
!    plot1d_limits -
!
!    plot1d_xaxis - 
! 
!    plot1d_yaxis -
!
!    plot1d_end -     
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module     plot
use iso_c_binding

interface
   subroutine usleep(useconds) bind(C)
    use iso_c_binding
    implicit none
    integer(c_int32_t), value :: useconds
   end subroutine
end interface
contains


subroutine plot_rectangles(pyname, plotname, nrects, rects, icolors)
implicit double precision (a-h,o-z)
character(len=*)             :: pyname, plotname
double precision             :: rects(4,nrects)
integer, optional            :: icolors(:)
integer*4                    :: useconds
!
!  Produce a plot of a collection of rectangles in the plane.
!
!  Input parameters:
!    pyname - name of the python script to produce or "*" to
!    plotname - name of the PDF file to produce or "*" to skip the production of the
!      PDF file and simply create display the output graph
!    nrects - the number of rectangles
!    rects - a (4,nrects) array each column of which is of the form (x1,y1,x2,y2)
!      and gives the coordinates of the lower left corner (x1,y1) and upper right
!      corner (x2,y2) of one rectangle
!    icolors - an optional argument specifying the color to be used for each
!      rectangle:
!
!      icolor = 0  indicates white background with black edges
!      icolor = 1  indicates blue background with black edges
!      icolor = 2  indicates red background with black edges
!      icolor = 3  indicates yellow background with black edges
!      icolor = 4  indicates purple background with black edges
!      icolor = 5  indicates green background with black edges
!      icolor = 6  indicates orange background with black edges
!      icolor = 7  indicates brown background with black edges
!
!      If this argument is missing, then the default of 0 will be used for
!      each box.
!
!  Output parameters:
!    N/A
!


character(len=:), allocatable     :: scriptname
character(len=:), allocatable     :: command

if (pyname == "*") then
call random_number(dd)
ii = dd * (10000)
allocate(character(12) :: scriptname)
write(scriptname,"(A,I5.5,A)") "plot",ii,".py"
else
lp = len(pyname)
allocate(character(lp) :: scriptname)
scriptname = pyname
endif

iw = 201
open(iw,FILE=scriptname)

write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "import matplotlib.pyplot as plt"
write(iw,"(A)") "import matplotlib.patches as patches"


!
!  Find the maximum extents
!
xmax = -1d300
xmin =  1d300
ymax = -1d300
ymin =  1d300

do i=1,nrects
x1 = rects(1,i)
y1 = rects(2,i)
x2 = rects(3,i)
y2 = rects(4,i)

xmax = max(xmax,x1)
xmin = min(xmin,x1)
xmax = max(xmax,x2)
xmin = min(xmin,x2)
ymax = max(ymax,y1)
ymin = min(ymin,y1)
ymax = max(ymax,y2)
ymin = min(ymin,y2)
end do


write(iw,"(A)") "fig, ax = plt.subplots()"



do i=1,nrects
x1 = rects(1,i)
y1 = rects(2,i)
x2 = rects(3,i)
y2 = rects(4,i)
dx = x2-x1
dy = y2-y1



if (present(icolors)) then
if (icolors(i) .eq. 0) then
  write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'white', edgecolor = 'black'"
elseif (icolors(i) .eq. 1) then
  write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'blue', edgecolor = 'black'"
elseif (icolors(i) .eq. 2) then
  write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'red', edgecolor = 'black'"
elseif (icolors(i) .eq. 3) then
  write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'yellow', edgecolor = 'black'"
elseif (icolors(i) .eq. 4) then
  write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'purple', edgecolor = 'black'"
elseif (icolors(i) .eq. 5) then
  write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'green', edgecolor = 'black'"
elseif (icolors(i) .eq. 6) then
  write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'orange', edgecolor = 'black'"
elseif (icolors(i) .eq. 7) then
  write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'brown', edgecolor = 'black'"
endif
else
  write(iw,"('ax.add_patch(patches.Rectangle((',F24.15,',',F24.15,'),',F24.15,',',F24.15,',',A,'))')") &
  x1,y1,dx,dy,"facecolor = 'white', edgecolor = 'black'"
endif

end do

write(iw,"(A,F24.16,',',F24.16,A)") "plt.xlim(",xmin,xmax,")"
write(iw,"(A,F24.16,',',F24.16,A)") "plt.ylim(",ymin,ymax,")"

if (plotname /= "*") then

!write(iw,"(A,A,A)") 'fig.set_size_inches(10.5,10.5)'
write(iw,"(A,A,A)") 'fig.savefig("',plotname,'")'

lc = len(scriptname) + 7
allocate(character(lc) :: command )
write(command,"(A,A,A)")  "python ", scriptname
call system(command)
deallocate(command)

else
write(iw,"(A,A,A)") 'fig.set_size_inches(12,10)'
write(iw,"(A,A,A)") 'fig.canvas.set_window_title("',"",'")'
write(iw,"(A)") 'plt.show()'
lc = len(scriptname) + 9
allocate(character(lc) :: command )
write(command,"(A,A,A)")  "python ", scriptname," &"

call system(command)
deallocate(command)
useconds = 100000
call usleep(useconds)
endif

close(iw)


if (pyname == "*") then
lc = len(scriptname) + 6
allocate(character(lc) :: command )
write(command,"(A,A)")  "rm -f ", scriptname
call system(command)
endif


end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine plot_function3d(filename,n,xs,ys,zs)
implicit double precision (a-h,o-z)

character(len=*)             :: filename
double precision             :: xs(n),ys(n),zs(n)

!
!  Use Python and the pyplot package to produce a plot a function
!  of two variables.  The final output is a PDF file; a python
!  script which is executed using the system command is an intermediate
!
!  Input parameters:
!    filename - the name of the PDF file to produce
!    xlabel - a label for the x-axis
!    ylabel - a label for the y-axis
!    n - the number of point in the graph of the function to be specified
!    xs - a vector of length n speciying the x coordinate of each point to plot
!    ys - a vector of length n speciying the y coordinate of each point to plot
!
!  Output parameters:
!    N/A
!
!

character(len=12)          :: scriptname
character(len=21)          :: command

call random_number(dd)
ii = dd * (10000)

write(scriptname,"(A,I5.5,A)") "plot",ii,".py"
write(command,"(A,A,I5.5,A)")  "python ","plot",ii,".py &"

iw = 20
open(iw,FILE=scriptname)
write(iw,"(A)") "import numpy as np"
write(iw,"(A)") "import matplotlib as mpl"
write(iw,"(A)") "import matplotlib.pyplot as plt"
write(iw,"(A)") "import matplotlib.pyplot as plt"
write(iw,"(A)") "from mpl_toolkits.mplot3d import Axes3D"
! write(iw,"(A)") "mpl.use('Agg')"

! fig = plt.figure()
! ax = plt.axes(projection='3D')
! ax.plot_wireframe(X, Y, Z, color='black')
! ax.set_title('wireframe');

write(iw,"(A,I5,A)") "xs = np.zeros(",n,")"
write(iw,"(A,I5,A)") "ys = np.zeros(",n,")"
write(iw,"(A,I5,A)") "zs = np.zeros(",n,")"

do i=1,n
write(iw,"(A,I5,A,E24.16)") "xs[",i-1,"] = ",xs(i)
write(iw,"(A,I5,A,E24.16)") "ys[",i-1,"] = ",ys(i)
write(iw,"(A,I5,A,E24.16)") "zs[",i-1,"] = ",zs(i)
end do

!write(iw,"(A)") "X, Y = np.meshgrid(xs,ys)"

write(iw,"(A)") "fig = plt.figure()"
write(iw,"(A)") "ax  = plt.axes(projection='3d')"


write(iw,"(A)") "ax.scatter(xs,ys,zs)"

write(iw,"(A)")     'ax.grid()'
! write(iw,"(A,A,A)") 'fig.savefig("',filename,'")'


write(iw,"(A,A,A)") 'fig.canvas.set_window_title("',filename,'")'
write(iw,"(A)") 'plt.show()'

close(iw)

call system(command)
! call system("rm -f plotscript.py")

end subroutine

subroutine randperm(n,iperm)
implicit double precision (a-h,o-z)
integer     :: iperm(n)

do i=1,n
iperm(i) = i
end do

do i=1,n-1
call random_number(dd)
j = i+1 + (n-i-1)*dd
ival     = iperm(j)
iperm(j) = iperm(i)
iperm(i) = ival
end do

end subroutine


end module plot
