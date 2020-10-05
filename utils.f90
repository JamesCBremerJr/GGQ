!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This module contains various utility routines.
!
!  The following routines provide formatted output:
!
!    prin2 - output an array of doubles with 7 digits displayed 
!    prind - output an array of doubles with 15 digits displayed
!    princ - output an array of complex numbers with 7 digits displayed
!    prinz - output an array of complex numbers with 15 digits displayed
!    prini - output an array of integers
!    prinl - output an array of long integers
!    prina - output a string
!
!  The following are miscellaneous routines:
!
!    elapsed - return the wall clock time in seconds which has elapsed since some
!     arbitrary point in the past 
!    insort - sort an array of real numbers
!    insorti  - sort an array of integers
!    insort2 - sort an array a of real numbers and return an array which describes
!      the rearrangment of the the array a 
!    insorti2 - sort an array ia of integers and return an array which describes
!      the rearrangment of the array ia
!    iremove - remove from a list of integers all integers which occur in a second
!      list of integers
!    randperm - return a random permutation in S_n
!    iduplicates - remove duplicates from a sorted list of integers
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module utils

interface prin2
   module procedure prin2_0
   module procedure prin2_1
   module procedure prin2_2
   module procedure prin2_3
end interface prin2

interface prind
   module procedure prind_0
   module procedure prind_1
   module procedure prind_2
end interface prind

interface princ
   module procedure princ_0
   module procedure princ_1
   module procedure princ_2
end interface princ

interface prinz
   module procedure prinz_0
   module procedure prinz_1
   module procedure prinz_2
end interface prinz

interface prini
   module procedure prini_0
   module procedure prini_1
   module procedure prini_2
end interface prini

interface prinl
   module procedure prinl_0
   module procedure prinl_1
   module procedure prinl_2
end interface prinl

contains

subroutine prin2_0(str,a)
implicit double precision (a-h,o-z)

double precision a
character(len=*), intent(in) :: str

print *,str
print "(4(2x,e15.7))",a

write (13,*) str
write (13,"(4(2x,e15.7))") a

end subroutine

subroutine prin2_1(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(4(2x,e15.7))",a

write (13,*) str
write (13,"(4(2x,e15.7))") a

end subroutine

subroutine prin2_2(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(4(2x,e15.7))",a

write (13,*) str
write (13,"(4(2x,e15.7))") a

end subroutine


subroutine prin2_3(str,a)
implicit double precision (a-h,o-z)

double precision, intent (in) :: a(:,:,:)
character(len=*), intent(in) :: str

print *,str
print "(4(2x,e15.7))",a

write (13,*) str
write (13,"(4(2x,e15.7))") a

end subroutine

subroutine prind_0(str,a)
implicit double precision (a-h,o-z)

double precision :: a
character(len=*), intent(in) :: str

print *,str
print "(2(2x,e24.16))",a

write (13,*) str
write (13,"(2(2x,e24.16))") a

end subroutine

subroutine prind_1(str,a)
implicit double precision (a-h,o-z)

double precision,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(2(2x,e24.16))",a

write (13,*) str
write (13,"(2(2x,e24.16))") a

end subroutine


subroutine prind_2(str,a)
implicit double precision (a-h,o-z)

double precision,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(2(2x,e24.16))",a

write (13,*) str
write (13,"(2(2x,e24.16))") a

end subroutine



subroutine princ_0(str,a)
implicit double precision (a-h,o-z)
double complex               :: a
character(len=*), intent(in) :: str

print *,str
print "(2(d15.7,',',d15.7,2X))",a

write (13,*) str
write (13,"(2(d15.7,',',d15.7,2X))") a

end subroutine


subroutine princ_1(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(2(d15.7,',',d15.7,2X))",a

write (13,*) str
write (13,"(2(d15.7,',',d15.7,2X))") a

end subroutine


subroutine princ_2(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(2(d15.7,',',d15.7,2X))",a

write (13,*) str
write (13,"(2(d15.7,',',d15.7,2X))") a

end subroutine


subroutine prinz_0(str,a)
implicit double precision (a-h,o-z)
double complex               :: a
character(len=*), intent(in) :: str

print *,str
print "(2(d24.15,',',d24.15,2X))",a

write (13,*) str
write (13,"(2(d24.15,',',d24.15,2X))") a

end subroutine


subroutine prinz_1(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(2(d24.15,',',d24.15,2X))",a

write (13,*) str
write (13,"(2(d24.15,',',d24.15,2X))") a

end subroutine


subroutine prinz_2(str,a)
implicit double precision (a-h,o-z)

double complex,intent(in)    :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(2(d24.15,',',d24.15,2X))",a

write (13,*) str
write (13,"(2(d24.15,',',d24.15,2X))") a

end subroutine


subroutine prini_0(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I9))",a

write (13,*) str
write (13,"(8(2x,I9))") a

end subroutine

subroutine prini_1(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I8))",a

write (13,*) str
write (13,"(8(2x,I8))") a

end subroutine



subroutine prini_2(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(8(2x,I8))",a

write (13,*) str
write (13,"(8(2x,I8))") a

end subroutine

subroutine prina(str)
implicit double precision (a-h,o-z)

character(len=*), intent(in) :: str

print *,str
write (13,*) str

end subroutine


subroutine prinl_0(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine

subroutine prinl_1(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine

subroutine prinl_2(str,a)
implicit double precision (a-h,o-z)

integer,intent(in) :: a(:,:)
character(len=*), intent(in) :: str

print *,str
print "(6(2x,I13))",a

write (13,*) str
write (13,"(6(2x,I13))") a

end subroutine



subroutine elapsed(t)
implicit double precision (a-h,o-z)
integer*8 i,irate
real t1
call system_clock(i,irate)

dd = i
dd = dd/irate
t = dd
return
end subroutine




subroutine insort(k,a)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
double precision, intent (inout) :: a(k)
!
!  Sort an array a of k double precision numbers.
!
if (k .le. 1) return

do i=2,k
val=a(i)
j=i-1
do while (j .ge. 1 .AND. a(j) .gt. val) 
a(j+1)=a(j)
j=j-1
end do
a(j+1)=val
end do
end subroutine


subroutine insorti(k,ia)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k)
!
!  Sort an array ia of k integers.
!
if (k .le. 1) return

do i=2,k
ival=ia(i)
j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)=ia(j)
j=j-1
end do
ia(j+1)=ival
end do
end subroutine


subroutine insorti2(k,ia,idxs)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k),idxs(k)
!
!  Sort an integer array ia and return a second array idxs
!  such that the sorted array is ia(idxs), where ia is the 
!  original list of indices.
!

do i=1,k
idxs(i) = i
end do

if (k .le. 1) return

do i=2,k
ival   = ia(i)
idxval = idxs(i)
j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)   = ia(j)
idxs(j+1) = idxs(j)
j=j-1
end do
ia(j+1)   = ival
idxs(j+1) = idxval
end do

end subroutine



subroutine quicksorti2(n,ia,idxs)
implicit double precision (a-h,o-z)
dimension istack(2,10000),ia(n),idxs(n)
!
!  Sort an integer array ia and return a second array idxs
!  such that the sorted array is ia(idxs), where ia is the 
!  original list of unsorted integers.
!

maxstack = 10000
k        = 60

do i=1,n
idxs(i) = i
end do

if (n .lt. k) then
call insorti2(n,ia,idxs)
return
endif

!
nstack      = 1
istack(1,1) = 1
istack(2,1) = n
!
do while( nstack .gt. 0) 

i1 = istack(1,nstack)
i2 = istack(2,nstack)
nstack=nstack-1
!
l = i2-i1+1

if (l .le. k) then
call insorti20(l,ia(i1:i2),idxs(i1:i2))
cycle
endif
!
!  Otherwise perform quicksort step
!
call quicksorti20(ia,idxs,i1,i2,i3)
!
!  This should never happen, but just in case ...
!
if (nstack+2 .ge. maxstack) then
print *,"quicksort out of memory"
stop
endif
!
!       Make sure the smaller half is processed first to reduce storage
!       to O(logn).
!             
n1 = i3-i1+1
n2 = i2-(i3+1)+1
!
if (n2 .lt. n1) then
!
nstack = nstack+1
istack(1,nstack) = i1
istack(2,nstack) = i3
!
nstack = nstack+1
istack(1,nstack) = i3+1
istack(2,nstack) = i2
!
else
!
nstack=nstack+1
istack(1,nstack) = i3+1
istack(2,nstack) = i2
!
nstack=nstack+1
istack(1,nstack) = i1
istack(2,nstack) = i3
!
endif
!
end do

end subroutine


subroutine insorti20(k,ia,idxs)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k),idxs(k)

if (k .le. 1) return

do i=2,k
ival   = ia(i)
idxval = idxs(i)
j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)   = ia(j)
idxs(j+1) = idxs(j)
j=j-1
end do
ia(j+1)   = ival
idxs(j+1) = idxval
end do

end subroutine


subroutine quicksorti20(ivals,ivals2,i1,i2,i3)
implicit double precision (a-h,o-z)
dimension ivals(:),ivals2(:)
!
!  Randomly choose a pivot index.
!
call random_number(r)
ipiv = i1+floor((i2-i1)*r)
!
ival  = ivals(ipiv)
ival2 = ivals2(ipiv)
!
! Swap the pivot element and the last element.
!
ivals(ipiv)  = ivals(i2)
ivals2(ipiv) = ivals2(i2)
ivals(i2)    = ival
ivals2(i2)   = ival2
!
i3 = i1
!
do i=i1,i2-1
if( ivals(i) .lt. ival) then
id  = ivals(i)
id2 = ivals2(i)
!
ivals(i)   = ivals(i3)
ivals2(i)  = ivals2(i3)
ivals(i3)  = id
ivals2(i3) = id2
!
i3=i3+1
endif
end do
!
id         = ivals(i3)
id2        = ivals2(i3)
ivals(i3)  = ivals(i2)
ivals2(i3) = ivals2(i2)
ivals(i2)  = id
ivals2(i2) = id2        
!
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine quicksorti(n,ia)
implicit double precision (a-h,o-z)
dimension istack(2,10000),ia(n)
!
!  Sort an integer array ia
!

maxstack = 10000
k        = 60


if (n .lt. k) then
call insorti(n,ia)
return
endif

!
nstack      = 1
istack(1,1) = 1
istack(2,1) = n
!
do while( nstack .gt. 0) 

i1 = istack(1,nstack)
i2 = istack(2,nstack)
nstack=nstack-1

l = i2-i1+1

if (l .le. k) then
call insorti0(l,ia(i1:i2))
cycle
endif
!
!  Otherwise perform quicksort step
!
call quicksorti0(ia,i1,i2,i3)
!
!  This should never happen, but just in case ...
!
if (nstack+2 .ge. maxstack) then
print *,"quicksort out of memory"
stop
endif
!
!       Make sure the smaller half is processed first to reduce storage
!       to O(logn).
!             
n1 = i3-i1+1
n2 = i2-(i3+1)+1
!
if (n2 .lt. n1) then
!
nstack = nstack+1
istack(1,nstack) = i1
istack(2,nstack) = i3
!
nstack = nstack+1
istack(1,nstack) = i3+1
istack(2,nstack) = i2
!
else
!
nstack=nstack+1
istack(1,nstack) = i3+1
istack(2,nstack) = i2
!
nstack=nstack+1
istack(1,nstack) = i1
istack(2,nstack) = i3
!
endif
!
end do

end subroutine


subroutine insorti0(k,ia)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
integer, intent (inout)          :: ia(k)

if (k .le. 1) return

do i=2,k
ival   = ia(i)
j=i-1
do while (j .ge. 1 .AND. ia(j) .gt. ival) 
ia(j+1)   = ia(j)
j=j-1
end do
ia(j+1)   = ival
end do

end subroutine


subroutine quicksorti0(ivals,i1,i2,i3)
implicit double precision (a-h,o-z)
dimension ivals(:)
!
!  Randomly choose a pivot index.
!
call random_number(r)
ipiv = i1+floor((i2-i1)*r)
!
ival  = ivals(ipiv)
!
! Swap the pivot element and the last element.
!
ivals(ipiv)  = ivals(i2)
ivals(i2)    = ival
!
i3 = i1
!
do i=i1,i2-1
if( ivals(i) .lt. ival) then
id  = ivals(i)
!
ivals(i)   = ivals(i3)
ivals(i3)  = id
!
i3=i3+1
endif
end do
!
id         = ivals(i3)
ivals(i3)  = ivals(i2)
ivals(i2)  = id
!
 end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine insort2(k,a,idxs)
implicit double precision (a-h,o-z)
integer, intent(in)              :: k
double precision                 :: a(k)
integer, intent (inout)          :: idxs(k)
!
!  Sort a double precision array a and return a second array idxs
!  such that the sorted array is equal to a(idxs).
!


do i=1,k
idxs(i) = i
end do

if (k .le. 1) return

do i=2,k
val    = a(i)
idxval = idxs(i)
j=i-1
do while (j .ge. 1 .AND. a(j) .gt. val) 
a(j+1)    = a(j)
idxs(j+1) = idxs(j)
j=j-1
end do
a(j+1)    = val
idxs(j+1) = idxval
end do

end subroutine



subroutine iremove(n,ia,m,ib)
implicit double precision (a-h,o-z)
dimension ia(n),ib(m)
!
!       Remove from the list ia of length n all integers appearing in 
!       the list ib of length m.  Both the list ia and the list ib
!       must be sorted before this call is made.  The results will
!       also be sorted.
!

        isrc = 1
        itar = 1
        ii   = 1
 1000 continue

        if (ii .gt. m)   goto 2000
        if (isrc .gt. n) goto 3000
        if (ia(isrc) .gt. ib(ii)) then
        ii=ii+1
        goto 1000
        endif

        if (ia(isrc) .lt. ib(ii)) then         
        ia(itar) = ia(isrc)
        itar=itar+1
        isrc=isrc+1
        goto 1000
        endif
        isrc=isrc+1
        goto 1000

 2000 continue
        if (isrc .gt. n) goto 3000
        ia(itar) = ia(isrc)
        itar=itar+1
        isrc=isrc+1
        goto 2000

 3000 continue
        n = itar-1
        return        
        
end subroutine


subroutine randperm(n,iperm)
implicit double precision (a-h,o-z)
integer     :: iperm(n)
!
!  Return a random permuation of the integers 1,2,...,n.
!
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


subroutine permute_randomly(n,ilist,ilist2)
implicit double precision (a-h,o-z)
integer     :: ilist(n), ilist2(n)
!
!  Return a random permuatio nof the list ilist.
!
integer, allocatable :: iperm(:)

allocate(iperm(n))

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

ilist2 = ilist(iperm)

end subroutine



subroutine iduplicates(nin,ilistin,nout,ilistout)
implicit double precision (a-h,o-z)
integer                           :: ilistin(nin)
integer, allocatable, intent(out) :: ilistout(:)
!
!  Remove duplicates from a sorted list of integers.
!
!  Input parameters:
!
!  Output parameters:
!

integer, allocatable :: ilist(:)


if (nin .eq. 0) then
allocate(ilistout(0))
return
endif

allocate(ilist(nin))
nout = 0
idx    = 1

i             = ilistin(idx)
nout        = 1
ilist(nout) = i

do while(idx .lt. nin)
idx = idx+1

if (ilistin(idx) .ne. i) then
i             = ilistin(idx)
nout        = nout+1
ilist(nout) = i
endif

end do

allocate(ilistout(nout))
ilistout = ilist(1:nout)

end subroutine


subroutine quicksort(n,vals)
implicit double precision (a-h,o-z)
dimension istack(2,20000)
dimension vals(1),idxs(1)
!
!       Sort a list of double precision numbers.
!
k        = 60

if (n .lt. k) then
call insort(n,vals)
return
endif

maxstack = 10000

m = 1
istack(1,1) = 1
istack(2,1) = n
!
 1000 continue
if (m .eq. 0) goto 1100
i1 = istack(1,m)
i2 = istack(2,m)
m=m-1
!
l = i2-i1+1
if (l .le. k) then
call insort(l,vals(i1))
goto 1000
endif
!
! Otherwise perform quicksort.
!
call quicksort01(vals,i1,i2,i3)
!
!       This should never happen, but just in case ...
!
! if (m+2 .ge. maxstack) then
! print *,"quicksort out of memory"
! stop
! endif
!
!  Make sure the smaller half is processed first to reduce storage
!  to O(logn).
!             
n1 = i3-i1+1
n2 = i2-i3
!
if (n2 .lt. n1) then
!
m = m+1
istack(1,m) = i1
istack(2,m) = i3
!
m = m+1
istack(1,m) = i3+1
istack(2,m) = i2
!
else
!
m = m+1
istack(1,m) = i3+1
istack(2,m) = i2
!
m = m+1
istack(1,m) = i1
istack(2,m) = i3
!
endif
!
goto 1000
 1100 continue 
end subroutine


        subroutine quicksort01(vals,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension vals(1)
!
!       Randomly choose a pivot index.
!
!        call corrand3(1,r)
        call random_number(r)
        ipiv = i1+floor((i2-i1)*r)
!
!        ipiv = i1+(i2-i1)/2
!
        val  = vals(ipiv)
!
!       Swap the pivot element and the last element.
!
        vals(ipiv) = vals(i2)
        vals(i2)   = val
!
       i3 = i1
!
        do 1000 i=i1,i2-1
!
        if( vals(i) .lt. val) then
        d  = vals(i)
!
        vals(i)  = vals(i3)
        vals(i3) = d       
!
        i3=i3+1
        endif
!
 1000 continue
!
        dd = vals(i3)
        vals(i3) = vals(i2)
        vals(i2) = dd
!
        end subroutine

end module
