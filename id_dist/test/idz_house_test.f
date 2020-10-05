c
c
c       dependencies: prini
c
c
        implicit none
c
        integer len
        parameter(len = 1 000 000)
c
        integer n,k,ifrescal
        real*8 scal,diffmax,difffrob,diff1,diffmax2ton,r1,rss
        complex*16 x(len),vn(len),h(len),c(len),y(len),
     1             hadjoint(len),ci,csshh
c
        r1 = 1
        ci = (0,1)
c
c
        call prini(6,13)
c
c
        print *,'Enter n (the length of the vector to reflect '
     1        //'into its first component): '
        read *,n
        call prinf('n = *',n,1)
c
c
c       Fill x with something.
c
        do k = 1,n
          x(k) = sqrt(r1*k)-ci*k
        enddo ! k
        call prin2('x = *',x,2*n)
c
c
c       Calculate the normalized Householder vector vn
c       corresponding to x.
c
        call idz_house(n,x,csshh,vn,scal)
        call prin2('csshh = *',csshh,2)
c
c
c       Build the Householder transformation matrix h from vn.
c
        call idz_housemat(n,vn,scal,h)
c
c
c       Calculate the root-sum-square of the entries of x.
c
        call ccalcrss(n,x,rss)
        call prin2('rss = *',rss,1)
c
c
c       Apply the Householder matrix for vector vn and scalar scal
c       to x, yielding y.
c
        ifrescal = 1
        call idz_houseapp(n,vn,x,ifrescal,scal,y)
        call prin2('y = *',y,2*n)
c
c
c       Check that abs(y(1)) = rss.
c
        diff1 = abs( rss-abs(y(1)) )
        diff1 = diff1/rss
        call prin2('diff1 = *',diff1,1)
c
c
c       Check that y(2) = 0, ..., y(n) = 0.
c
        diffmax2ton = 0
c
        do k = 2,n
          if(abs(y(k)) .gt. diffmax2ton) diffmax2ton = abs(y(k))
        enddo ! k
c
        diffmax2ton = diffmax2ton/rss
        call prin2('diffmax2ton = *',diffmax2ton,1)
c
c
c       Check that h adjoint(h) = _1_
c       (h adjoint(h) = _1_ because h is both symmetric and unitary).
c
        call cmatadjoint(n,h,hadjoint)
        call cmultiply(n,h,hadjoint,c)
        call ccheckid(n,c,diffmax,difffrob)
        call prin2('diffmax = *',diffmax,1)
        call prin2('difffrob = *',difffrob,1)
c
c
        stop
        end
c
c
c
c
        subroutine cmatadjoint(n,h,hadjoint)
c
c       forms the adjoint hadjoint of h.
c
c       input:
c       n -- first and second dimensions of h and hadjoint
c       h -- matrix for which the adjoint is taken
c
c       output:
c       hadjoint -- adjoint of h
c
        implicit none
        integer n,j,k
        complex*16 h(n,n),hadjoint(n,n)
c
c
        do k = 1,n
          do j = 1,n
            hadjoint(j,k) = conjg(h(k,j))
          enddo ! j
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine cmultiply(n,a,b,c)
c
c       multiplies a and b to get c.
c
c       input:
c       n -- size of a, b, and c
c       a -- n x n matrix to be applied to b
c       b -- n x n matrix to which a is applied
c
c       output:
c       c -- matrix resulting from applying a to b
c
        implicit none
        integer n,j,k,l
        complex*16 a(n,n),b(n,n),c(n,n)
c
c
        do j = 1,n
          do l = 1,n
            c(j,l) = 0
          enddo ! l
        enddo ! j
c
        do l = 1,n
          do j = 1,n
c
            do k = 1,n
              c(j,l) = c(j,l)+a(k,j)*b(l,k)
            enddo ! k
c
          enddo ! j
        enddo ! l
c
c
        return
        end
c
c
c
c
        subroutine ccheckid(n,c,diffmax,difffrob)
c
c       calculates the difference between c and the identity matrix.
c
c       input:
c       n -- size of c
c       c -- matrix that is supposed to be close to the identity
c
c       output:
c       diffmax -- maximum entrywise difference
c                  between c and the identity
c       difffrob -- root-sum-square of the entries
c                   of the matrix identity_matrix - c
c
        implicit none
        integer n,j,k
        real*8 diffmax,difffrob,diff
        complex*16 c(n,n)
c
c
        diffmax = 0
        difffrob = 0
c
        do j = 1,n
          do k = 1,n
c
            if(k .eq. j) diff = abs(1-c(k,j))
            if(k .ne. j) diff = abs(c(k,j))
c
            if(diff .gt. diffmax) diffmax = diff
            difffrob = difffrob+diff**2
c
          enddo ! k
        enddo ! j
c
        difffrob = sqrt(difffrob)
c
c
        return
        end
c
c
c
c
        subroutine disp(n,a)
c
c       displays the n x n matrix a.
c
c       input:
c       n -- size of a
c       a -- n x n matrix to be written to the output stream
c
        implicit none
        integer n,k
        real*8 a(n,n)
c
c
        do k = 1,n
          call prin2('*',a(1,k),n)
        enddo ! j
c
c
        return
        end
c
c
c
c
        subroutine ccalcrss(n,v,rss)
c
c       calculates the root-sum-square of the entries of v.
c
c       input:
c       n -- size of v
c       v -- vector whose entries are to be root-sum-squared
c
c       output:
c       rss -- root-sum-square of the entries of v
c
        implicit none
        integer n,k
        real*8 rss
        complex*16 v(n)
c
c
        rss = 0
        do k = 1,n
          rss = rss+v(k)*conjg(v(k))
        enddo ! k
        rss = sqrt(rss)
c
c
        return
        end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
