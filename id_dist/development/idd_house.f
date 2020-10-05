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
        real*8 x(len),vn(len),h(len),c(len),diffmax,difffrob,scal,
     1         y(len),rss,diff1,diffmax2ton,rsshh
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
          x(k) = -k
        enddo ! k
        call prin2('x = *',x,n)
c
c
c       Calculate the normalized Householder vector vn
c       corresponding to x.
c
        call idd_house(n,x,rsshh,vn,scal)
        call prin2('rsshh = *',rsshh,1)
c
c
c       Build the Householder transformation matrix h from vn.
c
        call idd_housemat(n,vn,scal,h)
c
c
c       Calculate the root-sum-square of the entries of x.
c
        call calcrss(n,x,rss)
        call prin2('rss = *',rss,1)
c
c
c       Apply the Householder matrix for vector vn and scalar scal
c       to x, yielding y.
c
        ifrescal = 1
        call idd_houseapp(n,vn,x,ifrescal,scal,y)
        call prin2('y = *',y,n)
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
c       Check that h^2 = _1_
c       (h^2 = _1_ because h is both symmetric and orthogonal).
c
        call multiply(n,h,h,c)
        call checkid(n,c,diffmax,difffrob)
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
        subroutine multiply(n,a,b,c)
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
        real*8 a(n,n),b(n,n),c(n,n)
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
        subroutine checkid(n,c,diffmax,difffrob)
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
        real*8 c(n,n),diffmax,difffrob,diff
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
        subroutine calcrss(n,v,rss)
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
        real*8 v(n),rss
c
c
        rss = 0
        do k = 1,n
          rss = rss+v(k)**2
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
c
c
c       routine idd_house calculates the vector and scalar
c       needed to apply the Householder tranformation reflecting
c       a given vector into its first component.
c
c       routine idd_houseapp applies a Householder matrix to a vector.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idd_houseapp(n,vn,u,ifrescal,scal,v)
c
c       applies the Householder matrix
c       identity_matrix - scal * vn * transpose(vn)
c       to the vector u, yielding the vector v;
c
c       scal = 2/(1 + vn(2)^2 + ... + vn(n)^2)
c       when vn(2), ..., vn(n) don't all vanish;
c
c       scal = 0
c       when vn(2), ..., vn(n) do all vanish
c       (including when n = 1).
c
c       input:
c       n -- size of vn, u, and v, though the indexing on vn goes
c            from 2 to n
c       vn -- components 2 to n of the Householder vector vn;
c             vn(1) is assumed to be 1 
c       u -- vector to be transformed
c       ifrescal -- set to 1 to recompute scal from vn(2), ..., vn(n);
c                   set to 0 to use scal as input
c       scal -- see the entry for ifrescal in the decription
c               of the input
c
c       output:
c       scal -- see the entry for ifrescal in the decription
c               of the input
c       v -- result of applying the Householder matrix to u;
c            it's O.K. to have v be the same as u
c            in order to apply the matrix to the vector in place
c
c       reference:
c       Golub and Van Loan, "Matrix Computations," 3rd edition,
c            Johns Hopkins University Press, 1996, Chapter 5.
c
        implicit none
        save
        integer n,k,ifrescal
        real*8 vn(2:*),scal,u(n),v(n),fact,sum
c
c
c       Get out of this routine if n = 1.
c
        if(n .eq. 1) then
          v(1) = u(1)
          return
        endif
c
c
        if(ifrescal .eq. 1) then
c
c
c         Calculate (vn(2))^2 + ... + (vn(n))^2.
c
          sum = 0
          do k = 2,n
            sum = sum+vn(k)**2
          enddo ! k
c
c
c         Calculate scal.
c
          if(sum .eq. 0) scal = 0
          if(sum .ne. 0) scal = 2/(1+sum)
c
c
        endif
c
c
c       Calculate fact = scal * transpose(vn) * u.
c
        fact = u(1)
c
        do k = 2,n
          fact = fact+vn(k)*u(k)
        enddo ! k
c
        fact = fact*scal
c
c
c       Subtract fact*vn from u, yielding v.
c      
        v(1) = u(1) - fact
c
        do k = 2,n
          v(k) = u(k) - fact*vn(k)
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_house(n,x,rss,vn,scal)
c
c       constructs the vector vn with vn(1) = 1
c       and the scalar scal such that
c       H := identity_matrix - scal * vn * transpose(vn) is orthogonal
c       and Hx = +/- e_1 * the root-sum-square of the entries of x
c       (H is the Householder matrix corresponding to x).
c
c       input:
c       n -- size of x and vn, though the indexing on vn goes
c            from 2 to n
c       x -- vector to reflect into its first component
c
c       output:
c       rss -- first entry of the vector resulting from the application
c              of the Householder matrix to x;
c              its absolute value is the root-sum-square
c              of the entries of x
c       vn -- entries 2 to n of the Householder vector vn;
c             vn(1) is assumed to be 1
c       scal -- scalar multiplying vn * transpose(vn);
c
c               scal = 2/(1 + vn(2)^2 + ... + vn(n)^2)
c               when vn(2), ..., vn(n) don't all vanish;
c
c               scal = 0
c               when vn(2), ..., vn(n) do all vanish
c               (including when n = 1)
c
c       reference:
c       Golub and Van Loan, "Matrix Computations," 3rd edition,
c            Johns Hopkins University Press, 1996, Chapter 5.
c
        implicit none
        save
        integer n,k
        real*8 x(n),rss,sum,v1,scal,vn(2:*),x1
c
c
        x1 = x(1)
c
c
c       Get out of this routine if n = 1.
c
        if(n .eq. 1) then
          rss = x1
          scal = 0
          return
        endif
c
c
c       Calculate (x(2))^2 + ... (x(n))^2
c       and the root-sum-square value of the entries in x.
c
c
        sum = 0
        do k = 2,n
          sum = sum+x(k)**2
        enddo ! k
c
c
c       Get out of this routine if sum = 0;
c       flag this case as such by setting v(2), ..., v(n) all to 0.
c
        if(sum .eq. 0) then
c
          rss = x1
          do k = 2,n
            vn(k) = 0
          enddo ! k
          scal = 0
c
          return
c
        endif
c
c
        rss = x1**2 + sum
        rss = sqrt(rss)
c
c
c       Determine the first component v1
c       of the unnormalized Householder vector
c       v = x - rss * (1 0 0 ... 0 0)^T.
c
c       If x1 <= 0, then form x1-rss directly,
c       since that expression cannot involve any cancellation.
c
        if(x1 .le. 0) v1 = x1-rss
c
c       If x1 > 0, then use the fact that
c       x1-rss = -sum / (x1+rss),
c       in order to avoid potential cancellation.
c
        if(x1 .gt. 0) v1 = -sum / (x1+rss)
c
c
c       Compute the vector vn and the scalar scal such that vn(1) = 1
c       in the Householder transformation
c       identity_matrix - scal * vn * transpose(vn).
c
        do k = 2,n
          vn(k) = x(k)/v1
        enddo ! k
c
c       scal = 2
c            / ( vn(1)^2 + vn(2)^2 + ... + vn(n)^2 )
c
c            = 2
c            / ( 1 + vn(2)^2 + ... + vn(n)^2 )
c
c            = 2*v(1)^2
c            / ( v(1)^2 + (v(1)*vn(2))^2 + ... + (v(1)*vn(n))^2 )
c
c            = 2*v(1)^2
c            / ( v(1)^2 + (v(2)^2 + ... + v(n)^2) )
c
        scal = 2*v1**2 / (v1**2+sum)
c
c
        return
        end
c
c
c
c
        subroutine idd_housemat(n,vn,scal,h)
c
c       fills h with the Householder matrix
c       identity_matrix - scal * vn * transpose(vn).
c
c       input:
c       n -- size of vn and h, though the indexing of vn goes
c            from 2 to n
c       vn -- entries 2 to n of the vector vn;
c             vn(1) is assumed to be 1
c       scal -- scalar multiplying vn * transpose(vn)
c
c       output:
c       h -- identity_matrix - scal * vn * transpose(vn)
c
        implicit none
        save
        integer n,j,k
        real*8 vn(2:*),h(n,n),scal,factor1,factor2
c
c
c       Fill h with the identity matrix.
c
        do j = 1,n
          do k = 1,n
c
            if(j .eq. k) h(k,j) = 1
            if(j .ne. k) h(k,j) = 0
c
          enddo ! k
        enddo ! j
c
c
c       Subtract from h the matrix scal*vn*transpose(vn).
c
        do j = 1,n
          do k = 1,n
c
            if(j .eq. 1) factor1 = 1
            if(j .ne. 1) factor1 = vn(j)
c
            if(k .eq. 1) factor2 = 1
            if(k .ne. 1) factor2 = vn(k)
c
            h(k,j) = h(k,j) - scal*factor1*factor2
c
          enddo ! k
        enddo ! j
c
c
        return
        end
