c
c
c       dependencies: prini, idd_house, idd_qrpiv, lapack.a, blas.a,
c                     and (for the debugging code) id_rand
c
c
        implicit none
c
        integer len
        parameter(len = 4 000 000)
c
        integer m,n,k,krank,ier,ifdisp,iterations,its,n2,iu,iv,is,
     1          lwork,loop
        real*8 specnorm(10 000),eps,s(len),a(len),work(len),
     1         v(len),u(len),a2(len),a0(len)
c
c
        ifdisp = 0
c
c
        call prini(6,13)
c
c
        print *,'Enter m:'
        read *,m
        call prinf('m = *',m,1)
c
        print *,'Enter n:'
        read *,n
        call prinf('n = *',n,1)
c
c
c       Fill the matrix to be SVD'd.
c
        krank = 11
        call fill(krank,m,n,a0,s)
        call prin2('s = *',s,krank+1)
c
        if(ifdisp .eq. 1) then
c
          call prin2('a0 = *',a0,0)
          do k = 1,n
            call prin2('*',a0(n*(k-1)+1),2*m)
          enddo ! k
c
        endif
c
c
        do loop = 1,2
c
c
          do k = 1,m*n
            a(k) = a0(k)
          enddo ! k
c
c
          if(loop .eq. 1) then
c
c           SVD a.
c
            call iddr_svd(m,n,a,krank,u,v,s,ier,work)
            call prin2('s = *',s,krank)
c
          endif ! loop .eq. 1
c
c
          if(loop .eq. 2) then
c
c           SVD a.
c
            eps = .1d-8
            lwork = len
c
            call iddp_svd(lwork,eps,m,n,a,krank,iu,iv,is,work,ier)
c
c           Copy u, v, and s from work.
c
            do k = 1,krank*m
              u(k) = work(iu+k-1)
            enddo ! k
c
            do k = 1,krank*n
              v(k) = work(iv+k-1)
            enddo ! k
c
            do k = 1,krank
              s(k) = work(is+k-1)
            enddo ! k
            call prin2('s = *',s,krank)
c
          endif ! loop .eq. 2
c
c
c         Form a2 = U S V^T.
c
          call svdrecon(m,krank,u,s,n,v,a2)
c
          if(ifdisp .eq. 1) then
c
            call prin2('a2 = *',a2,0)
            do k = 1,n
              call prin2('*',a2(n*(k-1)+1),2*m)
            enddo ! k
c
          endif
c
c
c         Calculate the spectral norm of a0-a2.
c
          n2 = 2*n
          call id_srand(n2,u)
c
          eps = .1d-1
          iterations = 100
          call powerchk(eps,m,n,a0,a2,u,iterations,
     1                  its,specnorm,work)
c
          call prin2('specnorm(its) = *',specnorm(its),1)
c
c
        enddo ! loop
c
c
        stop
        end
c
c
c
c
        subroutine powerchk(eps,m,n,a,b,u,iterations,its,specnorm,t)
c
c       estimates the spectral norm || a - b ||,
c       using the power method.
c
c       input:
c       eps -- desired accuracy of the norm estimate
c       m -- first dimension of a and b
c       n -- second dimension of a and b
c       u -- vector to which a-b is to be applied;
c            u is destroyed by this routine
c       iterations -- length of specnorm; number of power method
c                     iterations to conduct
c
c       output:
c       its -- number of iterations actually conducted
c              (the lesser of iterations and the number required
c               to attain convergence to within eps)
c       specnorm -- estimates of the spectral norm || a - b ||,
c                   as provided by power method iterations
c
c       work:
c       t -- must be at least m real*8 elements long
c
        implicit none
        integer m,n,ktest,iteration,j,k,iterations,its
        real*8 eps,specnorm(iterations),rss1,rss2,a(m,n),b(m,n),
     1         u(n),t(m)
c
c
        its = iterations
c
        do iteration = 1,iterations+1
c
c
c         Note that the first iteration normalizes u,
c         but that rss1 is meaningless during the first iteration,
c         unless the input u was normalized.
c
c
          ktest = 1
          do j = 1,m
c
            t(j) = 0
c
            do k = 1,n
              t(j) = t(j)+a(j,k)*u(k)-b(j,k)*u(k)
            enddo ! j
c
          enddo ! k
c
c
          rss1 = 0
c
          do j = 1,m
            rss1 = rss1 + t(j)**2
          enddo ! j
c
          rss1 = sqrt(rss1)
c
c
          if(rss1 .gt. 0) then
            do j = 1,m
              t(j) = t(j) / rss1
            enddo ! j
          endif
c
c
          ktest = 1
          do k = 1,n
c
            u(k) = 0
c
            do j = 1,m
              u(k) = u(k)+a(j,k)*t(j)-b(j,k)*t(j)
            enddo ! j
c
          enddo ! k
c
c
          rss2 = 0
c
          do k = 1,n
            rss2 = rss2 + u(k)**2
          enddo ! k
c
          rss2 = sqrt(rss2)
c
c
          if(rss2 .gt. 0) then
            do k = 1,n
              u(k) = u(k) / rss2
            enddo ! k
          endif
c
c
          if(iteration .gt. 1) then
            specnorm(iteration-1) = sqrt(rss1*rss2)
          endif
c
          if(iteration .gt. 2) then
            if(abs(specnorm(iteration-2)-specnorm(iteration-1))
     1       .le. eps*specnorm(iteration-1)) then
             its = iteration-1
             goto 1000
            endif
          endif
c
c
        enddo ! iteration
c
 1000   continue
c
c
        return
        end
c
c
c
c
        subroutine svdrecon(m,krank,u,s,n,v,a)
c
c       forms a = u diag(s) v^T.
c
c       input:
c       m -- first dimension of u and a
c       krank -- size of s, second dimension of u,
c                and second dimension of v
c       u -- leftmost matrix in the product a = u diag(s) v^T
c       s -- entries on the diagonal in the middle matrix
c            in the product a = u diag(s) v^T
c       n -- second dimension of v^T and a
c       v -- rightmost matrix in the product a = u diag(s) v^T
c
c       output:
c       a -- matrix product u diag(s) v^T
c
        implicit none
        integer m,n,krank,j,k,l
        real*8 s(krank),u(m,krank),v(n,krank),a(m,n),sum
c
c
        do k = 1,n
          do j = 1,m
c
            sum = 0
c
            do l = 1,krank
              sum = sum+u(j,l)*s(l)*v(k,l)
            enddo ! l
c
            a(j,k) = sum
c
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
        subroutine fill(krank,m,n,a,s)
c
c       fills an m x n matrix with suitably decaying singular values,
c       and left and right singular vectors taken from the DCT-IV.
c
c       input:
c       krank -- one less than the rank of the matrix to be constructed
c       m -- first dimension of a
c       n -- second dimension of a
c
c       output:
c       a -- filled matrix
c       s -- singular values of a
c
        implicit none
        integer krank,j,k,l,m,n
        real*8 r1,pi,a(m,n),sum,s(krank+1)
c
        r1 = 1
        pi = 4*atan(r1)
c
c
c       Specify the singular values.
c
        do k = 1,krank
          s(k) = exp(log(1d-10)*(k-1)/(krank-1))
        enddo ! k
c
        s(krank+1) = 1d-10
c
c
c       Construct a.
c
        do k = 1,n
          do j = 1,m
c
            sum = 0
c
            do l = 1,krank
              sum = sum+cos(pi*(j-r1/2)*(l-r1/2)/m)*sqrt(r1*2/m)
     1                 *cos(pi*(k-r1/2)*(l-r1/2)/n)*sqrt(r1*2/n)*s(l)
            enddo ! l
c
            l = krank+1
            sum = sum+cos(pi*(j-r1/2)*(l-r1/2)/m)*sqrt(r1*2/m)
     1               *cos(pi*(k-r1/2)*(l-r1/2)/n)*sqrt(r1*2/n)*s(l)
c
            a(j,k) = sum
c
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       The above code is for testing and debugging; the remainder of
c       this file contains the following user-callable routines:
c
c
c       routine iddr_svd computes an approximation of specified rank
c       to a given matrix, in the usual SVD form U S V^T,
c       where U has orthonormal columns, V has orthonormal columns,
c       and S is diagonal.
c
c       routine iddp_svd computes an approximation of specified
c       precision to a given matrix, in the usual SVD form U S V^T,
c       where U has orthonormal columns, V has orthonormal columns,
c       and S is diagonal.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine iddr_svd(m,n,a,krank,u,v,s,ier,r)
c
c       constructs a rank-krank SVD  u diag(s) v^T  approximating a,
c       where u is an m x krank matrix whose columns are orthonormal,
c       v is an n x krank matrix whose columns are orthonormal,
c       and diag(s) is a diagonal krank x krank matrix whose entries
c       are all nonnegative. This routine combines a QR code
c       (which is based on plane/Householder reflections)
c       with the LAPACK routine dgesdd.
c
c       input:
c       m -- first dimension of a and u
c       n -- second dimension of a, and first dimension of v
c       a -- matrix to be SVD'd
c       krank -- desired rank of the approximation to a
c
c       output:
c       u -- left singular vectors of a corresponding
c            to the k greatest singular values of a
c       v -- right singular vectors of a corresponding
c            to the k greatest singular values of a
c       s -- k greatest singular values of a
c       ier -- 0 when the routine terminates successfully;
c              nonzero when the routine encounters an error
c
c       work:
c       r -- must be at least
c            (krank+2)*n+8*min(m,n)+15*krank**2+8*krank
c            real*8 elements long
c
c       _N.B._: This routine destroys a. Also, please beware that
c               the source code for this routine could be clearer.
c
        implicit none
        character*1 jobz
        integer m,n,k,krank,iftranspose,ldr,ldu,ldvt,lwork,
     1          info,j,ier,io
        real*8 a(m,n),u(m,krank),v(n*krank),s(krank),r(*)
c
c
        io = 8*min(m,n)
c
c
        ier = 0
c
c
c       Compute a pivoted QR decomposition of a.
c
        call iddr_qrpiv(m,n,a,krank,r,r(io+1))
c
c
c       Extract R from the QR decomposition.
c
        call idd_retriever(m,n,a,krank,r(io+1))
c
c
c       Rearrange R according to ind (which is stored in r).
c
        call idd_permuter(krank,r,krank,n,r(io+1))
c
c
c       Use LAPACK to SVD R,
c       storing the krank (krank x 1) left singular vectors
c       in r(io+krank*n+1 : io+krank*n+krank*krank).
c
        jobz = 'S'
        ldr = krank
        lwork = 2*(3*krank**2+n+4*krank**2+4*krank)
        ldu = krank
        ldvt = krank
c
        call dgesdd(jobz,krank,n,r(io+1),ldr,s,r(io+krank*n+1),ldu,
     1              v,ldvt,r(io+krank*n+krank*krank+1),lwork,r,info)
c
        if(info .ne. 0) then
          ier = info
          return
        endif
c
c
c       Multiply the U from R from the left by Q to obtain the U
c       for A.
c
        do k = 1,krank
c
          do j = 1,krank
            u(j,k) = r(io+krank*n+j+krank*(k-1))
          enddo ! j
c
          do j = krank+1,m
            u(j,k) = 0
          enddo ! j
c
        enddo ! k
c
        iftranspose = 0
        call idd_qmatmat(iftranspose,m,n,a,krank,krank,u,r)
c
c
c       Transpose v to obtain r.
c
        call idd_transer(krank,n,v,r)
c
c
c       Copy r into v.
c
        do k = 1,n*krank
          v(k) = r(k)
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine iddp_svd(lw,eps,m,n,a,krank,iu,iv,is,w,ier)
c
c       constructs a rank-krank SVD  U Sigma V^T  approximating a
c       to precision eps, where U is an m x krank matrix whose
c       columns are orthonormal, V is an n x krank matrix whose
c       columns are orthonormal, and Sigma is a diagonal krank x krank
c       matrix whose entries are all nonnegative.
c       The entries of U are stored in w, starting at w(iu);
c       the entries of V are stored in w, starting at w(iv).
c       The diagonal entries of Sigma are stored in w,
c       starting at w(is). This routine combines a QR code
c       (which is based on plane/Householder reflections)
c       with the LAPACK routine dgesdd.
c
c       input:
c       lw -- maximum usable length of w (in real*8 elements)
c       eps -- precision to which the SVD approximates a
c       m -- first dimension of a and u
c       n -- second dimension of a, and first dimension of v
c       a -- matrix to be SVD'd
c
c       output:
c       krank -- rank of the approximation to a
c       iu -- index in w of the first entry of the matrix
c             of orthonormal left singular vectors of a
c       iv -- index in w of the first entry of the matrix
c             of orthonormal right singular vectors of a
c       is -- index in w of the first entry of the array
c             of singular values of a
c       w -- array containing the singular values and singular vectors
c            of a; w doubles as a work array, and so must be at least
c            (krank+1)*(m+2*n+9)+8*min(m,n)+15*krank**2
c            real*8 elements long, where krank is the rank
c            output by the present routine
c       ier -- 0 when the routine terminates successfully;
c              -1000 when lw is too small;
c              other nonzero values when dgesdd bombs
c
c       _N.B._: This routine destroys a. Also, please beware that
c               the source code for this routine could be clearer.
c               w must be at least 
c               (krank+1)*(m+2*n+9)+8*min(m,n)+15*krank**2
c               real*8 elements long, where krank is the rank
c               output by the present routine.
c
        implicit none
        character*1 jobz
        integer m,n,k,krank,iftranspose,ldr,ldu,ldvt,lwork,
     1          info,j,ier,io,iu,iv,is,ivi,isi,lw,lu,lv,ls
        real*8 a(m,n),w(*),eps
c
c
        io = 8*min(m,n)
c
c
        ier = 0
c
c
c       Compute a pivoted QR decomposition of a.
c
        call iddp_qrpiv(eps,m,n,a,krank,w,w(io+1))
c
c
        if(krank .gt. 0) then
c
c
c         Extract R from the QR decomposition.
c
          call idd_retriever(m,n,a,krank,w(io+1))
c
c
c         Rearrange R according to ind (which is stored in w).
c
          call idd_permuter(krank,w,krank,n,w(io+1))
c
c
c         Use LAPACK to SVD R,
c         storing the krank (krank x 1) left singular vectors
c         in w(io+krank*n+1 : io+krank*n+krank*krank).
c
          jobz = 'S'
          ldr = krank
          lwork = 2*(3*krank**2+n+4*krank**2+4*krank)
          ldu = krank
          ldvt = krank
c
          ivi = io+krank*n+krank*krank+lwork+1
          lv = n*krank
c
          isi = ivi+lv
          ls = krank
c
          if(lw .lt. isi+ls+m*krank-1) then
            ier = -1000
            return
          endif
c
          call dgesdd(jobz,krank,n,w(io+1),ldr,w(isi),w(io+krank*n+1),
     1                ldu,w(ivi),ldvt,w(io+krank*n+krank*krank+1),
     2                lwork,w,info)
c
          if(info .ne. 0) then
            ier = info
            return
          endif
c
c
c         Transpose w(ivi:ivi+lv-1) to obtain V.
c
          iv = 1
          call idd_transer(krank,n,w(ivi),w(iv))
c
c
c         Copy w(isi:isi+ls-1) into w(is:is+ls-1).
c
          is = iv+lv
c
          do k = 1,ls
            w(is+k-1) = w(isi+k-1)
          enddo ! k
c
c
c         Multiply the U from R from the left by Q to obtain the U
c         for A.
c
          iu = is+ls
          lu = m*krank
c
          do k = 1,krank
c
            do j = 1,krank
              w(iu-1+j+krank*(k-1)) = w(io+krank*n+j+krank*(k-1))
            enddo ! j
c
          enddo ! k
c
          do k = krank,1,-1
c
            do j = m,krank+1,-1
              w(iu-1+j+m*(k-1)) = 0
            enddo ! j
c
            do j = krank,1,-1
              w(iu-1+j+m*(k-1)) = w(iu-1+j+krank*(k-1))
            enddo ! j
c
          enddo ! k
c
          iftranspose = 0
          call idd_qmatmat(iftranspose,m,n,a,krank,krank,w(iu),
     1                     w(iu+lu+1))
c
c
        endif ! krank .gt. 0
c
c
        return
        end
c
c
c
c
        subroutine idd_permuter(krank,ind,m,n,a)
c
c       permutes the columns of a according to ind obtained
c       from routine iddr_qrpiv or iddp_qrpiv, assuming that
c       a = q r from iddr_qrpiv or iddp_qrpiv.
c
c       input:
c       krank -- rank specified to routine iddr_qrpiv
c                or obtained from routine iddp_qrpiv
c       ind -- indexing array obtained from routine iddr_qrpiv
c              or iddp_qrpiv
c       m -- first dimension of a
c       n -- second dimension of a
c       a -- matrix to be rearranged
c
c       output:
c       a -- rearranged matrix
c
        implicit none
        integer k,krank,m,n,j,ind(krank)
        real*8 rswap,a(m,n)
c
c
        do k = krank,1,-1
          do j = 1,m
c
            rswap = a(j,k)
            a(j,k) = a(j,ind(k))
            a(j,ind(k)) = rswap
c
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
        subroutine idd_retriever(m,n,a,krank,r)
c
c       extracts R in the QR decomposition specified by the output a
c       of the routine iddr_qrpiv or iddp_qrpiv
c
c       input:
c       m -- first dimension of a
c       n -- second dimension of a and r
c       a -- output of routine iddr_qrpiv or iddp_qrpiv
c       krank -- rank specified to routine iddr_qrpiv,
c                or output by routine iddp_qrpiv
c
c       output:
c       r -- triangular factor in the QR decomposition specified
c            by the output a of the routine iddr_qrpiv or iddp_qrpiv
c
        implicit none
        integer m,n,j,k,krank
        real*8 a(m,n),r(krank,n)
c
c
c       Copy a into r and zero out the appropriate
c       Householder vectors that are stored in one triangle of a.
c
        do k = 1,n
          do j = 1,krank
            r(j,k) = a(j,k)
          enddo ! j
        enddo ! k
c
        do k = 1,n
          if(k .lt. krank) then
            do j = k+1,krank
              r(j,k) = 0
            enddo ! j
          endif
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_transer(m,n,a,at)
c
c       forms the transpose at of a.
c
c       input:
c       m -- first dimension of a and second dimension of at
c       n -- second dimension of a and first dimension of at
c       a -- matrix to be transposed
c
c       output:
c       at -- transpose of a
c
        implicit none
        integer m,n,j,k
        real*8 a(m,n),at(n,m)
c
c
        do k = 1,n
          do j = 1,m
            at(k,j) = a(j,k)
          enddo ! j
        enddo ! k
c
c
        return
        end
