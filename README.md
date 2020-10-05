This is a package for constructing generalized Chebyshev and generalized Gaussian
quadrature rules for collections of univariate and bivariate functions.
A generalized  Chebyshev quadrature rule for a linearly independent collection of n 
functions f_1,...,f_n is an quadrature rule of the form

                             n
    \int f(x) dx  \approx  \sum  f(x_j) w_j                                                 (1)
       \Omega               j=1

which is exact for  the functions f_1, ..., f_n (or, at least, the formula (1) 
achieves near machine precision accuracy when one of the f_j is substituted for f). 

A generalized Gaussian quadrature rule for a linearly independent collection of 
functions f_1, ..., f_m functions in d dimensions with m = (d+1)*n is a quadrature rule 
of the form

                             n
    \int f(x) dx  \approx  \sum  f(x_j) w_j                                                 (1)
       \Omega               j=1

which is exact for each of the f_j's.

This package requires BLAS and LAPACK, as well as the interpolative decomposition
(ID) library of Martinsson, Tygert, Shkolnisky and Tygert.  The latter can be
downloaded from

   http://tygert.com/software.html

but it is included in this package in the directory id_dist.

The convenience routines in makequad.f90 are the most straightfoward to use,
and we suggest users begin by looking at the documentation there and at the
file test_makequad.f90.  

------------------------------------------------------------------------------------------

The package comprises the following files:

1.  The file utils.f90 contains some basic utility routines for timing, printing,
sorting and the like.

2.  The file adapquad.f90 contains primitive routines for adaptively evaluating
integrals; they are used only for testing other routines.

3.  The file linalg.f90 contains code for performing certain linear algebraic
operations.  These routines are principally wrappers around LAPACK and the ID
library of Tygert, et. al.

4.  The file legendre.f90 contains code for constructing Gauss-Legendre quadrature rules,
and for constructing and manipulating univariate Legendre expansions.

5.  The file legepw.f90 contains code for constructing and manipulating piecewise
Legendre expanions.

6.  The file chebquad.f90 contains code for constructing  ``generalized Chebyshev''  
quadrature rules for collections of functions given on intervals.

7. The file gaussquad.f90 contains code for constructing  ``generalized Gaussian'' 
quadrature rules for collections of functions given on intervals.

8.  The file makequad.f90 contains convenience routines which perform many necessary
steps for the user are easier to use than those found in gaussquad.f90 or chebquad.f90.

9.  The file logquads.f90 contains code for constructing a collection of quadrature rules
useful for discretizing one-dimensional integral operators whose kernels have logarithmic
singularities.

10.  The file sqquads.f90 contains a collection of quadrature rules for integrating polynomials
on squares; these quadratures are used in bilege.f90.

11.  The file bilege.f90 contains code for constructing and manipulating bivariate Legendre
expansions.

12.  The file bilegepw.f90 contains code for constructing and manipulating piecewise bivariate
Legendre expansions.

13.  The file gaussquad2d.f90 contains code for constructing generalized Chebyshev and
generalized Gaussian quadrature rules on two-dimensional domains.

14.  The file gausssq.f90 contains code for constructing generalized Gaussian quadrature
rules for polynomials given on the square [-1,1] x [-1,1].
