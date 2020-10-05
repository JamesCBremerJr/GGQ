.phony:  all PROGRAM clean

# Choose the compiler and set compiler options
OMP_NUM_THREADS      := $(shell nproc)
export OMP_NUM_THREADS

OMP_STACKSIZE        := 1024M
export OMP_STACKSIZE

# THE ID LIBRARY REQUIRES -fdefault-integer-8 if arithemtic is promoted to extended precision
FCOMP     =  gfortran
FOPTS     =  -fopenmp -Ofast -march=native -w -fexternal-blas -fallow-argument-mismatch 

#FOPTS    +=  -fdefault-real-8 -fdefault-integer-8
#FOPTS    +=  -freal-8-real-10 -fdefault-integer-8

# specify BLAS and LAPACK library
LDOPT     =  -lblas -llapack
#LDOPT     = OpenBLAS/libopenblas.a
#LDOPT     =  FLAME/libflame.a BLIS/libblis.a

FCOMPQUAD = $(FCOMP)
FC        = $(FCOMP)
FFLAGS    = $(FOPTS)

CC        = gcc
COPTS     = -O3  -I./include

export FC FFLAGS


# Set the list of programs to compile

PROGRAMS = test_adapquad test_linalg test_legendre test_legepw test_chebquad test_gaussquad \
   test_makequad logquads radquads test_bilege test_bilegepw test_gaussquad2d gausssq

# Compile all of the test programs and the library

all	      	             : clean $(PROGRAMS) 

# List the dependencies for each module's test program


GAUSSSQ_FILES                = utils.o                                                    \
                               plot.o                                                     \
                               legendre.o                                                 \
                               linalg.o                                                   \
                               sqquads.o                                                  \
                               bilege.o                                                   \
                               bilegepw.o                                                 \
                               gaussquad2d.o                                              \
                               id_lib.a


GAUSSQUAD2D_FILES            = utils.o                                                    \
                               plot.o                                                     \
                               legendre.o                                                 \
                               linalg.o                                                   \
                               sqquads.o                                                  \
                               bilege.o                                                   \
                               bilegepw.o                                                 \
                               gaussquad2d.o                                              \
                               id_lib.a


TEST_BILEGEPW_FILES          = utils.o                                                    \
                               plot.o                                                     \
                               legendre.o                                                 \
                               sqquads.o                                                  \
                               bilege.o                                                   \
                               bilegepw.o

TEST_BILEGE_FILES            = utils.o                                                    \
                               legendre.o                                                 \
                               sqquads.o                                                  \
                               bilege.o

RADQUADS_FILES               = utils.o                                                    \
                               legendre.o                                                 \
                               linalg.o                                                   \
                               legepw.o                                                   \
                               chebquad.o                                                 \
                               gaussquad.o                                                \
                               makequad.o                                                 \
                               id_lib.a

LOGQUADS_FILES               = utils.o                                                    \
                               legendre.o                                                 \
                               linalg.o                                                   \
                               legepw.o                                                   \
                               chebquad.o                                                 \
                               gaussquad.o                                                \
                               makequad.o                                                 \
                               id_lib.a

MAKEQUAD_FILES               = utils.o                                                    \
                               legendre.o                                                 \
                               adapquad.o                                                 \
                               linalg.o                                                   \
                               legepw.o                                                   \
                               chebquad.o                                                 \
                               gaussquad.o                                                \
                               makequad.o                                                 \
                               id_lib.a

GAUSSQUAD_FILES              = utils.o                                                    \
                               legendre.o                                                 \
                               linalg.o                                                   \
                               legepw.o                                                   \
                               chebquad.o                                                 \
                               gaussquad.o                                                \
                               id_lib.a

CHEBQUAD_FILES               = utils.o                                                    \
                               legendre.o                                                 \
                               linalg.o                                                   \
                               legepw.o                                                   \
                               chebquad.o                                                 \
                               id_lib.a

LEGEPW_FILES                 = utils.o                                                    \
                               legendre.o                                                 \
                               legepw.o

LEGENDRE_FILES               = utils.o                                                    \
                               legendre.o

LINALG_FILES                 = utils.o                                                    \
                               linalg.o                                                   \
                               id_lib.a

ADAPQUAD_FILES               = utils.o                                                    \
                               adapquad.o

id_lib.a	            :  
	cd id_dist; make clean; make; cp id_lib.a ..; cd ..

gausssq.o                   : $(GAUSSSQ_FILES) gausssq.f90
gausssq                     : $(GAUSSSQ_FILES) gausssq.o

test_gaussquad2d.o          : $(GAUSSQUAD2D_FILES) test_gaussquad2d.f90
test_gaussquad2d            : $(GAUSSQUAD2D_FILES) test_gaussquad2d.o

test_bilegepw.o             : $(TEST_BILEGEPW_FILES) test_bilegepw.f90
test_bilegepw               : $(TEST_BILEGEPW_FILES) test_bilegepw.o

test_bilege.o               : $(TEST_BILEGE_FILES) test_bilege.f90
test_bilege                 : $(TEST_BILEGE_FILES) test_bilege.o

radquads.o                  : $(RADQUADS_FILES) radquads.f90
radquads                    : $(RADQUADS_FILES) radquads.o

logquads.o                  : $(LOGQUADS_FILES) logquads.f90
logquads                    : $(LOGQUADS_FILES) logquads.o

test_makequad.o             : $(MAKEQUAD_FILES) test_makequad.f90
test_makequad               : $(MAKEQUAD_FILES) test_makequad.o

test_gaussquad.o            : $(GAUSSQUAD_FILES) test_gaussquad.f90
test_gaussquad              : $(GAUSSQUAD_FILES) test_gaussquad.o

test_chebquad.o             : $(CHEBQUAD_FILES) test_chebquad.f90
test_chebquad               : $(CHEBQUAD_FILES) test_chebquad.o

test_legepw.o               : $(LEGEPW_FILES) test_legepw.f90
test_legepw                 : $(LEGEPW_FILES) test_legepw.o

test_legendre.o             : $(LEGENDRE_FILES) test_legendre.f90
test_legendre               : $(LEGENDRE_FILES) test_legendre.o

test_linalg.o               : $(LINALG_FILES) test_linalg.f90
test_linalg                 : $(LINALG_FILES) test_linalg.o

test_adapquad.o             : $(ADAPQUAD_FILES) test_adapquad.f90
test_adapquad               : $(ADAPQUAD_FILES) test_adapquad.o


# Setup the general compilation rules

%		: %.o
	$(FCOMP) $(FOPTS) -o $@ $^ $(LDOPT)
	@echo  
	@echo 
	@echo "---------[ $@     ]--------------------------------------------"
	@echo 
	@./$@
	@echo 
	@echo "--------------------------------------------------------------------------"
	@echo 


%.o		: %.f90
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.f
	$(FCOMP) -c $(FOPTS)  $<

%.o		: %.cpp
	$(CPPCOMP) -c $(CPPOPTS)  $<

%.o		: %.c
	$(CC) -c $(COPTS)  $<

.PHONY: clean

clean:
	rm -f *.o *.mod *~ fort.* a.out *.a
	rm -f $(PROGRAMS)
	rm -f *.dat
	rm -f *.py
	rm -f *.pdf
	rm -f table*.tex

