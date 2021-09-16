################################################################################
LIBS  =
#LIBS += -L/usr/local/opt/lapack/lib/ -llapack 
#INCLUDES = -I/usr/local/opt/lapack/include

################################################################################
ifeq ($(intel),1)
  FC    = ifort
  FLAGS+= -fpp 
  LIBS +=  -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include
  LIBS += -lmkl_lapack95_lp64 -lmkl_intel_lp64
  LIBS += -lmkl_intel_thread -lmkl_core
  LIBS += -liomp5 -lpthread
else
  FC    = gfortran
  FLAGS+= -cpp -ffree-line-length-0
  LIBS += -L/usr/lib -llapack
endif



################################################################################

correlation: correlation.f90 
	$(FC) $< -o $@ $(FLAGS) ${INCLUDES} $(LIBS)

clean:
	rm -rf correlation 

################################################################################

