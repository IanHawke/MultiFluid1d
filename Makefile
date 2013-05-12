# gfortran options
FC = gfortran
FFLAGS = -g -O3
FPPFLAGS =

# Intel fortran options
#FC = ifort
#FPPFLAGS = -fpp
#FFLAGS = -g -O3 -heap-arrays

OBJS =  mf_evolve.o \
	real_type_mod.o \
	output_mod.o \
	variables_mod.o \
	parameters_mod.o \
	conversion_mod.o \
	evolve_mod.o \
	initial_mod.o \
	local_cons_mod.o \
	master_function_mod.o \
	system_mod.o \
	analysis_mod.o \
	compute_fft_mod.o

OBJS_MINPACK = \
	minpack/dogleg.o \
	minpack/dpmpar.o \
	minpack/enorm.o \
	minpack/fdjac1.o \
	minpack/hybrd.o \
	minpack/hybrd1.o \
	minpack/qform.o \
	minpack/qrfac.o \
	minpack/r1mpyq.o \
	minpack/r1updt.o

# Options building included minpack
LIBS = -llapack -lfftw3
# Options building system distributed minpack
#LIBS = -llapack -lminpack -lfftw3
INCS = -I/usr/include

mf_evolve: $(OBJS) $(OBJS_MINPACK)
	$(FC) $(FFLAGS) -o mf_evolve $(OBJS) $(OBJS_MINPACK) $(INCS) $(LIBS)

clean:
	rm -f *.o *.mod minpack/*.o

cleandata:
	rm -rf results*

mf_evolve.o: real_type_mod.o output_mod.o variables_mod.o parameters_mod.o evolve_mod.o conversion_mod.o initial_mod.o analysis_mod.o  system_mod.o
analysis_mod.o: real_type_mod.o output_mod.o variables_mod.o parameters_mod.o compute_fft_mod.o master_function_mod.o
output_mod.o: real_type_mod.o  system_mod.o parameters_mod.o
parameters_mod.o: real_type_mod.o variables_mod.o
conversion_mod.o: real_type_mod.o variables_mod.o parameters_mod.o master_function_mod.o local_cons_mod.o
evolve_mod.o: real_type_mod.o parameters_mod.o variables_mod.o conversion_mod.o master_function_mod.o
initial_mod.o: real_type_mod.o parameters_mod.o variables_mod.o conversion_mod.o
master_function_mod.o: real_type_mod.o variables_mod.o parameters_mod.o
local_cons_mod.o: real_type_mod.o variables_mod.o
compute_fft_mod.o: real_type_mod.o

%.o: %.f90
	$(FC) -c $(FFLAGS) $(INCS) $< -o $@
%.o: %.F90
	$(FC) -c $(FFLAGS) $(FPPFLAGS) $(INCS) $< -o $@
