KIND=GNU
## KIND shoule be GNU or IFORT

ifeq ($(KIND), GNU)
FC=gfortran
FLAG=-O3
endif

ifeq ($(KIND), IFORT)
FC=ifort
FLAG=-fPIC
endif

SRC=regrid.f90
LIB=$(SRC:.f90=.cpython-37m-x86_64-linux-gnu.so)

#all: $(LIB)
#	export FFLAGS=
#	f2py --fcompiler=FC -m $(basename SRC) -c --opt=FLAG regrid.f90
all: $(LIB)

%.cpython-37m-x86_64-linux-gnu.so: %.f90
	f2py --fcompiler=$(FC) -m $* -c --opt=$(FLAG) $<
