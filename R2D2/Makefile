#SERVER=OFP
SERVER=XC

ifeq ($(SERVER), OFP)
FC=intelem
FLAG=-fPIC
PRE_EXE=export FFLAGS=
endif

ifeq ($(SERVER), XC)
FC=ifort
FLAG=-fPIC
PRE_EXE=
endif

SRC=regrid.f90
LIB=$(SRC:.f90=.cpython-37m-x86_64-linux-gnu.so)

all: $(LIB)
%.cpython-37m-x86_64-linux-gnu.so: %.f90
	$(PRE_EXE)
	f2py --fcompiler=$(FC) -m $* -c --opt=$(FLAG) $<

clean:
	rm -rf *.cpython-37m-x86_64-linux-gnu.so
