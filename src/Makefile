FC=gfortran 
FFLAGS=-O2 -fopenmp -fdefault-real-8 
LIBFLAGS=-Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl
OUT=dragon
OBJ=densities.o operations.o print.o read.o

all: ${OBJ} main.f90
		${FC} ${FFLAGS} ${OBJ} main.f90 ${LIBFLAGS} -o ${OUT}
		mkdir -p ../bin
			mv ${OUT} ../bin
				rm *.o
densities.o : densities.f90
		${FC} ${FFLAGS} -c densities.f90
operations.o : operations.f90
		${FC} ${FFLAGS} -c operations.f90
print.o : print.f90
		${FC} ${FFLAGS} -c print.f90
read.o : read.f90
		${FC} ${FFLAGS} -c read.f90
clean:
		rm -rf *.mod *.o
