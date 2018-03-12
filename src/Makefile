RM = rm -f
OBJ=kink_rwh.o gen_data.o io.o rand.o hamiltonian.o kink.o mod_sparse_matrix.o checkpointing.o timing.o mod_random.o vartypes.o memory.o

HDF5_LIBFLAGS=$(shell pkg-config --libs hdf5)
HDF5_INCFLAGS=$(shell pkg-config --cflags-only-I hdf5)

hdf5 = -lhdf5_fortran $(HDF5_LIBFLAGS) -lz

#DEBUGFLAGS = -g -pg

CFLAGS = -O2 -Wall -fopenmp $(DEBUGFLAGS) -std=c99
FFLAGS = -O2 -Wall -fopenmp $(DEBUGFLAGS) -fbounds-check -m128bit-long-double -ffree-line-length-none $(COVERAGE_FLAGS) -std=gnu -DKINK_GSL

CC=mpicc $(CFLAGS)
F90=mpif90 $(FFLAGS) $(HDF5_INCFLAGS)

LDFLAGS= $(FFLAGS) -llapack -lgmp -lgsl -lblas -fbounds-check
EXENAME = kink_rwh

all: ${EXENAME}

${EXENAME}: ${OBJ} Makefile
	$(F90) -o ${EXENAME} ${OBJ} ${LDFLAGS} ${hdf5}

%.o %.mod: %.F90
	$(F90) -c -DTRAD_CPP $<

mod_sparse_matrix_int.mod mod_sparse_matrix_dp.mod mod_sparse_matrix_ham.mod: mod_sparse_matrix.F90 mod_sparse_matrix.o 

mod_sparse_matrix.o: checkpointing.mod
timing.o: checkpointing.mod
hamiltonian.o: kink.mod mod_sparse_matrix_int.mod mod_sparse_matrix_dp.mod mod_sparse_matrix_ham.mod checkpointing.mod
mod_random.o: vartypes.mod checkpointing.mod
checkpointing.o: vartypes.mod
gen_data.o: io.mod vartypes.mod kink.mod
kink_rwh.o: gen_data.o hamiltonian.o kink.o timing.o mod_random.o memory.o

rand.o: rand.c
	$(CC) -c -o rand.o rand.c



clean:
	@rm -f *.e *.o *.mod

rm: 
	@ rm -f *.dat fort.*


