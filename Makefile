#
# Bridges - PSC
#
# Intel Compilers are loaded by default; for other compilers please check the module list
#
CC = icpc
MPCC = mpicxx -cxx=icpc
OPENMP = -qopenmp
CFLAGS = -g -O3 -xHOST -Ofast -march=core-avx2 -funroll-loops
LIBS =


TARGETS = serial mpi autograder

all:	$(TARGETS)

serial: serial.o common.o
	$(CC) -o $@ $(LIBS) serial.o common.o
autograder: autograder.o common.o
	$(CC) -o $@ $(LIBS) autograder.o common.o
mpi: mpi.o common.o
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) mpi.o common.o

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(CFLAGS) mpi.cpp
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt *.optrpt
