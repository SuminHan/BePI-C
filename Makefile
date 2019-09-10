# gcc compiler
CC = g++
CCFLAGS = -pthread -fopenmp -fPIC -w
LFLAGS = -fopenmp

PARALUTION_LIB =/usr/local/paralution-1.1.0/build/lib/libparalution.so
PARALUTION_INC =-I/usr/local/paralution-1.1.0/build/inc

BEPI_INC = src/mergesort.c src/queue.c src/bmpcreator.c src/LU.c

all: bepi_prep bepi_wprep bepi_run

bepi_prep:
	$(CC) $(CCFLAGS) $(PARALUTION_INC) bepi_prep.cpp $(BEPI_INC) $(PARALUTION_LIB) $(LFLAGS) -o bepi_prep

bepi_wprep:
	$(CC) $(CCFLAGS) $(PARALUTION_INC) bepi_wprep.cpp $(BEPI_INC) $(PARALUTION_LIB) $(LFLAGS) -o bepi_wprep

bepi_run:
	$(CC) $(CCFLAGS) $(PARALUTION_INC) bepi_run.cpp $(BEPI_INC) $(PARALUTION_LIB) $(LFLAGS) -o bepi_run

clean:
	rm -rf bepi_prep bepi_wprep bepi_run
