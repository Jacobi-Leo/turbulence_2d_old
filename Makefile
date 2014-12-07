FC = gfortran
FCFLAGS = -ggdb -c \
		  -Wall -Wsurprising -Wextra -Wunderflow \
		  -fbacktrace \
		  -std=legacy \
		  #-O2 \
		  #fdefault-real-8 -fdefault-double-8 \
		  
FLFLAGS = -ggdb -fbacktrace

main: decaynd.o fft1d.o mpi0.o pickvor.o ndraw.o
	$(FC) $(FLFLAGS) -o main *.o

decaynd.o: decaynd.f fft1d.o mpi0.o pickvor.o ndraw.o
	$(FC) $(FCFLAGS) decaynd.f fft1d.o mpi0.o pickvor.o ndraw.o

pickvor.o: pickvor.f ndraw.o
	$(FC) $(FCFLAGS) pickvor.f ndraw.o

ndraw.o: ndraw.f
	$(FC) $(FCFLAGS) ndraw.f 

mpi0.o: mpi0.f
	$(FC) $(FCFLAGS) mpi0.f 

run: main decay.in
	if [ -e dat ]; then rm -r dat/; fi
	mkdir ./dat/
	./main

clean:
	rm *.o

purge:
	rm -r dat/
	rm main *.o


