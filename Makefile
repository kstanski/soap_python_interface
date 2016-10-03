CPP	= g++
CPPFLAGS	= -O3 -std=c++11 -Wall -fexceptions -msse4.2 -fopenmp
LFLAGS	= -lgomp -pthread
INCPATH	= -I../boost_1_61_0 -I../boost.simd/include

DEPS = setup.h molecule.h solver.h neighbourhood.h stats.h stratify.h structural_similarity.h descriptor.h power_spectrum.h local_similarity.h run.h
OBJ = setup.o molecule.o solver.o neighbourhood.o stats.o stratify.o structural_similarity.o descriptor.o power_spectrum.o local_similarity.o run.o

default: main.exe

%.o: %.cpp $(DEPS) $(CDEPS)
	$(CPP) $(CPPFLAGS) $(INCPATH) -c -o $@ $<

main.exe: $(OBJ) main.o
	$(CPP) -o $@ $^ $(LFLAGS)

clean:
	rm -rf *.o *.exe
