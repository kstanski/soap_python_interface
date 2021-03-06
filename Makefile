CPP	= g++
CPPFLAGS	= -O3 -fPIC -std=c++11 -Wall -fexceptions -msse4.2
LFLAGS	= -shared -pthread
INCPATH	= -I../boost_1_62_0 -I../boost.simd/include -I../pybind11/include

SRCDIR = src
OBJDIR = obj

OBJ = soap_kernel_function.o molecule.o solver.o neighbourhood.o structural_similarity.o descriptor.o power_spectrum.o local_similarity.o

default: soap.so

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR) $(SRCDIR)/%.h
	$(CPP) $(CPPFLAGS) $(INCPATH) -c -o $@ $^

$(OBJDIR):
	mkdir $@

soap.so: $(patsubst %, $(OBJDIR)/%, $(OBJ))
	$(CPP) $(CPPFLAGS) $(INCPATH) -o $@ $^ $(LFLAGS)

clean:
	rm -rf $(OBJDIR) *.exe *.so

