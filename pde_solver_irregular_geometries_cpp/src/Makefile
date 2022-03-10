PROGRAM		= ../bin/solver
CXX		= g++
CXXFLAGS	= -O3 -Wall -fopenmp

IDIR	= ../include
ODIR	= ../obj
LIBS	= -lgsl -lcblas

_DEPS	= solver.h
DEPS	= $(patsubst %,$(IDIR)/%,$(_DEPS))

CPPFILES	:= $(wildcard *.cpp)
OBJ		:= $(patsubst %.cpp,$(ODIR)/%.o,$(CPPFILES))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS) -I$(IDIR)

$(PROGRAM): $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) -I$(IDIR) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
	rm -f $(PROGRAM)
	rm -f *.txt
