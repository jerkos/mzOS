progs = emass

objs = emass.o parser.o formula.o getopt.o
headers = emass.h parser.h formula.h getopt.h

CXX = g++

libs = -lm
#CXXFLAGS = -Wall -W -ansi -pedantic -g -pg
CXXFLAGS = -O3

LDFLAGS =

emass : $(objs) $(headers)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(objs) -o $@ $(libs)

emass.o : emass.cpp $(headers)
	$(CXX) $(CXXFLAGS) -c $<

parser.o : parser.cpp $(headers)
	$(CXX) $(CXXFLAGS) -c $<

formula.o : formula.cpp $(headers)
	$(CXX) $(CXXFLAGS) -c $<

getopt.o : getopt.cpp $(headers)
	$(CXX) $(CXXFLAGS) -c $<


.PHONY: all
all: $(progs)

.PHONY: clean
clean:
	rm -f core *.o $(progs)

.PHONY: bare
bare: clean
	rm -f *~
