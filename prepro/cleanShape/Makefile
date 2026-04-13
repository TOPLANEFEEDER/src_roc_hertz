# Paths
COMMONPATH = ../../common
TCLAPPATH = ../../tclap
SRCPATH = ../../src

UNAME_S := $(shell uname -s)

# The compiler to be used
ifeq ($(UNAME_S),Darwin)
  CXX = clang++
  CXXFLAGS = -O3 -Wall -std=c++11 -I $(COMMONPATH) -I $(TCLAPPATH) -I $(SRCPATH)
else
  CXX = g++
  CXXFLAGS = -O3 -Wall -std=c++0x -I $(COMMONPATH) -I $(TCLAPPATH) -I $(SRCPATH)
endif

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<

all: cleanShape

clean:
	rm -f *.o cleanShape

Shape.o: $(SRCPATH)/Shape.cpp $(SRCPATH)/Shape.hpp
	$(CXX) $(CXXFLAGS) -c $<

cleanShape: cleanShape.cpp Shape.o
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(CXXFLAGS) -o $@ $@.o Shape.o
