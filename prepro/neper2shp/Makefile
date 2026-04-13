APP = neper2shp
CC =   g++
COMM_PATH = ../../common
ROCK_PATH = ../../src
TCLA_PATH = ../../tclap
OPTIONS = -DnQUAT_ACC -DnFT_CORR -DnROT_MATRIX -DnBREAK_ONCE

CFLAGS = -O3 -Wall -Wextra -std=c++0x -I $(COMM_PATH) -I $(TCLA_PATH) -I $(ROCK_PATH) $(OPTIONS)

OBJS = $(APP).o tessToshp.o $(LIB_PATH)/*.o

all: $(APP)

%.o:%.cpp
	$(CC) $(CFLAGS) -c $<

Shape.o: $(ROCK_PATH)/Rockable.cpp $(ROCK_PATH)/Rockable.hpp
	$(CXX) $(CXXFLAGS) -c $<

$(APP) : $(OBJS)
	$(CC) -o $@ $(CFLAGS) $(OBJS)

clean:
	rm -f *~ *.o $(APP)

depend:
	makedepend -- $(CFLAGS) -- *.cpp
# DO NOT DELETE

