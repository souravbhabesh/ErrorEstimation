SRC = $(wildcard *.cpp)
OBJS = $(SRC:.cpp=.o)
DEPS = variables.h
CXX = g++
DEBUG = -g
CXXFLAGS = -Wall -c $(DEBUG) -std=c++11 
LFLAGS = $(DEBUG) -O2 -Wall 

$JKHOOMDTHERMAL : $(OBJS)
	$(CXX) -o JKHOOMDTHERMAL $(OBJS) $(LFLAGS)

jack_knife.o : $(DEPS) jack_knife.h  jack_knife.cpp 
	$(CXX) $(CXXFLAGS) jack_knife.cpp

main.o : $(DEPS) jack_knife.h main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

clean:
	\rm *.o *~ JKHOOMDTHERMAL

