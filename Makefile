SRC = $(wildcard *.cpp)
OBJS = $(SRC:.cpp=.o)
DEPS = variables.h
CXX = g++
DEBUG = -g
CXXFLAGS = -Wall -c $(DEBUG)
LFLAGS = $(DEBUG) -O2 -Wall 

GSD : $(OBJS)
	$(CXX) -o GSD $(OBJS) $(LFLAGS)

analyze.o : $(DEPS) analyze.h  analyze.cpp 
	$(CXX) $(CXXFLAGS) analyze.cpp

main.o : $(DEPS) gsd_fn.h gsd_tools.h analyze.h main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

clean:
	\rm *.o *~ GSD


