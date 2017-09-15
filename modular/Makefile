SRC = $(wildcard *.c)
OBJS = $(SRC:.c=.o)
DEPS = variables.h
CXX = gcc
DEBUG = -g
CXXFLAGS = -Wall -c $(DEBUG)
LFLAGS = $(DEBUG) -O2 -Wall 

GSD : $(OBJS)
	$(CXX) -o GSD $(OBJS) $(LFLAGS) -lm

gsd.o : $(DEPS) gsd.h  gsd_fn.c
	$(CXX) $(CXXFLAGS) gsd_fn.c

analyze.o : $(DEPS) analyze.h  analyze.c 
	$(CXX) $(CXXFLAGS) analyze.c -lm

gsd_read.o : $(DEPS) gsd.h gsd_read.h  gsd_read.c
	$(CXX) $(CXXFLAGS) gsd_read.c -lm

main.o : $(DEPS) gsd.h gsd_read.h analyze.h main.c
	$(CXX) $(CXXFLAGS) main.c -lm

clean:
	\rm *.o *~ GSD

