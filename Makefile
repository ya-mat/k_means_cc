LINKER	     = g++
CC	     = g++
#LDFLAGS	     = -L/usr/lib/ -lblas -llapack
LDFLAGS	     = 
COPTS	     = -O3 -std=c++17 -Wall --pedantic-errors
INCLUDE		= -I/home/fafa/lib/eigen

OBJS          = main.o\

PROGRAM	      = a.out

all:		$(PROGRAM)

$(PROGRAM): $(OBJS)
		$(LINKER) $(COPTS) $(OBJS) -o $(PROGRAM) $(LDFLAGS)

clean:
		rm -f $(PROGRAM) *.o *~ ;\

.SUFFIXES: .o .cc

.cc.o :
		$(CC) $(COPTS) $(INCLUDE) -c -o $*.o $*.cc

