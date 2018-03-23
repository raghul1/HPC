CC = g++
CFLAGS = -std=c++11 -Wall -O2
HDRS = callfunctions.h
OBJS = main.o vtkwrite.o zeros_int.o zeros_double.o
LDLIBS = -llapack -lblas

default: c1
all: c1 c2 c3

%.o : %.cpp $(HDRS)
	$(CC) $(CFLAGS) -o $@ -c $<

compile: $(OBJS)
	$(CC) -o $@ $^ $(LDLIBS)

	# case_no, a, h1, h2, L, th, nelem_x, nelem_y, T0, q0, kx, ky, kxy
c1: compile
	./compile 1 0 1 1 2 0.2 10 5  10 2500 250 250 0

c2: compile
	./compile 2 0 1 1 2 0.2 10 5 3 1 10 2500 250 250 0

c3: compile
	./compile 3 0.25 1 1.3 3 0.2 15 8 -20 -5000 250 250 0

clean:
	rm -rf *o c1 c2 c3 compile
