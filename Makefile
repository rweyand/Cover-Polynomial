OBJS = main.o tools.o matrix.o poly.o remaindering.o
OBJS1 = small.o tools.o matrix.o poly.o remaindering.o
OBJS2 = fast.o tools.o matrix.o poly.o remaindering.o
CC = cc
CFLAGS = -Wall -lgmp -pg


small: $(OBJS1)
	cc $(OBJS1) -o small -lgmp
	
fast: $(OBJS2)
	cc $(OBJS2) -o fast -lgmp -fopenmp

fastdebug: $(OBJS2)
	 cc $(OBJS2) -pg -o fast -lgmp -fopenmp
	
smalldebug: $(OBJS1)
	cc $(OBJS1) -o small -g -lgmp -O0

cover: $(OBJS)
	cc $(OBJS) -o cover -lgmp

debug: $(OBJS)
	cc $(OBJS) -o cover -g -lgmp -O0

profile: main_profile.o tools.h
	cc -pg main.o -o cover -lgmp 


%.o: %.c
	cc -c  $? -o  $@ -fopenmp

clean:
	rm *.o cover small
