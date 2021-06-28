CC=gcc
CFLAGS=-g -Wall -I
LIBS=-lm
DEPS= gp.h gploc.h

gp:	gp.o gploc.o $(DEPS)
	$(CC) -g -o gp gp.o gploc.o $(LIBS)

gp.o: gp.c $(DEPS)
	$(CC) -g -c gp.c

gploc.o: gploc.c $(DEPS)
	$(CC)  -g -c gploc.c

clean:
	rm -rf gp.o gploc.o
