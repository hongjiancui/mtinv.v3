CC	= cc
INCLUDE = -I../../include
CFLAGS  = -g -m32 $(INCLUDE)
INSTALL = $(HOME)/bin

.c.o: ; $(CC) $(CFLAGS) -c $*.c

all : wiggins

wiggins :  wiggins.o wiggins_sub.o sacio.o
	$(CC) $(CFLAGS) -o $@ wiggins.o wiggins_sub.o sacio.o -lm -L../../lib -lget

clean :
	rm -rf *.o wiggins
