CC = gcc
INCLUDE = ../../include
CFLAGS = -I$(INCLUDE)
LIBS = -L../../lib -lget -lm
INSTALL = ../../bin

all : sac2xy

.c.o :
	$(CC) $(CFLAGS) -c $< -o $@

sac2xy : sac2xy.o
	$(CC) $(CFLAGS) -o sac2xy sac2xy.o $(LIBS)
	mv $@ $(INSTALL)

clean :
	rm -f sac2xy.o sac2xy
