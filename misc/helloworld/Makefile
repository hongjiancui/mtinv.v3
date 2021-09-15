CC = cc

#CFLAGS = -m32
CFLAGS = 

#LIBS = -L/Users/ichinose/Work/mtinv.v3.0.3/lib -lget
LIBS = -L../../lib -lget

.c.o :
	$(CC) $(CFLAGS) -c $< -o $@

all : helloworld 

helloworld.o : helloworld.c 

helloworld : helloworld.o
	cc $(CFLAGS) -o $@ helloworld.o $(LIBS)
	mv $@ ../../bin

clean : 
	rm -f  helloworld.o 
