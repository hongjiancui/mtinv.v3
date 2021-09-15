CC = cc

#CFLAGS = -m32
CFLAGS = 

#LIBS = -L/Users/ichinose/Work/mtinv.v3.0.3/lib -lget
LIBS = -L../../lib -lget

.c.o :
	$(CC) $(CFLAGS) -c $< -o $@

all : multithread_mkgrnlib 

multithread_mkgrnlib.o : multithread_mkgrnlib.c 

multithread_mkgrnlib : multithread_mkgrnlib.o
	cc $(CFLAGS) -o $@ multithread_mkgrnlib.o $(LIBS)
	mv $@ ../../bin

clean : 
	rm -f  multithread_mkgrnlib.o 
