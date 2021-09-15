CC=gcc
INSTALL=../../../bin
#CFLAGS = -m32
CFLAGS = 

all : sacswapbytes

sacswapbytes : sacswapbytes.o
	$(CC)  $(CFLAGS) -o $@ sacswapbytes.o
	mv $@ $(INSTALL)

clean :
	rm -f *.o
