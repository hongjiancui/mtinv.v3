CC = cc
INCLUDE = -I../../../include
CFLAGS = $(INCLUDE)
LIBDIRS = -L../../../lib
LIBS = -lglib -lm
INSTALL_BIN = ../../../bin

all : sacmerge

sacmerge : sacmerge.o
	$(CC) $(CFLAGS) -o $@ sacmerge.o $(LIBDIRS) $(LIBS)
	mv $@ $(INSTALL_BIN)

clean : 
	rm -f *.o sacmerge
