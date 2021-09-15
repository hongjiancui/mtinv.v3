CC = cc
INCLUDE = -I.
#CFLAGS = -m32 $(INCLUDE)
CFLAGS = $(INCLUDE)
LIBS= -lm -L../../lib -lget
# INSTALL = $(HOME)/bin
INSTALL = ../../bin

all : hudson makehudson ellipse

hudson : hudson.o
	$(CC) $(CFLAGS) -o $@ hudson.o $(LIBS)
	mv $@ $(INSTALL)

makehudson : makehudson.o
	$(CC) $(CFLAGS) -o $@ makehudson.o $(LIBS)
	mv $@ $(INSTALL)

ellipse : ellipse.o
	$(CC) $(CFLAGS) -o $@ ellipse.o -lm $(LIBS)
	mv $@ $(INSTALL)

clean : 
	rm -rf *.o hudson makehudson ellipse

