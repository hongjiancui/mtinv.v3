CC = clang
CFLAGS = -I../../include

LIBS = -lget -lm -lglib
LIBDIR = -L../../lib

all : mtdecomp

SUBS = ../mtinv_subs.o ../tred2.o ../tqli.o ../trans1.o ../tpdss.o \
       ../mtinv_gmtsubs.o ../indexx.o ../eig2lune.o ../ascii_mech.o 

mtdecomp : mtdecomp.o $(SUBS)
	$(CC) $(CFLAGS) -o $@ mtdecomp.o $(SUBS) $(LIBDIR) $(LIBS)

clean :
	/bin/rm -f *.o mtdecomp
