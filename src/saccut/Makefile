CC	= cc
FC	= f77
INCLUDE	= -I../../include
#CFLAGS	= -mcpu=7450 -pipe -maltivec -mabi=altivec -mpowerpc-gfxopt $(INCLUDE) -O 
#FFLAGS	= -mcpu=7450 -pipe -maltivec -mabi=altivec -mpowerpc-gfxopt $(INCLUDE) -ff90
CFLAGS  =
FFLAGS  =
LIBS	= -lm -lget
LIBDIR	=
INSTALL	= $(HOME)/bin

EXEC = saccut

all : $(EXEC)

OBJS = saccut.o saccut_subs.o ../timesubs.o ../sacextrema/sacextrema.o
saccut : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS) $(LIBDIR)
	mv $@ $(INSTALL)

clean :
	rm -f *.o $(EXEC)
