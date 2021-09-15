CC	= gcc
FC 	= g77
INCLUDE = -I../../include
FFLAGS  = -g $(INCLUDE)
CFLAGS  = -g $(INCLUDE)
INSTALL = $(HOME)/bin

.f.o: ; $(FC) $(FFLAGS) -c $*.f
.c.o: ; $(CC) $(CFLAGS) -c $*.c

EXEC 	= interpolate_drv
OBJS    = interpolate_drv.o interpolate_subs.o

all : $(EXEC)

$(EXEC) : $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) -lm -lget
	mv $(EXEC) $(INSTALL)

clean :
	rm -rf *.o $(EXEC)
