CC = cc
INCLUDE = -I.
INSTALL = ../bin
CFLAGS = -g -m32 $(INCLUDE)

.c.o :
	$(CC) $(CFLAGS) -c $< -o $@

all : envelope_driv

envelope_driv.o : envelope_driv.c 
envelope_sub.o : envelope_sub.c complex.h sac.h 
fftsub.o : fftsub.c complex.h
sacio.o : sacio.c sac.h

SUBS = envelope_driv.o envelope_sub.o fftsub.o sacio.o

envelope_driv : $(SUBS)
	$(CC) $(CFLAGS) -o $@ $(SUBS) -lm -L../../lib -lget


clean : 
	rm -f *.o envelope_driv
