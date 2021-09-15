CC = cc
CFLAGS = -m32

all :  driver


OBJS1 = sactransfer.o transfersubs.o fftsub.o

driver : driver.o fftsub.o
	cc -m32 -o driver  driver.o fftsub.o -lm -L/Users/ichinose/Work/mtinv.v3.0/lib -lget

sactransfer : $(OBJS1)
	cc m32 -o $@ $(OBJS1) -lm -lget

clean : 
	rm -f *.o main
