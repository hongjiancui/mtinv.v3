CC = gcc
#CC = clang

INCLUDE = -I../include -I.
#INCLUDE = -I../include -I/Users/ichinose/Work/libresponse/include
#INCLUDE = -I../include -I/home/gichinos/libresponse/include
#INCLUDE = -I../include -I/home/gichinose/libresponse/include
CFLAGS =  $(INCLUDE)
LIBDIRS = -L../../../lib

#LIBDIRS = -L/Users/ichinose/Work/libresponse/lib
#LIBDIRS = -L/home/gichinos/libresponse/lib
#LIBDIRS = -L/home/gichinose/libresponse/lib

LIBS = -lget -lm
#INSTALL_BIN = $(HOME)/bin
INSTALL_BIN = ../../../bin

.c.o :
	$(CC) $(CFLAGS) -c $< -o $@

all : sacqc hist switch_color

hist : hist.o histsub.o
	$(CC) $(CFLAGS) -o $@ hist.o histsub.o $(LIBDIRS) $(LIBS)
	mv $@ $(INSTALL_BIN)

OBJS = sacqc.o sacio.o timesubs.o remove_mean.o histsub.o shorten_path.o stats_subs.o rtrend.o

HDRS = ../include/sacfile.h ../include/sac.h ../include/mytime.h ../include/hist.h

sacqc : $(OBJS) $(HDRS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBDIRS) $(LIBS)
	mv $@ $(INSTALL_BIN)

switch_color : switch_color.o
	$(CC) $(CFLAGS) -o $@ switch_color.o 
	mv $@ $(INSTALL_BIN)

clean :
	rm -rf *.o sacqc
