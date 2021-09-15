CC = cc

INCLUDE = -I../../include 

CFLAGS = -m32 $(INCLUDE)

LIBS = -lm -lget

all : driver 

driver : driver.o rtrend.o ../sacio.o ../timesubs.o ../math.o ../rmean.o
	$(CC) $(CFLAGS) -o $@ driver.o rtrend.o ../sacio.o ../timesubs.o ../math.o ../rmean.o $(LIBS)

clean :
	rm -f *.o

