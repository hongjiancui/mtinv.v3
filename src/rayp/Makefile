CC = cc
FC = g77
INCLUDE = -I../include
GRN_LIBS = -lm -L../lib -lglib -lget
PLOT_LIBS = -L../lib -lcgraph
INSTALL = ../bin
CFLAGS = -g -m32 $(INCLUDE)
FFLAGS = -g -m32

.c.o :
	$(CC) $(CFLAGS) -c $< -o $@

.f.o :
	$(FC) $(FFLAGS) -c $< -o $@

rayp :	rayp.o rayp_subs.o ../stasubs.o ../modsubs.o ../misc_tools/distaz.o
	cc -m32 -o rayp rayp.o rayp_subs.o ../stasubs.o ../modsubs.o ../misc_tools/distaz.o -lm -L../../lib -lglib -lget

clean : 
	rm -f *.o rayp
