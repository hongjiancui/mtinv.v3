CC = gcc
#FC = f77
#CC = clang
FC = gfortran
#
INCLUDE = -I../include
GRN_LIBS = -lm -L../lib -lglib -lget
PLOT_LIBS = -L../lib -lcgraph
INSTALL = ../bin

#CFLAGS = -Wall -g -m32 $(INCLUDE)
#FFLAGS = -m32
#
CFLAGS = $(INCLUDE)
FFLAGS = 

.c.o :
	$(CC) $(CFLAGS) -c $< -o $@
.f.o : 
	$(FC) $(FFLAGS) -c $< -o $@

EXEC = libglib mkgrnlib grnlib2sac glib2inv sacdata2inv mtinv multi_mtinv makepar \
       setupMT pltmod mtsim pea whatshere mtbestfit grn2Mxy

all : $(EXEC)

### now in sacio.o

OBJS0 = rtrend/rtrend.o \
	rmean.o \
	transfer/transfersubs.o transfer/fftsub.o \
	interpolate/interpolate_subs.o \
	filter/filtersubs.o filter/complex.o \
	saccut/saccut_subs.o \
	source/source_subs.o \
	timesubs.o \
	modsubs.o \
	misc_tools/distaz.o misc_tools/ampshift.o misc_tools/fmul.o \
	taper/taper_subs.o \
	rotate/rotate.o \
	math/math.o \
	crosscorrelation/cross_correlation.o \
	tdif/Differentiates.o \
	shift/shift.o \
	random/ransubs.o \
	envelope/envelope_sub.o \
	envelope/fftsub.o \
	wiggins/wiggins_sub.o \
	sacio/sacio.o \
	writesacfile.o

libglib : $(OBJS0)
	ar -rv libglib.a $(OBJS0)
	mv libglib.a ../lib
	ranlib ../lib/libglib.a

OBJS1 = mkgrnlib.o wrtgrn2sac.o stasubs.o greensf.o besasy.o bessel.o cinv2.o fork.o \
	genrefl.o timesubs.o modsubs.o misc_tools/distaz.o rayp/rayp_subs.o shorten_path.o

mkgrnlib : $(OBJS1)
	$(FC) $(FFLAGS) -o $@ $(OBJS1) -lm -L../lib -lget
	mv $@ $(INSTALL)

OBJS2 = grnlib2sac.o wrtgrn2sac.o
grnlib2sac : $(OBJS2)
	$(CC) $(CFLAGS) -o $@ $(OBJS2) -lm -L../lib -lglib -lget
	mv $@ $(INSTALL)

grn2Mxy : grn2Mxy.o timesubs.o
	$(CC) $(CFLAGS) -o $@ grn2Mxy.o timesubs.o -lm -L../lib -lglib -lget
	mv $@ $(INSTALL)

OBJS3 = sacdata2inv.o sacdata2inv_subs.o glib2inv_subs.o getrespfile_sub.o compute_Peak_to_Peak.o \
         shorten_path.o sacdata2inv_serial.o sacdata2inv_parallel.o sacdata2inv_loaddata.o

sacdata2inv : $(OBJS3) libglib
	$(CC) $(CFLAGS) -pthread -o $@ $(OBJS3) -lm -L../lib -lglib -lget -lpthread
	mv $@ $(INSTALL)

OBJS4 = glib2inv.o glib2inv_subs.o glib2inv_serial.o glib2inv_parallel.o wrtgrn2sac.o shorten_path.o
glib2inv : $(OBJS4)
	$(CC) $(CFLAGS) -pthread -o $@ $(OBJS4) -lm -L../lib -lglib -lget -lpthread
	mv $@ $(INSTALL)

MTGS_SUBS = mtgs.o glib2inv_subs.o svbksb.o svdcmp.o tred2.o tqli.o trans1.o tpdss.o \
            mtinv_subs.o mtinv_gmtsubs.o ascii_mech.o psplot.o svdvar.o psmech.o make_amatrix.o indexx.o \
            check_depths.o 

mtgs : $(MTGS_SUBS)
	$(CC) $(CFLAGS)  -o $@ $(MTGS_SUBS) $(GRN_LIBS) $(PLOT_LIBS)
	mv $@ $(INSTALL)

MTINV_SUBS= mtinv.o glib2inv_subs.o svbksb.o svdcmp.o tred2.o tqli.o trans1.o tpdss.o \
	mtinv_subs.o mtinv_gmtsubs.o ascii_mech.o psplot.o svdvar.o psmech.o make_amatrix.o indexx.o \
        check_depths.o dbsubs.o eig2lune.o shorten_path.o find_best_shift.o realtime_subs.o dumpSAC.o

mtinv : $(MTINV_SUBS) $(OBJS0)
	$(CC) $(CFLAGS)  -o $@ $(MTINV_SUBS) $(GRN_LIBS) $(PLOT_LIBS)
	mv $@ $(INSTALL)

multi_mtinv : multi_mtinv.o
	$(CC) $(CFLAGS)  -o $@ multi_mtinv.o $(GRN_LIBS) $(PLOT_LIBS)
	mv $@ $(INSTALL)

MTSIM_SUBS = mtsim.o mtsim_subs.o glib2inv_subs.o svbksb.o svdcmp.o tred2.o tqli.o trans1.o tpdss.o \
        mtinv_subs.o mtinv_gmtsubs.o ascii_mech.o psplot.o svdvar.o psmech.o make_amatrix.o indexx.o \
        check_depths.o dbsubs.o eig2lune.o shorten_path.o find_best_shift.o mtsim_invert_sim_parallel.o \
        dumpSAC.o

mtsim : $(MTSIM_SUBS) 
	$(CC) $(CFLAGS) -pthread -o $@ $(MTSIM_SUBS) $(GRN_LIBS) $(PLOT_LIBS)
	mv $@ $(INSTALL)

SETUPMT_OBJS = setupMT.o indexx.o shorten_path.o sacio/sacio.o timesubs.o misc_tools/distaz.o 
setupMT : $(SETUPMT_OBJS)
	$(CC) $(CFLAGS) -o $@ $(SETUPMT_OBJS) -lm -L../lib -lget
	mv $@ $(INSTALL)

makepar : makepar.o indexx.o shorten_path.o
	$(CC) $(CFLAGS) -o $@ makepar.o indexx.o -lm -L../lib -lget
	mv $@ $(INSTALL)

pltmod : pltmod.o modsubs.o
	$(CC) $(CFLAGS) -o $@ pltmod.o modsubs.o -lm -L../lib -lget
	mv $@ $(INSTALL)

pea  : pea.o tred2.o tqli.o eigsrt.o hudson_subs.o
	$(CC) $(CFLAGS) -o $@ pea.o tred2.o tqli.o eigsrt.o hudson_subs.o -lm -L../lib -lget $(GRN_LIBS)
	mv $@ $(INSTALL)

whatshere : whatshere.o
	$(CC) $(CFLAGS) -o $@ whatshere.o -lm -L../lib -lget $(GRN_LIBS)
	mv $@ $(INSTALL)

testc2 : testc2.o math/math.o
	$(CC) $(CFLAGS) -o $@ testc2.o math/math.o -L../lib -lget
	mv $@ $(INSTALL)

mtbestfit : mtbestfit.o timesubs.o 
	$(CC) $(CFLAGS) -o $@ mtbestfit.o timesubs.o -L../lib -lget -lm
	mv $@ $(INSTALL)

clean :
	rm -f *.o */*.o
