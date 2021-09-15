#
# Top Level Makefile for mtinv package
#

# removed misc/rdseedv5.3.1  - compile seperately when needed

all :
	for dir in misc/Cgraph misc/getpar; do cd $$dir ; make all ; cd ../../ ; done
	for dir in src ; do cd $$dir ; make all ; cd .. ; done
	for dir in misc/sacqc misc/Hudson_Plot misc/renamesac misc/sac2xy misc/sacmerge misc/sacswapbytes misc/multithread_makeglib ; do cd $$dir ; make all ; cd ../../ ; done
clean :
	for dir in misc/Cgraph misc/getpar; do cd $$dir ; make clean ; cd ../.. ; done
	for dir in src ; do cd $$dir ; make clean ; cd .. ; done
	for dir in misc/sacqc misc/Hudson_Plot misc/renamesac misc/sac2xy misc/sacmerge misc/sacswapbytes misc/multithread_makeglib ; do cd $$dir ; make clean ; cd ../../ ; done
