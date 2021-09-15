#!/bin/csh

cat >! kyrg.par << EOF
velmod=kyrg
zrange=3,3,54
evla=39.395
evlo=72.812
dt=0.1
nt=2048
fmax=0.2
t0=-1.0
redv=14.
damp=1.
kmax=9999999
eps=0.00001
smin=0.00001
modeldb=/Users/ichinose1/Work/mtinv.v3.0.4/data/modeldb
stadb=rdseed.stations
noverbose
nodump
EOF
                                                                                                                                                  
#mkgrnlib par=kyrg.par stnm=AAK   net=II dt=0.1
#mkgrnlib par=kyrg.par stnm=KKAR  net=KZ dt=0.1
#mkgrnlib par=kyrg.par stnm=MKAR  net=KZ dt=0.25
#mkgrnlib par=kyrg.par stnm=ABKAR net=KZ dt=0.25
#mkgrnlib par=kyrg.par stnm=MAKZ  net=IU dt=0.25
#mkgrnlib par=kyrg.par stnm=EKS2  net=KN dt=0.1
#mkgrnlib par=kyrg.par stnm=USP   net=KN dt=0.1
#mkgrnlib par=kyrg.par stnm=KBK   net=KN dt=0.1
#mkgrnlib par=kyrg.par stnm=ULHL  net=KN dt=0.1
#mkgrnlib par=kyrg.par stnm=UCH   net=KN dt=0.1
#mkgrnlib par=kyrg.par stnm=TKM2  net=KN dt=0.1
#mkgrnlib par=kyrg.par stnm=CHM   net=KN dt=0.1
#mkgrnlib par=kyrg.par stnm=AML   net=KN dt=0.1
                                                                                                                                                  

cat >! mkgrnlib.par << EOF
AAK   II kyrg.par 0.1
KKAR  KZ kyrg.par 0.1
MKAR  KZ kyrg.par 0.25
ABKAR KZ kyrg.par 0.25
MAKZ  IU kyrg.par 0.25
EKS2  KN kyrg.par 0.1
USP   KN kyrg.par 0.1
KBK   KN kyrg.par 0.1
ULHL  KN kyrg.par 0.1
UCH   KN kyrg.par 0.1
AML   KN kyrg.par 0.1
EOF

cat >! mkgrnlib2.par << EOF
TES01 XX kyrg.par 0.1
TES02 XX kyrg.par 0.1
TES03 XX kyrg.par 0.1
TES04 XX kyrg.par 0.1
TES05 XX kyrg.par 0.1
TES06 XX kyrg.par 0.1
TES07 XX kyrg.par 0.1
TES08 XX kyrg.par 0.1
TES09 XX kyrg.par 0.1
TES10 XX kyrg.par 0.1
TES11 XX kyrg.par 0.1
TES12 XX kyrg.par 0.1
TES13 XX kyrg.par 0.1
TES14 XX kyrg.par 0.1
TES15 XX kyrg.par 0.1
TES16 XX kyrg.par 0.1
TES17 XX kyrg.par 0.1
TES18 XX kyrg.par 0.1
TES19 XX kyrg.par 0.1
TES20 XX kyrg.par 0.1
TES21 XX kyrg.par 0.1
TES22 XX kyrg.par 0.1
TES23 XX kyrg.par 0.1
TES24 XX kyrg.par 0.1
TES25 XX kyrg.par 0.1
TES26 XX kyrg.par 0.1
TES27 XX kyrg.par 0.1
TES28 XX kyrg.par 0.1
TES29 XX kyrg.par 0.1
TES30 XX kyrg.par 0.1
EOF

multithread_mkgrnlib parfile=mkgrnlib.par \
                     executable_pathname=/Users/ichinose1/Work/mtinv.v3.0.4/bin/mkgrnlib > test.out
