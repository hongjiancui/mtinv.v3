#!/bin/csh

rm -f x.sac x.sac.int xint.sac

sac << EOF
fg random 1 1 npts 2048 delta 0.01
write x.sac
read x.sac
lh npts delta
interpolate delta 0.2
write append .int
quit
EOF

sacsun2linux x.sac
wiggins f=x.sac nt=103 dt=0.2

cat >! sac.mac << EOF
read xint.sac x.sac.int
qdp off
color on inc on
lh npts delta
p1
EOF

sac sac.mac
