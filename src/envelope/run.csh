#!/bin/csh

rm -f env.sac

sac << EOF
fg sine 1 0 npts 1024 delta 0.01
taper type cosine w 0.5
taper type hann w 0.25
write sine.sac
read sine.sac
envelope
write sine.sac.env
quit
EOF

sacsun2linux sine.sac
envelope_driv sacf=sine.sac


cat >! sac.mac << EOF
window x 0.1 0.5 y 0.1 0.9
bd x
qdp off
read sine.sac.env env.sac
color on inc on
p2
pause
quit
EOF

sac sac.mac
