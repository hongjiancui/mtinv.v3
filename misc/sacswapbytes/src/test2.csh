#!/bin/csh

#             01     02    03   04  05  06  07  08   09   10   11    12     13 
set dt=( 10000.0 1000.0 100.0 10.0 5.0 2.0 1.0 0.1 0.01 0.02 0.05 0.001 0.0001 )

foreach i ( 01 02 03 04 05 06 07 08 09 10 11 12 13 )

echo $i $dt[$i]

sac << EOF
read test.sac.littleendian
ch delta $dt[$i]
write append .$i
quit
EOF
sacsun2linux test.sac.littleendian.$i
end


foreach i ( 01 02 03 04 05 06 07 08 09 10 11 12 13 )

dumpsac if=test.sac.littleendian.$i

end
