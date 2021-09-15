#!/bin/csh

#### example mixed
./mtdecomp verbose Mo=1.0E+23 str=45 dip=45 rak=90 pdc=0.1 clvd_type=+v pclvd=0.2 piso=0.7

#### pure dc
#./mtdecomp verbose Mo=1.0E+23 str=45 dip=45 rak=90 pdc=1.0 clvd_type=+v pclvd=0.0 piso=0.0

### pure clvd
#./mtdecomp verbose Mo=1.0E+23 str=45 dip=45 rak=90 pdc=0.0 clvd_type=+v pclvd=1.0 piso=0.0

### pure ex
#./mtdecomp verbose Mo=1.0E+23 str=45 dip=45 rak=90 pdc=0.0 clvd_type=+v pclvd=0.0 piso=1.0

