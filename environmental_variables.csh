#!/bin/csh -x

set workingdir = `pwd`

append_path PATH ${workingdir}/bin
append_path MANPATH ${workingdir}/man

setenv MTINV_GMT_CPT_FILE  MYTOPO.cpt
setenv MTINV_GMT_GRID_FILE  /usr/local/gmt/share/dbase/etopo5.grd
setenv MTINV_GMT_INT_FILE  /usr/local/gmt/share/dbase/etopo5.int
setenv MTINV_PATH /Users/ichinose/mtinv.v3.0.5
setenv MT_DATABASE_FILE /Users/ichinose/mtinv.v3.0.5/mt.db
