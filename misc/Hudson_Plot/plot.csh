#!/bin/csh

makehudson R="-R-1.5/1.5/-1.5/1.5" J="-JX6i" X="-X1i" Y="-Y1i" W="-W1.8p/0" init orientation="-P" ps=plot.ps 

#ellipse hudson x0=0.0 y0=0.0 theta=350 a=0.2 b=0.05 | psxy -R -JX -M -W1p/blue -G200/200/255 -O -K >> plot.ps
ellipse       x0=0.0 y0=0.0 theta=350 a=0.2 b=0.05 | psxy -R -JX -M -W1p/red -O -K >> plot.ps

#ellipse hudson x0=0.6 y0=0.6 theta=135 a=0.2 b=0.1 | psxy -R -JX -M -W1p/blue -G200/200/255 -O -K >> plot.ps
ellipse        x0=0.6 y0=0.6 theta=135 a=0.2 b=0.1 | psxy -R -JX -M -W1p/red -O -K >> plot.ps

#ellipse hudson x0=-0.6 y0=0.6 theta=45 a=0.2 b=0.1 | psxy -R -JX -M -W1p/blue -G200/200/255 -O -K >> plot.ps
ellipse        x0=-0.3 y0=0.6 theta=45 a=0.2 b=0.1 | psxy -R -JX -M -W1p/red -O -K >> plot.ps

#ellipse hudson x0=0.3 y0=-0.7 theta=85 a=0.2 b=0.1 | psxy -R -JX -M -W1p/blue -G200/200/255 -O -K >> plot.ps
ellipse   x0=0.3 y0=-0.7 theta=85 a=0.2 b=0.1 | psxy -R -JX -M -W1p/red -O -K >> plot.ps

echo " 0.0 0.0 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps

echo " -0.1 0.1 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " -0.2 0.2 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " -0.3 0.3 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " -0.4 0.4 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " -0.5 0.5 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " -0.6 0.6 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " -0.7 0.7 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " -0.8 0.8 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " -0.9 0.9 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " -0.99 0.99 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps

echo " +0.1 -0.1 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " +0.2 -0.2 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " +0.3 -0.3 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " +0.4 -0.4 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " +0.5 -0.5 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " +0.6 -0.6 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " +0.7 -0.7 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " +0.8 -0.8 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " +0.9 -0.9 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps
echo " +0.99 -0.99 " | hudson | psxy -R -JX -Sc0.1i -W1p/black -Gred -O -K >> plot.ps

#gs plot.ps
ps2pdf plot.ps
open plot.pdf
