c
c takes values of tk as i/p and writes x,y coords for source type plot
c of Hudson et al. 1989
c
      real t,k,x,y
c
      write(0,*) 'Enter values of T and k'
      read(*,*) t,k
c  protection
      if (abs(t).gt.1.0.or.abs(k).gt.1.0) then
         write(0,*) 'T or k out of range try again'
         read(*,*) t,k
      endif
c  crashes if k = +/- 1 and t = +/-1
      if(abs(k).eq.1.0.and.abs(t).eq.1.0) t=0.0
c
      call centre(t,k,x,y)
      write(*,99) x,y
 99   format(1x,f9.6,2x,f9.6)
      stop
      end

      SUBROUTINE CENTRE (T,K,XCEN,YCEN)
C
      IMPLICIT DOUBLE PRECISION (A-H)
      REAL K
C
      DT=T
      DK=K
C
      IF ((DK .GE. 0.D0 .AND. DT .LE. 0.D0) .OR.
     +    (DK .LE. 0.D0 .AND. DT .GE. 0.D0)) THEN
       DX=DT*(1.D0-DABS(DK))
       DY=DK
      ELSE
       DX = DSIGN(1.D0,DK)/(1/(DABS(DT)*(1.D0-DABS(DK)))-0.5D0)
       DY = DK*(1.D0+(DABS(DX)/2.D0))
       IF ((DY/DX) .LT. 0.25D0)THEN
        IF (DK .NE. 0.D0) THEN
         DY = (1.D0/(1.D0/DK -(2.D0*DSIGN(1.D0,DK))))
        ELSE
         DY=0.D0
        ENDIF
        DX = DT*(1.D0+DABS(DY))
       ENDIF
      ENDIF
C
      XCEN=DX
      YCEN=DY
C
      RETURN
C
      END

