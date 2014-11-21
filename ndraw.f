!----------------------------------------------------------------------------
            subroutine draw0(xd,yd,n,xs,ns)
            real,dimension (n)::xd,yd
            real,dimension (ns)::xs

            pi2 = 8.*atan(1.)
            dx = pi2/float(n)
            do i=1,n
               xd(i) = float(i-1)*dx
               yd(i) = xd(i)
            enddo
            do i=1,ns
               xs(i) = float(i)
            enddo

            return
            end 
!----------------------------------------------------------------------------
            subroutine draw(xd,yd,vor,n,xs,ys,ns)
            real xd(n),yd(n),vor(n,n),xs(ns),ys(ns)

            return
            end
!------------------------------------
         subroutine plotp2d(xx,yy,dx)
         
         return
         end
!------------------------------------
         subroutine newpen(nn)
         
         return
         end
!--------------------------------------------------------------
         subroutine plot_vortex(n,iv,jv,nc,ic)
         dimension iv(nc),jv(nc)

         return
         end
!*************************************************
         SUBROUTINE MOVESPA(FIN,NN)
         CHARACTER*(*) FIN
         
!.. 1 FIND LAST NO SPACE CHARACTER NUMBER
         KL=0
         DO I=NN,1,-1
          IF(FIN(I:I).NE.' ') THEN
            KL=I
            GOTO 11
          ENDIF
         ENDDO
11       CONTINUE
         IF(KL.LE.1) RETURN

!.. 2 DELETE SPACE FROM NN1 TO 1
         NN1=KL-1
         DO 2 I=NN1,1,-1
          IF(FIN(I:I).EQ.' ') THEN
           DO 3 J=I,NN1
            JP1=J+1
            FIN(J:J)=FIN(JP1:JP1)
 3         CONTINUE
           FIN(KL:KL)=' '
          ENDIF            
 2       CONTINUE

         RETURN
         END
! ------------
         SUBROUTINE ANG_ABC(XA,YA,XB,YB,XC,YC,ANGLE)
        
         PARAMETER(DPI=2.0*3.14159265)
         PARAMETER(U1=-1.+1.E-10,U2=1.-1.E-10,Z0=1.E-10)
         C=SQRT((XA-XB)**2+(YA-YB)**2)
         B2=(XA-XC)**2+(YA-YC)**2
         A=SQRT((XC-XB)**2+(YC-YB)**2)
         IF(A.LT.Z0 .OR. C.LT. Z0) THEN
          ANGLE = 0.0
          RETURN 
         ENDIF 
         COSABC=(A*A+C*C-B2)/(2.*A*C)
         IF(COSABC.LE.U1) THEN
            ANGLE=3.14159265
         ELSEIF(COSABC.GE.U2) THEN
            ANGLE=0.0
         ELSE    
            ANGLE=ACOS(COSABC)
         ENDIF
         CALL AREAM(XA,YA,XB,YB,XC,YC,S)
         IF(S.LT.Z0) ANGLE=DPI-ANGLE
         RETURN
         END
! -------------------------
!...AREA COORDINATES AT M
           SUBROUTINE AREAM(X1,Y1,X2,Y2,XM,YM,AR3)
           AR3=((X1-XM)*Y2-(Y1-YM)*X2
     *     -(X1-XM)*YM+(Y1-YM)*XM)*0.5
           RETURN
           END
! -------------------------
!...AREA COORDINATES AT M
           SUBROUTINE run_time_msec(ii)
           ier = ii
           RETURN
           END
