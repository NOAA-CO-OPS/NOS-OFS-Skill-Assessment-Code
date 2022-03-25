       DIMENSION SIGMA(11),ZSIGMA(11)
!       data sigma/-0.05,-0.15,-0.25,-0.35,-0.45,-0.55,-0.65,
!     1  -0.75,-0.85, -0.95/
!       data sigma/0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85, 0.95/

       KB = 11
       H = 50.0
       ELE = 0.0
       THETAS = 4.5
       THETAB = 0.95
       TC = 10.0
       HC = 6.0
       HMIN = 6.0
       P_SIGMA = 2

!--------  SET SIGMA LEVELS  FVCOM--------------------------------------------------!  
      IF(P_SIGMA == 1) THEN
         DO K = 1, KB
           SIGMA(K) = -((K-1)/FLOAT(KB-1))**P_SIGMA 
         END DO
      ELSE
         DO K = 1,(KB+1)/2
           SIGMA(K) = -((K-1)/FLOAT((KB+1)/2-1))**P_SIGMA/2 
         END DO
         DO K = (KB+1)/2+1,KB
           SIGMA(K) = ((KB-K)/FLOAT((KB+1)/2-1))**P_SIGMA/2-1.0
         END DO
      END IF

!       CALL SIGMA2Z_ROMS(SIGMA,H,ELE,KB,ZSIGMA)
       CALL SIGMA2Z_POM(SIGMA,H,ELE,KB,ZSIGMA)
       DO K = 1, KB
         PRINT *,K,SIGMA(K),ZSIGMA(K)
       ENDDO

       END


       SUBROUTINE SIGMA2Z_ROMS(SIGMA,H,ELE,KB,ZSIGMA)
       DIMENSION SIGMA(KB),ZSIGMA(KB)
       THETAS = 4.5
       THETAB = 0.95
       TC = 10.0
       HC = 6.0
       HMIN = 6.0
       DO K = 1,KB
!         SIGMA(K)=-1+(K-0.5)/KB
         PTHETA = SINH(THETAS*SIGMA(K))/SINH(THETAS)
         RTHETA = TANH(THETAS*(SIGMA(K)+0.5))/(2.0*TANH(THETAS*0.5))-0.5
         CSIGMA = (1.0-THETAB)*PTHETA+THETAB*RTHETA
         Z0 = HC*SIGMA(K)+(H-HC)*CSIGMA
         ZSIGMA(K) = Z0+(1.+Z0/H)*ELE
         IF (ZSIGMA(K) .LT. 0.0) ZSIGMA(K) = -ZSIGMA(K)
!         PRINT *,K,SIGMA(K),Z
       ENDDO

       RETURN
       END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       SUBROUTINE SIGMA2Z_POM(SIGMA,H,ELE,KB,ZSIGMA)
       DIMENSION SIGMA(KB),ZSIGMA(KB)
       DO K = 1,KB
         ZSIGMA(K) = SIGMA(K)*(H+ELE)+ELE !FOR POM
         IF (ZSIGMA(K) .LT. 0.0) ZSIGMA(K) = -ZSIGMA(K)
       ENDDO

       RETURN
       END
