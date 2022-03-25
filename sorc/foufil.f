C     The subroutine FOUFIL_OLD does not work if the data points is 88559,88558,88555,...
C     the filtered output is not correct, the subroutine FOURT does not work.
C     so the data series is splited into different segments each one is 4096 (2^12)
C     Aijun Zhang, Nov. 18, 2004     
       SUBROUTINE FOUFIL(NMAX,DELMIN,TCUT,U,AU)
       dimension u(NMAX),au(NMAX),xtmp(NMAX),ytmp(NMAX)
       NTWO=4096
       NTCUT=2*INT(TCUT*60/DELMIN)
       NF=INT(NMAX/NTWO)
       IF (NF .LT. 1) THEN
         CALL FOUFIL_OLD(NMAX,DELMIN,TCUT,U,AU)
         GOTO 777
       ENDIF
       DO K=1,NF
         JB=(K-1)*NTWO-NTCUT+1
         JE=K*NTWO+NTCUT
         J0=0
         IF (K .EQ. 1) THEN
           DO J=1,MIN0(NMAX,JE)
             J0=J0+1
             XTMP(J0)=U(J)
           ENDDO
           CALL FOUFIL_OLD(J0,DELMIN,TCUT,XTMP,YTMP)
           DO J=1,NTWO
             AU(J)=YTMP(J)
           ENDDO
         ELSE IF (K .GT.1 ) THEN
           DO J=MAX0(1,JB),MIN0(NMAX,JE)
             J0=J0+1
             XTMP(J0)=U(J)
           ENDDO
           CALL FOUFIL_OLD(J0,DELMIN,TCUT,XTMP,YTMP)
           JB=(K-1)*NTWO+1
           JE=K*NTWO
           J0=1
           DO J=JB,JE
            AU(J)=YTMP(J0+NTCUT)
            J0=J0+1
           ENDDO
         ENDIF
       ENDDO
       NRES=NMAX-NF*NTWO
       IF (NRES .GT. 40) THEN
         JB=NF*NTWO-NTCUT+1
         JE=NMAX
         J0=0
         DO J=JB,JE
            J0=J0+1
            XTMP(J0)=U(J)
         ENDDO
         CALL FOUFIL_OLD(J0,DELMIN,TCUT,XTMP,YTMP)
         JB=NF*NTWO+1
         JE=NMAX
         J0=1
         DO J=JB,JE
           AU(J)=YTMP(J0+NTCUT)
           J0=J0+1
         ENDDO
       ELSE
         JB=NF*NTWO+1
         JE=NMAX
         DO J=JB,JE
           AU(J)=U(J)
         ENDDO
       ENDIF
777    CONTINUE      

       RETURN
       END


       SUBROUTINE FOUFIL0(IMX,NMAX,DELMIN,TCUT,U,AU)
       DIMENSION U(IMX),AU(IMX),XTMP(NMAX),YTMP(NMAX)
       DIMENSION UTMP(NMAX),AUTMP(NMAX)

       NTWO = 4096
       NTCUT = 2*INT(TCUT*60/DELMIN)
       NF = INT(NMAX/NTWO)
       IF(NF .LT. 1) THEN
         DO K = 1, NMAX
           UTMP(K) = U(K)
         ENDDO
         CALL FOUFIL_OLD(NMAX,DELMIN,TCUT,UTMP,AUTMP)
         DO K = 1, NMAX
           AU(K) = AUTMP(K)
         END DO
         RETURN
       ENDIF

       DO K = 1, NF
         JB = (K-1)*NTWO-NTCUT + 1
         JE = K*NTWO+NTCUT
         J0 = 0
         IF(K .EQ. 1) THEN
           DO J = 1, MIN0(NMAX,JE)
             J0 = J0 + 1
             XTMP(J0) = U(J)
           ENDDO
           CALL FOUFIL_OLD(J0,DELMIN,TCUT,XTMP,YTMP)
           DO J = 1, NTWO
             AU(J) = YTMP(J)
           ENDDO
         ELSE IF(K .GT. 1) THEN
           DO J = MAX0(1,JB), MIN0(NMAX,JE)
             J0 = J0 + 1
             XTMP(J0) = U(J)
           ENDDO
           CALL FOUFIL_OLD(J0,DELMIN,TCUT,XTMP,YTMP)
           JB = (K-1)*NTWO + 1
           JE = K * NTWO
           J0=1
           DO J = JB, JE
            AU(J) = YTMP(J0+NTCUT)
            J0 = J0 + 1
           ENDDO
         ENDIF
       ENDDO

       NRES = NMAX - NF*NTWO
       IF(NRES .GT. 40) THEN
         JB = NF*NTWO - NTCUT + 1
         JE = NMAX
         J0 = 0
         DO J = JB, JE
            J0 = J0+1
            XTMP(J0) = U(J)
         ENDDO
         CALL FOUFIL_OLD(J0,DELMIN,TCUT,XTMP,YTMP)
         JB = NF*NTWO + 1
         JE = NMAX
         J0 = 1
         DO J = JB, JE
           AU(J) = YTMP(J0+NTCUT)
           J0 = J0 + 1
         ENDDO
       ELSE
         JB = NF*NTWO + 1
         JE = NMAX
         DO J = JB, JE
           AU(J) = U(J)
         ENDDO
       ENDIF

       RETURN
       END


C  This program is used as low-pass fourier filter with IOPT=1 and ITYPE=3 
C-----------------------------------------------------------------------
C
      SUBROUTINE FOUFIL_OLD(length,DELMIN,TCUT,U,AU)
c       length = number of points 
c       delmin = interval between points (minutes)
c       tcut = cutoff period (hrs)
c         u = unfiltered series
c        au = filtered series
C       FOURIER FILTER

      DIMENSION U(LENGTH),AU(LENGTH)
      REAL AMEAN,UMEAN,LMEAN,ASD,USD,LSD

      LENGTH0 = LENGTH
      NRECS = LENGTH
      IUNIT = 10
      IFILTER = 1
      ISORT = 1
      ITYPE = 3 
      IOPT = 1             ! 1=lo pass, 2=hi pass, 3=band pass
      NDATAS = 0
      CPERIOD = TCUT*60.0
      RADPDEG = ACOS(-1.0)/180.0
      DEGPRAD = 1.0/RADPDEG
c
C   fname  = input file name 
C   ndatas = data records to be skipped
c   length = length of file
C   isort  = output intervals vs input interval 
C   IFILTER = 1 TO FILTER
c   itype  = 3 for Fourier filter
c   iopt   = 1 for low pass
c   cperiod = TCUT * 60 (min)
C   path4  = output file name 
C   path4  = log file name 
C
C=====================================================
C NSPH0 is the number of sampling per hour (use CEOB's 
C     default setting nsph0=6
      NSPH0 = NINT(60.0/DELMIN)
C=====================================================
      CALL CMSV(U,LENGTH,AMEAN,UMEAN,LMEAN,ASD,USD,LSD)
      IF(LENGTH .LE. 10) THEN
        DO I = 1,LENGTH
          AU(I) = U(I)
        ENDDO
        RETURN
      ENDIF

      IF(MOD(LENGTH,2).NE.0) LENGTH = LENGTH-1
      PERIODS = FLOAT(LENGTH)*DELMIN
      IF(IOPT.EQ.1) THEN
        IBEGIN = 1
        IEND = NINT(PERIODS/CPERIOD)
      END IF
      IF(IOPT.EQ.2) THEN
        IBEGIN = NINT(PERIODS/CPERIOD)
        IEND = LENGTH/2+1
      END IF
      CALL TFTFTR(U,AU,LENGTH,IBEGIN,IEND,IOPT)

      NLOST = 0
      NDEC = isort
      NEND = LENGTH-2*NLOST
      LENGTH = LENGTH0
      DO I = 1, 5
        AU(I) = U(I)
      ENDDO
      DO I = LENGTH-5, LENGTH
        AU(I) = U(I)
      ENDDO

      RETURN
      END


c-------------------------------------------------------------------------
      SUBROUTINE TFTFTR(XR,XI,NRLPTS,IBEGIN,IEND,IOPT)
      REAL XR(NRLPTS),XI(NRLPTS),WORK(4*NRLPTS),FACTOR(NRLPTS)
      COMPLEX HCDAT(NRLPTS)
      INTEGER TENPCT

      LWRK=4*NRLPTS
C --- COMPUTE AND REMOVE MEAN FORM THE DATA 
C
      SUM = 0.0
      DO I = 1, NRLPTS
        SUM = SUM + XR(I)
      ENDDO
      AVG = SUM/FLOAT(NRLPTS)
      DO I = 1, NRLPTS
        XR(I) = XR(I)-AVG
      ENDDO

      SIGNEX = 1.0
      NHCPTS = NRLPTS/2 + 1

      CALL FFTRC(XR,NRLPTS,SIGNEX,HCDAT,NHCPTS,WORK,LWRK,IERR)
      DO I = 1, NHCPTS
        HCDAT(I) = HCDAT(I)/FLOAT(NRLPTS)
      END DO
  
      NBAND = IEND - IBEGIN + 1
      CALL TENCOS(FACTOR,NBAND,TENPCT)

      GO TO (10,20,30) IOPT
10    CONTINUE   !low-passed filtering
      II = NBAND - TENPCT
      DO I = IEND-TENPCT+1, IEND
        II = II + 1
        HCDAT(I) = FACTOR(II)*HCDAT(I)
      END DO
      DO I = IEND+1, NHCPTS
        HCDAT(I) = CMPLX(0.0,0.0)
      END DO
      GOTO 40

20    CONTINUE   !high-passed filtering
      DO I = 1, IBEGIN-1
        HCDAT(I) = CMPLX(0.0,0.0)
      END DO
      II = 0
      DO I = IBEGIN, IBEGIN+TENPCT
        II = II + 1
        HCDAT(I) = FACTOR(II)*HCDAT(I)
      END DO
      GOTO 40

30    CONTINUE   !band-passed filtering
      DO I = 1, IBEGIN-1
        HCDAT(I) = CMPLX(0.0,0.0)
      END DO
      II=0
      DO I =I BEGIN, IEND
        II = II + 1
        HCDAT(I) = FACTOR(II)*HCDAT(I)
      END DO
      DO I = IEND+1, NHCPTS
        HCDAT(I) = CMPLX(0.0,0.0)
      END DO
40    CONTINUE

C
C  COMPUTE THE INVERSE TRANSFORM(S) OF THE FOURIER COEFFICIENTS AND
C  OUTPUT THE NEW TIME SERIES AFTER SCALING.
C
      SIGNEX = -1.0
      CALL FFTCR(HCDAT,NHCPTS,SIGNEX,XI,NRLPTS,WORK,LWRK,IERR)
C
C ADD MEAN TO THE FILTERED DATA
C
      IF(IOPT .NE. 1) RETURN
      DO I = 1, NRLPTS
        XI(I) = XI(I) + AVG
        XR(I) = XR(I) + AVG
      END DO

      RETURN
      END


c-------------------------------------------------------------------------
      SUBROUTINE TENCOS(A,LENGTH,TENPCT)
      REAL A(LENGTH),SCALE(LENGTH)
      INTEGER TENPCT
C
C  10% COSINE WINDOW.
C
      PI = ACOS(-1.0)
      DO 100 I = 1, LENGTH
        A(I) = 1.0
100   CONTINUE
101   TENPCT = LENGTH/10 + 1

      FACTOR = 10.0*PI/(FLOAT(LENGTH)-1.0)
      DO 110 I = 1, TENPCT
        SCALE(I) = 0.5*(1.0-COS(FACTOR*FLOAT(I-1)))
110   CONTINUE
      J = LENGTH - TENPCT
      K = TENPCT + 1
      DO 111 I = 1, TENPCT
         A(J+I)=A(J+I)*SCALE(K-I)
111   CONTINUE

      RETURN
      END


C PACKAGE FFT            DESCRIPTION OF INDIVIDUAL USER ENTRIES         
C                        FOLLOWS THE PACKAGE DESCRIPTION BELOW.         
C                                                                       
C LATEST REVISION        JANUARY 1985                                   
C                                                                       
C PURPOSE                FAST FOURIER TRANSFORMS FOR DATA OF ARBITRARY  
C                        LENGTH.  THE FILE FFT CONTAINS THREE ROUTINES  
C                        TO HANDLE VARIOUS FORMS OF INPUT DATA AS       
C                        TABULATED BELOW:                               
C                          ROUTINE NAME  INPUT FORM     OUTPUT FORM     
C                          ------------  ----------     -----------     
C                             FFTRC         REAL        HALF COMPLEX    
C                             FFTCR      HALF COMPLEX      REAL         
C                             FFTCC        COMPLEX        COMPLEX       
C                                                                       
C                        HALF COMPLEX REFERS HERE TO THE FIRST N/2+1    
C                        COMPLEX VALUES OF A CONJUGATE SYMMETRIC ARRAY  
C                        OF LENGTH N.                                   
C                                                                       
C SPECIAL CONDITIONS     THE EFFICIENCY OF THESE ROUTINES IS GREATLY    
C                        AFFECTED BY THE NUMBER OF POINTS TO BE         
C                        TRANSFORMED.  IF N IS THE LENGTH OF THE        
C                        TRANSFORM (NRLPTS IN FFTRC OR FFTCR; NCPTS IN  
C                        FFTCC) THE EXECUTION TIME IS ROUGHLY           
C                        PROPORTIONAL TO N*SUMPF WHERE SUMPF IS THE SUM 
C                        OF THE PRIME FACTORS OF N.  CLEARLY, NUMBERS   
C                        WITH LARGE PRIME FACTORS SHOULD BE AVOIDED.    
C                        IF N IS A POWER OF TWO, THE PACKAGE FFTPOW2    
C                        PROVIDES GREATER EFFICIENCY.                   
C                                                                       
C I/O                    ERROR MESSAGES ARE PRINTED BY ROUTINE ULIBER.  
C                                                                       
C PRECISION              SINGLE                                         
C                                                                       
C REQUIRED LIBRARY       ULIBER AND Q8QST4, WHICH ARE LOADED BY         
C FILES                  DEFAULT ON NCAR'S CRAY MACHINES.               
C                                                                       
C LANGUAGE               FORTRAN                                        
C                                                                       
C HISTORY                DEVELOPED AT NCAR BY DAVE FULKER OF THE        
C                        SCIENTIFIC COMPUTING DIVISION IN THE EARLY     
C                        1970'S.                                        
C                                                                       
C PORTABILITY            FORTRAN 66                                     
C***********************************************************************
C                                                                       
C SUBROUTINE FFTRC(RLDAT,NRLPTS,SIGNEX,HCTRN,NHCPTS,WORK,LWRK,IERR)     
C                                                                       
C DIMENSION OF           REAL RLDAT (NRLPTS)                            
C ARGUMENTS              COMPLEX HCTRN (NHCPTS)                         
C                          WHERE NHCPTS = NRLPTS/2+1                    
C                        REAL WORK (LWRK)                               
C                                                                       
C PURPOSE                FOURIER TRANSFORM (REAL TO HALF COMPLEX).  (IN 
C                        CASE SIGNEX = +1., THE COMPUTATION PERFORMED BY
C                        THIS ROUTINE IS SOMETIMES CALLED A FORWARD     
C                        TRANSFORM OR FOURIER ANALYSIS.)  THE DISCRETE  
C                        FOURIER TRANSFORM OF AN ARRAY RLDAT, CONTAINING
C                        NRLPTS REAL VALUES, IS A SET CF OF NRLPTS      
C                        COMPLEX VALUES SATISFYING THE CONJUGATE        
C                        SYMMETRY RELATION CF(NRLPTS+2-K) = CONJG(CF(K))
C                        (K = 2,NRLPTS).  DUE TO THIS SYMMETRY RELATION,
C                        IT IS ONLY NECESSARY TO COMPUTE THE FIRST      
C                        NHCPTS = NRLPTS/2+1 COMPLEX VALUES, AND THESE  
C                        ARE RETURNED IN THE COMPLEX ARRAY HCTRN.       
C                                                                       
C USAGE                  CALL FFTRC (RLDAT,NRLPTS,SIGNEX,HCTRN,NHCPTS,  
C                                    WORK,LWRK,IERR)                    
C                          THE ORIGINAL VALUES OF RLDAT MAY BE          
C                          REGENERATED FROM HCTRN BY FIRST DIVIDING ALL 
C                          VALUES OF HCTRN BY NRLPTS AND THEN CALLING   
C                          FFTCR (HCTRN,NHCPTS,-SIGNEX,RLDAT,NRLPTS,    
C                                 WORK,LWRK,IERR).                      
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON INPUT               RLDAT                                          
C                          A REAL ARRAY CONTAINING THE NRLPTS DATA      
C                          VALUES TO BE TRANSFORMED.                    
C                                                                       
C                        NRLPTS                                         
C                          THE NUMBER OF REAL DATA VALUES TO BE         
C                          TRANSFORMED.  FOR THE GREATEST EFFICIENCY, IT
C                          SHOULD BE A PRODUCT OF SMALL PRIMES.         
C                                                                       
C                        SIGNEX                                         
C                          A VARIABLE WHOSE SIGN DETERMINES THE SIGN OF 
C                          THE ARGUMENT OF THE COMPLEX EXPONENTIAL USED 
C                          IN THE TRANSFORM COMPUTATIONS.  FOR          
C                          CONVENIENCE, WE ASSUME IN THESE COMMENTS THAT
C                          SIGNEX IS +1. OR -1., BUT THE ROUTINE IN FACT
C                          ONLY USES THE SIGN OF ITS VALUE.             
C                                                                       
C                        NHCPTS                                         
C                          THE NUMBER OF COMPLEX VALUES TO BE RETURNED  
C                          AS THE HALF COMPLEX TRANSFORM RESULT IN ARRAY
C                          HCTRN.  IT MUST BE = NRLPTS/2+1 OR A FATAL   
C                          ERROR IS FLAGGED.                            
C                          NOTE:  NHCPTS IS NOT AN OUTPUT PARAMETER.  IT
C                                 MUST BE SET TO NRLPTS/2+1 BY THE USER,
C                                 AND NHCPTS COMPLEX LOCATIONS MUST BE  
C                                 PROVIDED FOR THE OUTPUT ARRAY HCTRN.  
C                                                                       
C                        WORK                                           
C                          A WORKSPACE OF LENGTH LWRK (.GE. 4*NRLPTS)   
C                          FOR USE BY THE ROUTINE.                      
C                                                                       
C                        LWRK                                           
C                          THE LENGTH OF ARRAY WORK.  IT MUST BE        
C                          .GE. 4*NRLPTS OR A FATAL ERROR IS FLAGGED.   
C                                                                       
C ON OUTPUT              HCTRN                                          
C                          THE COMPLEX ARRAY CONTAINING THE             
C                          FIRST HALF OF THE CONJUGATE SYMMETRIC        
C                          TRANSFORM RESULT.  THESE VALUES ARE          
C                          REFERRED TO AS FOURIER COEFFICIENTS.         
C                                                                       
C                          ASSUME SIGNEX = +1. OR -1. AND CEX(X)        
C                          (FOR REAL X) IS THE COMPLEX EXPONENTIAL      
C                          OF SIGNEX*2*PI*I*X/NRLPTS  WHERE PI = 3.14...
C                          AND I = SQRT(-1.).                           
C                                                                       
C                          THEN (FOR K = 1,NHCPTS)                      
C                            HCTRN(K) =  SUM   RLDAT(J)*CEX((J-1)*(K-1) 
C                          WHERE THE SUM INDEX IS J = 1,....NRLPTS.     
C                                                                       
C                        WORK                                           
C                          WORKSPACE CONTAINING INTERMEDIATE RESULTS.   
C                                                                       
C                        IERR                                           
C                          AN ERROR FLAG WITH THE FOLLOWING MEANINGS:   
C                          =   0  NO ERROR.                             
C                          = 101  NRLPTS IS LESS THAN 1.                
C                          = 102  NHCPTS IS NOT NRLPTS/2+1.             
C                          = 103  INSUFFICIENT WORKSPACE HAS BEEN       
C                                 PROVIDED; LWRK IS LESS THAN 4*NRLPTS. 
C***********************************************************************
C                                                                       
C SUBROUTINE FFTCR (HCDAT,NHCPTS,SIGNEX,RLTRN,NRLPTS,WORK,LWRK,IERR)    
C                                                                       
C DIMENSION OF           COMPLEX HCDAT (NHCPTS)                         
C ARGUMENTS              REAL RLTRN (NRLPTS)   WHERE NRLPTS/2+1 = NHCPTS
C                        REAL WORK (LWRK)                               
C                                                                       
C PURPOSE                FOURIER TRANSFORM (HALF COMPLEX TO REAL).  (IN 
C                        CASE SIGNEX = -1., THE COMPUTATION PERFORMED BY
C                        THIS ROUTINE IS SOMETIMES CALLED A BACKWARD    
C                        TRANSFORM OR FOURIER SYNTHESIS.) HCDAT IS      
C                        ASSUMED TO BE THE FIRST NHCPTS = NRLPTS/2+1    
C                        COMPLEX VALUES OF A SET CF CONTAINING NRLPTS   
C                        COMPLEX VALUES WHICH SATISFY THE CONJUGATE     
C                        SYMMETRY RELATION CF(NRLPTS+2-K) = CONJG(CF(K))
C                        (K = 2,NRLPTS).  IN ADDITION, HCTRN(1) IS      
C                        ASSUMED TO BE REAL.  THE DISCRETE FOURIER      
C                        TRANSFORM OF SUCH A CONJUGATE SYMMETRIC ARRAY  
C                        IS A SET OF NRLPTS REAL VALUES WHICH ARE       
C                        RETURNED IN THE REAL ARRAY RLTRN.              
C                                                                       
C USAGE                  CALL FFTCR (HCDAT,NHCPTS,SIGNEX,RLTRN,NRLPTS,  
C                                    WORK,LWRK,IERR)                    
C                          THE ORIGINAL VALUES OF HCDAT MAY BE          
C                          REGENERATED FROM RLTRN BY FIRST DIVIDING ALL 
C                          VALUES OF RLTRN BY NRLPTS AND THEN CALLING   
C                          FFTRC (RLTRN,NRLPTS,-SIGNEX,HCDAT,NHCPTS,    
C                                 WORK,LWRK,IERR).                      
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON INPUT               HCDAT                                          
C                          A COMPLEX ARRAY CONTAINING NHCPTS COMPLEX    
C                          VALUES WHICH COMPRISE ESSENTIALLY THE FIRST  
C                          HALF OF THE CONJUGATE SYMMETRIC DATA TO BE   
C                          TRANSFORMED.                                 
C                                                                       
C                        NHCPTS                                         
C                          THE NUMBER OF COMPLEX VALUES ENTERED IN      
C                          HCDAT.                                       
C                                                                       
C                        SIGNEX                                         
C                          A VARIABLE WHOSE SIGN DETERMINES THE SIGN OF 
C                          THE ARGUMENT OF THE COMPLEX EXPONENTIAL USED 
C                          IN THE TRANSFORM COMPUTATIONS.  FOR          
C                          CONVENIENCE WE ASSUME IN THESE COMMENTS THAT 
C                          SIGNEX IS +1. OR -1., BUT THE ROUTINE IN FACT
C                          ONLY USES THE SIGN OF ITS VALUE.             
C                                                                       
C                        NRLPTS                                         
C                          THE NUMBER OF REAL TRANSFORM VALUES TO BE    
C                          RETURNED.  NHCPTS MUST BE NRLPTS/2+1 OR A    
C                          FATAL ERROR IS FLAGGED.                      
C                          NOTE:  NRLPTS IS NOT AN OUTPUT PARAMETER.  IT
C                                 MUST SATISFY NRLPTS/2+1 = NHCPTS, AND 
C                                 NRLPTS REAL LOCATIONS MUST BE PROVIDED
C                                 FOR THE OUTPUT ARRAY RLTRN.           
C                                                                       
C                        WORK                                           
C                          A WORKSPACE OF LENGTH LWRK (.GE. 4*NRLPTS)   
C                          FOR USE BY THE ROUTINE.                      
C                                                                       
C                        LWRK                                           
C                          THE LENGTH OF ARRAY WORK.  IT MUST BE        
C                          .GE. 4*NRLPTS OR A FATAL ERROR IS FLAGGED.   
C                                                                       
C ON OUTPUT              RLTRN                                          
C                          REAL ARRAY CONTAINING TRANSFORM RESULT.      
C                                                                       
C                          ASSUME SIGNEX = +1. OR -1. AND CEX(X)        
C                          (FOR REAL X) IS THE COMPLEX EXPONENTIAL      
C                          OF SIGNEX*2*PI*I*X/NRLPTS  WHERE PI = 3.14...
C                          AND I = SQRT(-1.).                           
C                                                                       
C                          THEN FOR J = 1,...,NRLPTS                    
C                            RLTRN(J) =  SUM   CF(K)*CEX(J-1)*(K-1))    
C                          WHERE THE SUM INDEX IS K = 1,...,NRLPTS      
C                          AND WHERE                                    
C                            CF(K) = HCDAT(K)/NRLPTS                    
C                                    FOR 1 .LE. K .LE. NHCPTS           
C                                  = CONJG(CF(NRLPTS+2-K))              
C                                    OTHERWISE                          
C                                                                       
C                        WORK                                           
C                          WORKSPACE CONTAINING INTERMEDIATE RESULTS.   
C                                                                       
C                        IERR                                           
C                          AN ERROR FLAG WITH THE FOLLOWING MEANINGS:   
C                            0  NO ERROR.                               
C                          101  NRLPTS IS LESS THAN 1.                  
C                          102  NHCPTS IS NOT NRLPTS/2+1.               
C                          103  INSUFFICIENT WORKSPACE HAS BEEN         
C                               PROVIDED:  LWRK IS LESS THAN 4*NRLPTS.  
C***********************************************************************
C                                                                       
C SUBROUTINE FFTCC (CDATA,NCPTS,SIGNEX,CTRAN,WORK,IERR)                 
C                                                                       
C DIMENSION OF           COMPLEX CDATA (NCPTS),CTRAN (NCPTS)            
C ARGUMENTS              REAL WORK(2*NCPTS)                             
C                                                                       
C PURPOSE                FOURIER TRANSFORM (COMPLEX TO COMPLEX).  (IN   
C                        CASE SIGNEX = +1., THE COMPUTATION PERFORMED BY
C                        THIS ROUTINE IS SOMETIMES CALLED A FORWARD     
C                        TRANSFORM OR FOURIER ANALYSIS.  FOR            
C                        SIGNEX = -1., IT IS CALLED A BACKWARD TRANSFORM
C                        OR FOURIER SYNTHESIS.)  THE DISCRETE FOURIER   
C                        TRANSFORM OF AN ARRAY CDATA, CONTAINING NCPTS  
C                        COMPLEX VALUES, IS A SET OF NCPTS COMPLEX      
C                        VALUES WHICH ARE RETURNED IN THE COMPLEX ARRAY 
C                        CTRAN.                                         
C                                                                       
C USAGE                  CALL FFTCC (CDATA,NCPTS,SIGNEX,CTRAN,WORK,IERR)
C                          THE ORIGINAL VALUES OF CDATA MAY BE          
C                          REGENERATED FROM CTRAN BY FIRST DIVIDING ALL 
C                          VALUES OF CTRAN BY NCPTS AND THEN CALLING    
C                          FFTCC (CTRAN,NCPTS,-SIGNEX,CDATA,WORK,IERR). 
C                                                                       
C ARGUMENTS                                                             
C                                                                       
C ON INPUT               CDATA                                          
C                          A COMPLEX ARRAY CONTAINING THE NCPTS COMPLEX 
C                          DATA VALUES TO BE TRANSFORMED.               
C                                                                       
C                        NCPTS                                          
C                          THE NUMBER OF COMPLEX DATA VALUES TO BE      
C                          TRANSFORMED.  IT MUST BE A POSITIVE POWER OF 
C                          2 OR A FATAL ERROR IS FLAGGED.               
C                                                                       
C                        SIGNEX                                         
C                          A VARIABLE WHOSE SIGN DETERMINES THE SIGN OF 
C                          THE ARGUMENT OF THE COMPLEX EXPONENTIAL USED 
C                          IN THE TRANSFORM COMPUTATIONS.  FOR          
C                          CONVENIENCE, WE ASSUME IN THESE COMMENTS THAT
C                          SIGNEX IS +1. OR -1., BUT THE ROUTINE IN FACT
C                          ONLY USES THE SIGN OF ITS VALUE.             
C                                                                       
C                        WORK                                           
C                          A WORKSPACE OF LENGTH 2*NCPTS FOR USE BY THE 
C                          ROUTINE.                                     
C                                                                       
C ON OUTPUT              CTRAN                                          
C                          THE COMPLEX ARRAY OF LENGTH NCPTS            
C                          CONTAINING THE TRANSFORMED RESULTS.          
C                          THESE VALUES ARE ALSO REFERRED TO AS THE     
C                          FOURIER COEFFICIENTS.                        
C                                                                       
C                          ASSUME SIGNEX = +1. OR -1. AND CEX(X)        
C                          (FOR REAL X) IS THE COMPLEX EXPONENTIAL      
C                          OF SIGNEX*2*PI*I*X/NRLPTS  WHERE PI = 3.14...
C                          AND I = SQRT(-1.).                           
C                                                                       
C                          FOR K = 1,NCPTS                              
C                            CTRAN(K) =  SUM  CDATA(J)*CEX(J-1)*(K-1))  
C                          WHERE THE SUM INDEX IS J = 1,...,NCPTS.      
C                                                                       
C                        WORK                                           
C                          WORKSPACE CONTAINING INTERMEDIATE RESULTS.   
C                                                                       
C                        IERR                                           
C                          AN ERROR FLAG WITH THE FOLLOWING MEANINGS:   
C                            0  NO ERROR.                               
C                          101  NCPTS IS LESS THAN 1.                   
C***********************************************************************
      SUBROUTINE FFTRC(RLDAT,NRLPTS,SIGNEX,HCTRN,NHCPTS,WORK,LWRK,IERR)
      REAL      RLDAT(NRLPTS), WORK(LWRK)
      COMPLEX   HCTRN(NHCPTS)
      DIMENSION NSCRT(1)

C THE FOLLOWING CALL IS FOR GATHERING STATISTICS ON LIBRARY USE AT NCAR 
c     LOGICAL Q8Q4
c     SAVE Q8Q4
c     DATA Q8Q4 /.TRUE./
c     IF (Q8Q4) THEN
c         CALL Q8QST4('LOCLIB','FFT','FFTRC','VERSION 08')
c         Q8Q4 = .FALSE.
c     ENDIF

      IERR = 0
      IF(NRLPTS .LT. 2) GO TO 103
      IF(NHCPTS .NE. NRLPTS/2+1) GO TO 105
      IF(LWRK .LT. 4*NRLPTS) GO TO 106
      DO 101 J = 1, NRLPTS
        WORK(2*J-1) = RLDAT(J)
        WORK(2*J) = 0.0
101   CONTINUE
      ISIGN = SIGNEX
      NSCRT(1) = NRLPTS
      CALL FOURT (WORK,NSCRT,1,ISIGN,0,WORK(2*NRLPTS+1))
      DO 102 K = 1, NHCPTS
        HCTRN(K) = CMPLX(WORK(2*K-1),WORK(2*K))
102   CONTINUE
      RETURN

103   IF(NRLPTS .LT. 1) GO TO 104
      IF(NHCPTS .NE. 1) GO TO 105
      HCTRN(1) = CMPLX(RLDAT(1),0.0)
      RETURN                                                            

104   IERR = 101
      CALL ULIBER(IERR,' FFTRC   NRLPTS IS .LT. 1',25)
      RETURN

105   IERR = 102
      CALL ULIBER(IERR,' FFTRC   NHCPTS IS NOT NRLPTS/2+1',33)
      RETURN

106   IERR = 103
      CALL ULIBER(IERR,  ' FFTRC   INSUFFICIENT WORKSPACE - LWRK IS .LT
     1. 4*NRLPTS   ',55)

      RETURN
      END


c-------------------------------------------------------------------------
      SUBROUTINE FFTCR(HCDAT,NHCPTS,SIGNEX,RLTRN,NRLPTS,WORK,LWRK,IERR)
      COMPLEX   HCDAT(NHCPTS)
      REAL      RLTRN(NRLPTS), WORK(LWRK)
      DIMENSION NSCRT(1)

C THE FOLLOWING CALL IS FOR GATHERING STATISTICS ON LIBRARY USE AT NCAR
C      LOGICAL Q8Q4
C      SAVE Q8Q4
C      DATA Q8Q4 /.TRUE./
C      IF (Q8Q4) THEN
C          CALL Q8QST4('LOCLIB','FFT','FFTCR','VERSION 08')
C          Q8Q4 = .FALSE.
C      ENDIF

      IERR = 0
      IF(NRLPTS .LT. 2) GO TO 103
      IF(NHCPTS .NE. NRLPTS/2+1) GO TO 105
      IF(LWRK .LT. 4*NRLPTS) GO TO 106
      NC = 2*NRLPTS + 4
      WORK(1) = REAL(HCDAT(1))
      WORK(2) = 0.0
      DO 101 K = 2, NHCPTS
        WORK(2*K-1) = REAL(HCDAT(K))
      NC2 = NC - 2*K - 1
      WORK(NC2) = REAL(HCDAT(K))
      WORK(2*K) = AIMAG(HCDAT(K))
      NC2 = NC - 2*K
      WORK(NC2) = -AIMAG(HCDAT(K))
101   CONTINUE
      ISIGN = SIGNEX
      NSCRT(1) = NRLPTS
      CALL FOURT(WORK,NSCRT,1,ISIGN,1,WORK(2*NRLPTS+1))
      DO 102 J = 1, NRLPTS
        RLTRN(J) = WORK(2*J-1)
102   CONTINUE
      RETURN

103   IF(NRLPTS .LT. 1) GO TO 104
      IF(NHCPTS .NE. 1) GO TO 105
      RLTRN(1) = REAL(HCDAT(1))
      RETURN

104   IERR = 101
      CALL ULIBER (IERR,25H FFTCR   NRLPTS IS .LT. 1,25)
      RETURN

105   IERR = 102
      CALL ULIBER (IERR,33H FFTCR   NHCPTS IS NOT NRLPTS/2+1,33)
      RETURN

106   IERR = 103
      CALL ULIBER (IERR,55H FFTCR   INSUFFICIENT WORKSPACE - LWRK IS .LT
     1. 4*NRLPTS   ,55)

      RETURN
      END


c-------------------------------------------------------------------------
      SUBROUTINE FFTCC(CDATA,NCPTS,SIGNEX,CTRAN,WORK,IERR)
      COMPLEX  CDATA(NCPTS), CTRAN(NCPTS), WORK(NCPTS)         
      DIMENSION  NSCRT(1)

C THE FOLLOWING CALL IS FOR GATHERING STATISTICS ON LIBRARY USE AT NCAR
C      LOGICAL Q8Q4
C      SAVE Q8Q4
C      DATA Q8Q4 /.TRUE./
C      IF (Q8Q4) THEN
C          CALL Q8QST4('LOCLIB','FFT','FFTCC','VERSION 08')
C          Q8Q4 = .FALSE.
C      ENDIF

      IERR = 0
      IF(NCPTS .LT. 1) GO TO 102
      DO 101 J = 1, NCPTS
        CTRAN(J) = CDATA(J)
101   CONTINUE

      ISIGN = SIGNEX
      NSCRT(1) = NCPTS
      CALL FOURT(CTRAN,NSCRT,1,ISIGN,1,WORK)
      RETURN

102   IERR = 101 
      CALL ULIBER (IERR,25H FFTCC    NCPTS IS .LT. 1,25)

      RETURN
      END


C     THE COOLEY-TUKEY FAST FOURIER TRANSFORM IN USASI BASIC FORTRAN    
C                                                                       
C     TRANSFORM(J1,J2,,,,) = SUM(DATA(I1,I2,,,,)*W1**((I2-1)*(J2-1))    
C                                 *W2**((I2-1)*(J2-1))*,,,),            
C     WHERE I1 AND J1 RUN FROM 1 TO NN(1) AND W1=EXP(ISIGN*2*PI=        
C     SQRT(-1)/NN(1)), ETC.  THERE IS NO LIMIT ON THE DIMENSIONALITY    
C     (NUMBER OF SUBSCRIPTS) OF THE DATA ARRAY.  IF AN INVERSE          
C     TRANSFORM (ISIGN=+1) IS PERFORMED UPON AN ARRAY OF TRANSFORMED    
C     (ISIGN=-1) DATA, THE ORIGINAL DATA WILL REAPPEAR.                 
C     MULTIPLIED BY NN(1)*NN(2)*,,,  THE ARRAY OF INPUT DATA MUST BE    
C     IN COMPLEX FORMAT.  HOWEVER, IF ALL IMAGINARY PARTS ARE ZERO (I.E.
C     THE DATA ARE DISGUISED REAL) RUNNING TIME IS CUT UP TO FORTY PER- 
C     CENT.  (FOR FASTEST TRANSFORM OF REAL DATA, NN(1) SHOULD BE EVEN.)
C     THE TRANSFORM VALUES ARE ALWAYS COMPLEX AND ARE RETURNED IN THE   
C     ORIGINAL ARRAY OF DATA, REPLACING THE INPUT DATA.  THE LENGTH     
C     OF EACH DIMENSION OF THE DATA ARRAY MAY BE ANY INTEGER.  THE      
C     PROGRAM RUNS FASTER ON COMPOSITE INTEGERS THAN ON PRIMES, AND IS  
C     PARTICULARLY FAST ON NUMBERS RICH IN FACTORS OF TWO.              
C                                                                       
C     TIMING IS IN FACT GIVEN BY THE FOLLOWING FORMULA.  LET NTOT BE THE
C     TOTAL NUMBER OF POINTS (REAL OR COMPLEX) IN THE DATA ARRAY, THAT  
C     IS, NTOT=NN(1)*NN(2)*...  DECOMPOSE NTOT INTO ITS PRIME FACTORS,  
C     SUCH AS 2**K2 * 3**K3 * 5**K5 * ...  LET SUM2 BE THE SUM OF ALL   
C     THE FACTORS OF TWO IN NTOT, THAT IS, SUM2 = 2*K2.  LET SUMF BE    
C     THE SUM OF ALL OTHER FACTORS OF NTOT, THAT IS, SUMF = 3*K3*5*K5*..
C     THE TIME TAKEN BY A MULTIDIMENSIONAL TRANSFORM ON THESE NTOT DATA 
C     IS T = T0 + NTOT*(T1+T2*SUM2+T3*SUMF).  ON THE CDC 3300 (FLOATING 
C     POINT ADD TIME = SIX MICROSECONDS), T = 3000 + NTOT*(600+40*SUM2+ 
C     175*SUMF) MICROSECONDS ON COMPLEX DATA.                           
C                                                                       
C     IMPLEMENTATION OF THE DEFINITION BY SUMMATION WILL RUN IN A TIME  
C     PROPORTIONAL TO NTOT*(NN(1)+NN(2)+...).  FOR HIGHLY COMPOSITE NTOT
C     THE SAVINGS OFFERED BY THIS PROGRAM CAN BE DRAMATIC.  A ONE-DIMEN-
C     SIONAL ARRAY 4000 IN LENGTH WILL BE TRANSFORMED IN 4000*(600+     
C     40*(2+2+2+2+2)+175*(5+5+5)) = 14.5 SECONDS VERSUS ABOUT 4000*     
C     4000*175 = 2800 SECONDS FOR THE STRAIGHTFORWARD TECHNIQUE.        
C                                                                       
C     THE FAST FOURIER TRANSFORM PLACES THREE RESTRICTIONS UPON THE     
C     DATA.                                                             
C     1.  THE NUMBER OF INPUT DATA AND THE NUMBER OF TRANSFORM VALUES   
C     MUST BE THE SAME.                                                 
C     2.  BOTH THE INPUT DATA AND THE TRANSFORM VALUES MUST REPRESENT   
C     EQUISPACED POINTS IN THEIR RESPECTIVE DOMAINS OF TIME AND         
C     FREQUENCY.  CALLING THESE SPACINGS DELTAT AND DELTAF, IT MUST BE  
C     TRUE THAT DELTAF=2*PI/(NN(I)*DELTAT).  OF COURSE, DELTAT NEED NOT 
C     BE THE SAME FOR EVERY DIMENSION.                                  
C     3.  CONCEPTUALLY AT LEAST, THE INPUT DATA AND THE TRANSFORM OUTPUT
C     REPRESENT SINGLE CYCLES OF PERIODIC FUNCTIONS.                    
C                                                                       
C     THE CALLING SEQUENCE IS--                                         
C     CALL FOURT(DATA,NN,NDIM,ISIGN,IFORM,WORK)                         
C                                                                       
C     DATA IS THE ARRAY USED TO HOLD THE REAL AND IMAGINARY PARTS       
C     OF THE DATA ON INPUT AND THE TRANSFORM VALUES ON OUTPUT.  IT      
C     IS A MULTIDIMENSIONAL FLOATING POINT ARRAY, WITH THE REAL AND     
C     IMAGINARY PARTS OF A DATUM STORED IMMEDIATELY ADJACENT IN STORAGE 
C     (SUCH AS FORTRAN IV PLACES THEM).  NORMAL FORTRAN ORDERING IS     
C     EXPECTED, THE FIRST SUBSCRIPT CHANGING FASTEST.  THE DIMENSIONS   
C     ARE GIVEN IN THE INTEGER ARRAY NN, OF LENGTH NDIM.  ISIGN IS -1   
C     TO INDICATE A FORWARD TRANSFORM (EXPONENTIAL SIGN IS -) AND +1    
C     FOR AN INVERSE TRANSFORM (SIGN IS +).  IFORM IS +1 IF THE DATA ARE
C     COMPLEX, 0 IF THE DATA ARE REAL.  IF IT IS 0, THE IMAGINARY       
C     PARTS OF THE DATA MUST BE SET TO ZERO.  AS EXPLAINED ABOVE, THE   
C     TRANSFORM VALUES ARE ALWAYS COMPLEX AND ARE STORED IN ARRAY DATA. 
C     WORK IS AN ARRAY USED FOR WORKING STORAGE.  IT IS FLOATING POINT  
C     REAL, ONE DIMENSIONAL OF LENGTH EQUAL TO TWICE THE LARGEST ARRAY  
C     DIMENSION NN(I) THAT IS NOT A POWER OF TWO.  IF ALL NN(I) ARE     
C     POWERS OF TWO, IT IS NOT NEEDED AND MAY BE REPLACED BY ZERO IN THE
C     CALLING SEQUENCE.  THUS, FOR A ONE-DIMENSIONAL ARRAY, NN(1) ODD,  
C     WORK OCCUPIES AS MANY STORAGE LOCATIONS AS DATA.  IF SUPPLIED,    
C     WORK MUST NOT BE THE SAME ARRAY AS DATA.  ALL SUBSCRIPTS OF ALL   
C     ARRAYS BEGIN AT ONE.                                              
C                                                                       
C     EXAMPLE 1.  THREE-DIMENSIONAL FORWARD FOURIER TRANSFORM OF A      
C     COMPLEX ARRAY DIMENSIONED 32 BY 25 BY 13 IN FORTRAN IV.           
C     DIMENSION DATA(32,25,13),WORK(50),NN(3)                           
C     COMPLEX DATA                                                      
C     DATA NN/32,25,13/                                                 
C     DO 1 I=1,32                                                       
C     DO 1 J=1,25                                                       
C     DO 1 K=1,13                                                       
C  1  DATA(I,J,K)=COMPLEX VALUE                                         
C     CALL FOURT(DATA,NN,3,-1,1,WORK)                                   
C                                                                       
C     EXAMPLE 2.  ONE-DIMENSIONAL FORWARD TRANSFORM OF A REAL ARRAY OF  
C     LENGTH 64 IN FORTRAN II,                                          
C     DIMENSION DATA(2,64)                                              
C     DO 2 I=1,64                                                       
C     DATA(1,I)=REAL PART                                               
C  2  DATA(2,I)=0.                                                      
C     CALL FOURT(DATA,64,1,-1,0,0)                                      
C                                                                       
C     THERE ARE NO ERROR MESSAGES OR ERROR HALTS IN THIS PROGRAM.  THE  
C     PROGRAM RETURNS IMMEDIATELY IF NDIM OR ANY NN(I) IS LESS THAN ONE.
C                                                                       
C     PROGRAM BY NORMAN BRENNER FROM THE BASIC PROGRAM BY CHARLES       
C     RADER,  JUNE 1967.  THE IDEA FOR THE DIGIT REVERSAL WAS           
C     SUGGESTED BY RALPH ALTER.                                         
C                                                                       
C     THIS IS THE FASTEST AND MOST VERSATILE VERSION OF THE FFT KNOWN   
C     TO THE AUTHOR.  A PROGRAM CALLED FOUR2 IS AVAILABLE THAT ALSO     
C     PERFORMS THE FAST FOURIER TRANSFORM AND IS WRITTEN IN USASI BASIC 
C     FORTRAN.  IT IS ABOUT ONE THIRD AS LONG AND RESTRICTS THE         
C     DIMENSIONS OF THE INPUT ARRAY (WHICH MUST BE COMPLEX) TO BE POWERS
C     OF TWO.  ANOTHER PROGRAM, CALLED FOUR1, IS ONE TENTH AS LONG AND  
C     RUNS TWO THIRDS AS FAST ON A ONE-DIMENSIONAL COMPLEX ARRAY WHOSE  
C     LENGTH IS A POWER OF TWO.                                         
C                                                                       
C     REFERENCE--                                                       
C     IEEE AUDIO TRANSACTIONS (JUNE 1967), SPECIAL ISSUE ON THE FFT.    
                                                                       
      SUBROUTINE FOURT(DATA,NN,NDIM,ISIGN,IFORM,WORK)
      DIMENSION DATA(1), NN(1), IFACT(32), WORK(1)
      DATA NP0/0/,NPREV/0/
      DATA TWOPI/6.2831853071796/, RTHLF/0.70710678118655/

      IF(NDIM-1) 232, 101, 101
101   NTOT = 2
      DO 103 IDIM = 1, NDIM
        IF(NN(IDIM)) 232, 232, 102
102     NTOT = NTOT*NN(IDIM)
103   CONTINUE
C                                                                       
C     MAIN LOOP FOR EACH DIMENSION                                      
C                                                                       
      NP1 = 2
      DO 231 IDIM = 1, NDIM
        N = NN(IDIM)
        NP2 = NP1*N
        IF(N-1) 232, 230, 104
C                                                                       
C     IS N A POWER OF TWO AND IF NOT, WHAT ARE ITS FACTORS              
C                                                                       
104     M = N
        NTWO = NP1
        IF2 = 1
        IDIV = 2

105     IQUOT = M/IDIV
        IREM = M-IDIV*IQUOT
        IF(IQUOT-IDIV) 113, 106, 106
106     IF(IREM) 108, 107, 108
107     NTWO = NTWO+NTWO
        IFACT(IF2) = IDIV
        IF2 = IF2+1
        M = IQUOT
        GO TO 105
108     IDIV = 3
        INON2 = IF2
109     IQUOT = M/IDIV
        IREM = M-IDIV*IQUOT
        IF(IQUOT-IDIV) 115, 110, 110
110     IF(IREM) 112, 111, 112
111     IFACT(IF2) = IDIV
        IF2 = IF2+1   
        M = IQUOT
        GO TO 109
112     IDIV = IDIV+2
        GO TO 109
113     INON2 = IF2
        IF(IREM) 115, 114, 115
114     NTWO = NTWO+NTWO
        GO TO 116
115     IFACT(IF2) = M

C                                                                       
C     SEPARATE FOUR CASES--                                             
C        1. COMPLEX TRANSFORM OR REAL TRANSFORM FOR THE 4TH, 9TH,ETC.   
C           DIMENSIONS.                                                 
C        2. REAL TRANSFORM FOR THE 2ND OR 3RD DIMENSION.  METHOD--      
C           TRANSFORM HALF THE DATA, SUPPLYING THE OTHER HALF BY CON-   
C           JUGATE SYMMETRY.                                            
C        3. REAL TRANSFORM FOR THE 1ST DIMENSION, N ODD.  METHOD--      
C           SET THE IMAGINARY PARTS TO ZERO.                            
C        4. REAL TRANSFORM FOR THE 1ST DIMENSION, N EVEN.  METHOD--     
C           TRANSFORM A COMPLEX ARRAY OF LENGTH N/2 WHOSE REAL PARTS    
C           ARE THE EVEN NUMBERED REAL VALUES AND WHOSE IMAGINARY PARTS 
C           ARE THE ODD NUMBERED REAL VALUES.  SEPARATE AND SUPPLY      
C           THE SECOND HALF BY CONJUGATE SYMMETRY.                      
C                                                                       
116     ICASE = 1
        IFMIN = 1
        I1RNG = NP1
        IF(IDIM-4) 117, 122, 122                                        
117     IF(IFORM) 118, 118, 122                                         
118     ICASE = 2                                                      
        I1RNG = NP0*(1+NPREV/2)                                        
        IF(IDIM-1) 119, 119, 122                                        
119     ICASE = 3                                                      
        I1RNG = NP1                                                    
        IF(NTWO-NP1) 122, 122, 120                                      
120     ICASE = 4                                                      
        IFMIN = 2                                                      
        NTWO = NTWO/2                                                  
        N = N/2                                                        
        NP2 = NP2/2                                                    
        NTOT = NTOT/2                                                  
        I = 1                                                          
        DO 121 J = 1, NTOT                                                
          DATA(J) = DATA(I)                                           
          I = I + 2                                                     
121     CONTINUE                                                       

C                                                                       
C     SHUFFLE DATA BY BIT REVERSAL, SINCE N=2**K.  AS THE SHUFFLING     
C     CAN BE DONE BY SIMPLE INTERCHANGE, NO WORKING ARRAY IS NEEDED     
C                                                                       
122     IF(NTWO-NP2) 132, 123, 123                                      
123     NP2HF = NP2/2                                                  
        J = 1                                                          
        DO 131 I2 = 1, NP2, NP1                                            
          IF(J-I2) 124, 127, 127                                       
124       I1MAX = I2+NP1-2                                            
          DO 126 I1 = I2, I1MAX, 2                                        
            DO 125 I3 = I1, NTOT, NP2                                    
              J3 = J+I3-I2                                          
              TEMPR = DATA(I3)                                      
              TEMPI = DATA(I3+1)                                    
              DATA(I3) = DATA(J3)                                   
              DATA(I3+1) = DATA(J3+1)                               
              DATA(J3) = TEMPR                                      
              DATA(J3+1) = TEMPI                                    
125         CONTINUE                                                 
126       CONTINUE                                                    
127       M = NP2HF                                                   
128       IF(J-M) 130, 130, 129                                        
129       J = J-M                                                     
          M = M/2                                                     
          IF(M-NP1) 130, 128, 128                                      
130       J = J+M                                                     
131     CONTINUE                                                       
        GO TO 142                                                      
C                                                                       
C     SHUFFLE DATA BY DIGIT REVERSAL FOR GENERAL N                      
C                                                                       
132     NWORK = 2*N                                                    
        DO 141 I1 = 1, NP1, 2                                              
          DO 140 I3 = I1, NTOT, NP2                                       
            J = I3                                                   
            DO 138 I = 1, NWORK, 2                                       
              IF(ICASE-3) 133, 134, 133                              
133           WORK(I) = DATA(J)                                     
              WORK(I+1) = DATA(J+1)                                 
              GO TO 135                                             
134           WORK(I) = DATA(J)                                     
              WORK(I+1) = 0.0                                        
135           IFP2 = NP2                                            
              IF2 = IFMIN                                            
136           IFP1 = IFP2/IFACT(IF2)                                 
              J = J+IFP1                                            
              IF(J-I3-IFP2) 138, 137, 137                            
137           J = J-IFP2                                            
              IFP2 = IFP1                                           
              IF2 = IF2+1                                             
              IF(IFP2-NP1) 138, 138, 136                             
138         CONTINUE                                                 
            I2MAX = I3+NP2-NP1                                       
            I = 1                                                    
            DO 139 I2 = I3, I2MAX, NP1                                   
              DATA(I2) = WORK(I)                                    
              DATA(I2+1) = WORK(I+1)                                
              I = I+2                                               
139         CONTINUE                                                 
140       CONTINUE                                                    
141     CONTINUE                                                       
C                                                                       
C     MAIN LOOP FOR FACTORS OF TWO.  PERFORM FOURIER TRANSFORMS OF      
C     LENGTH FOUR, WITH ONE OF LENGTH TWO IF NEEDED.  THE TWIDDLE FACTOR
C     W=EXP(ISIGN*2*PI*SQRT(-1)*M/(4*MMAX)).  CHECK FOR W=ISIGN*SQRT(-1)
C     AND REPEAT FOR W=W*(1+ISIGN*SQRT(-1))/SQRT(2).                    
C                                                                       
142     IF(NTWO-NP1) 174, 174, 143                                      
143     NP1TW = NP1+NP1                                                
        IPAR = NTWO/NP1                                                
144     IF (IPAR-2) 149, 146, 145                                        
145     IPAR = IPAR/4                                                  
        GO TO 144                                                      
146     DO 148 I1 = 1, I1RNG, 2                                            
          DO 147 K1 = I1, NTOT, NP1TW                                     
            K2 = K1+NP1                                              
            TEMPR = DATA(K2)                                         
            TEMPI = DATA(K2+1)                                       
            DATA(K2) = DATA(K1)-TEMPR                                
            DATA(K2+1) = DATA(K1+1)-TEMPI                            
            DATA(K1) = DATA(K1)+TEMPR                                
            DATA(K1+1) = DATA(K1+1)+TEMPI                            
147       CONTINUE                                                    
148     CONTINUE                                                       

149     MMAX = NP1                                                     
150     IF(MMAX-NTWO/2) 151, 174, 174                                   
151     LMAX = MAX0(NP1TW,MMAX/2)                                      
        DO 173 L = NP1, LMAX, NP1TW                                        
          M = L                                                       
          IF(MMAX-NP1) 156, 156, 152                                   
152       THETA = -TWOPI*FLOAT(L)/FLOAT(4*MMAX)                       
          IF(ISIGN) 154, 153, 153                                      
153       THETA = -THETA                                              
154       WR = COS(THETA)                                             
          WI = SIN(THETA)                                             
155       W2R = WR*WR-WI*WI                                           
          W2I = 2.0*WR*WI                                              
          W3R = W2R*WR-W2I*WI                                         
          W3I = W2R*WI+W2I*WR                                         
156       DO 169 I1 = 1, I1RNG, 2                                         
            KMIN = I1+IPAR*M                                         
            IF(MMAX-NP1) 157, 157, 158                                
157         KMIN = I1                                                
158         KDIF = IPAR*MMAX                                         
159         KSTEP = 4*KDIF                                           
            IF(KSTEP-NTWO) 160, 160, 169                              
160         DO 168 K1 = KMIN, NTOT, KSTEP                                
              K2 = K1+KDIF                                          
              K3 = K2+KDIF                                          
              K4 = K3+KDIF                                          
              IF(MMAX-NP1) 161, 161, 164                             
161           U1R = DATA(K1)+DATA(K2)                               
              U1I = DATA(K1+1)+DATA(K2+1)                           
              U2R = DATA(K3)+DATA(K4)                               
              U2I = DATA(K3+1)+DATA(K4+1)                           
              U3R = DATA(K1)-DATA(K2)                               
              U3I = DATA(K1+1)-DATA(K2+1)                           
              IF(ISIGN) 162, 163, 163                                
162           U4R = DATA(K3+1)-DATA(K4+1)                           
              U4I = DATA(K4)-DATA(K3)                               
              GO TO 167                                             
163           U4R = DATA(K4+1)-DATA(K3+1)                           
              U4I = DATA(K3)-DATA(K4)                               
              GO TO 167                                             
164           T2R = W2R*DATA(K2)-W2I*DATA(K2+1)                     
              T2I = W2R*DATA(K2+1)+W2I*DATA(K2)                     
              T3R = WR*DATA(K3)-WI*DATA(K3+1)                       
              T3I = WR*DATA(K3+1)+WI*DATA(K3)                       
              T4R = W3R*DATA(K4)-W3I*DATA(K4+1)                     
              T4I = W3R*DATA(K4+1)+W3I*DATA(K4)                     
              U1R = DATA(K1)+T2R                                    
              U1I = DATA(K1+1)+T2I                                  
              U2R = T3R+T4R                                         
              U2I = T3I+T4I                                         
              U3R = DATA(K1)-T2R                                    
              U3I = DATA(K1+1)-T2I                                  
              IF(ISIGN) 165, 166, 166                                
165           U4R = T3I-T4I                                         
              U4I = T4R-T3R                                         
              GO TO 167                                             
166           U4R = T4I-T3I                                         
              U4I = T3R-T4R                                         
167           DATA(K1) = U1R+U2R                                    
              DATA(K1+1) = U1I+U2I                                  
              DATA(K2) = U3R+U4R                                    
              DATA(K2+1) = U3I+U4I                                  
              DATA(K3) = U1R-U2R                                    
              DATA(K3+1) = U1I-U2I                                  
              DATA(K4) = U3R-U4R                                    
              DATA(K4+1) = U3I-U4I                                  
168         CONTINUE                                                 
            KDIF = KSTEP                                             
            KMIN = 4*(KMIN-I1)+I1                                    
            GO TO 159                                                
169       CONTINUE                                                    
          M = M+LMAX                                                  
          IF(M-MMAX) 170, 170, 173                                     
170       IF(ISIGN) 171, 172, 172                                      
171       TEMPR = WR                                                  
          WR = (WR+WI)*RTHLF                                          
          WI = (WI-TEMPR)*RTHLF                                       
          GO TO 155                                                   
172       TEMPR = WR                                                  
          WR = (WR-WI)*RTHLF                                          
          WI = (TEMPR+WI)*RTHLF                                       
          GO TO 155                                                   
173     CONTINUE                                                       
        IPAR = 3-IPAR                                                  
        MMAX = MMAX+MMAX                                               
        GO TO 150                                                      
C                                                                       
C     MAIN LOOP FOR FACTORS NOT EQUAL TO TWO.  APPLY THE TWIDDLE FACTOR 
C     W=EXP(ISIGN*2*PI*SQRT(-1)*(J1-1)*(J2-J1)/(IFP1+IFP2)), THEN       
C     PERFORM A FOURIER TRANSFORM OF LENGTH IFACT(IF), MAKING USE OF    
C     CONJUGATE SYMMETRIES.                                             
C                                                                       
174     IF(NTWO-NP2) 175, 201, 201                                      
175     IFP1 = NTWO                                                    
        IF2 = INON2                                                     
        NP1HF = NP1/2                                                  
176     IFP2 = IFACT(IF2)*IFP1                                          
        J1MIN = NP1+1                                                  
        IF(J1MIN-IFP1) 177, 177, 184                                    
177     DO 183 J1 = J1MIN, IFP1, NP1                                       
          THETA = -TWOPI*FLOAT(J1-1)/FLOAT(IFP2)                      
          IF (ISIGN) 179, 178, 178                                      
178       THETA = -THETA                                              
179       WSTPR = COS(THETA)                                          
          WSTPI = SIN(THETA)                                          
          WR = WSTPR                                                  
          WI = WSTPI                                                  
          J2MIN = J1+IFP1                                             
          J2MAX = J1+IFP2-IFP1                                        
          DO 182 J2 = J2MIN, J2MAX, IFP1                                  
            I1MAX = J2+I1RNG-2                                       
            DO 181 I1 = J2, I1MAX, 2                                     
              DO 180 J3 = I1, NTOT, IFP2                                
                TEMPR = DATA(J3)                                   
                DATA(J3) = DATA(J3)*WR-DATA(J3+1)*WI               
                DATA(J3+1) = TEMPR*WI+DATA(J3+1)*WR                
180           CONTINUE                                              
181         CONTINUE                                                 
            TEMPR = WR                                               
            WR = WR*WSTPR-WI*WSTPI                                   
            WI = TEMPR*WSTPI+WI*WSTPR                                
182       CONTINUE                                                    
183     CONTINUE                                                       
184     THETA = -TWOPI/FLOAT(IFACT(IF2))                                
        IF (ISIGN) 186, 185, 185                                         
185     THETA = -THETA                                                 
186     WSTPR = COS(THETA)                                             
        WSTPI = SIN(THETA)                                             
        J2RNG = IFP1*(1+IFACT(IF2)/2)                                   
        DO 200 I1=1,I1RNG,2                                            
          DO 199 I3=I1,NTOT,NP2                                       
            J2MAX = I3+J2RNG-IFP1                                    
            DO 197 J2=I3,J2MAX,IFP1                                  
              J1MAX = J2+IFP1-NP1                                   
              DO 193 J1=J2,J1MAX,NP1                                
                J3MAX = J1+NP2-IFP2                                
                DO 192 J3=J1,J3MAX,IFP2                            
                  JMIN = J3-J2+I3                                 
                  JMAX = JMIN+IFP2-IFP1                           
                  I = 1+(J3-I3)/NP1HF                             
                  IF (J2-I3) 187, 187, 189                          
187               SUMR = 0.0
                  SUMI = 0.0
                  DO 188 J=JMIN,JMAX,IFP1                         
                    SUMR = SUMR+DATA(J)                          
                    SUMI = SUMI+DATA(J+1)                        
188               CONTINUE                                        
                  WORK(I) = SUMR                                  
                  WORK(I+1) = SUMI                                
                  GO TO 192                                       
189               ICONJ = 1+(IFP2-2*J2+I3+J3)/NP1HF               
                  J = JMAX                                        
                  SUMR = DATA(J)                                  
                  SUMI = DATA(J+1)                                
                  OLDSR = 0.0                                     
                  OLDSI = 0.0                                     
                  J = J-IFP1                                      
190               TEMPR = SUMR                                    
                  TEMPI = SUMI                                    
                  SUMR = TWOWR*SUMR-OLDSR+DATA(J)                 
                  SUMI = TWOWR*SUMI-OLDSI+DATA(J+1)               
                  OLDSR = TEMPR                                   
                  OLDSI = TEMPI                                   
                  J = J-IFP1                                      
                  IF (J-JMIN) 191, 191, 190                         
191               TEMPR = WR*SUMR-OLDSR+DATA(J)                   
                  TEMPI = WI*SUMI                                 
                  WORK(I) = TEMPR-TEMPI                           
                  WORK(ICONJ) = TEMPR+TEMPI                       
                  TEMPR = WR*SUMI-OLDSI+DATA(J+1)                 
                  TEMPI = WI*SUMR                                 
                  WORK(I+1) = TEMPR+TEMPI                         
                  WORK(ICONJ+1) = TEMPR-TEMPI                     
192             CONTINUE                                           
193           CONTINUE                                              
              IF (J2-I3) 194, 194, 195                                
194           WR = WSTPR                                            
              WI = WSTPI                                            
              GO TO 196                                             
195           TEMPR = WR                                            
              WR = WR*WSTPR-WI*WSTPI                                
              WI = TEMPR*WSTPI+WI*WSTPR                             
196           TWOWR = WR+WR                                         
197         CONTINUE                                                 
            I = 1                                                    
            I2MAX = I3+NP2-NP1                                       
            DO 198 I2=I3,I2MAX,NP1                                   
              DATA(I2) = WORK(I)                                    
              DATA(I2+1) = WORK(I+1)                                
              I = I+2                                               
198         CONTINUE                                                 
199       CONTINUE                                                    
200     CONTINUE                                                       
        IF2 = IF2+1                                                      
        IFP1 = IFP2                                                    
        IF (IFP1-NP2) 176, 201, 201                                      
C                                                                       
C     COMPLETE A REAL TRANSFORM IN THE 1ST DIMENSION, N EVEN, BY CON-   
C     JUGATE SYMMETRIES.                                                
C                                                                       
201     GO TO (230,220,230,202),ICASE                                  
202     NHALF = N                                                      
        N = N+N                                                        
        THETA = -TWOPI/FLOAT(N)                                        
        IF (ISIGN) 204, 203, 203                                         
203     THETA = -THETA                                                 
204     WSTPR = COS(THETA)                                             
        WSTPI = SIN(THETA)                                             
        WR = WSTPR                                                     
        WI = WSTPI                                                     
        IMIN = 3                                                       
        JMIN = 2*NHALF-1                                               
        GO TO 207                                                      
205     J = JMIN                                                       
        DO 206 I=IMIN,NTOT,NP2                                         
          SUMR = (DATA(I)+DATA(J))/2.0
          SUMI = (DATA(I+1)+DATA(J+1))/2.0
          DIFR = (DATA(I)-DATA(J))/2.                                 
          DIFI = (DATA(I+1)-DATA(J+1))/2.                             
          TEMPR = WR*SUMI+WI*DIFR                                     
          TEMPI = WI*SUMI-WR*DIFR                                     
          DATA(I) = SUMR+TEMPR                                        
          DATA(I+1) = DIFI+TEMPI                                      
          DATA(J) = SUMR-TEMPR                                        
          DATA(J+1) = -DIFI+TEMPI                                     
          J = J+NP2                                                   
206     CONTINUE                                                       
        IMIN = IMIN+2                                                  
        JMIN = JMIN-2                                                  
        TEMPR = WR                                                     
        WR = WR*WSTPR-WI*WSTPI                                         
        WI = TEMPR*WSTPI+WI*WSTPR                                      
207     IF (IMIN-JMIN) 205, 208, 211                                     
208     IF (ISIGN) 209, 211, 211                                         
209     DO 210 I=IMIN,NTOT,NP2                                         
          DATA(I+1) = -DATA(I+1)                                      
210     CONTINUE                                                       
211     NP2 = NP2+NP2                                                  
        NTOT = NTOT+NTOT                                               
        J = NTOT+1                                                     
        IMAX = NTOT/2+1                                                
212     IMIN = IMAX-2*NHALF                                            
        I = IMIN                                                       
        GO TO 214                                                      
213     DATA(J) = DATA(I)                                              
        DATA(J+1) = -DATA(I+1)                                         
214     I = I+2                                                        
        J = J-2                                                        
        IF (I-IMAX) 213, 215, 215                                        
215     DATA(J) = DATA(IMIN)-DATA(IMIN+1)                              
        DATA(J+1) = 0.0                                             
        IF (I-J) 217, 219, 219                                           
216     DATA(J) = DATA(I)                                              
        DATA(J+1) = DATA(I+1)                                          
217     I = I-2                                                        
        J = J-2                                                        
        IF (I-IMIN) 218, 218, 216                                        
218     DATA(J) = DATA(IMIN)+DATA(IMIN+1)                              
        DATA(J+1) = 0.0                                                
        IMAX = IMIN                                                    
        GO TO 212                                                      
219     DATA(1) = DATA(1)+DATA(2)                                      
        DATA(2) = 0.0                                                 
        GO TO 230                                                      
C                                                                       
C     COMPLETE A REAL TRANSFORM FOR THE 2ND OR 3RD DIMENSION BY         
C     CONJUGATE SYMMETRIES.                                             
C                                                                       
220     IF (I1RNG-NP1) 221, 230, 230                                     
221     DO 229 I3 = 1, NTOT, NP2                                           
          I2MAX = I3+NP2-NP1                                          
          DO 228 I2 = I3, I2MAX, NP1                                      
            IMIN = I2+I1RNG                                          
            IMAX = I2+NP1-2                                          
            JMAX = 2*I3+NP1-IMIN                                     
            IF (I2-I3) 223, 223, 222                                   
222         JMAX = JMAX+NP2                                          
223         IF (IDIM-2) 226, 226, 224                                  
224         J = JMAX+NP0                                             
            DO 225 I = IMIN, IMAX, 2                                     
              DATA(I) = DATA(J)                                     
              DATA(I+1) = -DATA(J+1)                                
              J = J-2                                               
225         CONTINUE                                                 
226         J = JMAX                                                 
            DO 227 I = IMIN, IMAX, NP0                                   
              DATA(I) = DATA(J)                                     
              DATA(I+1) = -DATA(J+1)                                
              J = J-NP0                                             
227         CONTINUE                                                 
228       CONTINUE                                                    
229     CONTINUE                                                       
C                                                                       
C     END OF LOOP ON EACH DIMENSION                                     
C                                                                       
230     NP0 = NP1                                                      
        NP1 = NP2                                                      
        NPREV = N                                                      
231   CONTINUE

232   RETURN                                                            
      END                                                               


      SUBROUTINE ULIBER(IERR,MESSAGES,NMESS)
      INTEGER IERR
      CHARACTER MESSAGES(NMESS)
   
      WRITE(6,100) IERR,MESSAGES
100   FORMAT(' IERR:',I3,80A)

      RETURN
      END


C --- THIS ROUTINE COMPUTES AND SUBTRACTS FROM THE SERIES X EITHER
C ---   THE SERIES MEAN OR THE LEAST SQUARES STRAIGHT LINE.
C ---   INPUT PARAMETERS ARE
C --- X(I) CONTAINS THE TIME SERIES AND ON RETURN THE DETRENDED 
C ---   TIME SERIES
C --- N IS THE SERIES LENGTH

      SUBROUTINE DETREND(X,N)
      DIMENSION X(N)
      SUMX = 0.0
      DO 20 I = 1, N
20    SUMX = SUMX + X(I)
      XBAR = SUMX / FLOAT(N)

      NUE = N/3
      SUM1 = 0.0
      DO I = N-NUE, N
        SUM1 = SUM1 + X(I)
      END DO

      SUM2 = 0.0
      DO I = 1, NUE
        SUM2 = SUM2 + X(I)
      END DO

      ALPHA = (SUM1-SUM2)/FLOAT(NUE)/FLOAT(N-NUE)
      TD2 = FLOAT(N)/2.0

      DO I = 1, N
        X(I) = X(I) - XBAR - ALPHA*(FLOAT(I)-TD2)
      END DO

      RETURN
      END


      SUBROUTINE CMSV(A,NOBS,AMEAN,UMEAN,LMEAN,ASD,USD,LSD)
      REAL A(NOBS)
      REAL AMEAN,UMEAN,LMEAN,ASD,USD,LSD
      EXTERNAL TDIS

C  COMPUTE THE MEAN VALUE
      RNOBS=FLOAT(NOBS)
      SUM = 0.0
      DO 40 I = 1, NOBS
        SUM = SUM + A(I)
40    CONTINUE
      AMEAN = SUM/RNOBS

C  DE-MEAN TIME SERIES
      DO 50 I = 1, NOBS
        A(I) = A(I) - AMEAN
50    CONTINUE

C  CALCULATION OF THE STANDARD DEVIATION
      SUM = 0.0
      DO 60 I = 1, NOBS
        SUM = SUM + A(I)*A(I)
60    CONTINUE
      ASD = SQRT(SUM/FLOAT(NOBS - 1))
      AVAR = ASD*ASD
      N = NOBS-1
      CALL CONINT(N,CUPPER,CLOWER)
      USD = SQRT(CUPPER*AVAR)
      LSD = SQRT(CLOWER*AVAR)
      RN  = FLOAT(N)

      UMEAN = AMEAN+ASD*TDIS(RN,1)/SQRT(RNOBS)
      LMEAN = AMEAN-ASD*TDIS(RN,1)/SQRT(RNOBS)

C RESTORE INPUT TIME SERIES
      DO 70 I = 1, NOBS
        A(I) = A(I) + AMEAN
70    CONTINUE

      RETURN
      END


C --- PERCENTILE FOR T-DISTRIBUUTION WITH EF DEGREES OF FREEDOM
C ---    IAL (1,2) FOR 95% OR 80% PROBABILITY LEVEL. FROM
C ---    ABRAMOWITZ P949
      FUNCTION TDIS(EF,IAL)
      X=1.96
      IF(IAL .EQ. 2) X = 1.282
      X3 = X**3
      X5 = X**5
      X7 = X**7
      G1 = (X3+X)/4
      G2 = (5*X5+16*X3+3*X)/96
      G3 = (3*X7+19*X5+17*X3-15*X)/384
      TDIS = X+(G1+(G2+G3/EF)/EF)/EF

      RETURN
      END


C  THIS SUBROUTINE CALCULATES A 95% CONFIDENCE INTERVAL FOR SPECTRAL 
C  ESTIMATES GIVEN THE NUMBER OF DEGREES OF FREEDOM. THE MULTIPLICATIVE 
C  CONFIDENCE INTERVAL FACTORS ARE CALCULATED FROM VALUES OF THE CHI 
C  SQUARE DISTRIBUTION.
C
C  DEGREES OF FREEDOM MUST BE LESS THAN OR EQUAL TO 240.
C
C      VARIABLE LIST
C
C  NAVG          THE NUMBER OF AVERAGING OVER SPECTRAL WINDOW
C  NDOF          THE NUMBER OF DEGREES OF FREEDOM.
C  CUPPER        THE MULTIPLICATIVE FACTOR USED IN DETERMINING THE UPPER
C                VALUE OF THE CONFIDENCE INTERVAL.
C  CLOWER        THE MULTIPLICATIVE FACTOR USED IN DETERMINING THE LOWER
C                VALUE OF THE CONFIDENCE INTERVAL.
C  RNDOF         SAME AS NDOF BUT REAL.
C  CP025         ARRAY CONTAINING VALUES OF THE CHI SQUARE DISTRIBUTION
C                USED IN CALCULATING THE LOWER VALUE OF THE CONFIDENCE
C                INTERVAL.
C  CP975         ARRAY CONTAINING VALUES OF THE CHI SQUARE DISTRIBUTION
C                USED IN CALCULATING THE UPPER VALUE OF THE CONFIDENCE
C                INTERVAL.
C
      SUBROUTINE CONINT(NDOF,CUPPER,CLOWER)
      REAL CP025(33),CP975(33)
      DATA CP975/0.00098,0.0506,0.216,0.484,0.831,1.24, 1.69, 2.18,
     *        2.70, 3.25, 3.82, 4.40, 5.01, 5.63, 6.26, 6.91, 7.56,
     *        8.23, 8.91, 9.59,10.28,10.98,11.69,12.40,13.12,13.84,
     *       14.57,15.31,16.05,16.79,24.43,40.48,91.58/
      DATA CP025/5.02, 7.38, 9.35,11.14,12.83, 14.45,16.01,17.53,
     *    19.02,20.48,21.92,23.34,24.47,26.12, 27.49,28.85,30.19,
     *    31.53,32.85,34.17,35.48,36.78,38.08, 39.36,40.65,41.92,
     *    43.19,44.46,45.72,46.98,59.34,83.30,152.21/ 
C
C --- COMPUTE THE NUMBER OF DENSITY SPECTRA AVERAGED OVER THE SPECTRAL 
C     WIDTH
C
      NAVG = NDOF
C
C  IF LESS THAN 61 DEGREES OF FREEDOM, SET UPPER AND
C  LOWER CONFIDENCE INTERVAL FACTORS USING VALUES OF THE CHI SQUARE DIS-
C  TRIBUTION.
C
      RNAVG = FLOAT(NAVG)
      IF(NAVG .GT. 30) GO TO 100
      CUPPER = RNAVG/CP975(NAVG)
      CLOWER = RNAVG/CP025(NAVG)
      GO TO 300
C
C  IF LESS THAN 80 DEGREES OF FREEDOM(BUT MORE THAN 60),SET UPPER AND
C  LOWER CONFIDENCE INTERVAL FACTORS BY INTERPOLATING VALUES OF THE CHI
C  SQUARE DISTRIBUTION.
C
100   IF(NAVG .GT. 40) GO TO 200
      J = 30
      CUPPER = RNAVG/
     1        (CP975(J)+AMOD(RNAVG,10.0)/10.0*(CP975(J+1)-CP975(J)))
      CLOWER = RNAVG/
     1        (CP025(J)+AMOD(RNAVG,10.0)/10.0*(CP025(J+1)-CP025(J)))
      GO TO 300
C
C  IF LESS THAN 120 DEGREES OF FREEDOM(BUT MORE THAN 80),SET UPPER AND
C  LOWER CONFIDENCE INTERVAL FACTORS BY INTERPOLATING VALUES OF THE CHI
C  SQUARE DISTRIBUTION.
C
200   IF(NAVG .GT. 60) GO TO 210
      J = 31
      CUPPER = RNAVG/
     1        (CP975(J)+AMOD(RNAVG,20.0)/20.0*(CP975(J+1)-CP975(J)))
      CLOWER = RNAVG/
     1        (CP025(J)+AMOD(RNAVG,20.0)/20.0*(CP025(J+1)-CP025(J)))
      GO TO 300
C
C  IF LESS THAN 240 DEGREES OF FREEDOM(BUT MORE THAN 120),SET UPPER AND
C  LOWER CONFIDENCE INTERVAL FACTORS BY INTERPOLATING VALUES OF THE CHI
C  SQUARE DISTRIBUTION.
C
210   IF(NAVG .GT. 120) GO TO 220
      J = 32
      CUPPER = RNAVG/
     1        (CP975(J)+AMOD(RNAVG,60.0)/60.0*(CP975(J+1)-CP975(J)))
      CLOWER = RNAVG/
     1        (CP025(J)+AMOD(RNAVG,60.0)/60.0*(CP025(J+1)-CP025(J)))
      GO TO 300
C
C  IF GREATER THAN 240 DEGREES OF FREEDOM, USE THE FORMULA LISTED ON
C  P.388, BENDAT & PIERSOL, 1971 
C
220   Zalpha = 1.96
      TERM = 2.0/9.0/RNAVG
      TERMSQ = SQRT(TERM)
      CUPPER = RNAVG/(RNAVG*(1.0-TERM-ZALPHA*TERMSQ)**3)
      CLOWER = RNAVG/(RNAVG*(1.0-TERM+ZALPHA*TERMSQ)**3)
300   CONTINUE

      RETURN
      END

