      CHARACTER*10   ALIST
      CHARACTER*80   HEAD(2)
      DIMENSION ISTA(6),NO(6),AMP1(37),EPOC1(37),AMP2(37),EPOC2(37)
      DIMENSION   ALIST(37)
      DATA (ALIST(I),I=1,37) /'M(2)      ','S(2)      ','N(2)      ',
     1                        'K(1)      ','M(4)      ','O(1)      ',
     2                        'M(6)      ','MK(3)     ','S(4)      ',
     3                        'MN(4)     ','NU(2)     ','S(6)      ',
     4                        'MU(2)     ','2N(2)     ','OO(1)     ',
     5                        'LAMBDA(2)  ','S(1)      ','M(1)      ',
     6                        'J(1)      ','MM        ','SSA       ',
     7                        'SA        ','MSF       ','MF        ',
     8                        'RHO(1)    ','Q(1)      ','T(2)      ',
     9                        'R(2)      ','2Q(1)     ','P(1)      ',
     1                        '2SM(2)    ','M(3)      ','L(2)      ',
     2                        '2MK(3)   ','K(2)      ','M(8)      ',
     3                        'MS(4)     '/
      open(30,file='cons.out' )

111   READ (30,550)  HEAD(1),HEAD(2)
      READ (30,532)DATUM,ISTA(1),NO(1),(AMP1(J),EPOC1(J),J=1,7),
     1 ISTA(2),NO(2),
     2 (AMP1(J),EPOC1(J),J=8,14),ISTA(3),NO(3),
     3 (AMP1(J),EPOC1(J),J=15,21),ISTA(4),NO(4),
     4 (AMP1(J),EPOC1(J),J=22,28),ISTA(5),NO(5),(AMP1(J),
     5 EPOC1(J),J=29,35),ISTA(6),NO(6),(AMP1(J),EPOC1(J),J=36,37)
      Ind=INDEX(HEAD(2),'along')+5
      read(HEAD(2)(ind:ind+4),'(I4)')IPCD
      XMAJOR=IPCD
      READ (30,550)  HEAD(1),HEAD(2)
      READ (30,532)DATUM,ISTA(1),NO(1),(AMP2(J),EPOC2(J),J=1,7),
     1 ISTA(2),NO(2),(AMP2(J),EPOC2(J),J=8,14),ISTA(3),NO(3),
     2 (AMP2(J),EPOC2(J),J=15,21),ISTA(4),NO(4),
     3 (AMP2(J),EPOC2(J),J=22,28),ISTA(5),NO(5),(AMP2(J),
     4 EPOC2(J),J=29,35),ISTA(6),NO(6),(AMP2(J),EPOC2(J),J=36,37)

  532 FORMAT (F6.3,6(/2I4,7(F5.3,F4.1)))
  550 FORMAT (A80)
      N=0
      DO J=1,37
        PRINT 5000, ALIST(J),AMP1(J),EPOC1(J),AMP2(J),EPOC2(J)
        IF(AMP1(J) .GT. 0.0001)THEN
          N=N+1
          AMP1(n)=AMP1(J)
          AMP2(n)=AMP2(J)
          EPOC1(n)=EPOC1(j)
          EPOC2(n)=EPOC2(j)
        ENDIF
      ENDDO
      NCON=N
      CALL CELIPSE(NCON,ALIST,AMP2,AMP1,EPOC2,EPOC1,XMAJOR,0.)

 5000 FORMAT (1X,A10,2X,F7.3,3X,F6.2,2X,F7.3,3X,F6.2)
       
      end
      SUBROUTINE CELIPSE (NCON,LABLSV,D,SAVH,W,SAVK,AZI,STM)                                 
C                                                                               
C     DETERMINES AND PRINTS ELLIPSE PARAMETERS                                  
C     DATE: 01/31/2005
C     BY: RICHARD PATCHEN
C                                                                               
C     INPUT PARAMETERS:                                                         
C        NCON = THE NUMBER OF CONSTITUENTS TO PROCESS
C        LABLSV = THE CONSTITUENT NAME, I.E. M(2)
C                 (SEE SUBROUTINE NNAME FOR THE AVAILABLE
C                  CONSTITUENTS)
C        MAJOR OR NORTH AXIS --- SAVH: AMPLITUDE                                
C                                SAVK: PHASE (BASED ON STM)                     
C        MINOR OR EAST AXIS --- D: AMPLITUDE                                    
C                               W: PHASE (BASED ON STM)                         
C        AZI --- THE AXIS OF THE MAJOR COMPONENT of the harmonic constants (Principle Current Direction)                                 
C        TM = THE TIME MERIDIAN OF THE PHASES                                   
C                                                                               
c     HISTORY: BASED ON SUBROUTINE ELIPSE IN LSQ ANALYSIS
c
      PARAMETER (MX=180)                                                        
C                 *                                                              
C                 ******** THE MAXIMUM NO. OF ALLOW CONSTITUENTS                 
C                          CHECK THE DIMENSION OF <LABLEL> IN OTHER              
C                          SUBROUTINES                                           
C                                                                               
      DIMENSION D(*),SAVH(*),W(*),SAVK(*),C(5,2)                                
      DIMENSION A(MX),AP(MX),GSAVE(MX),HSAVE(MX),ISAVE(MX)                      
      REAL   SPEED                                                            
      CHARACTER*10   LABLSV(NCON),C*4                                                     
C
      PARAMETER (MXDRC=9000,NPR=3,MXDESC=3,MXFILES=1000)
C
      DATA(C(I,1),I=1,5)/4HCLOC,4HKWIS,4HE   ,4H    ,4H    /                    
      DATA(C(I,2),I=1,5)/4HCOUN,4HTER-,4HCLOC,4HKWIS,4HE   /                    
      LOUT = 6
      AZ = AZI*0.017453292                                                      
C
      DO 95444 I=1,NCON                                                         
      ISAVE(I)=0                                                                
      GSAVE(I)=0.0                                                              
      HSAVE(I)=0.0                                                              
95444 CONTINUE                                                                  
      WRITE (LOUT,40)                                                           
      WRITE(LOUT,42)STM                                                         
      AZJ =AZI + 90.                                                            
      IF(AZJ.GE.360.)AZJ=AZJ-360.                                               
      WRITE (LOUT,43)AZI,AZJ                                                    
      WRITE (LOUT,44)                                                           
      WRITE (LOUT,45)                                                           
C
      DO 88 I=1,NCON                                                            
      ZETA1 = W(I)*0.017453292                                                  
      ZETA2 = SAVK(I)*0.017453292                                               
      A1 = D(I)*COS(ZETA1)*COS(AZ)+SAVH(I)*COS(ZETA2)*SIN(AZ)                   
      B1 = D(I)*SIN(ZETA1)*COS(AZ)+SAVH(I)*SIN(ZETA2)*SIN(AZ)                   
      A2 = SAVH(I)*COS(ZETA2)*COS(AZ)-D(I)*COS(ZETA1)*SIN(AZ)                   
      B2 = SAVH(I)*SIN(ZETA2)*COS(AZ)-D(I)*SIN(ZETA1)*SIN(AZ)                   
      CC = A1*A1 + A2*A2                                                        
      DD = B1*B1 + B2*B2                                                        
      E  = A1*B1 + A2*B2                                                        
      F  = (CC-DD)/2.0                                                          
      G  = (CC+DD)/2.0                                                          
      IF (F.EQ.0.0.AND.E.EQ.0.0)   THEN
           AN = 0.0
      ELSE IF (F.NE.0.0)   THEN
           AN = ATAN2(E,F)                                              
      ELSE IF (F.EQ.0.0.AND.E.GT.0.0 )   THEN 
C
C          SET AN = PI/2
C
           AN =  1.570796                                 
      ELSE IF (F.EQ.0.0.AND.E.LT.0.0 )    THEN
C
C          SET AN = -PI/2
C
           AN = -1.570796                                 
      END IF
C
      IF (AN.LT.0.0)AN = 6.283185 + AN                                          
      IF(AN.NE.0.0)GO TO 15                                                     
      H = F/COS(AN)                                                             
      GO TO 14                                                                  
   15 H = E/SIN(AN)                                                             
   14 AN = AN/2.0                                                               
C                                                                               
      W1 = SQRT(G+H)                                                            
      TEST=G-H                                                                  
      IF(TEST.GT.0.0) GO TO 83                                                  
      W2=0.00001                                                                
      ISAVE(I)=1                                                               
      GSAVE(I)=G                                                               
      HSAVE(I)=H                                                               
      GO TO 84                                                                  
   83 W2 = SQRT(G-H)                                                            
   84 WS1 = A2*COS(AN) + B2  *SIN(AN)                                           
      WC1 = A1*COS(AN) + B1* SIN(AN)                                            
      THETA1 = ATAN2(WS1,WC1)*57.2957812                                        
      WS2 = B2*COS(AN) - A2 * SIN(AN)                                           
      WC2 = B1*COS(AN) - A1  *SIN(AN)                                           
      THETA2 = ATAN2(WS2,WC2)*57.2957812                                        
      IF(THETA1.LT.0.0)THETA1=THETA1+360.0                                      
      IF(THETA2.LT.0.0)THETA2=THETA2+360.0                                      
      TEST1 = THETA1 - 89.0                                                     
      IF(TEST1.LT.0.0) TEST1 = TEST1 +360.0                                     
      TEST1A= THETA1 - 91.0                                                     
      IF(TEST1A.LT.0.0)TEST1A= TEST1A+360.0                                     
      IROT=2                                                                    
      IF(TEST1.GE.THETA2.AND.TEST1A.LE.THETA2)IROT=1                            
      THETA1 = 90.0 - THETA1                                                    
      IF(THETA1.LT.0.0)THETA1 = THETA1 + 360.0                                  
      THETA2 = 90.0 - THETA2                                                    
      IF(THETA2.LT.0.0)THETA2 = THETA2 + 360.0                                  
      IF(IROT.EQ.2) THETA2=THETA2 + 180.                                        
      IF(THETA2.GE.360.) THETA2 = THETA2 - 360.                                 
      AN = AN*57.2957812                                                        
C                                                                               
C     **** RETRIEVE THE SPEED OF THE CONSTITUENT                                
C                                                                               
      CALL NNAME (SPEED,LABLSV(I),IDUM,IDUM,2)                                  
C                                                                               
      HOUR = AN/SPEED                                                           
      IF(IROT.EQ.1)ANN = AN + 90.0                                              
      IF(IROT.EQ.2) ANN = AN - 90.0                                             
      IF(ANN.GT.360.) ANN = ANN - 360.0                                         
      IF(ANN.LT.0.0)  ANN = ANN + 360.0                                         
      HOUR1 = ANN/SPEED                                                         
      ECC = SQRT(W1**2-W2**2)/W1                                                
      IF(IROT.EQ.1) ECC = -1.0*ECC                                                 
      WRITE (LOUT,66)  LABLSV(I),THETA1,W1,AN,HOUR,THETA2,W2,ANN,HOUR1,           
     1                 (C(J,IROT),J=1,3),ECC,SAVH(I),SAVK(I),D(I),W(I)             
   88 CONTINUE                                                                  
C                                                                               
      WRITE (LOUT,95)                                                           
      WRITE (LOUT,96)                                                           
      WRITE (LOUT,'(1X)')
      DO 85 I2=1,NCON
      IF(ISAVE(I2).EQ.0) GO TO 85                                               
      WRITE (LOUT,'(1X)')
      WRITE (LOUT,'(''ATTENTION -- AXIS RESET'')')                                      
      WRITE (LOUT,'(''CONSTITUENT'',2X,A10,2X,''W1 SET'',5X,''G='',
     1              F12.9,5X,''H='',F12.9)')  LABLSV(I2),GSAVE(I2),
     2              HSAVE(I2)
      WRITE (LOUT,'(1X)')
   85 CONTINUE                                                                  
      RETURN                                                                    
C                                                                               
   40 FORMAT(1X,'ELLIPSE PARAMETERS'/' (RIGHT-HANDED)')                        
   42 FORMAT( 80X,' TIME MERIDIAN OF RESULTS =',F5.0/)                          
   43 FORMAT(108X,F5.0,' AXIS',4X,F5.0,' AXIS')                                 
   44 FORMAT(22X,'MAJOR ELLIPSE AXIS',18X,'MINOR ELLIPSE AXIS',38X,'KAPP        
     *A',9X,'KAPPA')                                                            
   45 FORMAT(' CONSTITUENT DIR(TRUE) AMPLITUDE  PHASE  HOUR     DIR(TRUE        
     *) AMLITUDE  PHASE  HOUR     ROTATION',6X,'ECC',5X,'H(A)  PRIME   H        
     *(A)  PRIME'/)                                                             
   66 FORMAT(1X,A8 ,1X,2(3X,F7.3,2X,F7.3,2X,F7.3,1X,F7.1,2X),3A4,F8.2,2(        
     *F8.3,F6.1))                                                               
   95 FORMAT( /20X,'NOTE --- MINOR AXIS DIRECTION IS CONSIDERED AS THE M        
     1AJOR AXIS PLUS 90'/29X,'DEGREES IN ALL CASES.  IF ROTATION IS CLOC        
     2KWISE, PHASE WILL'/29X,'BE MAJOR AXIS PHASE PLUS 90 DEGREES.  IF C        
     3OUNTERCLOCKWISE, MINUS 90 DEGREES.')                                      
   96 FORMAT(  20X,'NOTE --- DIRECTION OF MAJOR AXIS IS USUALLY FLOOD DI        
     *RECTION' /29X,'BUT MAY BE EBB SOMETIMES. IF THIS OCCURS, TO GET TH        
     *E'/29X,'RESULTS FOR FLOOD DIRECTION, ADD 180 DEGREES TO ALL'/29X,'        
     *DIRECTIONS AND PHASES, AND (0.5 X CONSTITUENT PERIOD)TO ALL HOUR'/        
     *29X,'VALUES'/20X,'NOTE --- ECC = ECCENTRICITY OF ELLIPSE ( + IF CC        
     * ROTATION)')                                                              
      END                                                                       
      SUBROUTINE NNAME(SPDD,ITAG,ISUB,INUM,ICODE)                               
C     THIS SUBROUTINE IDENTIFIES THE CONSTITUENT BY ITS SPEED               1618
C     AND MAKES IT AVAILABLE FOR LABELING                                   1619
C     OR IDENTIFIES THE CONSTITUENT BY LABEL AND MAKES AVAILABLE ITS        1620
C     CONSTITUENT SPEED                                                     1621
C     IT ALSO DETERMINE THE SUBSCRIPT OF THE CONSTITUENT                    1622
C     SUBROUTINE PREPARED BY -- E. E. LONG  (NOAA/NOS) 1979                     
C        ORDER OF CONSTITUENT SPEEDS***  M(2),N(2),S(2),O(1),K(1),K(2)      1623
C      L(2),2N(2)R(2),T(2),LAMBDA(2),MU(2),NU(2),J(1),M(1),OO(1),P(1)       1624
C      Q(1),2Q(1),RHO(1),M(4),M(6),M(8),S(4),S(6),M(3),S(1),MK(3),2MK(3)    1625
C      MN(4),MS(4),2SM(2),MF,MSF,MM,SA,SSA                                  1626
C                                                                               
C     ICODE = 1 --- GIVEN A SPEED <SPDD>, RETURN A 10 CHARACTER LABEL           
C             2 --- GIVEN A 10 CHARACTER NAME <ITAG>, RETURN A                  
C                   CONSTITUENT SPEED                                           
C             3 --- GIVEN A CONSTITUENT NO., RETURN BOTH A 10 CHARACTER         
C                   NAME AND A CONSTITUENT SPEED                                
C                                                                               
      COMMON /IO/ LIN,LOUT                                                      
      REAL SPD,SPDD                                                           
      CHARACTER*10 LABLE(180),ITAG                                              
      DIMENSION SPD(180),IP(180)                                                
      DATA(IP(N),N = 1,37)/3*2,2*1,8*2,7*1,4,6,8,4,6,3,1,2*3,2*4,2,5*0/         
      DATA IP(38)/3/                                                            
      DATA(IP(N),N=39,114)/5*1,13*2,2*3,7*4,8*5,12*6,4*7,10*8,4*9,6*10,         
     111,4*12/                                                                  
      DATA(IP(N),N=115,120)/3,4,2*6,2*8/                                        
      DATA(IP(N),N=121,128)/4*1,2*2,3,4/                                    1633
      DATA(IP(N),N=129,140)/3*2,4,3*6,8,10,12,2*1/                          1634
      DATA(IP(N),N=141,150)/2*1,4*2,4*4/                                        
      DATA(IP(N),N=151,164)/2,11*0,3,5/                                         
      DATA(IP(N),N=165,175)/2*2,3*3,3*4,7,8,10/                                 
      DATA(SPD(I),I=  1, 37)/ 28.9841042, 28.4397295, 30.0000000,
     1    13.9430356,   15.0410686,   30.0821373,   29.5284789,       
     2    27.8953548,   30.0410667,   29.9589333,   29.4556253,       
     3    27.9682084,   28.5125831,   15.5854433,   14.4966939,       
     4    16.1391017,   14.9589314,   13.3986609,   12.8542862,       
     5    13.4715145,   57.9682084,   86.9523127,  115.9364169,       
     6    60.0000000,   90.0000000,   43.4761563,   15.0000000,       
     7    44.0251729,   42.9271398,   57.4238337,   58.9841042,       
     8    31.0158958,   01.0980331,   01.0158958,   00.5443747,       
     9    00.0410686,   00.0821373/                                   
      DATA SPD(38)/ 45.0410686/                                       
      DATA(SPD(I),I= 39, 77)/ 12.9271398, 14.0251729, 14.5695476,
     A    15.9748272,   16.0569644,   26.8794590,   26.9615963,       
     B    27.4238337,   27.5059710,   27.8039338,   28.6040041,       
     C    28.9019669,   29.0662415,   29.1483788,   29.3734880,       
     D    30.5443747,   30.7086493,   31.0980331,   42.3827651,       
     E    43.9430356,   56.8794590,   56.9523127,   57.5059710,       
     F    58.4397295,   58.5218668,   59.0662415,   59.5284789,       
     G    71.3668693,   71.9112440,   71.9933813,   72.4649023,       
     H    72.9271398,   73.0092770,   74.0251728,   74.1073100,       
     I    85.4013258,   85.8635632,   85.9457005,   86.4079380/       
      DATA(SPD(I),I= 78,116)/ 86.4900752, 87.4238337, 87.5059710,
     J    87.9682084,   88.0503457,   88.5218668,   88.9841042,       
     K    89.0662415,  100.3509735,  100.9046318,  101.9112440,       
     L   103.0092771,  114.8476674,  115.3920422,  115.4741794,       
     M   116.4079380,  116.4900752,  116.9523127,  117.0344500,       
     N   117.5059710,  117.9682084,  118.0503457,  129.8887360,       
     O   130.4331108,  130.9774855,  131.9933813,  144.3761464,       
     P   144.9205211,  145.3920422,  145.9364169,  146.4900752,       
     Q   146.9523127,  160.9774855,  174.3761464,  174.9205211,       
     R   175.4741794,  175.9364169,   41.9205276,   58.5125831/       
      DATA(SPD(I),I=117,140)/ 87.4966873, 88.5125831,117.4966873,
     S   116.4807916,   14.9178647,   15.0821353,   15.1232059,       
     T    15.5125897,   30.6265119,   27.3416965,   42.9271397,       
     U    60.0821373,   26.9523126,   27.4966873,   28.5947204,       
     V    57.4966873,   85.3920421,   85.9364168,   86.4807916,       
     W   115.4648958,  146.4807916,  175.4648958,   16.1391016,       
     X    12.8450026/                                           
      DATA(SPD(I),I=141,150)/ 15.1232058, 14.8767942, 30.0000001,         
     Y    29.9178627,   30.1642746,   29.9178666,   59.9589333,         
     Z    59.9178627,   60.2464119,   59.8767999/                         
      DATA(SPD(I),I=151,164)/ 28.9430356, 01.0569644, 00.5490165,         
     1          00.5079479, 00.0410667, 00.1232059, 00.1642746,         
     2          00.2464118, 00.3285841, 00.4106864, 00.4928237,         
     3          00.9856473, 45.0000000, 75.0000000/                       
      DATA(SPD(I),I=165,175)/ 27.8860712, 30.0410686, 43.4807981,         
     4          44.9589314, 45.1232059, 56.3258007, 56.8701754,         
     5          57.8860712,105.0000000,120.0000000,150.0000000/         
      DATA(LABLE(M),M = 1,37)/'M(2)'     ,'N(2)'     ,'S(2)'     ,              
     1          'O(1)'     ,'K(1)'     ,'K(2)'     ,'L(2)'     ,                
     2          '2N(2)'    ,'R(2)'     ,'T(2)'     ,'LAMBDA(2)',                
     3          'MU(2)'    ,'NU(2)'    ,'J(1)'     ,'M(1)'     ,                
     4          'OO(1)'    ,'P(1)'     ,'Q(1)'     ,'2Q(1)'    ,                
     5          'RHO(1)'   ,'M(4)'     ,'M(6)'     ,'M(8)'     ,                
     6          'S(4)'      ,'S(6)'      ,'M(3)'     ,'S(1)'     ,              
     7          'MK(3)'    ,'2MK(3)'   ,'MN(4)'    ,'MS(4)'    ,                
     8          '2SM(2)'    ,'MF'        ,'MSF'      ,'MM'        ,             
     9          'SA'        ,'SSA'       /                                      
      DATA LABLE(38)/'SK(3)'     /                                              
      DATA(LABLE(M),M =39,77)/'SIGMA(1)'  ,'MP(1)'     ,'CHI(1)'    ,           
     A          '2PO(1)'    ,'SO(1)'     ,'2NS(2)'    ,'2NK2S(2)'  ,            
     B          'MNS(2)'    ,'MNK2S(2)'  ,'2MS2K(2)'  ,'2KN2S(2)'  ,            
     C          'OP(2)'     ,'MKS(2)'    ,'M2(KS)(2)' ,'2SN(MK)(2)',            
     D          'MSN(2)'    ,'2KM(SN)(2)','SKM(2)'    ,'NO(3)'     ,            
     E          'SO(3)'     ,'N(4)'      ,'3MS(4)'    ,'MNKS(4)'   ,            
     F          'SN(4)'     ,'KN(4)'     ,'MK(4)'     ,'SL(4)'     ,            
     G          'MNO(5)'    ,'2MO(5)'    ,'3MP(5)'    ,'MNK(5)'    ,            
     H          '2MP(5)'    ,'2MK(5)'    ,'MSK(5)'    ,'3KM(5)'    ,            
     I          '3NKS(6)'   ,'2NM(6)'    ,'2NMKS(6)'  ,'2MN(6)'    /            
      DATA(LABLE(M),M=78,114)/'2MNKS(6)'  ,'MSN(6)'    ,'MKN(6)'    ,           
     J          '2MS(6)'    ,'2MK(6)'    ,'NSK(6)'    ,'2SM(6)'    ,            
     K          'MSK(6)'    ,'2MNO(7)'   ,'2NMK(7)'   ,'2MSO(7)'   ,            
     L          'MSKO(7)'   ,'2(MN)(8)'  ,'3MN(8)'    ,'3MNKS(8)'  ,            
     M          '2MSN(8)'   ,'2MNK(8)'   ,'3MS(8)'    ,'3MK(8)'    ,            
     N          'MSNK(8)'   ,'2(MS)(8)'  ,'2MSK(8)'   ,'2M2NK(9)'  ,            
     O          '3MNK(9)'   ,'4MK(9)'    ,'3MSK(9)'   ,'4MN(10)'   ,            
     P          'M(10)'     ,'3MNS(10)'  ,'4MS(10)'   ,'2MNSK(10)' ,            
     Q          '3M2S(10)'  ,'4MSK(11)'  ,'4MNS(12)'  ,'5MS(12)'   ,            
     R          '3MNKS(12)' ,'4M2S(12)'  /                                      
      DATA(LABLE(M),M=115,120)/'2NP(3)'    ,'ML(4)'     ,'2ML(6)'               
     1,          'MSL(6)'    ,'2MSL(8)'   ,'3ML(8)'    /                        
      DATA(LABLE(M),M=121,128)/'TK(1)'     ,'RP(1)'     ,'KP(1)'                
     S,          'THETA(1)'  ,'KJ(2)'     ,'OO(2)'     ,'MO(3)'                 
     T,          'SK(4)'     /                                             1701 
      DATA(LABLE(M),M=129,138)/'MLN2S(2)'  ,'2ML2S(2)'  ,'MKL2S(2)'             
     U,          '2MLS(4)'   ,'2NMLS(6)'  ,'2MLNS(6)'  ,'3MLS(6)'               
     V,          '4MLS(8)'   ,'3MSL(10)'  ,'4MSL(12)'  /                        
      DATA(LABLE(M),M=139,140)/'2KO(1)'    ,'2OK(1)'    /                 1705  
      DATA(LABLE(M),M=141,150)/'2KP(1)'   ,'2PK(1)'   ,'KP(2)'                  
     W,          '2SK(2)'   ,'2KS(2)'   ,'2TS(2)'   ,'ST(4)'                    
     X,          '3SK(4)'    ,'3KS(4)'    ,'3TS(4)'     /                       
      DATA(LABLE(M),M=151,164)/'SO(2)'    ,'SO(0)'    ,'.5MF'                   
     1,          '.5MSF'     ,'ST(0)'     ,'3SA'       ,'4SA'                   
     2,          '6SA'       ,'8SA'       ,'10SA'      ,'12SA'                  
     3,          '24SA'      ,'HS(3)'      ,'HS(5)'      /                      
      DATA(LABLE(M),M=165,175)/'O(2)'     ,'SK(2)'    ,'NK(3)'                  
     4,          'SP(3)'    ,'K(3)'     ,'NO(4)'    ,'MO(4)'                    
     5,         'SO(4)'     ,'S(7)'      ,'S(8)'      ,'S(10)'     /            
    1 FORMAT(1H1)                                                           1706
    2 FORMAT( 10X,' CONSTITUENT OF SPEED: ',F12.7,' NOT IN LIST') 
    3 FORMAT(/// 10X, 30H**** EXECUTION TERMINATED ****)                    1708
    4 FORMAT(10X,13HCONSTITUENT  , A10, 15HNOT IN THE LIST  )               1709
    5 FORMAT(10X,15HCONSTITUENT NO., I5, 13H  NOT IN LIST )                 1710
      GO TO (20,30,40),ICODE                                                1711
   20 DO 100 J = 1,175                                                      1712
      IF(SPDD.NE.SPD(J)) GO TO 100                                          1713
      ITAG = LABLE(J)                                                       1714
      ISUB = IP(J)                                                          1715
      INUM = J                                                              1716
      GO TO 101                                                             1717
  100 CONTINUE                                                              1718
      WRITE (LOUT,1)
      WRITE (LOUT,2) SPDD
      WRITE (LOUT,3)
      STOP
   30 DO 200 I = 1,175                                                      1723
      IF(ITAG.NE.LABLE(I)) GO TO 200                                        1724
      SPDD = SPD(I)                                                         1725
      ISUB = IP(I)                                                          1726
      INUM = I                                                              1727
      GO TO 101                                                             1728
  200 CONTINUE                                                              1729
      WRITE (LOUT,1)
      WRITE (LOUT,4) ITAG
      WRITE (LOUT,3)
      STOP
   40 DO 300 K = 1,175                                                      1734
      IF(INUM.NE.K) GO TO 300                                               1735
      ITAG = LABLE(K)                                                       1736
      SPDD = SPD(K)                                                         1737
      ISUB = IP(K)                                                          1738
      GO TO 101                                                             1739
  300 CONTINUE                                                              1740
      WRITE (LOUT,1)
      WRITE (LOUT,5) INUM 
      WRITE (LOUT,3)
      STOP
  101 RETURN                                                                1745
      END                                                                   1746
