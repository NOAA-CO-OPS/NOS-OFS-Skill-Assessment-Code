CC       This program will find maximum flood and ebb current speed,
C        start and end time of slack current before flood and ebb
C       t:  time in hours
C       h:  current speed
C       dirr: current direction in degrees
C       thighs:  time of maximum speed, obtained by calling subroutine extremes
C       hhighs:  maximum speed, obtained by calling subroutine extremes
C       idx:     =1 for maximum, and -1 for minimum, obtained by calling subroutine extremes
C       dirflood:  flood direction given by user
C      delt:     time interval of t()
C      NSMAX:    number of extremas from subroutine extremes
C      nmax:     number of data points

      SUBROUTINE SLACK_TIME(t,h,DIRR,thighs,hhighs,idx,NSMAX,nmax
     1  ,dirflood,delt,AFC,AEC,TFC,TEC,DFC,DEC,TSF,TEF,TSE,TEE,
     2  nafc,naec,NTSF,NTSE,NTSEE,nmx,nmx2)
      dimension t(nmx),h(nmx),dirr(nmx),thighs(nmx2)
     2  ,hhighs(nmx2),idx(nmx2)
     3  ,AFC(nmx2),AEC(nmx2),TFC(nmx2),TEC(nmx2),DFC(nmx2),DEC(nmx2)
     4  ,TSF(nmx2),TEF(nmx2),TSE(nmx2),TEE(nmx2)   

       nperhr=int(1./delt+0.001)
CC       fine Max current speed of flood and ebb
       nafc=0
       naec=0
       DO n=1,NSMAX
         IF (idx(n) .eq. 1)THEN
            NM=0
            DO I=1,nmax
              if (abs(thighs(n)-t(i)) .lt. 0.001)then
                 NM=I
                 goto 500
              endif
            ENDDO
500         CONTINUE
!        write(83,'(3f12.4,I3)')thighs(N),hhighs(N),DIRR(NM),idx(n)
            IF (NM .GT. 0)THEN
              IF( (abs(DIRR(NM)-dirflood) .LE. 90.) .OR. 
     1            (abs(DIRR(NM)-dirflood) .GT. 270.) )THEN
                nafc=nafc+1  
                AFC(nafc)=hhighs(N)
                TFC(nafc)=thighs(N)
                DFC(nafc)=DIRR(NM)
!        write(84,'(3f12.4)')thighs(N),hhighs(N),DIRR(NM)
              ELSE
                naec=naec+1  
                AEC(naec)=Hhighs(N)
                TEC(naec)=thighs(N)
                DEC(naec)=DIRR(NM)
!        write(85,'(3f12.4)')thighs(N),hhighs(N),DIRR(NM)
              ENDIF
            ELSE
              WRITE(6,*)'do not find match point, t=',thighs(n)
            ENDIF 
         ENDIF 
       ENDDO
       WRITE(6,*)'NAFC=',Nafc,'  NAEC=',Naec              
CC       find start and end time before flood and ebb around minimum within a window
       nwin=nperhr
       NTSF=1
       NTEF=1
       NTSE=1
       NTSEE=1
       NTEE=1
       DO N=1,NSMAX-1
!          write(6,*)'ntest=',n,'hhighs=',hhighs(n),idx(n),idx(n+1)
         IF ((idx(n) .EQ. -1) .and.(idx(n+1) .EQ. 1) )THEN
           IF (hhighs(n) .lt. 0.26)then
CC    find the slack position (index) in the original time series
            NM=0
            DO I=1,nmax
               if (abs(thighs(n)-t(i)) .lt. 0.1)then
                 NM=I
                 goto 550
               endif
            ENDDO
550         CONTINUE
            NSLACK=NM
CC  find next flood or ebb position (index) in the original time series
            NM=0
            DO I=1,nmax
               if (abs(thighs(n+1)-t(i)) .lt. 0.1)then
                 NM=I
                 goto 560
               endif
            ENDDO
560         CONTINUE
            NFLOOD=NM
!            write(6,*)'NSLACK=',NSLACK,'NFLOOD=',NFLOOD
CC          check if it is flood or ebb 
            IF( (abs(DIRR(NFLOOD)-dirflood) .LE. 90.) .OR. 
     1       (abs(DIRR(NFLOOD)-dirflood) .GT. 270.) )THEN
CC      find start and end slack time: slack is defined as if speed < 0.5 knots/s or 0.26 m/s
CC     find the point h>0.26 and point h<0.26, then using linear interpolation to get time where h=0.26 
              ifirst=0
              nwin=nperhr
570           hlocalmax=-1.0e5
              do nx=NSLACK,max0(1,NSLACK-nwin),-1
                   IF (h(nx) .gt. hlocalmax)then
                     hlocalmax=h(nx)
                     NBEGIN=nx
                   ENDIF
              enddo
              IF (hlocalmax .lt. 0.26)THEN
                 nwin=nwin+1
                 if (nwin .ge. NSLACK)goto 620
                 GOTO 570
              ENDIF
              do nx=max0(1,NBEGIN),min0(nmax,NSLACK)
                 if( (h(nx) .lt. 0.26 ) .and. (ifirst .eq. 0) )then
C  Zheng
                   if(n.ge.2) then
                     hh0=h(n-1)
                   else
                     hh0=0
                   end if
C  Zheng
                   hh1=h(nx)
                   tt0=t(nx-1)
                   tt1=t(nx)
                   t26=tt0+(0.26-hh0)*(tt0-tt1)/(hh0-hh1)   
!                   NTSF=NTSF+1
!                   tsf(NTSF)=t(nx)
                   tsf(NTSF)=t26
                   hstart=h(nx)
                   ifirst=1
                   goto 572
                 endif
              enddo
572           continue
              ifirst=0
              nwin=nperhr
580           hlocalmax=-1.0e5
              if (NSLACK+nwin .gt. NMAX)goto 620
              do nx=max0(1,NSLACK),min0(nmax,NSLACK+nwin)
                   IF (h(nx) .gt. hlocalmax)then
                     hlocalmax=h(nx)
                     NBEGIN=nx
                   ENDIF
              enddo
              iF (hlocalmax .lt. 0.26)THEN
                 nwin=nwin+1
                 GOTO 580
              ENDIF
              do nx=min0(nmax,NBEGIN),max0(1,NSLACK),-1
                 if( (h(nx) .lt. 0.26 ) .and. (ifirst .eq. 0) )then
                   hh0=h(nx+1)
                   hh1=h(nx)
                   tt0=t(nx+1)
                   tt1=t(nx)
                   t26=tt0+(0.26-hh0)*(tt0-tt1)/(hh0-hh1)   
!                   tef(NTEF)=t(nx)
                   tef(NTSF)=t26
                   hlast=h(nx)
                   NTSF=NTSF+1
                   ifirst=1
                   goto 582
                 endif
              enddo
582           continue  
!              write(33,'(4f12.4)')tsf(NTSF),tef(NTSF),0.26
            ELSE
              nwin=nperhr
              ifirst=0
590           hlocalmax=-1.0e5
              if (nwin .ge. NSLACK)goto 620
              do nx=min0(nmax,NSLACK),max0(1,NSLACK-nwin),-1
                   IF (h(nx) .gt. hlocalmax)then
                     hlocalmax=h(nx)
                     NBEGIN=nx
                   ENDIF
              enddo
              iF (hlocalmax .lt. 0.26)THEN
                 nwin=nwin+1
                 GOTO 590
              ENDIF
              do nx=max0(1,NBEGIN),min0(nmax,NSLACK)
                 if( (h(nx) .lt. 0.26 ) .and. (ifirst .eq. 0) )then
                   if(nx.gt.1) then
                    hh0=h(nx-1)
                    tt0=t(nx-1)
                   else
                    hh0=0.0
                    tt0=0.0
                   end if
                   hh1=h(nx)
                   tt1=t(nx)
                   t26=tt0+(0.26-hh0)*(tt0-tt1)/(hh0-hh1)
C                   WRITE(*,*) 'NTSEE=',NTSEE,'  TSE=',t26
                   tse(NTSEE)=t26
                   hstart=h(nx)
                   ifirst=1
                   NTSEE=NTSEE+1
                   goto 592
                 endif
              enddo
592           nwin=nperhr
              ifirst=0
595           hlocalmax=-1.0e5
              if (NSLACK+nwin .gt. NMAX)goto 620
              do nx=max0(1,NSLACK),min0(nmax,NSLACK+nwin)
                   IF (h(nx) .gt. hlocalmax)then
                     hlocalmax=h(nx)
                     NBEGIN=nx
                   ENDIF
              enddo
              iF (hlocalmax .lt. 0.26)THEN
                 nwin=nwin+1
                 GOTO 595
              ENDIF
              do nx=min0(nmax,NBEGIN),max0(1,NSLACK),-1
                 if( (h(nx) .lt. 0.26 ) .and. (ifirst .eq. 0) )then
                   hh0=h(nx+1)
                   hh1=h(nx)
                   tt0=t(nx+1)
                   tt1=t(nx)
                   t26=tt0+(0.26-hh0)*(tt0-tt1)/(hh0-hh1)   
C                   WRITE(*,*) 'NTSE=',NTSE,'  TEE=',t26
                   tee(NTSE)=t26
                   hlast=h(nx)
                   ifirst=1
                   NTSE=NTSE+1
!                   write(6,'(5f10.4)')tt0,hh0,tt1,hh1,t26
                   goto 600
                 endif
              enddo
600           continue
!              write(34,'(4f12.4)')tse(NTSE),tee(NTSE),0.26
            ENDIF               
           ENDIF               
         ENDIF
620     continue
!       WRITE(6,*)'NTSF=',NTSF,'  NTSE=',NTSE              
       ENDDO
       NTSF=NTSF-1
       NTSE=NTSE-1
       NTSEE=NTSEE-1
       WRITE(6,*) 'NTSF=',NTSF,'  NTSE=',NTSE,'  NTSEE=',NTSEE            
!       open(10,file='test_afc.dat')
!       open(20,file='test_aec.dat')
!       do i=1,nafc
!       write(10,"(f12.5,2f8.3)")tFC(i),AFC(i),DFC(i)
!       enddo
!       do i=1,naec
!       write(20,"(f12.5,2f8.3)")TEC(i),AEC(i),DEC(i)
!       enddo
!       close(10)
!       close(20)
!       open(10,file='test_tf.dat')
!       open(20,file='test_te.dat')
!       do i=1,nTSF
!       write(10,"(2f12.5,2f8.3)")tSF(i),TEF(i)
!       enddo
!       do i=1,nTSE
!       write(20,"(2f12.5,2f8.3)")TSE(i),TEE(i)
!       enddo
!       close(10)
!       close(20)
      write(6,*)'end of slack.f'
      return
      end

