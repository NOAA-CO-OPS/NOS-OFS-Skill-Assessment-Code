C  Write out a table
      Subroutine table_nontide(fileshort,stationname,KINDAT)
      include 'skills.inc'
      character*120 stationname,fileshort
      character*20 cmethod(2),cgap(2),cbdate,cedate
      character title(6)*36,c1*1,c2*1,ttl*9,ttl2*15
      dimension b(imx),e(imx),r(imx),p(imx),tz(imx),
     1          ifrst(imx),ilast(imx),
     2          hhwp0(nmx),thwp0(nmx),hlwp0(nmx),tlwp0(nmx),
     3   thighsr(nmx),hhighsr(nmx),idxr(nmx),
     4   thighsp(nmx),hhighsp(nmx),idxp(nmx),
     5   thighs0(nmx),hhighs0(nmx),idx0(nmx)
      dimension dire(imx),dirr(imx),dirp(imx),dirop(imx),dirap(imx),
     1  AFCR(nmx),AECR(nmx),TFCR(nmx),TECR(nmx),DFCR(nmx),DECR(nmx),
     2  TSFR(nmx),TEFR(nmx),TSER(nmx),TEER(nmx),   
     3  AFCP(nmx),AECP(nmx),TFCP(nmx),TECP(nmx),DFCP(nmx),DECP(nmx),
     4  TSFP(nmx),TEFP(nmx),TSEP(nmx),TEEP(nmx)   
      dimension op(imx),ap(imx),tp(imx),p2(imx),xtmp(imx),
     1  ytmp(imx),xtmp1(imx),ytmp1(imx)
      real*8 jday,jday0,jday1,jbase_date
      real*8 JULIAN,yearb,monthb,dayb,hourb

      data title/'SCENARIO: TIDAL SIMULATION ONLY     ',
     1           'SCENARIO: HINDCAST                  ',
     2           'SCENARIO: SEMI-OPERATIONAL NOWCAST  ',
     3           'SCENARIO: SEMI-OPERATIONAL FORECAST ',
     4           'COMPARISON: PERSISTENCE FORECAST    ',
     5           'COMPARISON: ASTRONOMICAL TIDE ONLY  '/
      data cmethod/'Cubic Spline','SVD         '/
      data cgap/'Gap is not filled','Gap is filled    '/ 

      WRITE(*,*) TRIM(FILESHORT),' ',TRIM(STATIONNAME),KINDAT
      Write(*,"(/,' <tab> ')")
C  Calculate start and end date
      ii = 0
      do i = 1,imaxB
        IF(of(i) .GT. -900.0) then
          II = II+1
          tp(ii) = t(i)
        ENDIF
      enddo

      nmax1 = II
      gap = 1.5*delt
      call continuous(tp,nmax1,gap,Nsegments,ifrst,ilast)
      difm = -999.0
      Write(*,*) 'Nsegments = ',Nsegments
      do NG = 1,Nsegments
        Istart = ifrst(NG)
        IEND = ilast(NG)
        DIF = tp(IEND)-tp(Istart)
        IF(DIF .gt. difm) THEN
          difm = dif
          II = NG
        ENDIF  
      ENDDO
      Write(*,*) 'II=',II,'the longest segment=',difm, ' hours'

      IYR = iyear
      Istart = ifrst(II)
      IEND = ilast(II)
      yearb = iyear
      monthb = 1.0
      dayb = 1.0
      hourb = 0.0
      jday0 = JULIAN(yearb,monthb,dayb,hourb)
      jday = tp(Istart)/24.0+jday0-1     !! tp start from 1 on Jan. 1
      call GREGORIAN(jday,yearb,monthb,dayb,hourb)

      IYRS = INT(yearb)
      ICMS = int(monthb+0.001)
      ICDS = INT(dayb+0.001)
      IHRS = INT(hourb+0.001)
      IMNS = INT((hourb-IHRS)*60+0.1)

      jday = tp(IEND)/24.0+jday0-1
      call GREGORIAN(jday,yearb,monthb,dayb,hourb)
      IYRE = INT(yearb)
      ICME = int(monthb+0.001)
      ICDE = INT(dayb+0.001)
      IHRE = INT(hourb+0.001)
      IMNE = INT((hourb-IHRE)*60+0.1)

      lu = 10
      lu1 = 12
      open(lu,file=trim(fileshort)//'_table.out',form='formatted')
C  Write column header
      write(*,"(/,'Station: ',a40)")trim(stationname)
      write(*,"(/,'Variable   X   N   Imax    SM    RMSE    SD     NOF'
     1  ,'   CF    POF   MDNO  MDPO  WOF   CORR  SKILL',
     2  /'Criterion  -   -    -      -  ',
     3  '    -       -     <1%  >90%   <1%    <N    <N  <.5%')")
      write(lu,"(/,'Station: ',a40)") trim(stationname)
      write(lu,"('Observed data time period from: /'
     1  I2,'/',I2,'/',I4,'  to /',
     2  I2,'/',I2,'/',I4)") ICMS,ICDS,IYRS,ICME,ICDE,IYRE

      IF(KINDAT .EQ. 1) THEN
        open(lu1,file=trim(fileshort)//'phase_table.out')
        write(lu1,"(/,'Station: ',a40)")trim(stationname)
        write(lu1,"('Observed data time period from: /',
     1     I2,'/',I2,'/',I4,'  to /',
     2     I2,'/',I2,'/',I4)") ICMS,ICDS,IYRS,ICME,ICDE,IYRE
      ENDIF

      IF(Igapfill .eq. 0) then
        write(lu,"('Data gap is not filled')")
        IF(KINDAT .EQ. 1) THEN
          write(lu1,"('Data gap is not filled')")
        endif
      ELSE
        IF(method .eq. 0)then
          write(lu,"('Data gap is filled using cubic spline method')")
          IF(KINDAT .EQ. 1) THEN
          write(lu1,"('Data gap is filled using cubic spline method')")
          endif
        else if (method .eq. 1)then
          IF(KINDAT .EQ. 1) THEN
            write(lu1,"('Data gap is filled using SVD method')")
          endif
          write(lu,"('Data gap is filled using SVD method')")
        endif
      ENDIF  

      IF(tcut .gt. 0.0) THEN
        write(lu,"('Data are filtered using ',F5.1,
     1    ' Hour Fourier Filter')") tcut
        IF(KINDAT .EQ. 1) THEN
          write(lu1,"('Data are filtered using ',F5.1,
     1      ' Hour Fourier Filter')") tcut
        endif
      else
        write(lu,"('Data are not filtered')")
        IF(KINDAT .EQ. 1) THEN
          write(lu1,"('Data are not filtered')")
        endif
      endif

      write(lu,"('------------------------------------------',
     1 '-------------------------------------------------')")
      write(lu,"('VARIABLE    X     N   IMAX    SM    RMSE    SD',
     1  '     NOF   CF    POF   MDNO  MDPO  WOF   CORR  SKILL',
     2  /'CRITERION   -     -     -      -  ',
     3  '    -      -     <1%  >90%   <1%    <N    <N  <.5%')")
      write(lu,"('------------------------------------------',
     1  '-----------------------------------------',/)")

      IF(KINDAT .EQ. 1) THEN
        write(lu1,"('------------------------------------------',
     1    '-------------------------------------------------')")
        write(lu1,"('VARIABLE    X     N   IMAX    SM    RMSE    SD',
     1    '     NOF   CF    POF   MDNO  MDPO  WOF   CORR   SKILL',
     2    /'CRITERION   -     -     -      -  ',
     3    '    -      -     <1%  >90%   <1%    <N    <N  <.5%')")
        write(lu1,"('------------------------------------------',
     1    '-------------------------------------------------',/)")
      endif

C  Loop thru scenarios
      do 220 is = 1,6
        if(ISWITCH(is) .eq. 0) goto 220
C  Write header
        write(*,"(/,5x,a36)") title(is)
        write(lu,"(5x,a36)") title(is)
        IF(nmax1 .LT. 2*IPDAY) THEN
          write(lu,"(5x,a36)")'no observational data'
          goto 220
        ENDIF
        IF(KINDAT .EQ. 1) THEN
          write(lu1,"(5x,a36)")title(is)
        endif

C  Compute error/difference series
        if(is .eq. 4 .or. is .eq. 5) goto 130
        IF(is .eq. 1) THEN  !!! tidal simulation only
          nmax = imaxB
          II = 0
          do i = 1,nmax
            IF((a(i) .GT. -900.0) .and. (s(i) .GT. -900.0)) then
              II = II+1
              r(ii) = a(i)
              p(ii) = s(i)
              ytmp(ii) = of(i)
              e(ii) = p(ii)-r(ii)
              IF(KINDAT .EQ. 1) THEN
                IF((a(i) .GT. 0.26) .and. (s(i) .GT. 0.26) ) then
                  dirr(ii) = dira(i)
                  dirp(ii) = dirs(i)
                  dire(ii) = dirp(ii)-dirr(ii)
                  if(dire(ii) .GT. 180.0) dire(ii) = dire(ii)-360.0
                  if(dire(ii) .LT. -180.0) dire(ii) = dire(ii)+360.0
                ELSE 
                  dirr(ii) = dira(i)
                  dirp(ii) = dirs(i)
                  dire(ii) = 0.0
                ENDIF
              ENDIF                
              tp(ii) = t(i)
            ENDIF
            if(i.le.10) write(*,"('  is=1  i=',i5,' t=',f11.4,'  p=',
     1      f9.4,' r=',f9.4,'  e=',f9.4)")ii,tp(ii),p(ii),r(ii),e(ii)
          enddo
          nmax = II
          wcof = wof(p,r,r,imx,nmax,2.0*X1)

        ELSE IF(is .eq. 2) THEN  !!! model hindcast
          nmax = imaxb
          II = 0
          do i = 1,nmax
            IF((of(i) .GT. -900.0) .and. (hi(i) .GT. -900.0)) then
              II = II+1
              r(ii) = of(i)
              p(ii) = hi(i)
              e(ii) = p(ii)-r(ii)
              tp(ii) = t(i)
              IF(KINDAT .EQ. 1) THEN
                IF((of(i) .GT. 0.26) .and. (hi(i) .GT. 0.26)) then
                  dirr(ii) = diro(i)
                  dirp(ii) = dirhi(i)
                  dire(ii) = dirp(ii)-dirr(ii)
                  if(dire(ii) .GT. 180.0) dire(ii) = dire(ii)-360.0
                  if(dire(ii) .LT. -180.0) dire(ii) = dire(ii)+360.0
                ELSE 
                  dirr(ii) = diro(i)
                  dirp(ii) = dirhi(i)
                  dire(ii) = 0.0
                ENDIF
              ENDIF                
            ENDIF
          enddo
          nmax = II

        ELSE IF(is .eq. 3)THEN  !!! model nowcast
          nmax = imaxb
          Write(*,*) 'first nowcast=',tc(1),c(1)
          II = 0
          do i = 1,nmax
            IF((of(i) .GT. -900.0) .and. (c(i) .GT. -900.0)) then
              II=II+1
              r(ii)=of(i)
              p(ii)=c(i)
              e(ii)=p(ii)-r(ii)
              tp(ii)=t(i)
              IF(KINDAT .EQ. 1) THEN
                IF((of(i) .GT. 0.26) .and. (c(i) .GT. 0.26)) then
                  dirr(ii) = diro(i)
                  dirp(ii) = dirc(i)
                  dire(ii) = dirp(ii)-dirr(ii)
                  if(dire(ii) .GT. 180.0) dire(ii) = dire(ii)-360.0
                  if(dire(ii) .LT. -180.0) dire(ii) = dire(ii)+360.0
                ELSE 
                  dirr(ii) = diro(i)
                  dirp(ii) = dirc(i)
                  dire(ii) = 0.0
                ENDIF
              ENDIF                
            ENDIF
          enddo
          nmax = II

        ELSE IF(is .eq. 6) THEN
          nmax = imaxb
          II = 0
          do i = 1,nmax
            IF(of(i) .GT. -900.0) then
              II = II+1
              r(ii) = of(i)
              p(ii) = a(i)
              e(ii) = p(ii)-r(ii)
              tp(ii) = t(i)
              IF(KINDAT .EQ. 1) THEN
                IF((of(i) .GT. 0.26) .and. (a(i) .GT. 0.26)) then
                  dirr(ii) = diro(i)
                  dirp(ii) = dira(i)
                  dire(ii) = dirp(ii)-dirr(ii)
                  if(dire(ii) .GT. 180.0) dire(ii) = dire(ii)-360.0
                  if(dire(ii) .LT. -180.0) dire(ii) = dire(ii)+360.0
                ELSE 
                  dirr(ii) = diro(i)
                  dirp(ii) = dira(i)
                  dire(ii) = 0.0
                ENDIF
              ENDIF                
            ENDIF
          enddo
          nmax = II
          wcof = wof(p,p,r,imx,nmax,2.*X1)
        ENDIF

C  Create series
        if(is .lt. 5) then   !!! for is=1:tide only; =2:hindcast; =3:nowcast
          ttl2 = '         '
          gap = 1.5*delt
          if(KINDAT .eq. 1) THEN
            call prtline1('U        ',ttl2,lu,p,tp,imx,nmax,X1,1,gap)
            call prtline1('u        ',ttl2,lu,r,tp,imx,nmax,X1,1,gap)
            call prtline1('D        ',ttl2,lu1,dirp,tp,imx,nmax,X11,
     1        1,gap)
            call prtline1('d        ',ttl2,lu1,dirr,tp,imx,nmax,X11,
     1        1,gap)

          else if(KINDAT .eq. 2) THEN
            call prtline1('H        ',ttl2,lu,p,tp,imx,nmax,X1,1,gap)
            call prtline1('h        ',ttl2,lu,r,tp,imx,nmax,X1,1,gap)
          else if(KINDAT .eq. 3) THEN
            call prtline1('T        ',ttl2,lu,p,tp,imx,nmax,X1,1,gap)
            call prtline1('t        ',ttl2,lu,r,tp,imx,nmax,X1,1,gap)
          else if(KINDAT .eq. 4) THEN
            call prtline1('S        ',ttl2,lu,p,tp,imx,nmax,X1,1,gap)
            call prtline1('s        ',ttl2,lu,r,tp,imx,nmax,X1,1,gap)
          endif
        endif
        imax = nmax
        gap = 1.5*delt
        if(KINDAT .eq. 1) THEN
          skillv = skilla(p,r,imx,nmax)
          CORR_C = CORRELATION(P,R,IMX,NMAX)
          IX1 = INT(X1*100+0.1)
          Write(ttl2,'(a3,I3,a9)')'   ',IX1,' cm/s 24h'
          Call prtline1('U-u      ',ttl2,lu,e,tp,imx,nmax,X1,6,gap)
          ttl2 = ' 22.5 dg 24h'
          Write(ttl2,'(a2,F4.1,a9)') '  ',X11,'    dg 24h'
          skillv = skilla(dirp,dirr,imx,nmax)
          CORR_C=CORRELATION(DIRP,DIRR,IMX,NMAX)
          call prtline1('D-d      ',ttl2,lu1,dire,tp,imx,nmax,X11,
     1      6,gap)

        else if(KINDAT .eq. 2) THEN
          IX1 = INT(X1*100+0.1)
          Write(ttl2,'(a3,I3,a9)') '   ',IX1,'   cm 24h'
          skillv = skilla(p,r,imx,nmax)
          CORR_C = CORRELATION(P,R,IMX,NMAX)
          Call prtline1('H-h      ',ttl2,lu,e,tp,imx,nmax,X1,5,gap)
        else if(KINDAT .eq. 3) THEN
          Write(ttl2,'(a3,F3.1,a9)') '   ',X1,'    c 24h'
          skillv = skilla(p,r,imx,nmax)
          CORR_C = CORRELATION(P,R,IMX,NMAX)
          Call prtline1('T-t      ',ttl2,lu,e,tp,imx,nmax,X1,6,gap)
        else if(KINDAT .eq. 4) THEN
          write(ttl2,'(a3,F3.1,a9)') '   ',X1,'  psu 24h'
          skillv = skilla(p,r,imx,nmax)
          CORR_C = CORRELATION(P,R,IMX,NMAX)
          Call prtline1('S-s      ',ttl2,lu,e,tp,imx,nmax,X1,6,gap)
        endif 
        goto 170

 130  continue
C  Add on Dec. 13, 2004 
        if(is .eq. 4) indx1 = 2  ! fcst
        if(is .eq. 5) indx1 = 3  ! persistence
        NLOOP = INT(NFDURATION/6+0.001) 
        NLOOPS = 1
        if(FINCLUDE0) NLOOPS = 0

C  Extend forecast hours to 30z, 36z,etc. for the skill assessment 
        do 160 jj = 0,NLOOP  !! do 00z,06z, 12z, 18z and 24z
          Call addpred(indx1,2,ytmp,dirp,xtmp,ipmax,jj)
          Write(*,*) 'ipmax in addpred=',ipmax
          Call collect_obs(op,ap,dirop,dirap,xtmp,tz,ipmax,nmax)
          if(nmax. ne. ipmax) then
            write(*,*) ' nmax=',nmax,' is unequal to ipmax=',ipmax
            DO I = 1,nmax
              KK = 0
              DO i1=1,IPMAX
                IF(abs(tz(i)-xtmp(i1)) .le. 0.001) then
                  ytmp(i) = ytmp(i1)
                  xtmp(i) = xtmp(i1)
                  KK = I1
                  goto 135
                ENDIF
              ENDDO

135           Continue
              IF(KK .EQ. 0) THEN
                WRITE(*,*) 'No match data found in obs. at time ',
     1            tz(i),jj
                WRITE(*,*) 'program stop at table.f!'
                stop
              endif
            ENDDO
            IPMAX = nmax
          endif

          II = 0
          do i = 1,nmax
            IF((op(i) .GT. -900.0) .and. (ytmp(i) .GT. -900.0)) then
              II = II+1
              r(ii) = op(i)
              p(ii) = ytmp(i)
              tp(ii) = xtmp(i)
              e(ii) = p(ii)-r(ii)
              IF(KINDAT .EQ. 1) THEN
                IF((op(i) .GT. 0.26) .and. (ytmp(i) .GT. 0.26)) then
                  dirp(ii) = dirp(i)
                  dirop(ii) = dirop(i)
                  dirap(ii) = dirap(i)
                  dire(ii) = dirp(ii)-dirop(ii)
                  if(dire(ii) .GT. 180.0) dire(ii) = dire(ii)-360.0
                  if(dire(ii) .LT. -180.0) dire(ii) = dire(ii)+360.0
                ELSE 
                  dirp(ii) = dirp(i)
                  dirop(ii) = dirop(i)
                  dirap(ii) = dirap(i)
                  dire(ii) = 0.0
                ENDIF
              ENDIF                
    
              ap(ii) = ap(i)
            endif
          enddo

          nmax = II
          c1 = char(ichar('0')+(jj*6)/10)
          c2 = char(ichar('0')+(jj*6)-10*((jj*6)/10))
          if(KINDAT .eq. 1) THEN
            IX1 = INT(X1*100+0.1)
            Write(ttl2,'(a3,I3,a9)') '   ',IX1,' cm/s 24h'
            ttl = 'U'//c1//c2//'-'//'u'//c1//c2
            Write(ttl,3234) 'U',jj*6,'-u',jj*6

            wcof = wof(p,ap,r,imx,nmax,2.*X1)
            skillv = skilla(p,r,imx,nmax)
            CORR_C = CORRELATION(P,R,IMX,NMAX)
            gap = 25.0
            call prtline1(ttl,ttl2,lu,e,tp,imx,nmax,X1,3,gap)
            write(ttl2,'(a2,F4.1,a9)') '  ',X11,'   dg 24h'
            ttl = 'D'//c1//c2//'-'//'d'//c1//c2
            write(ttl,3234) 'D',jj*6,'-d',jj*6
            skillv = skilla(dirp,dirop,imx,nmax)
            CORR_C=CORRELATION(DIRP,DIRop,IMX,NMAX)
            call prtline1(ttl,ttl2,lu1,dire,tp,imx,nmax,X1,3,gap)

          else if(KINDAT .eq. 2) THEN
            IX1 = INT(X1*100+0.1)
            write(ttl2,'(a3,I3,a9)') '   ',IX1,'   cm 24h'
            ttl = 'H'//c1//c2//'-'//'h'//c1//c2
            write(ttl,3234) 'H',jj*6,'-h',jj*6
3234	    format(a1,I3.3,a2,I3.3)
            wcof = wof(p,ap,r,imx,nmax,2.*X1)
            skillv = skilla(p,r,imx,nmax)
            CORR_C = CORRELATION(P,R,IMX,NMAX)
            gap = 25.0
            call prtline1(ttl,ttl2,lu,e,tp,imx,nmax,X1,3,gap)

          else if(KINDAT .eq. 3) THEN
            write(ttl2,'(a3,F3.1,a9)') '   ',X1,'    c 24h'
            ttl = 'T'//c1//c2//'-'//'t'//c1//c2
            write(ttl,3234) 'T',jj*6,'-t',jj*6
            skillv = skilla(p,r,imx,nmax)
            CORR_C = CORRELATION(P,R,IMX,NMAX)
            gap = 25.0
            call prtline1(ttl,ttl2,lu,e,tp,imx,nmax,X1,3,gap)

          else if(KINDAT .eq. 4) THEN
            write(ttl2,'(a3,F3.1,a9)') '   ',X1,'  psu 24h'
            ttl = 'S'//c1//c2//'-'//'s'//c1//c2
            write(ttl,3234) 'S',jj*6,'-s',jj*6
            skillv = skilla(p,r,imx,nmax)
            CORR_C = CORRELATION(P,R,IMX,NMAX)
            gap = 25.0
            call prtline1(ttl,ttl2,lu,e,tp,imx,nmax,X1,3,gap)
          endif
 160    continue

        if((KINDAT .eq. 3) .or. (KINDAT .eq. 4)) GOTO 220
        Call addpred(indx1,1,ytmp,dirp,xtmp,ipmax,0)
        Write(*,*) 'ipmax in addpred=',ipmax
        Call collect_obs(op,ap,dirop,dirap,xtmp,tz,ipmax,nmax)

        if(nmax .ne. ipmax) then
          write(*,*) ' nmax=',nmax,' is unequal to ipmax=',ipmax
          write(*,*) 'program stop here!'
          DO I = 1,nmax
            DO i1 = 1,IPMAX
              IF(abs(tz(i)-xtmp(i1)) .le. 0.001) then
                ytmp(i) = ytmp(i1)
                xtmp(i) = xtmp(i1)
                dirp(i) = dirp(i1)
                goto 165
              ENDIF
            ENDDO
165         continue
          ENDDO  
        endif

        ii = 0
        do i = 1,nmax
          IF((op(i) .GT. -900.0) .and. (ytmp(i) .GT. -900.0)) then
            II = II+1
            r(ii) = op(i)
            p(ii) = ytmp(i)
            tp(ii) = xtmp(i)
            e(ii) = p(ii)-r(ii)
            IF(KINDAT .EQ. 1) THEN
              IF((op(i) .GT. 0.26) .and. (ytmp(i) .GT. 0.26)) then
                dirr(ii) = dirop(i)
                dirp(ii) = dirp(i)
                dire(ii) = dirp(ii)-dirr(ii)
                if(dire(ii) .GT. 180.0) dire(ii) = dire(ii)-360.0
                if(dire(ii) .LT. -180.0) dire(ii) = dire(ii)+360.0
              ELSE 
                dirr(ii) = dirop(i)
                dirp(ii) = dirp(i)
                dire(ii) = 0.0
              ENDIF
            ENDIF                
          endif
        enddo
        nmax=ii

170     continue
        if((KINDAT .eq. 3) .or. (KINDAT .eq. 4)) GOTO 220
        sm0 = sm(r,imx,nmax)
        sd0 = sd(r,imx,nmax)
        rmse0 = rmse(r,imx,nmax)
        hupper = sm0+2*sd0    
        hlower = sm0-2*sd0
C  The hupper and hlower are re-computed inside the subroutine of 
C  events according to the parameter "factor"
        WRITE(*,*) 'mean=',sm0,'SD=',sd0,'RMSE=',rmse0
        WRITE(*,*) 'Lower and Upper limits =',hlower,hupper
C  Using continuous time series to get extrema
        ipr = 0
        if(is .eq. 0) ipr = 1
        write(*,"(/,'   reference series')")
        gap = 1.5*delt
        call continuous(tp,nmax,gap,Nsegments,ifrst,ilast)
        write(*,*) 'Nsegments=',Nsegments

C        rmdo = 0
        nmaxrh9 = 0
        nmaxrl9 = 0
        nmaxph9 = 0
        nmaxpl9 = 0
        nsmaxr = 0
        nsmaxp = 0
        DO NG = 1,Nsegments
          Istart = ifrst(NG)
          IEND = ilast(NG)
          DIF = tp(IEND)-tp(Istart)
          IF(DIF .gt. 48.0) THEN
            numb = IEND-Istart+1  
            DO I = Istart,IEND
              I0 = I-istart+1 
              xtmp(i0) = tp(i)
              ytmp(i0) = r(i)
            ENDDO

            Call extremes(xtmp,ytmp,numb,ipr1,delhr,delamp,delpct,
     1        iopta,thighs0,hhighs0,idx0,nsmax0,tcut,delt,hhwp0,
     2        thwp0,nmaxrh0,hlwp0,tlwp0,nmaxrl0,NTYPE)

            DO J = 1,nsmax0
              nsmaxr = nsmaxr+1
              thighsr(nsmaxr) = thighs0(j)
              hhighsr(nsmaxr) = hhighs0(j)
              idxr(nsmaxr) = idx0(j)
              IF(idx0(j) .LT. 0) THEN
                nmaxrl9 = nmaxrl9+1
                hlwr(nmaxrl9) = hhighs0(j)
                tlwr(nmaxrl9) = thighs0(j)
              ELSEIF(idx0(j) .GT. 0) THEN
                nmaxrh9 = nmaxrh9+1
                hhwr(nmaxrh9) = hhighs0(j)
                thwr(nmaxrh9) = thighs0(j)
              ENDIF
            ENDDO

C  For prediction
            DO I = Istart,IEND
              I0 = I-istart+1 
              ytmp(i0) = p(i)
            ENDDO

            Call extremes(xtmp,ytmp,numb,ipr1,delhr,delamp,delpct,
     1        iopta,thighs0,hhighs0,idx0,nsmax0,tcut,delt,hhwp0,
     2        thwp0,nmaxph0,hlwp0,tlwp0,nmaxpl0,NTYPE)

            DO J = 1,nsmax0
              nsmaxp = nsmaxp+1
              thighsp(nsmaxp) = thighs0(j)
              hhighsp(nsmaxp) = hhighs0(j)
              idxp(nsmaxp) = idx0(j)
              IF(idx0(j) .LT. 0) THEN
                nmaxpl9 = nmaxpl9+1
                tlwp(nmaxpl9) = thighs0(j)
                hlwp(nmaxpl9) = hhighs0(j)
              ELSEIF(idx0(j) .GT. 0) THEN
                nmaxph9 = nmaxph9+1
                thwp(nmaxph9) = thighs0(j)
                hhwp(nmaxph9) = hhighs0(j)
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        nmaxrh = nmaxrh9
        nmaxrl = nmaxrl9
        nmaxph = nmaxph9
        nmaxpl = nmaxpl9

C  Select events from peaks
C  Calculate standard deviation within 7 days = 168hours
        CALL events(tp,r,nmax,thighsr,hhighsr,idxr,nsmaxr,
     1    hupper,hlower,FACTOR)
        call collect_events_prd(thighsr,hhighsr,idxr,nsmaxr,
     1    p,tp,nmax,thighsp,hhighsp,idxp,nsmaxp)
        nmaxrh9 = 0
        nmaxrl9 = 0
        nmaxph9 = 0
        nmaxpl9 = 0
        DO J = 1,nsmaxr
          IF(idxr(j) .LT. 0) THEN
            nmaxrl9 = nmaxrl9+1
            hlwr(nmaxrl9) = hhighsr(j)
            tlwr(nmaxrl9) = thighsr(j)
          ELSEIF(idxr(j) .GT. 0) THEN
            nmaxrh9 = nmaxrh9+1
            hhwr(nmaxrh9) = hhighsr(j)
            thwr(nmaxrh9) = thighsr(j)
          ENDIF
        ENDDO

        DO J = 1,nsmaxp
          IF(idxp(j) .LT. 0) THEN
            nmaxpl9 = nmaxpl9+1
            tlwp(nmaxpl9) = thighsp(j)
            hlwp(nmaxpl9) = hhighsp(j)
          ELSEIF(idx0(j) .GT. 0)THEN
            nmaxph9 = nmaxph9+1
            thwp(nmaxph9) = thighsp(j)
            hhwp(nmaxph9) = hhighsp(j)
          ENDIF
        ENDDO
        nmaxrh = nmaxrh9
        nmaxrl = nmaxrl9
        nmaxph = nmaxph9
        nmaxpl = nmaxpl9
        Write(*,*) 'nsmaxr=',nsmaxr,'nsmaxp=',nsmaxp

        if(is .eq. 2) then
          open(50,file=trim(fileshort)//'_timeseries_obsandhind.dat')
          open(51,file=trim(fileshort)//'_extreme_obsandhind.dat')
          do i = 1,nmax
            write(50,'(4f10.4)') tp(i)/24.0,r(i),p(i)
          enddo

          m9 = 0
          do i = 1,nsmaxr
            nc = 0
            timemin = 999.0
            do n2 = 1,nsmaxp
              if(idxr(i)*idxp(n2) .GT. 0) then
                difff = abs(thighsr(i)-thighsp(n2))
                if(difff .lt. timemin) then
                  timemin = difff
                  nc = n2
                endif
              endif
            enddo  

            if((timemin .lt. 15) .and. (nc.gt.0)) then 
              write(51,'(8f12.5)') thighsr(i)/24.0,hhighsr(i),
     1          thighsp(nc)/24.0,hhighsp(nc)
            endif
          enddo
          close(50)
          close(51)
C          close(52)

        elseif (is .eq. 3) then
          open(50,file=trim(fileshort)//'_timeseries_obsandnow.dat')
          open(51,file=trim(fileshort)//'_extreme_obsandnow.dat')
          do i = 1,nmax
            write(50,'(4f10.4)') tp(i)/24.0,r(i),p(i)
          enddo

          m9 = 0
          do i = 1,nsmaxr
            nc = 0
            timemin = 999.0
            do n2 = 1,nsmaxp
              if(idxr(i)*idxp(n2) .GT. 0) then
                difff = abs(thighsr(i)-thighsp(n2))
                if(difff .lt. timemin) then
                  timemin = difff
                  nc = n2
                endif
              endif
            enddo  

            if((timemin .lt. 15) .and. (nc.gt.0)) then 
              write(51,'(8f12.5)') thighsr(i)/24.0,hhighsr(i),
     1          thighsp(nc)/24.0,hhighsp(nc)
            endif
          enddo
          close(50)
          close(51)
C          close(52)

        elseif(is .eq. 4) then
          open(50,file=trim(fileshort)//'_timeseries_obsandfr.dat')
          open(51,file=trim(fileshort)//'_extreme_obsandfr.dat')
          do i = 1,nmax
             write(50,'(4f10.4)') tp(i)/24.,r(i),p(i)
          enddo

          m9 = 0
          do i = 1,nsmaxr
            nc = 0
            timemin = 999.0
            do n2 = 1,nsmaxp
              if(idxr(i)*idxp(n2) .GT. 0) then
                difff = abs(thighsr(i)-thighsp(n2))
                if(difff .lt. timemin) then
                  timemin = difff
                  nc = n2
                endif
              endif
            enddo  

            if((timemin .lt. 15) .and. (nc.gt.0)) then 
              write(51,'(8f12.5)') thighsr(i)/24.0,hhighsr(i),
     1          thighsp(nc)/24.0,hhighsp(nc)
            endif
          enddo
          close(50)
          close(51)
        endif

        IF(KINDAT .eq. 1) THEN
C  Derive current extrema time series of reference (observation)
          write(*,*) 'Calculating slack time of observation'
          CALL SLACK_TIME(tp,r,dirr,thighsr,hhighsr,idxr,NSMAXr,
     1      nmax,dirflood,delt,AFCr,AECr,TFCr,TECr,DFCr,DECr,TSFr,
     2      TEFr,TSEr,TEEr,nafcr,naecr,NTSFr,NTSEr,NTSEEr,imx,nmx)

          write(*,*) 'Calculating slack time of model'
          CALL SLACK_TIME(tp,p,dirp,thighsp,hhighsp,idxp,NSMAXp,
     1      nmax,dirflood,delt,AFCp,AECp,TFCp,TECp,DFCp,DECp,TSFp,
     2      TEFp,TSEp,TEEp,nafcp,naecp,NTSFp,NTSEp,NTSEEp,imx,nmx)

          gap = 25.0
          ttl2 = ' 26 cm/s 25h'
          IX1 = INT(X1*100+0.1)
          write(ttl2,'(a1,I2,a9)') ' ',IX1,' cm/s 24h'
          call extems(AFCR,TFCR,nafcr,afcp,tfcp,nafcp,e,xtmp,nmaxb,1)
          if(nmaxb .GT. 1) Then
            call prtline1('AFC-afc',ttl2,lu,e,xtmp,imx,nmaxb,X1,3,gap)
          end if

          call extems(AECR,TECR,naecr,aecp,tecp,naecp,e,xtmp,nmaxb,1)
          if(nmaxb .GT. 1) then
           call prtline1('AEC-aec',ttl2,lu,e,xtmp,imx,nmaxb,X1,3,gap)
          end if

          write(ttl2,'(a1,f4.2,a7)') ' ',X2,'h   25h'
          call extems(AFCR,TFCR,nafcr,afcp,tfcp,nafcp,e,xtmp,nmaxb,2)
          if(nmaxb .GT. 1) then
            call prtline1('TFC-tfc',ttl2,lu,e,xtmp,imx,nmaxb,X2,3,gap)
          end if

          call extems(AECR,TECR,naecr,aecp,tecp,naecp,e,xtmp,nmaxb,2)
          if(nmaxb .GT. 1) then
            call prtline1('TEC-tec',ttl2,lu,e,xtmp,imx,nmaxb,X2,3,gap)
          end if

          X3 = X2
          write(ttl2,'(a1,f4.2,a7)') ' ',X3/2,'h   25h'
          call extems(AFCR,TSFR,ntsfr,afcp,tsfp,ntsfp,e,xtmp,nmaxb,2)
          if(nmaxb .GT. 1) then
            call prtline1('TSF-tsf',ttl2,lu,e,xtmp,imx,nmaxb,X3,3,gap)
          end if

          call extems(AFCR,TEFR,ntsfr,afcp,tefp,ntsfp,e,xtmp,nmaxb,2)
          if(nmaxb .GT. 1) then
            call prtline1('TEF-tef',ttl2,lu,e,xtmp,imx,nmaxb,X3,3,gap)
          end if

          call extems(AFCR,TSER,ntser,afcp,tsep,ntseep,e,xtmp,nmaxb,2)
          if(nmaxb .GT. 1) then
            call prtline1('TSE-tse',ttl2,lu,e,xtmp,imx,nmaxb,X3,3,gap)
          end if

          call extems(AFCR,TEER,ntser,afcp,teep,ntsep,e,xtmp,nmaxb,2)
          if(nmaxb .GT. 1) then
            call prtline1('TEE-tee',ttl2,lu,e,xtmp,imx,nmaxb,X3,3,gap)
          end if

          ttl2 = ' 22.5 dg 24h'
          write(ttl2,'(a1,F4.1,a7)') ' ',X11,' dg 24h'
          call extems(DFCR,TFCR,nafcr,dfcp,tfcp,nafcp,e,xtmp,nmaxb,1)
          DO I = 1,nmaxb
            if(e(i) .GT. 180.0) e(i) = e(i)-360.0
            if(e(i) .LT. -180.0) e(i) = e(i)+360.0
          ENDDO
          if(nmaxb .GT. 1) then
            call prtline1('DFC-dfc',ttl2,lu1,e,xtmp,imx,nmaxb,X11,
     1        3,gap)
          end if

          call extems(DECR,TECR,naecr,decp,tecp,naecp,e,xtmp,nmaxb,1)
          DO I = 1,nmaxb
            if(e(i) .GT. 180.0) e(i) = e(i)-360.0
            if(e(i) .LT. -180.0) e(i) = e(i)+360.0
          ENDDO
          if(nmaxb .GT. 1) then
            call prtline1('DEC-dec',ttl2,lu1,e,xtmp,imx,nmaxb,X11,
     1        3,gap)
          end if
        ENDIF

        IF(KINDAT .eq. 2) THEN
          IX1 = INT(X1*100+0.1)
          write(ttl2,'(a3,I2,a7)') '   ',IX1,' cm 24h'
          gap = 25.0
     
C  Create extrema series for amplitudes
          call extems(hhwr,thwr,nmaxrh,hhwp,thwp,nmaxph,e,xtmp,nmaxb,1)
          IF((IS .eq. 4) .or. (IS .eq. 5)) then
            if(nmaxb .GT. 1) then
              call prtline1('AHW-ahw  ',ttl2,lu,e,xtmp,imx,nmaxb,X1,
     1          4,gap)
            endif
          else
            if(nmaxb .GT. 1) then
              call prtline1('AHW-ahw  ',ttl2,lu,e,xtmp,imx,nmaxb,X1,
     1          3,gap)
            endif
          endif

          call extems(hlwr,tlwr,nmaxrl,hlwp,tlwp,nmaxpl,e,xtmp,nmaxb,1)
          IF((IS .eq. 4) .or. (IS .eq. 5)) then
            if(nmaxb .GT. 1)then
              call prtline1('ALW-alw  ',ttl2,lu,e,xtmp,imx,nmaxb,X1,
     1          4,gap)
            endif
          else
            if(nmaxb .GT. 1) then
              call prtline1('ALW-alw  ',ttl2,lu,e,xtmp,imx,nmaxb,X1,
     1          3,gap)
            endif
          endif

C  Create extrema series for times 
          ttl2 = '    .5h  25h'
          write(ttl2,'(a1,f4.2,a7)') ' ',X2,' hr 25h'
          call extems(hhwr,thwr,nmaxrh,hhwp,thwp,nmaxph,e,xtmp,nmaxb,2)
          IF((IS .eq. 4) .or. (IS .eq. 5)) then
            if(nmaxb .GT. 1) then
              call prtline1('THW-thw  ',ttl2,lu,e,xtmp,imx,nmaxb,X2,
     1          4,gap)
            endif
          else
            if(nmaxb .GT. 1) then
              call prtline1('THW-thw  ',ttl2,lu,e,xtmp,imx,nmaxb,X2,
     1          3,gap)
            endif
          endif

          call extems(hlwr,tlwr,nmaxrl,hlwp,tlwp,nmaxpl,e,xtmp,nmaxb,2)
          IF((IS .eq. 4) .or. (IS .eq. 5)) then
            if(nmaxb .GT. 1) then
              call prtline1('TLW-tlw  ',ttl2,lu,e,xtmp,imx,nmaxb,X2,
     1          4,gap)
            endif
          else  
            if(nmaxb .GT. 1) then
              call prtline1('TLW-tlw  ',ttl2,lu,e,xtmp,imx,nmaxb,X2,
     1          3,gap)
            endif
          endif
        ENDIF
220   Continue
      CLOSE(lu)
      IF(KINDAT .EQ. 1) THEN
        CLOSE(lu1)
      END IF

      RETURN
      END


C  Construct event time series from model prediction corresponding 
C  to the obs. event time series
      Subroutine collect_events_prd(thighsr,hhighsr,idxr,nsmaxr,
     1   p,tp,nmax,thighsp,hhighsp,idxp,nsmaxp)

      dimension thighsr(nsmaxr),hhighsr(nsmaxr),idxr(nsmaxr)
      dimension thighsp(nsmaxr),hhighsp(nsmaxr),idxp(nsmaxr)
      dimension tp(nmax),p(nmax)
      window = 6.0
      do i = 1,nsmaxr
        t1 = max(tp(1),thighsr(i)-window)
        t2 = min(tp(nmax),thighsr(i)+window)
        aminwl = 999.0
        amaxwl = -999.0
        do ii = 1,nmax
          if((tp(ii) .ge. t1) .and. (tp(ii).lt.t2)) then
            if(p(ii) .lt. aminwl) then
              aminwl = p(ii)
              timemin = tp(ii)
            endif

            if(p(ii) .gt. amaxwl) then
              amaxwl = p(ii)
              timemax = tp(ii)
            endif
          endif
        enddo

        if(idxr(i) .GT. 0) then
          thighsp(i) = timemax
          hhighsp(i) = amaxwl
          idxp(i) = idxr(i)
        elseif(idxr(i) .LT. 0) then
          thighsp(i) = timemin
          hhighsp(i) = aminwl
          idxp(i) = idxr(i)
        endif
      enddo
      nsmaxp = nsmaxr

      RETURN
      END

