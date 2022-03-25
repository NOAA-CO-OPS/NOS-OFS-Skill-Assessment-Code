	CHARACTER*80 BUFFER,staID,shortname,stationname,c1,c2
	CHARACTER*80 OFS,FIN
	OFS='ngofs'
	FIN='/opt/comf/ohms/'//trim(OFS)
	FIN=trim(FIN)//'/info/stationdata.dat'
	print *,'fin=',trim(FIN)
	OPEN(1,file=trim(FIN))
1	READ(1,*,end=99)staID,shortname,stationname
!	print *,staID,shortname,stationname
	READ(1,*)alat,a,c1,alon,a,c2
	FIN=trim(staID)//'  '//trim(staID)//'W_'//trim(OFS)
	FIN=trim(FIN)//'  "'//trim(stationname)//'"'
	WRITE(2,*)trim(FIN)
	WRITE(2,100)alat,-alon,0.0,0.0,0.0
        GOTO 1
100     FORMAT(2F12.5,2x,f3.1,2x,f3.1,2x,f3.1)
99	STOP
	END
