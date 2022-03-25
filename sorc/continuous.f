      parameter(imx=2000)
      dimension time(imx),m(imx),ifrst(imx),ilast(imx)
      do I=1,10
        time(i)=1.+(i-1)*0.1
      enddo
      do I=11,15
        time(i)=1.+(i-1)*0.3
      enddo
      do I=16,35
        time(i)=1.+(i-1)*0.1
      enddo
      do I=36,45
        time(i)=1.+(i-1)*0.4
      enddo
      do I=46,100
        time(i)=1.+(i-1)*0.1
      enddo
!        time(99)=1.0+98*0.2
!        time(100)=1+99*0.1
      do I=1,100
        time(i)=1.+(i-1)*0.1
      enddo
      gap=0.15
      call continuous(time,100,gap,mdo,Istart,IEND)
      end
      subroutine continuous(z,imx,gap,mdo,Istart,IEND)
c        this program to search maximum continuous time duration in a time series
      dimension z(imx),m(imx),ifrst(imx),ilast(imx)
      m(1)=1
      ipr=1
      do i=2,imx
        gp=z(i)-z(i-1)
        m(i)=1
        if(gp.ge.gap)m(i)=0     
      enddo
      do i=1,imx
        write(*,*)'time,m=',I,z(i),m(i)
      enddo
c        look for start of a series of '1's.  ifrst is i of first '1'
      is=0
      do i=1,imx-1
        if(i.eq.1.and.m(i).eq.1)then
          is=is+1
          ifrst(is)=i
        endif
        if(i.ge.1.and.m(i).eq.0.and.m(i+1).eq.1)then
          is=is+1
          ifrst(is)=i+1
        endif
      enddo
      if(ipr.eq.1)then
        write(6,*)' is=',is
        if(is.gt.0)then
          do i=1,is
            write(6,*)' i=',i,' ifrst=',ifrst(i)
          enddo
        endif
      endif

c        look for end of a series of '1's.  ilast is i of last '1'
      ie=0
      do i=1,imx-1
        if(m(i).eq.1.and.m(i+1).eq.0)then
          ie=ie+1
          ilast(ie)=i
        endif
      enddo
      if(m(imx).eq.1)then
        ie=ie+1
        ilast(ie)=i
      endif

      if(ipr.eq.1)then
        write(6,*)' ie=',ie
        if(ie.gt.0)then
          do i=1,ie
            write(6,*)' i=',i,' ilast=',ilast(i)
          enddo
        endif
      endif

c        find longest duration
      if(is.ne.ie)then
        write(6,"(/,3x,'**function tmdo. starts not equal to ends**')")
        write(6,"(5x,'is=',i4,' ie=',i4)")is,ie
        if(is.gt.0)then
          do i=1,is
            write(6,*)' i=',i,' ifrst=',ifrst(i)
          enddo
        endif
        if(ie.gt.0)then
          do i=1,ie
            write(6,*)' i=',i,' ilast=',ilast(i)
          enddo
        endif
        do i=1,imx
          write(6,"(' i=',i4,' m=',i2)")i,m(i)
        enddo
        stop
      else if(is.gt.0)then
        mdo=0
        do i=1,is
          NDIF=ilast(i)-ifrst(i)
          IF (NDIF .GT. mdo)then
            mdo=NDIF
            IMAX=I
          ENDIF  
        enddo
      endif
      Istart=ifrst(imax)
      IEND=ilast(imax)
      write(*,*)'mdo=',mdo,IMAX,Istart,IEND 
      return
      end
