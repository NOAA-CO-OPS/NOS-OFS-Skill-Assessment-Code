
      parameter (kmx=735)              ! Max number of days                  
      parameter (imx=kmx*240)          ! Length of time series of tides and observations
      parameter (nmx=imx/10)           ! Maximum number of high/lows in a series (up to 5/day)   
      parameter (jmx=4)                ! Number of forecasts per day (includes  hour 00 06 12 18)       
      parameter (MAXDURATION=190*10+1) ! Maximum number of forecasts per cycle 
                                       ! (e.g. forecast duration (in hours) X data number per hour (1/DELT)

      logical FINCLUDE0
      real*8 delt,delt_o,delt_m,delt_t
      real*8 t,tc,tstart,tfinish,tf
      common/a0/tstart,tfinish,jmax,nday,delt,delt_o,delt_m,delt_t,
     1  ipday,tcut,ISWITCH(6),criteria1,criteria2,method,Igapfill,
     2  dirflood,imaxa,IMAXB,IMAXF,NTYPE,FACTOR,NFDURATION,ISURGE,
     3  FINCLUDE0
      common/a1/t(imx),a(imx),er(imx),o(imx),of(imx),s(imx),
     1  hi(imx),sgo(imx),sgm(imx)
      common/a2/c(imx),tc(imx),f(MAXDURATION,jmx,kmx),
     1  tf(MAXDURATION,jmx,kmx),pf(MAXDURATION,jmx,kmx)
      common/a3/wcof,skillv,corr_c
      common/a4/hhwr(nmx),thwr(nmx),hlwr(nmx),tlwr(nmx),nmaxrh,nmaxrl,
     1          hhwp(nmx),thwp(nmx),hlwp(nmx),tlwp(nmx),nmaxph,nmaxpl
      common/a5/ipr1,delhr,delamp,delpct,iopta,iyear,x1,x2,x11
      common/a6/dira(imx),diro(imx),dirof(imx),dirs(imx),dirhi(imx),
     1  dirc(imx),dirf(MAXDURATION,jmx,kmx),dirpf(MAXDURATION,jmx,kmx)

