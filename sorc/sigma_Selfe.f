       
!       SUBROUTINE sigma2Z_SELFE(sigma,H,ele,KB,zsigma)
       parameter(KB=10)
       dimension sigma(KB),zsigma(KB),cs(Kb),dcs(kb)
       data sigma/-1,-0.888888,-0.777777,-0.666666,-0.555555,-0.444444
     1   ,-0.333333, -0.222222, -0.111111, 0/
       H=12.1850004
       ele=0.8452
       NVRT=10
       KZ=1
       h_s=4000.0
       thetas=4.5
       thetab=0.95
       tc=10.0
       h_c=10.0
       CALL sigma2Z_SELFE(sigma,H,ele,KB,zsigma
     1            ,h_s,h_c,thetas,thetab)
       do k=1,kb
          write(6,*)k,zsigma(K),sigma(k)
       enddo !k
       stop
       END
       SUBROUTINE sigma2Z_SELFE(sigma,H,ele,KB,zsigma
     1            ,h_s,h_c,thetas,thetab)
       dimension sigma(KB),zsigma(KB),cs(kb),dcs(kb)
       NVRT=kb
       KZ=1
!       h_s=4000.0
!       thetas=4.5
!       thetab=0.95
!       hc=10.0
!       ztot(kz)=-h_s
       nsig=nvrt-kz+1 !# of S levels (including "bottom" & f.s.)
       s_con1=sinh(thetas)
!       sigma(1)=-1 !bottom
!       sigma(nsig)=0 !surface
!     Compute C(s) and C'(s)
       do k=1,nsig
         cs(k)=(1-thetab)*sinh(thetas*sigma(k))/sinh(thetas)+
     1     thetab*(tanh(thetas*(sigma(k)+0.5))-tanh(thetas*0.5))
     2     /2/tanh(thetas*0.5)
         dcs(k)=(1-thetab)*thetas*cosh(thetas*sigma(k))/
     1      sinh(thetas)+ thetab*thetas/2/tanh(thetas*0.5)
     2      /cosh(thetas*(sigma(k)+0.5))**2
       enddo !k=1,nvrt
       do k=kz,nvrt
         kin=k-kz+1
         hmod2=amin1(H,h_s)
!         z0=h_c*sigma(k)+(H-h_c)*cs(kin)
!         zsigma(k)=z0+(1.+z0/H)*ele
         if(hmod2<=h_c) then
           zsigma(K)=sigma(kin)*(hmod2+ele)+ele
         else
           zsigma(K)=ele*(1+sigma(kin))+h_c*sigma(kin)
     1       +(hmod2-h_c)*cs(kin)
         endif
         IF (zsigma(k) .LT. 0.0)zsigma(k)=-zsigma(k)
 !        print *,cs(kin),sigma(k),zsigma(k),hmod2,z0
       enddo !k
       RETURN
       END
