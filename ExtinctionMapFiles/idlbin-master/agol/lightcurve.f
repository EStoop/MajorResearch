      program lightcurve
      implicit none
      integer i,nb
      parameter(nb=2000)
      real*8 b0(nb),muom(nb),muop(nb),mu(nb)
      real*8 rs,rl
      real*8 c1,c2,c3,c4,mulimb(nb,5),mulimb0(nb)
      read(5,*) rs,rl,c1,c2,c3,c4
      do i=1,nb
        b0(i)=4.d0*dble(i-1)/dble(nb-1)
      enddo
      if(c1.eq.0.d0.and.c2.eq.0.d0.and.c3.eq.0.d0.and.c4.eq.0.d0) then
c        write(6,*) 'Calling extsource'
        call extsource(rs,rl,b0,muom,muop,mu,nb)
        do i=1,nb
          mulimb0(i)=muom(i)+muop(i)
        enddo
      else
c        write(6,*) 'Calling microccultnl'
        call microccultnl(rs,rl,c1,c2,c3,c4,b0,mulimb0,mulimb,nb)
      endif
      do i=1,nb
        write(6,*) b0(i),mulimb0(i)
      enddo
      stop
      end


      subroutine extsource(rs,rl,b0,muom,muop,mu,nb)
      implicit none
      integer i,nb
      real*8 pi,b0(nb),muom(nb),muop(nb),mu(nb)
      real*8 a1,a2,be,bl,bvec,k,n,rs,rl,phi0,phi1,phi2,
     &       psi0,psi1,psi2,u0,u1,u2,u3,v1,v2,gofphi
      pi=acos(-1.d0)
C This routine computes the lightcurve for microlensing & occultation of a uniform
C extended source by a point lens (Agol 2002).
C
CInput:
C  rs      source size in units of Einstein radius projected onto source plane
C  rl      radius of the lens in units of the Einstein radius
C  bvec    position of source in source-radius units
C  bmin    minimum impact parameter (if bvec is unspecified)
C  bmax    maximum impact parameter (if bvec is unspecified)
C
COutput:
C  b0      impact parameter grid (units of source radius)
C  muom    magnification of negative (inner) image
C  muop    magnification of positive (outer) image
C  mu      1/2 of the magnification of an unocculted extended uniform source
C 
C If the call doesn't specify a list of impact parameters, then 
C generate a list:
      do i=1,nb
        bvec=b0(i)
C position of source in Einstein-radius units (in the paper
C this is the parameter \zeta_0):
        be=bvec*rs
C Point-source lensing formula:      
        if(rs.eq.0.d0) then
          if(rl.lt.1.d0) then
            muop(i)=0.5d0*( 1.d0+(be**2+2.d0)/be/sqrt(be**2+4.d0))
            if(be*rl.lt.1.d0-rl**2) then
              muom(i)=0.5d0*(-1.d0+(be**2+2.d0)/be/sqrt(be**2+4.d0))
            else
              muom(i)=0.d0
            endif
          else
            if(be*rl.gt.rl**2-1.d0) then
              muop(i)=0.5d0*( 1.d0+(be**2+2.d0)/be/sqrt(be**2+4.d0))
            else
              muop(i)=0.d0
            endif
            muom(i)=0.d0
          endif
          goto 10
        endif
      
C Now, compute expression from Witt & Mao (1994):
        n=4.d0*be*rs/(be+rs)**2
        k=sqrt(4.d0*n/(4.d0+(be-rs)**2))
        if(abs(be-rs).gt.1.d-7) then
          mu(i)=gofphi(pi/2.d0,be,rs)/4.d0/pi/rs**2
        else
C Use special formula when the impact parameter equals the source star radius:
          mu(i)=(rs+(1.d0+rs**2)*atan(rs))/pi/rs**2
        endif

C Now, compute the occulted magnification:
C First, set the magnifications equal to the unocculted magnification:
        muom(i)=mu(i)-0.5d0
        muop(i)=mu(i)+0.5d0
C There are 3 cases - rl=0, 0<rl<1, rl>1:
C (1) rl = 0
C There is no occultation - same as Witt & Mao (1995) expression:
        if(rl.eq.0.d0) goto 10
      
        if(rl.eq.1.d0) then
C Inner image is occulted, outer is not:
          muom(i)=0.d0
          goto 10
        endif
C (2) 0 < rl < 1:
C     In this case, the inner (negative parity) image can be occulted, 
C     while the outer (positive parity) image is never occulted:
        if(rl.lt.1.d0) then
C Project lens radius onto source plane:
          bl=1.d0/rl-rl
C (a) Inner image is completely occulted by lens star:
          if(be.gt.bl+rs) then begin
            muom(i)=0.d0 
            goto 5
          endif
C (b) When the edge of the star crosses the lens singularity, we need a 
C special expression:
          if(abs(be-rs).lt.1.d-7) then begin
C Partial occultation (equation 19):
            if(rs.gt.0.5d0*bl) then begin
              v1=sqrt(4.d0+bl**2)
              v2=sqrt(4.d0*rs**2-bl**2)
              a1=(0.25d0*v2*(bl-v1)+(1.d0+rs**2)*
     &           (acos(0.5d0*bl/rs)-atan(v2/v1)))/pi/rs**2
              a2=rl**2/pi/rs**2*acos(0.5d0*bl/rs)
              muom(i)=muom(i)+a1-a2
            endif
            goto 5
          endif
C (c) Inner image is unocculted - use Witt & Mao expression:
          if(be.le.bl-rs) goto 5
C (d) When the source star covers the unocculted window, then the magnification
C is a constant:
          if(be.lt.rs-bl) then
            muom(i)=(1.d0-rl**2)/rs**2 
            goto 5
          endif
C (e) The inner image is partially occulted by the lens star:
C (see 4/18/02 notes) equation (16)
          u0=bl**2
          u1=(be-rs)**2
          u2=(be+rs)**2
          u3=be**2-rs**2
C Now, use notation from equation (18) in paper:
          phi0=acos(sqrt(u1*(u2-u0)/u0/(u2-u1)))
          phi1=acos((u1+u2-2.d0*u0)/(u2-u1))
          phi2=acos((u3+u0)/2.d0/be/bl)
          if(be.lt.rs) then
            muom(i)=1.d0/rs**2
          else
            muom(i)=0.d0
          endif
          muom(i)=muom(i)-(-4.d0*u3/abs(u3)*phi0+2.d0*(1.d0+rs**2)
     &       *phi1+4.d0*rl**2*phi2+sqrt((u2-u0)*(u0-u1))
     &       *(sqrt(1.d0+4.d0/u0)-1.d0)-gofphi(phi0,be,rs))
     &       /4.d0/pi/rs**2
 5        continue
C The outer (positive image) is unocculted.
          muop(i)=muop(i)
          goto 10
        endif 
      
C Now for rl gt 1 case:
C Compute the size of the lens star projected onto the source plane:
        bl=rl-1.d0/rl
C Image is unocculted - use usual formula:
        if(be.gt.bl+rs) goto 15
C Image is completely occulted:
        if(rs.le.bl.and.be.le.bl-rs) then
          muop(i)=0.d0
          goto 15
        endif
C Image has a hole in the center:
        if(rs.gt.bl.and.abs(be).lt.rs-bl) then begin
          muop(i)=mu(i)+0.5d0-(rl**2-1.d0)/rs**2 
          goto 15
        endif
C Image is partially occulted:  use equations (22) & (23):
        u1=(be-rs)**2
        u2=(be+rs)**2
        u3=be**2-rs**2
        u0=bl**2
        psi0=acos(sqrt((u0-u1)*(4.d0+u2)/(4.d0+u0)/(u2-u1)))
        phi0=acos(sqrt(u1*(u2-u0)/u0/(u2-u1)))
        phi1=acos((u1+u2-2.d0*u0)/(u2-u1))
        phi2=acos((u3+u0)/2.d0/be/bl)
        psi1=pi/2.d0-phi0
        psi2=pi-phi1+2.d0*acos(sqrt(u0*(4.d0+u0)
     &       /(u0*(4.d0+u1+u2)-u1*u2)))
        if(u3.ne.0.d0) then
          muop(i)=(-4.d0*u3/abs(u3)*psi1+2.d0*(1.d0+rs**2)*psi2
     &            -4.d0*rl**2*phi2
     &        +sqrt((u2-u0)*(u0-u1))*(sqrt(u0/(4.d0+u0))+1.d0)
     &        +gofphi(psi0,be,rs))/4.d0/pi/rs**2
        else
          muop(i)=(-4.d0*psi1+2.d0*(1.d0+rs**2)*psi2-4.d0*rl**2*phi2
     &        +sqrt((u2-u0)*(u0-u1))*(sqrt(u0/(4.d0+u0))+1.d0)
     &        +gofphi(psi0,be,rs))/4.d0/pi/rs**2
        endif
 15    continue
C Negative (inner) image is always occulted:
        muom(i)=0.d0
 10     continue
      enddo
      return
      end
      
      function gofphi(phi,z0,rs)
      implicit none
      real*8  gofphi,phi,z0,rs,u1,u2,u3,en,k,ee,ff,pp,gf,
     &        cp2,sp,yy,rf,rd,rj
      u1=(z0-rs)**2
      u2=(z0+rs)**2
      u3=z0**2-rs**2
      en=1.d0-u1/u2
      k=sqrt(4.d0*(u2-u1)/u2/(4.d0+u1))
      sp=sin(phi)
      cp2=1.d0-sp*sp
      yy=1.d0-k*k*sp*sp
      ff=sp*rf(cp2,yy,1.d0)
      ee=ff-k*k*sp*sp*sp/3.d0*rd(cp2,yy,1.d0)
      pp=ff+en*sp*sp*sp/3.d0*rj(cp2,yy,1.d0,1.d0-en*sp*sp)
      gf=(u2*(4.d0+u1)*ee-(u1*u2+8.d0*u3)*ff+4.d0*u1*(1.d0+rs**2)*pp)
     &   /sqrt(u2*(4.d0+u1))
      gofphi=gf
      end

      FUNCTION rf(x,y,z)
      REAL*8 rf,x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08d0,TINY=1.5d-38,BIG=3.d37,THIRD=1.d0/3.d0,
     *C1=1.d0/24.d0,C2=.1d0,C3=3.d0/44.d0,C4=1.d0/14.d0)
      REAL*8 alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if(min(x,y,z).lt.0.d0.or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        ave=THIRD*(xt+yt+zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0NL&WR2.

      FUNCTION rd(x,y,z)
      REAL*8 rd,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
      PARAMETER (ERRTOL=.05d0,TINY=1.d-25,BIG=4.5d21,C1=3.d0/14.d0,
     *C2=1.d0/6.d0,C3=9.d0/22.d0,C4=3.d0/26.d0,C5=.25d0*C3,C6=1.5d0*C4)
      REAL*8 alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,
     *sqrtz,sum,xt,yt,zt
      if(min(x,y).lt.0.d0.or.min(x+y,z).lt.TINY.or.max(x,y,
     *z).gt.BIG)pause 'invalid arguments in rd'
      xt=x
      yt=y
      zt=z
      sum=0.d0
      fac=1.d0
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=.25d0*fac
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        ave=.2d0*(xt+yt+3.d0*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6.d0*eb
      ee=ed+ec+ec
      rd=3.d0*sum+fac*(1.d0+ed*(-C1+C5*ed-C6*delz*ee)+delz*
     *(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      FUNCTION rj(x,y,z,p)
      REAL*8 rj,p,x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6,C7,C8
      PARAMETER (ERRTOL=.05d0,TINY=2.5d-13,BIG=9.d11,C1=3.d0/14.d0,
     *C2=1.d0/3.d0,C3=3.d0/22.d0,C4=3.d0/26.d0,C5=.75d0*C3,
     *C6=1.5d0*C4,C7=.5d0*C2,C8=C3+C3)
CU    USES rc,rf
      REAL*8 a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,
     *ed,ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt,rc,rf
      if(min(x,y,z).lt.0.d0.or.min(x+y,x+z,y+z,abs(p)).lt.TINY.or.
     *max(x,y,z,abs(p)).gt.BIG)pause 'invalid arguments in rj'
      sum=0.d0
      fac=1.d0
      if(p.gt.0.d0)then
        xt=x
        yt=y
        zt=z
        pt=p
      else
        xt=min(x,y,z)
        zt=max(x,y,z)
        yt=x+y+z-xt-zt
        a=1.d0/(yt-p)
        b=a*(zt-yt)*(yt-xt)
        pt=yt+b
        rho=xt*zt/yt
        tau=p*pt/yt
        rcx=rc(rho,tau)
      endif
1     continue
        sqrtx=sqrt(xt)
        sqrty=sqrt(yt)
        sqrtz=sqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
        beta=pt*(pt+alamb)**2
        sum=sum+fac*rc(alpha,beta)
        fac=.25d0*fac
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        zt=.25d0*(zt+alamb)
        pt=.25d0*(pt+alamb)
        ave=.2d0*(xt+yt+zt+pt+pt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
        delp=(ave-pt)/ave
      if(max(abs(delx),abs(dely),abs(delz),abs(delp)).gt.ERRTOL)goto 1
      ea=delx*(dely+delz)+dely*delz
      eb=delx*dely*delz
      ec=delp**2
      ed=ea-3.d0*ec
      ee=eb+2.d0*delp*(ea-ec)
      rj=3.d0*sum+fac*(1.d0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*
     *(-C8+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave))
      if (p.le.0.d0) rj=a*(b*rj+3.d0*(rcx-rf(xt,yt,zt)))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software

      FUNCTION rc(x,y)
      REAL*8 rc,x,y,ERRTOL,TINY,SQRTNY,BIG,TNBG,COMP1,COMP2,THIRD,
     *C1,C2,C3,C4
      PARAMETER (ERRTOL=.04d0,TINY=1.69d-38,SQRTNY=1.3d-19,BIG=3.d37,
     *TNBG=TINY*BIG,COMP1=2.236d0/SQRTNY,COMP2=TNBG*TNBG/25.d0,
     *THIRD=1.d0/3.d0,C1=.3d0,C2=1.d0/7.d0,C3=.375d0,C4=9.d0/22.d0)
      REAL*8 alamb,ave,s,w,xt,yt
      if(x.lt.0.d0.or.y.eq.0.d0.or.(x+abs(y)).lt.TINY.or.(x+
     *abs(y)).gt.BIG.or.(y.lt.-COMP1.and.x.gt.0.d0.and.x.lt.COMP2))pause 
     *'invalid arguments in rc'
      if(y.gt.0.d0)then
        xt=x
        yt=y
        w=1.d0
      else
        xt=x-y
        yt=-y
        w=sqrt(x)/sqrt(xt)
      endif
1     continue
        alamb=2.d0*sqrt(xt)*sqrt(yt)+yt
        xt=.25d0*(xt+alamb)
        yt=.25d0*(yt+alamb)
        ave=THIRD*(xt+yt+yt)
        s=(yt-ave)/ave
      if(abs(s).gt.ERRTOL)goto 1
      rc=w*(1.d0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software
      subroutine microccultnl(rse,rle,c1,c2,c3,c4,b0,mulimb0,mulimbf,nb)
c Please cite Agol (2002) if making use of this routine.
      implicit none
      integer i,j,nb,nr,i1,i2,nmax
      parameter (nmax=513)
      real*8 mulimbf(nb,5),pi,c1,c2,c3,c4,rl,bt0(nb),b0(nb),mulimb0(nb),
     &       mulimb(nb),mulimbp(nb),dt,t(nmax),th(nmax),r(nmax),sig,
     &       mulimb1(nb),mulimbhalf(nb),mulimb3half(nb),mulimb2(nb),
     &       sig1,sig2,omega,dmumax,fac,mu(nb),f1,f2,rs,rse,rle,
     &       muop(nb),muom(nb)
      pi=acos(-1.d0)
C  This routine uses the results for a uniform source to
C  compute the lightcurve for a limb-darkened source
C  (5-1-02 notes)
C Input:
C   rle       radius of the lens   in units of the Einstein radius
C   rse       radius of the source in units of the Einstein radius
C   c1-c4     limb-darkening coefficients
C   b0        impact parameter normalized to source radius
C Output:
C  mulimb0 limb-darkened magnification
C  mulimbf lightcurves for each component
C  
C  First, make grid in radius:
C  Call magnification of uniform source:
      call extsource(rse,rle,b0,muop,muom,mu,nb-1)
      do i=1,nb
        mulimb0(i)=muom(i)+muop(i)
      enddo
c      call occultuniform(b0,rl,mulimb0,nb)
      i1=nb
      i2=1
      fac=0.d0
      do i=1,nb
        bt0(i)=b0(i)
        mulimbf(i,1)=1.d0
        mulimbf(i,2)=0.8d0
        mulimbf(i,3)=2.d0/3.d0
        mulimbf(i,4)=4.d0/7.d0
        mulimbf(i,5)=0.5d0
        mulimb(i)=mulimb0(i)
        if(mulimb0(i).ne.1.d0) then
          i1=min(i1,i)
          i2=max(i2,i)
        endif
        fac=max(fac,abs(mulimb0(i)-1.d0))
      enddo
C print,rl
      omega=4.*((1.d0-c1-c2-c3-c4)/4.+c1/5.+c2/6.+c3/7.+c4/8.)
      nr=2
      dmumax=1.d0
c      write(6,*) 'i1,i2 ',i1,i2
      do while (dmumax.gt.fac*1.d-3)
        do i=i1,i2
          mulimbp(i)=mulimb(i)
        enddo
        nr=nr*2
c        write(6,*) 'nr ',nr
        dt=0.5d0*pi/dble(nr)
        do j=1,nr+1
          t(j) =dt*dble(j-1)
          th(j)=t(j)+0.5d0*dt
          r(j)=sin(t(j))
        enddo
        sig=sqrt(cos(th(nr)))
        do i=i1,i2
          mulimbhalf(i) =sig**3*mulimb0(i)/(1.d0-r(nr))
          mulimb1(i)    =sig**4*mulimb0(i)/(1.d0-r(nr))
          mulimb3half(i)=sig**5*mulimb0(i)/(1.d0-r(nr))
          mulimb2(i)    =sig**6*mulimb0(i)/(1.d0-r(nr))
        enddo
        do j=2,nr
          do i=1,nb
            b0(i)=bt0(i)/r(j)
          enddo
C  Calculate uniform magnification at intermediate radii:
          call extsource(rse*r(j),rle,b0,muop,muom,mu,nb-1)
          do i=1,nb
            mu(i)=muom(i)+muop(i)
          enddo
c          call occultuniform(b0,rl/r(j),mu,nb)
C  Equation (29):
          sig1=sqrt(cos(th(j-1)))
          sig2=sqrt(cos(th(j)))
          dmumax=0.d0
          do i=i1,i2
            f1=r(j)*r(j)*mu(i)/(r(j)-r(j-1))
            f2=r(j)*r(j)*mu(i)/(r(j+1)-r(j))
            mulimbhalf(i) =mulimbhalf(i) +f1*sig1**3-f2*sig2**3
            mulimb1(i)    =mulimb1(i)    +f1*sig1**4-f2*sig2**4
            mulimb3half(i)=mulimb3half(i)+f1*sig1**5-f2*sig2**5
            mulimb2(i)    =mulimb2(i)    +f1*sig1**6-f2*sig2**6
            mulimb(i)=((1.d0-c1-c2-c3-c4)*mulimb0(i)+c1*mulimbhalf(i)*dt
     &        +c2*mulimb1(i)*dt+c3*mulimb3half(i)*dt+c4*mulimb2(i)*dt)
     &        /omega
            if(mulimb(i)+mulimbp(i).ne.0.d0) then 
              dmumax=max(dmumax,abs(mulimb(i)-mulimbp(i))/(mulimb(i)+
     &               mulimbp(i)))
            endif
          enddo
        enddo
      enddo
      do i=i1,i2
        mulimbf(i,1)=mulimb0(i)
        mulimbf(i,2)=mulimbhalf(i)*dt
        mulimbf(i,3)=mulimb1(i)*dt
        mulimbf(i,4)=mulimb3half(i)*dt
        mulimbf(i,5)=mulimb2(i)*dt
        mulimb0(i)=mulimb(i)
      enddo
      do i=1,nb
        b0(i)=bt0(i)
      enddo
      return
      end
