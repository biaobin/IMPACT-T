!201908, this module is from Kilean
!newly added:
!LagrangeInterp, interp, which is from M. Borland's Elegant source code
module MathModule
  implicit none
  private
  real*8, parameter :: pi    = 3.141592653589793d0
  real*8, parameter :: twopi = 6.283185307179586d0
  type :: mathTemplate
    contains
      procedure :: hamsl,iErf,randn
      procedure :: LagrangeInterp,interp,get_vhphi
      procedure :: bessk,bessk0,bessk1,bessi0,bessi1
  end type
  type(mathTemplate), public :: math
  ! === example use ==================================
  ! x = math%randn()             ! normal distribution
  ! x_std = math%std(x,size(x))  ! get s.t.d. of x_std
  ! ==================================================
contains

real*8 function hamsl(self,j,at)
!==================================================================
! same with hammersley except that the sequence now start from "at"
!------------------------------------------------------------------
  implicit none
  class(mathTemplate) :: self
  integer, optional, intent(in) :: j,at
  integer, parameter            :: jmax=10
  real*8                        :: xs(jmax),xsi(jmax)
  integer                       :: i(jmax),nbase(jmax),i1(jmax),i2(jmax)
  data nbase/2,3,5,7,11,13,17,19,23,29/
  data i/jmax*1/
  if((.not. present(j)) .or. (.not. present(at))) then
    call random_number(hamsl) 
    return
  endif
  if(j<1) then
    call random_number(hamsl) 
    return
  endif
  if(at<=0) then
    call random_number(hamsl) 
    return
  endif
  if(at>i(j)) i(j) = at  ! starting index
  xs (j)=0.d0
  xsi(j)=1.0d0
  i2(j)= i(j)
  do
    xsi(j)=xsi(j)/float(nbase(j))
    i1 (j)= i2(j)/nbase(j)
    xs (j)= xs(j)+(i2(j)-nbase(j)*i1(j))*xsi(j)
    i2 (j)= i1(j)
    if(i2(j)<=0) exit
  enddo
  hamsl=xs(j)
  i(j)= i(j)+1
  return
end function 

real*8 function iErf(self,x)
! ==========================================================
! inverted error function
! original author:  Mark Law
! reference: https://github.com/markthelaw/GoStatHelper
! ------------------------------------------------------------      
  implicit none
  class(mathTemplate) :: self
  real*8, intent(in) :: x
  real*8 :: y, w, p
  
  y = x
  if(x>1d0-1d-16) y=1d0-1d-16    ! divergence correction 
  if(x<-1d0+1d-16) y=-1d0+1d-16  ! divergence correction 
  w = -log((1d0-y)*(1d0+y))      ! 1.0 - x * x would induce rounding errors near the boundaries +/-1
  if(w<6.25d0) then
    w = w - 3.125
    p =  -3.6444120640178196996d-21
    p =   -1.685059138182016589d-19 + p * w
    p =   1.2858480715256400167d-18 + p * w
    p =    1.115787767802518096d-17 + p * w
    p =   -1.333171662854620906d-16 + p * w
    p =   2.0972767875968561637d-17 + p * w
    p =   6.6376381343583238325d-15 + p * w
    p =  -4.0545662729752068639d-14 + p * w
    p =  -8.1519341976054721522d-14 + p * w
    p =   2.6335093153082322977d-12 + p * w
    p =  -1.2975133253453532498d-11 + p * w
    p =  -5.4154120542946279317d-11 + p * w
    p =    1.051212273321532285d-09 + p * w
    p =  -4.1126339803469836976d-09 + p * w
    p =  -2.9070369957882005086d-08 + p * w
    p =   4.2347877827932403518d-07 + p * w
    p =  -1.3654692000834678645d-06 + p * w
    p =  -1.3882523362786468719d-05 + p * w
    p =    0.0001867342080340571352 + p * w
    p =  -0.00074070253416626697512 + p * w
    p =   -0.0060336708714301490533 + p * w
    p =      0.24015818242558961693 + p * w
    p =       1.6536545626831027356 + p * w 
  elseif(w<16d0) then
    w = sqrt(w) - 3.25d0
    p =   2.2137376921775787049d-09
    p =   9.0756561938885390979d-08 + p * w
    p =  -2.7517406297064545428d-07 + p * w
    p =   1.8239629214389227755d-08 + p * w
    p =   1.5027403968909827627d-06 + p * w
    p =   -4.013867526981545969d-06 + p * w
    p =   2.9234449089955446044d-06 + p * w
    p =   1.2475304481671778723d-05 + p * w
    p =  -4.7318229009055733981d-05 + p * w
    p =   6.8284851459573175448d-05 + p * w
    p =   2.4031110387097893999d-05 + p * w
    p =   -0.0003550375203628474796 + p * w
    p =   0.00095328937973738049703 + p * w
    p =   -0.0016882755560235047313 + p * w
    p =    0.0024914420961078508066 + p * w
    p =   -0.0037512085075692412107 + p * w
    p =     0.005370914553590063617 + p * w
    p =       1.0052589676941592334 + p * w
    p =       3.0838856104922207635 + p * w      
  else
    w =  sqrt(w) - 5d0
    p =  -2.7109920616438573243d-11
    p =  -2.5556418169965252055d-10 + p * w
    p =   1.5076572693500548083d-09 + p * w
    p =  -3.7894654401267369937d-09 + p * w
    p =   7.6157012080783393804d-09 + p * w
    p =  -1.4960026627149240478d-08 + p * w
    p =   2.9147953450901080826d-08 + p * w
    p =  -6.7711997758452339498d-08 + p * w
    p =   2.2900482228026654717d-07 + p * w
    p =  -9.9298272942317002539d-07 + p * w
    p =   4.5260625972231537039d-06 + p * w
    p =  -1.9681778105531670567d-05 + p * w
    p =   7.5995277030017761139d-05 + p * w
    p =  -0.00021503011930044477347 + p * w
    p =  -0.00013871931833623122026 + p * w
    p =       1.0103004648645343977 + p * w
    p =       4.8499064014085844221 + p * w
  endif
  iErf = p*x
end function      

real*8 function cdfn(self,low,up,lowCut,upCut)
  implicit none
  class(mathTemplate) :: self
  real*8, intent(in) :: low,up
  real*8, optional, intent(in) :: lowCut,upCut
  real*8 :: Erf_lc,Erf_hc,ErfL,ErfH
  if(present(lowCut)) then
    Erf_lc = Erf(lowCut/sqrt(2d0))
  else
    Erf_lc = -1d0
  endif
  if(present(upCut)) then
    Erf_hc = Erf(upCut/sqrt(2d0))
  else
    Erf_hc =  1d0
  endif
  ErfL = Erf(low/sqrt(2d0))
  ErfH = Erf(up/sqrt(2d0))
  cdfn = (ErfH-ErfL)/(Erf_hc-Erf_lc)
end function cdfn

real*8 function randn(self, jHamm, at, low, up)
!==========================================
! hammersley sampling of local gaussian distribution 
! starting at "at" using inverse error function
!------------------------------------------
  implicit none
  class(mathTemplate) :: self
  integer,optional, intent(in) :: jHamm,at
  real*8, optional, intent(in) :: low,up
  real*8                       :: erfL,erfH
  if(present(low)) then
    erfL = Erf(low/sqrt(2d0))
  else
    erfL = -1d0
  endif
  if(present(up)) then
    erfH = Erf(up/sqrt(2d0))
  else
    erfH = 1d0
  endif  
  randn = sqrt(2d0)*self%iErf((erfH-erfL)*self%hamsl(jHamm,at)+erfL)
end function randn

! code form M. Borland's Elegant source C code
real*8 function interp(self,x,Fx,Uj,Nbin1)
  implicit none
  class(mathTemplate) :: self
  real*8, intent(in), dimension(Nbin1) :: x,Fx
  integer,intent(in) :: Nbin1
  real*8,intent(in) :: Uj
  integer :: lo,hi,mid

  lo=1
  hi=Nbin1
  do while ((hi-lo)>1)
      mid = (lo+hi)/2
      if (Uj<Fx(mid)) then
          hi=mid
      else
          lo=mid
      endif
  enddo
  interp = self%LagrangeInterp(x,Fx,lo,Uj,Nbin1)
end function interp

real*8 function LagrangeInterp(self,x,Fx,lo,Uj,Nbin1)
implicit none
  class(mathTemplate) :: self
  real*8, intent(in), dimension(Nbin1) :: x,Fx
  integer,intent(in) :: lo,Nbin1
  real*8,intent(in) :: Uj
  real*8 :: sum,denom,numer
  integer :: i,j
  
  sum=0.0d0
  do i=0,1
    denom=1.0d0
    numer=1.0d0
    do j=0,1
      if (i.ne.j) then
        denom=denom*(Fx(lo+i)-Fx(lo+j))
        numer=numer*(Uj-Fx(lo+j))
        if (numer.eq.0) then
          LagrangeInterp=x(lo+j)
          return
        end if
      end if
    end do
    if (denom.eq.0) then
      LagrangeInterp=0.0d0
      return
    end if
    sum=sum+x(lo+i)*numer/denom
  end do
  LagrangeInterp=sum
end function LagrangeInterp

function get_vhphi(self,tj)
  implicit none
  class(mathTemplate) :: self
  real*8 :: tj
  real*8, allocatable :: t(:),volt(:),phase(:)
  real*8 :: tmp1,tmp2,tmp3,tmp4
  integer :: i,io,linenum,harmnum
  real*8,dimension(3) :: get_vhphi
  
  open(unit=100, file='rfdata_ac.in')
  
  !get the total line number and harmonic number
  read(100,*)linenum,harmnum
  !print*,linenum,harmnum
  
  allocate(t(linenum))
  !allocate(brho(linenum)) !no need to get brho
  allocate(volt(linenum))
  allocate(phase(linenum))
  
  do i=1,linenum
      read(100,*,iostat=io)tmp1,tmp2,tmp3,tmp4
      t(i) = tmp1             !s
      !brho(i) = tmp2
      volt(i) = tmp3*1.0d9    !PyORBIT uses GV, change to V
      phase(i) = tmp4-90      !degree, IMPACT-Z use cos() function
  end do
  close(100)

  if(tj>t(linenum)) then
    print*,"ERROR, turn number too large, out of rf_curve time range."
    stop
  end if
  
  !get the value based on t0 using linear interp,
  !pay attention to the order
  get_vhphi(1) = self%interp(volt,t,tj,linenum)
  get_vhphi(2) = harmnum
  get_vhphi(3) = self%interp(phase,t,tj,linenum)
 
  !print*,"volt=",get_vhphi(1),"h=",get_vhphi(2),"phi=",get_vhphi(3)
  !deallocate the array
  deallocate(t)
  deallocate(volt)
  deallocate(phase)

end function get_vhphi

!-------------------------
! Besselk function
!-------------------------
      FUNCTION BESSK(self,N,X)
      IMPLICIT NONE
      class(mathTemplate) :: self
      INTEGER N,J
      REAL *8 X,BESSK,BESSK0,BESSK1,TOX,BK,BKM,BKP
! ------------------------------------------------------------------------
!     CE SOUS-PROGRAMME CALCULE LA FONCTION BESSEL MODIFIFIEE 3E ESPECE
!     D'ORDRE N ENTIER POUR TOUT X REEL POSITIF > 0.  ON UTILISE ICI LA
!     FORMULE DE RECURRENCE CLASSIQUE EN PARTANT DE BESSK0 ET BESSK1.
!
!     THIS ROUTINE CALCULATES THE MODIFIED BESSEL FUNCTION OF THE THIRD
!     KIND OF INTEGER ORDER, N FOR ANY POSITIVE REAL ARGUMENT, X. THE
!     CLASSICAL RECURSION FORMULA IS USED, STARTING FROM BESSK0 AND BESSK1.
! ------------------------------------------------------------------------
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
! ------------------------------------------------------------------------
      IF (N.EQ.0) THEN
      BESSK = self%BESSK0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSK = self%BESSK1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.D0) THEN
      BESSK = 1.D30
      RETURN
      ENDIF
      TOX = 2.D0/X
      BK  = self%BESSK1(X)
      BKM = self%BESSK0(X)
      DO 11 J=1,N-1
      BKP = BKM+DFLOAT(J)*TOX*BK
      BKM = BK
      BK  = BKP
   11 CONTINUE
      BESSK = BK
      RETURN
      END
! ----------------------------------------------------------------------
      FUNCTION BESSK0(self,X)
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DU 3EME ESPECE D'ORDRE 0
!     POUR TOUT X REEL NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF
!     ORDER ZERO FOR ANY POSITIVE REAL ARGUMENT, X.
! ----------------------------------------------------------------------
      IMPLICIT NONE
      class(mathTemplate) :: self
      REAL*8 X,BESSK0,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,    &
      BESSI0
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0, &
      0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1, &
      -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF(X.EQ.0.D0) THEN
      BESSK0=1.D30
      RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
      Y=X*X/4.D0
      AX=-LOG(X/2.D0)*self%BESSI0(X)
      BESSK0=AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      Y=(2.D0/X)
      AX=EXP(-X)/DSQRT(X)
      BESSK0=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
! ----------------------------------------------------------------------
      FUNCTION BESSK1(self,X)
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DE 3EME ESPECE D'ORDRE 1
!     POUR TOUT X REEL POSITF NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF
!     ORDER ONE FOR ANY POSITIVE REAL ARGUMENT, X.
! ----------------------------------------------------------------------
      IMPLICIT NONE
      class(mathTemplate) :: self
      REAL*8 X,BESSK1,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,BESSI1
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,0.15443144D0,-0.67278579D0,  &
      -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1, &
      0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
      IF(X.EQ.0.D0) THEN
      BESSK1=1.D32
      RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
      Y=X*X/4.D0
      AX=LOG(X/2.D0)*self%BESSI1(X)
      BESSK1=AX+(1.D0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      Y=(2.D0/X)
      AX=EXP(-X)/DSQRT(X)
      BESSK1=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
!
!     Bessel Function of the 1st kind of order zero.
!
      FUNCTION BESSI0(self,X)
      IMPLICIT NONE
      class(mathTemplate) :: self
      REAL *8 X,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067429D0,  &
      0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,  &
      0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
      0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/DSQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
      ENDIF
      RETURN
      END
!
!     Bessel Function of the 1st kind of order one.
!
      FUNCTION BESSI1(self,X)
      IMPLICIT NONE
      class(mathTemplate) :: self
      REAL *8 X,BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,  &
      0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1, &
      -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,        &
      -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/DSQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
      ENDIF
      RETURN
      END

end module MathModule
