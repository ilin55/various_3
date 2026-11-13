*     .         Programme C-M_SPH (version 5.1q) *** 28.01.17  ***
*     .          C-M_SPH (ver.4.2qb + hom6_6q, Peterhof)
*     .         Programme C-M_SPH (version 4.2q...b benchmarks) *** 21.08.04  ***
*     .         Programme C-M_SPH (version 4.2q)   *** 21.08.04  ***
*     .         (Peterhof, functions for prolate + DIFF)
*     .         Programme C-M_SPH (version 4.2q)   *** 23.07.04  ***
*     .         ('Krupeli-Tolmachevo')
*     ,           from cm4_1qwo + complex functons cdrsf212, cdrsf20
*     .         Programme C-M_SPH (version 4.1q)   *** 11.07.04  ***
*     .         ('Krupeli-Tolmachevo', +  Jaffe's functions + ...
*     .          from hom5_2qw, 23.02.04, Peterhof)
*     .         Programme C-M_SPH (version 3.3q)   *** 04.12.01  ***
*     .         (Bielefeld, new subroutines BESSJJ0, CESSEL0, CDRB12aa
*     .          cinvm01)
*     .........................................................................
*     .  Light Scattering by Core-mantle Spheroidal Particles
*     .  by N.V.Voshchinnikov, V.V.Somsikov and V.G.Farafonov
*     .  Copyright (c) 1994/97/2004/2017 Astronomical Institute,
*     .                        St. Petersburg University
*     ........................................................................
*     .........................................................................
*     .  Calculations of extinction, scattering and absorption efficiency
*     .  factors for prolate or oblate core-mantle spheroids.
*     .........................................................................
*     .  Description of method, equations and some results are in
*     .  V.G.Farafonov, N.V.Voshchinnikov and V.V.Somsikov,
*     .  "Light scattering by a core-mantle spheroidal particle",
*     .  Applied Optics, 1996, V.35, N 27, P.5412-5426.
*     .........................................................................
*     .........................................................................
*     .  Parameters:
*     .      nterms: number of terms in Eqs. for efficiency factors,
*     .     nalphas: number of angles for incident radiation,
*     .........................................................................
*     .    Input data:           file:   cmnn.dat
*     .             k: parameter describing type of particle,
*     .                           = 0 for prolate spheroids,
*     .                           = 1 for oblate spheroids,
*     .           a/b: outer (mantle) aspect ratio for spheroid,
*     . rii = n + k*i: complex refraction index of core,
*     .  ri = n + k*i: complex refraction index of mantle,
*     .          eps: relative error for calculations of spheroidal functions,
*     .         eps1: relative error for calculations of efficiency factors,
*     .          wvl: wavelength of radiation (in mkm),
*     .          icx: parameter describing outer size (size of mantle)
*     .                     of spheroids to be input,
*     .       = 1 - Rv (outer radius of equivolume sphere, in mkm),
*     .       = 2 - Xv (outer size parameter of equivolume sphere),
*     .       = 3 - c  (outer size parameter proportional to focus distance),
*     .       = 4 - 2*pi*a/wvl (=c*ksi1 for prolate spheroids),
*     .                        (=c*ksi1*(a/b) for oblate spheroids),
*     .       = 5 - Xc (outer size parameter for cylinders),
*     .          NRv: number of sizes of spheroids to be calculated,
*     .         aRv0: initial size,
*     .          dRv: step over sizes,
*     .        nalph: number of angles for incident radiation 'alpha',
*     .         ial0: initial angle 'alpha',
*     .         idal: step over angles 'alpha',
*     .      iswitch: parameter describing core of spheroids to be input,
*     .       = 0 - a/b of core (aspect ratio for core),
*     .       = 1 - Vol_rat (ratio of core volume,
*     .                      to total volume of particle),
*     .       = 2 - Xv of core (size parameter of equivolume sphere for core),
*     .       = 3 - Rv of core (radius of equivolume sphere for core, in mkm).
*     .........................................................................
*     . IMPORTANT: in all cases must be: (a/b)_core > (a/b)_mantle
*     .........................................................................
*     .  Output data:         files cm.res, cm.nrm
*     .........................................................................
*     .   file:   cm.res
*     .   Qext : extinction efficiency factors,
*     .   Pext : extinction efficiency factors (imaginary part),
*     .   Qsca : scattering efficiency factors,
*     .   Qabs : absorption efficiency factors.
*     .........................................................................
*     .   file:   cm.nrm
*     .   Cext/(Pi*Rv**2) : normalized extinction cross-sections,
*     .   Csca/(Pi*Rv**2) : normalized scattering cross-sections,
*     .   Cabs/(Pi*Rv**2) : normalized absorption cross-sections,
*     .........................................................................
*     .  Intermediate data:   file:  cm.out
*     .   Qext: extinction efficiency factors,
*     .   Qsca: scattering efficiency factors,
*     .   Qabs: absorption efficiency factors.
*     .........................................................................
*     .   Functions and subroutines used (alphabetically):
*...............................................................
*     .      functions: AKAPPA, COFF, fact2, fact3,
*     .                 OMEGA, ss1, TAU,
*     .                 P01AAF, X02AAF, A02ABF,
*     .      subroutines: asigma11,  CMLM1,    cmlm2,
*     .                   CXSC11b,   DASWF1,   delta21,
*     .                   cFEHL,     funlegnn, funlegnn1,
*     .                   gamkap1,   gj11,     homfunq,
*     .                   INT3,      INVMAT,   cdcof3,
*     .                   CMLM0,     sphtm,    sphte,
*     .                   CDRB12aa,  CDRf12,   DUGOLM1,
*     .                   DLEGF00,   DLEGF1,   DLEGF2,
*     .                   BESSJJ,    CESSEL0,  CDRG1bn
*     .                   cDRFP0,    cDRFP2,   cDRF222,
*     .                   cDRF2k02,  cDRKF45,  cDRKFS,
*     .                   cDRSF20,   cDRSF212, SEGV1,
*     .                   CDLAMn,    CDLAM_m,  sf02,
*     .                   F02AJF,    M01AAF,   F01AMZ,
*     .                   F02ANF,    A02ACF,   F01AMF,
*     .                   F01AMY,    A02AAF,   X04AAF,
*     .                   INVM01
c     .........................................................................

*     .      functions: AKAPPA, COFF,
*     .                 OMEGA,  TAU,
*     .      subroutines: BESSJJ,    cdcof4,
*     .                   CDLAMo_l,  CDLAMo_n,   CDLAMo_s,
*     .                   cdlamp,   CDRB12c,
*     .                   CDRf12c,    CESSEL,
*     .                (*)CGAMMA,    CMLM,
*     .                   CXSC10a,    DASWF1a, (*)DELSIG,
*     .                   DLEGF00,   DLEGF1,     DLEGF2,
*     .                   DLEGF22,   DRKF45,     DRKFS,
*     .                   DRFP0,     DRFP2,      DRF22,
*     .      ----?        DRF2K0,    DRSF20,     DRSF21,
*     .                   FEHL,      FUNLEG,     FUNLEG1,
*     .                (*)GAMKAP,    GJ,         INT3,
*     .                  (INVMAT,)   ks_new,
*     .                  (te,        te0,        te00,
*     .                   te2,       te_II,      te_III,)
*     .                  (tm,        tm0,        tm00,
*     .                   tm2,       tm_II,      tm_III,)
*     .                   cmtm,      cmte,
*     .                   Lambda,    Cmfunq,
*     .                   DELcc,     GAMKAP2,    sigma22
*     .                   GAMc2c1
*     .........................................................................
*     .........................................................................

      parameter (nalphas=19,nterms=330)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16(R-S)
      REAL*8 ksi1,ksi2
      COMPLEX*16 C2,c3

      DIMENSION Ialpha(nalphas), Jalpha(nalphas), Malpha(nalphas),
     *          ASIN(nalphas), ACOS(nalphas),
c     *          QXM(nalphas), PXM(nalphas), QSM(nalphas),
c     *          QXMM(nalphas), PXMM(nalphas), QSMM(nalphas)
c     *         ,QXE(nalphas), PXE(nalphas), QSE(nalphas),
c     *          PXEM(nalphas),QXEM(nalphas), QSEM(nalphas),
     *          QXM(nalphas), QSM(nalphas),
     *          QXMM(nalphas), QSMM(nalphas)
     *         ,QXE(nalphas), QSE(nalphas),
     *          QXEM(nalphas), QSEM(nalphas),
     *          QBKM(nalphas),  QBKMM(nalphas),
     *          QBKE(nalphas),  QBKEM(nalphas)

      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /K2/ Ialpha, Jalpha
      COMMON /EPS1/ EPS1
      COMMON /EPS3/ EPS3
      COMMON /EPS33/ EPS33
      COMMON /FF/ ASIN,ACOS
      COMMON /FACT/ FACT(300)
      COMMON /PI/ PI
      COMMON /Print/ iprint
  100 FORMAT(2F8.4)
  101 FORMAT(I4/(F8.4))
  102 FORMAT(I4)
  201 FORMAT(37X,4(D16.8,2X),D10.3,/)
  202 FORMAT(D12.4)
  204 FORMAT(1X,78('*'))
  205 FORMAT(1X,78('='))
  206 FORMAT(11X,'c',4X,'cksi1',2X,'Nmax',2X,'alpha',5X,'Qext',
     *       10X,'Qsca',10X,'Qabs')
  207 FORMAT(1X,'M=',I2,2X,2(F6.3,2X),2(I3,2X),1X,3(1pD13.6,1x),
     *       1X,'TM')
  277 FORMAT(34X,3(1pD13.6,1x))
  209 FORMAT(1X,78('.'))
  213 FORMAT(1X,78('-'))
  214 FORMAT(1X,'M=',I2,2X,2(F6.3,2X),2(I3,2X),1X,3(1pD13.6,1x),
     *       1X,'TE')
  220 FORMAT(2X,' ******* Nalph > nalphas *** Nalph =',I8)
c  224 FORMAT(1X,'alpha',3x,'m',6X,'Qext',11X,'Pext',
c     *       11X,'Qsca',11X,'Qabs')
  224 FORMAT(1X,'alpha',3x,'m',6X,'Qext',11X,'Qsca',
     *       11X,'Qabs',11X,'Qbk')

  225 FORMAT(1X,I3,4X,I2,4(1pD14.6,1X),1X,'TM')
 2259 FORMAT(1X,I3,4X,I2,4(1pD14.6,2X),1X,'TM')
22598 FORMAT(1X,f6.3,6(1pD12.5,1x))
22599 FORMAT(1X,f6.3,6X,4(1pD14.6,3X),1X,'TM')
  226 FORMAT(1X,I3,4X,I2,4(1pD14.6,1X),1X,'TE')
 2269 FORMAT(1X,I3,4X,I2,4(1pD14.6,2X),1X,'TE')
  227 FORMAT(10X,I4,10X,5(1pD12.3,1X),3X,2(D12.3,7X))
* -----------------------------------------------------------------
      print *,'C-M_SPH v.5_1q'
* -----------------------------------------------------------------
      open (unit=05,file='cmnn.dat',status='old',
     *      access='sequential')
      open (unit=07,file='cm.out', access='append',STATUS= 'UNKNOWN')
      open (unit=09,file='cm.res', access='append',STATUS= 'UNKNOWN')
      open (unit=77,file='cm.nrm', access='append',STATUS= 'UNKNOWN')
      open (unit=13,file='h3.ben',status='unknown',access='append')
* -----------------------------------------------------------------
*     setting program constants
* -----------------------------------------------------------------
      S=(0D0,0D0)
      S1=(0D0,1D0)
      PI = 4d0 * datan(1d0)
      DPI=2D0*PI
      PI180=PI/180D0
*       eps3 = 1d-80
*       eps3 = 1d-100
       eps3 = 1d-200
       eps33 = eps3
c factorial
       FACT(1)=1D0
       FACT(2)=1D0
       DO 130 I=3,170
  130  FACT(I)=FACT(I-1)*(I-1D0)
       FACT(170)=FACT(170)*(1.D-300)
       DO J = 170, 297
        FACT(J+1)=J*FACT(J)
       end do
* .................................................................
*     input data
* .................................................................
      READ(5,102) K
      write (*,*) k
      READ(5,100) AB
      write (*,*) ab
      READ(5,100) RII
      write (*,*) rii
      READ(5,100) RI
      write (*,*) ri
      READ(5,202) EPS,EPS1,eps2
      write (*,*)  EPS,EPS1,eps2
c!!
      if(k.eq.0) eps1 = 1d-10
      READ(5,202) eps22
      READ (5,102) iprint
      READ(5,100) WVL
      write (*,*)  wvl
      WRITE(7,202) EPS,EPS1,eps2
      WRITE(7,202) eps22
      READ (5,102) ICX
      GO TO (300,301,302,303,404),ICX
  300       READ (5,101) NC,ARV0,DRV
            WRITE (7,101) NC,ARV0,DRV,WVL
            XV0=DPI*ARV0/WVL
       GO TO 304
  301       READ (5,101) NC, XV0,DXV
            WRITE (7,101) NC, XV0,DXV,WVL
            ARV0=XV0*WVL/DPI
       GO TO 304
  302       READ (5,101) NC,  C0,DC
            WRITE (7,101) NC,  C0,DC ,WVL
            if(K.eq.0) ksi1=DSQRT(1.0D0+1.0D0/(AB**2-1.0D0))
            if(K.eq.1) ksi1=1D0/DSQRT(AB**2-1D0)
            if(K.eq.0) then
            Cksi1=C0*ksi1
            end if
            if(K.eq.1) then
            Cksi1=C0*ksi1*ab
            end if
            IF(K.EQ.0)   XV0=Cksi1/(AB)**(2D0/3D0)
            IF(K.EQ.1)   XV0=Cksi1*(AB)**(2D0/3D0)
            ARV0=XV0*WVL/DPI
       GO TO 304
  303       READ (5,101) NC,Cksi1,DCKSI
            WRITE (7,101) NC,Cksi1,DCKSI,WVL
c            WRITE (*,101) NC,Cksi1,DCKSI,WVL
            IF(K.EQ.0)   XV0=Cksi1/(AB)**(2D0/3D0)
            IF(K.EQ.1)   XV0=Cksi1*(AB)**(2D0/3D0)
            IF(K.EQ.1)   XV0 = Xv0 / AB
            ARV0=XV0*WVL/DPI
c            print *, xv0, arv0
c            pause
       GO TO 304
  404       READ (5,101) NC,XC0,DXC
            WRITE (7,101) NC,  XC0,DXC  ,WVL
  304 READ(5,102) nalph,IAL0,IDAL
      WRITE(7,102) nalph,IAL0,IDAL
      NAL=nalph
      IF(IAL0.EQ.0) NAL= nalph-1
  999 continue
      read(5,102,end=1000,err=1000) ISWITCH
      WRITE(7,*) 'K,AB',K,AB

 9000 continue
* ....................................
      if(iswitch.eq.0) then
      READ(5,*,end=1000) AB2
      if(AB.gt.AB2) then
      WRITE(*,199) ab2,ab
  199 format(1x,'a/b of core:',F8.4/1x,
     * 'a/b of mantle:',f8.4/1x,70('-'))
      stop 'a/b of core must be greater than a/b of mantle!'
      end if
      end if
* ....................................
      if(iswitch.eq.1) then
      read(5,*,end=1000) Vol_rat
      ARVcore=ARV0*Vol_rat**(1.D0/3.D0)
      x_v_core=XV0*Vol_rat**(1.D0/3.D0)
      if(Vol_rat.gt.1.) then
      write(*,245) ARVcore,ARV0,x_v_core,xv0,Vol_rat
  245 format(1x,'Rv of core:',f6.2/1x,
     * 'Rv of total particle:',f6.2/1x,
     * 'Xv of core:',f6.2/1x,
     * 'Xv of total particle:',f6.2/1x,
     * 'Core volume/Total volume:',f6.2/1x,70('-'))
      stop 'Core must be smaller than total particle!'
      end if
      end if
* ....................................
      if(iswitch.eq.2) then
      read(5,*,end=1000) x_v_core
      write(7,*) 'Xv of core',x_v_core
      ARVcore=x_v_core*WVL/DPI
      Vol_rat=(ARVcore/ARV0)**3.0D0
      if(Vol_rat.gt.1.) then
      write(*,245) ARVcore,ARV0,x_v_core,xv0,Vol_rat
      stop 'Core must be smaller than total particle!'
      end if
      end if
* ....................................
      if(iswitch.eq.3) then
      read(5,*,end=1000) ARVcore
      Vol_rat=(ARVcore/ARV0)**3.D0
      x_v_core=DPI*ARVcore/WVL
      if(ARVcore.gt.ARV0) then
      write(*,245) ARVcore,ARV0,x_v_core,xv0,Vol_rat
      stop 'Core must be smaller than total particle!'
      end if
      end if
* .................................................................
*     calculation of program constants
* .................................................................
      IF(K) 2,3,2
    2 CONTINUE
* oblate
      ksi1=1D0/DSQRT(AB**2-1D0)
      AKSI=ksi1**2+1D0
      AK=-1D0
* Vol_rat - The ratio of V_core to V_total
      if(iswitch.ne.0) then
         if(Vol_rat.eq.1) then
             ab2=ab
             ksi2=ksi1
             write(7,*) 'ab2',ab2
           else
             call ks_new(K,ksi1,Vol_rat,ksi2)
             ab2=dsqrt(1d0/ksi2/ksi2+1d0)
             write(7,*) 'ab2',ab2
         end if
      else
      KSI2=1D0/DSQRT(AB2**2-1D0)
      Vol_rat=sngl(ksi2*(ksi2*ksi2+1d0)/ksi1/(ksi1*ksi1+1d0))
      ARVcore=ARV0*Vol_rat**(1.D0/3.D0)
      x_v_core=XV0*Vol_rat**(1.D0/3.D0)
      end if
      GO TO 4
    3 CONTINUE
* prolate
      ksi1=DSQRT(1.0D0+1.0D0/(AB**2-1.0D0))

c          ksi1 = 1.1d0

      AKSI=ksi1**2-1D0
      AK= 1D0
* Vol_rat - The ratio of V_core to V_total
      if(iswitch.eq.0) then
      KSI2=DSQRT(1.0D0+1.0D0/(AB2**2-1.0D0))
      Vol_rat=sngl(ksi2*(ksi2*ksi2-1d0)/ksi1/(ksi1*ksi1-1d0))
      ARVcore=ARV0*Vol_rat**(1.D0/3.D0)
      x_v_core=XV0*Vol_rat**(1.D0/3.D0)
      else

         if(Vol_rat.eq.1) then
             ab2=ab
             ksi2=ksi1
             write(7,*) 'ab2',ab2
          else
             call ks_new(K,ksi1,Vol_rat,ksi2)
             ab2=dsqrt(1d0/(ksi2*ksi2-1d0)+1d0)
             write(7,*) 'ab2',ab2
         end if

      end if
    4 CONTINUE
      WRITE(7,204)
          IF(NALph.GT.nalphas) then
          WRITE(7,220) NALph
          WRITE(*,220) NALph
          STOP 'Nalph exceeds nalphas'
          end if
      IAL=IAL0
      IF(NAL.EQ.0) GO TO 666
****************      IF(NAL.EQ.0) stop 'NAL=0'
      IF(IAL.EQ.0) GO TO 66
      IAL=IAL0-IDAL
   66 CONTINUE
      DO 6 I=1,NAL
      IAL=IAL+IDAL
      ASIN(I)=DSIN(PI180*IAL)
      ACOS(I)=DCOS(PI180*IAL)
      IF(IAL.EQ.90) ACOS(I)=0D0
    6 Ialpha(I)=IAL
  666 CONTINUE

* .................................................................
*     start of calculations over size parameter
* .................................................................
      GO TO (305,306,307,308,409),ICX
  305       ARV=ARV0-DRV
            GO TO 309
  306        XV= XV0-DXV
            GO TO 309
  307        C1=  C0-DC
            GO TO 309
  308      CKSI=Cksi1-DCKSI
           dC=dCKSI/ksi1
           IF(K.EQ.1) dC=dC/AB
            GO TO 309
  409      XC  =  XC0-DXC
  309 CONTINUE

**
      IF(K.EQ.0) WRITE(13,203)
      IF(K.EQ.1) WRITE(13,208)
 8200 FORMAT(1X,'RI=',2F7.4,'i','/',2F7.4,'i,',
     *      1x,'a/b=',F6.2,1X,'/',F6.2,',',
     *       1x,'ksi=',F6.3,'/',F6.3,',',
     *       1X,'Vcore/Vtot=',F5.3)
      WRITE(13,8200) RII,RI,ab2,AB,ksi2,ksi1,Vol_rat
      WRITE(13,4206)
 2224 FORMAT(3X,'Xv',2x,'alpha',6x,'Cext/Srv(TM)',13X,'Cabs/Srv(TM)',
     *       13x,'Cext/Srv(TE)',13X,'Cabs/Srv(TE)')
      WRITE(13,2224)
      WRITE(13,4206)
 4206 FORMAT(1X,109('='))
**

      DO 1 I1=1,Nc

* -----------------------------------------------------------------
*      close(7)
*      close(9)
*      close(77)
      close(13)

*      open (unit=07,file='cm.out',status='old',
*     * access='append')
*      open (unit=09,file='cm.res',status='old',
*     * access='append')
*      open (unit=77,file='cm.nrm',status='old',
*     * access='append')
       open (unit=13,file='h3.ben',status='unknown',access='append')

  516 FORMAT(1X,78('*'))
      GO TO (310,311,312,313,414),ICX
  310       ARV=ARV+DRV
            XV=DPI*ARV/WVL
            if(iswitch.eq.2.or.iswitch.eq.3) then
            Vol_rat=(ARVcore/ARV)**3.D0
            call ks_new(K,ksi1,Vol_rat,ksi2)
            if(k.eq.0) ab2=dsqrt(1d0/(ksi2*ksi2-1d0)+1d0)
            if(k.eq.1) ab2=dsqrt(1d0/ksi2/ksi2+1d0)
            end if
            if(iswitch.eq.0.or.iswitch.eq.1) then
* that is Vol_rat is constant => we should recalculate x_v_core and ARVcore
            ARVcore=ARV*Vol_rat**(1.D0/3.D0)
            x_v_core=XV*Vol_rat**(1.D0/3.D0)
            end if
            GO TO 314
  311       XV=XV+DXV
*****************
c             if(i1.eq.3) Xv = 1d0
c             if(i1.eq.4) Xv = 3d0
c             if(i1.eq.5) Xv = 10d0
c             if(i1.eq.6) Xv = 30d0
*****************
            ARV=XV*WVL/DPI
            if(iswitch.eq.2.or.iswitch.eq.3) then
            Vol_rat=(ARVcore/ARV)**3.D0
            call ks_new(K,ksi1,Vol_rat,ksi2)
            if(k.eq.0) ab2=dsqrt(1d0/(ksi2*ksi2-1d0)+1d0)
            if(k.eq.1) ab2=dsqrt(1d0/ksi2/ksi2+1d0)
            end if
            if(iswitch.eq.1) then
* that is Vol_rat is constant => we should recalculate x_v_core and ARVcore
            ARVcore=ARV*Vol_rat**(1.D0/3.D0)
            x_v_core=XV*Vol_rat**(1.D0/3.D0)
            end if
  314       CONTINUE
            IF(K.EQ.0) then
            CKSI=XV*(AB)**(2D0/3D0)
            CKSI_core=x_v_core*AB2**(2D0/3D0)
c            xc_core=cksi_core/ab2
            end if
            IF(K.EQ.1) then
            CKSI=XV/(AB)**(2D0/3D0)
            CKSI_core=x_v_core/AB2**(2D0/3D0)
c            xc_core=cksi_core/ab2
            end if
            C1=CKSI/ksi1
            XC=CKSI/AB
            GO TO 315
  312       C1 =C1 +DC
            CKSI=C1*ksi1
            XC=CKSI/AB
            GO TO 316
  313 CKSI=CKSI+DCksi
      C1=CKSI/ksi1
      IF(K.EQ.1) C1=C1/AB
            XC=CKSI/AB
            GO TO 316
  414   XC=XC+DXC
        CKSI=XC*AB
      C1=CKSI/ksi1
      IF(K.EQ.1) C1=C1/AB
  316       CONTINUE
            IF(K.EQ.0)   XV=CKSI/(AB)**(2D0/3D0)
            IF(K.EQ.1)   XV=CKSI*(AB)**(2D0/3D0)
            IF(K.EQ.1)   XV=XV/AB
            x_v_core=XV*Vol_rat**(1.D0/3.D0)
            ARVcore=x_v_core*WVL/DPI
            IF(K.EQ.0) then
            CKSI_core=x_v_core*AB2**(2D0/3D0)
c            xc_core=cksi_core/ab2
            end if
            IF(K.EQ.1) then
            CKSI_core=x_v_core/AB2**(2D0/3D0)
            end if
            ARV=XV*WVL/DPI
  315 continue
      C2=RI*C1
      C3=RII*C1
      CCC=dmax1(cdabs(c2),cdabs(c3))
      if(k.eq.0) bcksi = cksi
      if(k.eq.1) bcksi = c1 * ksi1 * ab
      NMAX1 = 1.1*DABS(CCC) + 3 + 2 * k
c      NMAX1 = 1.2*DABS(CCC) + 3 + 2 * k
c      NMAX = idint(bcksi) + 5
      NMAX = idint(bcksi) + 3
c        print *, nmax, nmax1
       if(nmax1.gt.nmax) nmax = nmax1
c      NMAX=NMAX+6
      IF(KSI1.GT.2.5D0) NMAX= 2 * idint(xv) + 4
c
      IF(k.eq.0.and.ab.Gt.5D0.and.c1.ge.30d0) NMAX = nmax + 8
      IF(k.eq.0.and.ab.Ge.10D0.and.c1.gt.50d0) NMAX = nmax + 10
      IF(k.eq.0.and.ab.Ge.10D0.and.c1.gt.100d0) NMAX = nmax + 10
cc      IF(k.eq.0.and.ab.Ge.10D0.and.c1.gt.150d0) NMAX = nmax + 10
      IF(k.eq.1.and.ksi1.gt.0.5D0.and.c1.gt.2d0) NMAX=nmax+6
                              NmAX=nmax+4
*                              NmAX=nmax+10
*                                 NmAX=nmax+10
*                                 NmAX=nmax+10
      IF(MOD(NMAX,2)) 33,34,33
   33 NMAX=NMAX+1
   34 CONTINUE

c        print *, nmax
c        pause

* ...............................
* ??????
c       IF(k.eq.0.and.ksi2.lt.1.0001.or.k.eq.1.and.ksi2.lt.0.04d0) then
c       WRITE(*,*) 'ksi2 !!', ksi2
c       go to 1
c       end if
* ...............................

* ...............................
          IF(Nmax.GT.nterms) then
          WRITE(*,*) 'Nterms ! ', Nmax, nterms
c          STOP 'Nmax exceeds nterms'
          end if
* ...............................

      IF(K.EQ.0) WRITE(7,203)
  203 FORMAT(1X,'CORE / MANTLE SPHEROID: Prolate')
      IF(K.EQ.1) WRITE(7,208)
  208 FORMAT(1X,'CORE / MANTLE SPHEROID: Oblate')
      WRITE(7,200) RII,RI,wvl,ARVcore,ARV,
     *             x_v_core,XV,ab2,AB,ksi2,ksi1,Vol_rat
  200 FORMAT(1X,'Refractive index:',2F8.4,'i',1x,'/',2F8.4,'i'
     *      /1x,'Lambda=',F6.3,1x,' mkm,',
     *       2x,'Rv=  ',F6.3,1X,'/',F6.3,1x,'mkm,',
     *       2x,'Xv=',F6.3,1X,'/',F6.3,','
     *      /1x,'a/b=',F6.2,1X,'/',F6.2,',',
     *       2x,'ksi =',F6.3,1x,'/',F6.3,',',
     *       2X,' Volume: core / total= ',F5.3)
      WRITE(7,205)
      WRITE(7,206)
      WRITE(7,213)
**
      IF(K.EQ.0) WRITE(*,203)
      IF(K.EQ.1) WRITE(*,208)
      WRITE(*,200) RII,RI,wvl,ARVcore,ARV,
     *             x_v_core,XV,ab2,AB,ksi2,ksi1,Vol_rat
      WRITE(*,205)
      WRITE(*,206)
      WRITE(*,213)
* ...............................

      M=1
      NK=1.5*NMAX+10
      NK=NMAX+40
      IF(NK.LT.60) NK=60
      IF(ksi1.gt.1.5D0) NK=NK+40
      IF(NK.GT.4*nterms) NK=4*nterms
      IF(NMAX.GT.nterms) NMAX=nterms
      NMXE = NMAX + 6
      IF(DABS(CCC).LT.4D0) NMXE=NMXE+2
c      IF(AB.GT.4D0) NMXE=NMXE+4
c      IF(AB.GT.4D0.and.k.eq.1) NMXE=NMXE+4
*************************************************************
c      IF(C1.lT.4D0.and.k.eq.1) NMXE = NMXE + 8
c      IF(C1.lT.4D0.and.k.eq.1) NMXE = NMXE + 8
c      IF(C1.lT.4D0.and.k.eq.0.and.ksi1.gt.2d0) NMXE = NMXE + 10
      IF(NMXE.GT.nterms) NMXE=nterms
      if(cksi.lt.3d0.and.ksi1.gt.2d0) nmxe=nmax
      if(ksi1.gt.2.5d0) nmxe=nmax

      DO 7 I=1,NAL
      QbkE(I)=0D0
      Qbkm(I)=0D0
      QXE(I)=0D0
      QXM(I)=0D0
      QSE(I)=0D0
      QSM(I)=0D0
c      PXE(I)=0D0
c      PXM(I)=0D0
      Jalpha(I)=0
    7 CONTINUE
*-------------
      IF(IAL0.EQ.0) GO TO 50
*-------------

   54 CONTINUE
      JAM=0
* .................................................................
*     alpha is NOT equal 0 degr.:
*     calculation of efficiency factors
* .................................................................

   11 CONTINUE

*-------------------------------------------------------------
* <TM-mode>

c      CALL Cmfunq2(1,m,c1,c2,c3,KSI1,ksi2,EPS,NMXE,*12)
      CALL Cmfunq3(1,m,c1,c2,c3,KSI1,ksi2,NMXE,*12)
      call cmtm(1,m,rii,ri,c1,ksi1,KSI2,nmax,nmxe)
c      call CXSC10a(1,m,C1,NMXE,NMAX,QXMM,PXMM,QSMM,QEXT,QSCA)
      CALL CXSC11b(1,M,C1,NMXE,NMAX,QXMM,QSMM,qbkmm,QEXT,QSCA,qbk)

      JM=JAM+1

      DO   I = JM, NAL
        QXM(I)=QXM(I)+QXMM(I)
        QbkM(I)=QbkM(I)+QbkMM(I)
        QSM(I)=QSM(I)+QSMM(I)
c        PXM(I)=PXM(I)+PXMM(I)
        QABS=QXM(I)-QSM(I)
        QABSm=QXMm(I)-QSMm(I)
c        WRITE(7,207) M,C1,CKSI,NMAX,Ialpha(I),QXM(I),QSM(I),QABS,QbkM(I)
c        WRITE(7,277) QXmM(I),QSmM(I),QABSm,QbkmM(I)
        WRITE(7,207) M,C1,CKSI,NMAX,Ialpha(I),QXM(I),QSM(I),QABS
        WRITE(7,277) QXmM(I),QSmM(I),QABSm
        WRITE(*,277) QXmM(I),QSmM(I),QABSm
        WRITE(*,207) M,C1,CKSI,NMAX,Ialpha(I),QXM(I),QSM(I),QABS
      end do

*-------------------------------------------------------------
* <TE-mode>

      call cmte(m,c1,ksi1,KSI2,nmxe)
c      call CXSC10a(2,m,C1,NMXE,NMXe,QXeM,PXeM,QSeM,QEXT,QSCA)
      CALL CXSC11b(2,M,C1,NMXE,NMXE,QXEM,QSEM,QbkEM,QEXT,QSCA,Qbk)

      DO   I = JM, NAL
        QXE(I) = QXE(I) + QXEM(I)
        QBKE(I) = QBKE(I) + QBKEM(I)
        QSE(I) = QSE(I) + QSEM(I)
c        PXE(I) = PXE(I) + PXEM(I)
        QABS =   QXE(I) - QSE(I)
        QABSm = QXEm(I) - QSEm(I)
c        WRITE(7,214) M,C1,CKSI,NMXE,Ialpha(I),QXE(I),QSE(I),QABS,QBKE(I)
c        WRITE(7,277) QXeM(I),QSeM(I),QABSm,QBKEm(I)
        WRITE(7,214) M,C1,CKSI,NMXE,Ialpha(I),QXE(I),QSE(I),QABS
        WRITE(7,277) QXeM(I),QSeM(I),QABSm
        WRITE(*,277) QXeM(I),QSeM(I),QABSm
        WRITE(*,214) M,C1,CKSI,NMXE,Ialpha(I),QXE(I),QSE(I),QABS

        IF(DABS(QXEM(I)/QXE(I)).LE.EPS2.and.
     *     DABS(QXMM(I)/QXM(I)).LE.EPS2) Jalpha(I)=1
        JAM=JAM+Jalpha(I)
        Malpha(I)=M
        IF(JAM.EQ.NAL) GO TO 13
      end do
*-------------------------------------------------------------
      WRITE(7,209)
      WRITE(*,209)
      M=M+1

      GO TO 11
   13 WRITE(7,213)
      WRITE(*,213)
      GO TO 53

   50 CONTINUE
* -----------------------------------------------------------------
*      alpha is equal to 0 degrees
* -----------------------------------------------------------------
      QEX1=1D0
      inv = 0

c      IF(NMAX.gt.nterms) NMAX=nterms
      if(nal.eq.0) NMXE = NMAX

* -----------------------------------------------------------------
*     calculation of efficiency factors
* -----------------------------------------------------------------
   55 CONTINUE
c      CALL Cmfunq2(0,1,c1,c2,c3,KSI1,ksi2,EPS,NMXE,*12)
      CALL Cmfunq3(0,1,c1,c2,c3,KSI1,ksi2,NMXE,*12)
      call cmtm(0,1,rii,ri,c1,ksi1,KSI2,nmax,nmxe)
c      call CXSC10a(0,1,C1,NMXE,NMAX,QXMM,PXMM,QSMM,QEXT,QSCA)
      CALL CXSC11b(0,1,C1,NMXE,NMAX,QXMM,QSMM,Qbkm,QEXT,QSCA,qbk)
      QAB=QEXT-QSCA
      qQAB = QSCA / QEXT - 1d0
c      PEXT=0D0

c      WRITE(7,207) M,C1,CKSI,NMAX,IAL0,QEXT,QSCA,QAB,qbk
      WRITE(7,207) M,C1,CKSI,NMAX,IAL0,QEXT,QSCA,QAB
      WRITE(*,207) M,C1,CKSI,NMAX,IAL0,QEXT,QSCA,QAB
cccccccccccccccccccccccccccccccccccccccc
c      IF(DABS(QEXT/QEX1-1D0).LE.EPS2.and.qAB.gt.-0.1d0) GO TO 52
      IF(DABS(QEXT/QEX1-1D0).LE.eps22.and.qAB.gt.-0.1d0.or.
     *        DABS(qQAB).lt.eps22.or.inv.ge.12) GO TO 52
      inv = inv + 1
      IF(DABS(QAB).lt.eps22) GO TO 52
      IF(DABS(QAB).lt.eps2) GO TO 52
*      IF(DABS(qAB).gt.1d3.and.nmax.gt.nterms) GO TO 1
      QEX1=QEXT
      NMAX = NMAX + 2
c ** DRAINE+N! ???
       IF(qAB.lt.-1d0) NMAX=NMAX+4
       IF(qAB.lt.-10d0) NMAX=NMAX+4
       IF(qAB.lt.-20d0) NMAX=NMAX+10
      NMXe = NMXe + 4
c      IF(c1.lT.4D0.and.ab.le.5d0.and.k.eq.1) NMXe=NMXe+2*ab
c      IF(c1.lT.4D0.and.ab.gt.5d0.and.k.eq.1) NMXe=NMXe+30
      if(cksi.lt.5d0.and.ksi1.gt.1.5d0)  nmxe=nmax
c      IF(NMXE.GT.nterms) NMXE=nterms
      IF(NMaX.GT.nterms) go to 52
      if(nal.eq.0) NMXE=NMAX
      GO TO 55
   52 CONTINUE
      WRITE(7,205)
      WRITE(*,205)

      IF(NAL.GT.0) GO TO 54
   53 CONTINUE
* -----------------------------------------------------------------
*     print of final results
* -----------------------------------------------------------------
      IF(K.EQ.0) write(9,203)
      IF(K.EQ.1) write(9,208)
      WRITE(9,200) RII,RI,wvl,ARVcore,ARV,
     *             x_v_core,XV,ab2,AB,ksi2,ksi1,Vol_rat
      write(9,209)
      write(9,223) C1,C2,c3,CKSI,cksi_core,NMAX,NMXE
  223 FORMAT(1X,'C1=',f5.2,1X,'C2=',2f5.2,1X,'C3=',2f5.2,1X,
     *       'Cksi=',f5.2,1X,'Cksi2=',f5.2,1X,
     *       'Nmax=',I3,1X,'Nmxe=',I3)

      write(9,213)
      write(9,224)
      write(9,213)

      IF(K.EQ.0) write(77,203)
      IF(K.EQ.1) write(77,208)
      WRITE(77,200) RII,RI,wvl,ARVcore,ARV,
     *             x_v_core,XV,ab2,AB,ksi2,ksi1,Vol_rat

      write(77,209)
      write(77,223) C1,C2,c3,CKSI,cksi_core,NMAX,NMXE

      write(77,213)
      write(77,2249)
      write(77,209)

 2249 FORMAT(1X,'alpha',3x,'m',1X,'Cext/(Pi*Rv**2)',
c     *       1X,'Csca/(Pi*Rv**2)',1X,'Cabs/(Pi*Rv**2)')
     *  1X,'Csca/(Pi*Rv**2)',1X,'Cabs/(Pi*Rv**2)'1X,'Cbk/(Pi*Rv**2)')

      M0=1
c      IF(IAL0.EQ.0) write(9,225) IAL0,M0,QEXT,PEXT,QSCA,QAB
      IF(IAL0.EQ.0) write(9,225) IAL0,M0,QEXT,QSCA,QAB,qbk
      IF(IAL0.EQ.0) then
       yyy=ab**(2d0/3d0)
       if(k.eq.0) then
       g1=qext/yyy
       g2=qsca/yyy
       g3=qab/yyy
       gbk=qbk/yyy
       else
       g1=qext*yyy
       g2=qsca*yyy
       g3=qab*yyy
       gbk=qbk*yyy
       end if
       write(77,2259) IAL0,M0,g1,g2,g3,gbk
       WRITE(13,22252) xv,ial0,g1,g3
22252 FORMAT(f6.2,i4,6(1pD25.15))
***********       write(77,22599) xv,g1,g2,g3
*****************       write(77,22599) cksi,QEXT,QSCA,QAB
      end if

      IF(nal.eq.0) GO TO 114

      DO 112 IA=1,NAL
      QABSm=QXM(IA)-QSM(IA)

* normalization
      if(k.eq.0) then
      ttt=dsqrt(ab*ab*asin(ia)*asin(ia)+acos(ia)*acos(ia))/
     * (ab)**(2d0/3d0)
      else
      ttt=dsqrt(ab*ab*acos(ia)*acos(ia)+asin(ia)*asin(ia))/
     * (ab)**(1d0/3d0)
      end if
      u1=QXM(ia)*ttt
      u2=QSM(ia)*ttt
      u3=u1-u2
      ubk=qbkm(ia)*ttt
      write(77,2259) Ialpha(IA),Malpha(IA),u1,u2,u3,ubk
********       write(77,22598) cksi,QXM(IA),QSM(IA),QABS,
********     *                 QXE(IA),QSE(IA),QXE(IA)-QSE(IA)
c      write(9,225) Ialpha(IA),Malpha(IA),QXM(IA),PXM(IA),QSM(IA),QABS
      write(9,225) Ialpha(IA),Malpha(IA),QXM(IA),QSM(IA),QABS,qbkm(ia)
  112 CONTINUE

      write(9,209)
      write(77,209)

      DO 412 IA=1,NAL
      QABSe=QXE(IA)-QSE(IA)
* normalization
      if(k.eq.0) then
      ttt=dsqrt(ab*ab*asin(ia)*asin(ia)+acos(ia)*acos(ia))/
     * (ab)**(2d0/3d0)
      else
      ttt=dsqrt(ab*ab*acos(ia)*acos(ia)+asin(ia)*asin(ia))/
     * (ab)**(1d0/3d0)
      end if
      w1=QXE(ia)*ttt
      w2=QSE(ia)*ttt
      w3=w1-w2
      wbk=qbke(ia)*ttt
      u1=QXM(ia)*ttt
      u2=QSM(ia)*ttt
      u3=u1-u2
      write(77,2269) Ialpha(IA),Malpha(IA),w1,w2,w3,wbk
      WRITE(13,22252) xv,Ialpha(IA),u1, u3, w1, w3
c      write(9,226) Ialpha(IA),Malpha(IA),QXE(IA),PXE(IA),QSE(IA),QABSe
      write(9,226) Ialpha(IA),Malpha(IA),QXE(IA),QSE(IA),QABSe,qbke(ia)
  412 CONTINUE

  114 write(9,205)
      write(77,205)
      WRITE(13,4206)
      GO TO 5
   12 continue
    5 CONTINUE

*      close(7)
*      close(9)
*       close(77)
* -----------------------------------------------------------------
    1 CONTINUE
      GO TO 999
 1000 continue
      close(7)
      close(9)
      close(77)
      STOP
      END

c================================================
c                 FUNCTIONS
c================================================
c************************************************
C AKAPPA
C
      FUNCTION AKAPPA(M,L,N,NN,A,B,C,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NN),B(NN)
      COMMON /EPS3/ EPS3
      IF(MOD(IABS(L-N),2)) 20,20,9
   20 AKAPPA =0D0
      RETURN
    9 CONTINUE
      M2=2*M
      IF(MOD(L-M,2))1,1,2
    1 I1=0
      I2=1
      K=M2
      E=A(1)*(-B(1)*M/(M2+3D0))*2D0
      GO TO 3
    2 I1=1
      I2=0
      K=M2+1
      E=A(1)*(B(1)*(M+2D0)/(M2+1D0)-B(2)*(M+1D0)*(M2+2D0)/
     *  (M2+5D0))*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E=E*F
      DO 5 I=2,NN-1
      J=2*I-1-I2
      JM=J+M
      JM2=2*(J+M)
c      F=1D0
      F=a(i)
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
      IF(M2.EQ.0) F=1D0
      G=(B(I-I2)*J*(JM+1D0)/(JM2-1D0)-B(I+I1)*JM*(M2+J+1D0)/
     *  (JM2+3D0))*2D0/(JM2+1D0)*F
      IF(DABS(G).LE.eps3.and.i.gt.l+n) GO TO 6
    5 E=E+G
    6 AKAPPA =E/C/D
      RETURN
      END

c************************************************
C COFF
C
      FUNCTION COFF(L,N,A)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N)
      IF(MOD(L-1,2))1,1,2
    1 I2=1
      GO TO 3
    2 I2=0
    3 D=0D0
      DO 4 I=1,N
      J=2*I-1-I2
    4 D=D+(J+1D0)*(J+2D0)*A(I)
      COFF=D
      RETURN
      END

c************************************************
C OMEGA
C
      FUNCTION OMEGA(M,L,N,NN,A,B,C,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NN),B(NN)
      COMMON /EPS3/ EPS3
      IF(MOD(IABS(L-N),2)) 9,9,20
   20 OMEGA = 0D0
      RETURN
    9 CONTINUE
      M2=2*M
      IF(MOD(L-M,2))1,1,2
    1 I2=1
      K=M2
      E=A(1)*(B(2)*(M2+2D0)*(M2+1D0)/(M2+5D0)/(M2+3D0)-
     *  B(1)*2D0*(M+1D0)*(M2-1D0)/(M2+3D0)/(M2-1D0))*2D0/(M2+1D0)
      GO TO 3
    2 I2=0
      K=M2+1
      E=A(1)*(B(2)*(M2+3D0)*(M2+2D0)/(M2+7D0)/(M2+5D0)-B(1)*2D0*
     *  (M+1)*(M2+1D0)/(M2+5D0)/(M2+1D0))*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E=E*F
      DO 5 I=2,NN-1
      J=2*I-1-I2
      JM=J+M
      JM2=2*(J+M)
c      F=1D0
      F=a(i)
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
      IF(M2.EQ.0) F=1D0
      G = (B(I-1)*J*(J-1D0)/(JM2-1D0)/(JM2-3D0)+B(I+1)*(J+M2+2D0)*
     *  (J+M2+1D0)/(2D0*JM+5D0)/(2D0*JM+3D0)-B(I)*2D0*(JM*(JM+1D0)+M*M
     *  -1D0)/(2D0*JM+3D0)/(2D0*JM-1D0))*2D0/(JM2+1D0)*F
      IF(DABS(G).LE.eps3.and.i.gt.l+n) GO TO 6
    5 E=E+G
    6 OMEGA =-E/C/D
      RETURN
      END

c************************************************
C TAU
C
      FUNCTION TAU(M,L,N,NN,A,B,C,D)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NN),B(NN)
      COMMON /EPS3/ EPS3
      IF(MOD(IABS(L-N),2)) 9,9,20
   20 TAU    =0D0
      RETURN
    9 CONTINUE
      M2=2*M
      IF(MOD(L-M,2))1,1,2
    1 I2=1
      K=M2
      E=A(1)*B(1)*M*(M+1D0)*2D0/(M2+1D0)
      GO TO 3
    2 I2=0
      K=M2+1
      E=A(1)*B(1)*(M+1D0)*(M+2D0)*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E=E*F
      DO 5 I=2,NN
      J=2*I-1-I2
      JM=J+M
c      F=1D0
      F=a(i)
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
      IF(M2.EQ.0) F=1D0
      G = B(I)*JM*(JM+1D0)*2D0/(2D0*JM+1D0)*F
      IF(DABS(G).LE.eps3.and.i.gt.l+n) GO TO 6
    5 E=E+G
    6 TAU    =E/C/D
      RETURN
      END

c************************************************
C fact2
C
      real FUNCTION fact2*8(n,k)
      REAL*8 f
      f = 1d0
       do i = k + 1, n
        f = f * i
       end do
      fact2 = f
      RETURN
      END

c************************************************
C fact3
C
      real FUNCTION fact3*8(n,k,f1)
      REAL*8 f, f1
      f = f1
       do i = k + 1, n
        f = f * i
       end do
      fact3 = f
      RETURN
      END

c************************************************
C ss1
C
      complex*16 FUNCTION ss1(n)
      complex*16 s, s1, ss0
      REAL*8 aksi, ak
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      Ss0=(1D0,0D0)
      n1 = mod(iabs(n),4)
      go to (1, 2, 3, 4) n1+1
    1  ss1 = ss0
       RETURN
    2  ss1 = s1
       if(n.lt.0) ss1 = -ss1
       RETURN
    3  ss1 = -ss0
       RETURN
    4  ss1 = -s1
       if(n.lt.0) ss1 = -ss1
       RETURN
      END

c================================================
c                 SUBROUTINES
c================================================

c************************************************
C BESSJJ0
c
      SUBROUTINE BESSJJ0(A,NUM,BESJ)
c+&      parameter (nterm=200)
      IMPLICIT COMPLEX*16 (A-H,O-Q,T-Z)
      DIMENSION BESJ(NUM+1)
   11 continue
c+&
c      BESJ(NUM+1)=0D0
          do i =1,num+1
            BESJ(i)=(0D0,0D0)
          end do
      BESJ(NUM)=(1D-300,1d-300)
      N=2*NUM+1
      NUM1=NUM-1
      DO I = 1, NUM1
       N=N-2
       I1=NUM-I
       BESj(I1) = n * a * BESj(I1+1) - BESj(I1+2)
         if(cdabs(besJ(i1)).gt.1d300) then
           num = i - 1
           go to 11
         end if
      end do
      N=2*(NUM/2)
      B=1.2533141373155002D0*CDSQRT(A)
      C=B*BESJ(1)
      DO 12 I=3,N,2
      B=B*(I-0.5D0)*(I-2.0D0)/(I-2.5D0)/(I-1.0D0)
   12 C=C+B*BESJ(I)
      C=1.0D0/C
      DO 13 I=1,NUM
   13 BESJ(I)=C*BESJ(I)
      RETURN
      END

c************************************************
C cdcof4a
C
      SUBROUTINE cdcof4a(VL,AD,IA,BD,IB,ID,VN,IN,M,N,C,IER)
      parameter (nterms=330)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 VL,AD,VN,C,bd
      DIMENSION AD(IA),BD(IB),RL(4*nterms)
      COMMON /EPS3/ EPS3
      COMMON /EPS33/ EPS33
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      common /FACT/ FACT(300)
      AC=C
      IER=0
c      eps33 = eps3
c      eps4 = 1d-20
                 eps4 = 1d-40
C* èPOBEPKA BBOÑàMõX áHAóEHàâ èAPAMETPOB
      IF(M.LT.0.OR.N.LT.M) IER=1
      IF(AC.LT.0D0) IER=3
      IF(ID.EQ.1.AND.IA.LE.((N-M)/2+2)) IER = 4
      IF(ID.EQ.0.AND.IN.NE.0) IER = 6
      IF(IER.NE.0) GO TO 900
      aM = M
      aMM = 2d0 * M
      S2=(1D0,0D0)
      bM = 4D0 * am**2 - 1D0
      RC2=C*C
      RC4=RC2*RC2
      N1=N-M
      IF(ID.EQ.0) GO TO 900
C**** BõóàCãEHàE K-TOB PAáãOÜEHàü
      IF(MOD(N1+1,2)) 21,20,21
   20 I11=1
      GO TO 22
   21 I11=0
   22 CONTINUE
      L1=(N1-I11)/2
      L11=L1+1
      RL(L11)=S2
      IF(L1.EQ.0) GO TO 24
C***  BõóàCãEHàE K-TOB èPà R<N-M
      DO 23 J=1,L1
      JJ1=N1-2*J
      RU=S2
      RB=(am+JJ1)*(am+JJ1+1D0)+AK*RC2/2D0*(1D0-bm/(amm+2D0*JJ1-1D0)/
     *   (amm+2D0*JJ1+3D0))        -VL
      RV=S2/RB
      RW=RV
      JJ2=JJ1
      IF(JJ2.LE.1) GO TO 103
  102 CONTINUE
      II=2*M+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(am+JJ2-2D0)*(am+JJ2-1D0)+AK*RC2/2D0*(1D0-bm/(II1-5D0)/
     *   (II1-1D0))        -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV
      IF(CDABS(RV/rw).LE.eps4.and.j.ge.l1/2) GO TO 103
c      IF(CDABS(RV/rw).LE.eps4) GO TO 103
      JJ2=JJ2-2
      IF(JJ2.LE.1) GO TO 103
      GO TO 102
  103 RL(L11-J)=-(amm+JJ1+1D0)*(amm+JJ1+2D0)*AK*RC2/(amm+2D0*JJ1+3D0)/
     *           (amm+2D0*JJ1+5D0)*RW
   23 RL(L11-J)=RL(L11-J)*RL(L11+1-J)
   24 CONTINUE
C***  BõóàCãEHàE K-TOB èPà R>N-M
      IA1=IA-1
      DO 25 J=L11,IA1
      JJ1=N1+2*(J-L11)+2
      II=2*M+JJ1
      II1=II+JJ1
      RU=S2
      RB=(am+JJ1)*(am+JJ1+1D0)+AK*RC2/2D0*(1D0-bm/(II1-1D0)/(II1+3D0))
     *   -VL
      RV=-JJ1*(JJ1-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RW=RV
      JJ2=JJ1
  104 CONTINUE
      JJ2=JJ2+2
      II=2*m+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(am+JJ2)*(am+JJ2+1D0)+AK*RC2/2D0*(1D0-bm/(II1-1D0)/(II1+3D0))
     *          -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV

c       IF(CDABS(RV/rw).LE.eps4.and.j.ge.ia1/2) GO TO 105
       IF(CDABS(RV/rw).LE.eps4) GO TO 105
c       IF(CDABS(RV/rw).LE.eps3) GO TO 105
      GO TO 104
  105 RL(J+1)=RW*(amm+2D0*JJ1-1D0)*(amm+2D0*JJ1+1D0)/((amm+JJ1-1D0)*
     *      (amm+JJ1)*AK*RC2)
      RL(J+1)=RL(J)*RL(J+1)
       IF(J.GT.L11) RL(J+1)=RL(J+1)/RL(L11)
       IF(CDABS(RL(J+1)).LT.eps3) GO TO  125
   25 CONTINUE
      GO TO 126
  125 IJ=J+2
         DO 127 I=IJ,IA
  127    RL(I)=S
  126 CONTINUE

C *                  BõóàCãEHàE K-TA D(N-M)
      L4=(N+M+I11)/2
c      if(l1+1.gt.170) stop 1701
      F2=1D0/Fact(l1+1)
      DO  I = 1, L4
       F2 = F2 * (L4 + I)
       if(dabs(f2).gt.1d300) go to 27
      end do
      iff = 0
      go to 270
   27 continue
      f2 = f2 * 1d-300
      iff = 1
c Draine!
c      DO Ij = i, l4
      DO Ij = i+1, l4
c       F2 = F2 * (J3 + Ij)
       F2 = F2 * (L4 + Ij)
      end do
 270  continue
      il1 = 1
      if(mod(l1,2).ne.0) il1 = -1
      D1 = iL1 * F2 / 2D0**N1
       if(iff.eq.1) d1 = d1 * 1d300
      SNM=S

      DO 28 J=1,IA
      J1=2*(J-1)+I11
      J2=(J1-I11)/2
      J3=(J1+2*M+I11)/2
      rF2=RL(J)/2D0**J1
      DO I=1,J3
       rF2 = rF2 * (J3 + I)
       if(cdabs(rf2).gt.1d300) go to 30
      end do
      iff = 0
      go to 330
   30 continue
      rf2 = rf2 * 1d-300
      iff = 1
c DRAINE+!
c      DO Ij = i, J3
      DO Ij = i+1, J3
       rF2 = rF2 * (J3 + Ij)
      end do
 330  continue
      rF1=rF2
      DO 29 I=1,J2
   29 rF1=rF1/I
       SG=(-1D0)**J2*rF1
       if(iff.eq.1) SG = sg * 1d300
       IF(CDABS(SG).LT.eps3.AND.J.GT.ia/2) GO TO 128
   28 SNM=SNM+SG
  128 AD(L11)=D1/SNM

      DO 31 I=1,IA
   31 AD(I)=RL(I)*AD(L11)
      IF(ID.NE.2) GO TO 901
C*       BõóàCãEHàE KOùîîàñàEHTOB PAáãOÜEHàü èPà  R<0
      BD(1)=AD(1)
      DO 223 J=2,IB
      JJ1=I11-2*(J-1)
      RU=S2
      RB=(am+JJ1)*(am+JJ1+1D0)+AK*RC2/2D0*(1D0-bm/(amm+2D0*JJ1-1D0)/
     *   (amm+2D0*JJ1+3D0))        -VL
      RV=1D0/RB
      RW=RV
      JJ2=JJ1
  224 CONTINUE
      II=2*M+JJ2
      II1=II+JJ2
      RR=-JJ2*(JJ2-1D0)*II*(II-1D0)*RC4/((II1-1D0)**2*(II1-3D0)*
     *   (II1+1D0))/RB
      RB=(am+JJ2-2D0)*(am+JJ2-1D0)+AK*RC2/2D0*(1D0-bm/(II1-5D0)/
     *   (II1-1D0))        -VL
      RR=RR/RB
      RU=1D0/(1D0+RR*RU)
      RV=RV*(RU-1D0)
      RW=RW+RV
c          print *, j, RV/rw
c          pause
c      IF(CDABS(RV/rw).LE.eps4.and.j.ge.ib/2) GO TO 225
      IF(CDABS(RV/rw).LE.eps4) GO TO 225
      JJ2=JJ2-2
      GO TO 224
  225 CONTINUE
      IF(JJ1.EQ.-2*M-2.OR.JJ1.EQ.-2*M-1) GO TO 226
      BD(J)=-(amm+JJ1+1D0)*(amm+JJ1+2D0)*AK*RC2/(amm+2D0*JJ1+3D0)/
     *           (amm+2D0*JJ1+5D0)*RW
      GO TO 227
  226 BD(J)= AK*RC2/(amm+2D0*JJ1+3D0)/(amm+2D0*JJ1+5D0)*RW
      IF(I11.EQ.1) BD(J)=-BD(J)
  227 BD(J)=BD(J)*BD(J-1)
          IF(cDABS(BD(J)).LT.eps33) GO TO 1223
  223 CONTINUE
        GO TO 901
 1223 IJ=J+1
        DO 1224 I=IJ,IB
 1224   BD(I)=s
  901 IF(IN.EQ.0) GO TO 900
C*****  BõóàCãEHàE HOPMàPYûôEÉO MHOÜàTEãü *
      VN=S
         iff = 0
         if(cdabs(ad(1)).gt.1d150) iff = 1

      DO 32 J=1,IA
      JJ1=2*(J-1)+I11
         rf1 = 1d0
         if(iff.eq.1) rf1 = 1d-300
      rF1 = rf1 * AD(J) / (amm+2D0*JJ1+1D0) * ad(j)
      IF(M.EQ.0) GO TO 333
      DO 33 I=1,2*M
   33 rF1=rF1*(JJ1+I)
  333 RVN=rF1
      VN=VN+RVN
      IF(CDABS(RVN/(vn+1d-50)).LT.eps4.and.j.gt.ia/2) GO TO 34
c      IF(CDABS(RVN/(vn+1d-50)).LT.eps4) GO TO 34
   32 continue
   34 VN=2D0*VN
         if(iff.eq.1) vn = vn * 1d300
  900 RETURN
      END

c************************************************
c cdrb12cc
c
      SUBROUTINE CDRB12cc(R,DR,AR,ADR,ID,M,N,C,X,AD,IA,IER)
c      parameter (nterm=200)
      parameter (nterm=500)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 AD,C,S1,BESJ,BESY,AC,SIG,
     *           RK,R,DR,F2,RS,RD,S,BR,ABR,ar,adr
      COMPLEX*16 sS1
      DIMENSION AD(IA),BESJ(4*nterm),BESY(4*nterm)
      DIMENSION iiy(4*nterm)
      COMMON /EPS3/ EPS3
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /pi/ pi
      factor=1d250
      IER=0
      IF(X.LT.1D-10) IER=10
      IF(IER.NE.0) GO TO 900
      M2=2*M
      IF(MOD((N-M),2)) 1,2,1
    1 JK=1
      GO TO 5
    2 JK=0
    5 CONTINUE
      AM=M
      AM=AM/2D0
      PI2 = pi / 2d0
      X2=X*x
      AC=C*X
      SIG=S
      DO 7 I=1,IA
      J=2*(I-1)+JK
      F2=AD(I)
      DO 8 I1=1,M2
    8 F2=F2*(J+I1)
    7 SIG=SIG+F2
      RK=CDSQRT(PI2/AC)/SIG*((X2-AK)/X2)**AM
       NUM=3*CDABS(AC)+10
c       if(x.gt.2d0) NUM = 3 * CDABS(AC) + 10 * (m + x)
       if(x.gt.2d0) NUM = 3 * CDABS(AC) + 6 * (m + x)
      if(NUM.lt.60) num=60
      if(NUM-1.ge.4*nterm) write(*,*) num, 4*nterm
      if(NUM-1.ge.4*nterm) num = 4*nterm
      IF(ID.EQ.1) CALL CESSEL0(1D0/AC,NUM,BESJ,BESY,iiy)
      IF(ID.EQ.0) CALL BESSJJ0(1D0/AC,NUM,BESJ)
      R=S
      DR=S
      AR=s
      ADR=s
      F1=AK/X*M/(X2-AK)
      NM2=MIN0(IA,NUM/2-2)

      DO 13 I=1,NM2
      J=2*(I-1)+JK
c      F2=AD(I)*S1**(J+M-N)
      F2=AD(I)*sS1(J+M-N)
      DO 18 I1=1,M2
   18 F2=F2*(J+I1)
      IF(I.EQ.1) GO TO 21
c+&      IF(CDABS(RD).LT.eps3.AND.CDABS(RS).LT.eps3.AND.I.GT.(M+5)
c+&     *   .AND.ID.EQ.0)   GO TO 19
   21 RS=F2*BESJ(J+M+1)
      RD=F2*((F1+(M+J)/X)*BESJ(J+M+1)-C*BESJ(J+M+2))
      R=R+RS
      DR=DR+RD
   20 CONTINUE
      IF(ID.EQ.0) GO TO 13
      IF(I.EQ.1) GO TO 22
c+&      IF(CDABS(BR).LT.eps3.AND.CDABS(ABR).LT.eps3.AND.I.GT.(M+5))
c+&     *   GO TO 14
   22 BR = BESY(J+M+1) * F2 * factor**(iiy(j+m+1))
                     AR  = AR  + BR
      ABR=(F1 +(M+J)/X)*BESY(J+M+1) * f2* factor**(iiy(j+m+1))
     *    -C*BESY(J+M+2) * f2 * factor**(iiy(j+m+1))
                     ADR = ADR + ABR
   13 CONTINUE
   14 AR= AR* RK
      ADR= ADR* RK
   19 R=R*RK
      DR=DR*RK
  900 RETURN
      END

c************************************************
c CDRF12cc
c
      SUBROUTINE CDRF12cc(R,DR,AR,ADR,I12,M,N,C,AD,IA,BD,IB,IER)
      parameter (nterms=330)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 AD,C,DR,ar,adr,bd
      DIMENSION AD(IA ),BD(IB+1)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /LEG/ P(2*nterms),PD(2*nterms),Q(2*nterms),QD(2*nterms)
      COMMON /FACT/ FACT(300)
      COMMON /eps3/ eps3
c?       eps4 = 1d-20
       eps4 = 1d-40
      ier=0
      M2=2*M
      aM2 = 2d0 * M
      nm1=n-m
      nm2=n+m
      IF(MOD((N-M),2)) 1,2,1
    1  CONTINUE
      JK=1
      GO TO 36
    2  CONTINUE
      JK=0
   36 CONTINUE
      IF(MOD(N,2)) 3,4,3
    3  CONTINUE
      Jn=1
      GO TO 5
    4  CONTINUE
      Jn=0
    5 CONTINUE
      IF(MOD(m,2)) 33,34,33
   33  CONTINUE
      Jm=1
      GO TO 35
   34  CONTINUE
      Jm=0
   35 CONTINUE
      SIG=S
      DO 7 I=1,IA
      J=2*(I-1)+JK
       rf2=AD(I)*fact2(m2+j,j)
      SIG=SIG+rf2
      if(cdabs(rf2/sig).lt.eps4) go to 8
    7  CONTINUE
    8  CONTINUE
      if((nm1-jk)/2+1.gt.170) stop 1702
      ff3 =  fact3(nm2+jk,(nm2+jk)/2,1d0/fact((nm1-jk)/2+1))
      RK1=(am2+1d0+jk*2d0)/2d0**nm2 * ff3 /
     *    AD(1)/fact(m+1)*sig/C**(M+JK)
      if(k.eq.1) rk1=rk1/(-s1)**(M+JK)*(-s1)**jn
      R=S
      DR=S
      IF(I12.EQ.0) GO TO 70
      RK2=2d0**nm1/(am2-1d0)*fact(m2+1)/
     *    fact(m+1) / ff3
     *    *bd(m+1)*sig/C**(M-1-JK)
      if(jk.eq.1.and.k.eq.0) rk2=-rk2/(am2-3d0)
      if(k.eq.0) go to 339
      if(jk.eq.1.and.jm.eq.0) rk2=rk2/(am2-3d0)
      if(jk.eq.1.and.jm.eq.1) rk2=-rk2/(am2-3d0)*(-s1)
      if(jk.eq.0.and.jm.eq.0) rk2=rk2*s1
      if(k.eq.1) rk2=rk2/(-s1)**(M-1-JK)
  339 continue
      AR=s
      ADR=s
      DO 12 I=1,M
        J2=2*(M-I)+1+JK
      ADR=ADR+BD(I+1)*QD(J2)
   12 AR=AR+BD(I+1)*Q(J2)           
   70 CONTINUE
      DO 13 I=1,IA
        J1=2*(I-1)+1+JK
        J2=2*(M+I)-1+JK
      RS=AD(I)*P(J1)
      RD=AD(I)*PD(J1)
      R=R+RS
      DR=DR+RD
      IF(I12.EQ.0) GO TO 17
      rBR=AD(I)*Q(J2)
      rABR=AD(I)*QD(J2)
      AR=AR+rBR
      ADR=ADR+rABR
c      IF(cDABS(rBR/ar).LT.eps4.AND.cDABS(rABR/adr).LT.eps3) GO TO 14
   17 CONTINUE
      IF(CDABS(RS/r).LT.eps4.AND.CDABS(RD/dr).LT.eps4) GO TO 14
   13 CONTINUE
   14 CONTINUE
      IF(I12.EQ.0) GO TO 71
      MR=M+1
      DO 15 I=MR,IB
        J2=2*(I-M)-JK
      if (k.eq.1.and.mod(m,2).eq.0) go to 116
      rBR=BD(I+1)*P(J2)
      rABR=BD(I+1)*PD(J2)
      go to 115
  116 CONTINUE
      rBR=-BD(I+1)*P(J2)
      rABR=-BD(I+1)*PD(J2)
  115 CONTINUE
      AR=AR+rBR
      ADR=ADR+rABR
c      IF(cDABS(rBR/ar).LT.eps4.AND.cDABS(rABR/adr).LT.eps4) GO TO 16
   15 CONTINUE
   16 AR=AR/RK2
      ADR=ADR/RK2
   71 R=R/RK1
      DR=DR/RK1
  900 RETURN
      END

c************************************************
C CESSEL0
c
      SUBROUTINE CESSEL0(A,NUM,BESJ,BESY,iiy)
      IMPLICIT COMPLEX*16 (A-H,O-Q,T-Z)
      DIMENSION BESJ(NUM+1),BESY(NUM+1),iiy(num+1)
      factor=1d250
   11 continue
          do i =1,num+1
            iiy(i) = 0
c+&
            BESJ(i)=(0D0,0D0)
            BESY(i)=(0D0,0D0)
          end do
c      BESJ(NUM+1)=(0D0,0D0)
      BESJ(NUM)=(1D-300,1D-300)
      N=2*NUM+1
      NUM1=NUM-1
      DO I=1,NUM1
       N=N-2
       I1=NUM-I
       BESj(I1) = n * a * BESj(I1+1) - BESj(I1+2)
        if(cdabs(besJ(i1)).gt.1d300) then
           num = i - 1
           go to 11
        end if
      end do
      B1=A*BESJ(1)-BESJ(2)
      B2=-A*B1-BESJ(1)
      N=2*(NUM/2)
      B=1.2533141373155002D0*CDSQRT(A)
      C=B*BESJ(1)
      DO 12 I=3,N,2
      B=B*(I-0.5D0)*(I-2.0D0)/(I-2.5D0)/(I-1.0D0)
   12 C=C+B*BESJ(I)
      C=1.0D0/C
      DO 13 I=1,NUM
      BESJ(I)=C*BESJ(I)
   13 continue
c Y
      BESY(1)=-C*B1
      BESY(2)=C*B2
      DO 14 I=3,NUM
      I2=I-2
         if(iiy(i-1).eq.iiy(i2)) then
          BESY(I) = (2.0D0*I2+1.0D0) * a * BESY(I-1) - BESY(I2)
         else
          BESY(I) = (2.0D0*I2+1.0D0) * a * BESY(I-1) - BESY(I2)/factor
         end if
         iiy(i) = iiy(i-1)
         if(cdabs(besY(i)).gt.1d300) then
           besY(i) = besY(i) / factor
           iiy(i) = iiy(i-1) + 1
         end if
   14 continue
      RETURN
      END

c************************************************
C Cmfunq3
C
c      SUBROUTINE Cmfunq2(II,M,C1,C2,c3,ksi1,ksi2,EPS,NE,*)
      SUBROUTINE Cmfunq3(II,M,C1,C2,c3,ksi1,ksi2,NE,*)
      parameter (nterms=330, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      REAL*8 ksi1,ksi2
      COMPLEX*16 bdc1
      COMPLEX*16 C2, c3, bdc2, BCOF, VV2, ccof, vv3, p1, p3,
     *           DEL21, del23, GAM21, GAM22, akap22
      DIMENSION
     *          RLC1(nterms), RLC2(nterms), RLC3(nterms),
     *          ADC1(4*nterms), bDC1(4*nterms),
     *          rDC2(4*nterms), bDC2(4*nterms),

     *          ACOF(nterms,4*nterms), VV1(nterms),
     *          BCOF(nterms,4*nterms), VV2(nterms),
     *          cCOF(nterms,4*nterms), VV3(nterms),

     *          BD(nalphas,2*nterms), B0(nterms),
     *          Ialpha(nalphas),Jalpha(nalphas),
     *          ASIN(nalphas), ACOS(nalphas)

  202 FORMAT(1X,'Cmfunq3 NN>nterms',5X,'NN=',I5)
  210 FORMAT(1X,'L=',I4,5X,'IER1,2,3,4,5=',5I5)
  212 FORMAT(1X,'I,L=',2I5,5X,'IER=',I5)
      COMMON /Kc1/  ACOF, VV1
      COMMON /Kc2/  BCOF, VV2
      COMMON /Kc3/  cCOF, VV3
      COMMON /Kfun/ ar1_11(2*nterms), r3_11(2*nterms),
     *              aar1_11(nterms), rr3_11(nterms),
     *              r1_21(2*nterms), r3_21(2*nterms),
     *              rr1_21(nterms), rr3_21(nterms),

     *              r1_22(2*nterms), r3_22(2*nterms),
     *              rr1_22(nterms), rr3_22(nterms),
     *              r1_32(2*nterms), rr1_32(nterms),
     *              p1(nterms), p3(nterms)

      COMMON /K1/ S, S1, AKSI, AK, kK, NK, nal
      COMMON /K2/  Ialpha, Jalpha
      COMMON /K10/ B0
      COMMON /K11/ BD
      COMMON /FF/  ASIN, ACOS
      COMMON /INT/ DEL21(nterms,nterms), del23(nterms,nterms),
     *             GAM21(nterms,nterms), GAM22(nterms,nterms),
     *             akap22(nterms,nterms), sig22(nterms,nterms)
      COMMON /EPS1/ EPS1
      COMMON /EPS33/ EPS33
      COMMON /Print/ iprint
      k = kk
      W1=1D0/(C1*(ksi1**2-1D0+2*K))
      W11 = 1D0/(C2*(ksi1**2-1D0+2*K))
      W111 = 1D0/(C2*(ksi2**2-1D0+2*K))
c * DRAINE+!
c      ncc = c1
c ** DRAINE+N! ???
      ncc0 = c1
      ncc2 = cdabs(c2)
      ncc3 = cdabs(c3)
       ncc = 2*ne/3
      IF(NE-nterms)  40,40,41
   41 WRITE(7,202) NE
      WRITE(*,202) NE
      RETURN 1
   40 CONTINUE
      RC1=DCMPLX(C1,0D0)

c++++++++
c     Lambda ---> c1
c++++++++
          if(k.eq.0) kd = 1
          if(k.eq.1) kd = -1
            do jj = 1,ne
            CALL SEGV1(M,jj,rC1,KD,CV,rlc1)
c            write(7,*) jj, rlc1(jj)
c            write(*,*) jj, rlc1(jj)
c            pause
            end do
c-------

c++++++++
c     Lambda ---> c2
c++++++++
      CALL CDLAM_m(RLC2,M,Ne,C2,IE)
c            write(7,*) (rlc2(jj), jj=1,ne)
         if(ie.ne.0) return 1
c-------

c++++++++
c     Lambda ---> c3
c++++++++
c      call lambda(K,M,ne,C3,EPS,rlc3)
      CALL CDLAM_m(RLC3,M,Ne,C3,IE)
c-------
       if(ii.eq.1) CALL funlegnn1(M,2*nk+2)

      inum1=nk
      if(nk.lt.inum1) inum1=nk+10

c       if(ii.eq.1) CALL funlegnn1(M,2*nk+2)
       if(ii.eq.1) CALL funlegnn1(M,inum1)

      if(k.eq.0) then
      ifun1 = 44
c ** DRAINE+N! ???
cccccccccccccccccccccccccccc      IF(ksi0.GE.1.07D0) IFUN1 = 11
c      ifun1 = 11
c--back
       IF(ksi1.GE.1.5D0) IFUN1 = 22
c+&      IF(ksi0.GE.1.5D0) IFUN1 = 22
cccccccccccccccccc      IF(ksi1.GE.1.3D0) IFUN1 = 22
      end if

      if(k.eq.1) then
      ifun1 = 33
c      ifun1 = 11
c      if(l.gt.ncc.and.ksi0.gt.0.5d0) ifun1 = 11
      IF(ksi1.GE.1.5D0) IFUN1 = 22
c      IF(ksi0.GE.1.3D0) IFUN1 = 22
c+&      IF(ksi0.GE.0.6D0) IFUN1 = 22
      end if
c      IF(ifun1.ne.22.or.ifun1.ne.44) CALL funlegnn(m,ksi1,inum1)

c         print *, c1, c2, c3
c         pause

c
c-------------- functions: c1,c2/ksi1
c
      DO 1 L=1,NE
      LM=L+M-1
      L1=L+NE
      ier1 = 0
      ier2 = 0
      ier3 = 0
      ier4 = 0
      ier5 = 0
      ier0 = 0
      ier = 0
c----------------------------------------------------
c  Functions for `c1'
c----------------------------------------------------
      RL=RLC1(L)

      CALL cdcof4a(RL,RDC2,NK,BDC1,NK,2,rvn,1,M,LM,RC1,IER1)
      if(cdabs(bdc1(m+1)).le.eps33*1d20) eps33 = eps33 / 1d20

c *****  functions ************************************
c        ifun = 11: LEG/LEG
c        ifun = 22: BES/BES
c        ifun = 33: DIFF/DIFF
c        ifun = 44: JAF/JAF

c      if(k.eq.1.and.l.gt.ncc.and.ksi0.gt.0.5d0) ifun1 = 11
c * DRAINE+!
       if(k.eq.1.and.l.gt.ncc.and.ksi1.gt.0.25d0.
     *    and.ksi1.lt.0.8d0) ifun1 = 11
c            if(k.eq.1.and.l.gt.ncc.and.ksi0.gt.0.25d0.
c     *    and.ksi0.lt.1.1d0) ifun1 = 11
c      if(k.eq.1.and.l.gt.ncc) ifun1 = 11
c+&
c       print *, ifun1
c * DRAINE+!
c       if(ifun1.eq.22.or.ifun1.eq.33.or.ifun1.eq.44) go to 553
c ** DRAINE+!
c       if(ifun1.eq.22.or.ifun1.eq.44) go to 553
c         CALL funlegnn(m,ksi0,inum1)
c  553  continue
c ****************************************************
c         print *, l, ifun1
c         pause

      if(ifun1.eq.11) go to 11
      if(ifun1.eq.22) go to 22
      if(ifun1.eq.33) go to 33
      if(ifun1.eq.44) go to 44

c----------------------------------------------------
c    LEG/LEG    ifun1 = 11
c----------------------------------------------------
   11  CONTINUE
c      ifun1 = 11
c ** DRAINE+!
         CALL funlegnn(m,ksi1,inum1)
      CALL CDRF12cc(R1,R2,r3,r4,1,M,LM,RC1,RDC2,NK,
     *            BDC1,NK,IER2)
* -----
      W = r1 * r4 - r2 * r3
      W2a = dabs(W / W1 - 1d0)
      if(W2a.le.eps1) go to 5
      if(W2a.gt.eps1.and.iprint.ne.0) then
        write(*,99999) ifun1,l,W,W1,W2a
       end if
* -----
      ar1_11(L)  = R1
      ar1_11(L1) = R2
      r3_11(l) = r1 + s1 * r3
      r3_11(l1) = r2 + s1 * r4
* -----
c Draine!
c      if(W2a.gt.W0.and.k.eq.0) go to 13
      if(W2a.gt.W0.and.k.eq.0.and.ifun1.eq.44) go to 13
c         print *, l, ifun1
c         pause
c      if(W2a.lt.W0) then
      go to 5
c      end if

c99999 format(1x,'ifun=',i3,2x,'L=',i3,2x,'W=',1pd13.6,2x,
c     *       'W1=',d13.6,2x,'W2=',d13.6)
99999 format(1x,'c1/ksi1  ifun=',i3,1x,'L=',i3,1x,'W=',1pd9.2,1x,
     *       'W1=',d9.2,1x,'W2=',d9.2)

c----------------------------------------------------
c    LEG/DIFF (JAF/DIFF)  ifun1 = 13
c----------------------------------------------------
   13  CONTINUE
      IF(K.EQ.0) CALL cDRSF20(r3,r4,ksi1,R1,R2,m,rc1,rl,IER4)
      W = r1 * r4 - r2 * r3
      W2b = dabs(W / W1 - 1d0)

      ifu = 13
      if(ifun1.eq.22) ifu = 23
      if(ifun1.eq.44) ifu = 43
      if(W2b.gt.eps1.and.iprint.ne.0)  write(*,99999) ifu,l,W,W1,W2b
c*** DRAINE+N!
         IF(W2b.lt.w2a) then
c        IF(W2b.lt.w0) then
        r3_11(l) = r1 + s1 * r3
        r3_11(l1) = r2 + s1 * r4
        end if
      go to 5

c----------------------------------------------------
c    BES/BES    ifun1 = 22
c----------------------------------------------------
   22  CONTINUE
      CALL CDRB12cc(R1,R2,r3,r4,1,M,LM,RC1,ksi1,RDC2,NK,IER3)
* -----
      ar1_11(L)  = R1
      ar1_11(L1) = R2
      r3_11(l) = r1 + s1 * r3
      r3_11(l1) = r2 + s1 * r4
* -----
      W = r1 * r4 - r2 * r3
      W2a = dabs(W - W1)
      if(W2a.le.eps1) go to 5

c * DRAINE+!
      if(W2a.gt.eps1) then
       if(iprint.ne.0) write(*,99999) ifun1,l,W,W1,W2a
      end if

c ** DRAINE+N!
c      if(k.eq.0) go to 13
c      if(k.eq.0.and.l.gt.ncc) go to 11
c      if(k.eq.0.and.l.gt.ncc0) go to 11
c--NEW
      if(k.eq.0.and.l.gt.ncc0.and.c1.gt.20d0) go to 11
c      if(k.eq.0) go to 11
c ** DRAINE+!
c      go to 5

c----------------------------------------------------
c    DIFF/DIFF   ifun1 = 33
c----------------------------------------------------
   33  CONTINUE

c       DO I = 1,NK
c       bDC1(I) = bDC2(I)
c       end do

c       print *, bdc1(m+1)
c?
       do i = 1, m+1
       bdc2(i) = bdc1(i)
       end do

      CALL cDRSF212(R1,R2,r3,r4,2,m,LM,rc1,ksi1,rl,
     *             RDC2,NK,BDC2,m+1,IER4)
c?     *             RDC2,NK,BDC1,m+1,IER4)
* -----
      ar1_11(L)  = R1
      ar1_11(L1) = R2
      r3_11(l) = r1 + s1 * r3
      r3_11(l1) = r2 + s1 * r4
* -----
c         print *, l, ifun1
c         pause

      W = r1 * r4 - r2 * r3
      W2a = dabs(W / W1 -1d0)
c ** DRAINE+!
      if(W2a.gt.eps1) then
       if(iprint.ne.0) write(*,99999) ifun1,l,W,W1,W2a
c DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
c      if(l.gt.ncc0) go to 11
c--NEW
      if(l.gt.ncc0.and.c1.gt.20d0) go to 11
cccccccccc        go to 11
c      go to 21
c      go to 31
      end if
      go to 5

c----------------------------------------------------
c    JAF/JAF    ifun1 = 44
c----------------------------------------------------
 44   CONTINUE
c........
        CALL CDRG1cn(r1,r2,raR3,rAr4,1,M,lm,rC1,ksi1,
     *              RDC2,NK,bdc1,nk,rl,IER0)
c........
* -----
      ar1_11(L)  = R1
      ar1_11(L1) = R2
      r3_11(l) = r1 + s1 * rar3
      r3_11(l1) = r2 + s1 * rar4
* -----
      W = r1 * rar4 - r2 * rar3
      W0 = dabs(W / W1 - 1d0)
      ifun1 = 44

c ** DRAINE+N!
      if(W0.gt.eps1) then
      if(iprint.ne.0) write(*,99999) ifun1,l,W,W1,W0
c      if(l.gt.ncc0) go to 11
c--NEW
      if(l.gt.ncc0.and.c1.gt.20d0) go to 11
      end if
c ---------------------final results
    5   continue
c ---------------------final results

c      IF((IER1+IER2+IER3+IER4).NE.0) WRITE(7,210) L,
c     *    IER1,IER2,IER3,IER4

       DO I = 1,NK
       ADC1(I) = RDC2(I)
       ACOF(L,i) = ADC1(i)
       end do

      VC1L2 = RVN
      VV1(L) = DSQRT(VC1L2)

      IF(II.EQ.0) then
         B0(L) = COFF(L,NK,ADC1)
         GO TO 88
      end if

      DO 8 I = 1, NAL
      IF(Jalpha(I).EQ.1) GO TO 8
      CALL DASWF1a(BD(I,L),BD(I,L1),1,M,LM,ACOS(I),i,ADC1,NK,IER)
      IF(IER.NE.0) WRITE(7,212) I,L,IER
      BD(I,L) = BD(I,L) / VC1L2
    8 CONTINUE

   88 CONTINUE

      ar1_11(L) = ar1_11(L) * VV1(L)
      ar1_11(L1) = ar1_11(L1) * VV1(L)
      r3_11(L) = r3_11(L) * VV1(L)
      r3_11(L1) = r3_11(L1) * VV1(L)
      aar1_11(L) = ar1_11(L1) / ar1_11(L)
      rr3_11(L) = r3_11(L1) / r3_11(L)

c----------------------------------------------------
c  Functions for `c2'/ `ksi1'
c----------------------------------------------------
      RP=RLC2(L)
      CALL cdcof4a(RP,RDC2,NK,BDC2,NK,2,RVN2,1,M,LM,C2,IER5)

      if(ifun1.eq.11) go to 111
      if(ifun1.eq.22) go to 222
      if(ifun1.eq.33) go to 333
      if(ifun1.eq.44) go to 444
c----------------------------------------------------
c    LEG/LEG    ifun1 = 11
c----------------------------------------------------
  111  CONTINUE
      CALL CDRF12cc(r1,r2,r3,r4,1,M,LM,C2,RDC2,NK,
     *             BDC2,NK,IER1)
* -----
      r1_21(L)  = R1
      r1_21(L1) = R2
      r3_21(l) = r1 + s1 * r3
      r3_21(l1) = r2 + s1 * r4
* -----
      W = r1 * r4 - r2 * r3
      W0 = dabs(W / W11 -1d0)
      if(W0.gt.eps1.and.iprint.ne.0) write(*,98) ifun1,l,W,W11,W0
ccccccccccccccc      if(W0.gt.eps1.and.k.eq.0) go to 113
      GO TO 7
c----------------------------------------------------
c    LEG/DIFF (JAF/DIFF)  ifun1 = 13
c----------------------------------------------------
  113  CONTINUE
      IF(K.EQ.0) CALL cDRSF20(r3,r4,ksi1,R1,R2,m,c2,rp,IER2)
      W = r1 * r4 - r2 * r3
      W2b = dabs(W / W11 - 1d0)

      ifu = 13
      if(ifun1.eq.22) ifu = 23
      if(ifun1.eq.22) ifu = 43
      if(W2b.gt.eps1.and.iprint.ne.0) write(*,98) ifu,l,W,W11,W2b

        IF(W2b.lt.w0) then
        r3_21(l) = r1 + s1 * r3
        r3_21(l1) = r2 + s1 * r4
        end if
      go to 7
c----------------------------------------------------
c    BES/BES    ifun1 = 22
c----------------------------------------------------
  222  CONTINUE
      CALL CDRB12cc(r1,r2,r3,r4,1,M,LM,C2,ksi1,RDC2,NK,IER3)
* -----
      r1_21(L)  = R1
      r1_21(L1) = R2
      r3_21(l) = r1 + s1 * r3
      r3_21(l1) = r2 + s1 * r4
* -----
      W = r1 * r4 - r2 * r3
      W0 = dabs(W / W11 -1d0)
      if(W0.le.eps1) go to 7
      if(W0.gt.eps1) then
        if(iprint.ne.0) write(*,98) ifun1,l,W,W11,W0
c        if(k.eq.0) go to 113
c        go to 111
c--New
      if(l.gt.ncc2.and.cdabs(c2).gt.20d0) go to 111
      end if
      go to 7
c----------------------------------------------------
c    DIFF/DIFF   ifun1 = 33
c----------------------------------------------------
  333  CONTINUE
      CALL cdRSF212(r1,r2,r3,r4,2,m,LM,
     *              c2,ksi1,rp,RDC2,NK,BDC2,m+1,IER4)
* -----
      r1_21(L)  = R1
      r1_21(L1) = R2
      r3_21(l) = r1 + s1 * r3
      r3_21(l1) = r2 + s1 * r4
* -----
      W = r1 * r4 - r2 * r3
      W2a = dabs(W / W11 -1d0)
      if(W2a.gt.eps1) then
       if(iprint.ne.0) write(*,98) ifun1,l,W,W11,W2a
c DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
c--New
      if(l.gt.ncc2.and.cdabs(c2).gt.20d0) go to 111
c      if(l.gt.ncc2) go to 111
c        go to 111
      end if
      GO TO 7
c----------------------------------------------------
c    JAF/JAF    ifun1 = 44
c----------------------------------------------------
 444  CONTINUE
c........
        CALL CDRG1cn(r1,r2,r3,r4,1,M,lm,C2,ksi1,
     *              RDC2,NK,bdc2,nk,rp,IER0)
c........
      r1_21(L)  = R1
      r1_21(L1) = R2
      r3_21(l) = r1 + s1 * r3
      r3_21(l1) = r2 + s1 * r4
* -----
      W = r1 * r4 - r2 * r3
      W0 = dabs(W / W11 - 1d0)
      if(W0.gt.eps1.and.iprint.ne.0) write(*,98) ifun1,l,W,W11,W0
      if(l.gt.ncc2.and.cdabs(c2).gt.20d0) go to 111
cccccccccccccccccc      if(W0.gt.eps1)  go to 113

c ---------------------final results
    7 CONTINUE
c=====
   98 format(1x,'c2/ksi1 L=',2i4,1x,3(2d10.3,2x))
c=====
c      IF((IER0+IER1+IER2+IER3+IER4).NE.0) WRITE(7,210) L,
c     *    IER1,IER2,IER3,IER4,IER5
c.......................
           if(cdabs(rvn2).lt.1d20) vv2(l) = cdsqrt(rvn2)
           if(cdabs(rvn2).gt.1d20.and.cdabs(rvn2).lt.1d80)
     *        vv2(l) = cdsqrt(rvn2/1d40)*1d20
           if(cdabs(rvn2).gt.1d80) vv2(l) = cdsqrt(rvn2/1d80)*1d40
c.......................

      r1_21(L ) = r1_21(L ) * VV2(L)
      r1_21(L1) = r1_21(L1) * VV2(L)
      r3_21(L ) = r3_21(L ) * VV2(L)
      r3_21(L1) = r3_21(L1) * VV2(L)
      rr1_21(L) = r1_21(L1) / r1_21(L)
      rr3_21(L) = r3_21(L1) / r3_21(L)

      DO N = 1, NK
      BCOF(L,N)=rDC2(N)
      cCOF(L,N)=bDC2(N)
      end do

    1 CONTINUE

        CALL DELcc(21,M,NE,DEL21)
        CALL GAMKAP2(M,NE,GAM22,akap22)
        CALL GAMc2c1(M,NE,GAM21)
        CALL sigma22(M,NE,sig22)

c
c-------------- functions: c2,c3/ksi2
c
* ----------------------------------
      if(k.eq.0) then
      ifun3 = 44
c      ifun3 = 22
      IF(ksi2.GE.1.5D0) IFUN3 = 22
c     IF(ksi2.GE.1.3D0) IFUN3 = 22
c      IF(ksi2.GE.2D0) IFUN3 = 22
c      IF(ksi2.lt.1.2D0) IFUN3 = 44
      end if

      if(k.eq.1) then
c      ifun3 = 11
      ifun3 = 33
      IF(ksi2.GE.1.5D0) IFUN3 = 22
      end if

c      print *, 'ifun3', ifun3
c      pause
* ----------------------------------

      DO 10 L=1,NE
      LM=L+M-1
      L1=L+NE
      ier1 = 0
      ier2 = 0
      ier3 = 0
      ier4 = 0
      ier5 = 0
      ier0 = 0
      ier = 0
      RP=RLC2(L)
c      CALL cdcof4a(RP,RDC2,NK,BDC2,NK,2,RVN2,1,M,LM,C2,IER5)
c----------------------------------------------------
c  Functions for `c2' / `ksi2'
c----------------------------------------------------
      DO N = 1, NK
      rDC2(N) = BCOF(L,N)
      bDC2(N) = cCOF(L,N)
      end do

      if(ifun3.eq.11) go to 1111
      if(ifun3.eq.22) go to 2222
      if(ifun3.eq.33) go to 3333
      if(ifun3.eq.44) go to 4444
c----------------------------------------------------
c    LEG/LEG    ifun3 = 11
c----------------------------------------------------
 1111  CONTINUE
c ** DRAINE+!
         CALL funlegnn(m,ksi2,inum1)
      CALL CDRF12cc(r1,r2,r3,r4,1,M,LM,C2,RDC2,NK,
     *             BDC2,NK,IER4)
* -----
      r1_22(L)  = R1
      r1_22(L1) = R2
      r3_22(l) = r1 + s1 * r3
      r3_22(l1) = r2 + s1 * r4
* -----
      W = r1 * r4 - r2 * r3
      W0 = dabs(W / W111 -1d0)
      if(W0.gt.eps1.and.iprint.ne.0) then
        write(*,97) ifun3,l,W,W111,W0
cccccccccccccccccccc        if(k.eq.0) go to 1113
      end if

      GO TO 77
c----------------------------------------------------
c    LEG/DIFF (JAF/DIFF)  ifun3 = 13
c----------------------------------------------------
 1113  CONTINUE
      IF(K.EQ.0) CALL cDRSF20(r3,r4,ksi2,R1,R2,m,c2,rp,IER4)
      W = r1 * r4 - r2 * r3
      W2b = dabs(W / w111 - 1d0)
      ifu = 13
      if(ifun3.eq.22) ifu = 23
      if(ifun3.eq.22) ifu = 43
      if(W2b.gt.eps1.and.iprint.ne.0) write(*,97) ifu,l,W,w111,W2b

        IF(W2b.lt.w0) then
        r3_22(l) = r1 + s1 * r3
        r3_22(l1) = r2 + s1 * r4
        end if
      go to 77
c----------------------------------------------------
c    BES/BES    ifun3 = 22
c----------------------------------------------------
 2222  CONTINUE
      CALL CDRB12cc(r1,r2,r3,r4,1,M,LM,C2,ksi2,RDC2,NK,IER3)
* -----
      r1_22(L)  = R1
      r1_22(L1) = R2
      r3_22(l) = r1 + s1 * r3
      r3_22(l1) = r2 + s1 * r4
* -----
      W = r1 * r4 - r2 * r3
      W0 = dabs(W / w111 -1d0)
      if(W0.gt.eps1) then
       if(iprint.ne.0) write(*,97) ifun3,l,W,w111,W0
c        if(k.eq.0) go to 1113
c--New
      if(l.gt.ncc2.and.cdabs(c2).gt.20d0) go to 1111
c        if(l.gt.ncc2) go to 1111
c        go to 1111
      end if
      go to 77

c----------------------------------------------------
c    DIFF/DIFF   ifun3 = 33
c----------------------------------------------------
 3333  CONTINUE
      CALL cDRSF212(r1,r2,r3,r4,2,m,LM,
     *               c2,ksi2,rp,RDC2,NK,BDC2,m+1,IER4)
* -----
      r1_22(L)  = R1
      r1_22(L1) = R2
      r3_22(l) = r1 + s1 * r3
      r3_22(l1) = r2 + s1 * r4
* -----
      W = r1 * r4 - r2 * r3
      W2a = dabs(W / w111 -1d0)
      if(W2a.gt.eps1) then
        if(iprint.ne.0) write(*,97) ifun3,l,W,w111,W2a
c--New
      if(l.gt.ncc2.and.cdabs(c2).gt.20d0) go to 1111
c        if(l.gt.ncc2) go to 1111
c        go to 1111
      end if
      go to 77

c----------------------------------------------------
c    JAF/JAF    ifun3 = 44
c----------------------------------------------------
4444  CONTINUE
c........
        CALL CDRG1cn(r1,r2,r3,r4,1,M,lm,C2,ksi2,
     *              RDC2,NK,bdc2,nk,rp,IER0)
c........
      r1_22(L)  = R1
      r1_22(L1) = R2
      r3_22(l) = r1 + s1 * r3
      r3_22(l1) = r2 + s1 * r4
* -----
      W = r1 * r4 - r2 * r3
      W0 = dabs(W / w111 - 1d0)
      if(W0.gt.eps1) then
          if(iprint.ne.0) write(*,97) ifun3,l,W,W111,W0
      if(l.gt.ncc2.and.cdabs(c2).gt.20d0) go to 1111
cccccccccccccccccccccccc        go to 1113
      end if
c ---------------------final results
   77 CONTINUE
c=====
   97 format(1x,'c2/ksi2  L=',2i4,2x,3(2d10.3,2x))
c=====
      r1_22(L ) = r1_22(L ) * VV2(L)
      r1_22(L1) = r1_22(L1) * VV2(L)
      r3_22(L ) = r3_22(L ) * VV2(L)
      r3_22(L1) = r3_22(L1) * VV2(L)
      rr1_22(L) = r1_22(L1) / r1_22(L)
      rr3_22(L) = r3_22(L1) / r3_22(L)

      p1(L) = r1_22(L) / r1_21(L)
      p3(L) = r3_22(L) / r3_21(L)

   10 CONTINUE
*--------------------------------------

      DO 110 L=1,NE
      LM=L+M-1
      L1=L+NE
      ier1 = 0
      ier2 = 0
      ier3 = 0
      ier4 = 0
      ier5 = 0
      ier0 = 0
      ier = 0
c----------------------------------------------------
c  Functions for `c3' / `ksi2'
c----------------------------------------------------
      Rq=RLC3(L)
      CALL cdcof4a(Rq,RDC2,NK,BDC2,NK,2,RVN3,1,M,LM,C3,IER1)

      if(IFUN3.eq.11)
     *   CALL CDRF12cc(r1_32(l),r1_32(l1),r3,r4,0,M,LM,
     *                             C3,RDC2,NK,BDC2,NK,IER2)
      if(IFUN3.eq.22) CALL CDRB12cc(r1_32(l),r1_32(l1),r3,r4,0,M,LM,
     *                             C3,ksi2,RDC2,NK,IER3)
      if(ifun3.eq.33) CALL cDRSF212(r1_32(l),r1_32(l1),r3,r4,0,m,LM,
     *                             c3,ksi2,rq,RDC2,NK,BDC2,m+1,IER4)
      if(ifun3.eq.44)
     *   CALL CDRG1cn(r1_32(l),r1_32(l1),r3,r4,0,M,lm,
     *                            C3,ksi2,RDC2,NK,bdc2,nk,rq,IER0)
c       print *,c3, ksi2
c       print *,l, r1_32(l)
c       print *,l, r1_32(l1)
c       pause
c.......................
           if(cdabs(rvn3).lt.1d20) vv3(l) = cdsqrt(rvn3)
           if(cdabs(rvn3).gt.1d20.and.cdabs(rvn3).lt.1d80)
     *        vv3(l) = cdsqrt(rvn3/1d40)*1d20
           if(cdabs(rvn3).gt.1d80) vv3(l) = cdsqrt(rvn3/1d80)*1d40
c.......................
      r1_32(L ) = r1_32(L ) * VV3(L)
      r1_32(L1) = r1_32(L1) * VV3(L)
      rr1_32(L) = r1_32(L1) / r1_32(L)

      DO N = 1, NK
      cCOF(L,N) = rDC2(N)
      end do

  110 CONTINUE

      CALL DELcc(23,M,NE,DEL23)

      RETURN
      END

c************************************************
c cmlm
c
      SUBROUTINE CMLM(N,A1,A2,A3)
      parameter (nalphas=19,nterms=330)
      IMPLICIT COMPLEX*16(A-H,O-Z)
      DIMENSION A1(nterms,nterms),A2(nterms,nterms),A3(nterms,nterms),
     * A4(nterms,nterms)

* multiplication of two complex matrices

      DO 1 I=1,N
      DO 1 J=1,N
      B=(0D0,0D0)
      DO 2 K=1,N
    2 B=B+A1(I,K)*A2(K,J)
      A4(I,J)=B
    1 CONTINUE
      DO 3 I=1,N
      DO 3 J=1,N
    3 A3(I,J)=A4(I,J)
      RETURN
      END

c************************************************
c  daswf1a
C
      subroutine DASWF1a(ar,adr,is,M,n,x,icos,AD,ia,IER)
      parameter (nterms=330, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16 (R-S)
      DIMENSION AD(IA),P(nalphas,2*nterms),PD(nalphas,2*nterms)
      COMMON /LEG1/ P,PD
      ier=0
      ar=0d0
      aDR=0d0
      eps4 = 1d-40
c      eps4 = 1d-20
      if(dabs(x-1d0).lt.1d-6) return
      if(dabs(x).lt.1d-6) go to 14
      IF(MOD((N-M),2)) 1,2,1
    1 JK=1
      GO TO 5
    2 JK=0
    5 CONTINUE
      DO 13 I=1,IA
      J1=2*(I-1)+1+JK
      aRS=AD(I)*P(icos,J1)
      aRD=AD(I)*PD(icos,J1)
      aR=aR+aRS
      aDR=aDR+aRD
      IF(DABS(aRS/ar).LT.eps4.AND.DABS(aRD/adr).LT.eps4) GO TO 900
   13 CONTINUE
      go to 900
c                          x = 0
   14 CONTINUE
      CALL DLEGF00(ar,adr,n,M,is,IER1)
  900 continue
      RETURN
      END

c************************************************
c DELcc
c
      SUBROUTINE DELcc(ic,M,NN,DEL)
      parameter (nalphas=19,nterms=330)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 BCOF,B,A,VV2,DEL,C,E1,G1,D,ccof, VV3
      COMPLEX*16 s, s1
      DIMENSION A(4*nterms),B(4*nterms),
     *          cCOF(nterms,4*nterms), VV3(nterms),
     *          BCOF(nterms,4*nterms), VV2(nterms),
     *          ACOF(nterms,4*nterms), VV1(nterms),
     *          DEL(nterms,nterms)
      COMMON /EPS3/ EPS3
      COMMON /K1/ S, S1, AKSI, AK, kK, NK, nal
      COMMON /Kc1/  ACOF, VV1
      COMMON /Kc2/  BCOF, VV2
      COMMON /Kc3/  cCOF, VV3
c      eps4 = 1d-20
       M2=2*M
      DO 10 N=1,NN
        NM=N+M-1
         DO I = 1, NK
          B(I) = BCOF(N,I)
         end do
        C = VV2(N)
      DO 11 L=1,NN
        LM=L+M-1
        IF(MOD(IABS(LM-NM),2)) 9,9,20
    9 CONTINUE

      if(ic.eq.21) then
      DO I = 1, NK
      a(I) = aCOF(l,I)
      end do
      d = VV1(l)
      end if

      if(ic.eq.23) then
      DO I = 1,NK
      a(I) = cCOF(l,I)
      end do
      d = VV3(l)
      end if

      IF(MOD(LM-M,2))1,1,2
    1 I2=1
      K=M2
      E1=A(1)*B(1)*2D0/(M2+1D0)
      GO TO 3
    2 I2=0
      K=M2+1
      E1=A(1)*B(1)*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E1=E1*F
      DO 5 I=2,NK
      J=2*I-1-I2
      JM=J+M
      JM2=2*JM
c      F=1D0
      F=a(i)      
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
c      IF(M2.EQ.0) F=1D0
      G1=B(I)*2D0/(JM2+1D0)*F
      IF(CDABS(G1).LE.1D-140.and.i.gt.nn) GO TO 6
c      IF(CDABS(G1).LE.eps3.and.i.gt.nn) GO TO 6
c      IF(CDABS(G1/e1).LE.eps4.and.i.gt.nn) GO TO 6
      E1=E1+G1
    5 continue
    6 DEL(N,L)=E1/C/D
      GO TO 11
   20 DEL(N,L)=s
   11 CONTINUE
   10 CONTINUE
      RETURN
      END

c************************************************
C DLEGf00
C
      subroutine DLEGF00(res,red,n,M,id,iER)
      IMPLICIT REAL*8(A-H,O-Z)                                        
      COMMON  /FACT/ FACT(300)
      ier=0
      res=0d0
      red=0d0
      if(n.eq.0.and.m.eq.0) res=1d0
      if(n.eq.0.and.m.eq.0)  go to 100
      if(n.le.0.or.m.lt.0) ier=1
      if(ier.ne.0) go to 100
c**************************************   << x1 = 0 >>
      IF(MOD(iabs(N-M),2)) 50,1,50
    1 continue
      k1=(m+n)/2
c+&
      e=1d0/2d0**n
      do 2 kk=1,k1
    2 e=e*(k1+kk)
      k2=(n-m)/2
      if(k2+1.gt.170) stop 1703
c+&      res=(-1)**k2*e/fact(k2+1)/2d0**n
      res=(-1)**k2*e/fact(k2+1)
c           print *,m,n, e, (-1)**k2, 2d0**n
c           print *, fact(k2+1), res
      if(id) 50,100,50
   50 continue
c=========================             <<     Ø‡Æ®ß¢Æ§≠†Ô    >>
      if(id.eq.0) go to 100
c************************************  << x1 = 0 >>
      IF(MOD(iabs(N-M),2)) 71,100,71
   71 continue
      k1=(m+n+1)/2
c+&      e=1d0
      e=1d0/2d0**n
      do 72 kk=1,k1
   72 e=e*(k1+kk)
      k2=(n-m-1)/2
      red=(-1)**k2*e/fact(k2+1)
c+&      red=(-1)**k2*e/fact(k2+1)/2d0**n  
  100 CONTINUE                                                      
      RETURN                                                        
      END      

c************************************************
C DLEGf1
C   
      subroutine DLEGF1(rs,rd,X1,x2,x3,n,M,ic,id,iER)
      IMPLICIT REAL*8(A-H,O-Z) 
      complex*16 s, s1
      COMMON /FACT/ FACT(300)
      COMMON /K1/ S, S1, AKSI, AK, Kkk, NK, nal
      ier=0
      res=0d0
      red=0d0
c      eps4 = 1d-20
      xx2=dabs(x2)
      nm=n-m
      j=mod(n,4)
      jm=mod(m,4)
      if(n.eq.0.and.m.eq.0) res=1d0
      if(n.eq.0.and.m.eq.0)  go to 100
      if(n.le.0.or.m.lt.0) ier=1
      if(id.ne.0.and.dabs(x1-1d0).lt.1d-6.and.ic.eq.0) ier=2
      if(ier.ne.0) go to 100
      if(dabs(x1).gt.1d-6)  go to 20
c**************************************   << x1 = 0 >>
      IF(MOD(iabs(nm),2)) 100,1,100
    1 continue
      k1=(m+n)/2
      e=1d0
      do 2 k=1,k1
    2 e=e*(k1+k)
      if(nm/2+1.gt.170) stop 1703
      res=e/fact(nm/2+1)/2d0**n
      if(id) 50,102,50
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   << x1 /= 0 >>
   20 continue
      if(m.gt.n) go to 100
      if(nm.ge.2.and.x1.gt.1d0.and.ic.eq.0) go to 48
      z=1d0/(2d0**n)
      k1=(nm+2)/2
      d=0d0
      if(ic.ne.0) go to 21
c.....................................   << ic = 0 >>
c      if(nm.gt.169) go to 291
      if(nm.gt.169) stop 170

c       print *, m, n, k1
c       pause

      do 22 k=1,k1
      k3=nm-2*(k-1)
      if(k3+1.gt.170) stop 1704
      h=1d0/fact(k3+1)
      k4=n-k+1

c       print *,  k, k4
c       pause

      do 27 i=1,k4
      if(dabs(h).gt.1d300) print *,  n, h

   27 h=h*(k4+i)
      id = 1
       if(mod(k-1,2).ne.0) id = -1
       dd = id * h / fact(k) * x1**k3
c       print *,  dd, d
      d = d + dd
c+&      if(dabs(dd/d).lt.eps4) go to 23
   22 continue
   23 continue
       res=d*z*dsqrt(dabs(1d0-x1**2)**m)
      if(id) 50,102,50     
c<<<<<<<<<++++++++++++++++++++++++++++   << ic /= 0 >>
   21 continue
c      if(nm.gt.169) go to 29
      if(nm.gt.169) stop 1705
      do 28 k=1,k1
      k3=nm-2*(k-1)
      h=1d0/fact(k3+1)
      k4=n-k+1
c      if(k4.gt.169) stop 17055
      do 33 i=1,k4
   33 h=h*(k4+i)
   28 d=d+h/fact(k)*x1**k3
      res=d*z*dsqrt((1+x1**2)**m)
      if(id) 50,102,50     
   48 continue
      res=((2d0*n-1d0)*x1*x2-(n+m-1d0)*x3)/nm
      if(id) 50,102,50     
   50 continue
c=========================             <<  derivature  >>
c*****************************************   << x1 /= + - 1 >>
      if(ic.ne.0) go to 101
      red=((nm)*res-x1*n*x2)/(x1**2-1d0)
       rs=res
       rd=red
       return
  101 continue
      red=((nm)*res-x1*n*xx2)/(x1**2-ak)
  100 CONTINUE        
      if(jm.eq.2.or.jm.eq.3) go to 103
      rs=res
      if(j.eq.2.or.j.eq.3) rs=-dabs(res)
      go to 104
  103 continue
      rs=-dabs(res)
      if(j.eq.0.or.j.eq.1) rs=-rs
  104 continue
      j=mod(n-1,4)                          
      rd=red
      if(j.eq.2.or.j.eq.3) rd=-dabs(red)
      RETURN      
  102 continue
       rs=res
      if(ic.eq.0.and.jm.eq.3.and.x1.gt.1d0) rs=dabs(res)
      if(ic.ne.0.and.jm.eq.2) rs=-dabs(res)
      if(ic.ne.0.and.jm.eq.3) rs=-dabs(res)
      RETURN                                               
      END

c************************************************
C DLEGf2
C   
      subroutine DLEGF2(res,red,X1,x2,n,M,ic,id,iER) 
      IMPLICIT REAL*8(A-H,O-Z)  
      complex*16 s1, s
c + &
      dimension r(3000),q(3000)
      COMMON /K1/ S, S1, AKSI, AK, Kkk, NK, nal
      COMMON /FACT/ FACT(300)
      COMMON /PI/ PI
      ier=0
      res=0d0
      red=0d0
      if(n+m.lt.0) ier=1
      if(ic.eq.0.and.dabs(x1).le.1d0) ier=2
      if(ic.ne.0.and.dabs(x1).lt.0d0) ier=2
      if(n+m-1.lt.0.and.id.ne.0) ier=3
      if(ier.ne.0) go to 100
      if(dabs(x1).gt.1d-6.or.ic.eq.0)  go to 20
c**************************************   << x1 = 0 & ic /= 0 >>
      if(m.gt.n) go to 3
c                             m  < = n
      k1=n-m
      k2=n+m-1
      d=1d0
      e=1d0
      IF(MOD(k1,2)) 2,1,2
    1 continue
      if(k1.eq.0) go to 6
      do 4 i=1,k2,2
    4 d=d*i
      do 5 i=2,k1,2
    5 e=e*i
    6 res=-pi/2d0*d/e
      if(id) 50,100,50
    2 continue
      do 7 i=2,k2,2
    7 d=d*i
      do 8 i=1,k1,2
    8 e=e*i
      res=d/e
      if(id) 50,100,50
c                                m  >  n 
    3 continue
      IF(MOD(iabs(N-M),2)) 10,9,10
    9 go to 100
   10 continue
      k1=(n-m-1)/2
      k2=(m+n-1)/2
      d=1d0
      do 11 i=1,k1
   11 d=d*(k1+i)
      if(k1+1.gt.170) stop 1706
      res=2d0**n*d*fact(k1+1)
      if(id) 50,100,50
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$   << x1 /= 0 >>
   20 continue
      if(m.eq.0)  go to 61
      if(m.eq.1)  go to 71
c *********************************     m  > 1
      if(ic.eq.0) res=dlog((x1+1d0)/(x1-1d0))/2d0
      if(ic.ne.0.and.x1.lt.0.5d0) res=(datan(x1)-pi/2d0)
      if(ic.ne.0.and.x1.ge.0.5d0) res=-datan(1d0/x1)
      r(1)=-ak/(dsqrt(x1**2-ak))
      r(2)=-2d0*x1*ak/(dsqrt(x1**2-ak))*r(1)
      if(m.eq.2) go to 221
      do 21 i=3,m
   21 r(i)=-2d0*(i-1d0)*x1*ak/(dsqrt(x1**2-ak))*r(i-1)+
     *     (2d0-i)*(i-1d0)*r(i-2)
  221 continue
      q0=r(m)
      if(n.eq.0) res=q0
      if(n.eq.0) go to 26
      j=n
      if(n.lt.0) j=-n
      nmax=30*x1+j+50
c+&      if(nmax.gt.500) nmax=500
c      if(nmax.gt.1500) print *, 'dlegf2', nmax
c      if(nmax.gt.1500) nmax=1500
      if(nmax.gt.3000) print *, 'dlegf2', nmax
      if(nmax.gt.3000) nmax=3000
      r(nmax)=x1-dsqrt(x1**2-ak)
      xer=dabs(x1-1d0/dsqrt(3d0))
      do 22 i=nmax-1,1,-1
   22 r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
      q(1)=r(1)*q0*ak
      if(xer.lt.1d-10.and.m.eq.3) q(1)=3d0
      if(xer.lt.1d-10.and.m.eq.6) q(1)=-144d0
      if(xer.lt.1d-10.and.m.eq.9) q(1) = 4.536d4
      if(xer.lt.1d-10.and.m.eq.12) q(1) = -4.35456d+07
      if(xer.lt.1d-10.and.m.eq.15) q(1) = 9.340529643126767d+10
      if(xer.lt.1d-10.and.m.eq.18) q(1) = -3.76610155207262d+14
      if(xer.lt.1d-10.and.m.eq.21) q(1) = 2.554546682742534d+18
      if(xer.lt.1d-10.and.m.eq.24) q(1) = -2.697601737669382d+022
      if(xer.lt.1d-10.and.m.eq.27) q(1) =  4.188026680200701d+026
      if(xer.lt.1d-10.and.m.eq.30) q(1) = -9.146650294330210d+030
      if(xer.lt.1d-10.and.m.eq.33) q(1) =  2.713536749598336d+035
      if(xer.lt.1d-10.and.m.eq.36) q(1) = -1.062838075962123d+040
      if(xer.lt.1d-10.and.m.eq.39) q(1) =  5.367863692727776d+044
      if(xer.lt.1d-10.and.m.eq.42) stop 42
      if(xer.lt.1d-10.and.mod(m,3).eq.0) q0=0d0
      if(n.eq.1) go to 24
      if(n.lt.0) go to 225
      do 23 i=1,j-1
   23 q(i+1)=r(i+1)*q(i)*ak**(i+1)
       go to 24
c******************************************      n < 0
  225 continue 
      if(n.eq.-1) res=(ak*x1*q0+(m-1d0)*q(1))/m
      if(n.eq.-1) go to 26
      r(1)=q0
      r(2)=(ak*x1*q0+(m-1d0)*q(1))/m

c         print *, m,n,q0,q(1)
c         write(7,*) m,q(1)

      do 25 i=-2,n,-1
      r(1-i)=(ak**i*(2d0*i+3d0)*x1*r(-i)-
     *       (i-m+2d0)*r(-i-1))/(i+m+1d0)
   25 continue 
      res=r(-n+1)

c         print *,res

      go to 26
   24 continue
      res=q(n)
   26 continue
      if(id) 50,100,50
c***********************************    m = 0
   61 continue
      if(ic.eq.0) res=dlog((x1+1d0)/(x1-1d0))/2d0
      if(ic.ne.0.and.x1.lt.0.5d0) res=(datan(x1)-pi/2d0)
      if(ic.ne.0.and.x1.ge.0.5d0) res=-datan(1d0/x1)
      if(n.eq.0) go to 65
      q0=res
      nmax=1.5*x1+n+10
      r(nmax)=x1-dsqrt(x1**2-ak)
      do 62 i=nmax-1,1,-1
      r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
   62 continue
      q(1)=r(1)*q0*ak
      if(n.eq.1) GO TO 64
      do 63 i=1,n-1
      q(i+1)=r(i+1)*q(i)*ak**(i+1)
   63 continue
   64 continue
      res=q(n)
   65 if(id) 50,100,50
c***********************************    m = 1
   71 continue
      res=-ak/(dsqrt(x1**2-ak))
      if(n.eq.0) go to 75
      q0=res
      if(n.eq.-1) res=ak*x1*q0
      if(n.eq.-1) go to 75
      nmax=30*x1+n+50
c+&
      if(nmax.gt.3000) print *, 'dlegf2', nmax
      if(nmax.gt.3000) nmax=3000
      r(nmax)=x1-dsqrt(x1**2-ak)
      do 72 i=nmax-1,1,-1
      r(i)=((i+m)*ak)/((2d0*i+1d0)*x1-(i-m+1d0)*r(i+1))
   72 continue
      q(1)=r(1)*q0*ak
      if(n.eq.1) GO TO 74
      do 73 i=1,n-1
   73 q(i+1)=r(i+1)*q(i)*ak**(i+1)
   74 continue
      res=q(n)
   75 if(id) 50,100,50
   50 continue
c=========================             <<     Ø‡Æ®ß¢Æ§≠†Ô    >>
      if(id.eq.0) go to 100
      if(ic.eq.0.or.mod(n,2).eq.0)  go to 174
      red=-((n-m)*res+x1*n*x2)/(x1**2-ak)
      go to 100
  174 continue
      red=((n-m)*res-x1*n*x2)/(x1**2-ak)
  100 CONTINUE                                                      
      RETURN                                                        
      END

c************************************************
c GAMc2c1
c
      SUBROUTINE GAMc2c1(M,NN,GAM)
      parameter (nalphas=19,nterms=330)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 s, s1
      COMPLEX*16 BCOF,A,B,GAM,VV2,C,E,G
      DIMENSION A(4*nterms),B(4*nterms),
     *          BCOF(nterms,4*nterms), VV2(nterms),
     *          ACOF(nterms,4*nterms), VV1(nterms),
     *          GAM(nterms,nterms)
      COMMON /EPS3/ EPS3
      COMMON /K1/ S, S1, AKSI, AK, kK, NK, nal
      COMMON /Kc1/  ACOF, VV1
      COMMON /Kc2/  BCOF, VV2

      eps4 = 1d-40
      M2=2*M
       DO 10 N=1,NN
             NM=N+M-1
       DO 12 I=1,NK
   12 B(I)=BCOF(N,I)
      C=VV2(N)
       DO 11 L=1,NN
             LM=L+M-1
      IF(MOD(IABS(LM-NM),2)) 20,20,9
    9 CONTINUE
       DO 13 I=1,NK
   13 A(I)=aCOF(L,I)
      D=VV1(L)
      IF(MOD(LM-M,2))1,1,2
    1 I1=0
      I2=1
      K=M2
      E=A(1)*(B(1)/(M2+3D0))*2D0
      GO TO 3
    2 I1=1
      I2=0
      K=M2+1
      E=A(1)*(B(2)*(M2+2D0)/(M2+5D0)+B(1)/(M2+1D0))*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E=E*F
      DO 5 I=2,NK
      J=2*I-1-I2
c      AJ=2d0*I-1d0-I2
      JM2=2*(J+M)
c      AJM2=2d0*(aJ+M)
c      F=1D0
      F=a(i)
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
c    7 F=F*(aJ+JJ)
      IF(M2.EQ.0) F=1D0

c        print *,kk,m,nn,nk
c        print *,i, A(I)
c        print *,i, B(I+I1)
c        print *,i, B(I-I2)
c        print *,i, f
c        print *,i, i1, I+I1
c        print *,i, i2, I-I2
c        print *,e, g
c        print *,aj, aJM2-1D0
cc        print *,j, J+M2+1D0,JM2+3D0
cc        print *,(J+M2+1D0)/(JM2+3D0)
c        pause

      G=(B(I+I1)*(J+M2+1D0)/(JM2+3D0)+B(I-I2)*J/(JM2-1D0))*2D0/
     *(JM2+1D0)*F
c      G=A(I)*(B(I+I1)*(aJ+M2+1D0)/(aJM2+3D0)+B(I-I2)*aJ/(aJM2-1D0))
c     *  * 2D0 / (aJM2+1D0)*F
c      IF(CDABS(G).LE.1D-140.and.i.gt.nn)  GO TO 6
c      IF(DABS(G).LE.eps3.and.i.gt.nn) go to 6
      IF(cDABS(G/e).LE.eps4.and.i.gt.nn) go to 6

    5 E=E+G
    6 GAM(N,L) =E/C/D
      GO TO 11
   20 GAM(N,L) = s
   11 CONTINUE
   10 CONTINUE
      RETURN
      end

c************************************************
c GAMKAP2
c
      SUBROUTINE GAMKAP2(M,NN,GAM,AKAP)
      parameter (nalphas=19,nterms=330)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 BCOF,B,VV2,GAM,AKAP,C,E,G,E1,G1,A,D
      COMPLEX*16 s, s1
      DIMENSION A(4*nterms), B(4*nterms),
     *          BCOF(nterms,4*nterms), VV2(nterms),
     *          GAM(nterms,nterms), AKAP(nterms,nterms)
      COMMON /K1/ S, S1, AKSI, AK, kK, NK, nal
      COMMON /Kc2/  BCOF, VV2

c      eps4 = 1d-20
      eps4 = 1d-40
        M2=2*M
       DO 10 N=1,NN
          NM=N+M-1
       DO 12 I=1,NK
   12 B(I)=BCOF(N,I)
      C = VV2(N)
       DO 11 L=1,NN
          LM=L+M-1
      IF(MOD(IABS(LM-NM),2)) 20,20,9
    9 CONTINUE
       DO 13 I=1,NK
   13 A(I) = bCOF(L,I)
      D = VV2(L)
      IF(MOD(LM-M,2))1,1,2
    1 I1=0
      I2=1
      K=M2
      E=A(1)*(-B(1)*M/(M2+3D0))*2D0
      E1=A(1)*(B(1)/(M2+3D0))*2D0
      GO TO 3
    2 I1=1
      I2=0
      K=M2+1
      E=A(1)*(B(1)*(M+2D0)/(M2+1D0)-B(2)*(M+1D0)*(M2+2D0)/
     *(M2+5D0))*2D0/(M2+3D0)
      E1=A(1)*(B(1)/(M2+1D0)+B(2)*(M2+2D0)/(M2+5D0))*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E=E*F
      E1=E1*F
      DO 5 I=2,NK
      J=2*I-1-I2
      JM=J+M
      JM2=2*JM
c+&
c      F=1D0
      F=a(i)
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
      IF(M2.EQ.0) F=1D0
      G=(B(I-I2)*J*(JM+1D0)/(JM2-1D0)-B(I+I1)*JM*(M2+J+1D0)/
     *(JM2+3D0))*2D0/(JM2+1D0)*F
      G1=(B(I+I1)*(J+M2+1D0)/(JM2+3D0)+B(I-I2)*J/(JM2-1D0))
     *   *2D0/(JM2+1D0)*F
c      IF(CDABS(G).LE.1D-140.AND.CDABS(G1).LE.1D-140.and.i.gt.nn)
      IF(CDABS(G/e).LE.eps4.AND.CDABS(G1/e1).LE.eps4.and.i.gt.nn)
     *   GO TO 6
      E1=E1+G1
    5 E=E+G
    6 GAM(N,L)=E1/C/D
      AKAP(N,L)=E/C/D
      GO TO 11
   20 GAM(N,L)=s
      AKAP(N,L)=s
   11  CONTINUE
   10  CONTINUE
      RETURN
      END

c************************************************
C GJ1
c
      SUBROUTINE GJ1(A,N,NP,B,M,MP)
      implicit complex*16(a-h,o-z)
      PARAMETER (nterms=330, nalphas=19)
      DIMENSION A(NP,NP),B(mP,nP),IPIV(2*nterms),
     *          INDXR(2*nterms),INDXC(2*nterms)
      COMMON /K2/ Ialpha(nalphas),Jalpha(nalphas)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=(0.0d0,0.0d0)
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (cdABS(A(J,K)).GE.cdabs(BIG))THEN
                  BIG=cdABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
           IF(Jalpha(l).EQ.1) GO TO 15
            DUM=B(L,IROW)
            B(L,IROW)=B(L,ICOL)
            B(L,ICOL)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (cdabs(A(ICOL,ICOL)).EQ.(0.0d0,0.0d0))
     *  PAUSE 'Singular matrix.'
        PIVINV=(1.0d0,0.0d0)/A(ICOL,ICOL)
        A(ICOL,ICOL)=(1.0d0,0.0d0)
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
           IF(Jalpha(l).EQ.1) GO TO 17
          B(L,ICOL)=B(L,ICOL)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=(0.0d0,0.0d0)
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
             IF(Jalpha(l).EQ.1) GO TO 19
              B(L,LL)=B(L,LL)-B(L,ICOL)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END

c************************************************
C INT3
C
      SUBROUTINE INT3(M,NN,KAP,OMG,TU)
      parameter (nalphas=19,nterms=330)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KAP
      COMPLEX*16 s, s1
      DIMENSION A(4*nterms),B(4*nterms),ACOF(nterms,4*nterms),
     *          VV(nterms),OMG(nterms,nterms),TU(nterms,nterms),
     *          KAP(nterms,nterms)
       COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
       COMMON /kc1/ ACOF,VV
      DO 6 L=1,NN
      LM=L+M-1
       V1=VV(L)
      DO 12 I=1,NK
   12 A(I)=ACOF(L,I)
      DO 5 N =1,NN
      NM=N+M-1
       V2=VV(N)
      DO 13 I=1,NK
   13 B(I)=ACOF(N ,I)
          OMG(L,N)=OMEGA (M,LM,NM,NK,A,B,V1,V2)
          TU (L,N)=TAU   (M,LM,NM,NK,A,B,V1,V2)
          KAP(L,N)=AKAPPA(M,LM,NM,NK,A,B,V1,V2)
    5 CONTINUE
    6 CONTINUE
       RETURN
        END

c************************************************
c ks_new
c
      subroutine ks_new(K,ksi1,Vol_rat,ksi2)
      implicit real*8 (a-h,o-z)
      real*8 ksi1,ksi2,ksi20

      if(k.eq.0) then
* prolate
      ksi20=ksi1
* iter
101   ksi2=(ksi20+ksi1*(ksi1*ksi1-1d0)*Vol_rat)**(1d0/3d0)

      if(dabs(ksi2-ksi20).gt.1d-12) then
          ksi20=ksi2
          goto 101
      end if

      V_approx=ksi2*(ksi2*ksi2-1d0)/ksi1/(ksi1*ksi1-1d0)

         if(dabs(V_approx-Vol_rat).gt.1d-6) then
             print *,'ksi2 is bad',ksi2,ksi20,Vol_rat,V_approx
             stop
         end if

      else
* oblate
         if(Vol_rat.lt.1d0) then
           ksi20=ksi1
         else
           ksi20=(ksi1*(ksi1*ksi1+1d0)*Vol_rat)**(1d0/3d0)
         end if
* iter
102   ksi2=ksi20-(ksi20*(ksi20*ksi20+1d0)-ksi1*(ksi1*ksi1+1d0)*
     *    Vol_rat)/(3d0*ksi20*ksi20+1d0)

      if(dabs(ksi20-ksi2).gt.1d-12) then
      ksi20=ksi2
      goto 102
      end if

* check
      V_approx=ksi2*(ksi2*ksi2+1d0)/ksi1/(ksi1*ksi1+1d0)
        if(dabs(V_approx-Vol_rat).gt.1d-6) then
             print *,'ksi2 is bad',ksi2,Vol_rat,V_approx
             stop
        end if

      end if

       return
       end

c************************************************
c SIGma22
c
      SUBROUTINE SIGma22(M,NN,SIG)
      parameter (nterms=330)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 BCOF,B,A,VV2,SIG,C,E,G,D
      COMPLEX*16 s, s1
      DIMENSION A(4*nterms), B(4*nterms),
     *          BCOF(nterms,4*nterms), VV2(nterms), SIG(nterms,nterms)
      COMMON /EPS3/ EPS3
      COMMON /K1/ S, S1, AKSI, AK, kK, NK, nal
      COMMON /Kc2/ BCOF,VV2
c      eps4 = 1d-12
       M2=2*M
      DO 10 N=1,NN
        NM=N+M-1
      DO 12 I=1,NK
   12 B(I)=BCOF(N,I)
      C=VV2(N)
      DO 11 L=1,NN
        LM=L+M-1
      IF(MOD(IABS(LM-NM),2)) 9,9,20
    9 CONTINUE
      DO 13 I=1,NK
   13 A(I)=bCOF(L,I)
      D=VV2(L)
      IF(MOD(LM-M,2))1,1,2
    1 I2=1
      K=M2
      E=A(1)*(-B(2)*(M+1D0)*(M2+1D0)*(M2+2D0)/((M2+3D0)*(M2+5D0))
     *+B(1)*(2D0*M*M+3D0*M-2D0)/((M2-1D0)*(M2+3D0)))*2D0/(M2+1D0)
      GO TO 3
    2 I2=0
      K=M2+1
      E=A(1)*(-B(2)*(M+2D0)*(M2+2D0)*(M2+3D0)/((M2+5D0)*(M2+7D0))+
     *B(1)*(2D0*M*M+9D0*M+4D0)/((M2+1D0)*(M2+5D0)))*2D0/(M2+3D0)
    3 CONTINUE
      F=1D0
      DO 4 I=2,K
    4 F=F*I
      E=E*F
      DO 5 I=2,NK
      J=2*I-1-I2
      JM=J+M
      JM2=2*JM
c      F=1D0
c+&
      F=a(i)
      DO 7 JJ=1,M2
    7 F=F*(J+JJ)
      IF(M2.EQ.0) F=1D0
      G=(B(I-1)*J*(J-1D0)*JM/(JM2-3D0)/(JM2-1D0)-B(I+1)*(JM+1D0)*
     *(J+M2+1D0)*(J+M2+2D0)/(JM2+3D0)/(JM2+5D0)+B(I)*(JM*3D0*(JM+1D0)-
     *M*M-2D0)/(JM2-1D0)/(JM2+3D0))*2D0/(JM2+1D0)*F
      IF(CDABS(G).LE.1D-140.and.i.gt.nn) GO TO 6
c      IF(DABS(G).LE.eps3.and.i.gt.nn) GO TO 6
c      IF(cDABS(G/e).LE.eps4.and.i.gt.nn) GO TO 6
c      IF(cDABS(G/e).LE.eps4.and.i.gt.nn/2) GO TO 6
c      IF(cDABS(G/e).LE.eps4) GO TO 6
    5 E=E+G
    6 SIG(N,L)=E/C/D
      GO TO 11
   20 SIG(N,L)=s
   11 CONTINUE
   10 CONTINUE
      RETURN
      END

c************************************************
c crr
c
      SUBROUTINE Crr(N,A1,A2,A3)
      parameter (nalphas=19,nterms=330)
      IMPLICIT COMPLEX*16(A-H,O-Z)
      DIMENSION A1(nterms,nterms),A2(nterms),A3(nterms,nterms),
     *          A4(nterms,nterms)

* multiplication of complex matrix (a1), vector (a2) and a1^*

      DO I = 1, N
        DO J = 1, N
         a4(i,j) = a1(i,j) * a2(j)
        end do
      end do

      DO I = 1, N
        DO J = 1, N
         B=(0D0,0D0)
          DO 2 K=1,N
    2     B = B + A4(I,K) * A1(j,K)
         A3(I,J) = B
        end do
      end do

      RETURN
      END
c************************************************
C CXSC11b
C
      SUBROUTINE CXSC11b(II,M,C,NE,NN,QX,QS,QB,QEXT,QSCA,qbk)
      parameter (nterms=330, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z), COMPLEX*16 (R-S)
      COMPLEX*16 g, g1
      REAL*8 KAP
      DIMENSION B0(nterms),ACOF(nterms,4*nterms),VV(nterms),
     *          OMG(nterms,nterms),TU(nterms,nterms),
     *          BD(nalphas,2*nterms),KAP(nterms,nterms),IALFA(nalphas),
     *          JALFA(nalphas),RB0(nalphas,nterms),
     *          RA(nalphas,nterms),RB (nalphas,nterms),ASIN(nalphas),
     *          QX(nalphas),QS(nalphas),QB(nalphas)
      COMMON /EPS3/ EPS3
       COMMON /kc1/ ACOF,VV
       COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
       COMMON /K2/ IALFA, JALFa
       COMMON /K10/ B0
       COMMON /K11/ BD
       COMMON /FF/ ASIN,acos(nalphas)
       COMMON /RES/ RB0, RA, RB
        eps4 = 1d-100
c         eps4 = 1d-40
c         eps4 = 1d-12
       CALL INT3(M,NE,KAP,OMG,TU)

       IF(II.EQ.0) GO TO 101
       CK=C**2*DSQRT(AKSI)
       DO 20 IA=1,NAL
       DK=CK*DSQRT(AKSI+AK*ASIN(IA)**2)
       R=S
       Rr=S
       D=0D0
       IF(JALFA(IA).EQ.1) GO TO 20
       IF(M.GT.1) GO TO 23
       DO 21 L=1,NN
       V=VV(L)**2
        if(ii.eq.1) rb0(ia,l) = -rb0(ia,l)
       R1=sS1(-L)*RB0(IA,L)*BD(IA,L)*V
       Rr1=sS1(L+1)*RB0(IA,L)*BD(IA,L)*V
       D1=RB0(IA,L)*DCONJG(RB0(IA,L))*V
       R=R+R1
       Rr=Rr+Rr1
       D=D+D1
      IF(CDABS(R1/r).LT.eps4.AND.DABS(D1/d).LT.eps4.
     *   AND.IALFA(IA).NE.90)    GO TO 22
   21 CONTINUE
   22 D=2D0*D

ccccccccccccccccccccccccccccc       go to 26

   23 CONTINUE
       DO 24 L=1,NN
       LM=L+M-1
       R1=-sS1(-LM)*(S1*RA(IA,L)*BD(IA,L)*VV(L)**2-
     *    RB(IA,L)*BD(IA,L+NE))*ASIN(IA)
       Rr1=-sS1(LM)*(RA(IA,L)*BD(IA,L)*VV(L)**2-
     *    S1*RB(IA,L)*BD(IA,L+NE))*ASIN(IA)
       R=R+R1
       Rr=Rr+Rr1
       DO 25 N=1,NN
       NM=N+M-1
        RAN=DCONJG(RA(IA,N))
        RBN=DCONJG(RB(IA,N))
      D1=sS1(NM-LM)*(RA(IA,L)*RAN*OMG(L,N)+
     *  S1*(RB(IA,L)*RAN* KAP(L,N)-
     *  RA(IA,L)*RBN* KAP(N ,L))+
     *  RB(IA,L)*RBN*TU(L,N))*VV(N)*VV(L)
       D=D+D1
       if(dabs(d1).gt.1d300) print *,'cxsc11 ', l,n,vv(n),vv(l),d1,d
c       IF(DABS(D1/d).LT.eps4.AND.IALFA(IA).NE.90) GO TO 27
       IF(DABS(D1/d).LT.eps4.AND.n.gt.nn/2.and.d/dk.gt.1d-5) GO TO 27
c       IF(DABS(D1).LT.eps3.AND.IALFA(IA).NE.90) GO TO 27
   25 CONTINUE
c   27  IF(CDABS(R1).LT.eps3.AND.IALFA(IA).NE.90) GO TO 26
c   27  IF(CDABS(R1/r).LT.eps4.AND.IALFA(IA).NE.90) GO TO 26
  27  IF(CDABS(R1/r).LT.eps4.AND.l.gt.nn/2.and.cdabs(4d0*r/dk).gt.1d-5)
     * GO TO 26
   24 CONTINUE
   26 CONTINUE
      R = R * 4D0/DK
cc-----------------
      Rr = Rr * dconjg(rr)
      Rr = Rr * 4D0/DK
      QX(IA)=R
      QB(IA)=Rr
      QS(IA)=D/DK
   20 CONTINUE
        IF(II.nE.0) GO TO 100
c.................................................................
c     alpha = 0 deg.
c.................................................................
  101 CONTINUE
      D=0d0
      r=s
      rr=s
      CK=C**2*AKSI
      DO 6 L=1,NN
      DO 5 N =1,NN
        RAN=DCONJG(RA(1,N))
        RBN=DCONJG(RB(1,N))
      G=sS1(N-L)*(RA(1,L)*RAN*OMG(L,N)+S1*(RB(1,L)*RAN*KAP(L,N)-
     *  RA(1,L)*RBN*KAP(N,L))+RB(1,L)*RBN*TU(L,N))*VV(L)*VV(N)
      d=d+G
      IF(cDABS(g/d).LT.eps4.AND.n.gt.nn/2.and.d/ck.gt.1d-5) GO TO 7
cccccccccccccccccccccc      IF(cDABS(G).LE.eps3.AND.L.GT.3.and.d/ck.gt.1d-5) GO TO 7
c      IF(cDABS(G/d).LE.eps4.and.cdabs(g).lt.1d-150) GO TO 7
    5 continue
    7 G1=sS1(-L)*RB(1,L)*B0(L)
      rg1=sS1(L)*RB(1,L)*B0(L)
c       write(7,*) L,rg1
c       pause
      r=r+G1
      rr=rr+rG1
      IF(CDABS(g1/r).LT.eps4.AND.l.gt.nn/2.and.cdabs(4d0*r/ck).gt.1d-5)
     * GO TO 11
cccccccccccccccccccccc      IF(cDABS(G1).LE.eps3.AND.L.GT.3.and.cdabs(r/ck).gt.1d-5) GO TO 11
c      if(cDABS(G1/r).LE.eps4) GO TO 11
    6 CONTINUE
   11 QEXT=-2D0*r/CK
c      Qbk=2D0*rr*dconjg(rr)/CK
      Qbk=rr*dconjg(rr)/CK
      QSCA=d/CK
  100 RETURN
      END

c************************************************
C CmTM
C   TM mode, core-mantle spheroid
C................................................
C   Matrix, Eq.(36)
C ⁄ƒƒƒ                                                                            ƒƒƒƒƒø
C ≥ ksi1*rw1             ksi1*rw1             ksi1*rw1             ksi1*rw1            ≥
C ≥ G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3 ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ¿ƒƒƒƒ                                                                           ƒƒƒƒƒŸ
C................................................

      SUBROUTINE CMtm(II,M,rii,RI,C1,KSI1,ksi2,NN,ne)
      parameter (nterms=330, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16(R-S)
      REAL*8 ksi1, ksi2
      COMPLEX*16 p1, p3, DEL21, del23, GAM21, GAM22, akap22,
     *           AAA(2*nterms,2*nterms),BBB(nalphas,2*nterms),
     *           CCC(nalphas,2*nterms)
      DIMENSION
     *          BD(nalphas,2*nterms), B0(nterms),
     *          Ialpha(nalphas),Jalpha(nalphas),
     *          ASIN(nalphas), ACOS(nalphas),

     *          Rw3(nterms,nterms), Rw4(nterms,nterms),
     *          Rw7(nterms,nterms), Rw8(nterms,nterms),

     *          rfi(nalphas,2*nterms)

      COMMON /Kc1/  ACOF(nterms,4*nterms), VV1(nterms)
      COMMON /Kfun/ ar1_11(2*nterms), r3_11(2*nterms),
     *              aar1_11(nterms), rr3_11(nterms),
     *              r1_21(2*nterms), r3_21(2*nterms),
     *              rr1_21(nterms), rr3_21(nterms),

     *              r1_22(2*nterms), r3_22(2*nterms),
     *              rr1_22(nterms), rr3_22(nterms),
     *              r1_32(2*nterms), rr1_32(nterms),
     *              p1(nterms), p3(nterms)

      common /bbb/ bbb, Ra00(nalphas,nterms), Ra0(nalphas,nterms)
     *             ,rfi

      common /ires/ Rw1(nterms,nterms), Rw2(nterms,nterms),
     *              Rw5(nterms,nterms), Rw6(nterms,nterms),
     *              RA1(nterms,nterms), RA3(nterms,nterms),
     *              Rg1(nterms,nterms), Rg3(nterms,nterms),
     *              Rf1, Rf2, Rf3

      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /K2/ Ialpha, Jalpha
      COMMON /K10/ B0
      COMMON /K11/ BD
      COMMON /FF/  ASIN, ACOS
      COMMON /INT/ DEL21(nterms,nterms), del23(nterms,nterms),
     *             GAM21(nterms,nterms), GAM22(nterms,nterms),
     *             akap22(nterms,nterms), sig22(nterms,nterms)
      COMMON /RES/ Rz0(nalphas,nterms),
     *             Rz1(nalphas,nterms),Rz2(nalphas,nterms)

      NDIM = 2*nterms
      NP2=NN/2
      rf1 = ri**2 - 1d0
      rf3 = ri**2 / rii**2
      rf2 = rf3 - 1d0

c ra1 --> ~R_1 (Eq.(20))
c ra3 --> ~R_3 (Eq.(20))

c       print *, rr1_32(1)
c       print *, rr1_32(2)
c       print *, del23(1,1)
c       print *, del23(1,2)
c       print *, del23(1,3)
c       write(7,*) rr1_32(1)
c       write(7,*) rr1_32(2)
c       write(7,*) del23(1,1)
c       write(7,*) del23(1,2)
c       write(7,*) del23(1,3)
c       pause

      CALL Crr(Ne,del23,rr1_32,RA1)
      CALL Crr(Ne,del21,rr3_11,RA3)

         IF(M.GT.1.OR.II.EQ.0) GO TO 20
*
*  Axisymmeric part

      f3 = ksi2**2 - ak

c**
c** matrix Eq.(20)
c**
      DO  N = 1, NN
      n1 = n + nn
        DO  L = 1, NN
        l1 = l + nn
c--------
         AAA(N ,L ) = aksi * (rf1 + 1d0) * ra3(n,l)
       if(n.eq.l) AAA(N,L) = AAA(N,L) - aksi *  rr1_21(l) + rf1 * ksi1
         AAA(N ,L1) = aksi * (rf1 + 1d0) * ra3(n,l)
       if(n.eq.l) AAA(N,L1) = AAA(N,L1) - aksi * rr3_21(l) + rf1 * ksi1
c--------
         AAA(N1,L ) =  f3 * rf3 * ra1(n,l) * p1(l)
       if(n.eq.l) AAA(N1,L) =  AAA(N1,L) +
     *            (rf2 * ksi2 - f3 * rr1_22(l)) * p1(l)
         AAA(N1,L1) =  f3 * rf3 * ra1(n,l) * p3(l)
       if(n.eq.l) AAA(N1,L1) =  AAA(N1,L1) +
     *            (rf2 * ksi2 - f3 * rr3_22(l)) * p3(l)
c--------
        end do
      end do

       DO 23 IA = 1, NAL
         DO N = 1, Ne
           RA00(ia,n) = 2D0 * S1**N * BD(IA,N) * ar1_11(n)
         end do

        DO N = 1, NN
          S2=S
           DO i = 1, NN
            S2 = S2 + DEL21(n,i) * (Rr3_11(I) - aar1_11(i)) * RA00(ia,i)
           end do
          ccc(ia,N) = ri**2 * aksi * S2
          ccc(ia,N+nn) = S
        end do
   23  CONTINUE

          CALL GJ1(AAA,2*nn,NDIM,ccc,nal,nalphas)

* calculation of Z1 (Eq.(22))

      do ia = 1, nal
      do i = 1, nn
         s2 = s
         do jj = 1, nn
           s2 = s2 + DEL21(jj,i) * (ccc(ia,jj) + ccc(ia,jj+nn))
         end do
         rz0(ia,i) = (s2 - RA00(ia,i)) / R3_11(i)
      end do
      end do

c**
*  Non-axisymmeric part
c**
   20     CONTINUE

      CALL Cmlm(Ne,gam22,ra1,rg1)
      CALL Cmlm(Ne,gam22,ra3,rg3)

      f1 = ksi1 / (ksi1**2 - ak)
      f2 = ak / (ksi1**2 - ak)
      f3 = ksi2 / (ksi2**2 - ak)
      f4 = ak / (ksi2**2 - ak)

      do i = 1, ne
         do j = 1, ne
            rw1(i,j) = ra3(i,j)
            if(i.eq.j) rw1(i,j) = rw1(i,j) - rr1_21(i)
            rw2(i,j) = ra3(i,j)
            if(i.eq.j) rw2(i,j) = rw2(i,j) - rr3_21(i)
            rw3(i,j) = rg3(i,j) + f1 * akap22(i,j)
            rw4(i,j) = ksi1 * ra3(i,j) + f2 * sig22(i,j)
            if(i.eq.j) rw4(i,j) = rw4(i,j) + 1d0

            rw5(i,j) = ra1(i,j)
            if(i.eq.j) rw5(i,j) = rw5(i,j) - rr1_22(i)
            rw6(i,j) = ra1(i,j)
            if(i.eq.j) rw6(i,j) = rw6(i,j) - rr3_22(i)
            rw7(i,j) = rg1(i,j) + f3 * akap22(i,j)
            rw8(i,j) = ksi2 * ra1(i,j) + f4 * sig22(i,j)
            if(i.eq.j) rw8(i,j) = rw8(i,j) + 1d0
         end do
      end do

      IF(II.EQ.0) then
         DO  L = 1, Nn
            RA0(1,L) = 2D0*S1**(L-1)*B0(L)*ar1_11(L)/VV1(L)**2
         end do
      end if

      IF(II.nE.0) then
         DO IA = 1, NAl
           IF(Jalpha(IA).EQ.1) GO TO 15
            DO N = 1, Ne
              RA0(ia,n) = 4D0*S1**(N+M-2)*BD(IA,N)*ar1_11(N)/ASIN(IA)
            end do
   15      CONTINUE
         end do
      end if

c**
c** matrix Eq.(36)
c**
      DO 5 IV=1,2
      IG=IV-1
      DO 8 N=1,NP2
      N1=N  + NP2
      N2=N1 + NP2
      N3=N2 + NP2
      ND=2*(N-1)+IV
      NG=2*N-IG
      DO 9 L=1,NP2
      L1=L + NP2
      L2=L1+ NP2
      L3=L2+ NP2
      LD=2*(L-1)+IV
      LG=2*L-IG
c--------
      AAA(N ,L ) = ksi1 * rw1(nd,ld)
      AAA(N ,L1) = ak * (rg3(nd,lg) - gam22(nd,lg) * rr1_21(lg))
      AAA(N ,L2) = ksi1 * rw2(nd,ld)
      AAA(N ,L3) = ak * (rg3(nd,lg) - gam22(nd,lg) * rr3_21(lg))
c--------
      AAA(N1,L ) = (rg3(ng,ld) - gam22(ng,ld) * rr1_21(ld))
     *             + rf1 * rw3(ng,ld)
      AAA(N1,L1) = ksi1 * rw1(ng,lg) + rf1 * rw4(ng,lg)
      AAA(N1,L2) = (rg3(ng,ld) - gam22(ng,ld) * rr3_21(ld))
     *             + rf1 * rw3(ng,ld)
      AAA(N1,L3) = ksi1 * rw2(ng,lg) + rf1 * rw4(ng,lg)
c--------
      AAA(N2,L ) = ksi2 * rw5(nd,ld) * p1(ld)
      AAA(N2,L1) = ak * (rg1(nd,lg) - gam22(nd,lg) * rr1_22(lg))
     *             * p1(lg)
      AAA(N2,L2) = ksi2 * rw6(nd,ld) * p3(ld)
      AAA(N2,L3) = ak * (rg1(nd,lg) - gam22(nd,lg) * rr3_22(lg))
     *             * p3(lg)
c--------
      AAA(N3,L ) = ((rg1(ng,ld) - gam22(ng,ld) * rr1_22(ld))
     *             + rf2 * rw7(ng,ld)) * p1(ld)
      AAA(N3,L1) = (ksi2 * rw5(ng,lg) + rf2 * rw8(ng,lg)) * p1(lg)
      AAA(N3,L2) = ((rg1(ng,ld) - gam22(ng,ld) * rr3_22(ld))
     *             + rf2 * rw7(ng,ld))* p3(ld)
      AAA(N3,L3) = (ksi2 * rw6(ng,lg) + rf2 * rw8(ng,lg)) * p3(lg)
c--------
    9 CONTINUE
    8 CONTINUE

      NA=NAL
      IF(II.EQ.0) na = 1
      DO 11 IA = 1, NA
      IF(Jalpha(IA).EQ.1) GO TO 11
      DO N = 1, Ne
       S2=S
       S3=S
         DO I = 1, Ne
          S2 = S2 + DEL21(n,i) * (Rr3_11(I) - aar1_11(i)) * RA0(ia,i)
          S3 = S3 + gam21(n,i) * (Rr3_11(I) - aar1_11(i)) * RA0(ia,i)
         end do
       BBB(ia,N) = ksi1 * S2
       BBB(ia,N+Ne) = (rf1 + 1d0) * S3
      end do

      DO 13 N=1,NP2
        ND=2*(N-1)+IV
        NG=2*N-IG
        CCC(ia,N    )   = BBB(ia,ND)
        CCC(ia,N+NP2)   = BBB(ia,NG+Ne)
        CCC(ia,N+2*NP2) = s
        CCC(ia,N+3*NP2) = s
   13 continue

   11 CONTINUE
c---
          CALL GJ1(AAA,2*nn,NDIM,ccC,na,nalphas)
c---
      if(iv.eq.1) then
        do ia = 1, na
         do n = 1, 2*nn
           rfi(ia,n) = ccc(ia,n)
         end do
        end do
      end if

    5 CONTINUE

* calculations of Z1 / Z2

      do ia = 1, na
       IF(Jalpha(IA).EQ.1) GO TO 111
       do i = 1, nn
        s2 = s
        s3 = s
        do jj = 1, nn
         if(MOD(jj,2).eq.0) then
        s2 = s2 + DEL21(jj,i) * (ccc(ia,jj/2) + ccc(ia,jj/2+2*np2))
        s3 = s3 + DEL21(jj,i) *
     *            (rfi(ia,jj/2+np2) + rfi(ia,jj/2+3*np2))
       else
        s2 = s2 + DEL21(jj,i) * (rfi(ia,jj/2+1) + rfi(ia,jj/2+1+2*np2))
        s3 = s3 + DEL21(jj,i) *
     *            (ccc(ia,jj/2+1+np2) + ccc(ia,jj/2+1+3*np2))
         end if
        end do
        rz1(ia,i) = (s2-RA0(ia,i)) / R3_11(i)
        rz2(ia,i) = s3 / R3_11(i) / C1
       end do
 111   continue
      end do

      RETURN
      END


c******************************************
      SUBROUTINE DUGOLM1(X,M,LMAX,PI,TAU)
      REAL*8 X,PI(LMAX),TAU(LMAX),A
      M2=2*M
      A=1.D0
      DO I=3,M2-1
       IF(MOD(I,2).EQ.1) A=A*I
      ENDDO
      PI(1)=A
      PI(2)=A*X*(M2+1.D0)
      TAU(1)=A*M*X
      TAU(2)=A*(M2+1.D0)*((M+1.D0)*X*X-1.D0)
      DO L1=3,LMAX
        PI(L1)=((2.D0*L1+M2-3.D0)*X*PI(L1-1)-(L1+M2-2.D0)*PI(L1-2))
     *          /(L1-1.D0)
      TAU(L1)=(L1+M-1.D0)*X*PI(L1)-(L1+M2-1.D0)*PI(L1-1)
      ENDDO
      RETURN
      END

c************************************************
C funlegnn
C
      SUBROUTINE funlegnn(M,X,NF)
      parameter (nterms=330)
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 s, s1
      DIMENSION P(2*nterms),PD(2*nterms),Q(2*nterms),QD(2*nterms)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /LEG/ P, PD, Q, QD
      IF(NF.GT.2*nterms) WRITE(7,1) NF
    1 FORMAT(1X,'funlegnn  NF > 2*nterms, NF=',2I10)
      IF(NF.GT.2*nterms) GO TO 900
         DO I  = 1, 2*nterms
         P (I) = 0D0
         PD(I) = 0D0
         Q (I) = 0D0
         QD(I) = 0D0
         end do
      if(k.eq.0) then
        CALL DUGOLM1(X,M,Nf,P,PD)
      IF(X.LT.1.D0) THEN
        A1=DSQRT(1.D0-X*X)
      ELSE
        A1=DSQRT(X*X-1.D0)
      ENDIF
      DO I = 1, Nf
        P(I)=P(I)*A1**M
        PD(I)=-PD(I)*A1**(M-2)
        IF(X.GT.1.D0) PD(I)=-PD(I)
      ENDDO
      end if
       ier1 = 0
      if(k.eq.1) CALL DLEGF1(P(1),D,X,D,d,M,M,K,0,IER1)
      CALL DLEGF2(Q(1),D,X,D,-M,M,K,0,IER2)
       IF((IER1+IER2).NE.0) WRITE(7,2) M,IER1,IER2
c       IF(IER2.NE.0) WRITE(7,2) M,IER2
    2 FORMAT(1X,'funlegnn  M,I,IER1,2,3,4=',6I10)
       ier3 = 0
       DO 3 I=2,NF
      if(k.eq.1) then
       if(i.gt.2) d=P(I-2)
      CALL DLEGF1(P(I),PD(I-1),X,P(I-1),d,M+I-1,M,K,1,IER3)
      end if
      CALL DLEGF2(Q(I),QD(I-1),X,Q(I-1),-M+I-1,M,K,1,IER4)
c       IF((IER3+IER4).NE.0) WRITE(7,2) M,I,IER1,IER2,IER3,IER4
       IF(IER4.NE.0) WRITE(7,2) M,I,IER4
    3 CONTINUE
  900 RETURN
      END
c************************************************
C funlegnn1
C
      SUBROUTINE funlegnn1(M,NF)
      parameter (nterms=330,nalphas=19)
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 s, s1
      DIMENSION P(nalphas,2*nterms),PD(nalphas,2*nterms),
     *          Pp(2*nterms),PDp(2*nterms)
      COMMON /LEG1/ P,PD
      COMMON /K1/ S, S1, AKSI, AK, kkK, NK, nal
      COMMON /FF/ ASIN(nalphas),ACOS(nalphas)
      IF(NF.GT.2*nterms) WRITE(7,1) NF
      IF(NF.GT.2*nterms) WRITE(*,1) NF
    1 FORMAT(1X,'funlegnn1 NF > 250, NF=',2I10)
      IF(NF.GT.2*nterms) stop 666
         DO j  = 1, nalphas
         DO I  = 1, 2*nterms
         P (j,I) = 0D0
         PD(j,I) = 0D0
         end do
         end do
      do 6 ia=1,nal
      x=acos(ia)
        CALL DUGOLM1(X,M,Nf,Pp,PDp)
        A1=DSQRT(1.D0-X*X)
      DO I = 1, Nf
        P(ia,I)=Pp(I)*A1**M
        PD(ia,I)=-PDp(I)*A1**(M-2)
      ENDDO
    6 CONTINUE
  900 RETURN
      END


c************************************************
c CDRG1cn
c
      SUBROUTINE CDRG1cn(R,DR,AR,ADR,I12,M,n,C,x,AD,IA,BD,IB,rrl,IER)
      IMPLICIT REAL*8(A-H,O-Q,T-Z),COMPLEX*16(R-S)
      COMPLEX*16 C, dr, drr, drr2, dr2, AD(IA), drr3
      COMPLEX*16 ar, adr, bd
      DIMENSION SG(1000), sa(1000), sb(1000), BD(IB+1),
     *          sam(100)
      COMMON /K1/ S, S1, AKSIi, AK, K, NK, nal
      COMMON /FACT/ FACT(300)
      COMMON /eps3/ eps3
      ier=0
c      if(i12.eq.0) go to 900
c! lambda according Slavyanov
      rl = rrl - c**2
c! aksi!
      if(k.eq.0) aksi = x**2 - 1d0
      if(k.eq.1) aksi = x**2 + 1d0
c+&      eps4 = 1d-20
      eps4 = 1d-100
      aM = 0.5d0 * M
      if(k.eq.0) xx0 = ((x - 1d0) / (x + 1d0)) ** am
      if(k.eq.0) xx =  (x - 1d0) / (x + 1d0)

c      print *,k
c      print *,x, xx
c      pause

      rp = - S1 * c
      re = cdexp(-rp * (x - 1d0))
      rp2 = 2d0 * rp
      rrp2l = (rp2 + 1d0) * (m + 1d0) + rl

c coefficients
      M2=2*M
      nm1=n-m
      nm2=n+m
      IF(MOD((N-M),2)) 1,2,1
    1  CONTINUE
      JK=1
      GO TO 36
    2  CONTINUE
      JK=0
   36 CONTINUE

c       if((nm1-jk)/2+1.gt.170) stop 1707
      ff3 = fact3(nm2+jk,(nm2+jk)/2,1d0/fact((nm1-jk)/2+1))
      rg0=(2d0*m+1d0+jk*2d0)/2d0**(nm2+1d0) * ff3 / AD(1)/C**(M+JK)
      rg0 = 1d0 / rg0

c 2===
      IF(I12.ne.0) then
      IF(MOD(m,2)) 33,34,33
   33  CONTINUE
      Jm=1
      GO TO 35
   34  CONTINUE
      Jm=0
   35 CONTINUE
      SIG=S
      DO 7 I=1,IA
      J=2*(I-1)+JK
      rf2 = AD(I) * fact2(m2+j,j)
      SIG=SIG+rf2
      if(cdabs(rf2/sig).lt.eps4) go to 8
    7  CONTINUE
    8  CONTINUE
      RK2=2d0**nm1/(m2-1d0)*fact(m2+1)/fact(m+1) / ff3
     *    * bd(m+1)*sig/C**(M-1-JK)
      if(jk.eq.1.and.k.eq.0) rk2=-rk2/(m2-3d0)
      if(k.eq.0) go to 339
      if(jk.eq.1.and.jm.eq.0) rk2=rk2/(m2-3d0)
      if(jk.eq.1.and.jm.eq.1) rk2=-rk2/(m2-3d0)*(-s1)
      if(jk.eq.0.and.jm.eq.0) rk2=rk2*s1
      if(k.eq.1) rk2=rk2/(-s1)**(M-1-JK)
  339 continue
      end if
c 2===
c------------
      Sssa0 = S
c. 1s
      DO I = 1, M
        j = - 2 * (M - I + 1) + JK
        if(j.ge.-m) go to 351
      sssa0 = sssa0 + BD(m-I+2) * (-1)**(-j-1) * fact(-j)
     *        * fact(m2+j+1)
      end do

351   continue

c. 2s
      rf1 = s
      DO I = 1, IA
      ir = 2* (i - 1) + jk
      rf3 = s
c=
         do j = 1, m
           ff2 = fact(m+1) / fact(m-j+1)
           rf2 = s
           do jj = 1, j+1
            jjj = jj - 1
            rf2 = rf2 + fact2(m,m-jjj)**2 / fact(jjj+1) /
     *            fact(j-jjj+1) * fact2(ir+m2-jjj,ir+jjj)
           end do
           rf2 = rf2 * (-1)**j / (2d0 * j)  * ff2
           rf3 = rf3 + (-1)**(j-1) / (2d0 * j)  * ff2**2 *
     *           fact2(ir+m2-j,ir+j) + rf2
         end do
c=
         irm = (ir + m - 1 + jk) / 2
         rf4 = s
         if(ir.eq.0) go to 295
         do j = 1, irm
           jj = j - 1
           if(2*jj.gt.ir-1) go to 295
           rf4 = rf4 + (2d0 * ir + m2 - 4d0 * jj - 1d0) / (ir + m - jj)
     *           / (2d0 * jj + 1d0) * fact2(ir+m2-2*jj-1,ir-2*jj-1)
         end do
 295   CONTINUE
      rf11 =  AD(I) * (rf3 - rf4)
      rf1 = rf1 + rf11
      if(cdabs(rf11/rf1).lt.eps4) go to 297
      end do
 297   CONTINUE
      Sssa0 = Sssa0 + rf1
      rs2 = rf1

c. 3s
      rf1 = S
      MR = M + 1
      DO I = MR, IB
        J = 2 * I  - JK
        rf11 = bD(I+1) * fact2(j-1,j-m2-1)
        rf1 = rf1 + rf11
      if(cdabs(rf11/rf1).lt.eps4) go to 298
      end do

  298 continue

      Sssa0 = Sssa0 + rf1

c 1--
      SG0 = (1d0,0d0)
      SG(1) = rrp2l / (m + 1d0) * sg0
      SG(2) = ((2d0 * (m + 2d0 + rp2) + rrp2l) * sg(1) -
     *        (m + 1d0) * sg0) / (2d0 * (m + 2d0))
      r = sg0 + sg(1) * xx + sg(2) * xx**2
c      r = (sg0 + sg(1) * xx + sg(2) * xx**2) * re
      drr = sg(1) * 1d0 +  sg(2) * 2d0 * xx
c      drr = (sg(1) * 1d0 +  sg(2) * 2d0 * xx) * re

c 2--
      if(i12.ne.0) then

      do i = 1, 100
         Sam(i) = s
      end do

      Sam(1) = - 2d0 / m / c / rg0
      if(-m+1.eq.0) go to 352
      Sam(2) = (2d0 * m * (rp2 + 1d0) - rrp2l) / (m - 1d0) * sam(1)
      if(-m+2.eq.0) go to 352
       do i = 1, m-2
         j = - m + i
         Sam(i+2) = ((2d0 * j * (j + rp2 + m + 1d0) + rrp2l)
     *              * Sam(i+1) -
     *           (j * (j + m) * sam(i))) / ((j + 1d0) * (j + m + 1d0))
       end do

  352 continue

c a0--
      s2ic = 2d0 * s1 * c
      sssa1 = s
       do i = 1, m
       j = i - 1
       rf = s
         do ii = 1, m - j
           rf = rf + s2ic**ii / fact(ii+1) * fact(m-j) / fact (ii)
     *          / fact(m-ii-j+1)
         end do
       sssa1 = sssa1 + sam(i) * rf
       end do

      Sa00 = 2d0 / fact(m+1) * sssa0 / rk2 - sssa1

c a0--
c!!!!!!!!          Sa0 = sa00

      Sb0 = (-2d0 * (rp2 + m) + rrp2l) * Sam(m) / m
      if(m.gt.1) sb0 = sb0 + (m - 1d0) * sam(m-1) / m

      sb(1) = rrp2l / (m + 1d0) * sb0
      sb(2) = ((2d0 * (m + 2d0 + rp2) + rrp2l) * sb(1) -
     *        (m + 1d0) * sb0) / (2d0 * (m + 2d0))
      Sa(1) = (rrp2l * sa00 - (m + 2d0) * sb(1)
     *        + 2d0 * (rp2 + m + 1d0) * sb0) / (m + 1d0)
      Sa(2) = ((2d0 * (m + 2d0 + rp2) + rrp2l) * sa(1)
     *        - (m + 1d0) * sa00 - ((m + 4d0) * sb(2)
     *        - 2d0 * (rp2 + m + 3d0) * sb(1)
     *        + (m + 2d0) * sb0)) / (2d0 * (m + 2d0))

      r2 = s
      drr2 = s
       do i = 1, m
         r2 = r2 + sam(m-i+1) / xx**i
         drr2 = drr2 - i * sam(m-i+1) / xx**(i+1)
       end do

      r2 = r2 + sa00 + sa(1) * xx + sa(2) * xx**2
      rr2 = sb0 + sb(1) * xx + sb(2) * xx**2
      drr2 = drr2 + sa(1) * 1d0 +
     *       2d0 * sa(2) * xx
      drr3 = sb(1) * 1d0 + 2d0 * sb(2) * xx
c      r2 = r2 + sa00 + sa(1) * xx + sa(2) * xx**2
c      r2 = r2 * re
c      rr2 = (sb0 + sb(1) * xx + sb(2) * xx**2) * re
c      drr2 = drr2 + sa(1) * 1d0 +
c     *       2d0 * sa(2) * xx
c      drr2 = drr2 * re
c      drr3 = (sb(1) * 1d0 + 2d0 * sb(2) * xx) * re
      end if

c      cosn = cdcos(c*(x-1d0))
c      csin = cdsin(c*(x-1d0))

      DO  I = 2, 998
c 1--
      a = (i + 1d0) * (i + m + 1d0)
      rb = 2d0 * i * (i + m + 1d0 + rp2) + rrp2l
      g = i * (i + m)

      sg(i+1) = (rb * sg(i) - g * sg(i-1)) / a
      rs = sg(i+1) * xx**(i+1d0)
c      drs = dreal(rs)
c      airs = dimag(rs)
c      rsn = drs * cosn - airs * csin + s1 * (drs * csin + airs * cosn)

c      if(dreal(rs).ge.0d0) Rr01 = Rr01 + RS
c      if(dreal(rs).lt.0d0) Rr02 = Rr02 + RS
      R = R + RS
c       R = R + RSn

      rd = sg(i+1) * (i + 1d0) * xx**i
c      drn = dreal(rd)
c      airn = dimag(rd)
c      rdn = drn * cosn - airn * csin + s1 * (drn * csin + airn * cosn)

      dRr = dRr + Rd
c      dRr = dRr + Rdn
c 2--
      rs2 = s
      rrs2 = s
      rd21 = s
      rd23 = s

      if(i12.ne.0) then

        if(cdabs(sb(i)).gt.1d300) return

      sb(i+1) = (rb * sb(i) - g * sb(i-1)) / a
      sa(i+1) = (rb * sa(i) - g * sa(i-1) - ((2d0 * i + m + 2d0) *
     *          sb(i+1) - 2d0 * (2d0 * i + rp2 + m + 1d0) * sb(i) +
     *          (2d0 * i + m) * sb(i-1))) / a
      rs2 = sa(i+1) * xx**(i+1d0)
      rrs2 = sb(i+1) * xx**(i+1d0)

c      drs2 = dreal(rs2)
c      airs2 = dimag(rs2)
c      rs2n = drs2 * cosn - airs2 * csin + s1 * (drs2 * csin
c     *       + airs2 * cosn)

c      drs22 = dreal(rrs2)
c      airs22 = dimag(rrs2)
c      rrs2n = drs22 * cosn - airs22 * csin + s1 * (drs22 * csin
c     *       + airs22 * cosn)

      R2 = R2 + RS2
c      R2 = R2 + RS2n
      Rr2 = Rr2 + RrS2
c      Rr2 = Rr2 + RrS2n

      rd21 = sa(i+1) * (i + 1d0) * xx**i
      rd23 = sb(i+1) * (i + 1d0) * xx**i
c      drd21 = dreal(rd21)
c      aird21 = dimag(rd21)
c      rd21n = drd21 * cosn - aird21 * csin + s1 * (drd21 * csin
c     *       + aird21 * cosn)
c      drd23 = dreal(rd23)
c      aird23 = dimag(rd23)
c      rd23n = drd23 * cosn - aird23 * csin + s1 * (drd23 * csin
c     *       + aird23 * cosn)

      dRr2 = dRr2 + rd21
      dRr3 = dRr3 + rd23
c      dRr2 = dRr2 + rd21n
c      dRr3 = dRr3 + rd23n

c      if(CDABS(RS2n/r2).LT.eps4.AND.CDABS(Rrs2n/rr2).LT.eps4) GO TO 15
      if(CDABS(RS2/r2).LT.eps4.AND.CDABS(Rrs2/rr2).LT.eps4) GO TO 15
      end if
   15 CONTINUE
c      IF(CDABS(RS/r).LT.eps4.AND.CDABS(RDn/drr).LT.eps4) GO TO 14
      IF(CDABS(RS/r).LT.eps4.AND.CDABS(RD/drr).LT.eps4) GO TO 14
      end do
   14 CONTINUE
c 1--

c       R = R + rr01 + rr02

      R = R * xx0 / (x + 1d0) * re * rg0

c      R = R * xx0 / (x + 1d0) * rg0

c      aaR = R
c      R = aaR
      dr = r * (-rp  - 1d0/(x + 1d0) + m / aksi) +
     *     DRr * (2d0 / (x + 1d0)**3) * xx0 * re * rg0
c      dr = r * (-rp  - 1d0/(x + 1d0) + m / aksi) +
c     *     DRr * (2d0 / (x + 1d0)**3) * xx0 * rg0
c      aaR = dR
c      dR = aaR
c 2--
      if(i12.ne.0) then
      R2 = (rr2 * dlog(xx) + R2) * xx0 / (x + 1d0) * re
c      R2 = (rr2 * dlog(xx) + R2) * xx0 / (x + 1d0)
      dr2 = r2 * (-rp -1d0/(x + 1d0) + m / aksi) +
     *     (DRr2 + rr2 / xx + drr3 * dlog(xx)) *
     *     (2d0 / (x + 1d0)**3) * xx0 * re
c      dr2 = r2 * (-rp -1d0/(x + 1d0) + m / aksi) +
c     *     (DRr2 + rr2 / xx + drr3 * dlog(xx)) *
c     *     (2d0 / (x + 1d0)**3) * xx0
      AR = r2
      ADR = dr2
      end if

      if(dimag(c).eq.0d0) then
      r = dreal(r)
      dr = dreal(dr)
      ar = dreal(ar)
      adr = dreal(adr)
      end if

  900 RETURN
      END
c************************************************
C cdrkf45a
c
      SUBROUTINE cDRKF45a(F,NEQN,Y,T,TOUT,
c     1   RELERR,ABSERR,IFLAG,WORK,IWORK)
     1   RELERR,ABSERR,IFLAG)
c      INTEGER IWORK(5),NEQN
      INTEGER NEQN
c      REAL*8 Y(NEQN),T,TOUT,RELERR,ABSERR,
c     1   WORK(27)
      REAL*8 T,TOUT,RELERR,ABSERR
c      complex*16 Y(NEQN), WORK(27)
      complex*16 Y(NEQN)
      EXTERNAL F
c      K1M = NEQN + 1
c      K1 = K1M +1
c      K2 = K1 + NEQN
c      K3 = K2 + NEQN
c      K4 = K3 + NEQN
c      K5 = K4 + NEQN
c      K6=K5+NEQN
c      CALL cDRKFS (F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK(1),
c     1   WORK(K1M),WORK(K1),WORK(K2),WORK(K3),WORK(K4),WORK(K5),
c     2   WORK(K6),WORK(K6+1),IWORK(1),IWORK(2),IWORK(3),
c     3   IWORK(4),IWORK(5))
      CALL cDRKFSa (F,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG)
      RETURN
      END

c************************************************
C cdrkfsa
c
      SUBROUTINE cDRKFSa(F,NEQN,Y,T,TOUT,RELERR,ABSERR,
     1   IFLAG)
      LOGICAL  HFAILD, OUTPUT
c      complex*16 Y(NEQN),YP(NEQN),F1(NEQN),
c     *           F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN)
      complex*16 Y(NEQN),YP(4),F1(4),
     *           F2(4),F3(4),F4(4),F5(4)
      REAL*8 T,TOUT,RELERR,ABSERR,h
c      REAL*8 Y(NEQN),T,TOUT,RELERR,ABSERR,H,YP(NEQN),F1(NEQN)
     1   ,SAVRE,SAVAE
      EXTERNAL F
      REAL*8  A,AE,DT,EE,EEOET,ESTTOL,ET,HMIN,REMIN,
     1   RER,S,SCALE,TOL,TOLN,U26,EPSP1,EPS,YPK
c   * added (Draine!)
      SAVE EPS,MAXNFE,REMIN,U26
      DATA EPS/1.D0/
c   *
      DATA REMIN /1.0D-16/
      DATA MAXNFE /500000/
c  *
      IF(EPS.GE.1.D0)THEN
    5    EPS=EPS/2.0D0
         EPSP1=EPS+1.0D0
         IF(EPSP1.GT.1.0D0)GOTO 5
         U26=26.0D0*EPS
      ENDIF
c  **************************************** (Draine!)
      IF(NEQN.LT.1) GO TO 10
      IF((RELERR.LT.0D0).OR.(ABSERR.LT.0D0)) GO TO 10
      MFLAG = IABS(IFLAG)
      IF((MFLAG.EQ.0).OR.(MFLAG.GT.8))  GO TO 10
      IF(MFLAG.NE.1) GO TO 20
c      EPS = 1.0D0
c    5 EPS = EPS/2.0D0
c      EPSP1 = EPS + 1.0D0
c      IF(EPSP1.GT.1.0D0) GO TO 5
c      U26 = 26.0D0 * EPS
      GO TO 50
   10 IFLAG =8
      RETURN
   20 IF((T.EQ.TOUT).AND.(KFLAG.NE.3)) GO TO 10
      IF(MFLAG.NE.2) GO TO 25
      IF((KFLAG.EQ.3).OR.(INIT.EQ.0)) GO TO 45
      IF(KFLAG.EQ.4) GO TO 40
      IF((KFLAG.EQ.5).AND.(ABSERR.EQ.0.0D0)) GO TO 30
      IF((KFLAG.EQ.6).AND.(RELERR.LE.SAVRE).AND.
     1   (ABSERR.LE.SAVAE)) GO TO 30
      GO TO 50
   25 IF(IFLAG.EQ.3) GO TO 45
      IF(IFLAG.EQ.4) GO TO 40
      IF((IFLAG.EQ.5).AND.(ABSERR.GT.0.0D0))  GO TO 45
   30 STOP
   40 NFE = 0
      IF(MFLAG.EQ.2) GO TO 50
   45 IFLAG = JFLAG
      IF(KFLAG.EQ.3)  MFLAG = IABS(IFLAG)
   50 JFLAG = IFLAG
      KFLAG = 0
      SAVRE = RELERR
      SAVAE = ABSERR
      RER = 2.0D0*EPS + REMIN
      IF(RELERR.GE.RER)  GO TO 55
      RELERR = RER
      IFLAG = 3
      KFLAG = 3
      RETURN
   55 DT = TOUT -T
      IF(MFLAG.EQ.1)  GO TO 60
      IF(INIT.EQ.0)  GO TO 65
      GO TO 80
   60 INIT = 0
      KOP = 0
      A = T
      CALL F(A,Y,YP)
      NFE = 1
      IF(T.NE.TOUT) GO TO 65
      IFLAG = 2
      RETURN
   65 INIT = 1
      H = DABS(DT)
      TOLN = 0.0D0
      DO 70 K = 1,NEQN
c         TOL=RELERR*DABS(Y(K)) + ABSERR
         TOL=RELERR*cDABS(Y(K)) + ABSERR
          IF(TOL.LE.0.0D0)  GO TO 70
         TOLN = TOL
c         YPK = DABS(YP(K))
         YPK = cDABS(YP(K))
         IF(YPK *H**5.GT.TOL)  H=(TOL/YPK)**0.2D0
   70 CONTINUE
      IF(TOLN.LE.0.0D0)  H = 0.0D0
      H=1D-4
      H= DMAX1(H, U26*DMAX1(DABS(T),DABS(DT)))
      JFLAG = ISIGN(2,IFLAG)
   80 H =DSIGN(H,DT)
      IF(DABS(H).GE.2.0D0*DABS(DT))  KOP = KOP+1
      IF(KOP.NE.100) GO TO 85
      KOP = 0
      IFLAG = 7
      RETURN
   85 IF(DABS(DT).GT.U26*DABS(T))  GO TO 95
      DO 90  K = 1,NEQN
   90 Y(K) = Y(K) + DT*YP(K)
      A = TOUT
      CALL F(A,Y,YP)
      NFE = NFE+ 1
      GO TO  300
   95 OUTPUT=.FALSE.
      SCALE= 2D0/RELERR
      AE = SCALE * ABSERR
  100 HFAILD = .FALSE.
      HMIN = U26 * DABS(T)
      DT = TOUT -T
      IF(DABS(DT).GE.2.0D0*DABS(H))  GO TO 200
      IF(DABS(DT).GT.DABS(H))  GO TO 150
      OUTPUT= .TRUE.
      H = DT
      GO TO  200
  150 H = 0.5D0 * DT
  200 IF(NFE.LE.MAXNFE)  GO TO 220
      IFLAG = 4
      KFLAG = 4
      RETURN
  220 CALL cFEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,F1)
      NFE = NFE + 5
      EEOET=0D0
      DO 250  K = 1,NEQN
c         ET = DABS(Y(K))+ DABS(F1(K)) + AE
         ET = cDABS(Y(K))+ cDABS(F1(K)) + AE
         IF(ET.GT.0.0D0)  GO TO 240
               IFLAG =5
               RETURN
c  240    EE = DABS((-2.09D3 * YP(K) + (2.197D4 * F3(K) -
c     1       1.5048D4 * F4(K))) + (2.2528D4 *F2(K) -2.736D4 *F5(K)))
  240    EE = cDABS((-2.09D3 * YP(K) + (2.197D4 * F3(K) -
     1       1.5048D4 * F4(K))) + (2.2528D4 *F2(K) -2.736D4 *F5(K)))
  250    EEOET = DMAX1(EEOET,EE/ET)
      ESTTOL = DABS(H)*EEOET*SCALE/7.524D5
      IF (ESTTOL.LE.1.0D0)  GO TO 260
      HFAILD = .TRUE.
      OUTPUT = .FALSE.
      S = 0.1D0
      IF(ESTTOL.LT.5.9049D4)  S = 0.9D0/ESTTOL**0.2D0
      H = S * H
      IF(DABS(H).GT.HMIN)  GO TO  200
      IFLAG = 6
      KFLAG = 6
      RETURN
  260 T=T+H
      DO 270  K=1,NEQN
  270    Y(K) =F1(K)
      A = T
      CALL F(A,Y,YP)
      NFE = NFE + 1
      S = 5.0D0
      IF(ESTTOL.GT.1.889568D-4) S = 0.9D0/ESTTOL**0.2D0
      IF (HFAILD)  S = DMIN1(S,1.0D0)
      H= DSIGN( DMAX1(S*DABS(H),HMIN),H)
      IF(OUTPUT)  GO TO 300
      IF(IFLAG.GT.0) GO TO 100
      IFLAG = -2
      RETURN
  300 T = TOUT
      IFLAG = 2
      RETURN
      END

c************************************************
C cfehl
c
      SUBROUTINE cFEHL(F,NEQN,Y,T,H,YP,F1,F2,F3,F4,F5,S)
c      REAL*8  Y(NEQN),T,H,YP(NEQN),F1(NEQN),F2(NEQN),
c     1   F3(NEQN),F4(NEQN),F5(NEQN),S(NEQN),CH
      REAL*8  T,H,CH
      complex*16  Y(NEQN),YP(NEQN),F1(NEQN),S(NEQN),
     *            F2(NEQN),F3(NEQN),F4(NEQN),F5(NEQN)
      CH = H/4.0D0
      DO 221  K = 1,NEQN
  221    F5(K) = Y(K) + CH * YP(K)
      CALL F(T+CH,F5,F1)
      CH = 3.0D0 * H/3.2D1
      DO 222  K = 1,NEQN
  222    F5(K) = Y(K) + CH*(YP(K)+3.0D0*F1(K))
      CALL F(T+3.0D0*H/8.0D0,F5,F2)
      CH = H/2.197D3
      DO 223  K =1,NEQN
  223    F5(K) = Y(K)+CH*(1.932D3*YP(K)+(7.296D3*F2(K)-7.2D3*F1(K)))
      CALL F(T+1.2D1*H/1.3D1,F5,F3)
      CH = H/4.104D3
      DO 224  K= 1,NEQN
  224    F5(K)=Y(K)+CH*((8.341D3*YP(K)-8.45D2*F3(K))
     1         +(2.944D4*F2(K) - 3.2832D4*F1(K)))
      CALL F(T+H,F5,F4)
      CH = H/2.052D4
      DO 225 K= 1,NEQN
  225    F1(K)=Y(K)+CH*((-6.08D3*YP(K)+(9.295D3 *F3(K)-
     1         5.643D3*F4(K)))+(4.104D4*F1(K)-2.8352D4 *F2(K)))
      CALL F(T+H/2.0D0,F1,F5)
      CH = H/7.61805D6
      DO 230  K= 1,NEQN
  230    S(K)=Y(K)+CH*((9.0288D5*YP(K)+(3.855735D6 *F3(K) -
     1        1.371249D6*F4(K)))+(3.953664D6*F2(K)+2.7702D5*F5(K)))
      RETURN
      end
c************************************************
C cDRSF20
C
      SUBROUTINE cDRSF20(R,DR,X,R1,R2,m,c,rl,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C,S1,S,RRA,RRB,R1,R2,R3,RL,cc,rll
c      COMPLEX*16 r, dr, WORK(27),Y(4)
      COMPLEX*16 r, dr, Y(4)
c      DIMENSION IWORK(5)
      COMMON /K1/ S, S1, AKSI, AK, kkk, nk, nal
      COMMON /F2/ RLl, Cc, Mm, i12a
      EXTERNAL  cDRFP2
   11 FORMAT(1X,I10,2D10.3,4D20.8)
   31 FORMAT(17H TOLERANCES RESET,2D12.3)
   41 FORMAT(11H MANY STEPS)
   71 FORMAT(12H MUCH OUTPUT)
   81 FORMAT(14H IMPROPER CALL)
          mm = m
          cc = c
          rll = rl
      NEQN = 4
      PL = dreal(RL)
      AC = dreal(C)


c      print *, rL, C
c      print *, 1D0+PL/AC**2, 1D0-PL/AC**2
c      pause

      TKSI=DSQRT((1D0+PL/AC**2)/2D0+DSQRT((1D0-PL/AC**2)**2/4D0+
     *     (M**2-1D0)/AC**2))
      IF(M.EQ.1) TKSI=1D0
       T=TKSI+X/AC*10D0
*       ELERR=1D-5
      ELERR=1D-12
c      ELERR=1D-08
      ABSERR=0D0
      TFINAL=X
      IFLAG=1
      TOUT=T
      CALL cDRFP0(M,C,RL,T,Y)
c  10 CALL cDRKF45(cDRFP2,NEQN,Y,T,TOUT,ELERR,ABSERR,IFLAG,WORK,IWORK)
   10 CALL cDRKF45a(cDRFP2,NEQN,Y,T,TOUT,ELERR,ABSERR,IFLAG)
c      RRA=DCMPLX(Y(1),Y(2))
c      RRB=DCMPLX(Y(3),Y(4))
      RRA = Y(1) + Y(2)
      RRB = Y(3) + Y(4)
      R3=-S1/(C*AKSI*((RRB/RRA+X/AKSI)*R1+R2))
       R=-S1*R3
       DR=-(X/AKSI+RRB/RRA)*R+R1*S1*RRB/RRA
      GO TO (80,20,30,40,50,60,70,80),IFLAG
   20 TOUT=X
      IF(T .GT.TFINAL) GO TO 10
      IER=0
      GO TO  111
   30 WRITE(7,31)  ELERR,ABSERR
      IER=3
      GO TO 10
   40 WRITE(7,41)
      IER = 4
      GO TO 10
   50 ABSERR=1D-9
      WRITE(7,31)  ELERR,ABSERR
      IER=5
      GO TO 10
   60  ELERR=10D0* ELERR
      WRITE(7,31)  ELERR,ABSERR
      IFLAG=2
      IER=6
      GO TO 10
   70 WRITE(7,71)
      IFLAG=2
      IER=7
      GO TO 10
   80 WRITE(7,81)
      IER=8
      GO TO 1000
  111 CONTINUE
 1000 RETURN
      END
c************************************************
C cDRFP0
C
      SUBROUTINE cDRFP0(M,C,RL,KSI,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 RL,C, y, c2, gksi, al, an
      REAL*8 KSI
      DIMENSION Y(4)
      aM2 = M * M
      C2=C**2
      AKSI=KSI**2-1D0
      GKSI=C2-(RL-C2)/AKSI+(1D0-aM2)/AKSI**2
      AL=KSI*(2D0*(1D0-aM2)-(RL-C2)*AKSI)/(2D0*GKSI*AKSI**3)
      AN=cDSQRT(1D0+cdabs(GKSI)+AL**2)
      Y(1)=1D0/AN
      Y(2)=0D0
      Y(3)=-AL/AN
      Y(4)=-DSQRT(cdabs(GKSI))/AN
      RETURN
      END
c************************************************
C cDRFP2
c
      SUBROUTINE cDRFP2(X,Y,YP)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C,RL, Y(4),YP(4), c2, gksi, akk
      COMMON /F2/ RL,C,M, i1
      AKSI=X**2-1D0
      C2=C**2
      GKSI=C2-(RL-C2)/AKSI+(1D0-M**2)/AKSI**2
      AKK=(Y(1)*Y(3)+Y(2)*Y(4)*(1D0-GKSI)/(Y(1)**2+Y(2)**2+Y(3)**2+
     *     Y(4)**2))
       YP(1)=AKK*Y(1)-Y(3)
       YP(2)=AKK*Y(2)-Y(4)
       YP(3)=AKK*Y(3)+GKSI*Y(1)
       YP(4)=AKK*Y(4)+GKSI*Y(2)
      RETURN
      END
c************************************************
C cDRSF212
C
      SUBROUTINE cDRSF212(R,DR,ar,adr,i12,m,N,c,X,rl,AD,IA,BD,IB,IER)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C,S1,S,AD(IA),RL,cc,rll
c      COMPLEX*16 WORK(27),Y(4),BD(IB), ar, adr, r, dr
      COMPLEX*16 Y(4),BD(IB), ar, adr, r, dr
c      DIMENSION IWORK(5)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /F2/ RLl,Cc,Mm,i12a
      EXTERNAL  cDRF222
   11 FORMAT(1X,I10,2D10.3,4D20.8)
   31 FORMAT(17H TOLERANCES RESET,2D12.3)
   41 FORMAT(11H MANY STEPS)
   71 FORMAT(12H MUCH OUTPUT)
   81 FORMAT(14H IMPROPER CALL)

          mm = m
          cc = c
          rll = rl
          i12a = i12

      NEQN=2
      if(i12.eq.2) NEQN=4
       T=0D0
cccc       ELERR=1D-5
       ELERR=1D-15
      ABSERR=0D0
      TFINAL=X
      IFLAG=1
      TOUT=T
      CALL cDRF2K02(M,N,C,AD,IA,BD,IB,i12,Y)
c   10 CALL cDRKF45(cDRF222,NEQN,Y,T,TOUT,ELERR,ABSERR,IFLAG,WORK,IWORK)
   10 CALL cDRKF45a(cDRF222,NEQN,Y,T,TOUT,ELERR,ABSERR,IFLAG)

      if(i12.eq.0) then
      R=Y(1)
      DR=Y(2)
      end if

      if(i12.eq.1) then
      aR=Y(1)
      aDR=Y(2)
      end if

      if(i12.eq.2) then
      R=Y(1)
      DR=Y(2)
      aR=Y(3)
      aDR=Y(4)
      end if

      GO TO (80,20,30,40,50,60,70,80),IFLAG
   20 TOUT=X
      IF(T.LT.TFINAL) GO TO 10
      IER=0
      GO TO  111
   30 WRITE(7,31)  ELERR,ABSERR
      IER=3
      GO TO 10
   40 WRITE(7,41)
      IER = 4
      GO TO 10
   50 ABSERR=1D-9
      WRITE(7,31)  ELERR,ABSERR
      IER=5
      GO TO 10
   60  ELERR=10D0 * ELERR
      WRITE(7,31)  ELERR,ABSERR
      IFLAG=2
      IER=6
      GO TO 10
   70 WRITE(7,71)
      IFLAG=2
      IER=7
      GO TO 10
   80 WRITE(7,81)
      IER=8
      GO TO 1000
  111 CONTINUE
 1000 RETURN
      END

c************************************************
C cDRF222
c
      SUBROUTINE cDRF222(X,Y,YP)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 C, RL, Y(4),YP(4), ar
      COMMON /F2/ RL,C,M,i12
       AKSI=X**2+1D0
       am = m * m
       ar = RL - (C*X)**2 - aM / AKSI

      if(i12.eq.0.or.i12.eq.1) then
       YP(1)=Y(2)
       YP(2)=(-2D0*X*Y(2)+ aR * Y(1))/AKSI
      end if

      if(i12.eq.2) then
       YP(1)=Y(2)
       YP(2)=(-2D0*X*Y(2) + aR * Y(1))/AKSI
       YP(3)=Y(4)
       YP(4)=(-2D0*X*Y(4) + aR * Y(3))/AKSI
      end if

      RETURN
      END

c************************************************
C cDRF2K02
C
      SUBROUTINE cDRF2K02(M,N,C,AD,IA,BD,IB,i12,Y)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 AD,C,R1,S1,r,s,s2,t,r2,ss,ss2, bd, y
      DIMENSION AD(IA),BD(IB),FACT(300),Y(4)
      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /FACT/ FACT
      COMMON /eps3/ eps3
      COMMON /pi/ pi
       eps4 = 1d-40
c       eps4 = 1d-20
      M2=2*M
      IF(MOD((N-M),2)) 1,2,1
    1 JK=1
      GO TO 3
    2 JK=0
    3 CONTINUE
      S2 = s
      DO 4 I=1,IA
      J=2*(I-1)+JK
c      F=1D0
      F=ad(i)
      DO 5 JJ=1,M2
    5 F=F*(J+JJ)
c      T=AD(I)*F
      T=F
      S2=S2+T
      IF(cDABS(T/s2).LT.eps4) GO TO 6
    4 CONTINUE
    6 CONTINUE
      MF=N+M+1+JK
      F=FACT(MF)

c+&
c      IF(MF.LE.169) F=FACT(MF)
c        IF(MF.ge.168) print *, m,n,mf,f,s2
            f=f/s2
c      IF(MF.LE.170) F=FACT(MF)
c      IF(MF.gt.170) stop 1706

c+&      if(i12.eq.1.or.i12.eq.2) then
c      if(i12.eq.1.or.i12.eq.2) then
c----
       jmn = mod(N-m-jk,4)
       if(jmn.eq.0) ss = 1d0
       if(jmn.eq.1) ss = s1
       if(jmn.eq.2) ss = -1d0
       if(jmn.eq.3) ss = -s1
c+&
       if(i12.eq.1.or.i12.eq.2) then
c       jmn2 = dabs(mod(m-n-jk,4))
       jmn2 = iabs(mod(m-n-jk,4))
       if(jmn2.eq.0) ss2 = 1d0
       if(jmn2.eq.1) ss2 = s1
       if(jmn2.eq.2) ss2 = -1d0
       if(jmn2.eq.3) ss2 = -s1
c-----
      if(mf.le.169)
     * R=Ss*(M2-1D0)*FACT(M+1)*PI/2D0**(2*N-M+1)*C**(M-1-JK)/
     *  FACT(M2+1)*
     *  (F/FACT((N-M-JK)/2+1)/FACT((N+M+JK)/2+1))**2
     *  *s2/BD(M+1)
           if(mf.gt.169) then
       R=Ss*(M2-1D0)*FACT(M+1)*PI/2D0**(2*N-M+1)*C**(M-1-JK)/
     *  FACT(M2+1)*
     *  (F/FACT((N-M-JK)/2+1)/FACT((N+M+JK)/2+1)*1d300)
     *  *s2/BD(M+1)
c          print *,r
        r=r*(F/FACT((N-M-JK)/2+1)/FACT((N+M+JK)/2+1)*1d300)
       end if
c-----------------------------------------------------------------
      R1=Ss2*S2/2D0**M/C**(M+1+JK)/AD(1)/FACT(M+1)
c        if(m.eq.64.and.n.gt.207)  write(*,*) n,mf,
c     *  FACT(M2+1),BD(M+1),s2, r
      end if

      if(i12.eq.0.or.i12.eq.2) then
          R2=Ss*2D0**M*C**(M+JK)*AD(1)*FACT(M+1)/s2
      end if

      if(i12.eq.0.or.i12.eq.2) then
          R2=Ss*2D0**M*C**(M+JK)*AD(1)*FACT(M+1)/s2
      end if

      IF(JK.EQ.1) GO TO 10

      if(i12.eq.0) then
      Y(1)=R2/(m2+1d0)
      Y(2)=0D0
      end if

      if(i12.eq.1) then
      Y(1)=R
      Y(2)=R1*(M2+1D0)
      end if

      if(i12.eq.2) then
      Y(1)=R2/(m2+1d0)
      Y(2)=0D0
      Y(3)=R
      Y(4)=R1*(M2+1D0)
      end if

      GO TO 900

   10 CONTINUE

      if(i12.eq.0) then
      Y(1)=0D0
      Y(2)=R2/(M2+3D0)
      end if

      if(i12.eq.1) then
      Y(1)=R1*(M2+3D0)
      Y(2)=R*(M2-3D0)
      end if

      if(i12.eq.2) then
      Y(1)=0D0
      Y(2)=R2/(M2+3D0)
      Y(3)=R1*(M2+3D0)
      Y(4)=R*(M2-3D0)
      end if

  900 CONTINUE
      RETURN
      END
c************************************************
C Cmte
C   TE mode, core-mantle spheroid
C................................................
C   Matrix, Eq.(34)
C ⁄ƒƒƒ                                                                            ƒƒƒƒƒø
C ≥ ksi1*rw1             ksi1*rw1             ksi1*rw1             ksi1*rw1            ≥
C ≥ G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3  G*rw1+(e2/e1-1)*rw3 ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ≥                                                                                    ≥
C ¿ƒƒƒƒ                                                                           ƒƒƒƒƒŸ
C................................................

      SUBROUTINE CMte(M,C1,KSI1,ksi2,NN)
      parameter (nterms=330, nalphas=19)
      IMPLICIT REAL*8 (A-H,O-Q,T-Z),COMPLEX*16(R-S)
      REAL*8 ksi1, ksi2
      COMPLEX*16 p1, p3, DEL21, del23, GAM21, GAM22, akap22,
     *           AAA(2*nterms,2*nterms),BBB(nalphas,2*nterms),
     *           CCC(nalphas,2*nterms), BCOF, VV2,

     *          aa3(nterms,nterms), bb3(nterms,nterms),
     *          aa1(nterms,nterms), bb1(nterms,nterms),
     *          ag3(nterms,nterms), bg3(nterms,nterms),
     *          ag1(nterms,nterms), bg1(nterms,nterms),
     *          q21(nterms,nterms), q22(nterms,nterms)

      DIMENSION
     *          Ialpha(nalphas),Jalpha(nalphas),
c     *          Ainv(200,200),
     *          rfi(nalphas,2*nterms)

      COMMON /Kc2/  bCOF(nterms,4*nterms), VV2(nterms)
      COMMON /Kfun/ ar1_11(2*nterms), r3_11(2*nterms),
     *              aar1_11(nterms), rr3_11(nterms),
     *              r1_21(2*nterms), r3_21(2*nterms),
     *              rr1_21(nterms), rr3_21(nterms),

     *              r1_22(2*nterms), r3_22(2*nterms),
     *              rr1_22(nterms), rr3_22(nterms),
     *              r1_32(2*nterms), rr1_32(nterms),
     *              p1(nterms), p3(nterms)

      common /ires/ Rw1(nterms,nterms), Rw2(nterms,nterms),
     *              Rw5(nterms,nterms), Rw6(nterms,nterms),
     *              RA1(nterms,nterms), RA3(nterms,nterms),
     *              Rg1(nterms,nterms), Rg3(nterms,nterms),
     *              Rf1, Rf2, Rf3

      common /bbb/ bbb, Ra00(nalphas,nterms), Ra0(nalphas,nterms)
     *             ,rfi

      COMMON /K1/ S, S1, AKSI, AK, K, NK, nal
      COMMON /K2/ Ialpha, Jalpha
      COMMON /INT/ DEL21(nterms,nterms), del23(nterms,nterms),
     *             GAM21(nterms,nterms), GAM22(nterms,nterms),
     *             akap22(nterms,nterms), sig22(nterms,nterms)
      COMMON /RES/ Rz0(nalphas,nterms),
     *             Rz1(nalphas,nterms),Rz2(nalphas,nterms)

      NDIM = 2*nterms
      NP2=NN/2

         IF(M.GT.1) GO TO 20
*
*  Axisymmeric part

c**
c** matrix Eq.(20)
c**
      DO  N = 1, NN
      n1 = n + nn
        DO  L = 1, NN
        l1 = l + nn
c--------
         AAA(N ,L ) = ra3(n,l)
       if(n.eq.l) AAA(N,L) = AAA(N,L) -  rr1_21(l)
         AAA(N ,L1) = ra3(n,l)
       if(n.eq.l) AAA(N,L1) = AAA(N,L1) - rr3_21(l)
c--------
         AAA(N1,L ) =  ra1(n,l) * p1(l)
       if(n.eq.l) AAA(N1,L) =  AAA(N1,L) - rr1_22(l) * p1(l)
         AAA(N1,L1) =  ra1(n,l) * p3(l)
       if(n.eq.l) AAA(N1,L1) =  AAA(N1,L1) - rr3_22(l) * p3(l)
c--------
        end do
      end do

       DO 23 IA = 1, NAL

       DO N = 1, NN
         S2=S
          DO i = 1, NN
        S2 = S2 + DEL21(n,i) * (Rr3_11(I) - aar1_11(i)) * (-RA00(ia,i))
          end do
          ccc(ia,N) = S2
          ccc(ia,N+nn) = S
       end do
   23  CONTINUE

          CALL GJ1(AAA,2*nn,NDIM,ccc,nal,nalphas)

* calculation of Z1 (Eq.(22))

      do ia = 1, nal
      do i = 1, nn
         s2 = s
         do jj = 1, nn
           s2 = s2 + DEL21(jj,i) * (ccc(ia,jj) + ccc(ia,jj+nn))
         end do
         rz0(ia,i) = (s2 + RA00(ia,i)) / R3_11(i)
      end do
      end do

c**
*  Non-axisymmeric part
c**
   20     CONTINUE

      CALL Cmlm(NN,rg1,gam22,ag1)
      CALL Cmlm(NN,rg3,gam22,ag3)

          do i=1,nn
          do j=1,nn
          bi = 0d0
          if(i.eq.j) bi = 1d0
          bg1(i,j) = bi + ksi2 * ra1(i,j)
          bg3(i,j) = bi + ksi1 * ra3(i,j)
          end do
          end do

      CALL Cmlm(NN,bg1,gam22,bg1)
      CALL Cmlm(NN,bg3,gam22,bg3)

c--------
      do i = 1, nn
       do j = 1, nn
       bi = 0d0
       if(i.eq.j) bi = 1d0
       aa1(i,j) = ksi2 * (bi + ksi2 * ra1(i,j)) - ak * ag1(i,j)
       aa3(i,j) = ksi1 * (bi + ksi1 * ra3(i,j)) - ak * ag3(i,j)
       bb1(i,j) = ksi2 * rg1(i,j) - bg1(i,j)
       bb3(i,j) = ksi1 * rg3(i,j) - bg3(i,j)
       end do
      end do

c--------
c  Inverse matrices
c--------
* Q21 -- q(c_2,ksi1)
c--------
      CALL Cmlm(NN,gam22,gam22,q22)
      DO  N = 1, nn
        DO  L = 1, nn
         q21(n,l) = - ak * q22(n,l)
         if(n.eq.l) q21(n,l) =  q21(n,l) + ksi1**2
c         print *, n,l, q22(n,l), q21(n,l), gam22(n,l)
c         pause
        end do
      end do
c      CALL invm01(q21,NN,nterms)
      CALL cinvm01(q21,NN,nterms)

* Q22 -- q(c_2,ksi2)
c--------
      DO  N = 1, nn
        DO  L = 1, nn
         q22(n,l) = - ak * q22(n,l)
         if(n.eq.l) q22(n,l) =  q22(n,l) + ksi2**2
        end do
      end do
c      CALL invm01(q22,NN,nterms)
      CALL cinvm01(q22,NN,nterms)
c--------
c--------

      CALL Cmlm(NN,aa1,q22,aa1)
      CALL Cmlm(NN,aa3,q21,aa3)
      CALL Cmlm(NN,bb1,q22,bb1)
      CALL Cmlm(NN,bb3,q21,bb3)

      do i = 1, nn
       do j = 1, nn
        bb1(i,j) = bb1(i,j) +  akap22(i,j) / (ksi2**2 - ak)
        bb3(i,j) = bb3(i,j) +  akap22(i,j) / (ksi1**2 - ak)
       end do
      end do

c--------
      CALL Cmlm(NN,aa3,gam22,ag3)
      CALL Cmlm(NN,bb3,gam22,bg3)
      CALL Cmlm(NN,aa1,gam22,ag1)
      CALL Cmlm(NN,bb1,gam22,bg1)

c--------

c**
c** matrix Eq.(34)
c**
      DO 5 IV=1,2
      IG=IV-1
      DO 8 N=1,NP2
      N1=N  + NP2
      N2=N1 + NP2
      N3=N2 + NP2
      ND=2*(N-1)+IV
      NG=2*N-IG
      DO 9 L=1,NP2
      L1=L + NP2
      L2=L1+ NP2
      L3=L2+ NP2
      LD=2*(L-1)+IV
      LG=2*L-IG
c--------
      AAA(N ,L ) = ksi1 * (rw1(nd,ld) + rf1 * aa3(nd,ld))
      AAA(N ,L1) = ak * (rg3(nd,lg) - gam22(nd,lg) * rr1_21(lg)
     *             + rf1 * ag3(nd,lg))
      AAA(N ,L2) = ksi1 * (rw2(nd,ld) + rf1 * aa3(nd,ld))
      AAA(N ,L3) = ak * (rg3(nd,lg) - gam22(nd,lg) * rr3_21(lg)
     *             + rf1 * ag3(nd,lg))
c--------
      AAA(N1,L ) = (rg3(ng,ld) - gam22(ng,ld) * rr1_21(ld))
     *             + ksi1 * rf1 * bb3(ng,ld)
      AAA(N1,L1) = ksi1 * rw1(ng,lg) + ak * rf1 * bg3(ng,lg)
      AAA(N1,L2) = (rg3(ng,ld) - gam22(ng,ld) * rr3_21(ld))
     *             + ksi1 * rf1 * bb3(ng,ld)
      AAA(N1,L3) = ksi1 * rw2(ng,lg) + ak * rf1 * bg3(ng,lg)
c--------
      AAA(N2,L ) = ksi2 * (rw5(nd,ld) + rf2 * aa1(nd,ld)) * p1(ld)
      AAA(N2,L1) = ak * ((rg1(nd,lg) - gam22(nd,lg) * rr1_22(lg)) +
     *             rf2 * ag1(nd,lg)) * p1(lg)
      AAA(N2,L2) = ksi2 * (rw6(nd,ld) + rf2 * aa1(nd,ld)) * p3(ld)
      AAA(N2,L3) = ak * ((rg1(nd,lg) - gam22(nd,lg) * rr3_22(lg)) +
     *             rf2 * ag1(nd,lg)) * p3(lg)
c--------
      AAA(N3,L ) = ((rg1(ng,ld) - gam22(ng,ld) * rr1_22(ld))
     *             + ksi2 * rf2 * bb1(ng,ld)) * p1(ld)
      AAA(N3,L1) = (ksi2 * rw5(ng,lg) + ak * rf2 * bg1(ng,lg)) * p1(lg)
      AAA(N3,L2) = ((rg1(ng,ld) - gam22(ng,ld) * rr3_22(ld))
     *             + ksi2 * rf2 * bb1(ng,ld)) * p3(ld)
      AAA(N3,L3) = (ksi2 * rw6(ng,lg) + ak * rf2 * bg1(ng,lg)) * p3(lg)
c--------
    9 CONTINUE
    8 CONTINUE

      DO 11 IA = 1, nal

      DO 13 N=1,NP2
        ND=2*(N-1)+IV
        NG=2*N-IG
        CCC(ia,N    )   = BBB(ia,ND)
        CCC(ia,N+NP2)   = BBB(ia,NG+NN) / (rf1 + 1d0)
        CCC(ia,N+2*NP2) = s
        CCC(ia,N+3*NP2) = s
   13 continue

   11 CONTINUE
c---
          CALL GJ1(AAA,2*nn,NDIM,ccC,nal,nalphas)
c---
      if(iv.eq.1) then
        do ia = 1, nal
         do n = 1, 2*nn
           rfi(ia,n) = ccc(ia,n)
         end do
        end do
      end if

    5 CONTINUE

* calculations of Z1 / Z2

      CALL Cmlm(NN,gam22,q21,ra1)

c         do i = 1, nn
c         do j = 1, nn
c           write(7,*) i,j, gam22(i,j)
c         end do
c         end do
c         do i = 1, nn
c         do j = 1, nn
c           write(7,*) i,j, q21(i,j)
c         end do
c         end do


      do ia = 1, nal
       IF(Jalpha(IA).EQ.1) GO TO 111
       do i = 1, nn
        s2 = s
        s3 = s
        do jj = 1, nn

         bi = 0d0
         if(i.eq.jj) bi = 1d0

         if(MOD(jj,2).eq.0) then
        sx1y1 = ccc(ia,jj/2) + ccc(ia,jj/2+2*np2)
        sx2y2 = rfi(ia,jj/2+np2) + rfi(ia,jj/2+3*np2)

       else

        sx1y1 = rfi(ia,jj/2+1) + rfi(ia,jj/2+1+2*np2)
        sx2y2 = ccc(ia,jj/2+1+np2) + ccc(ia,jj/2+1+3*np2)

         end if

       s2 = s2 + (bi + rf1 * ksi1**2 * q21(i,jj)) * sx1y1
     *            + ak * rf1 * ksi1 * ra1(i,jj) * sx2y2

c          write(*,*) i, (bi + rf1 * ksi1**2 * q21(i,jj)) * sx1y1
c          write(*,*) i, rf1,  ksi1,  ra1(i,jj), sx2y2
c          write(7,*) i, jj,   ra1(i,jj)
c          write(*,*) i, ak * rf1 * ksi1 * ra1(i,jj) * sx2y2
c          write(*,*) s2
c          pause

       s3 = s3 + ((rf1 + 1d0) * bi - rf1 * ksi1**2 * q21(i,jj)) * sx2y2
     *            - rf1 * ksi1 * ra1(i,jj) * sx1y1
        end do
        p1(i) = s2
        p3(i) = s3
       end do

       do i = 1, nn
        s2 = s
        s3 = s
        do jj = 1, nn
        s2 = s2 + DEL21(jj,i) * p1(jj)
c          write(*,*) jj, DEL21(jj,i), p1(jj)
c          pause
        s3 = s3 + DEL21(jj,i) * p3(jj)
        end do
        rz1(ia,i) = (s2-RA0(ia,i)) / R3_11(i)
        rz2(ia,i) = s3 / R3_11(i) / C1
c          write(*,*) i, s2-RA0(ia,i)
c          write(*,*) R3_11(i)
c          write(*,*) i, rz1(1,i)
c
c          write(7,*) i, s2-RA0(ia,i)
c          write(7,*) R3_11(i)
c          write(7,*) i, rz1(1,i)
c          pause
       end do
 111   continue
      end do

c-------------------------------------------------

      RETURN
      END
C    1------------------------------------------------------------------
C
c      SUBROUTINE INVM01(a,n,np)
      SUBROUTINE CINVM01(a,n,np)
C
C **** complex matrix!
C     Computes the inverse of a matrix using the Gauss-Jordan
C     elimination procedure with full pivoting.
C     A is an input matrix of N by N elements, stored in an array
C     of physical dimensions NP by NP.
C     On output, A is replaced by its inverse matrix.
C     The integer arrays IPIV, INDXR, and INDXC are used for
C     bookkeeping on the pivoting.
C
      INTEGER n,np,Nterms
c      DOUBLE PRECISION a(np,np)
      complex*16 a(np,np)
      PARAMETER (nterms=330)
      INTEGER i,icol,irow,j,k,l,ll,indxc(nterms),indxr(nterms),
     *        ipiv(nterms)
c      DOUBLE PRECISION big,dum,pivinv, d0,d1
      complex*16   big,dum,pivinv, d0,d1
c * Draine!
      SAVE d0,d1
      DATA d0,d1/0.D+0,1.D+0/
c      d0 = (0d0,0d0)
c      d1 = (1d0,0d0)
C
      do 11 j=1,n
        ipiv(j)=0
11    continue
C
C     - The main loop over the columns to be reduced.
C
      do 22 i=1,n
        big=d0
C     - - - The outer loop of the search for a pivot element.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (cdabs(a(j,k)).ge.cdabs(big))then
                  big=cdabs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj. ipiv(k) > 1'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
C
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.d0) pause 'a(icol,icol) = 0'
        pivinv=d1/a(icol,icol)
        a(icol,icol)=d1
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=d0
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END
c
c        SUBROUTINE SEGV(M,N,C,KD,CV,EG)
        SUBROUTINE SEGV1(M,jN,rC1,KD,CV,rlc1)
C
C       =========================================================
C       Purpose: Compute the characteristic values of spheroidal
C                wave functions
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                KD --- Function code
C                       KD=1 for Prolate; KD=-1 for Oblate
C       Output:  CV --- Characteristic value for given m, n and c
C                EG(L) --- Characteristic value for mode m and n'
C                          ( L = n' - m + 1 )
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION B(1000),H(1000),D(3000),E(3000),F(3000),CV0(1000),
     &            A(3000),G(3000),EG(2000)
        complex*16 rc1, rlc1(2000)
c        DIMENSION B(100),H(100),D(300),E(300),F(300),CV0(100),
c     &            A(300),G(300),EG(200)
        c = rc1
        n = jn + m - 1

c            write(*,*) m,jn, n,kd, c
c            pause

        IF (C.LT.1.0D-10) THEN
           DO 5 I=1,N
5             EG(I)=(I+M)*(I+M-1.0D0)
           GO TO 70
        ENDIF                                           
        ICM=(N-M+2)/2
        NM=10+INT(0.5*(N-M)+C)
        CS=C*C*KD
        DO 60 L=0,1
           DO 10 I=1,NM
              IF (L.EQ.0) K=2*(I-1)
              IF (L.EQ.1) K=2*I-1
              DK0=M+K
              DK1=M+K+1
              DK2=2*(M+K)
              D2K=2*M+K
              A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
              D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &             *(DK2+3.0))*CS
10            G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
           DO 15 K=2,NM
              E(K)=DSQRT(A(K-1)*G(K))
15            F(K)=E(K)*E(K)
           F(1)=0.0D0
           E(1)=0.0D0
           XA=D(NM)+DABS(E(NM))
           XB=D(NM)-DABS(E(NM))
           NM1=NM-1
           DO 20 I=1,NM1
              T=DABS(E(I))+DABS(E(I+1))
              T1=D(I)+T
              IF (XA.LT.T1) XA=T1
              T1=D(I)-T
              IF (T1.LT.XB) XB=T1
20         CONTINUE
           DO 25 I=1,ICM
              B(I)=XA
25            H(I)=XB
           DO 55 K=1,ICM
              DO 30 K1=K,ICM
                 IF (B(K1).LT.B(K)) THEN
                    B(K)=B(K1)
                    GO TO 35
                 ENDIF
30            CONTINUE
c35            IF (K.NE.1.AND.H(K).LT.H(K-1)) H(K)=H(K-1)
35            IF (K.NE.1) then
                 IF(H(K).LT.H(K-1)) H(K)=H(K-1)
              ENDIF
40            X1=(B(K)+H(K))/2.0D0
              CV0(K)=X1
c+&              IF (DABS((B(K)-H(K))/X1).LT.1.0D-14) GO TO 50
              IF (DABS((B(K)-H(K))/X1).LT.1.0D-15) GO TO 50
              J=0
              S=1.0D0
              DO 45 I=1,NM
                 IF (S.EQ.0.0D0) S=S+1.0D-30
                 T=F(I)/S
                 S=D(I)-T-X1
                 IF (S.LT.0.0D0) J=J+1
45            CONTINUE
              IF (J.LT.K) THEN
                 H(K)=X1
              ELSE
                 B(K)=X1
                 IF (J.GE.ICM) THEN
                    B(ICM)=X1
                 ELSE
                    IF (H(J+1).LT.X1) H(J+1)=X1
                    IF (X1.LT.B(J)) B(J)=X1
                 ENDIF
              ENDIF
              GO TO 40
50            CV0(K)=X1
              IF (L.EQ.0) EG(2*K-1)=CV0(K)
              IF (L.EQ.1) EG(2*K)=CV0(K)
55         CONTINUE
60      CONTINUE
70      CONTINUE
c70      CV=EG(N-M+1)
        CV=EG(N-M+1)
        rlc1(n-m+1) = CV

c            write(*,*) m,n, cv
c            write(*,*) n-m+1
c            write(*,*)  rlc1(n-m+1)
c            pause

        RETURN
        END
c
c      subroutine CDLAM_m(RLC2,M,NN,Cin,EPS,IEr)
      subroutine CDLAM_m(RLC2,M,NN,Cin,IEr)
c complex 'lambda' according Martin....
c vector of eigenvalues found from diagonalizing real symmetric tri-diag
c matrix(tql1).finds egv(m,m),egv(m+1,m)... as xlam(1),xlam(2)... out
c and in -side spheroid.accuracy is improved and d(r,mn) made by bowkmp
      parameter (nterms=330)
      implicit real*8(a-h,o-z)
c ** s1
      complex*16 cin, cdiag,csubd,cc,ss,rs1,rlc2(nterms)
      COMMON /K1/ sS, rS1, AKSI, AK, K, NK, nal
      dimension cdiag(nterms),csubd(nterms)
      integer r,s
      maxdim=nterms
c      maxdim=100
      msq=m**2
      m2=2*m
      sgn=ak
c      nl=kmax1
      nl=nn
      numax=nl-1
cm* now do inside
      cc=cin
c      print *, M,NN, cin
c      pause
c find eigenvalues of even parity first
      s=0d0
c determine matrix size and check for overflow and other errors
      n=numax+3*idint(dsqrt(cdabs(cc)))+2
      if((nl.le.0).or.(n.gt.maxdim)) goto 100
c set the diagonal elements.note that we write r=2*q+s
 1005 r=s
      do 1001 i=1,n
         xa=dfloat((m+r)*(m+r+1))
         cdiag(i)=(2.0d0*xa-dfloat(2*msq+1))*sgn*cc*cc
     $        /dfloat((m2+2*r-1)*(m2+2*r+3))+xa
 1001 r=r+2
c now set the subdiagonal elements
      r=s+2d0
      do 1002 i=2,n
         csubd(i)=dsqrt(dfloat((m2+r)*(m2+r-1))*dfloat(r*(r-1))
     $        /dfloat((m2+2*r-3)*(m2+2*r+1)))*cc*cc/dfloat(m2+2*r-1)
 1002 r=r+2
      ib=0
      call sf02(n,cdiag,csubd,ib)
cm* this goes to setup for and execution of nag routine f02ajf and sorts
      if(ib.ne.0) goto 102
c sf02 returns eigenvalues of same parity in ascending order,
cm* with respect to the size of the real part (see onaka, asano's letter)
cm? is this completely general?
      j=0
      ls=s+1
      do 1030 i=ls,nl,2
         rlc2(i)=cdiag(j+1)
c      write(6,*) i,rlc2(i)
c      print *, i,rlc2(i)
c      pause
 1030 j=j+1
      if(s.eq.1) return
c find eigenvalues of odd parity
      s=1
      goto 1005
cm*
cm* error returns
 100  if(nl.le.0) write(6,600)
      if (n.gt.maxdim) write(6,602) n
c+&
      if (n.gt.maxdim) n=maxdim
c      if (n.gt.maxdim) stop 100
 600  format(20h-nl negative or zero)
 602  format(' -allowable array size exceeded. increase dimension 
     $   statements in sphwf to at least ',i4)
      return
 102  write(6,604) ib
 604  format(20h-f02ajf error number,i4)
      ier=ib
      return
      end
c
      subroutine sf02(n,cdiag,csubd,ifail)
      parameter (nterms=330)
      implicit real*8(a-h,o-z)
      complex*16 cdiag(n),csubd(n),dcmplx
c      dimension ar(100,100),ai(100,100),rr(100),ri(100)
c      dimension intger(100),irank(100)
      dimension ar(nterms,nterms),ai(nterms,nterms),rr(nterms),
     *          ri(nterms),intger(nterms),irank(nterms)
      do 1 i=1,n
         do 1 j=1,n
            ar(i,j)=0.
            ai(i,j)=0.
 1    continue
cm* diagonal
      do 2 i=1,n
         ar(i,i)=dreal(cdiag(i))
         ai(i,i)=dimag(cdiag(i))
 2    continue
cm* subdiagonal, and superdiagonal since symmetric
      do 3 i=2,n
         arr=dreal(csubd(i))
         brr=dimag(csubd(i))
         j=i-1
         ar(i,j) = arr
         ar(j,i) = arr
         ai(i,j) = brr
         ai(j,i) = brr
 3    continue
      ifail = 1
      call f02ajf(ar,nterms,ai,nterms,n,rr,ri,intger,ifail)
      if(ifail.eq.0) go to 20
      write(6,*) ' f02ajf',ifail
      stop
cm* rank the real parts, then rearrange, using mark 12 nag routines
c
c 20   continue
c      ifail=0
c      call m01daf(rr,1,n,'a',irank,ifail)
ccm* rank the real parts rr using a mark 1 nag routine
  20   ifail=1
       call m01aaf(rr,1,n,irank,intger,ifail)
      if(ifail.eq.0) go to 21
c      write(6,*) ' m01daf ',ifail
      write(6,*) ' m01aaf ',ifail
      stop
cm* sort into complex number
c in sub. sf02
c!!! ordering - here! cdiag(irank(i))-YES,cdiag(i)-NO
 21   do 9 i=1,n
c         cdiag(irank(i)) =dcmplx(rr(i),ri(i))
c          print *, i,irank(i)
c          pause
         cdiag(i) =dcmplx(rr(i),ri(i))
 9    continue
c
c      if(.not.wanted) return
c      write(6,99994)
c99994 format(1x,42hcomplex eigenvalues - even then odd parity)
c      write(6,99995) (cdiag(i),i=1,n)
c      do i=1,n
c        write(6,*) i, cdiag(i)
c      end do
c99995 format(1x,3(2e12.4,2x))
      return
      end
c
      SUBROUTINE F02AJF(AR, IAR, AI, IAI, N, RR, RI, INTGER, LFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C
C     EIGENVALUES OF A COMPLEX*16 UNSYMMETRIC MATRIX
C     1ST AUGUST 1971
C
      INTEGER P01AAF, ISAVE, LFAIL, N, IAR, IAI, INTGER(N)
C$P 1
c * Draine!
      CHARACTER*8 SRNAME
      SAVE SRNAME
c      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION ACC, AR(IAR,N), AI(IAI,N),	RR(N), RI(N),
     * X02AAF
      DATA SRNAME /8H F02AJF /
      ISAVE = LFAIL
      LFAIL = 1
      ACC = X02AAF(1.0D0)
      CALL F01AMF(N, 1,	N, AR, IAR, AI,	IAI, INTGER)
      CALL F02ANF(N, ACC, AR, IAR, AI, IAI, RR,	RI, LFAIL)
      IF (LFAIL.NE.0) LFAIL = P01AAF(ISAVE,LFAIL,SRNAME)
      RETURN
      END
c
      SUBROUTINE M01AAF(A, MM, NN, IP, IST, IFAIL)
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C     MARK 10B REVISED.	IER-409	(JAN 1983).
C     DETACHED KEY SORT, TO SORT AN INDEX ARRAY	TO GIVE	THE
C     ASCENDING	ORDER OF A REAL	ARRAY.
      INTEGER P01AAF, M, MM, N,	NN, I, IFAIL, LIMIT, JS, L, J, M2,K,
     * IP(NN), IST(NN)
C$P 1
c * Draine!
      CHARACTER*8 SRNAME
      SAVE SRNAME
c      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION A(NN)
      DATA SRNAME /8H M01AAF /
      M	= MM
      N	= NN
      IF (MM) 440, 440,	20
   20 IF (NN) 420, 420,	40
   40 IF (NN-MM) 460, 400, 60
   60 DO 80 I=M,N
	 IP(I) = I
	 IST(I)	= I
   80 CONTINUE
      IFAIL = 0
  100 LIMIT = (N+M-1)/2
      JS = M - 2
      L	= M - 1
  120 L	= L + 1
      JS = JS +	2
      I	= IST(JS)
      J	= IST(JS+1)
  140 IF (A(I).GT.A(J))	GO TO 160
      M2 = I
      I	= IP(I)
      GO TO 180
  160 M2 = J
      J	= IP(J)
  180 K	= M2
  200 IF (A(I).GT.A(J))	GO TO 220
      IP(K) = I
      K	= I
      I	= IP(I)
      IF (IP(K)-K) 200,	260, 200
  220 IP(K) = J
      K	= J
      J	= IP(J)
      IF (IP(K)-K) 200,	240, 200
  240 J	= I
  260 IP(K) = J
      K	= J
      IF (IP(J)-J) 280,	300, 280
  280 J	= IP(J)
      GO TO 260
  300 IST(L) = M2
      IF (L.LT.LIMIT) GO TO 120
      IF (JS+1-N) 320, 340, 340
  320 JS = JS +	2
      I	= IST(L)
      J	= IST(JS)
      GO TO 140
  340 N	= LIMIT
      IF (LIMIT.GT.M) GO TO 100
      I	= IST(M)
      J	= 0
  360 IF (IP(I).EQ.I) GO TO 380
      J	= J + 1
      K	= IP(I)
      IP(I) = J
      I	= K
      GO TO 360
  380 IP(I) = J	+ 1
      GO TO 480
  400 IP(MM) = 1
      IFAIL = 0
      GO TO 480
  420 IFAIL = P01AAF(IFAIL,3,SRNAME)
      GO TO 480
  440 IFAIL = P01AAF(IFAIL,2,SRNAME)
      GO TO 480
  460 IFAIL = P01AAF(IFAIL,1,SRNAME)
  480 RETURN
      END
c
      SUBROUTINE F01AMZ(AR, IAR, AI, IAI, M, N,	BR, IBR, BI, IBI,CR, CI)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C
C     COMPUTES	C = C +	 A*B  (COMPLEX)	WHERE
C     A	IS RECTANGULAR M BY N.
C     C	MUST BE	DISTINCT FROM B.
C     THE ELEMENTS OF B	MAY BE NON-CONSECUTIVE,	WITH OFFSETS IBR
C     AND IBI.
C
C     .. SCALAR	ARGUMENTS ..
      INTEGER IAR, IAI,	IBR, IBI, M, N
C     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION AR(IAR,N), AI(IAI,N), BR(IBR,N),	BI(IBI,N), CR(M)
     *,	CI(M)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I, J
C     ..
C
      DO 60 J=1,N
	 DO 40 I=1,M
	    CR(I) = (CR(I)+AR(I,J)*BR(1,J)) - AI(I,J)*BI(1,J)
	    CI(I) = (CI(I)+AR(I,J)*BI(1,J)) + AI(I,J)*BR(1,J)
   40	 CONTINUE
   60 CONTINUE
      RETURN
      END
c
      SUBROUTINE F02ANF(NN, ACC, HR, IHR, HI, IHI, WR, WI, LFAIL)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 3 REVISED.
C     MARK 4.5 REVISED
C     MARK 9 REVISED. IER-326 (SEP 1981).
C
C     COMLR
C     FINDS THE	EIGENVALUES OF A COMPLEX UPPER-HESSENBERG MATRIX,
C     H, STORED	IN THE ARRAYS HR(N,N) AND HI(N,N), AND
C     STORES THE REAL PARTS IN THE ARRAY WR(N) AND THE COMPLEX
C     PARTS IN THE ARRAY WI(N).	ACC IS THE RELATIVE MACHINE
C     PRECISION. THE SUBROUTINE	FAILS IF ALL EIGENVALUES TAKE MORE
C     THAN 30*N	ITERATIONS.
C     1ST DECEMBER 1971
C
      INTEGER P01AAF, ISAVE, LFAIL, N, NN, ITS,	L, LL, N1, I, J, M,MM,
     * M1, I1, IHR, IHI, ITN
C$P 1
c * Draine!
      CHARACTER*8 SRNAME
      SAVE SRNAME
c      DOUBLE PRECISION SRNAME
      DOUBLE PRECISION TR, TI, ACC, SR,	SI, XR,	XI, YR,	YI, ZR,	ZI, ZXR,
     * ZXI,HR(IHR,NN), HI(IHI,NN), WR(NN), WI(NN)
      DATA SRNAME /8H F02ANF /
      ISAVE = LFAIL
      LFAIL = 0
      N	= NN
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
   20 IF (N.EQ.0) GO TO	480
      ITS = 0
C     LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
   40 L	= N + 1
      IF (N.EQ.1) GO TO	80
      DO 60 LL=2,N
	 L = L - 1
	 IF ((DABS(HR(L,L-1))+DABS(HI(L,L-1))).LE.(ACC*(DABS(HR(L-1,L-1)
     *	 )+DABS(HI(L-1,L-1))+DABS(HR(L,L))+DABS(HI(L,L))))) GO TO100
   60 CONTINUE
   80 L	= 1
  100 IF (L.EQ.N) GO TO	460
      IF (ITN.GT.0) GO TO 120
      LFAIL = P01AAF(ISAVE,1,SRNAME)
      RETURN
C     FORM SHIFT
  120 IF ((ITS.NE.10) .AND. (ITS.NE.20)) GO TO 140
      SR = DABS(HR(N,N-1)) + DABS(HR(N-1,N-2))
      SI = DABS(HI(N,N-1)) + DABS(HI(N-1,N-2))
      GO TO 180
  140 SR = HR(N,N)
      SI = HI(N,N)
      N1 = N - 1
      XR = HR(N1,N)*HR(N,N1) - HI(N1,N)*HI(N,N1)
      XI = HR(N1,N)*HI(N,N1) + HI(N1,N)*HR(N,N1)
      IF ((XR.EQ.0.0D0)	.AND. (XI.EQ.0.0D0)) GO	TO 180
      YR = (HR(N-1,N-1)-SR)/2.0D0
      YI = (HI(N-1,N-1)-SI)/2.0D0
      CALL A02AAF(YR**2-YI**2+XR, 2.0D0*YR*YI+XI, ZR, ZI)
      IF ((YR*ZR+YI*ZI).GE.0.0D0) GO TO	160
      ZR = -ZR
      ZI = -ZI
  160 CALL A02ACF(XR, XI, YR+ZR, YI+ZI,	ZXR, ZXI)
      XR = ZXR
      XI = ZXI
      SR = SR -	XR
      SI = SI -	XI
  180 DO 200 I=1,N
	 HR(I,I) = HR(I,I) - SR
	 HI(I,I) = HI(I,I) - SI
  200 CONTINUE
      TR = TR +	SR
      TI = TI +	SI
      ITS = ITS	+ 1
      ITN = ITN	- 1
      J	= L + 1
C     LOOK FOR TWO CONSECUTIVE SMALL SUB-DIAGONAL ELEMENTS
      XR = DABS(HR(N-1,N-1)) + DABS(HI(N-1,N-1))
      YR = DABS(HR(N,N-1)) + DABS(HI(N,N-1))
      ZR = DABS(HR(N,N)) + DABS(HI(N,N))
      M	= N
      N1 = N - 1
      IF (J.GT.N1) GO TO 240
      DO 220 MM=J,N1
	 M = M - 1
	 YI = YR
	 YR = DABS(HR(M,M-1)) +	DABS(HI(M,M-1))
	 XI = ZR
	 ZR = XR
	 XR = DABS(HR(M-1,M-1))	+ DABS(HI(M-1,M-1))
	 IF (YR.LE.(ACC*ZR/YI*(ZR+XR+XI))) GO TO 260
  220 CONTINUE
  240 M	= L
C     TRIANGULAR DECOMPOSITION H = L*R
  260 M1 = M + 1
      DO 360 I=M1,N
	 XR = HR(I-1,I-1)
	 XI = HI(I-1,I-1)
	 YR = HR(I,I-1)
	 YI = HI(I,I-1)
	 IF ((DABS(XR)+DABS(XI)).GE.(DABS(YR)+DABS(YI))) GO TO 300
C     INTERCHANGE ROWS OF HR AND HI
	 I1 = I	- 1
	 DO 280	J=I1,N
	    ZR = HR(I-1,J)
	    HR(I-1,J) =	HR(I,J)
	    HR(I,J) = ZR
	    ZI = HI(I-1,J)
	    HI(I-1,J) =	HI(I,J)
	    HI(I,J) = ZI
  280	 CONTINUE
	 CALL A02ACF(XR, XI, YR, YI, ZR, ZI)
	 WR(I) = 1.0D0
	 GO TO 320
  300	 CALL A02ACF(YR, YI, XR, XI, ZR, ZI)
	 WR(I) = -1.0D0
  320	 HR(I,I-1) = ZR
	 HI(I,I-1) = ZI
	 DO 340	J=I,N
	    HR(I,J) = HR(I,J) -	ZR*HR(I-1,J) + ZI*HI(I-1,J)
	    HI(I,J) = HI(I,J) -	ZR*HI(I-1,J) - ZI*HR(I-1,J)
  340	 CONTINUE
  360 CONTINUE
C     COMPOSITION R*L =	H
      M1 = M + 1
      DO 440 J=M1,N
	 XR = HR(J,J-1)
	 XI = HI(J,J-1)
	 HR(J,J-1) = 0.0D0
	 HI(J,J-1) = 0.0D0
C     INTERCHANGE COLUMNS OF HR	AND HI,	IF NECESSARY
	 IF (WR(J).LE.0.0D0) GO	TO 400
	 DO 380	I=L,J
	    ZR = HR(I,J-1)
	    HR(I,J-1) =	HR(I,J)
	    HR(I,J) = ZR
	    ZI = HI(I,J-1)
	    HI(I,J-1) =	HI(I,J)
	    HI(I,J) = ZI
  380	 CONTINUE
  400	 DO 420	I=L,J
	    HR(I,J-1) =	HR(I,J-1) + XR*HR(I,J) - XI*HI(I,J)
	    HI(I,J-1) =	HI(I,J-1) + XR*HI(I,J) + XI*HR(I,J)
  420	 CONTINUE
  440 CONTINUE
      GO TO 40
C     A	ROOT FOUND
  460 WR(N) = HR(N,N) +	TR
      WI(N) = HI(N,N) +	TI
      N	= N - 1
      GO TO 20
  480 RETURN
      END
c
      SUBROUTINE A02ACF(XXR, XXI, YYR, YYI, ZR,	ZI)
C     MARK 2A RELEASE.	NAG COPYRIGHT 1973
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C
C     DIVIDES ONE COMPLEX NUMBER BY A SECOND
C
      DOUBLE PRECISION ONE, YYR, YYI, H, A, ZR,	XXR, XXI, ZI
c   * added (Draine!)
      SAVE ONE
      DATA ONE /1.0D0/
C
      IF (DABS(YYR).LE.DABS(YYI)) GO TO	20
      H	= YYI/YYR
      A	= ONE/(H*YYI+YYR)
      ZR = (XXR+H*XXI)*A
      ZI = (XXI-H*XXR)*A
      RETURN
   20 H	= YYR/YYI
      A	= ONE/(H*YYR+YYI)
      ZR = (H*XXR+XXI)*A
      ZI = (H*XXI-XXR)*A
      RETURN
      END
c
      SUBROUTINE F01AMF(N, K, L, AR, IAR, AI, IAI, INTGER)
C     MARK 2 RELEASE. NAG COPYRIGHT 1972
C     MARK 4.5 REVISED
C     MARK 11 REVISED. VECTORISATION (JAN 1984).
C
C     COMHES
C     GIVEN THE	COMPLEX	UNSYMMETRIC MATRIX,A, STORED IN	THE ARRAYS
C     AR(N,N) AND AI(N,N), THIS	SUBROUTINE REDUCES THE
C     SUB-MATRIX OF ORDER L - K	+ 1, WHICH STARTS AT THE ELEMENTS
C     AR(K,K) AND AI(K,K) AND FINISHES AT THE ELEMENTS AR(L,L)
C     AND AI(L,L), TO HESSENBERG FORM, H, BY NON-ORTHOGONAL
C     ELEMENTARY TRANSFORMATIONS. THE MATRIX H IS OVERWRITTEN ON A
C     WITH DETAILS OF THE TRANSFORMATIONS STORED IN THE	REMAINING
C     TRIANGLE UNDER H AND IN ELEMENTS K TO L OF THE ARRAY
C     INTGER(N).
C     1ST AUGUST 1971
C
      INTEGER LA, L, K1, K, M, I, J, M1, N, IAR, IAI, INTGER(N), MM
      DOUBLE PRECISION XR, XI, YR, YI, AR(IAR,N), AI(IAI,N)
      LA = L - 1
      K1 = K + 1
      IF (K1.GT.LA) RETURN
      DO 160 M=K1,N
	 I = M
	 XR = 0.0D0
	 XI = 0.0D0
	 IF (M.GT.L) GO	TO 120
	 DO 20 J=M,L
	    IF ((DABS(AR(J,M-1))+DABS(AI(J,M-1))).LE.(DABS(XR)+DABS(XI))
     *	    ) GO TO 20
	    XR = AR(J,M-1)
	    XI = AI(J,M-1)
	    I =	J
   20	 CONTINUE
	 INTGER(M) = I
	 IF (I.EQ.M) GO	TO 80
C     INTERCHANGE ROWS AND COLUMNS OF ARRAYS AR	AND AI.
	 DO 40 J=K,N
	    YR = AR(I,J)
	    AR(I,J) = AR(M,J)
	    AR(M,J) = YR
	    YI = AI(I,J)
	    AI(I,J) = AI(M,J)
	    AI(M,J) = YI
   40	 CONTINUE
	 DO 60 J=1,L
	    YR = AR(J,I)
	    AR(J,I) = AR(J,M)
	    AR(J,M) = YR
	    YI = AI(J,I)
	    AI(J,I) = AI(J,M)
	    AI(J,M) = YI
   60	 CONTINUE
   80	 IF (XR.EQ.0.0D0 .AND. XI.EQ.0.0D0) GO TO 120
	 IF (M.GE.L) GO	TO 120
	 M1 = M	+ 1
	 DO 100	I=M1,L
	    YR = AR(I,M-1)
	    YI = AI(I,M-1)
	    CALL A02ACF(YR, YI,	XR, XI,	AR(I,M-1), AI(I,M-1))
  100	 CONTINUE
	 CALL F01AMZ(AR(1,M+1),	IAR, AI(1,M+1),	IAI, L,	L-M,AR(M+1,M-1),
     *	  1, AI(M+1,M-1), 1, AR(1,M), AI(1,M))
  120	 CALL F01AMY(AR(K+1,K),	IAR, AI(K+1,K),	IAI, L-K, M-K,AR(K+1,M),
     *	  AI(K+1,M))
  160 CONTINUE
      DO 200 MM=K1,LA
	 M = LA	- MM + K1
	 M1 = M	- 2
	 I = INTGER(M)
	 IF (I.EQ.M) GO	TO 200
	 IF (K.GT.M1) GO TO 200
	 DO 180	J=K,M1
	    YR = AR(I,J)
	    AR(I,J) = AR(M,J)
	    AR(M,J) = YR
	    YI = AI(I,J)
	    AI(I,J) = AI(M,J)
	    AI(M,J) = YI
  180	 CONTINUE
  200 CONTINUE
      RETURN
      END
c
      SUBROUTINE F01AMY(AR, IAR, AI, IAI, M, N,	BR, BI)
C     MARK 11 RELEASE. NAG COPYRIGHT 1983.
C
C     COMPUTES	B = L**(-1) * B	 (COMPLEX) WHERE
C     L	IS A LOWER TRIANGULAR MATRIX OF	THE SPECIAL FORM
C     (	L11  0 )
C     (	L21  I )
C     L11 IS LOWER TRIANGULAR OF ORDER N,
C     L21 IS RECTANGULAR (M-N) BY N.
C     A	HOLDS THE SUBDIAGONAL ELEMENTS OF THE FIRST N COLUMNS OF L
C     AND THE DIAGONAL ELEMENTS	OF L ARE TAKEN TO BE 1.0 .
C
C     .. SCALAR	ARGUMENTS ..
      INTEGER IAR, IAI,	M, N
C     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION AR(IAR,N), AI(IAI,N), BR(M), BI(M)
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE PRECISION XR, XI
      INTEGER I, J, JP1
C     ..
C
      DO 40 J=1,N
	 JP1 = J + 1
	 IF (JP1.GT.M) GO TO 40
	 XR = BR(J)
	 XI = BI(J)
	 DO 20 I=JP1,M
	    BR(I) = (BR(I)-AR(I,J)*XR) + AI(I,J)*XI
	    BI(I) = (BI(I)-AR(I,J)*XI) - AI(I,J)*XR
   20	 CONTINUE
   40 CONTINUE
      RETURN
      END
c
      INTEGER FUNCTION P01AAF(IFAIL, ERROR, SRNAME)
C     MARK 1 RELEASE.  NAG COPYRIGHT 1971
C     MARK 3 REVISED
C     MARK 4A REVISED, IER-45
C     MARK 4.5 REVISED
C     MARK 7 REVISED (DEC 1978)
C     MARK 11 REVISED (FEB 1984)
C     RETURNS THE VALUE	OF ERROR OR TERMINATES THE PROGRAM.
      INTEGER ERROR, IFAIL, NOUT
C$P 1
c * Draine!
      CHARACTER*8 SRNAME
c      DOUBLE PRECISION SRNAME
      CHARACTER REPLAY
C     TEST IF NO ERROR DETECTED
      IF (ERROR.EQ.0) GO TO 20
C     DETERMINE	OUTPUT UNIT FOR	MESSAGE
      CALL X04AAF (0,NOUT)
C     TEST FOR SOFT FAILURE
      IF (MOD(IFAIL,10).EQ.1) GO TO 10
C     HARD FAILURE
****   !!!!!!
****   STDOUT is chanched on '*'
****   V.Grikurov, 15/09/94      
c      WRITE (*,99999) SRNAME, ERROR
c      READ (*, 11111) REPLAY
      WRITE (NOUT,99999) SRNAME, ERROR
      READ (NOUT, 11111) REPLAY
C     ******************** IMPLEMENTATION NOTE ********************
C     THE FOLLOWING STOP STATEMENT MAY BE REPLACED BY A	CALL TO	AN
C     IMPLEMENTATION-DEPENDENT ROUTINE TO DISPLAY A MESSAGE AND/OR
C     ABORT THE	PROGRAM.
C     *************************************************************
      IF(REPLAY.EQ.'Y' .OR. REPLAY.EQ.'y') GOTO 20
      STOP
C     SOFT FAIL
C     TEST IF ERROR MESSAGES SUPPRESSED
   10 IF (MOD(IFAIL/10,10).EQ.0) GO TO 20
      WRITE (NOUT,99999) SRNAME, ERROR
   20 P01AAF = ERROR
      RETURN
99999 FORMAT(1H , 38HERROR DETECTED BY	NAG LIBRARY ROUTINE , A8,
     * 11H - IFAIL = , I5/1H , 22HCONTINUE ANYWAY (Y/N)? )
11111 FORMAT( A1 )
      END
c
      SUBROUTINE A02AAF(XXR, XXI, YR, YI)
C     MARK 2A RELEASE.	NAG COPYRIGHT 1973
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C     COMPUTES THE SQUARE ROOT OF A COMPLEX NUMBER
C
      DOUBLE PRECISION ZERO, HALF, XI, XXI, H, XXR, YR,	YI, A02ABF
c   * added Draine!
      SAVE HALF,ZERO
      DATA ZERO	/0.0D0/, HALF /0.5D0/
C
      XI = XXI
      H	= DSQRT((DABS(XXR)+A02ABF(XXR,XI))*HALF)
      IF (XI.NE.ZERO) XI = XI/(H+H)
      IF (XXR.LT.ZERO) GO TO 20
      YR = H
      YI = XI
      RETURN
   20 IF (XI.LT.ZERO) GO TO 40
      YR = XI
      YI = H
      RETURN
   40 YR = -XI
      YI = -H
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION	X02AAF(X)
C     NAG COPYRIGHT 1975
C     MARK 4.5 RELEASE
      DOUBLE PRECISION X
C     *	EPS *
C
C     IBM DOUBLE PRECISION VERSION
C
C     RETURNS THE VALUE	EPS WHERE EPS IS THE SMALLEST
C     POSITIVE
C     NUMBER SUCH THAT 1.0 + EPS > 1.0
C     THE X PARAMETER IS NOT USED
C     FOR ICL 1900
C     X02AAF = 2.0**(-37.0)
C     FOR IBM 360/370
      X02AAF = 2.0D0**(-52.0D0)
c      DOUBLE PRECISION Z
C     DATA Z/Z'3410000000000000'/
c      DATA Z/Z'0000000000000001'/
c      X02AAF = Z
      RETURN
      END
c
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     IF I = 0,	SETS NERR TO CURRENT ERROR MESSAGE UNIT	NUMBER
C     (STORED IN NERR1).
C     IF I = 1,	CHANGES	CURRENT	ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     *** NOTE ***
C     THIS ROUTINE ASSUMES THAT	THE VALUE OF NERR1 IS SAVED
C     BETWEEN CALLS.  IN SOME IMPLEMENTATIONS IT MAY BE
C     NECESSARY	TO STORE NERR1 IN A LABELLED COMMON
C     BLOCK /AX04AA/ TO	ACHIEVE	THIS.
C
C     .. SCALAR	ARGUMENTS ..
      INTEGER I, NERR
C     ..
C     .. LOCAL SCALARS ..
      INTEGER NERR1
C     ..
c * Draine!
      SAVE NERR1
      DATA NERR1 /6/
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1	= NERR
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION	A02ABF(XXR, XXI)
C     NAG COPYRIGHT 1975
C     MARK 4.5 REVISED
C     MARK 5C REVISED
C
C     RETURNS THE ABSOLUTE VALUE OF A COMPLEX NUMBER VIA ROUTINE
C     NAME
C
      DOUBLE PRECISION ZERO, ONE, XR, XXR, XI, XXI, H
c   * added (Draine!)
      SAVE ONE,ZERO
      DATA ZERO	/0.0D0/, ONE /1.0D0/
C
      XR = DABS(XXR)
      XI = DABS(XXI)
      IF (XI.LE.XR) GO TO 20
      H	= XR
      XR = XI
      XI = H
   20 IF (XI.NE.ZERO) GO TO 40
      A02ABF = XR
      RETURN
   40 H	= XR*DSQRT(ONE+(XI/XR)**2)
      A02ABF = H
      RETURN
      END
