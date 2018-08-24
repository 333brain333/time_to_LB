
      subroutine igrf_sub(xlat,xlong,year,height,
     &          xl,icode,dip,dec)
c----------------------------------------------------------------
c   INPUT:
c       xlat    geodatic latitude in degrees
c       xlong   geodatic longitude in degrees
c       year    decimal year (year+month/12.0-0.5 or year+day-of-year/365 
c               or 366 if leap year) 
c       height  height in km
c   OUTPUT:
c       xl      L value
c       icode   =1  L is correct; =2  L is not correct;
c               =3  an approximation is used
c       dip     geomagnetic inclination in degrees
c       dec     geomagnetic declination in degress
c----------------------------------------------------------------

      INTEGER           EGNR,AGNR,OGNR
      REAL              LATI,LONGI
      COMMON/GENER/     UMR,ERA,AQUAD,BQUAD
C
      CALL INITIZE
      ibbb=0
      ALOG2=ALOG(2.)
      ISTART=1
      lati=xlat
      longi=xlong
c
C----------------CALCULATE PROFILES-----------------------------------
c
      CALL FELDCOF(YEAR,DIMO)
      CALL FELDG(LATI,LONGI,HEIGHT,BNORTH,BEAST,BDOWN,BABS)
      CALL SHELLG(LATI,LONGI,HEIGHT,DIMO,XL,ICODE,BAB1)
      DIP=ASIN(BDOWN/BABS)/UMR
      DEC=ASIN(BEAST/SQRT(BEAST*BEAST+BNORTH*BNORTH))/UMR
      RETURN
      END
c
c
C SHELLIG.FOR, Version 2.0, January 1992
C
C 11/01/91-DKB- SHELLG: lowest starting point for B0 search is 2  
C  1/27/92-DKB- Adopted to IGRF-91 coeffcients model
C  2/05/92-DKB- Reduce variable-names: INTER(P)SHC,EXTRA(P)SHC,INITI(ALI)ZE
C  8/08/95-DKB- Updated to IGRF-45-95; new coeff. DGRF90, IGRF95, IGRF95S
C  5/31/00-DKB- Updated to IGRF-45-00; new coeff.: IGRF00, IGRF00s
C  3/24/05-DKB- Updated to IGRF-45-10; new coeff.: IGRF05, IGRF05s
C
C*********************************************************************
C  SUBROUTINES FINDB0, SHELLG, STOER, FELDG, FELDCOF, GETSHC,        *
C       INTERSHC, EXTRASHC, INITIZE                                  *
C*********************************************************************
C*********************************************************************
C
C
      SUBROUTINE FINDB0(STPS,BDEL,VALUE,BEQU,RR0)
C--------------------------------------------------------------------
C FINDS SMALLEST MAGNETIC FIELD STRENGTH ON FIELD LINE
C
C INPUT:   STPS   STEP SIZE FOR FIELD LINE TRACING
C       COMMON/FIDB0/
C          SP     DIPOLE ORIENTED COORDINATES FORM SHELLG; P(1,*),
C                 P(2,*), P(3,*) CLOSEST TO MAGNETIC EQUATOR 
C          BDEL   REQUIRED ACCURACY  = [ B(LAST) - BEQU ] / BEQU  
C                 B(LAST)  IS FIELD STRENGTH BEFORE BEQU
C
C OUTPUT:  VALUE  =.FALSE., IF BEQU IS NOT MINIMAL VALUE ON FIELD LINE
C          BEQU   MAGNETIC FIELD STRENGTH AT MAGNETIC EQUATOR
C          RR0    EQUATORIAL RADIUS NORMALIZED TO EARTH RADIUS
C          BDEL   FINAL ACHIEVED ACCURACY
C--------------------------------------------------------------------
      DIMENSION         P(8,4),SP(3)
      LOGICAL           VALUE
      COMMON/FIDB0/     SP
C
      STEP=STPS
      IRUN=0
7777  IRUN=IRUN+1
      IF(IRUN.GT.5) THEN
        VALUE=.FALSE.
        GOTO 8888
      ENDIF
C*********************FIRST THREE POINTS 
      P(1,2)=SP(1)
      P(2,2)=SP(2)
      P(3,2)=SP(3)
      STEP=-SIGN(STEP,P(3,2))
      CALL STOER(P(1,2),BQ2,R2)
      P(1,3)=P(1,2)+0.5*STEP*P(4,2)
      P(2,3)=P(2,2)+0.5*STEP*P(5,2)
      P(3,3)=P(3,2)+0.5*STEP
      CALL STOER(P(1,3),BQ3,R3)
      P(1,1)=P(1,2)-STEP*(2.*P(4,2)-P(4,3))
      P(2,1)=P(2,2)-STEP*(2.*P(5,2)-P(5,3))
      P(3,1)=P(3,2)-STEP
      CALL STOER(P(1,1),BQ1,R1)
      P(1,3)=P(1,2)+STEP*(20.*P(4,3)-3.*P(4,2)+P(4,1))/18.
      P(2,3)=P(2,2)+STEP*(20.*P(5,3)-3.*P(5,2)+P(5,1))/18.
      P(3,3)=P(3,2)+STEP
      CALL STOER(P(1,3),BQ3,R3)
C******************INVERT SENSE IF REQUIRED
      IF(BQ3.LE.BQ1) GOTO 2
      STEP=-STEP
      R3=R1
      BQ3=BQ1
      DO 1 I=1,5
           ZZ=P(I,1)
           P(I,1)=P(I,3)
1          P(I,3)=ZZ
C******************INITIALIZATION 
2     STEP12=STEP/12.
      VALUE=.TRUE.
      BMIN=1.E4
      BOLD=1.E4
C******************CORRECTOR (FIELD LINE TRACING)
      N=0
5555  P(1,3)=P(1,2)+STEP12*(5.*P(4,3)+8.*P(4,2)-P(4,1))
      N=N+1
      P(2,3)=P(2,2)+STEP12*(5.*P(5,3)+8.*P(5,2)-P(5,1))
C******************PREDICTOR (FIELD LINE TRACING)
      P(1,4)=P(1,3)+STEP12*(23.*P(4,3)-16.*P(4,2)+5.*P(4,1))
      P(2,4)=P(2,3)+STEP12*(23.*P(5,3)-16.*P(5,2)+5.*P(5,1))
      P(3,4)=P(3,3)+STEP
      CALL STOER(P(1,4),BQ3,R3)
      DO 1111 J=1,3
        DO 1111 I=1,7 
1111       P(I,J)=P(I,J+1)
      B=SQRT(BQ3)
      IF(B.LT.BMIN) BMIN=B
      IF(B.LE.BOLD) THEN
           BOLD=B
           ROLD=1./R3
           SP(1)=P(1,4)
           SP(2)=P(2,4)
           SP(3)=P(3,4)
           GOTO 5555
      ENDIF
      IF(BOLD.NE.BMIN) THEN
           VALUE=.FALSE.
      ENDIF
      BDELTA=(B-BOLD)/BOLD
      IF(BDELTA.GT.BDEL) THEN
           STEP=STEP/10.
           GOTO 7777
      ENDIF
8888  RR0=ROLD
      BEQU=BOLD
      BDEL=BDELTA
      RETURN
      END
C
C
      SUBROUTINE SHELLG(GLAT,GLON,ALT,DIMO,FL,ICODE,B0)
C--------------------------------------------------------------------
C CALCULATES L-VALUE FOR SPECIFIED GEODAETIC COORDINATES, ALTITUDE
C AND GEMAGNETIC FIELD MODEL.
C REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTER, INTERNAL NOTE 
C      NO. 67, 1970.
C      G. KLUGE, COMPUTER PHYSICS COMMUNICATIONS 3, 31-35, 1972
C--------------------------------------------------------------------
C CHANGES (D. BILITZA, NOV 87):
C   - USING CORRECT DIPOL MOMENT I.E.,DIFFERENT COMMON/MODEL/
C   - USING IGRF EARTH MAGNETIC FIELD MODELS FROM 1945 TO 1990
C--------------------------------------------------------------------
C  INPUT:  ENTRY POINT SHELLG
C               GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
C               GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
C               ALT   ALTITUDE IN KM ABOVE SEA LEVEL
C
C          ENTRY POINT SHELLC
C               V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
C                       X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
C                       Y-AXIS POINTING TO EQUATOR AT 90 LONG.
C                       Z-AXIS POINTING TO NORTH POLE
C
C          DIMO       DIPOL MOMENT IN GAUSS (NORMALIZED TO EARTH RADIUS) 
C
C          COMMON 
C               X(3)    NOT USED
C               H(144)  FIELD MODEL COEFFICIENTS ADJUSTED FOR SHELLG
C-----------------------------------------------------------------------
C  OUTPUT: FL           L-VALUE
C          ICODE        =1 NORMAL COMPLETION
C                       =2 UNPHYSICAL CONJUGATE POINT (FL MEANINGLESS)
C                       =3 SHELL PARAMETER GREATER THAN LIMIT UP TO
C                          WHICH ACCURATE CALCULATION IS REQUIRED;
C                          APPROXIMATION IS USED.
C          B0           MAGNETIC FIELD STRENGTH IN GAUSS
C-----------------------------------------------------------------------
      DIMENSION         V(3),U(3,3),P(8,100),SP(3)
      COMMON            X(3),H(144)
      COMMON/FIDB0/     SP
      COMMON/GENER/     UMR,ERA,AQUAD,BQUAD
C
C-- RMIN, RMAX ARE BOUNDARIES FOR IDENTIFICATION OF ICODE=2 AND 3
C-- STEP IS STEP SIZE FOR FIELD LINE TRACING
C-- STEQ IS STEP SIZE FOR INTEGRATION
C 
      DATA RMIN,RMAX    /0.05,1.01/
      DATA STEP,STEQ    /0.20,0.03/
      BEQU=1.E10
C*****ENTRY POINT  SHELLG  TO BE USED WITH GEODETIC CO-ORDINATES
      RLAT=GLAT*UMR
      CT=SIN(RLAT)                                              
      ST=COS(RLAT)                                              
      D=SQRT(AQUAD-(AQUAD-BQUAD)*CT*CT)
      X(1)=(ALT+AQUAD/D)*ST/ERA
      X(3)=(ALT+BQUAD/D)*CT/ERA
      RLON=GLON*UMR
      X(2)=X(1)*SIN(RLON)                                       
      X(1)=X(1)*COS(RLON)                                       
      GOTO9                                                     
      ENTRY SHELLC(V,FL,B0)                                     
C*****ENTRY POINT  SHELLC  TO BE USED WITH CARTESIAN CO-ORDINATES
      X(1)=V(1)                                                  
      X(2)=V(2)                                                  
      X(3)=V(3)                                                  
C*****CONVERT TO DIPOL-ORIENTED CO-ORDINATES                     
      DATA U/                 +0.3511737,-0.9148385,-0.1993679,  
     A                        +0.9335804,+0.3583680,+0.0000000,  
     B                        +0.0714471,-0.1861260,+0.9799247/  
9     RQ=1./(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
      R3H=SQRT(RQ*SQRT(RQ))                                      
      P(1,2)=(X(1)*U(1,1)+X(2)*U(2,1)+X(3)*U(3,1))*R3H           
      P(2,2)=(X(1)*U(1,2)+X(2)*U(2,2)            )*R3H           
      P(3,2)=(X(1)*U(1,3)+X(2)*U(2,3)+X(3)*U(3,3))*RQ            
C*****FIRST THREE POINTS OF FIELD LINE                           
      STEP=-SIGN(STEP,P(3,2))                                    
      CALL STOER(P(1,2),BQ2,R2)                                  
      B0=SQRT(BQ2)                                               
      P(1,3)=P(1,2)+0.5*STEP*P(4,2)                              
      P(2,3)=P(2,2)+0.5*STEP*P(5,2)                              
      P(3,3)=P(3,2)+0.5*STEP                                     
      CALL STOER(P(1,3),BQ3,R3)                                  
      P(1,1)=P(1,2)-STEP*(2.*P(4,2)-P(4,3))                      
      P(2,1)=P(2,2)-STEP*(2.*P(5,2)-P(5,3))                      
      P(3,1)=P(3,2)-STEP                                         
      CALL STOER(P(1,1),BQ1,R1)                                  
      P(1,3)=P(1,2)+STEP*(20.*P(4,3)-3.*P(4,2)+P(4,1))/18.       
      P(2,3)=P(2,2)+STEP*(20.*P(5,3)-3.*P(5,2)+P(5,1))/18.       
      P(3,3)=P(3,2)+STEP                                         
      CALL STOER(P(1,3),BQ3,R3)                                  
C*****INVERT SENSE IF REQUIRED                                   
      IF(BQ3.LE.BQ1)GOTO2                                        
      STEP=-STEP                                                 
      R3=R1                                                      
      BQ3=BQ1                                                    
      DO 1 I=1,7                                                 
      ZZ=P(I,1)                                                  
      P(I,1)=P(I,3)                                              
1     P(I,3)=ZZ                                                  
C*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
2     IF(BQ1.LT.BEQU) THEN
        BEQU=BQ1
        IEQU=1
      ENDIF
      IF(BQ2.LT.BEQU) THEN
        BEQU=BQ2
        IEQU=2
      ENDIF
      IF(BQ3.LT.BEQU) THEN
        BEQU=BQ3
        IEQU=3
      ENDIF
C*****INITIALIZATION OF INTEGRATION LOOPS                        
      STEP12=STEP/12.
      STEP2=STEP+STEP                                            
      STEQ=SIGN(STEQ,STEP)                                       
      FI=0.                                                      
      ICODE=1                                                    
      ORADIK=0.                                                  
      OTERM=0.                                                   
      STP=R2*STEQ                                                
      Z=P(3,2)+STP                                               
      STP=STP/0.75
      P(8,1)=STEP2*(P(1,1)*P(4,1)+P(2,1)*P(5,1))                 
      P(8,2)=STEP2*(P(1,2)*P(4,2)+P(2,2)*P(5,2))
C*****MAIN LOOP (FIELD LINE TRACING)                             
      DO 3 N=3,3333                                              
C*****CORRECTOR (FIELD LINE TRACING)                             
          P(1,N)=P(1,N-1)+STEP12*(5.*P(4,N)+8.*P(4,N-1)-P(4,N-2))    
          P(2,N)=P(2,N-1)+STEP12*(5.*P(5,N)+8.*P(5,N-1)-P(5,N-2))    
C*****PREPARE EXPANSION COEFFICIENTS FOR INTERPOLATION           
C*****OF SLOWLY VARYING QUANTITIES                               
          P(8,N)=STEP2*(P(1,N)*P(4,N)+P(2,N)*P(5,N))                 
          C0=P(1,N-1)**2+P(2,N-1)**2                                 
          C1=P(8,N-1)                                                
          C2=(P(8,N)-P(8,N-2))*0.25                                  
          C3=(P(8,N)+P(8,N-2)-C1-C1)/6.0
          D0=P(6,N-1)                                                
          D1=(P(6,N)-P(6,N-2))*0.5
          D2=(P(6,N)+P(6,N-2)-D0-D0)*0.5                             
          E0=P(7,N-1)
          E1=(P(7,N)-P(7,N-2))*0.5                                   
          E2=(P(7,N)+P(7,N-2)-E0-E0)*0.5                             
C*****INNER LOOP (FOR QUADRATURE)                                
4         T=(Z-P(3,N-1))/STEP                                        
          IF(T.GT.1.)GOTO5                                           
          HLI=0.5*(((C3*T+C2)*T+C1)*T+C0)                            
          ZQ=Z*Z
          R=HLI+SQRT(HLI*HLI+ZQ)
          IF(R.LE.RMIN)GOTO30                               
          RQ=R*R
          FF=SQRT(1.+3.*ZQ/RQ)                              
          RADIK=B0-((D2*T+D1)*T+D0)*R*RQ*FF                 
          IF(R-RMAX)44,44,45                                
45        ICODE=2                                           
          RADIK=RADIK-12.*(R-RMAX)**2                       
44        IF(RADIK+RADIK.LE.ORADIK) GOTO 10
          TERM=SQRT(RADIK)*FF*((E2*T+E1)*T+E0)/(RQ+ZQ)      
          FI=FI+STP*(OTERM+TERM)                            
          ORADIK=RADIK                                      
          OTERM=TERM                                        
          STP=R*STEQ                                        
          Z=Z+STP                                           
          GOTO4                                             
C*****PREDICTOR (FIELD LINE TRACING)                    
5         P(1,N+1)=P(1,N)+STEP12*(23.*P(4,N)-16.*P(4,N-1)+5.*P(4,N-2))  
          P(2,N+1)=P(2,N)+STEP12*(23.*P(5,N)-16.*P(5,N-1)+5.*P(5,N-2))  
          P(3,N+1)=P(3,N)+STEP                                          
          CALL STOER(P(1,N+1),BQ3,R3)                                   
C*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
          IF(BQ3.LT.BEQU) THEN
             IEQU=N+1
             BEQU=BQ3
          ENDIF
3     CONTINUE
10    IF(IEQU.lt.2) IEQU=2 
      SP(1)=P(1,IEQU-1)
      SP(2)=P(2,IEQU-1)
      SP(3)=P(3,IEQU-1)
      IF(ORADIK.LT.1E-15)GOTO11                                     
      FI=FI+STP/0.75*OTERM*ORADIK/(ORADIK-RADIK)              
C
C-- The minimal allowable value of FI was changed from 1E-15 to 1E-12,
C-- because 1E-38 is the minimal allowable arg. for ALOG in our envir.
C-- D. Bilitza, Nov 87.
C
11    FI=0.5*ABS(FI)/SQRT(B0)+1E-12                       
C
C*****COMPUTE L FROM B AND I.  SAME AS CARMEL IN INVAR.  
C
C-- Correct dipole moment is used here. D. Bilitza, Nov 87.
C
      DIMOB0=DIMO/B0
      arg1=alog(FI)
      arg2=alog(DIMOB0)
c       arg = FI*FI*FI/DIMOB0
c       if(abs(arg).gt.88.0) arg=88.0
      XX=3*arg1-arg2
      IF(XX.GT.23.0) GOTO 776   
      IF(XX.GT.11.7) GOTO 775  
      IF(XX.GT.+3.0) GOTO 774    
      IF(XX.GT.-3.0) GOTO 773   
      IF(XX.GT.-22.) GOTO 772  
771   GG=3.33338E-1*XX+3.0062102E-1                                 
      GOTO777                                                          
772   GG=((((((((-8.1537735E-14*XX+8.3232531E-13)*XX+1.0066362E-9)*XX+  
     18.1048663E-8)*XX+3.2916354E-6)*XX+8.2711096E-5)*XX+1.3714667E-3)* 
     2XX+1.5017245E-2)*XX+4.3432642E-1)*XX+6.2337691E-1                 
      GOTO777                                                           
773   GG=((((((((2.6047023E-10*XX+2.3028767E-9)*XX-2.1997983E-8)*XX-    
     15.3977642E-7)*XX-3.3408822E-6)*XX+3.8379917E-5)*XX+1.1784234E-3)* 
     2XX+1.4492441E-2)*XX+4.3352788E-1)*XX+6.228644E-1                  
      GOTO777                                                           
774   GG=((((((((6.3271665E-10*XX-3.958306E-8)*XX+9.9766148E-07)*XX-    
     11.2531932E-5)*XX+7.9451313E-5)*XX-3.2077032E-4)*XX+2.1680398E-3)* 
     2XX+1.2817956E-2)*XX+4.3510529E-1)*XX+6.222355E-1                  
      GOTO777                                                           
775   GG=(((((2.8212095E-8*XX-3.8049276E-6)*XX+2.170224E-4)*XX-6.7310339
     1E-3)*XX+1.2038224E-1)*XX-1.8461796E-1)*XX+2.0007187E0             
      GOTO777                                                           
776   GG=XX-3.0460681E0                                                 
777   FL=EXP(ALOG((1.+EXP(GG))*DIMOB0)/3.0)
      RETURN                                                            
C*****APPROXIMATION FOR HIGH VALUES OF L.                               
30    ICODE=3                                                           
      T=-P(3,N-1)/STEP                                                  
      FL=1./(ABS(((C3*T+C2)*T+C1)*T+C0)+1E-15)                          
      RETURN                                                            
      END                                                               
C
C
      SUBROUTINE STOER(P,BQ,R)                                          
C*******************************************************************
C* SUBROUTINE USED FOR FIELD LINE TRACING IN SHELLG                *
C* CALLS ENTRY POINT FELDI IN GEOMAGNETIC FIELD SUBROUTINE FELDG   *
C*******************************************************************
      DIMENSION         P(7),U(3,3)
      COMMON            XI(3),H(144)
C*****XM,YM,ZM  ARE GEOMAGNETIC CARTESIAN INVERSE CO-ORDINATES          
      ZM=P(3)                                                           
      FLI=P(1)*P(1)+P(2)*P(2)+1E-15
      R=0.5*(FLI+SQRT(FLI*FLI+(ZM+ZM)**2))
      RQ=R*R
      WR=SQRT(R)                                                        
      XM=P(1)*WR                                                        
      YM=P(2)*WR                                                        
C*****TRANSFORM TO GEOGRAPHIC CO-ORDINATE SYSTEM                        
      DATA U/                 +0.3511737,-0.9148385,-0.1993679,         
     A                        +0.9335804,+0.3583680,+0.0000000,         
     B                        +0.0714471,-0.1861260,+0.9799247/         
      XI(1)=XM*U(1,1)+YM*U(1,2)+ZM*U(1,3)                               
      XI(2)=XM*U(2,1)+YM*U(2,2)+ZM*U(2,3)                               
      XI(3)=XM*U(3,1)          +ZM*U(3,3)                               
C*****COMPUTE DERIVATIVES                                               
C      CALL FELDI(XI,H)                                                  
      CALL FELDI                                                  
      Q=H(1)/RQ                                                         
      DX=H(3)+H(3)+Q*XI(1)                                              
      DY=H(4)+H(4)+Q*XI(2)                                              
      DZ=H(2)+H(2)+Q*XI(3)                                              
C*****TRANSFORM BACK TO GEOMAGNETIC CO-ORDINATE SYSTEM                  
      DXM=U(1,1)*DX+U(2,1)*DY+U(3,1)*DZ                                 
      DYM=U(1,2)*DX+U(2,2)*DY                                           
      DZM=U(1,3)*DX+U(2,3)*DY+U(3,3)*DZ                                 
      DR=(XM*DXM+YM*DYM+ZM*DZM)/R                                       
C*****FORM SLOWLY VARYING EXPRESSIONS                                   
      P(4)=(WR*DXM-0.5*P(1)*DR)/(R*DZM)                                 
      P(5)=(WR*DYM-0.5*P(2)*DR)/(R*DZM)                                 
      DSQ=RQ*(DXM*DXM+DYM*DYM+DZM*DZM)
      BQ=DSQ*RQ*RQ
      P(6)=SQRT(DSQ/(RQ+3.*ZM*ZM))                                      
      P(7)=P(6)*(RQ+ZM*ZM)/(RQ*DZM)                                     
      RETURN                                                            
      END                                                               
C
C
      SUBROUTINE FELDG(GLAT,GLON,ALT,BNORTH,BEAST,BDOWN,BABS)           
C-------------------------------------------------------------------
C CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL
C REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61, 
C      1970.
C--------------------------------------------------------------------
C CHANGES (D. BILITZA, NOV 87):
C   - FIELD COEFFICIENTS IN BINARY DATA FILES INSTEAD OF BLOCK DATA
C   - CALCULATES DIPOL MOMENT
C--------------------------------------------------------------------
C  INPUT:  ENTRY POINT FELDG
C               GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
C               GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
C               ALT   ALTITUDE IN KM ABOVE SEA LEVEL
C
C          ENTRY POINT FELDC
C               V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
C                       X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
C                       Y-AXIS POINTING TO EQUATOR AT 90 LONG.
C                       Z-AXIS POINTING TO NORTH POLE
C
C          COMMON BLANK AND ENTRY POINT FELDI ARE NEEDED WHEN USED
C            IN CONNECTION WITH L-CALCULATION PROGRAM SHELLG.
C       
C          COMMON /MODEL/ AND /GENER/
C               UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
C               ERA     EARTH RADIUS FOR NORMALIZATION OF CARTESIAN 
C                       COORDINATES (6371.2 KM)
C               AQUAD, BQUAD   SQUARE OF MAJOR AND MINOR HALF AXIS FOR 
C                       EARTH ELLIPSOID AS RECOMMENDED BY INTERNATIONAL 
C                       ASTRONOMICAL UNION (6378.160, 6356.775 KM).
C               NMAX    MAXIMUM ORDER OF SPHERICAL HARMONICS
C               TIME    YEAR (DECIMAL: 1973.5) FOR WHICH MAGNETIC 
C                       FIELD IS TO BE CALCULATED
C               G(M)    NORMALIZED FIELD COEFFICIENTS (SEE FELDCOF)
C                       M=NMAX*(NMAX+2)
C------------------------------------------------------------------------
C  OUTPUT: BABS   MAGNETIC FIELD STRENGTH IN GAUSS
C          BNORTH, BEAST, BDOWN   COMPONENTS OF THE FIELD WITH RESPECT
C                 TO THE LOCAL GEODETIC COORDINATE SYSTEM, WITH AXIS
C                 POINTING IN THE TANGENTIAL PLANE TO THE NORTH, EAST
C                 AND DOWNWARD.   
C-----------------------------------------------------------------------
      DIMENSION         V(3),B(3)   
      CHARACTER*12      NAME
      COMMON            XI(3),H(144)
C      DIMENSION            XI(3),H(144)      
      COMMON/MODEL/     NAME,NMAX,TIME,G(144)  
      COMMON/GENER/     UMR,ERA,AQUAD,BQUAD
C
C-- IS RECORDS ENTRY POINT
C
C*****ENTRY POINT  FELDG  TO BE USED WITH GEODETIC CO-ORDINATES         
      IS=1                                                              
      RLAT=GLAT*UMR
      CT=SIN(RLAT)                                                      
      ST=COS(RLAT)                                                      
      D=SQRT(AQUAD-(AQUAD-BQUAD)*CT*CT)                                 
      RLON=GLON*UMR
      CP=COS(RLON)                                                      
      SP=SIN(RLON)                                                      
      ZZZ=(ALT+BQUAD/D)*CT/ERA
      RHO=(ALT+AQUAD/D)*ST/ERA
      XXX=RHO*CP                                                       
      YYY=RHO*SP                                                       
      GOTO10                                                            
      ENTRY FELDC(V,B)                                                  
C*****ENTRY POINT  FELDC  TO BE USED WITH CARTESIAN CO-ORDINATES        
      IS=2                                                              
      XXX=V(1)                                                          
      YYY=V(2)                                                          
      ZZZ=V(3)                                                          
10    RQ=1./(XXX*XXX+YYY*YYY+ZZZ*ZZZ) 
      XI(1)=XXX*RQ                                                      
      XI(2)=YYY*RQ                                                      
      XI(3)=ZZZ*RQ                                                      
      GOTO20                                                            
      ENTRY FELDI                                                      
C      ENTRY FELDI(ARG_XI,H)                                                       
C*****ENTRY POINT  FELDI  USED FOR L COMPUTATION                        
      IS=3                                                              
20    IHMAX=NMAX*NMAX+1                                                 
      LAST=IHMAX+NMAX+NMAX                                              
      IMAX=NMAX+NMAX-1                                                  
      DO 8 I=IHMAX,LAST                                                 
8     H(I)=G(I)                                                         
      DO 6 K=1,3,2                                                      
      I=IMAX                                                            
      IH=IHMAX                                                          
1     IL=IH-I                                                           
      F=2./FLOAT(I-K+2)                                                 
      X=XI(1)*F                                                         
      Y=XI(2)*F                                                         
      Z=XI(3)*(F+F)                                                     
      I=I-2                                                             
      IF(I-1)5,4,2                                                      
2     DO 3 M=3,I,2                                                      
      H(IL+M+1)=G(IL+M+1)+Z*H(IH+M+1)+X*(H(IH+M+3)-H(IH+M-1))           
     A                               -Y*(H(IH+M+2)+H(IH+M-2))           
3     H(IL+M)=G(IL+M)+Z*H(IH+M)+X*(H(IH+M+2)-H(IH+M-2))                 
     A                         +Y*(H(IH+M+3)+H(IH+M-1))                 
4     H(IL+2)=G(IL+2)+Z*H(IH+2)+X*H(IH+4)-Y*(H(IH+3)+H(IH))             
      H(IL+1)=G(IL+1)+Z*H(IH+1)+Y*H(IH+4)+X*(H(IH+3)-H(IH))             
5     H(IL)=G(IL)+Z*H(IH)+2.*(X*H(IH+1)+Y*H(IH+2))                      
      IH=IL                                                             
      IF(I.GE.K)GOTO1                                                   
6     CONTINUE                                                          
      IF(IS.EQ.3)RETURN                                                 
      S=.5*H(1)+2.*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))                   
      T=(RQ+RQ)*SQRT(RQ)                                                
      BXXX=T*(H(3)-S*XXX)                                               
      BYYY=T*(H(4)-S*YYY)                                               
      BZZZ=T*(H(2)-S*ZZZ)                                               
      IF(IS.EQ.2)GOTO7                                                  
      BABS=SQRT(BXXX*BXXX+BYYY*BYYY+BZZZ*BZZZ)
      BEAST=BYYY*CP-BXXX*SP                                             
      BRHO=BYYY*SP+BXXX*CP                                              
      BNORTH=BZZZ*ST-BRHO*CT                                            
      BDOWN=-BZZZ*CT-BRHO*ST                                            
      RETURN                                                            
7     B(1)=BXXX                                                         
      B(2)=BYYY                                                         
      B(3)=BZZZ                                                         
      RETURN                                                            
      END                                                               
C
C
      SUBROUTINE FELDCOF(YEAR,DIMO)
C------------------------------------------------------------------------
C  DETERMINES COEFFICIENTS AND DIPOL MOMENT FROM IGRF MODELS
C
C       INPUT:  YEAR    DECIMAL YEAR FOR WHICH GEOMAGNETIC FIELD IS TO
C                       BE CALCULATED
C       OUTPUT: DIMO    GEOMAGNETIC DIPOL MOMENT IN GAUSS (NORMALIZED 
C                       TO EARTH'S RADIUS) AT THE TIME (YEAR)
C  D. BILITZA, NSSDC, GSFC, CODE 633, GREENBELT, MD 20771, 
C       (301)286-9536   NOV 1987.
C  ### updated to IGRF-2000 version -dkb- 5/31/2000
C  ### updated to IGRF-2005 version -dkb- 3/24/2000
C-----------------------------------------------------------------------
      CHARACTER*12    FILMOD, FIL1, FIL2           
C ### FILMOD, DTEMOD arrays +1
      DIMENSION       GH1(144),GH2(120),GHA(144),FILMOD(14),DTEMOD(14)
      DOUBLE PRECISION X,F0,F 
      COMMON/MODEL/   FIL1,NMAX,TIME,GH1
      COMMON/GENER/   UMR,ERAD,AQUAD,BQUAD
C ### updated to 2005
      DATA            FILMOD /'dgrf45.dat', 'dgrf50.dat',            
     1                  'dgrf55.dat', 'dgrf60.dat', 'dgrf65.dat',      
     2                  'dgrf70.dat', 'dgrf75.dat', 'dgrf80.dat',      
     3                  'dgrf85.dat', 'dgrf90.dat', 'dgrf95.dat','dgrf00.dat','dgrf05.dat',
     4                  'igrf00.dat','igrf05.dat','igrf05s.dat','igrf10.dat','igrf10s.dat'/
      DATA   DTEMOD / 1945., 1950., 1955., 1960., 1965., 1970.,
     1       1975., 1980., 1985., 1990., 1995., 2000.,2005.,2010.,2015./      
C
C ### numye = numye + 1 ; is number of years represented by IGRF
C
      NUMYE=13
C
C  IS=0 FOR SCHMIDT NORMALIZATION   IS=1 GAUSS NORMALIZATION
C  IU  IS INPUT UNIT NUMBER FOR IGRF COEFFICIENT SETS
C
      IU = 10
      IS = 0
C-- DETERMINE IGRF-YEARS FOR INPUT-YEAR
      TIME = YEAR
      IYEA = INT(YEAR/5.)*5
      L = (IYEA - 1945)/5 + 1
      IF(L.LT.1) L=1
      IF(L.GT.NUMYE) L=NUMYE         
      DTE1 = DTEMOD(L)   
      FIL1 = FILMOD(L)   
      DTE2 = DTEMOD(L+1) 
      FIL2 = FILMOD(L+1) 
C-- GET IGRF COEFFICIENTS FOR THE BOUNDARY YEARS
      CALL GETSHC (IU, FIL1, NMAX1, ERAD, GH1, IER)  
      IF (IER .NE. 0) STOP                           
      CALL GETSHC (IU, FIL2, NMAX2, ERAD, GH2, IER)  
      IF (IER .NE. 0) STOP                    
C-- DETERMINE IGRF COEFFICIENTS FOR YEAR
      IF (L .LE. NUMYE-1) THEN                        
          CALL INTERSHC (YEAR, DTE1, NMAX1, GH1, DTE2, 
     1          NMAX2, GH2, NMAX, GHA)                        
      ELSE               
          CALL EXTRASHC (YEAR, DTE1, NMAX1, GH1, NMAX2,     
     1          GH2, NMAX, GHA)                                    
      ENDIF 
C-- DETERMINE MAGNETIC DIPOL MOMENT AND COEFFIECIENTS G
      F0=0.D0
      DO 1234 J=1,3
           F = GHA(J) * 1.D-5
           F0 = F0 + F * F
1234  CONTINUE
      DIMO = DSQRT(F0)

      GH1(1) =  0.0
      I=2          
      F0=1.D-5                
      IF(IS.EQ.0) F0=-F0 
      SQRT2=SQRT(2.)      

      DO 9 N=1,NMAX           
        X = N
        F0 = F0 * X * X / (4.D0 * X - 2.D0)               
        IF(IS.EQ.0) F0 = F0 * (2.D0 * X - 1.D0) / X
        F = F0 * 0.5D0                                    
        IF(IS.EQ.0) F = F * SQRT2
        GH1(I) = GHA(I-1) * F0
        I = I+1                                         
      DO 9 M=1,N                                    
        F = F * (X + M) / (X - M + 1.D0)                 
        IF(IS.EQ.0) F = F * DSQRT((X - M + 1.D0) / (X + M))             
        GH1(I) = GHA(I-1) * F
        GH1(I+1) = GHA(I) * F
        I=I+2
9     CONTINUE                                          
      RETURN
      END
C
C
      SUBROUTINE GETSHC (IU, FSPEC, NMAX, ERAD, GH, IER)           
                                                                                
C ===============================================================               
C                                                                               
C       Version 1.01                                                 
C                                                                               
C       Reads spherical harmonic coefficients from the specified     
C       file into an array.                                          
C                                                                               
C       Input:                                                       
C           IU    - Logical unit number                              
C           FSPEC - File specification                               
C                                                                               
C       Output:                                                      
C           NMAX  - Maximum degree and order of model                
C           ERAD  - Earth's radius associated with the spherical     
C                   harmonic coefficients, in the same units as      
C                   elevation                                        
C           GH    - Schmidt quasi-normal internal spherical          
C                   harmonic coefficients                            
C           IER   - Error number: =  0, no error                     
C                                 = -2, records out of order         
C                                 = FORTRAN run-time error number    
C                                                                    
C       A. Zunde                                                     
C       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225    
C                                                                               
C ===============================================================               
                                                                                
      CHARACTER  FSPEC*(*), FOUT*55                                    
      DIMENSION       GH(*)                                        
C ---------------------------------------------------------------               
C       Open coefficient file. Read past first header record.        
C       Read degree and order of model and Earth's radius.           
C ---------------------------------------------------------------               
      WRITE(FOUT,667) FSPEC
c 667  FORMAT('/usr/local/etc/httpd/cgi-bin/natasha/IRI/',A12)
667   FORMAT(A12)
      OPEN (IU, FILE=FOUT, STATUS='OLD', IOSTAT=IER, ERR=999)     
     
      READ (IU, *, IOSTAT=IER, ERR=999)                            
      READ (IU, *, IOSTAT=IER, ERR=999) NMAX, ERAD                 
C ---------------------------------------------------------------               
C       Read the coefficient file, arranged as follows:              
C                                                                               
C                                       N     M     G     H          
C                                       ----------------------       
C                                   /   1     0    GH(1)  -          
C                                  /    1     1    GH(2) GH(3)       
C                                 /     2     0    GH(4)  -          
C                                /      2     1    GH(5) GH(6)       
C           NMAX*(NMAX+3)/2     /       2     2    GH(7) GH(8)       
C              records          \       3     0    GH(9)  -          
C                                \      .     .     .     .          
C                                 \     .     .     .     .          
C           NMAX*(NMAX+2)          \    .     .     .     .          
C           elements in GH          \  NMAX  NMAX   .     .          
C                                                                               
C       N and M are, respectively, the degree and order of the       
C       coefficient.                                                 
C ---------------------------------------------------------------               
                                                                                
      I = 0                                                        
      DO 2211 NN = 1, NMAX                                              
          DO 2233 MM = 0, NN                                            
              READ (IU, *, IOSTAT=IER, ERR=999) N, M, G, H         
              IF (NN .NE. N .OR. MM .NE. M) THEN                   
                  IER = -2                                         
                  GOTO 999                                         
              ENDIF                                                
              I = I + 1                                            
              GH(I) = G                                            
              IF (M .NE. 0) THEN                                   
                  I = I + 1                                        
                  GH(I) = H                                        
              ENDIF                                                
2233      CONTINUE                                                    
2211  CONTINUE                                                        
                                                                                
999   CLOSE (IU)                                                   
                                                                                
      RETURN                                                       
      END                                                          
C
C
      SUBROUTINE INTERSHC (DATE, DTE1, NMAX1, GH1, DTE2,          
     1                        NMAX2, GH2, NMAX, GH)                  
                                                                                
C ===============================================================               
C                                                                               
C       Version 1.01                                                 
C                                                                               
C       Interpolates linearly, in time, between two spherical        
C       harmonic models.                                             
C                                                                               
C       Input:                                                       
C           DATE  - Date of resulting model (in decimal year)        
C           DTE1  - Date of earlier model                            
C           NMAX1 - Maximum degree and order of earlier model        
C           GH1   - Schmidt quasi-normal internal spherical          
C                   harmonic coefficients of earlier model           
C           DTE2  - Date of later model                              
C           NMAX2 - Maximum degree and order of later model          
C           GH2   - Schmidt quasi-normal internal spherical          
C                   harmonic coefficients of later model             
C                                                                               
C       Output:                                                      
C           GH    - Coefficients of resulting model                  
C           NMAX  - Maximum degree and order of resulting model      
C                                                                               
C       A. Zunde                                                     
C       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225    
C                                                                               
C ===============================================================               
                                                                                
      DIMENSION       GH1(*), GH2(*), GH(*)                        
                                                                                
C ---------------------------------------------------------------               
C       The coefficients (GH) of the resulting model, at date        
C       DATE, are computed by linearly interpolating between the     
C       coefficients of the earlier model (GH1), at date DTE1,       
C       and those of the later model (GH2), at date DTE2. If one     
C       model is smaller than the other, the interpolation is        
C       performed with the missing coefficients assumed to be 0.     
C ---------------------------------------------------------------               
                                                                                
      FACTOR = (DATE - DTE1) / (DTE2 - DTE1)                       
                                                                               
      IF (NMAX1 .EQ. NMAX2) THEN                                   
            K = NMAX1 * (NMAX1 + 2)                                  
            NMAX = NMAX1                                             
      ELSE IF (NMAX1 .GT. NMAX2) THEN                              
            K = NMAX2 * (NMAX2 + 2)                                  
            L = NMAX1 * (NMAX1 + 2)                                  
            DO 1122 I = K + 1, L                                          
1122            GH(I) = GH1(I) + FACTOR * (-GH1(I))                  
            NMAX = NMAX1                                             
      ELSE                                                         
            K = NMAX1 * (NMAX1 + 2)                                  
            L = NMAX2 * (NMAX2 + 2)                                  
            DO 1133 I = K + 1, L                                          
1133            GH(I) = FACTOR * GH2(I)                              
            NMAX = NMAX2                                             
      ENDIF                                                        
                                                                                
      DO 1144 I = 1, K                                                  
1144        GH(I) = GH1(I) + FACTOR * (GH2(I) - GH1(I))              
                                                                                
      RETURN                                                       
      END                                                          
C
C
      SUBROUTINE EXTRASHC (DATE, DTE1, NMAX1, GH1, NMAX2,           
     1                        GH2, NMAX, GH)                           
                                                                                
C ===============================================================               
C                                                                               
C       Version 1.01                                                   
C                                                                               
C       Extrapolates linearly a spherical harmonic model with a        
C       rate-of-change model.                                          
C                                                                               
C       Input:                                                         
C           DATE  - Date of resulting model (in decimal year)          
C           DTE1  - Date of base model                                 
C           NMAX1 - Maximum degree and order of base model             
C           GH1   - Schmidt quasi-normal internal spherical            
C                   harmonic coefficients of base model                
C           NMAX2 - Maximum degree and order of rate-of-change         
C                   model                                              
C           GH2   - Schmidt quasi-normal internal spherical            
C                   harmonic coefficients of rate-of-change model      
C                                                                               
C       Output:                                                        
C           GH    - Coefficients of resulting model                    
C           NMAX  - Maximum degree and order of resulting model        
C                                                                               
C       A. Zunde                                                       
C       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225      
C                                                                               
C ===============================================================               
                                                                                
      DIMENSION       GH1(*), GH2(*), GH(*)                          
                                                                                
C ---------------------------------------------------------------               
C       The coefficients (GH) of the resulting model, at date          
C       DATE, are computed by linearly extrapolating the coef-         
C       ficients of the base model (GH1), at date DTE1, using          
C       those of the rate-of-change model (GH2), at date DTE2. If      
C       one model is smaller than the other, the extrapolation is      
C       performed with the missing coefficients assumed to be 0.       
C ---------------------------------------------------------------               
                                                                                
      FACTOR = (DATE - DTE1)                                         
                                                                                
      IF (NMAX1 .EQ. NMAX2) THEN                                     
            K = NMAX1 * (NMAX1 + 2)                                    
            NMAX = NMAX1                                               
      ELSE IF (NMAX1 .GT. NMAX2) THEN                                
            K = NMAX2 * (NMAX2 + 2)                                    
            L = NMAX1 * (NMAX1 + 2)                                    
            DO 1155 I = K + 1, L                                            
1155            GH(I) = GH1(I)                                         
            NMAX = NMAX1                                               
      ELSE                                                           
            K = NMAX1 * (NMAX1 + 2)                                    
            L = NMAX2 * (NMAX2 + 2)                                    
            DO 1166 I = K + 1, L                                            
1166            GH(I) = FACTOR * GH2(I)                                
            NMAX = NMAX2                                               
      ENDIF                                                          
                                                                                
      DO 1177 I = 1, K                                                    
1177        GH(I) = GH1(I) + FACTOR * GH2(I)                           
                                                                                
      RETURN                                                         
      END                                                            
C
C
      SUBROUTINE INITIZE
C----------------------------------------------------------------
C Initializes the parameters in COMMON/GENER/
C
C       UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
C       ERA     EARTH RADIUS FOR NORMALIZATION OF CARTESIAN 
C                       COORDINATES (6371.2 KM) 
C       EREQU   MAJOR HALF AXIS FOR EARTH ELLIPSOID (6378.160 KM)
C       ERPOL   MINOR HALF AXIS FOR EARTH ELLIPSOID (6356.775 KM)
C       AQUAD   SQUARE OF MAJOR HALF AXIS FOR EARTH ELLIPSOID
C       BQUAD   SQUARE OF MINOR HALF AXIS FOR EARTH ELLIPSOID
C
C ERA, EREQU and ERPOL as recommended by the INTERNATIONAL 
C ASTRONOMICAL UNION .
C-----------------------------------------------------------------
      COMMON/GENER/   UMR,ERA,AQUAD,BQUAD
      ERA=6371.2
      EREQU=6378.16
      ERPOL=6356.775
      AQUAD=EREQU*EREQU
      BQUAD=ERPOL*ERPOL
      UMR=ATAN(1.0)*4./180.
      RETURN
      END

      SUBROUTINE RECALC(IYR,IDAY,IHOUR,MIN,ISEC)
      
C  PREPARES ELEMENTS OF ROTATION MATRICES FOR TRANSFORMATIONS OF VECTORS BETWEEN
C  SEVERAL COORDINATE SYSTEMS, MOST FREQUENTLY USED IN SPACE PHYSICS.
C
C  THIS SUBROUTINE SHOULD BE INVOKED BEFORE USING THE FOLLOWING SUBROUTINES:
C         GEOGSM, MAGSM, SMGSM, GSMGSE, GEIGEO.
C
C  THERE IS NO NEED TO REPEATEDLY INVOKE RECALC, IF MULTIPLE CALCULATIONS ARE MADE
C    FOR THE SAME DATE AND TIME.
C
C-----INPUT PARAMETERS:

C     IYR   -  YEAR NUMBER (FOUR DIGITS)
C     IDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
C     IHOUR -  HOUR OF DAY (00 TO 23)
C     MIN   -  MINUTE OF HOUR (00 TO 59)
C     ISEC  -  SECONDS OF MINUTE (00 TO 59)
C
C-----OUTPUT PARAMETERS:   NONE (ALL OUTPUT QUANTITIES ARE PLACED
C                                 INTO THE COMMON BLOCK /GEOPACK/)
C
C    OTHER SUBROUTINES CALLED BY THIS ONE: SUN
C
C   ################################################
C   #  WRITTEN BY  N.A. TSYGANENKO ON DEC.1, 1991  #
C   ################################################
c
c    Last modification:  January 5, 2001.
c    The code has been modified to accept dates through 2005.
c    Also, a "save" statement was added, to avoid potential problems
c    with some Fortran compilers.
C
      COMMON /GEOPACK/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS
      COMMON /GEOPACK/ CPS,SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3
      COMMON /GEOPACK/  K,IY,CGST,SGST,BA(6)
c
      SAVE IYE,IDE,IPR
      DATA IYE,IDE,IPR/3*0/
C
      IF (IYR.EQ.IYE.AND.IDAY.EQ.IDE) GOTO 5
C
C   IYE AND IDE ARE THE CURRENT VALUES OF YEAR AND DAY NUMBER
C
      IY=IYR
      IDE=IDAY
      IF(IY.LT.1965) IY=1965
      IF(IY.GT.2005) IY=2005
C
C  WE ARE RESTRICTED BY THE INTERVAL 1965-2005,
C  FOR WHICH THE IGRF COEFFICIENTS ARE KNOWN; IF IYR IS OUTSIDE THIS INTERVAL
C  THE SUBROUTINE PRINTS A WARNING (BUT DOES NOT REPEAT IT AT NEXT INVOCATIONS)
C
      IF(IY.NE.IYR.AND.IPR.EQ.0) WRITE (*,10) IYR,IY
      IF(IY.NE.IYR) IPR=1
      IYE=IY
C
C  LINEAR INTERPOLATION OF THE GEODIPOLE MOMENT COMPONENTS BETWEEN THE
C  VALUES FOR THE NEAREST EPOCHS:
C
      IF (IY.LT.1970) THEN                            !1965-1970
        F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1965.)/5.
        F1=1.D0-F2
           G10=30334.*F1+30220.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
        G11=-2119.*F1-2068.*F2
        H11=5776.*F1+5737.*F2
      ELSEIF (IY.LT.1975) THEN                        !1970-1975
        F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1970.)/5.
        F1=1.D0-F2
           G10=30220.*F1+30100.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
        G11=-2068.*F1-2013.*F2
        H11=5737.*F1+5675.*F2
      ELSEIF (IY.LT.1980) THEN                        !1975-1980
        F2=(DFLOAT(IY)+DFLOAT(IDAY)/365.-1975.)/5.
        F1=1.D0-F2
           G10=30100.*F1+29992.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
        G11=-2013.*F1-1956.*F2
        H11=5675.*F1+5604.*F2
      ELSEIF (IY.LT.1985) THEN                        !1980-1985
        F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1980.)/5.
        F1=1.D0-F2
           G10=29992.*F1+29873.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
        G11=-1956.*F1-1905.*F2
        H11=5604.*F1+5500.*F2
      ELSEIF (IY.LT.1990) THEN!1985-1990
        F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1985.)/5.
        F1=1.D0-F2
           G10=29873.*F1+29775.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
        G11=-1905.*F1-1848.*F2
        H11=5500.*F1+5406.*F2
      ELSEIF (IY.LT.1995) THEN!1990-1995
        F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1990.)/5.
        F1=1.D0-F2
           G10=29775.*F1+29682.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
        G11=-1848.*F1-1789.*F2
        H11=5406.*F1+5318.*F2
      ELSEIF (IY.LT.2000) THEN                        !1995-2000
        F2=(FLOAT(IY)+FLOAT(IDAY)/365.-1995.)/5.
        F1=1.D0-F2
        G10=29682.*F1+29615.*F2 ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
        G11=-1789.*F1-1728.*F2
        H11=5318.*F1+5186.*F2
      ELSE                                            !2000-2005
C
C   LINEAR EXTRAPOLATION BEYOND 2000 BY USING SECULAR VELOCITY COEFFICIENTS:
C
        DT=FLOAT(IY)+FLOAT(IDAY)/366.-2000.
        G10=29615.-14.6*DT      ! HERE G10 HAS OPPOSITE SIGN TO THAT IN IGRF TABLES
        G11=-1728.+10.7*DT
        H11=5186.-22.5*DT
      ENDIF
C
C  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EzMAG IN GEO COORD.SYSTEM:
C   SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
C         ST0 * CL0                ST0 * SL0                CT0
C
      SQ=G11**2+H11**2
      SQQ=SQRT(SQ)
      SQR=SQRT(G10**2+SQ)
      SL0=-H11/SQQ
      CL0=-G11/SQQ
      ST0=SQQ/SQR
      CT0=G10/SQR
      STCL=ST0*CL0
      STSL=ST0*SL0
      CTSL=CT0*SL0
      CTCL=CT0*CL0
C
C      THE CALCULATIONS ARE TERMINATED IF ONLY GEO-MAG TRANSFORMATION
C       IS TO BE DONE  (IHOUR>24 IS THE AGREED INDICATOR FOR THIS CASE):
C
5     IF (IHOUR.GT.24) RETURN
C
      CALL SUN(IY,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C
C  S1,S2, AND S3 ARE THE COMPONENTS OF THE UNIT VECTOR EXGSM=EXGSE IN THE
C   SYSTEM GEI POINTING FROM THE EARTH'S CENTER TO THE SUN:
C
      S1=COS(SRASN)*COS(SDEC)
      S2=SIN(SRASN)*COS(SDEC)
      S3=SIN(SDEC)
      CGST=COS(GST)
      SGST=SIN(GST)
C
C  DIP1, DIP2, AND DIP3 ARE THE COMPONENTS OF THE UNIT VECTOR EZSM=EZMAG
C   IN THE SYSTEM GEI:
C
      DIP1=STCL*CGST-STSL*SGST
      DIP2=STCL*SGST+STSL*CGST
      DIP3=CT0
C
C  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EYGSM IN THE SYSTEM GEI
C   BY TAKING THE VECTOR PRODUCT D x S AND NORMALIZING IT TO UNIT LENGTH:
C
      Y1=DIP2*S3-DIP3*S2
      Y2=DIP3*S1-DIP1*S3
      Y3=DIP1*S2-DIP2*S1
      Y=SQRT(Y1*Y1+Y2*Y2+Y3*Y3)
      Y1=Y1/Y
      Y2=Y2/Y
      Y3=Y3/Y
C
C   THEN IN THE GEI SYSTEM THE UNIT VECTOR Z = EZGSM = EXGSM x EYGSM = S x Y
C    HAS THE COMPONENTS:
C
      Z1=S2*Y3-S3*Y2
      Z2=S3*Y1-S1*Y3
      Z3=S1*Y2-S2*Y1
C
C    THE VECTOR EZGSE (HERE DZ) IN GEI HAS THE COMPONENTS (0,-SIN(DELTA),
C     COS(DELTA)) = (0.,-0.397823,0.917462); HERE DELTA = 23.44214 DEG FOR
C   THE EPOCH 1978 (SEE THE BOOK BY GUREVICH OR OTHER ASTRONOMICAL HANDBOOKS).
C    HERE THE MOST ACCURATE TIME-DEPENDENT FORMULA IS USED:
C
      DJ=FLOAT(365*(IY-1900)+(IY-1901)/4 +IDAY)-0.5+
     &   FLOAT(IHOUR*3600+MIN*60+ISEC)/86400.
      T=DJ/36525.
      OBLIQ=(23.45229-0.0130125*T)/57.2957795
      DZ1=0.
      DZ2=-SIN(OBLIQ)
      DZ3=COS(OBLIQ)
C
C  THEN THE UNIT VECTOR EYGSE IN GEI SYSTEM IS THE VECTOR PRODUCT DZ x S :
C
      DY1=DZ2*S3-DZ3*S2
      DY2=DZ3*S1-DZ1*S3
      DY3=DZ1*S2-DZ2*S1
C
C   THE ELEMENTS OF THE MATRIX GSE TO GSM ARE THE SCALAR PRODUCTS:
C   CHI=EM22=(EYGSM,EYGSE), SHI=EM23=(EYGSM,EZGSE), EM32=(EZGSM,EYGSE)=-EM23,
C     AND EM33=(EZGSM,EZGSE)=EM22
C
      CHI=Y1*DY1+Y2*DY2+Y3*DY3
      SHI=Y1*DZ1+Y2*DZ2+Y3*DZ3
      HI=ASIN(SHI)
C
C    TILT ANGLE: PSI=ARCSIN(DIP,EXGSM)
C
      SPS=DIP1*S1+DIP2*S2+DIP3*S3
      CPS=SQRT(1.-SPS**2)
      PSI=ASIN(SPS)
C
C    THE ELEMENTS OF THE MATRIX MAG TO SM ARE THE SCALAR PRODUCTS:
C CFI=GM22=(EYSM,EYMAG), SFI=GM23=(EYSM,EXMAG); THEY CAN BE DERIVED AS FOLLOWS:
C
C IN GEO THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS (CT0*CL0,CT0*SL0,-ST0)
C  AND (-SL0,CL0,0), RESPECTIVELY.    HENCE, IN GEI THE COMPONENTS ARE:
C  EXMAG:    CT0*CL0*COS(GST)-CT0*SL0*SIN(GST)
C            CT0*CL0*SIN(GST)+CT0*SL0*COS(GST)
C            -ST0
C  EYMAG:    -SL0*COS(GST)-CL0*SIN(GST)
C            -SL0*SIN(GST)+CL0*COS(GST)
C             0
C  THE COMPONENTS OF EYSM IN GEI WERE FOUND ABOVE AS Y1, Y2, AND Y3;
C  NOW WE ONLY HAVE TO COMBINE THE QUANTITIES INTO SCALAR PRODUCTS:
C
      EXMAGX=CT0*(CL0*CGST-SL0*SGST)
      EXMAGY=CT0*(CL0*SGST+SL0*CGST)
      EXMAGZ=-ST0
      EYMAGX=-(SL0*CGST+CL0*SGST)
      EYMAGY=-(SL0*SGST-CL0*CGST)
      CFI=Y1*EYMAGX+Y2*EYMAGY
      SFI=Y1*EXMAGX+Y2*EXMAGY+Y3*EXMAGZ
C
      XMUT=(ATAN2(SFI,CFI)+3.1415926536)*3.8197186342
C
C  THE ELEMENTS OF THE MATRIX GEO TO GSM ARE THE SCALAR PRODUCTS:
C
C   A11=(EXGEO,EXGSM), A12=(EYGEO,EXGSM), A13=(EZGEO,EXGSM),
C   A21=(EXGEO,EYGSM), A22=(EYGEO,EYGSM), A23=(EZGEO,EYGSM),
C   A31=(EXGEO,EZGSM), A32=(EYGEO,EZGSM), A33=(EZGEO,EZGSM),
C
C   ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:
C
C  EXGEO=(CGST,SGST,0), EYGEO=(-SGST,CGST,0), EZGEO=(0,0,1)
C  EXGSM=(S1,S2,S3),  EYGSM=(Y1,Y2,Y3),   EZGSM=(Z1,Z2,Z3)
C                                                           AND  THEREFORE:
C
      A11=S1*CGST+S2*SGST
      A12=-S1*SGST+S2*CGST
      A13=S3
      A21=Y1*CGST+Y2*SGST
      A22=-Y1*SGST+Y2*CGST
      A23=Y3
      A31=Z1*CGST+Z2*SGST
      A32=-Z1*SGST+Z2*CGST
      A33=Z3
C
10    FORMAT('*YEAR IS OUT OF INTERVAL 1965-2005: IYR=',
     & I4,'CALCULATIONS WILL BE DONE FOR IYR=',I4,/)
      RETURN
      END

C----------------------------------------------------------------------------
      SUBROUTINE GEOMAG(XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,J,IYR)
C
C    CONVERTS GEOGRAPHIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICA VERSA.
C    IYR IS YEAR NUMBER (FOUR DIGITS).
C
C                           J >  0                J < 0
C-----INPUT:  J,XGEO,YGEO,ZGEO,IYR   J,XMAG,YMAG,ZMAG,IYR
C-----OUTPUT:    XMAG,YMAG,ZMAG        XGEO,YGEO,ZGEO
C
C   OTHER SUBROUTINES CALLED BY THIS ONE:  RECALC
C
C   LAST MOFIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
C                         SOME REDUNDANT STATEMENTS TAKEN OUT)
C
C   WRITTEN BY:  N. A. TSYGANENKO
C
      REAL XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,IYR
      INTEGER J
      COMMON /GEOPACK/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,AB(19),K,IY,BB(8)
      SAVE II
      DATA II/1/
C
      IF(IYR.EQ.II) GOTO 1
      II=IYR
      CALL RECALC(II,0,25,0,0)
1     CONTINUE
      IF(J.LT.0) GOTO 2
      XMAG=XGEO*CTCL+YGEO*CTSL-ZGEO*ST0
      YMAG=YGEO*CL0-XGEO*SL0
      ZMAG=XGEO*STCL+YGEO*STSL+ZGEO*CT0
      RETURN
2     XGEO=XMAG*CTCL-YMAG*SL0+ZMAG*STCL
      YGEO=XMAG*CTSL+YMAG*CL0+ZMAG*STSL
      ZGEO=ZMAG*CT0-XMAG*ST0
      RETURN
      END
c

C---------------------------------------------------------------------------

      SUBROUTINE SUN(IYR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
C
C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
C
C-------  INPUT PARAMETERS:
C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
C
C-------  OUTPUT PARAMETERS:
C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
C
C     LAST MODIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
C     SOME REDUNDANT STATEMENTS TAKEN OUT FROM THE PREVIOUS VERSION)
C
C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
C
      DOUBLE PRECISION DJ,FDAY
      DATA RAD/57.295779513/
C
      IF(IYR.LT.1901.OR.IYR.GT.2099) RETURN
      FDAY=DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
      DJ=365*(IYR-1900)+(IYR-1901)/4+IDAY-0.5D0+FDAY
      T=DJ/36525.
      VL=DMOD(279.696678+0.9856473354*DJ,360.D0)
      GST=DMOD(279.690983+.9856473354*DJ+360.*FDAY+180.,360.D0)/RAD
      G=DMOD(358.475845+0.985600267*DJ,360.D0)/RAD
      SLONG=(VL+(1.91946-0.004789*T)*SIN(G)+0.020094*SIN(2.*G))/RAD
      IF(SLONG.GT.6.2831853) SLONG=SLONG-6.2831853
      IF (SLONG.LT.0.) SLONG=SLONG+6.2831853
      OBLIQ=(23.45229-0.0130125*T)/RAD
      SOB=SIN(OBLIQ)
      SLP=SLONG-9.924E-5
C
C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
C   THE ORBITAL MOTION OF THE EARTH
C
      SIND=SOB*SIN(SLP)
      COSD=SQRT(1.-SIND**2)
      SC=SIND/COSD
      SDEC=ATAN(SC)
      SRASN=3.141592654-ATAN2(COS(OBLIQ)/SOB*SC,-COS(SLP)/COSD)
      RETURN
      END
      
C------------------------------------------------------------------------
      SUBROUTINE SPHCAR(R,TETA,PHI,X,Y,Z,J)
C
C   CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICA VERSA
C    (TETA AND PHI IN RADIANS).
C
C                  J > 0            J < 0
C-----INPUT:   J,R,TETA,PHI     J,X,Y,Z
C----OUTPUT:      X,Y,Z        R,TETA,PHI
C
C   LAST MOFIFICATION:  JAN 5, 2001 (NO ESSENTIAL CHANGES, BUT
C                        SOME REDUNDANT STATEMENTS TAKEN OUT)
C
C   WRITTEN BY:  N. A. TSYGANENKO
C
      REAL R,TETA,PHI,X,Y,Z
      INTEGER J
      
      IF(J.GT.0) GOTO 3
      SQ=X**2+Y**2
      R=SQRT(SQ+Z**2)
      IF (SQ.NE.0.) GOTO 2
      PHI=0.
      IF (Z.LT.0.) GOTO 1
      TETA=0.
      RETURN
1     TETA=3.141592654
      RETURN
2     SQ=SQRT(SQ)
      PHI=ATAN2(Y,X)
      TETA=ATAN2(SQ,Z)
      IF (PHI.LT.0.) PHI=PHI+6.28318531
      RETURN
3     SQ=R*SIN(TETA)
      X=SQ*COS(PHI)
      Y=SQ*SIN(PHI)
      Z=R*COS(TETA)
      RETURN
      END

C Based on BILCAL, VERSION 3.0, AUGUST 1995
C
C Alexey Nikolayevich Petrov
C SINP MSU Russia
C gluk@srd.sinp.msu.ru
C
C*****************************************************************
      INTEGER       INFILE,OUTFILE
      REAL          YEAR,HEIGHT
      REAL          LATI,LONGI,MR,MLAT,MLONG
      REAL          XMAG,YMAG,ZMAG,COLAT,XLON
      CHARACTER*11  OUTFILENAME
      CHARACTER*11  INFILENAME      
      LOGICAL       VAL
      DIMENSION     VARE(4),VARB(4)
      COMMON/GENER/ UMR,ERA,AQUAD,BQUAD
C
      DATA LATI,LONGI,HEIGHT,YEAR
     &                  /45.1,293.1,100,1985.5/
c### year limit modified
      DATA VARB         /-90.0,-360.0,0.00000,1940.0/
      DATA VARE         /+90.0,+360.0,30000.0,2010.0/
C
      CALL INITIZE

      MONITOR=6
      OUTFILE=16
      INFILE=17

c      WRITE(MONITOR, *) 'TYPE INPUT FILE NAME: '
c      READ '(13A)', INFILENAME
c      WRITE(MONITOR, *) 'INFILENAME=', INFILENAME
      
c      WRITE(MONITOR, *) 'TYPE RESULT FILE NAME: '
c      READ '(13A)', OUTFILENAME      
c      WRITE(MONITOR, *) 'OUTFILENAME=', OUTFILENAME
    
      INFILENAME='Path.txt'
      OUTFILENAME='Result.txt'      
      OPEN(UNIT=OUTFILE,FILE=OUTFILENAME,
     &     STATUS='NEW',FORM='FORMATTED') 
      OPEN(UNIT=INFILE,FILE=INFILENAME,STATUS='OLD',
     &     FORM='FORMATTED',ERR=8000)     

      WRITE(OUTFILE,3910) 
3910  FORMAT(
     & '  YEAR      ALT    LON     LAT     DIMO   B/B0     B       B0',
     & '      B-NORTH  B-EAST',
     & '  B-DOWN    DIP    DEC  L-VALUE MLONG  MLAT  C')  
      
1000  READ(INFILE, *, END=9999) YEAR, HEIGHT, LATI, LONGI
     
      CALL FELDCOF(YEAR,DIMO)
      CALL FELDG(LATI,LONGI,HEIGHT,BNORTH,BEAST,BDOWN,BABS)
      CALL SHELLG(LATI,LONGI,HEIGHT,DIMO,XL,ICODE,BAB1)
C      CALL GGM(0,LONGI,LATI,MLONG,MLAT)
      COLAT=(90.0 - LATI) * 0.01745329
      XLON=LONGI * 0.01745329
      CALL SPHCAR(HEIGHT,COLAT,XLON,XGEO,YGEO,ZGEO,1)
      CALL GEOMAG(XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,1,YEAR)
      CALL SPHCAR(MR,MLAT,MLONG,XMAG,YMAG,ZMAG,-1)
      MLAT = 90 - MLAT/0.01745329
      MLONG = MLONG/0.01745329

      BEQU=DIMO/(XL*XL*XL)
      BDEL=1.E-3
      CALL FINDB0(0.05,BDEL,VAL,BEQ,RR0)
      IF(VAL) BEQU=BEQ
      DIP=ASIN(BDOWN/BABS)/UMR
      DEC=ASIN(BEAST/SQRT(BEAST*BEAST+BNORTH*BNORTH))/UMR
      BBX=BABS/BEQU
      IF(BBX.GT.9999.999) BBX=9999.999
      WRITE(OUTFILE,7177) YEAR,HEIGHT,LONGI,LATI,DIMO,BBX,BABS,BEQU,
     &      BNORTH,BEAST,BDOWN,DIP,DEC,XL,MLONG,MLAT,ICODE           
7177  FORMAT(F8.1,F8.1,F8.2,F8.2,F8.4,F9.3,
     &      (1X,F7.5),(1X,F7.5),3(1X,F7.5),
     &      2F7.1,F8.3,F8.2,F8.2,I3)
      GOTO 1000
8000  WRITE(MONITOR, *) 'INPUT FILE NOT FOUND'    
9999  STOP
      END
