C***********************************************************************
FUNCTION ALD(EX,SPIN,IP)
    C***********************************************************************
    C
    C     This subroutine provides cubic interpolation of values
    C     of level density - based on the function AICC.
    C     At the lowest excitations (very low density) the linear interpolation 
    c     is used while at higher energies a cubic interpolation is adopted 
    C                                                    Version from 16-MAY-09
    C
          INTEGER*4    MAXIS
          PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
         *              MAXBIN = 700,  !maximum # of bins in continuum
         *              MAXIS  = 76000000        )   
          DIMENSION XX(4),YY(4),A(4)
          COMMON /LDTAB/ TABENLD(0:270),TABLD(0:270,0:MAXJC,0:1),NLD
    C
          ALD = 0.0
    c
    c    Low level densities treated in a special way
    c
          NLMIN = 0
          DO I = 1, NLD
            IF (TABLD(I,ISUBSC(SPIN),IP).GT.0.0) THEN
              NLMIN = I - 1
              GOTO 9
            ENDIF
          ENDDO
          GOTO 11
       9  CONTINUE
          IF (EX.LE.TABENLD(NLMIN)) GOTO 11   ! rho = 0
          IF (EX.LE.TABENLD(NLMIN+1)) THEN    ! linear interpolation for low rho
            EXPOM = (EX - TABENLD(NLMIN)) / 
         *          (TABENLD(NLMIN + 1) - TABENLD(NLMIN))
            ALD =TABLD(NLMIN+1,ISUBSC(SPIN),IP)-TABLD(NLMIN,ISUBSC(SPIN),IP)
            ALD = ALD * EXPOM
            GOTO 11
          ENDIF
          IF (EX.LE.TABENLD(NLMIN+2)) THEN    
           EXPOM = (EX - TABENLD(NLMIN + 1)) / 
         *         (TABENLD(NLMIN + 2) - TABENLD(NLMIN + 1))
           ALD=TABLD(NLMIN+2,ISUBSC(SPIN),IP)-TABLD(NLMIN+1,ISUBSC(SPIN),IP)
           ALD = TABLD(NLMIN+1,ISUBSC(SPIN),IP) + ALD * EXPOM
           GOTO 11
          ENDIF
    c
          IF (EX.LE.TABENLD(2)) THEN
             K=0
             ELSE
             IF (EX.LE.TABENLD(NLD-1)) THEN
                DO I=3,NLD-1
                   IF (EX.LE.TABENLD(I)) THEN
                      K=I-3
                      GO TO 1
                   ENDIF
                ENDDO !I
                ELSE
                K=NLD-4
             ENDIF
          ENDIF
    C
        1 exlog=log(ex)
          DO J=1,4
            if ((TABENLD(K+J)).LT.0.0) WRITE(*,*)'A problem with NLD energy'
            if (TABLD(K+J,ISUBSC(SPIN),IP).LE.0.0) then
             TABLD(K+J,ISUBSC(SPIN),IP) = 1.0e-2
            endif 
            XX(J)=log(TABENLD(K+J))
            YY(J)=log(TABLD(K+J,ISUBSC(SPIN),IP))
          ENDDO !J
          
          DO I=1,4
             A(I)=YY(I)
             DO J=1,4
                IF (I.NE.J) A(I)=A(I)/(XX(J)-XX(I))
             ENDDO !J
             DO J=1,4
                IF (I.NE.J) A(I)=A(I)*(XX(J)-EXLOG)
             ENDDO !J
             ALD=ALD+A(I)
          ENDDO !I
    C
          ALD=exp(ALD)
    
       11 CONTINUE
          RETURN
          END
    C

    C***********************************************************************
      FUNCTION AICC(ETRA,TABEN,TABICC,MAEL,MUL,N)
C***********************************************************************
C
C     This subroutine provides cubic interpolation of values
C     Y=TABICC(MAEL,MUL,I) for fixed MAEL and MUL and variable I.
C     These values Y for various I are assumed to represent internal
C     conversion coefficients (ICC's) for generally non-equdistant
C     transition energies that are specified by TABEN(I). MAEL stands for
C     the type of a transition (0 for M1 and 1 for E1) while MUL means
C     multipolarity. N is the length of the table for a fixed MAEL and MUL.
C     Be careful, if ETRA (the energy at which it is desirable to get
C     the ICC-coefficient) is lower than TABEN(1) or higher than TABEN(N),
C     then the interpolation changes to extrapolation and difficulties may
C     start ...
C                                                    Version from 6-OCT-95
C
      DIMENSION TABEN(50),TABICC(0:1,5,50),XX(4),YY(4),A(4)
C
      AICC=0.0
      IF (ETRA.LT.TABEN(1)) THEN
        AICC = TABICC(MAEL,MUL,1)
        GOTO 2
      ENDIF             ! ETRA=TABEN(1)
      IF (ETRA.GT.TABEN(N)) ETRA=TABEN(N)
      IF (ETRA.LE.TABEN(2)) THEN
         K=0
         ELSE
         IF (ETRA.LE.TABEN(N-1)) THEN
            DO I=3,N-1
               IF (ETRA.LE.TABEN(I)) THEN
                  K=I-3
                  GO TO 1
               ENDIF
            ENDDO !I
            ELSE
            K=N-4
         ENDIF
      ENDIF
C
    1 etrapom=etra
      etra=log(etra)
      DO J=1,4
         XX(J)=log(TABEN(K+J))
         YY(J)=log(TABICC(MAEL,MUL,K+J))
      ENDDO !J
      AICC=0.
      DO I=1,4
         A(I)=YY(I)
         DO J=1,4
            IF (I.NE.J) A(I)=A(I)/(XX(J)-XX(I))
         ENDDO !J
         DO J=1,4
            IF (I.NE.J) A(I)=A(I)*(XX(J)-ETRA)
         ENDDO !J
         AICC=AICC+A(I)
      ENDDO !I
C
      aicc=exp(aicc)
      etra=etrapom
    2 CONTINUE
      RETURN
      END
      FUNCTION APSF(EGX,MTYP)
        C***********************************************************************
        C
        C     This subroutine provides cubic interpolation of values
        C     of PSFs - based on the function AICC.
        C                                                 Version from 10-NOV-09
        C
              DIMENSION XX(4),YY(4),A(4)
              COMMON /PSFTAB/TABENPSF(3,0:400),TABPSF(3,0:400),NPSF(3)
        C
              APSF = 0.0
        c
        c    Low level densities treated in a special way
        c
              IF (EGX.LT.TABENPSF(MTYP,1)) GOTO 11            ! f(XL) = 0
              IF (EGX.GT.TABENPSF(MTYP,NPSF(MTYP))) GOTO 11   ! f(XL) = 0
        
              IF (EGX.LE.TABENPSF(MTYP,2)) THEN
                 K=0
                 ELSE
                 IF (EGX.LE.TABENPSF(MTYP,NPSF(MTYP)-1)) THEN
                    DO I=3,NPSF(MTYP)-1
                       IF (EGX.LE.TABENPSF(MTYP,I)) THEN
                          K=I-3
                          GO TO 1
                       ENDIF
                    ENDDO !I
                    ELSE
                    K=NPSF(MTYP)-4
                 ENDIF
              ENDIF
        C
            1 eglog=log(egx)
              DO J=1,4
                 XX(J)=log(TABENPSF(MTYP,K+J))
                 YY(J)=log(TABPSF(MTYP,K+J))
              ENDDO !J
              DO I=1,4
                 A(I)=YY(I)
                 DO J=1,4
                    IF (I.NE.J) A(I)=A(I)/(XX(J)-XX(I))
                 ENDDO !J
                 DO J=1,4
                    IF (I.NE.J) A(I)=A(I)*(XX(J)-EGLOG)
                 ENDDO !J
                 APSF=APSF+A(I)
              ENDDO !I
        C
              APSF=exp(APSF)
        
           11 CONTINUE
              RETURN
              END
        C
        