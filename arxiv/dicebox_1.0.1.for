c     version 1.0.1                                             11.3.2020
c - SMLO model corrected
c - NOPTM1.EQ.22 corrected

c - in EVENTS "S" is supra-realization; "R" is realization within supra-realization
c - variables renamed (NREAL -> NSREAL, NSUB -> NRL, ISUB -> IRL) in the code
c - checks of decay scheme - done
c - checks of the number of levels with the same spin and parity below ecrit (STDISa,...) - done
c - ERRPRIM - not used now - if changed, change manual

c - write if there is a level above Ecrit with no possible decay
c - treatment of RADW for TNC

c***********************************************************************
c     version 1.0                                             15.1.2018
c
c      - changes applied to avoid warnings from different compilers
c      - COMMON block /PHOTO/ is now consistent
c      - souroutine LEVELS updated - option LMODE=-1,-2
c
c      - individual levels can be printed into files LEVELS.nnn
c      - IF ((SNGL(TOTCON(0))+TOTDIS(0)).LE.0.) THEN - now writes down the cascade (GOTO 44 vs. GOTO 5)
c
c      - electron conversion should be treated fully correctly in quasicontinuum
c
c     - merged version for resonances and primaries with given fraction
c     - line in input file for LD=8 and 9; models in function DENSITY updated
c
c
c
c     Toshihiko's NLD added (LDTAB_K.DAT)
c     New version of GLO model added T=sqrt((BN-EG)/a)
c
c     Changes with respect to 2.85b
c        - model NOPTM1.EQ.25 and 26 added
c
c     Changes with respect to 2.85 (model NOPTE1 not included)
c        - model NOPTM1.EQ.31 added
c
c     Changes with respect to 2.84
c        - models added "Oslo-type" E1 models
c
c     Changes with respect to 2.83
c        - SLOPE for NOPTLD=11 added
c
c     Changes with respect to 2.82
c        - several PSF models added (11,67-69)
c
c
c  !!! IRCON EXTRA LARGE (for 243Am)
c     Changes with respect to 2.81
c        - MLO models updated - error in FKS() fixed
c        - M1 = 11 added (SC only at low energies + SP)
c
c     Changes made:
c        - subroutine GERMS changed - different ascription of precursors 
c          for extremely large level densities (NTOTAL > 30 000 000)
c        - two lines added to function RAN0() in order to prevent
c          RAN0 = 1.000       
c
c        - new models of LD and PSF added
c        - large dimension of IRCON
c     Changes with respect to 2.80:
c        - tabulated level density read in a different way (using LDTAB.DAT file)
c
c     Changes with respect to 2.79:
c        - some loops over real changed to loops over integer (gfortran warnings)
c        - some of the loops over low-lying levels are still only for SP = 0,8
c        - tabulated level density (Gorielly) added (NOPTLD = 11) 
c        - variables MAXJC,MAXBIN,MAXIS as PARAMETERS
c
c     Changes with respect to 2.78:
c        - New models added into SGAM
c     Changes with respect to 2.77:
c        - Wigner distribution of level spacings with "long-range"
c          correlations is working for LMODE=1
c     Changes with respect to 2.76:
c        - subroutine GERMS is changed in order to ascribe precursors
c          also in case of extreme number of levels (above 20 000 000)
c        - ICC in input data are in a different format: 
c          E1, M1, E2, M2, ... instead of E1, E2, E3, E4, M1,...
c          in previous versions
c     Changes with respect to 2.75:
c        - ICC on different shells distinguinsted in output
c        - if Eg<TABICC(1) then alpha(Eg)=0 ! 
c        - in version 2.75 might be a problem with ICC 
c          (different array bounds in CONV and TABICC)
c     Changes with respect to 2.70.5:
c        - Correlation matrix for population and sidefeeding of 
c          low-lying levels added -> changed COMMON /intr/
c        - Range for ICC enlarged to 5
c     Changes with respect to 2.61:
c        - Changes made into subroutines WIDTHS and ONESTEP to allow 
c          "random" energy of a level in a bin
c     Changes with respect to the version 2.60:
c        - Subroutine DENSITY updated (3 pars now) - new "parity-dependent"
c          level densities are allowed now
c        - input data file changed and the common block /DENP/ added 
c     Changes with respect to the version 2.54:
c        - Subroutine LEVELS updated (+2 new subroutines added) - long-range 
c          correlations in the Wigner distribution of level spacings 
c          taken into account
c     Changes with respect to the version 2.53:
c	 - added "more sophisticated" models into SGAM (30ff) + 
c          corresponding line	in input data added
c        - corrected model of Kopecky (EGLO) in SGAM
c        - checked models of Plujko
c     Changes with respect to the version 2.52:
c        - enlarged number of final levels (from 15 to 18)
c	 - covariance matrix written using FORMAT (<NFILEV>...)
c        - in function PEAKEFF the "zero efficiency" set to 0.0
c          instead of 10^(-4)
c	 - it is recommended to use formated output for covariance matrix,
c	   see lines about 620 
c
c     Changes with respect to the version 2.50:
c        - fixed problem with possible intensity < 0 in READ_INT   
C        - secondary intensities fluctuate from realiz. to realiz.
c          (procedure READ_INT is called in each realization)
c     Changes with respect to the version 2.42:
C        - random number generater RAN0 (from Numerical Recipes) is used
c        - changed some conditions in IF statements (simplyfying)
c        - fixed an inconsistence in INCREM caused due to the EC effect
c        - fixed problem with low level density at the beginning of
c          subroutine WIDTH
c
      PROGRAM DICE_EVENT
c
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
c
      LOGICAL        LOGS,ponverze,lopopgs,sidlev,ponverk,logfil
c      CHARACTER*80   OUTNAME
      CHARACTER*12   NAME
      CHARACTER*3    EXT
      INTEGER*4      IR1,IR2,IR3,IR4,STEPS,pointers,IRCON,MULT,NTOTAL,
     &               IRCONc,ISCON,ISDIS,IRINIT
      INTEGER        DEKOD,DENUM,DELEV,DEPARITY
      REAL*8         TOTCON,GACON,STCON
      DIMENSION      TOTDIS(0:2),TOTCON(0:2),
     &               STCON(0:2,0:MAXBIN,-2:2,0:1),GACON(0:2,-2:2,0:1),
     &               GADIS(0:2,-2:2,0:1),STDIS(0:2,0:30,-2:2,0:1),
     &               ISCON(0:2,0:MAXBIN,-2:2,0:1),
     &               ISDIS(0:2,20,0:MAXJC,0:1)
      DIMENSION      POPULT(99),POPERT(99),COVAP(99,99),CONO(99,99),
     &               POPULS(99),POPERS(99),COVAS(99,99),OLDPOP(99)
c
C     Hereafter, article 'CON' is related to a 'CONtinuum' part of
C     the level system, while 'DIS' -- to the 'DIScrete' part.
C     Suffixes 'p' and 's' are related to primary and secondary
C     transitions, respectively.
C
      COMMON /GAU/   U,IFLAG
     &       /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd,dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
     &       /incr/  nfilev,STEPS,enrgfin(18),buff(251),
     &               ncum,cumwidth(14),edo(18),eup(18)
     &       /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,NTOTAL,iregi,
     &               NOPTCS,NLINc,CAPFR(2)
     &       /PHOTO/ NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
     &               SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
     &               DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
     &               PAR_E1(3),DIPSLP,DIPZER
     &       /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
     &       /kon/   ponverze,ponverk
     &       /intr/  last,ir1,ir2,ir3,ir4,radwid,radwdi,
     &               fnoav,fnodi,rnn(18,14),
     &               rmult(18,126),rovfl(18),avea(18,14),erra(18,14),
     &               puvea(18,14),cova(18,18,14),aves(1:18,0:MAXBIN),
     &               errs(1:18,0:MAXBIN),
     &               diss(1:18,0:MAXBIN),poptlev(99),popslev(99)
     &       /readi/ nevents,NSREAL,NRL,ecrit,xrayk,xrayl,numlev,
     &               elowlev(99),elowsp(99),ilowip(99),
     &               STDISa(0:2,0:30,-2:2,0:1),RADW
     &       /grnd/  logs,kgs,LOpopGS,KpopGS
     &       /angle/ spinc(2),ipinc,spintl(18),ipintl(18),
     &               smear,Fk(4,4,0:16,0:16),F4(0:16,0:16)
     &       /monit/ ELQQ(0:126),SPQQ(0:126),IPQQ(0:126),ICQQ(126),
     &               DMQQ(126),WIQQ(0:126)
      COMMON /GOE/   EIGENVAL(0:MAXBIN,0:MAXJC,0:1),
     *               NEIGENVAL(0:MAXJC,0:1),LMODE
     &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
     &               ISWIC,ISWMX,ISWWI,ISWLS
     &       /rnd/   IRINIT(4,999)
c
c     lpreru=.FALSE.
      LOGS=.FALSE.
      lopopgs=.false.
      ponverze=.FALSE.
      ponverk=.FALSE.
      NDEAD=0
      DO J=0,MAXJC
        DO K=0,1
          NDIS(J,K)=0
        ENDDO
      ENDDO
      DO K=0,1
        DO L=0,20
          DO M=-2,2
            DO MM=0,2
              STDIS(MM,L,M,K)=0.
              STDISa(MM,L,M,K)=0.
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO K=1,99
        DO L=0,20
          SALL(K,L)=0.0
        ENDDO
        POPULT(K)=0.0
        POPERT(K)=0.0
        POPULS(K)=0.0
        POPERS(K)=0.0
        DO L=1,99
          COVAP(K,L)=0.0
          COVAS(K,L)=0.0
        ENDDO
      ENDDO
      RADWID=0.0
      RADWDI=0.0
c

cc    ! initial reading
c      CALL get_command_argument(1, NAME)
c      IF (LEN_TRIM(NAME).EQ.0) THEN
c          STOP 'Invalid input: need an input name as the first argument'
c      ENDIF
c      NAME=TRIM(NAME)
c      INQUIRE(FILE=NAME, EXIST=logfil)
c      IF (.not.logfil) THEN
c        STOP 'Invalid input: file not found'
c      ENDIF
c      CALL READ_EV(NAME)

      NAME='DICE_EV.DAT'
      CALL READ_EV(NAME)
      CALL RANDINIT()
C
C     Generation of GOE eigenvalues
c
      IF (LMODE.EQ.1) THEN
       NDIM=700              !Maximum allowed value is NDIM=1000
       NDIM=300              !Maximum allowed value is NDIM=1000
       CALL GENERATE_GOE_EIGEN_VAL(IR1,NDIM)
      ENDIF
c
c     PSFs used in tabular form
c
      OPEN (UNIT=12,FILE='PSF_GS.DAT',STATUS='UNKNOWN')
      OPEN (UNIT=13,FILE='PSF_INI.DAT',STATUS='UNKNOWN')

      EG_MAX = 11.0
      EG_STEP = 0.05
      NEG = INT (EG_MAX/EG_STEP)
      DO I = 1, NEG
        eg = FLOAT(I) * EG_STEP
        sf1_ini = sgamma(eg,bn,1)/eg**3     !from initial level
        sf2_ini = sgamma(eg,bn,3)/eg**3
        sf3_ini = sgamma(eg,bn,4)/eg**5
c        sf4_ini = sgamma(eg,bn,4)/eg**3
        sf1_fin = sgamma(eg,eg,1)/eg**3     !to GS
        sf2_fin = sgamma(eg,eg,3)/eg**3
        sf3_fin = sgamma(eg,eg,4)/eg**5
c        sf4_fin = sgamma(eg,eg,4)/eg**3
        write(12,201) eg,sf1_fin,sf2_fin,sf3_fin
        write(13,201) eg,sf1_ini,sf2_ini,sf3_ini
  201   format(f7.3,4e13.5)
      ENDDO
      CLOSE (12)
      CLOSE (13)
c
c          The loop through number of realisation of nucleus.
c
      DO NUC=1,NSREAL
c
c       write(*,*) NUC,IR1, IR2, IR3, IR4
        IR1 = IRINIT(1,NUC) !level scheme
        IR2 = IRINIT(2,NUC) !radiation widths in form of precursors and low-lying intensity fluctuations
        IR3 = IRINIT(3,NUC) !MC of cascades - actual search for the final state
        IR4 = IRINIT(4,NUC) !MC of cascades - the seed for a) mixing of primaries b) "coin flip" of mixing ratio
c
        CALL LEVELS (SPINC(1),IR1) ! Generation of level scheme
c
        IF (NTOTAL.GT.MAXIS) THEN
          WRITE (*,*)
     &      'The level density is too high. Increase IRCON dimension'
          WRITE (*,*) 'The required dimension is ', NTOTAL
          STOP 'ERROR'
        ENDIF
c
c       Writting levels 
c
        IF (ISWLS.EQ.1) THEN   ! Writing generated levels if requested
          WRITE (EXT,211) NUC
          S0=SPINC(1)-FLOAT(INT(SPINC(1)))
          IF (S0.LE.1e-3) THEN  ! only for security reason
           S0 = 0.0
          ELSE
           S0 = 0.5
          ENDIF
          OPEN (UNIT=13,FILE='LEVELS.'//EXT,STATUS='UNKNOWN')
          DO IS = 0, MAXJC
           SPFI = FLOAT(IS) + S0
           DO IPFI = 0, 1
            DO IBIN = 1, NBIN
              NL = LEVCON(IBIN,IS,IPFI)-LEVCON(IBIN-1,IS,IPFI)
              DO IL = 1, NL
                Q=1.-FLOAT(2*IL-1)/FLOAT(2*NL)
                ENRG = BN - (FLOAT(IBIN)-Q)*DELTA
                WRITE(13,212) ENRG, SPFI, IPFI
              ENDDO
            ENDDO ! IBIN
            NL = ndis(ISUBSC(spfi),ipfi)
            DO IL = 1, NL
              ENRG = endis(IL,ISUBSC(spfi),ipfi)
              WRITE(13,212) ENRG, SPFI, IPFI
            ENDDO
           ENDDO !IPFI
          ENDDO !IS
          CLOSE(13)
        ENDIF ! ISWLS
  211   FORMAT (I3.3)
  212   FORMAT(F8.5,F5.1,I2)  ! End of the part writing generated levels
c
        CALL GERMS(IR2) ! Precursors assigned to each level
c
c       Intensities of low-lying transitions flucutuate
c
        NAME='DICE_EV.DAT'
        CALL READ_INT(NAME)

        DO IRL = 1, NRL      !TODO pocitani pozorovatelnych
c
          IF (ISWWR.EQ.1) CALL OPEN_IT(NUC,IRL)
          WRITE(*,*) ' Supra-Realisation : ',NUC,' (of',NSREAL,')'
          WRITE(*,*) ' Realisation : ',IRL,' (of',NRL,')'

C
C      The following DO loop serves for computing of mixing ratios
C      delta for primary transitions (E2 admixture is probably not
C      very important, so this is done in a very simple way - and
C      probably not fully correctly in the case NLINc=2)
C
          DO ilinc=1,nlinc
           DO ip=0,1
            DO is=0,MAXJC
             DO il=1,NDIS(is,ip)
              DO ipom=1,100
               dummy=ran0(IR4)
              ENDDO
              isdis(ilinc,il,is,ip)=IR4
             ENDDO
            ENDDO
           ENDDO
          ENDDO
c
c     Here, the procedure 'WIDTHS' is called in order to evaluate a total
c     radiative width and proper values of STDISp.
c
c     Label 19 - Security reason! If TOTDIS(one of capt. state)>1 
c     (possible due to the random distribution of intensities between c.s.)
c

   19	    CONTINUE
          IF (IPRIM.EQ.1) THEN
            NAME='DICE_EV.DAT'
            IF (NLINc.GT.1) CALL READ_AGAIN(NAME)
          ENDIF
 
          IREGI=0
          IBIN=0
          ILIN=1
          RADW=0.
c
c
          IF (IPRIM.EQ.1) THEN ! Known primaries
           DO ILINc=1,NLINc
            CALL WIDTHS (ILINc,IPINC,SPINC(ILINc),IBIN,ILIN,TOTCON,
     &                   STCON,GACON,ISCON,TOTDIS,STDISa,GADIS,ISDIS)
           ENDDO
           DO ILINc=1,NLINc
            IF(TOTDIS(ILINc).GE.1.) GOTO 19
           ENDDO	 

           DO ILINc=1, NLINc
            RADWp=sngl(TOTCON(ILINc))/(1-TOTDIS(ILINc))
            DO IPFI=0,1
             SP=SPINC(ILINc)-INT(SPINC(ILINc)+.25)
             DO IS = 0, MAXJC
              SPFI = FLOAT(IS) + SP
              ISBS=NINT(SPFI+.25)-NINT(SPINC(ILINc)+.25)
              IF ((ISBS.GE.-2).AND.(ISBS.LE.2)) THEN
               DO L=1,NDIS(ISUBSC(SPFI),IPFI)
                STDIS(ILINc,L,ISBS,IPFI)=STDISa(ILINc,L,ISBS,IPFI)*RADWp
               ENDDO
              ENDIF
             ENDDO ! IS / SPFI
            ENDDO ! IPFI
            RADW=RADW+RADWp*CAPFR(ILINc)
           ENDDO ! ILINc
           IREGI=0
           IBIN=0
           ILIN=1
C
           DO ILINc=1,NLINc
            CALL WIDTHS(ILINc,IPINC,SPINC(ILINc),IBIN,ILIN,
     &           TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS)
           ENDDO

          ELSE
           DO ILINc=1,NLINc   ! Uknown primary intensities (resonances)
            CALL WIDTHS(ILINc,IPINC,SPINC(ILINc),IBIN,ILIN,
     &                TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS)
           ENDDO

           DO ILINc=1, NLINc
             RADW=RADW+sngl(TOTCON(ILINc))+TOTDIS(ILINc)
           ENDDO
           IREGI=0
           IBIN=0
           ILIN=1
          ENDIF !IPRIM


c
c   Population of low-lying levels - setting initial (zero) values
c
         DO I=1,NUMLEV
           POPTLEV(I)=0.0
           POPSLEV(I)=0.0
         ENDDO
C
c**************************************
c       The master DO-loop
c
         WRITE(*,*) ' Event (N=',NEVENTS,') : '
         WRITE(*,*)
         DO IEV=1,NEVENTS
          IF (MOD(IEV,20000).EQ.0) WRITE(*,5410) IEV
 5410     FORMAT('+',21X,I6)
          ponverze=.FALSE.
          ponverk=.FALSE.
          sidlev=.FALSE.
          EIN=BN
          STEPS=1
          pointers=0
C
C         Determination of spin of capture state
C
          ILINc=1
          IF (IVER.GE.1) THEN    ! Thermal capturing - more capturing states
           IF (NOPTCS.NE.1) THEN
            IF (RAN0(IR1).GT.CAPFR(1)) ILINc=2
           ENDIF
          ENDIF ! IVER

          ELQQ(0)=BN
          SPQQ(0)=SPINc(ILINc)
          IPQQ(0)=IPINc
C
          IREGI=0
          IBIN=0
          ILIN=1
          CALL ONESTEP(ILINc,IPINC,SPINC(ILINc),IBIN,ILIN,TOTCON,
     &                  STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS,
     &                  IPFI,SPFI,IBFI,ILFI,DMIX2,sign,IR3,IR4)
          DO WHILE (EFI.GT.0.)
            pointers=pointers+1
            ELQQ(steps)=EFI
            SPQQ(steps)=SPFI
            IPQQ(steps)=IPFI
            DMQQ(steps)=sign*sqrt(dmix2)
            IF (steps.EQ.1) THEN
              WIQQ(steps) = TOTCON(ILINc) + TOTDIS(ILINc)
            ELSE
              WIQQ(steps) = TOTCON(0) + TOTDIS(0)
            ENDIF
            IF (ELQQ(steps-1).LT.ECRIT) WIQQ(steps)=0.0
c            WRITE(*,*) ponverze
            IF (IREGI.GT.0) THEN
              do k=1,numlev
                if (efi.eq.elowlev(k)) then
                 poptlev(k)=poptlev(k)+1.
                 if (.NOT.sidlev) then
                  popslev(k)=popslev(k)+1.
                  sidlev=.TRUE.
                 endif
                endif
              enddo
            ENDIF
c
            if (ponverze) then
              if (ponverk) then
                ICQQ(steps)=1
              else
                ICQQ(steps)=2
              endif
            else
              ICQQ(steps)=0
            endif
            ponverze=.FALSE.
            ponverk=.FALSE.

c            WRITE (*,*) steps,buff(pointers),EIN,EFI,iregi,ilfi,ipfi,spfi
            IPIN=IPFI
            SPIN=SPFI
            IBIN=IBFI
            ILIN=ILFI
            IF (IREGI.LT.2) THEN
c              IF (IVER.GE.1) THEN ! Thermal neutron capture
                CALL WIDTHS(0,IPIN,SPIN,IBIN,ILIN,
     &               TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS)
c              ELSE
c                CALL WIDTHS_R(0,IPIN,SPIN,IBIN,ILIN,
c     &               TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS)
c              ENDIF ! IVER
              IF ((SNGL(TOTCON(0))+TOTDIS(0)).LE.0.) THEN
                NDEAD=NDEAD+1
                GOTO 44
              ENDIF
            ENDIF

            if (iregi.eq.2) then 
              if (denum(dekod(ilfi,isubsc(spfi),ipfi)).eq.0) then
                ndead=ndead+1
                go to 44
              endif
            endif
            STEPS=STEPS+1
            CALL ONESTEP(0,IPIN,SPIN,IBIN,ILIN,TOTCON,STCON,
     &        GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS,
     &        IPFI,SPFI,IBFI,ILFI,dmix2,sign,IR3,IR4)
          ENDDO  !WHILE EFI
          ELQQ(steps)=EFI
          SPQQ(steps)=SPFI
          IPQQ(steps)=IPFI
          DMQQ(steps)=sign*sqrt(dmix2)
          IF (steps.EQ.1) THEN
            WIQQ(steps) = TOTCON(ILINc) + TOTDIS(ILINc)
          ELSE
            WIQQ(steps) = TOTCON(0) + TOTDIS(0)
          ENDIF
          IF (ELQQ(steps-1).LT.ECRIT) WIQQ(steps)=0.0
C   **** feeding of ground state ****
          if (lopopgs) then
            poptlev(kpopgs)=poptlev(kpopgs)+1.
            if (.NOT.sidlev) then
             popslev(kpopgs)=popslev(kpopgs)+1.
c             sidlev=.TRUE.
            endif
          endif
c
            if (ponverze) then
              if (ponverk) then
                ICQQ(steps)=1
              else
                ICQQ(steps)=2
              endif
            else
              ICQQ(steps)=0
            endif
            ponverze=.FALSE.
            ponverk=.FALSE.
c
  44      IF (ISWWR.EQ.1) CALL DO_IT(STEPS)
c
c         WRITE (*,*) steps,buff(pointers),EIN,EFI,iregi
c
   5      CONTINUE
         ENDDO ! IEV
c*************************************
c
c*******    Cummulative quantities collected through the number
c                     of nuclear realizations :
c
         RNUC=FLOAT(IRL+(NUC-1)*NRL)
         RNUC1=RNUC-1.
         RNUC2=RNUC*RNUC
         RNUC12=RNUC1*RNUC1
c
c***    The total radiation width of the capturing state
c
         OLDWID=RADWID
         RADWID=(RNUC1*RADWID+RADW)/RNUC
         RADWDI=(RNUC1*(RADWDI+OLDWID*OLDWID)+RADW*RADW)/RNUC-
     &          RADWID*RADWID
c
c***    Population of low-lying levels (originaly nlowlev(level,nuc))
c
         DO K=1,NUMLEV
          POPTLEV(K)=POPTLEV(K)/FLOAT(NEVENTS)
          OLDPOP(K)=POPULT(K)
          POPULT(K)=(RNUC1*POPULT(K)+POPTLEV(K))/RNUC
          POPERT(K)=(RNUC1*(POPERT(K)+OLDPOP(K)*OLDPOP(K))+
     *         POPTLEV(K)*POPTLEV(K))/RNUC-POPULT(K)*POPULT(K)
          DO L=1,K
            COVAP(K,L)=(RNUC1*(COVAP(K,L)+OLDPOP(K)*OLDPOP(L))+
     &                 POPTLEV(K)*POPTLEV(L))/RNUC-POPULT(K)*POPULT(L)
          ENDDO
         ENDDO
         DO K=1,NUMLEV
          POPSLEV(K)=POPSLEV(K)/FLOAT(NEVENTS)
          OLDPOP(K)=POPULS(K)
          POPULS(K)=(RNUC1*POPULS(K)+POPSLEV(K))/RNUC
          POPERS(K)=(RNUC1*(POPERS(K)+OLDPOP(K)*OLDPOP(K))+
     *         POPSLEV(K)*POPSLEV(K))/RNUC-POPULS(K)*POPULS(K)
          DO L=1,K
            COVAS(K,L)=(RNUC1*(COVAS(K,L)+OLDPOP(K)*OLDPOP(L))+
     &                 POPSLEV(K)*POPSLEV(L))/RNUC-POPULS(K)*POPULS(L)
          ENDDO
         ENDDO
c
c***    Up-dating of the output files
c
         NAME='DICE.PRO'
         CALL WRITE_PARAMS(NAME)
         OPEN(9,FILE=NAME,ACCESS='APPEND',STATUS='OLD')
          WRITE(9,102) NOPTDE,NOPTE1,NOPTM1,NOPTE2
  102     FORMAT(' NOPTDE:',I2,'   NOPTE1:',I2,'   NOPTM1:',I2,
     &         '   NOPTE2:',I2)
          WRITE(9,*)
          WRITE(9,*) ' Real #',NUC,' of',NSREAL,
     &    '      Subr #',IRL,' of',NRL
          WRITE(9,*) ' Number of events',NEVENTS
          WRITE(9,*)
          ERR=SQRT(RADWDI)
          WRITE(9,*) 'Capt.state tot.rad.width (MeV): ',
     &              RADWID,' +/-',ERR
          WRITE (9,*)
          WRITE (9,*)
     &    'Number of cascades terminating at a dead end: ',NDEAD
c     ***  feeding of low-lying states ***
          WRITE(9,*) 'Population of low-lying states'
          if (rnuc1.GT.0) then   
	      rn1=rnuc/rnuc1
          else
	      rn1=1.
	    endif
          do i=1,numlev
           if (popert(i).GT.0.0) then
            poper=sqrt(popert(i)*rn1)
           else
            poper=0.0  
           endif
           WRITE(9,198)i,elowlev(i),popult(i),poper,
     &                elowsp(i),ilowip(i)
          enddo
  198     format(I3,f10.5,F11.5,' +- ',F7.5,F11.1,I3)
          write(9,*)
          WRITE(9,*)
     *        'Direct population of low-lying states from continuum'
          do i=1,numlev
           if (popers(i).GT.0.0) then
            poper=sqrt(popers(i)*rn1)
           else
            poper=0.0  
           endif
           WRITE(9,198)i,elowlev(i),populs(i),poper,
     &                 elowsp(i),ilowip(i)
          enddo
          write(9,*)
c
          WRITE(9,*) ' Population - covariance matrix'
          DO K=1,NUMLEV
           DO L=1,K
            CON=COVAP(K,K)*COVAP(L,L)
            IF (CON.GT.0.) THEN
              CONO(K,L)=COVAP(K,L)/SQRT(CON)
            ELSE
              CONO(K,L)=1.0
            ENDIF
           ENDDO
          ENDDO
          DO K=1,NUMLEV
c          WRITE(9,105) (CONO(K,L),L=1,NUMLEV)
       !    WRITE(9,*) (CONO(K,L),L=1,NUMLEV)   ! No covariance matrix written in the present version 
c  105     FORMAT(<NUMLEV>F7.3)
          ENDDO
          WRITE(9,*)

          WRITE(9,*) ' Sidefeeding - covariance matrix'
          DO K=1,NUMLEV
           DO L=1,K
            CON=COVAS(K,K)*COVAS(L,L)
            IF (CON.GT.0.) THEN
              CONO(K,L)=COVAS(K,L)/SQRT(CON)
            ELSE
              CONO(K,L)=1.0
            ENDIF
           ENDDO
          ENDDO
          DO K=1,NUMLEV
c          WRITE(9,105) (CONO(K,L),L=1,NUMLEV)
      !     WRITE(9,*) (CONO(K,L),L=1,NUMLEV)   ! No covariance matrix written in the present version 
          ENDDO

c        WRITE(9,*)
c        WRITE(9,*) 'last generator seed IR1: ',IR1
c        WRITE(9,*) '                    IR2: ',IR2
c        WRITE(9,*) '                    IR3: ',IR3
c        WRITE(9,*) '                    IR4: ',IR4

         CLOSE(9)
         CALL CLOSE_IT
c
c       write(*,*) NUC, IR1, IR2, IR3, IR4, 'l'

         IRCONc(1) = IR3 ! Is the OK - reproducibility?
         IRCONc(2) = IR3 - 1 ! Is the OK - reproducibility?
         IF (IRCONc(2).LE.0)  IRCONc(2) = IR3 + 1

        ENDDO ! IRL
c   55 CONTINUE
      ENDDO ! NUC
c
c******************************************************************
      STOP 'Execution of DICEBOX ended O.K.'
      END
c
C**********************************************************************
      FUNCTION ALPH_TOT (EI,SPI,IPI,EF,SPF,IPF,DMISQ,NEN,ELEN,CONV)
C**********************************************************************
C
      DIMENSION ELEN(50),CONV(0:1,5,50)
C
C           DMISQ   - SQUARED (!) mixing amplitude
C
C           MAEL0   - the type of the lowest order of gamma radiation.
C                     for MAgnetic radiation MAEL0=1, for ELectric
C                     MAEL0=0
C
C           MUL0    - the lowest multipolarity contributing to the
C                     transition (say, 2 in case of E2+M3)
C
C           MUL1    - the next contributing multipolarity
C                     (3 in this case)
C
      MUL0=NINT(ABS(SPI-SPF))
      IF (MUL0.EQ.0) MUL0=1
      MUL1=MUL0+1
      IAUX=(-1)**(IPI+IPF+MUL0+1)
      IF (IAUX.EQ.-1) THEN
         MAEL0=0    ! The dominating type is electric
         ELSE
         MAEL0=1    ! ... magnetic
      ENDIF
C
C     MAEL1 - the type  of the next-order contributing radiation
C
      MAEL1=1-MAEL0
      EG=EI-EF
      CT0=AICC(EG,ELEN,CONV,MAEL0,MUL0,NEN)
      CT1=AICC(EG,ELEN,CONV,MAEL1,MUL1,NEN)
c
      ALPH_TOT=(CT0+CT1*DMISQ)/(1.+DMISQ)
      RETURN
      END
C

C***********************************************************************
      FUNCTION ITYPE(SPIN,IPIN,SPFI,IPFI)
C
C     The function ITYPE determines the type of gamma-transition.
C     The meaning of the variables is evident.
C      ITYPE=1: Pure E1-transition
C            2: Mixed (M1+E2)-transition
C            3: Pure M1-transition
C            4: Pure E2-transition
C     This function also checks wheather the spin of the final state
C     falls within the limits 0. to 50. or -- in case of odd product
C     nuclei -- within the limits 0.5 to 50.5.
C
C***********************************************************************
C
      JL=NINT(ABS(2.*SPIN-2.*SPFI))
      JU=NINT(2.*SPIN+2.*SPFI)
      IF((SPFI.LT.(-.25)).OR.(SPFI.GT.(50.75))) GOTO 3
      IF ((JL.LE.2).AND.(JU.GE.2)) GOTO 1
      IF ((JL.LE.4).AND.(JU.GE.4)) GOTO 2
    3 ITYPE=0
      RETURN
    2 IF (IPIN.NE.IPFI) GOTO 3
      ITYPE=4
      RETURN
    1 IF ((JL.LE.4).AND.(JU.GE.4)) GOTO 4
      IF (IPIN.EQ.IPFI) GOTO 5
    6 ITYPE=1
      RETURN
    5 ITYPE=3
      RETURN
    4 IF (IPIN.NE.IPFI) GOTO 6
      ITYPE=2
      RETURN
      END
C
C
C

C***********************************************************************
c     version 2.61                                         19.1.2005              
      FUNCTION DENSITY(EEXC,SPIN,IPAR)
C
C     Explicit expressions for level density
C     NOPTDE=0: CTF-model
C           =1: Bethe's level-density formula following formulation
C               of T.von Egidy et al., Nucl.Phys. A (1988)
c           =2: modified BSFG
c           =3: modified CTF
c         =4,5: BSFG with modified spin cut-off parameter
C
c     Changed a factor 1/2 in the row DENSITY=DENSITY*FJ*.5
c          => the parity-dependent level density allowed (PRC67,015803)
c
C***********************************************************************
C
      COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
     *       /denp/  LDENP,ZNUM,DENPC,DENPEN(3),IPDEN
C
      DENSITY=0.
      IF (NOPTDE.EQ.0) THEN                                 ! CTF
        EEFF=EEXC-EZERO
        IF (EEFF.LE.0.) RETURN
        DENSITY=EXP(EEFF/TEMPER)/TEMPER
        SIGSQ=(.98*AMASS**.29)**2.
c        SIGSQ=(2.*AMASS**.29)**2.
      ELSEIF (NOPTDE.EQ.1) THEN                             ! BSFG
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0888*SQRT(ASHELL*EEFF)*AMASS**.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)
      ELSEIF (NOPTDE.EQ.2) THEN                             ! modified BSFG
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0888*SQRT(ASHELL*EEFF)*AMASS**.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)
        IF ((EEXC.GE.DENLO).AND.(EEXC.LE.DENHI)) THEN
         FCTDEN=DENPAR(1)+DENPAR(2)*EEXC+
     +          DENPAR(3)*EEXC**2+DENPAR(4)*EEXC**3
         DENSITY=DENSITY*FCTDEN
        ENDIF
      ELSEIF (NOPTDE.EQ.3) THEN                              ! modified CTF
        EEFF=EEXC-EZERO
        IF (EEFF.LE.0.) RETURN
        DENSITY=EXP(EEFF/TEMPER)/TEMPER
        SIGSQ=(.98*AMASS**.29)**2.
        IF ((EEXC.GE.DENLO).AND.(EEXC.LE.DENHI)) THEN
         FCTDEN=DENPAR(1)+DENPAR(2)*EEXC+
     +          DENPAR(3)*EEXC**2+DENPAR(4)*EEXC**3
         DENSITY=DENSITY*FCTDEN
        ENDIF
      ELSEIF (NOPTDE.EQ.4) THEN                             ! BSFG, s=0.1446 (Paar,Al-Quraishi)
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.1446*SQRT(ASHELL*EEFF)*AMASS**.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)
      ELSEIF (NOPTDE.EQ.5) THEN                             ! BSFG, another s (Al-Quraishi)
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0145*0.8*SQRT(EEFF/ASHELL)*AMASS**1.66667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)

      ELSEIF (NOPTDE.EQ.6) THEN                             ! BSFG - Von Egidy cut-off
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0146*(1+SQRT(1+4*ASHELL*EEFF))/2./ASHELL*AMASS**1.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)

      ELSEIF (NOPTDE.EQ.8) THEN                                 ! CTF - Von Egidy 2009
        EEFF=EEXC-EZERO09
        IF (EEFF.LE.0.) RETURN
        IF (EEXC.LT.0.5*PAIRING09) RETURN
        SIGSQ=.391*AMASS**0.675*(EEXC-0.5*PAIRING09)**0.312 
        DENSITY=EXP(EEFF/TEMPER09)/TEMPER09
      ELSEIF (NOPTDE.EQ.9) THEN                             ! BSFG - Von Egidy 2009
        EEFF=EEXC-EONE09
        IF (EEFF.LE.0.) RETURN
        IF (EEXC.LT.0.5*PAIRING09) RETURN
        SIGSQ=.391*AMASS**0.675*(EEXC-0.5*PAIRING09)**0.312 
        DENSITY=EXP(2.*SQRT(ASHELL09*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL09**.25*EEFF**1.25)

      ELSEIF (NOPTDE.EQ.11) THEN                            ! Goriely
        IF (EEXC.LE.0.) RETURN
        DENSITY = ALD(EEXC,SPIN,IPAR)
        GOTO 31
      ELSEIF (NOPTDE.EQ.12) THEN                            ! Kawano
        IF (EEXC.LE.0.) RETURN
        DENSITY = ALD(EEXC,SPIN,IPAR)
        GOTO 31

      ELSEIF (NOPTDE.EQ.7) THEN                             ! Voinov Mo BSFG
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0146*SQRT(EEFF/ASHELL)*AMASS**1.666667
        DENSITY=EXP(2.*SQRT(ASHELL*EEFF))/(16.9706*SQRT(SIGSQ)*
     *    ASHELL**.25*EEFF**1.25)

      ELSEIF (NOPTDE.EQ.13) THEN                             ! BSFG - Von Egidy cut-off
        EEFF=EEXC-EONE
        IF (EEFF.LE.0.) RETURN
        SIGSQ=.0146*(1+SQRT(1+4*ASHELL*EEFF))/2./ASHELL*AMASS**1.666667
        EEFF=EEXC-EZERO
        DENSITY=EXP(EEFF/TEMPER)/TEMPER

      ENDIF
c
      FJ=(SPIN+.5)*EXP(-(SPIN+.5)**2/(2.*SIGSQ))/SIGSQ
      IF ((NOPTDE.EQ.8).OR.(NOPTDE.EQ.9)) THEN
        IF (MOD(INT(AMASS+0.25),2).EQ.0) THEN
          IF (MOD(INT(ZNUM+0.25),2).EQ.0) THEN
            STAG = (EEXC - DENLO) / (DENHI - DENLO)
            IF (STAG.LE.0.0) STAG = 0.0
            IF (STAG.GE.1.0) STAG = 1.0
            IF (SPIN.LT.0.25) THEN
              FJ = FJ * (1.0 + 1.02 * (1.0 - STAG) )
            ELSEIF (MOD(INT(SPIN+0.25),2).EQ.0) THEN
              FJ = FJ * (1.0 + 0.227 * (1.0 - STAG) )
            ELSE
              FJ = FJ * (1.0 - 0.227 * (1.0 - STAG) )
            ENDIF
          ENDIF
        ENDIF
      ENDIF

C
C     "Parity-dependence" term
C
      IF (LDENP.EQ.0) THEN
        PARDEP=0.5
      ELSEIF (LDENP.EQ.2) THEN
        PAIRS=DENPEN(1)+DENPEN(2)/AMASS**DENPEN(3)    !see PRC67, 015803
      ELSEIF (LDENP.EQ.1) THEN                        !see PRC67, 015803
        IF (MOD(INT(AMASS+0.25),2).EQ.0) THEN
          IF (MOD(INT(ZNUM+0.25),2).EQ.0) THEN
            PAIRS= 1.34+75.22/AMASS**0.89             !E-E nucleus
          ELSE
            PAIRS=-0.90+75.22/AMASS**0.89             !O-O nucleus
          ENDIF
        ELSE
          IF (MOD(INT(ZNUM+0.25),2).EQ.0) THEN
            PAIRS=-0.08+75.22/AMASS**0.89             !E-O nucleus
          ELSE
            PAIRS=-0.42+75.22/AMASS**0.89             !O-E nucleus
          ENDIF
        ENDIF
      ENDIF
      IF (LDENP.GT.0) THEN
        IF (IPDEN.EQ.0) THEN 
          PARDEP=0.5*(1+1/(1+EXP(DENPC*(EEXC-PAIRS))))
        ELSE
          PARDEP=0.5*(1-1/(1+EXP(DENPC*(EEXC-PAIRS))))
        ENDIF
      ENDIF           
      IF (IPAR.EQ.1) PARDEP=1.0-PARDEP
c
      DENSITY=DENSITY*FJ*PARDEP
C
   31 CONTINUE
      RETURN
      END
c


c This file contain subroutines
c   LABELS - convert numeric value of give model to string label
c   WRITE_PARAMS - write down (to the output file) parameters concerning
c                  giant resonances, SP strength, Fermi liq. param., and
c                  values Ek0 and Eg0 (according to Kopecky)
c
c************************************************************************
      SUBROUTINE LABELS(NOPTDE,NOPTE1,NOPTM1,NOPTE2,tdens,tsfe1,tsfm1,
     &                      tsfe2)
c************************************************************************
        CHARACTER*8 tdens,tsfe1,tsfm1,tsfe2
c
c    Conversion number of model to string label of the model
c
      tdens='???'
      tsfe1='???'
      tsfm1='???'
      tsfe2='???'

      IF (noptde.EQ.0) THEN
        tdens='CTF'
      ELSEIF (noptde.EQ.1) THEN
        tdens='BSFG'
      ENDIF

      IF (nopte1.EQ.0) THEN
        tsfe1='SP'
      ELSEIF (nopte1.EQ.1) THEN
        tsfe1='BA'
      ELSEIF (nopte1.EQ.2) THEN
        tsfe1='TD-BA'
      ELSEIF (nopte1.EQ.3) THEN
        tsfe1='KMF+Chr'
      ELSEIF (nopte1.EQ.4) THEN
        tsfe1='KMF'
      ELSEIF (nopte1.EQ.5) THEN
        tsfe1='Chr'
      ELSEIF (nopte1.EQ.6) THEN
        tsfe1='Chr-phD'
      ELSEIF (nopte1.EQ.7) THEN
        tsfe1='SMLO'
      ELSEIF (nopte1.EQ.8) THEN
        tsfe1='8'
      ELSEIF (nopte1.EQ.9) THEN
        tsfe1='9'
      ELSEIF (nopte1.EQ.10) THEN
        tsfe1='10'
      ELSEIF (nopte1.EQ.11) THEN
        tsfe1='11'
      ENDIF

      IF (noptm1.EQ.0) THEN
        tsfm1='SP'
      ELSEIF (noptm1.EQ.1) THEN
        tsfm1='BA'
      ELSEIF (noptm1.EQ.2) THEN
        tsfm1='BAonSP'
      ELSEIF (noptm1.EQ.3) THEN
        tsfm1='Pow'
      ELSEIF (noptm1.EQ.4) THEN
        tsfm1='4'
      ENDIF

      IF (nopte2.EQ.0) THEN
        tsfe2='SP'
      ELSEIF (nopte2.EQ.1) THEN
        tsfe2='BA'
      ENDIF

      END


C***********************************************************************
      SUBROUTINE LEVELS(SPC,IR)
C***********************************************************************
C
C    - no fluctuations at all LMODE=-1
C    - the Poisson distribution of neighbourhood level spacing is
C      assumed for LMODE=0 
C    - the Wigner distribution (with long-range correlations) is
C      assumed for LMODE=1 
c    - the "restricted Wigner distribution" - no long-range correlations
C
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
c
      integer   dekod,denum,delev,deparity
      INTEGER*4 IR,IRCON,IRCONc,NTOTAL
      REAL*4    EIN,EFI
C
      COMMON  /WID/  sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
     *       /MAINE/ eall,EIN,EFI,noptfl,BN,DELTA,NBIN,NTOTAL,iregi,
     *               noptcs,nlinc,capfr(2)
      COMMON /GOE/   EIGENVAL(0:MAXBIN,0:MAXJC,0:1),
     *               NEIGENVAL(0:MAXJC,0:1),LMODE
      COMMON  /GAU/   U,IFLAG
C
C     (Variables printed in lower case are not used in this subroutine)
C
      NTOTAL=0
      IFLAG=0

      IF (LMODE.EQ.0) THEN              !Poisson distribution
       DO IP=0,1
        SP=SPC-INT(SPC+.25)-1.
        DO J=0,MAXJC
         SP=SP+1.
         DO I=1,NBIN
          E=BN+DELTA/2.-FLOAT(I)*DELTA
          AVNL=DELTA*DENSITY(E,SP,IP)
          LEVCON(I,J,IP)=NPOISS(IR,AVNL)
          NTOTAL=NTOTAL+LEVCON(I,J,IP)
         ENDDO !I
        ENDDO !J
       ENDDO !IP
      ELSEIF (LMODE.EQ.1) THEN      !Wigner distribution with long-range correlations
       DO IP=0,1
        SP=SPC-FLOAT(INT(SPC+.25))-1.
        DO J=0,MAXJC
         SP=SP+1.
  	   X=10.*RAN0(IR)  ! Only to "randomize" the beginning of the sequence
	   K=1
    8	   IF (GOE_EIGEN_VAL(K,J,IP).LT.X) THEN
          K=K+1
          GOTO 8
	   ENDIF
         KLO=K
         DO I=NBIN,1,-1
          E=BN+DELTA/2.-FLOAT(I)*DELTA
          AVNL=DELTA*DENSITY(E,SP,IP)
	    X=X+AVNL
    9     IF (GOE_EIGEN_VAL(K,J,IP).LT.X) THEN 
           K=K+1
           GOTO 9
          ENDIF   
          KHI=K
          LEVCON(I,J,IP)=KHI-KLO        
          NTOTAL=NTOTAL+KHI-KLO
          KLO=KHI
         ENDDO !I
        ENDDO !J
       ENDDO !IP
      ELSEIF (LMODE.EQ.2) THEN              !Wigner distribution - no long-range correlations
       DO IP=0,1
        SP=SPC-INT(SPC+.25)-1.
         DO J=0,MAXJC
          DO I=0,NBIN
           LEVCON(I,J,IP)=0
          ENDDO  !I
          SP=SP+1.
          OMEGA=0.
          ALPHA=0.
          I=0
C
C         OMEGA is a random sample drawn from the Wigner distri-
C         bution; the expectation value of average distance between
C         neighbouring levels of a given spin and parity is assumed
C         to be equal to 1.
C
C           The constant 1.1283791671 is equal to "two divided by
C           square root of pi"
C
    1     RA=RAN0(IR)
          IF (RA.LE.0.) GO TO 1
          OMEGA=OMEGA+1.1283791671*SQRT(-ALOG(RA))
    3     IF (OMEGA.LT.ALPHA) THEN
            LEVCON(I,J,IP)=LEVCON(I,J,IP)+1
            NTOTAL=NTOTAL+1
            GO TO 1
          ELSE
            I=I+1
            IF (I.GT.NBIN) GO TO 2
            E=BN+DELTA/2.-FLOAT(I)*DELTA
            ALPHA=ALPHA+DENSITY(E,SP,IP)*DELTA
            GO TO 3
          ENDIF
    2     CONTINUE
         ENDDO !J
       ENDDO !IP
      ELSEIF (LMODE.EQ.-1) THEN             ! cumulative number of levels used; based on Wigner distribution but no randomnes
       DO IP=0,1
        SP=SPC-INT(SPC+.25)-1.
         DO J=0,MAXJC
          DO I=0,NBIN
           LEVCON(I,J,IP)=0
          ENDDO  !I
          SP=SP+1.
          OMEGA=0.
          ALPHA=0.
          I=0
C
   11     CONTINUE 
          OMEGA=OMEGA+1.0
   13     IF (OMEGA.LT.ALPHA) THEN
            LEVCON(I,J,IP)=LEVCON(I,J,IP)+1
            NTOTAL=NTOTAL+1
            GO TO 11
          ELSE
            I=I+1
            IF (I.GT.NBIN) GO TO 12
            E=BN+DELTA/2.-FLOAT(I)*DELTA
            ALPHA=ALPHA+DENSITY(E,SP,IP)*DELTA
            GO TO 13
          ENDIF
   12     CONTINUE
         ENDDO !J
       ENDDO !IP
      ELSEIF (LMODE.EQ.-2) THEN                  ! No fluctuations
       DO IP=0,1
        SP=SPC-INT(SPC+.25)-1.
        DO J=0,MAXJC
         SP=SP+1.
         DO I=1,NBIN
          E=BN+DELTA/2.-FLOAT(I)*DELTA
          AVNL=DELTA*DENSITY(E,SP,IP)
          LEVCON(I,J,IP)=NINT(AVNL)
          NTOTAL=NTOTAL+LEVCON(I,J,IP)
         ENDDO !I
        ENDDO !J
       ENDDO !IP


      ENDIF
C
C     NTOTAL is the total number of generated levels.
C
C     At this moment for each value of I the variable LEVCON(...) contains
C     the number of those levels of a particular spin and parity that fall
C     within the corresponding energy bin (whose width is DELTA).
C     The following DO-loops, however, convert this differential distribution
C     of level energies into a CUMULATIVE (i.e. integral) form. This simple
C     conversion leads to a significant increase of the speed of
C     function SEED.
C
    5 IAUX=0
      DO IP=0,1
       DO J=0,MAXJC
        LEVCON(0,J,IP)=IAUX
        DO I=1,NBIN
          LEVCON(I,J,IP)=LEVCON(I,J,IP)+LEVCON(I-1,J,IP)
        ENDDO 
        IAUX=LEVCON(NBIN,J,IP)
       ENDDO !J
      ENDDO !IP
      RETURN
      END
C


C***********************************************************************
      SUBROUTINE GERMS(IR)
C
C     This subroutine generates a random-generator random seed for
C     each individual level. The seeds obtained are stored in IRCON(K).
C     For a level of interest, characterized by {IP,SP,IB,IL}, the
C     corresponding seed can be fetched in a simple way. This is evident
C     from the body of the function SEED that performs such an operation.
C
C***********************************************************************
c
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      integer   dekod,denum,delev,deparity
      INTEGER*4 IR,IRCON,IRCONc,I,IG,K,N,NTOTAL,IAUX(8000000),
     &          NOLD,NDIF,NMGEN
      REAL*4  EIN,EFI
C
      COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
     *        /MAINE/ eall,EIN,EFI,noptfl,bn,delta,NBIN,NTOTAL,iregi,
     *                NOPTCS,NLINc,capfr(2)
C
C     Seed for 1st capturing state
C
      k=1
      NOLD = 0
      IRCONc(k)=IR
C
C     IRCON(.) is seeded in a random fashion by seeds IR produced
C     consecutively by repeated calls of RAN(IR), see the body of the
C     routine. If two or more seeds are falling to the same word of
C     IRCON(.), then only the last of them is kept there.
C
      N = NTOTAL
      DO 1 K=1,NTOTAL+nddd
    1 IRCON(K)=0

      IF (N.GT.30000000) THEN
        KLAST = 1
        NMAXR = 0
    9   KLAST = KLAST + NMAXR 
        NMAXR = INT( N / 2) - 2
        DO K = KLAST, KLAST + NMAXR 
          DUMMY = RAN0(IR)
          IRCON(K) = IR
        ENDDO
        N = N - NMAXR
        IF (N.GT.15000000) GOTO 9 
      ENDIF

    4 CONTINUE
      IF (N.GT.15000000) THEN
        NRNDN = 100000000
      ELSE
        NRNDN = 10000000
      ENDIF
      DO 2 IG=1,NRNDN
      K=int(float(NTOTAL+nddd)*RAN0(IR))+1
      DUMMY=RAN0(IR)
    2 IRCON(K)=IR
C
C     Determine the number of not seeded sites N
C
      N=0
      DO 3 K=1,NTOTAL+nddd
      IF (IRCON(K).EQ.0) N=N+1
    3 CONTINUE

      NDIF = ABS(NOLD - N)
      NOLD = N
      IF ((N.GT.2000).AND.(NDIF.GT.1000)) GOTO 4
      IF (N.GT.2000) THEN
        IF (N.GT.8000000) THEN
          WRITE (*,*)
     &      'The level density is too high.'
          WRITE (*,*) ' '
          WRITE (*,*) 'The required dimension in subr GERMS is ', N
          STOP 'ERROR'
        ENDIF
      ENDIF
C
C     If N is too high (>2000), the WHOLE field is to be randomly
C     seeded with another handful of 10**6 seeds.
C
    7 IF (N.EQ.0) THEN
       IF (NOPTCS.NE.1) THEN
        DUMMY=RAN0(IR)
        k=2
        IRCONc(k)=IR         ! 2nd capture state is now seeded
       ENDIF
       RETURN
      ENDIF
C     If every site of the useful part of IRCON(.) is seeded,
C     the RETURN follows. If not, ONLY EMPTY WORDS  are additionally
C     seeded. A smaller number (100000) of seeds is used. Seeds
C     are again being distributed randomly. This new seeding to
C     a restricted area is to be repeated up to the point when ALL
C     words are seeded. After that the RETURN follows.
C
      N=0
      DO 5 K=1,NTOTAL+nddd
      IF (IRCON(K).NE.0) GOTO 5
      N=N+1
      IAUX(N)=K
    5 CONTINUE

      IF (N.LT.10000) THEN
        NMGEN = 100000
      ELSE
        NMGEN = 1000000
      ENDIF

      DO 6 IG=1,NMGEN
      I=int(float(N)*RAN0(IR))+1
      DUMMY=RAN0(IR)
    6 IRCON(IAUX(I))=IR
      GOTO 7
      END
C
C
C***********************************************************************
      SUBROUTINE WIDTHS (MODE,IPIN,SPIN,IBIN,ILIN,
     1             TOTCON,STCON,GACON,ISCON,TOTDIS,STDIS,GADIS,ISDIS)
C***********************************************************************
c    Diferrent treatment of primary transitions to low-lying levels
C
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      LOGICAL    renergy
      integer    dekod,denum,delev,deparity
      INTEGER*4  IRCON,IRCONc,ISCON,ISEED,SEEDS,ISDIS,ntotal
      REAL       Im2Res
      REAL*4     EIN,EFI
      REAL*8     TOTCON,STCON,GACON
      DIMENSION  TOTDIS(0:2),TOTCON(0:2),STCON(0:2,0:MAXBIN,-2:2,0:1),
     *           ISCON(0:2,0:MAXBIN,-2:2,0:1),GACON(0:2,-2:2,0:1),
     *           STDIS(0:2,0:30,-2:2,0:1),GADIS(0:2,-2:2,0:1),
     *           ISDIS(0:2,20,0:MAXJC,0:1),S(2),GG(2)
      COMMON  /WID/  sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd,dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
     &       /GAU/   u,IFLAG,
     &       /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,ntotal,iregi,
     &               noptcs,nlinc,capfr(2)
     &       /SFCE/  SFCEE1,SFCEM1,SFCEE2,GE1,GM1,GE2,
     &               GE1DIS(2),GM1DIS(2),GE2DIS(2)
     &       /enh/   sumg,ibmin,ibmax
     &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)
     &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
     &               ISWIC,ISWMX,ISWWI,ISWLS
     &       /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
     &               CONVk(0:1,5,50)

C
C     Variables printed in lower case are not EXPLICITELY used in
C     this subroutine
C
C     If you use MODE=1 or 2 you have to set IBIN=1 and ILIN=1 !
C
      GE1=0.
      GM1=0.
      GE2=0.

      IF (IREGI.EQ.0) THEN
        IF (MODE.NE.0) THEN
          EIN=BN
          Q=0.0
        ELSE
          NLEV=LEVCON(IBIN,ISUBSC(SPIN),IPIN)- 
     *         LEVCON(IBIN-1,ISUBSC(SPIN),IPIN)
          Q=1.-FLOAT(2*ILIN-1)/FLOAT(2*NLEV)
          EIN=BN-(FLOAT(IBIN)-Q)*DELTA
        ENDIF
        IF (NOPTFL.GE.1) ISEED=SEEDS(MODE,SPIN,IPIN,ILIN,IBIN)
      ELSE
        EIN=EFI
      ENDIF
      SPAC=1./DENSITY(EIN,SPIN,IPIN)
      IF (SPAC.LE.0.0) SPAC=0.00001
      TOTCON(MODE)=0.D+0
      TOTDIS(MODE)=0.
c
      DO IS=-2,2
       DO IP=0,1
        DO I=0,NBIN
         STCON(MODE,I,IS,IP)=0.D+0
        ENDDO
        IF (MODE.EQ.0) THEN
         DO I=0,20
          STDIS(MODE,I,IS,IP)=0.
         ENDDO
        ENDIF
        GACON(MODE,IS,IP)=0.D+0
        GADIS(MODE,IS,IP)=0.
       ENDDO
      ENDDO
C
      SP  = SPIN - INT(SPIN + .25)
      ISP = INT(SPIN + .25)
      DO 1 IPFI=0,1
      DO 1 ISPFI= ISP-2, ISP+2
      SPFI = SP + FLOAT(ISPFI)
      ISBS=NINT(SPFI+.25)-NINT(SPIN+.25)
C
      IT=ITYPE(SPIN,IPIN,SPFI,IPFI)
      IF (IT.EQ.0) GOTO 1
      IF (IT.EQ.2) THEN
         IT1=3
         IT2=4
         ELSE
         IT1=IT
         IT2=IT
      ENDIF
      STCON(MODE,IBIN,ISBS,IPFI)=0.D+0
      IF (iregi.eq.1) then
       ISEED=SEEDS(mode,spin,ipin,ilin,ibin)
       GOTO 24
      ENDIF
C
      DO I=IBIN+1,NBIN
       NL=LEVCON(I,ISUBSC(SPFI),IPFI)-LEVCON(I-1,ISUBSC(SPFI),IPFI)
       EG=(FLOAT(I-IBIN-1)+Q+.5)*DELTA    ! Q=0.5 ... fixed; trans. between mid-bins
       Z=0.
       DO ITT=IT1,IT2
        S(ITT-IT1+1)=SGAMMA(EG,EIN,ITT)
       ENDDO

       IF (NOPTFL.NE.1) THEN  !No fluctuation
         ZZ = 0.0
         GG(2)=0.0
         DO ITT=IT1,IT2
          GG(ITT-IT1+1) = S(ITT-IT1+1)
          ZZ = ZZ + GG(ITT-IT1+1)
         ENDDO
         IF (GG(1).GT.0.0) THEN
          DMIX2 = GG(2) / GG(1)
         ELSE
          DMIX2 = 0.0
         ENDIF
         alpha = ALPH_TOT(EG,spfi,ipfi,EG-EG,
     *                    spin,ipin,dmix2,nent,elent,convt)
         Z = Z + ZZ * FLOAT(NL) * (1 + alpha)
       ELSE
        ISCON(MODE,I,ISBS,IPFI)=ISEED
        IFLAG=0
        IF (MODE.EQ.0) THEN
         DO IL=1,NL
          ZZ = 0.0
          GG(2)=0.0
          DO ITT=IT1,IT2
           G=GAUSS(ISEED)
           GG(ITT-IT1+1) = G * G * S(ITT-IT1+1)
           ZZ = ZZ + GG(ITT-IT1+1)
          ENDDO
          IF (GG(1).GT.0.0) THEN
           DMIX2 = GG(2) / GG(1)
          ELSE
           DMIX2 = 0.0
          ENDIF
          alpha = ALPH_TOT(EG,spfi,ipfi,0.0,
     *                    spin,ipin,dmix2,nent,elent,convt)
          Z = Z + ZZ * (1 + alpha)
         ENDDO  !IL

        ELSE                    ! Primary transitions (the same fluctuation)
         DO IL=1,NL
          ZZ = 0.0
          GG(2)=0.0
          DO ITT=IT1,IT2
           IF (IVER.EQ.2) THEN
            G1=GAUSS(ISEED)
            G2=GAUSS(ISEED)+CORRI(MODE)*G1
            GSQ=(Re2Res(mode)*G2*G2/(1+CORRI(mode)**2)+
     *           Im2Res(mode)*G1*G1)/(Re2Res(mode)+Im2Res(mode))
            GG(ITT-IT1+1) = GSQ * S(ITT-IT1+1)
           ELSE
            G=GAUSS(ISEED)
            GG(ITT-IT1+1) = G * G * S(ITT-IT1+1)
           ENDIF
           ZZ = ZZ + GG(ITT-IT1+1)
          ENDDO
          IF (GG(1).GT.0.0) THEN
           DMIX2 = GG(2) / GG(1)
          ELSE
           DMIX2 = 0.0
          ENDIF
          alpha = ALPH_TOT(EG,spfi,ipfi,0.0,
     *                     spin,ipin,dmix2,nent,elent,convt)
          Z = Z + ZZ * (1 + alpha)
         ENDDO  !IL
        ENDIF ! MODE
       ENDIF ! NOPTFL
c
       STCON(MODE,I,ISBS,IPFI)=STCON(MODE,I-1,ISBS,IPFI)+
     *                         DBLE(Z*SPAC)
       IF (IT.EQ.1) THEN
         GE1=GE1+Z
       ELSEIF (IT.EQ.3) THEN
         GM1=GM1+Z
       ELSEIF (IT.EQ.4) THEN
         GE2=GE2+Z
       ELSEIF (IT.EQ.2) THEN          ! This is not fully correct (strength fluctuate)
         GM1=GM1+Z*S(1)/(S(1)+S(2))   ! But it is acceptable approximation
         GE2=GE2+Z*S(2)/(S(1)+S(2))
       ENDIF
      ENDDO !I
C
C     For primary transitions it is assumed that MODE=0. In such
C     a case it is understood that the values of the ACTUAL
C     subscripted variable that replaces the FICTIVE variable STDIS
C     are simply derived from input data without the need of Monte
C     Carlo simulation.
C
      IF (IPRIM.EQ.1) THEN
        IF (MODE.NE.0) GOTO 22    !!!!!!!! Compare to WIDTHS_R()
      ENDIF

   24 CONTINUE
      DO I=1,NDIS(ISUBSC(SPFI),IPFI)
       EG=EIN-ENDIS(I,ISUBSC(SPFI),IPFI)
       ISDIS(MODE,I,ISUBSC(SPFI),IPFI)=ISEED   ! Needed only for mixing ratio
       IFLAG=0
       Z=0.
       DO ITT=IT1,IT2
        S(ITT-IT1+1)=SGAMMA(EG,EIN,ITT)
       ENDDO
       IF (NOPTFL.EQ.1) THEN
        IF (MODE.EQ.0) THEN
          ZZ = 0.0
          GG(2)=0.0
          DO ITT=IT1,IT2
           G=GAUSS(ISEED)
           GG(ITT-IT1+1) = G * G * S(ITT-IT1+1)
           ZZ = ZZ + GG(ITT-IT1+1)
          ENDDO
          IF (GG(1).GT.0.0) THEN
           DMIX2 = GG(2) / GG(1)
          ELSE
           DMIX2 = 0.0
          ENDIF
          alpha = ALPH_TOT(EG,spfi,ipfi,0.0,
     *                    spin,ipin,dmix2,nent,elent,convt)
          Z = Z + ZZ * (1 + alpha)
        ELSE                    ! Primary transitions (the same fluctuation)
         ZZ = 0.0
         GG(2)=0.0
         DO ITT=IT1,IT2
           IF (IVER.EQ.2) THEN
            G1=GAUSS(ISEED)
            G2=GAUSS(ISEED)+CORRI(MODE)*G1
            GSQ=(Re2Res(mode)*G2*G2/(1+CORRI(mode)**2)+
     *           Im2Res(mode)*G1*G1)/(Re2Res(mode)+Im2Res(mode))
            GG(ITT-IT1+1) = GSQ * S(ITT-IT1+1)
           ELSE
            G=GAUSS(ISEED)
            GG(ITT-IT1+1) = G * G * S(ITT-IT1+1)
           ENDIF
           ZZ = ZZ + GG(ITT-IT1+1)
         ENDDO
         IF (GG(1).GT.0.0) THEN
           DMIX2 = GG(2) / GG(1)
         ELSE
           DMIX2 = 0.0
         ENDIF
         alpha = ALPH_TOT(EG,spfi,ipfi,0.0,
     *                     spin,ipin,dmix2,nent,elent,convt)
         Z = Z + ZZ * (1 + alpha)
        ENDIF ! MODE
       ELSE
         ZZ = 0.0
         GG(2)=0.0
         DO ITT=IT1,IT2
          GG(ITT-IT1+1) = S(ITT-IT1+1)
          ZZ = ZZ + GG(ITT-IT1+1)
         ENDDO
         IF (GG(1).GT.0.0) THEN
          DMIX2 = GG(2) / GG(1)
         ELSE
          DMIX2 = 0.0
         ENDIF
         alpha = ALPH_TOT(EG,spfi,ipfi,0.0,
     *                    spin,ipin,dmix2,nent,elent,convt)
         Z = Z + ZZ * (1 + alpha)
       ENDIF ! NOPTFL
c
       IF (IT.EQ.1) THEN
         GE1=GE1+Z
       ELSEIF (IT.EQ.3) THEN
         GM1=GM1+Z
       ELSEIF (IT.EQ.4) THEN
         GE2=GE2+Z
       ELSEIF (IT.EQ.2) THEN           ! This is not fully correct (strength fluctuate)
         GM1=GM1+Z*S(1)/(S(1)+S(2))    ! But it is acceptable approximation
         GE2=GE2+Z*S(2)/(S(1)+S(2))
       ENDIF
c
       STDIS(MODE,I,ISBS,IPFI)=STDIS(MODE,I-1,ISBS,IPFI)+Z*SPAC     
      ENDDO

   22 TOTCON(MODE)=TOTCON(MODE)+STCON(MODE,NBIN,ISBS,IPFI)
      GACON(MODE,ISBS,IPFI)=TOTCON(MODE)
      TOTDIS(MODE)=TOTDIS(MODE)+
     &             STDIS(MODE,NDIS(ISUBSC(SPFI),IPFI),ISBS,IPFI)
      GADIS(MODE,ISBS,IPFI)=TOTDIS(MODE)
C
    1 CONTINUE
      RETURN
      END
C
C
C***********************************************************************
      SUBROUTINE ONESTEP
     *   (MODE,IPIN,SPIN,IBIN,ILIN,TOTCON,STCON,GACON,ISCON,TOTDIS,
     *    STDIS,GADIS,ISDIS,IPFI,SPFI,IBFI,ILFI,DMIX2,SIGN,IR,IRX)
C***********************************************************************
C
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      logical    ponverze,ponverk
      integer    dekod,denum,delev,deparity,deaux
      REAL       Im2Res
      REAL*4     EIN,EFI
      REAL*8     AUX,AUX0,TOTCON,DRN,STCON,GACON,GAC
      INTEGER*4  ISEED,IR,IRX,ISDIS,ircon,irconc,ntotal,ISEEDEN,SEEDS
      DIMENSION  TOTDIS(0:2),TOTCON(0:2),STCON(0:2,0:MAXBIN,-2:2,0:1),
     *           ISCON(0:2,0:MAXBIN,-2:2,0:1),GACON(0:2,-2:2,0:1),
     *           STDIS(0:2,0:30,-2:2,0:1),GADIS(0:2,-2:2,0:1),
     *           ISDIS(0:2,20,0:MAXJC,0:1),S(2),GG(2)
C
      COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
     &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
     &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
     &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
     &               denum(99),delev(99,20),despin(99,20),
     &               deparity(99,20),IRCONc(2),deltx(99,20)
     &       /GAU/   u,IFLAG
     &       /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,ntotal,iregi,
     *               noptcs,nlinc,capfr(2)
     &       /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
     &               CONVk(0:1,5,50)
     *       /kon/   ponverze,ponverk
     &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
     &               ISWIC,ISWMX,ISWWI,ISWLS
     &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)
C
      IF (iregi.eq.2) GOTO 23
      SPAC=1./DENSITY(EIN,SPIN,IPIN)
      IF (SPAC.LE.0.0) SPAC=0.00001
C
    5 DRN=(DBLE(INT(DBLE(RAN0(IR))/1.D-4))+DBLE(RAN0(IR)))*1.D-4*
     *     (TOTCON(MODE)+DBLE(TOTDIS(MODE)))
      IF (DRN-TOTCON(MODE)) 1,2,2
C
C     Label #1 means that the transition ends in the continuum; this
C     can be learnt from IREGI=0.
C
    1 IREGI=0
      SP  = SPIN - INT(SPIN + .25)
      ISP = INT(SPIN + .25)
      DO IPF=0,1
       DO ISPF = ISP-2, ISP+2
        SPF = SP + FLOAT(ISPF)
        ISBS=NINT(SPF+.25)-NINT(SPIN+.25)
        IT=ITYPE(SPIN,IPIN,SPF,IPF)
        IF (IT.NE.0) THEN
          GAC=GACON(MODE,ISBS,IPF)
          IF(DRN.LT.GAC) GOTO 4
        ENDIF
       ENDDO !SPF
      ENDDO !IPF
C
C     GOTO 5 statements are here for security reasons
C
      GOTO 5
C
C     Now IPFI,SPFI (& ISBS) are already determined:
C
    4 IPFI=IPF
      SPFI=SPF
      DRN=DRN-GAC+STCON(MODE,NBIN,ISBS,IPFI)
      IF (IT.EQ.2) THEN
        IT1=3
        IT2=4
      ELSE
        IT1=IT
        IT2=IT
      ENDIF
C
C     A DO-loop that follows should be later replaced by a faster
C     algorithm
C
      DO I=IBIN+1,NBIN
        IF (STCON(MODE,I,ISBS,IPFI).GT.DRN) GOTO 7
      ENDDO
      GOTO 5
    7 IBFI=I
C
C     Now even IBFI is known. It remains to determine ILFI.
C
      ISEED=ISCON(MODE,IBFI,ISBS,IPFI)
      IFLAG=0
C
C     The appropriate ISEED is now fetched and IFLAG reset.
C     The algorithm for finding ILFI can start.
C
      EG=EIN-BN+(FLOAT(IBFI)-0.5)*DELTA
c      EG=(FLOAT(IBFI-IBIN))*DELTA	       !Energies between mid-bins
      NL=LEVCON(IBFI,ISUBSC(SPFI),IPFI)-
     *   LEVCON(IBFI-1,ISUBSC(SPFI),IPFI)
      AUX0=STCON(MODE,IBFI-1,ISBS,IPFI)
      Z=0.
      DO ITT=IT1,IT2
        S(ITT-IT1+1)=SGAMMA(EG,EIN,ITT)
      ENDDO





      DO IL=1,NL
       ZZ = 0.0
       GG(2)=0.0
       IF (NOPTFL.EQ.0) THEN   ! no fluctuations
         DO ITT=IT1,IT2
          GG(ITT-IT1+1) = S(ITT-IT1+1)
          ZZ = ZZ + GG(ITT-IT1+1)
         ENDDO
         IF (GG(1).GT.0.0) THEN
          DMIX2 = GG(2) / GG(1)
         ELSE
          DMIX2 = 0.0
         ENDIF
         alpha = ALPH_TOT(EG,spfi,ipfi,0.0,
     *                    spin,ipin,dmix2,nent,elent,convt)
         Z = Z + ZZ * (1 + alpha)
       ELSE
        IF (MODE.EQ.0) THEN  ! secondary transitions
         DO ITT=IT1,IT2
          G=GAUSS(ISEED)
          GG(ITT-IT1+1) = G * G * S(ITT-IT1+1)
          ZZ = ZZ + GG(ITT-IT1+1)
         ENDDO
         IF (GG(1).GT.0.0) THEN
          DMIX2 = GG(2) / GG(1)
         ELSE
          DMIX2 = 0.0
         ENDIF
         alpha = ALPH_TOT(EG,spfi,ipfi,0.0,
     *                    spin,ipin,dmix2,nent,elent,convt)
         Z = Z + ZZ * (1 + alpha)
        ELSE
         IF (IVER.EQ.2) THEN  ! primary transitions - not exactly PT distribution
          DO ITT=IT1,IT2
           G1=GAUSS(ISEED)
           G2=GAUSS(ISEED)+CORRI(MODE)*G1
           GSQ=(Re2Res(mode)*G2*G2/(1+CORRI(mode)**2)+
     *          Im2Res(mode)*G1*G1)/(Re2Res(mode)+Im2Res(mode))
           GG(ITT-IT1+1) = GSQ * S(ITT-IT1+1)
           ZZ = ZZ + GG(ITT-IT1+1)
          ENDDO
          IF (GG(1).GT.0.0) THEN
           DMIX2 = GG(2) / GG(1)
          ELSE
           DMIX2 = 0.0
          ENDIF
          alpha = ALPH_TOT(EG,spfi,ipfi,0.0,
     *                     spin,ipin,dmix2,nent,elent,convt)
          Z = Z + ZZ * (1 + alpha)
         ELSE  ! primary transitions - PT distribution
          DO ITT=IT1,IT2
            G=GAUSS(ISEED)
            GG(ITT-IT1+1) = G * G * S(ITT-IT1+1)
            ZZ = ZZ + GG(ITT-IT1+1)
          ENDDO
          IF (GG(1).GT.0.0) THEN
           DMIX2 = GG(2) / GG(1)
          ELSE
           DMIX2 = 0.0
          ENDIF
          alpha = ALPH_TOT(EG,spfi,ipfi,0.0,
     *                     spin,ipin,dmix2,nent,elent,convt)
          Z = Z + ZZ * (1 + alpha)
         ENDIF  ! IVER
        ENDIF  ! MODE
       ENDIF ! NOPTFL
       AUX=AUX0+DBLE(Z*SPAC)
       IF (AUX.GT.DRN) GOTO 9
      ENDDO !IL
      GO TO 5
C
C     The resulting ILFI:
C
    9 ILFI=IL
C
      NLEV1=LEVCON(IBFI,ISUBSC(SPFI),IPFI)-
     *      LEVCON(IBFI-1,ISUBSC(SPFI),IPFI)
      Q1=1.-FLOAT(2*ILFI-1)/FLOAT(2*NLEV1)
c      ISEEDEN=SEEDS(0,SPFI,IPFI,ILFI,IBFI)
c      Q1=RAN0(ISEEDEN)
      EFI=BN-(FLOAT(IBFI)-Q1)*DELTA  ! Energy of the final level random in bin
      IF (efi.le.eall) IREGI=2
c
c       delta**2 and conversion - TODO - sign must be changed 
c
      IF (NOPTFL.NE.1) THEN
        GG(1)=1.
c        GG(2)=1.
      ENDIF
      IF (IT.EQ.2) THEN
        GG(1) = GAUSS(ISEED)
        IF (GG(1).NE.0.) THEN
c          DMIX2=GG(2)**2*S(2)/GG(1)**2/S(1)
          IF (GG(1).GT.0) THEN
            SIGN=1.
          ELSE
            SIGN=-1.
          ENDIF
c        ELSE
c          DMIX2=1.E+6
c          SIGN=1.
        ENDIF
      ELSE
        DMIX2=0.
        SIGN=1.
      ENDIF
      dummy=RAN0(IRX)
      ALPHA=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,
     *               NENT,ELENT,CONVT)
      IF (dummy.GT.(ALPHA/(1+ALPHA))) THEN
        ponverze=.FALSE.
        ponverk=.FALSE.
      ELSE
        ALPHAK=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,
     *                  NENK,ELENK,CONVK)
        ponverze=.TRUE.
        IF (dummy.LE.(ALPHAK/(1+ALPHAK))) THEN
          ponverk=.TRUE.
        ELSE
          ponverk=.FALSE.
        ENDIF
      ENDIF
      RETURN
C
    2 IREGI=1
      RN=SNGL(DRN-TOTCON(MODE))
      SP  = SPIN - INT(SPIN + .25)
      ISP = INT(SPIN + .25)
      DO IPF=0,1
       DO ISPF = ISP-2, ISP+2
        SPF = SP + FLOAT(ISPF)
        ISBS=NINT(SPF+.25)-NINT(SPIN+.25)
        IT=ITYPE(SPIN,IPIN,SPF,IPF)
        IF (IT.NE.0) THEN
         IF (RN.LT.GADIS(MODE,ISBS,IPF)) GOTO 11
        ENDIF
       ENDDO !SPF
      ENDDO !IPF
      GOTO 5
   11 SPFI=SPF
      IPFI=IPF
C
      RN=RN-GADIS(MODE,ISBS,IPFI)+
     *   STDIS(MODE,NDIS(ISUBSC(SPFI),IPFI),ISBS,IPFI)
      DO I=1,NDIS(ISUBSC(SPFI),IPFI)
        IF (RN.LT.STDIS(MODE,I,ISBS,IPFI)) GOTO 15
      ENDDO
      GOTO 5
   15 ILFI=I
      ibin=0                                        !!!!Why?
      EFI=ENDIS(ILFI,ISUBSC(SPFI),IPFI)
      IF (EFI.LE.EALL) IREGI=2
c
c     delta**2 and conversion
c
      ISEED=ISDIS(MODE,ILFI,ISUBSC(SPFI),IPFI)
      IFLAG=0
      IF (IT.EQ.2) THEN
         IT1=3
         IT2=4
         ELSE
         IT1=IT
         IT2=IT
      ENDIF
      EG=EIN-EFI
      DO ITT=IT1,IT2
        S(ITT-IT1+1)=SGAMMA(EG,EIN,ITT)
      ENDDO
      IF (NOPTFL.EQ.1) THEN
       IF (MODE.EQ.0) THEN
        DO ITT=IT1,IT2
         G=GAUSS(ISEED)
         GG(ITT-IT1+1)=G
        ENDDO
       ELSE                             !Primary transitions
        DO ITT=IT1,IT2
c         G1=GAUSS(ISEED)
c         G2=GAUSS(ISEED)+CORRI(MODE)*G1
c         GSQ=(Re2Res(mode)*G2*G2/(1+CORRI(mode)**2)+
c     *        Im2Res(mode)*G1*G1)/(Re2Res(mode)+Im2Res(mode))
c         GG(ITT-IT1+1)=sqrt(GSQ)
c         IF (G1.LT.0) GG(ITT-IT1+1)=-1.*GG(ITT-IT1+1)
         G=GAUSS(ISEED)
         GG(ITT-IT1+1)=G
        ENDDO
       ENDIF
      ENDIF
      IF (NOPTFL.NE.1) THEN
        GG(1)=1.
        GG(2)=1.
      ENDIF
      IF (IT.EQ.2) THEN
        GG(1) = GAUSS(ISEED)
        IF (GG(1).NE.0.) THEN
c          DMIX2=GG(2)**2*S(2)/GG(1)**2/S(1)
          IF (GG(1).GT.0) THEN
            SIGN=1.
          ELSE
            SIGN=-1.
          ENDIF
        ELSE
          DMIX2=1.E+6
          SIGN=1.
        ENDIF
      ELSE
        DMIX2=0.
        SIGN=1.
      ENDIF
      dummy=RAN0(IRX)
      ALPHA=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,
     *               NENT,ELENT,CONVT)
      IF (dummy.GT.(ALPHA/(1+ALPHA))) THEN
        ponverze=.FALSE.
        ponverk=.FALSE.
      ELSE
        ALPHAK=ALPH_TOT(EIN,SPIN,IPIN,EFI,SPFI,IPFI,DMIX2,
     *                  NENK,ELENK,CONVK)
        ponverze=.TRUE.
        IF (dummy.LE.(ALPHAK/(1+ALPHAK))) THEN
          ponverk=.TRUE.
        ELSE
          ponverk=.FALSE.
        ENDIF
      ENDIF
      RETURN
C
C
C     Note that this time ILFI means order number of an actual
C     level, not the order number of a level in an energy bin of
C     the energy continuum.
C
C     End in a part, where one knows every branching ratios
C
  23  EIN=EFI
      deaux=dekod(ilfi,isubsc(spfi),ipfi)
      rn=ran0(ir)*sall(deaux,denum(deaux))
      DO ia=1,denum(deaux)
        IF (rn.lt.sall(deaux,ia)) GOTO 24
      ENDDO
      GOTO 23
  24  ilfi=delev(deaux,ia)
      spfi=despin(deaux,ia)
      ipfi=deparity(deaux,ia)
      ibin=0
      efi=endis(ilfi,ISUBSC(spfi),ipfi)
      dmix2=deltx(deaux,ia)**2
      IF (deltx(deaux,ia).LT.0) THEN
        sign=-1.
      ELSE
        sign=1.
      ENDIF
c
C     conversion
c
      if (ponv(deaux,ia).gt.0.0) then
         dummy=ran0(ir)
         if (dummy.gt.ponv(deaux,ia)) then
           ponverze=.false.
           ponverk=.false.
         else
           ponverze=.true.
           if (dummy.le.ponvk(deaux,ia)) then
             ponverk=.true.
           else
             ponverk=.false.
           endif
         endif
      endif
      RETURN
      END
C

c***********************************************************************
c                                             version 2.54
c
c   - changed meaning of some of NOPTE1 (3,12-17) in this version
c
c Functions TERM(Excit_En)
c           SGAM(Egam,Eini,Ityp)
c
C***********************************************************************
      FUNCTION SGAMMA(EGAM,EINI,ITYP)
c
c   Photon strengths for E1, M1+E2, M1 and E2 transitions. "Photon
c   strength" does not mean "photon strength function" here, but
c   photon strength function multiplied by EGAM**(2*L+1) and, in
c   the case of M1+E2 transitions, summed over both XL-components.
c
c     ITYP is equal to 1, 2, 3 or 4 (see ITYPE function)
c     EINI is the initial state energy in MeV
c     EGAM is gamma-ray energy in MeV
c
c   E1:  NOPTE1= 0: Single-particle approximation
c                1: Classical Lorentzian GDER
c                2: GDER with an energy and temperature dependent
c                   damping width (J.Kopecky, R.Chrien, Nucl.Phys.
c                   A468,p.285)
c                3: correct EGLO model - see 6
c                4: Kadmenskij-Markushev-Furman original Strength function
c                   (no high energy approximation according to Chrien)
c                5: The Chrien's Strength function (Nucl. Phys. A468, 285
c                   (1987)) only. In 3: is this model used only for
c                   high energy region
c                6: The strength function according to Chrien with
c                   phenomenological temperature dependent damping
c                   proposed by Kopecky (Distribution of Radiative Strength
c                   in Gd-156, 157 and 158 Nuclei)
c                   - I have found an error, correct EGLO is 3:
c                7: GDER with phenomenological temperature dependent
c                   damping width proposed by Kopecky
c                8: 4: with the first resonance of Lorentz type
c                9: 6: with the first resonance of Lorentz type
c               10: KMF (4:) for EG<4 MeV; lin. combination of KMF and BA
c                   for 4 MeV<EG<8 MeV; BA (1:) for EG>8 MeV
c               12: MLO1 (agrees with RIPL2 code)
c               16: MLO2 (agrees with RIPL2 code)
c               17: MLO3 (agrees with RIPL2 code)
c               31-40: correspond to 1-10 for high EGAM; for low EGAM original
c                      values are multiplied by a factor (given in input data)
c                      in between the PSF is a linear combination ... 
c                      motivation comes from Au
c               51: KMF (4:) without temperature-dependent term in damping
c                   width
c               52: EGLO (6:) without temperature-dependent term in damping
c                   width; temperature is taken into account only in the
c                   "second term" - FK*...
c               53: EELO (7:) without temperature-dependent term in damping
c                   width - i.e. no termperature dependence assumed
c               41: KMF according to Oslo group 
c
c   M1:  NOPTM1= 0: Single-particle approximation
c                1: Classical lorentzian GDMR
c                3: Scissors (first) resonance is build up only on states
c                   with excitation energy lower than PAR_M1(1)
c                4: Classical lorentzian build on the "background"
c                   that is described by the SP (constant function)
c               12: Scissors (first) resonance is build up only on states
c                   with excitation energy lower than PAR_M1(1) + SP on all states (originally NOPTM1=11)
c
c                ?: Enery of scissors resonance depends linearly on
c                   the energy of final state (and is build up only on states
c                   below certain excitation energy)
c                ?: Scissors resonance is considered only for primary transitions
c
c                ?: power dependence
c
c   E2:  NOPTE2= 0: Single-particle approximation
c                1: Classical Lorentzian GQER
c
C***********************************************************************
c
      PARAMETER  (PIH=  8.673592583E-08, ! 1/(3*(pi*hbar*c)**2)
     &            PIHQ= 5.204155555E-08, ! 1/(5*(pi*hbar*c)**2)
     &            PI42=39.4784176)       ! 4*pi**2
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
c
      COMMON /PHOTO/  NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
     &                SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
     &                DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
     &                PAR_E1(3),DIPSLP,DIPZER
     &       /DEN/    NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *                DENLO,DENHI,DENPAR(4),
     *                ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
     &       /MAINE/  eall,EIN,efi,NOPTFL,BN,delta,nbin,ntotal,iregi,
     &                noptcs,nlinc,capfr(2)
     &       /SFCE/   SFCEE1,SFCEM1,SFCEE2,GE1,GM1,GE2,
     &                GE1DIS(2),GM1DIS(2),GE2DIS(2)
     &       /incr/  nfilev,STEPS,enrgfin(18),buff(251),
     &               ncum,cumwidth(14),edo(18),eup(18)
c
      SGAMMA=0.
      SFCEE1=0
      SFCEM1=0
      SFCEE2=0
      IF ((ITYP.GT.4).OR.(ITYP.LT.1).OR.(EGAM.LE.0.)) RETURN
c
c
c*****                       E1 component
c
      IF (ITYP.EQ.1) THEN
c
        IF     (NOPTE1.EQ.0) THEN   ! The single-particle approximation
          SGAMMA=DEG*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.1) THEN   ! Classical Lorentzian
          Q=0.
          DO I=1,NGIGE              ! loop over both GDR peaks
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.2) THEN   ! GDER with E,T-dependent damping
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
     &!     ^ .......... energy and temperature dependent width
            QQ=SIG(I)*W0(I)*
     &        (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.3) THEN  !Empirical generalization of temperature
c          dependent damping according to Kopecky in Chrien model (EGLO)
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            WPHENZ=EK0-(1.-EK0)*EGZERO/(ER(I)-EGZERO)
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5*WPHENZ
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.4) THEN   ! Pure Fermi liquid theory (Kadmenskij)
          TFIN=TERM(EINI-EGAM)          ! (no high energy approximation)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.5) THEN    !Pure Chrien model (similar to !OPT=3,
c          but no Kadmenskij for low EGAM) viz Kopecky
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.6) THEN  !Empirical generalization of temperature
c          dependent damping according to Kopecky in Chrien model (EGLO)
c          - an error found in the SLIM - corrected in NOPTE1=3
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.7) THEN   ! SMLO
          IF (EGAM.LE.EINI) THEN
            TFIN=SQRT((EINI-EGAM)/AMASS*10.0)
          ELSE
            TFIN=0.0
          ENDIF
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM*ER(I)+PI42*TFIN**2)/ER(I)**2
            QQ=SIG(I)*W0(I)*W*EGAM / (1.0-EXP(-EGAM/TFIN)) /
     &         ((EGAM**2-ER(I)**2)**2+(EGAM*W)**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.207) THEN  !Empirical generelization of temperature
c        dependent damping according to Kopecky aplied to TD model (EELO)
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
     &!     ^ .......... energy and temperature dependent width
            QQ=SIG(I)*W0(I)*
     &        (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.8) THEN   ! Fermi liquid theory (Kadmenskij)
          TFIN=TERM(EINI-EGAM)      ! with 1st resonance of Lorentz. shape
          Q=0.
          QQ=SIG(1)*
     &       (EGAM*W0(1)**2/((EGAM**2-ER(1)**2)**2+(EGAM*W0(1))**2))
          Q=Q+QQ
          DO I=2,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.9) THEN  !Empirical generalization of temperature
c          dependent damping according to Kopecky in Chrien model (EGLO)
c          with the first resonance of Lorentzian type
          TFIN=TERM(EINI-EGAM)
          Q=0.
          QQ=SIG(1)*
     &       (EGAM*W0(1)**2/((EGAM**2-ER(1)**2)**2+(EGAM*W0(1))**2))
          Q=Q+QQ
          DO I=2,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.10) THEN   ! KMF for low energies
          TFIN=TERM(EINI-EGAM)       ! Mix KMF and BA for higher energies
          Q=0.
          IF (EGAM.LE.PAR_E1(1)) THEN
           DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
           ENDDO
          ELSE
           x=(EGAM-PAR_E1(1))/(PAR_E1(2)-PAR_E1(1))         ! Admixture of BA to KMF
           IF (x.GT.1.) x=1.
           DO I=1,NGIGE
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+x*QQ
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+(1.-x)*QQ
           ENDDO
          ENDIF
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.11) THEN                            ! Goriely
          SGAMMA = APSF(EGAM,1)*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.12) THEN                            ! Goriely
          SGAMMA = APSF(EGAM,1)
          SGAMMA = SGAMMA + EINI * PAR_E1(1) / (1.0+EXP(EGAM-PAR_E1(2))) 
          SGAMMA = SGAMMA * EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.25) THEN    !Pure Chrien+Kopecky model but different T (similar to OPT=3)
          TFIN=TERM(BN-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.31) THEN  ! Classical Lor.; suppressed for small EGAM
          Q=0.
          DO I=1,NGIGE              ! loop over both GDR peaks
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.34) THEN   ! Pure Fermi liquid theory (Kadmenskij)
c                                      suppressed for small EGAM
          TFIN=TERM(EINI-EGAM)          
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.36) THEN  !EGLO (6 - incorrect); suppressed for small EGAM
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.37) THEN  !EELO; suppressed for small EGAM
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
     &!     ^ .......... energy and temperature dependent width
            QQ=SIG(I)*W0(I)*
     &        (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          x=DIPSLP*EGAM+DIPZER 
          if (x.LT.PAR_E1(3)) x=PAR_E1(3)
          if (x.GT.1.0)    x=1.0
          SGAMMA=PIH*Q*EGAM**3*x
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.74) THEN   ! KMF with constant T
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.76) THEN  ! EGLO(6) with constant T
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.80) THEN   ! KMF for low energies
          TFIN=PAR_E1(3)                ! Mix KMF and BA for higher energies
          Q=0.
          IF (EGAM.LE.PAR_E1(1)) THEN
           DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
           ENDDO
          ELSE
           x=(EGAM-PAR_E1(1))/(PAR_E1(2)-PAR_E1(1))         ! Admixture of BA to KMF
           IF (x.GT.1.) x=1.
           DO I=1,NGIGE
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+x*QQ
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+(1.-x)*QQ
           ENDDO
          ENDIF
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.83) THEN  ! EGLO(6) with constant T + Lorentzian pygmy (Oslo - actinides); PAR_E1(1) !!!
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=2,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          DO I=1,1              ! pygmy
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.82) THEN  ! EGLO(6) with constant T + Lorentzian pygmy (Oslo - actinides); PAR_E1(1) !!!
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=4,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          DO I=1,3              ! pygmy
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.84) THEN   ! GLO + gaussian PR (Oslo - Sn)
          Q=0.
          TFIN = PAR_E1(1)
          DO I = 2, NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
c          !     ^ .......... energy and temperature dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          DO I = 1, 1
            QQ = SIG(I) / W0(I) / SQRT(2.0 * 3.141592) *
     &           EXP(-(EGAM-ER(I))**2 / 2.0 / W0(I)**2)
            Q=Q+QQ
          ENDDO
          SGAMMA= PAR_E1(3) * PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN

c
        ELSEIF (NOPTE1.EQ.64) THEN   ! SMLO (RIPL-3)
c
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W = EINI * W0(I) / ER(I)
            IF (TFIN.GT.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.67) THEN   ! Mughabghab+Dunford - formulation from RIPL2 
c
          S2=217.16/AMASS/AMASS      !global parametrization, see RIPL2 
          CDQ=1.05
          TFIN=TERM(EINI-EGAM)
          BETA2 = PAR_E1(3)*PAR_E1(3)    !square of deformation
c
          Q=0.
          DO I=1,NGIGE
            WDQ0= CDQ*SQRT(BETA2*ER(I)*ER(I)+ER(I)*S2)
            WDQ = CDQ*SQRT(BETA2*EGAM*EGAM+EGAM*S2)
            WC  = (W0(I)-WDQ0)/ER(I)/ER(I)*(EGAM**2+PI42*TFIN**2)
            W   = WDQ + WC 
            QQ=SIG(I)*W0(I)*FERMC*ER(I)*W/
     &         ((EGAM**2-ER(I)**2)**2+FERMC*(EGAM*W)**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.68) THEN   ! Mughabghab+Dunford - original formulation (PLB487,155)  
c
          S2=217.16/AMASS/AMASS      !global parametrization, see RIPL2 
          CDQ=1.05
          TFIN=TERM(EINI-EGAM)
          BETA2 = PAR_E1(3)*PAR_E1(3)    !square of deformation
c
          Q=0.
          DO I=1,NGIGE
            WDQ0= CDQ*SQRT(BETA2*ER(I)*ER(I)+ER(I)*S2)
            WDQ = CDQ*SQRT(BETA2*EGAM*EGAM+EGAM*S2)
            WM  = (W0(I)-WDQ0)/ER(I)/ER(I)*(EGAM**2+PI42*TFIN**2)
            W   = WDQ + WM 
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.69) THEN   ! Goriely (PLB436,10) - in RIPL2
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*FERMC*(EGAM**2+PI42*TFIN**2)/ER(I)/EGAM
     &!     ^ .......... energy and temperature dependent width
            QQ=SIG(I)*W0(I)*
     &        (EGAM*W/((EGAM**2-ER(I)**2)**2+EGAM**2*W0(I)*W))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
c
        ELSEIF (NOPTE1.EQ.112) THEN   ! MLO1 (Original Plujko)
c         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
c         !ktere mohou podle 'Plujkovy teorie' nastat 
c         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
c         !hodnoty nekterych 'promennych' nastaveny zde
c
          TFIN=TERM(EINI-EGAM)
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          ER0PL=SMFREQ()
c
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)/ALPPL
            WDR=ER(I)**2/ALPPL
            WR=2*WDR*(ER(I)**2+ER0PL**2)/
     &         ((ER(I)**2-ER0PL**2)**2+4*(WDR*ER(I))**2)
            W=2*WD*W0(i)/WR*(ER(I)**2+ER0PL**2)/    
     &        ((ER(I)**2-ER0PL**2)**2+4*(WD*EGAM)**2)

            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.13) THEN   ! MLO1 (Original Plujko); Jina T!!!!
c         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
c         !ktere mohou podle 'Plujkovy teorie' nastat 
c         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
c         !hodnoty nekterych 'promennych' nastaveny zde
c
          TFIN=TERMDILG(EINI-EGAM)
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          ER0PL=SMFREQ()
c
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)/ALPPL
            WDR=ER(I)**2/ALPPL
            WR=2*WDR*(ER(I)**2+ER0PL**2)/
     &         ((ER(I)**2-ER0PL**2)**2+4*(WDR*ER(I))**2)
            W=2*WD*W0(i)/WR*(ER(I)**2+ER0PL**2)/    
     &        ((ER(I)**2-ER0PL**2)**2+4*(WD*EGAM)**2)

            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.14) THEN   ! MLO - RIPL2 
c         !it is not clear if there is a factor EINI or (EGAM+a*TFIN**2) 
c          or (EGAM+a*TINI**2) in the width WD 
c
          TFIN=TERM(EINI-EGAM)
          TINI=TERM(EINI)
c          TINI=TERM1(EINI)
          ALPPL=0.00542              !in RIPL2 error in the order (5.42e-4)
          FACTOR=1.0                 !All parameters recommended by Plujko
          ALPPL=ALPPL/FACTOR
          FKs0=0.3                          
          FNS=1.0                           
c
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)*ALPPL
c            WD=(EGAM+UINI)*ER(I)*ALPPL
            FKR=(W0(I)-ALPPL*ER(I)**2)/WWALL1()
            IF (EGAM.LE.(2*ER(I))) THEN
             FKS1=FKR+(FKs0-FKR)*ABS((EGAM-W0(I))/W0(I))**FNS             
            ELSE
             FKS1=FKs0
            ENDIF
            W=WD+FKS1*WWALL1()
            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.16) THEN   ! MLO2 (Original Plujko)
c         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
c         !ktere mohou podle 'Plujkovy teorie' nastat 
c         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
c         !hodnoty nekterych 'promennych' nastaveny zde
c
          TFIN=TERM(EINI-EGAM)
          BC=1.0                            ! Recommended by Plujko
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          EFERMI=37.                        ! Recommended by Plujko 
          WWALL=6.857764*sqrt(EFERMI)/R0()  ! origin of the const. is a puzzle for me
          FKs0=0.3                          ! Recommended by Plujko
          FNS=1.0                           ! Recommended by Plujko
c
          Q=0.
          DO I=1,NGIGE
            WD=EINI*ER(I)/ALPPL
            IF (NGIGE.EQ.1) THEN
              W=WD+FKS(EGAM,FNS,FKs0,WWALL,ALPPL)*WWALL
            ELSE
              W=WD+(W0(I)-ER(I)**2/ALPPL)*FKS(EGAM,FNS,FKs0,WWALL,ALPPL)
     &          /FKS(ER(I),FNS,FKs0,WWALL,ALPPL)
            ENDIF
            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.17) THEN   ! MLO3 (Original Plujko)
c         !tato procedura nepostihuje mozna uplne vsechny pripady (hodnoty parametru),
c         !ktere mohou podle 'Plujkovy teorie' nastat 
c         !urcite neni zakomponovana jeho moznost KEYSET=2 (je jine ALPPL) 
c         !hodnoty nekterych 'promennych' nastaveny zde
c
          TFIN=TERM(EINI-EGAM)
          BC=1.0                            ! Recommended by Plujko
          FACTOR=1.0                        ! Recommended by Plujko
          ALPPL=185.659                     ! 4*pi^2*alpha_free - from Plujko
          ALPPL=ALPPL/FACTOR
          EFERMI=37.                        ! Recommended by Plujko 
          WWALL=6.857764*sqrt(EFERMI)/R0()  ! origin of the const. is a puzzle for me
          FKs0=0.3                          ! Recommended by Plujko
          FNS=1.0                           ! Recommended by Plujko
c
          Q=0.
          DO I=1,NGIGE
            WD=(EGAM**2+PI42*TFIN**2)/ALPPL
            IF (NGIGE.EQ.1) THEN
              W=WD+FKS(EGAM,FNS,FKs0,WWALL,ALPPL)*WWALL
            ELSE
              W=WD+(W0(I)-ER(I)**2/ALPPL)*FKS(EGAM,FNS,FKs0,WWALL,ALPPL)
     &          /FKS(ER(I),FNS,FKs0,WWALL,ALPPL)
            ENDIF
            IF (TFIN.NE.0) THEN
             QQ=SIG(I)*W0(I)/(1-exp(-EGAM/TFIN))*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ELSE
             QQ=SIG(I)*W0(I)*
     &          (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.18) THEN   ! Sirotkin (1) - with KMF damping
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            O=FERMC*ER(I)
            IF (TFIN.NE.0) THEN
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/
     &         O**2/(1-exp(-EGAM/TFIN))
            ELSE
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/
     &         O**2
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.19) THEN   ! Sirotkin (2) - with damping from Sir
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
c            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            W=0.666*4./35.*(EGAM**2+PI42*TFIN**2)
            O=FERMC*ER(I)
            IF (TFIN.NE.0) THEN
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/
     &         O**2/(1-exp(-EGAM/TFIN))
            ELSE
            QQ=SIG(I)*W0(I)*W*EGAM*(EGAM**2+O**2)/(EGAM**2-ER(I)**2)**2/
     &         O**2
            ENDIF
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.20) THEN   ! Pure Fermi liquid theory (Kadmenskij)
          TFIN=TERM(EINI-EGAM)       ! with theoretical damping from Sir
          Q=0.
          DO I=1,NGIGE
c            W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
            W=0.666*4./35.*(EGAM**2+PI42*TFIN**2)
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
        ELSEIF (NOPTE1.EQ.51) THEN    ! Fermi liquid theory (Kadmenskij)
                                      ! without dependence on T !!!
          Q=0.
          DO I=1,NGIGE
            W=W0(I)*(EGAM**2)/ER(I)**2
            QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.52) THEN  !Modified empirical generalization
c          of temperature dependent damping according to Kopecky
c          in Chrien model (EGLO); Temperature not assumed in damping
c          width - see above
          TFIN=TERM(EINI-EGAM)
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2)/ER(I)**2
c          !     ^ .......... energy dependent width
            SLIM=FERMC*PI42*TFIN**2*W0(I)/ER(I)**5
c          !     ^ ...... the non-zero limit at Egam-->0
            QQ=SIG(I)*W0(I)*
     &          (SLIM+EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.53) THEN  !Modified empirical generelization
c        of temperature dependent damping according to Kopecky aplied
c        to TD model (EELO); temperature dependence removed - see above
          Q=0.
          DO I=1,NGIGE
            WPHEN=EK0+(1.-EK0)*(EGAM-EGZERO)/(ER(I)-EGZERO)
            W=WPHEN*W0(I)*(EGAM**2)/ER(I)**2
     &!     ^ .......... energy and temperature dependent width
            QQ=SIG(I)*W0(I)*
     &        (EGAM*W/((EGAM**2-ER(I)**2)**2+(EGAM*W)**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEE1=SGAMMA
          RETURN
c
c
c
        ELSEIF (NOPTE1.EQ.41) THEN   !!!!!! Oslo KMF for Yb
          TFIN=TERMOSLO(EINI-EGAM)        
c          TFIN=0.34        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.42) THEN   !!!!!! Oslo KMF - our temperature
          TFIN=TERM(EINI-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.43) THEN   !!!!!! Oslo KMF; including softpole
          TFIN=TERMOSLO(EINI-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
          RETURN
        ELSEIF (NOPTE1.EQ.44) THEN   !!!!!! Oslo KMF - our temperature; including softpole
          TFIN=TERM(EINI-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.45) THEN   !!!!!! Oslo KMF - const T 
          TFIN=TERMOSLO(BN-EGAM)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.46) THEN   !!!!!! Oslo KMF; including softpole
          TFIN=TERMOSLO(BN-EGAM)       !Temperature independent PSF  
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.47) THEN   !!!!!! Oslo KMF; including softpole
          TFIN=TERMOSLO(EINI-EGAM)     !!!! const below Etr   
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_M1(1)) EGAM=PAR_M1(1)
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          EGAM=EGAMOLD
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.48) THEN   !!!!!! Oslo KMF - T independent; including softpole
          TFIN=TERMOSLO(BN-EGAM)     !!!! const below Etr   
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_M1(1)) EGAM=PAR_M1(1)
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          EGAM=EGAMOLD
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN

        ELSEIF (NOPTE1.EQ.49) THEN   !!!!!! Oslo KMF - T const; including softpole
          TFIN=EK0        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN



        ELSEIF (NOPTE1.EQ.71) THEN   ! 34 (Kadmenskij)
c                                      with enhancement at high energies
          TFIN=TERM(EINI-EGAM)          
          Q=0.

          IF (EGAM.LT.7.35) THEN
           DO I=1,NGIGE
             W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
             QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
             Q=Q+QQ
           ENDDO
           x=DIPSLP*EGAM+DIPZER 
           if (x.LT.PAR_E1(3)) x=PAR_E1(3)
           if (x.GT.1.0)    x=1.0
           if ((EGAM.GT.6.0).AND.(EGAM.LT.7.35)) 
     &         x=x+x*2*(EGAM-6.0)
           SGAMMA=PIH*Q*EGAM**3*x
          ELSE
           DO I=1,NGIGE              ! BA
             QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
             Q=Q+QQ
           ENDDO
           SGAMMA=PIH*Q*EGAM**3
          ENDIF
          SFCEE1=SGAMMA
          RETURN

        ELSEIF (NOPTE1.EQ.78) THEN   !!!!!! Voinov Mo - no LLR
c          TFIN=0.015        
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=1,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
        ELSEIF (NOPTE1.EQ.79) THEN   !!!!!! Voinov Mo LLR
c          TFIN=0.015        
          TFIN=PAR_E1(3)        
          Q=0.
          DO I=1,1                   ! Lorentzian LLR
            QQ=SIG(I)*
     &        (EGAM*W0(I)**2/((EGAM**2-ER(I)**2)**2+(EGAM*W0(I))**2))
            Q=Q+QQ/PAR_E1(1)
          ENDDO
          DO I=2,NGIGE
           W=W0(I)*(EGAM**2+PI42*TFIN**2)/ER(I)**2
           QQ=FERMC*SIG(I)*W0(I)*W*ER(I)/(EGAM**2-ER(I)**2)**2
           Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3    
          SFCEE1=SGAMMA
        RETURN
c
        ELSEIF (NOPTE1.EQ.28) THEN  !Hokus - pokus
          SGAMMA=DEG*EGAM**8
          SFCEE1=SGAMMA
          RETURN
        ENDIF
      ENDIF
c
c
c*****                       M1 component
c
      IF ((ITYP.EQ.2).OR.( ITYP.EQ.3)) THEN
c
        IF (NOPTM1.EQ.0) THEN       ! The single-particle approximation
          SGAMMA=DMG*EGAM*EGAM*EGAM
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.1) THEN   ! Classical Lorentzian
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.3) THEN   ! Scissors mode is build up only on states
          Q=0.                      ! with Efinal lower than Etr
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF (Efinal.LE.PAR_M1(1)) THEN
            QQ=SIGM(1)*EGAM*WM0(1)**2/
     &        ((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
            Q=Q+QQ
          ENDIF
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.112) THEN   ! Scissors mode is build up only on states
          Q=0.                       ! with Efinal lower than Etr; SP everywhere ! Originally NOPTM1=11
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF (Efinal.LE.PAR_M1(1)) THEN
            QQ=SIGM(1)*EGAM*WM0(1)**2/
     &        ((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
            Q=Q+QQ
          ENDIF
          SGAMMA=(PIH*Q+DMG)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.4) THEN   ! Classical Lorentzian on backgroung
          Q=0.                      ! which is described by SP up to Etr
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          IF (EGAM.LE.PAR_M1(1)) THEN
            SGAMMA=(PIH*Q+DMG)*EGAM**3
          ELSE
            SGAMMA=PIH*Q*EGAM**3
          ENDIF
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.7) THEN   ! SMLO + "upbend"
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q
          SGAMMA = SGAMMA + 
     *             PAR_M1(1) * EXP(-PAR_M1(2)*EGAM) *
     *             (1.0+PAR_M1(3)*EGAM*EGAM*EGAM)
          SGAMMA = SGAMMA * EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN


        ELSEIF (NOPTM1.EQ.5) THEN  ! SP on states
          Q=0.                     ! with Efinal lower than Etr
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &         ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=EINI-EGAM
          IF (Efinal.LE.PAR_M1(1)) THEN
            SGAMMA=(PIH*Q+DMG)*EGAM**3
          ELSE
            SGAMMA=PIH*Q*EGAM**3
          ENDIF
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.6) THEN   ! Classical Lorentzian on backgroung
          Q=0.                      ! which is described by exponencial fcion
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=(PIH*Q+DMG*exp(-EGAM/PAR_M1(1)))*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.18) THEN   ! Classical Lorentzian on a background
          Q=0.                       ! given by exponencial fcion of Eg and Ei
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=EINI-EGAM
          SGAMMA=(PIH*Q+DMG*exp(-EINI/PAR_E1(3))
     *           * exp(-EGAM/PAR_M1(1)))*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
c

        ELSEIF (NOPTM1.EQ.31) THEN  ! SP on states
          Q=0.                      ! with Eini lower than Etr
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &         ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          IF (EINI.LE.PAR_M1(1)) THEN
            SGAMMA=(PIH*Q+DMG)*EGAM**3
          ELSE
            SGAMMA=PIH*Q*EGAM**3
          ENDIF
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.32) THEN   ! Classical Lorentzian ...
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=(PIH*Q+DMG*exp(-EINI/PAR_M1(1)))*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN


c
c
        ELSEIF (NOPTM1.EQ.8) THEN   ! SP + Double-bumped Scissors mode - higher bump only on states
          Q=0.                      ! with Efinal lower than Etr
          DO I=1,1
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          DO I=3,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF (Efinal.LE.PAR_M1(1)) THEN
            QQ=SIGM(2)*EGAM*WM0(2)**2/
     &        ((EGAM**2-ERM(2)**2)**2+(EGAM*WM0(2))**2)
            Q=Q+QQ
          ENDIF
c          SGAMMA=PIH*Q*EGAM**3
          SGAMMA=(PIH*Q+DMG)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
c
        ELSEIF (NOPTM1.EQ.9) THEN   ! Scissors mode is build up only on
          Q=0.                      ! some of the terminal states
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF ((Efinal.EQ.enrgfin(1)).OR.(Efinal.EQ.enrgfin(2)).OR.
c     &        (Efinal.EQ.enrgfin(3)).OR.(Efinal.EQ.enrgfin(4)).OR.
c     &        (Efinal.EQ.enrgfin(5)).OR.(Efinal.EQ.enrgfin(6)).OR.
     &        (Efinal.EQ.enrgfin(3))) THEN
            QQ=SIGM(1)*EGAM*WM0(1)**2/
     &        ((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
            Q=Q+QQ
          ENDIF
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.10) THEN   ! Scissors mode up to Etr 
          Q=0.                       ! on SP backgroung (everywhere)
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Efinal=Eini-Egam
          IF (Efinal.LE.PAR_M1(1)) THEN
            QQ=SIGM(1)*EGAM*WM0(1)**2/
     &        ((EGAM**2-ERM(1)**2)**2+(EGAM*WM0(1))**2)
            Q=Q+QQ
          ENDIF
          SGAMMA=(PIH*Q+DMG)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
c
        ELSEIF (NOPTM1.EQ.11) THEN                            ! Goriely
          SGAMMA = APSF(EGAM,2)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.12) THEN                            ! Goriely
          SGAMMA = APSF(EGAM,2)
          SGAMMA = SGAMMA + 
     *             PAR_M1(1) * EXP(-PAR_M1(2)*EGAM) *
     *             (1.0+PAR_M1(3)*EGAM*EGAM*EGAM)
          SGAMMA = SGAMMA * EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.13) THEN                            ! Goriely (incorrect Eg dep)
          SGAMMA = APSF(EGAM,2)
          SGAMMA = SGAMMA + 
     *          PAR_M1(1)*EXP(-PAR_M1(2)*EGAM)*(1.0+PAR_M1(3)*EGAM*EGAM)
          SGAMMA = SGAMMA * EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
c

        ELSEIF (NOPTM1.EQ.14) THEN   ! Width of Scissors resonance depends on Efin
          Efinal=Eini-Egam           ! => total strength increases
          WM0SC=WM0(1)
          WM0(1)=WM0(1)+PAR_M1(1)*Efinal
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          WM0(1)=WM0SC               !restoring the original value
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.15) THEN   ! Width of Scissors resonance depends on Efin
          Efinal=Eini-Egam           ! Total strength should remain 
          WM0SC=WM0(1)
          WM0(1)=WM0(1)+PAR_M1(1)*Efinal
          SIGMSC=SIGM(1)
          SIGM(1)=WM0SC*SIGMSC/WM0(1)
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          WM0(1)=WM0SC               !restoring the original values
          SIGM(1)=SIGMSC
          IF (ITYP.EQ.3) RETURN
c
        ELSEIF (NOPTM1.EQ.16) THEN   ! \sigma_0 of Scissors resonance depends on Efin
          Efinal=Eini-Egam            
          SIGMSC=SIGM(1)
          SIGM(1) = SIGMSC + PAR_M1(1) * Efinal
          IF (SIGM(1).LE.0.0) SIGM(1) = 0.0
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          SIGM(1)=SIGMSC               !restoring the original values
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.21) THEN   ! Oslo - factor k
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.22) THEN   ! Oslo - factor k + softpole
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Q=Q+PAR_M1(2)/(EGAM**PAR_M1(3))
          SGAMMA=PAR_M1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.23) THEN   ! Oslo - factor k + softpole with a const below Etr
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_M1(1)) EGAM=PAR_M1(1)
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          EGAM=EGAMOLD
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.24) THEN   ! SP - factor k + softpole 
          Q=0.
          EGAMOLD=EGAM
          IF (EGAM.LT.PAR_M1(1)) EGAM=PAR_M1(1)
          Q=Q+PAR_E1(2)/(EGAM**PAR_E1(3))
          EGAM=EGAMOLD
          SGAMMA=PAR_E1(1)*(PIH*Q+DMG)*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ELSEIF (NOPTM1.EQ.25) THEN   ! Oslo - factor k + softpole/rossendorf
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          QQ= PAR_M1(1) * EXP(-EGAM / PAR_E1(2)) 
          SGAMMA=PAR_E1(1) * (PIH*Q + QQ) * EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN
        ELSEIF (NOPTM1.EQ.26) THEN   ! Oslo - factor k + softpole/rossendorf (no PAR_E1(1))
          Q=0.
          DO I=1,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          QQ= PAR_M1(1) * EXP(-EGAM / PAR_E1(2)) 
          SGAMMA= (PIH*Q + QQ) * EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN


        ELSEIF (NOPTM1.EQ.27) THEN   ! Oslo - factor k not for LLR
          Q=0.
          DO I=1,1
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ/PAR_E1(1)
          ENDDO
          DO I=2,NGIGM
            QQ=SIGM(I)*EGAM*WM0(I)**2/
     &        ((EGAM**2-ERM(I)**2)**2+(EGAM*WM0(I))**2)
            Q=Q+QQ
          ENDDO
          SGAMMA=PAR_E1(1)*PIH*Q*EGAM**3
          SFCEM1=SGAMMA
          IF (ITYP.EQ.3) RETURN

        ENDIF
      ENDIF
c
c
c*****                       E2 component
c
      IF ((ITYP.EQ.2).OR.(ITYP.EQ.4)) THEN
        IF (NOPTE2.EQ.0) THEN
          SGAMMA=SGAMMA+QEL*EGAM**5
            SFCEE2=SGAMMA-SFCEM1
        ELSEIF (NOPTE2.EQ.11) THEN                            ! Goriely
          SGAMMA = SGAMMA + APSF(EGAM,3)*EGAM**5
          SFCEE2=SGAMMA-SFCEM1
          RETURN
        ELSEIF (NOPTE2.EQ.1) THEN                            ! Goriely
          Q=0.
          DO I=1,NGIGE2
            QQ=(SIGE(I)*WE0(I)**2)/(EGAM*((EGAM**2-ERE(I)**2)**2+
     &        (EGAM*WE0(I))**2))
            Q=Q+QQ
          ENDDO
          SGAMMA=SGAMMA+PIHQ*Q*EGAM**5
          SFCEE2=SGAMMA-SFCEM1
        ENDIF
        RETURN
      ENDIF
c
      END
c




