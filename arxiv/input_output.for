c***********************************************************************
SUBROUTINE READ_EV(NAME)
    c***********************************************************************
    
          INTEGER*4     MAXIS
          PARAMETER     (MAXJC  = 49,   !maximum spin generated in level scheme
         *               MAXBIN = 700,  !maximum # of bins in continuum
         *               MAXIS  = 76000000        )   
          LOGICAL        LOGS,lopopgs
          CHARACTER*80   TITLE1,TITLE2,TITLE3
          CHARACTER*12   NAME,NAME_RES
          INTEGER*4      IR1,IR2,IR3,IR4,STEPS,IRCON,IRCONc,NTOTAL,IRINIT
          INTEGER        DEKOD,DENUM,DELEV,DEPARITY
          REAL           Im2Res
    
          DIMENSION      ityp(2),qg(2),help(2),isbspin(2)      ! TNC
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
         &       /SFCE/  SFCEE1,SFCEM1,SFCEE2,GE1,GM1,GE2,
         &               GE1DIS(2),GM1DIS(2),GE2DIS(2)
         &       /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
         &               CONVk(0:1,5,50)
         &       /vars/  factnrm
         &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
         &               ISWIC,ISWMX,ISWWI,ISWLS
         &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)             ! TNC
          COMMON /GOE/   EIGENVAL(0:MAXBIN,0:MAXJC,0:1),
         *               NEIGENVAL(0:MAXJC,0:1),LMODE
         *       /denp/  LDENP,ZNUM,DENPC,DENPEN(3),IPDEN
          COMMON /LDTAB/ TABENLD(0:270),TABLD(0:270,0:MAXJC,0:1),NLD
          COMMON /PSFTAB/TABENPSF(3,0:400),TABPSF(3,0:400),NPSF(3)
         &       /rnd/   IRINIT(4,999)
    
          OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
    C     OPEN (UNIT=4,FILE='AUXI.DAT')
    C
    C     User's alphanumeric titles:
    C
          READ (5,100) TITLE1
          READ (5,100) TITLE2
          READ (5,100) TITLE3
      100 FORMAT (A80)
    C
    C     The regime of run:
    C
          READ (5,*)
          READ (5,*) IVER,IPRIM
          READ (5,*)
          READ (5,*) ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,ISWIC,ISWMX,ISWWI,ISWLS
          READ (5,*)
          READ (5,*) NBIN, (IRINIT(I,1),I=1,4)
          IR1 = IRINIT(1,1) 
          IR2 = IRINIT(2,1)                   ! not really necessary
          IR3 = IRINIT(3,1)
          IR4 = IRINIT(4,1)
          READ (5,*)
          READ (5,*) NOPTFL,NOPTE1,NOPTM1,NOPTE2,NOPTDE,LMODE,LDENP
          READ (5,*)
          Read (5,*) NSREAL, NEVENTS, NRL
    C
    C     Giant Resonaces:
    C
          READ (5,*)
          READ (5,*) NGIGE
          DO 15 I=1, NGIGE
            READ (5,*) ER(I),W0(I),SIG(I)
       15 CONTINUE
          READ (5,*)
          READ (5,*) NGIGM
          DO 16 I=1, NGIGM
            READ (5,*) ERM(I),WM0(I),SIGM(I)
       16 CONTINUE
          READ (5,*)
          READ (5,*) NGIGE2
          DO 17 I=1, NGIGE2
            READ (5,*) ERE(I),WE0(I),SIGE(I)
       17 CONTINUE
    C
    C     Other data needed for photon strength:
    C
          READ (5,*)
          READ (5,*) DEG,DMG,QEL
          READ (5,*)
          READ (5,*) FERMC, EK0, EGZERO 
          READ (5,*)
          READ (5,*) (PAR_E1(I),I=1,3)
            DIPSLP=(1.0-PAR_E1(3))/(PAR_E1(2)-PAR_E1(1))
            DIPZER=1.0-DIPSLP*PAR_E1(2)
          READ (5,*)
          READ (5,*) (PAR_M1(I),I=1,3)
    C
    C     Data needed for level density formulas:
    C
          READ (5,*)
          READ (5,*) ASHELL,EONE,TEMPER,EZERO,AMASS,ZNUM,PAIRING
          READ (5,*)
          READ (5,*) ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
          READ (5,*)
          READ (5,*) DENPC,(DENPEN(I),I=1,3)        !see PRC67, 015803 
          READ (5,*)
          READ (5,*) DENLO,DENHI,(DENPAR(I),I=1,4)
    C
    C     Data relating to the (thermal) neutron capturing state:
    C
          READ (5,*)
          IF (IVER.EQ.0) THEN
            READ (5,*) BN,SPINCS,IPINC 
          ELSEIF (IVER.EQ.1) THEN
            READ (5,*) BN,SPINt,IPINC,CAPFR(1)
          ELSEIF (IVER.EQ.2) THEN
            READ (5,*) BN,SPINt,IPINC, En
            NAME_RES = 'RES_PAR.DAT'
            CALL READ_RES_PAR(NAME_RES,En,SPINt)
          ENDIF
          IF (IVER.EQ.0) THEN        ! IVER - resonances
            SPINc(1)=SPINCS
            NOPTCS=1
            NLINc=1
            CAPFR(1)=1.
          ELSEIF (IVER.EQ.1) THEN    ! TNC with fraction / for resonances already calculated in READ_RES_PAR()
           IF (SPINt.NE.0.) THEN
             NLINc=2
             NOPTCS=2
             spinc(1)=spint-0.5
             spinc(2)=spint+0.5
             capfr(2)=1.-capfr(1)
           ELSE
             NLINc=1
             NOPTCS=1
             spinc(1)=spint+0.5
             capfr(1)=1.
             capfr(2)=0.
           ENDIF
          ENDIF
    C
    C     Now read the tables of ICC-coefficients and other information
    C     on electron conversion. The used notation:
    C
    C             XRAYK,XRAYL         -- electron binding energy for K-shell
    C                                    (and L-shell) expressed in MeV
    C             NENT,NENK           -- the number of electron energy points
    C             NMU                 -- the highest multipolarity, e.g. 3 for
    C                                    octupole (electric and magnetic)
    C             CONVK(I,MU,K)       -- K-shell ICC for I-th type of radiation
    C                                    (0 for electric, 1 for magnetic),
    C                                    multipolarity MU and K-th electron
    C                                    energy value
    C             CONVT(I,M,K)        -- total ICC ... (as the CONVK)
    C             ELENT(K),ELENK(K)   -- K-th value of electron energy
    C
          read (5,*)
          read (5,*) xrayk,xrayl
          READ (5,*) NMU
          READ (5,*) NENT
          READ (5,*) (ELENT(K),K=1,NENT)
    c      READ (5,*) (((CONVT(ITY,MU,K),MU=1,NMU),ITY=0,1),K=1,NENT)
          READ (5,*) (((CONVT(ITY,MU,K),ITY=0,1),MU=1,NMU),K=1,NENT)
          READ (5,*) NENK
          READ (5,*) (ELENK(K),K=1,NENK)
    c      READ (5,*) (((CONVK(ITY,MU,K),MU=1,NMU),ITY=0,1),K=1,NENK)
          READ (5,*) (((CONVK(ITY,MU,K),ITY=0,1),MU=1,NMU),K=1,NENK)
    C
    C     Data related to the discrete levels (J, pi, Eexc, primary intensities
    C     and branchings):
    C
          do i=1,2
             GE1DIS(i)=0.
             GM1DIS(i)=0.
             GE2DIS(i)=0.
          enddo
          IPDEN = 0
    c
          READ (5,*)
          READ (5,*) ECRIT
          EALL=ECRIT
          DELTA=(BN-ECRIT)/FLOAT(NBIN)
    c      IF (IPRIM.GE.1) THEN   ! IPRIM
            READ (5,*)
            READ (5,*) FACTNRM
    c      ENDIF
    
          auxi=0
          read (5,*)
          READ (5,*) numlev
          DO i=1,numlev
           IF (IPRIM.EQ.0) THEN
             READ (5,*) enrg,spfi,ipfi,denum(i)
           ELSE 
             READ (5,*) enrg,spfi,ipfi,denum(i),prim,errprim
           ENDIF
           IF (denum(i).GT.20) THEN
             WRITE (*,*) 'number of de-exciting transitions is higher than',
         &               ' 20 (for level)',enrg
             STOP 'ERROR'  
           ENDIF
           IF (ndis(ISUBSC(spfi),ipfi).EQ.30) THEN
             WRITE (*,*) 'number of discrete levels with one spin/parity ',
         &               'is higher than 30',SPFI,IPFI
             STOP 'ERROR'  
           ENDIF
           ndis(ISUBSC(spfi),ipfi)=ndis(ISUBSC(spfi),ipfi)+1
           endis(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=enrg
           dekod(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=i
           IF (I.LE.5) THEN
             IPDEN = IPDEN + ipfi
           ENDIF
           IF (denum(i).GT.0) then
            DO k=1,denum(i)
             READ (5,*) enrgf,sal,errsal,desp,ipar,dlt,alpha
             if (sal.LE.0.0) then
               WRITE(*,*)'intensity must be >0 (between levels)',enrg,enrgf
               STOP 'ERROR'
             endif 
             DO il=1,ndis(ISUBSC(desp),ipar)
               IF (ABS(endis(il,ISUBSC(desp),ipar)-enrgf).LE.1e-5) THEN
                delev(i,k)=il
                despin(i,k)=desp
                deparity(i,k)=ipar
                deltx(i,k)=dlt
                dmixsq=dlt**2
                if (alpha.LE.1e-6) alpha=ALPH_TOT(enrg,spfi,ipfi,enrgf,
         *                            desp,ipar,dmixsq,nent,elent,convt)
                alphak=ALPH_TOT(enrg,spfi,ipfi,enrgf,desp,ipar,dmixsq,
         *                     nenk,elenk,convk)
                if (alphak.GT.alpha) alphak=alpha
                ponv(i,k)=alpha/(1+alpha)
                ponvk(i,k)=alphak/(1+alpha)
                sall(i,k)=sall(i,k-1)+sal*(1+alpha)
               ENDIF
             ENDDO
             if (delev(i,k).EQ.0) then
               WRITE(*,*) 'Wrong final level (Eini, Efin)',enrg,enrgf
               STOP 'ERROR'
             endif
            ENDDO
           ENDIF
           elowlev(i)=enrg
           IF(elowlev(i).EQ.0) THEN
            LOpopGS=.TRUE.
            KpopGS=i
           ENDIF
           elowsp(i)=spfi
           ilowip(i)=ipfi
           IF (IPDEN.GT.2) THEN 
             IPDEN = 1
           ELSE
             IPDEN = 0
           ENDIF
    c
    c     Intensities of primary transitions
    c
           IF (IPRIM.GE.1) THEN
            DO ILINc=1,NLINc
              isbspin(ilinc)=NINT(spfi+.25)-NINT(spinc(ILINc)+.25)
              ityp(ilinc)=ITYPE(spinc(ilinc),ipinc,spfi,ipfi)
            ENDDO
            IF (NLINc.EQ.1) ityp(2)=0
            IF ((ityp(1).EQ.0).AND.(ityp(2).EQ.0)) THEN
              IF (prim.GT.0) THEN
                write(*,*) "No possible transition for primary feeding"
                write(*,*) "for level",i," with spin/parity",spfi,ipfi
                STOP
              ENDIF
            ELSEIF ((ityp(1).NE.0).AND.(ityp(2).EQ.0)) THEN
             STDISa(1,ndis(ISUBSC(spfi),ipfi),isbspin(1),ipfi)=
         &         STDISa(1,ndis(ISUBSC(spfi),ipfi)-1,isbspin(1),ipfi)
         &         +prim*factnrm/capfr(1)
             IF (ityp(1).EQ.1) ge1dis(1)=ge1dis(1)+prim*factnrm/capfr(1)
             IF ((ityp(1).EQ.3).OR.(ityp(1).EQ.2)) gm1dis(1)=gm1dis(1)
         &                                +prim*factnrm/capfr(1)
             IF (ityp(1).EQ.4) ge2dis(1)=ge2dis(1)+prim*factnrm/capfr(1)
            ELSEIF ((ityp(1).EQ.0).AND.(ityp(2).NE.0)) THEN
              STDISa(2,ndis(ISUBSC(spfi),ipfi),isbspin(2),ipfi)=
         &         STDISa(2,ndis(ISUBSC(spfi),ipfi)-1,isbspin(2),ipfi)
         &         +prim*factnrm/capfr(2)
             IF (ityp(2).EQ.1) ge1dis(2)=ge1dis(2)+prim*factnrm/capfr(2)
             IF ((ityp(2).EQ.3).OR.(ityp(2).EQ.2)) gm1dis(2)=gm1dis(2)
         &                                +prim*factnrm/capfr(2)
             IF (ityp(2).EQ.4) ge2dis(2)=ge2dis(2)+prim*factnrm/capfr(2)
            ELSEIF ((ityp(1).NE.0).AND.(ityp(2).NE.0)) THEN
             DO ILINc=1,NLINc
               IF (prim.GT.0) THEN
                 QG(ilinc)=sgamma(bn-enrg,bn,ityp(ilinc))/
         *                 density(bn,spinc(ilinc),ipinc)/(prim*factnrm)
               ELSE
                 QG(ilinc)=0.
               ENDIF
             ENDDO
             qgsum=0.
             DO ilinc=1,nlinc
               qgsum=qgsum+qg(ilinc)
             ENDDO
             DO ilinc=1,nlinc
                if (qgsum.GT.0) then
                  help(ilinc)=QG(ilinc)/qgsum*prim*factnrm/capfr(ilinc)
                else
                  help(ilinc)=0.
                endif
               ENDDO
             DO ilinc=1,nlinc
               STDISa(ilinc,ndis(ISUBSC(spfi),ipfi),isbspin(ilinc),ipfi)=
         &       STDISa(ilinc,ndis(ISUBSC(spfi),ipfi)-1,isbspin(ilinc),ipfi)
         &       +help(ilinc)
               IF (ityp(ilinc).EQ.1)
         &       ge1dis(ilinc)=ge1dis(ilinc)+help(ilinc)
               IF ((ityp(ilinc).EQ.3).OR.(ityp(ilinc).EQ.2))
         &       gm1dis(ilinc)=gm1dis(ilinc)+help(ilinc)
               IF (ityp(ilinc).EQ.4)
         &       ge2dis(ilinc)=ge2dis(ilinc)+help(ilinc)
             ENDDO
            ENDIF        ! ityp
    c      
           ENDIF ! IPRIM
    c
          ENDDO  ! i - NUMLEV
    c
          DO 31 ipfi=0,1
           sp=spinc(1)-int(spinc(1)+.25)
           DO 31 ispfi = 0, 8
            spfi = sp + float(ispfi)
            levdis(ISUBSC(spfi),ipfi)=auxi+ndis(ISUBSC(spfi),ipfi)
            auxi=levdis(ISUBSC(spfi),ipfi)
       31   CONTINUE
          nddd=auxi
    c
    c      KGS=0
    c      DO I=1,NFILEV
    c        READ (5,*) ENRGFIN(I)
    c        IF (ENRGFIN(I).EQ.0.) THEN
    c          LOGS=.TRUE.
    c          KGS=I
    c        ENDIF
    c        EDO(I)=ENRGFIN(I)-1.E-6
    c        EUP(I)=ENRGFIN(I)+1.E-6
    c      ENDDO
    c
          CLOSE (UNIT=5,STATUS='KEEP')
    c
    c     Tabulated level density (Goriely)
    c
          IF (NOPTDE.EQ.11) THEN
            NAME = "LDTAB.DAT"
            OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
            READ(5,*) 
            READ(5,*) NLD, SPACRES,SPINi,IPINi,corAlpLD,corDelLD
            READ(5,*) 
            DO IEN = 1, NLD
              READ(5,*) TABENLD(IEN),DUMMY,DUMMY,DUMMY,DUMMY,
         *              (TABLD(IEN,J,0),J=0,MAXJC) 
            ENDDO
            READ(5,*) 
            DO IEN = 1, NLD
              READ(5,*) TABENLD(IEN),DUMMY,DUMMY,DUMMY,DUMMY,
         *              (TABLD(IEN,J,1),J=0,MAXJC) 
            ENDDO
            CLOSE(5)
    c
    c     "Normalization to experimental level spacing" - SPINt/IPINt not needed at the moment
    c      PRC78,064307(2008): \rho(u,J,p)=\exp(alpha x \sqrt(U-delta)) * \rho(U-delta,J,p)
    c
            DO I = 1, NLD
              FSPAC = EXP(corAlpLD * sqrt(TABENLD(I)))
              DO J = 0, MAXJC
               DO K = 0, 1
                 TABLD(I,J,K) = TABLD(I,J,K) * FSPAC
               ENDDO
              ENDDO
              TABENLD(I)=TABENLD(I)+corDelLD
            ENDDO
          ENDIF 
    c
    c  Normalization is not "correct": PRC78,064307(2008) gives \rho(u,J,p)=\exp(alpha x \sqrt(U-delta)) * \rho(U-delta,J,p)
    c  parameters in RIPL correcton file are (likely) alpha = ctable and delta = ptable
    c
    
    c
    c     Tabulated level density (Kawano)
    c
          IF (NOPTDE.EQ.12) THEN
            NAME = "LDTAB_K.DAT"
            OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
            READ(5,*) 
            READ(5,*) NLD
            READ(5,*) 
            DO IEN = NLD, 1, -1
              READ(5,*) TABENLD(IEN),DUMMY,
         *              (TABLD(IEN,J,0),J=0,9) 
            ENDDO
            DO IEN = 1, NLD
              DO J = 0, 9
                TABLD(IEN,J,0) = TABLD(IEN,J,0) / 2.0 ! Toshihiko gives total NLD (p-independent)
                TABLD(IEN,J,1) = TABLD(IEN,J,0) 
              ENDDO
            ENDDO
          ENDIF
    c
    c     Tabulated PSF
    c
          IF ((NOPTE1.GE.11).AND.(NOPTE1.LE.15)) THEN
            NAME = "PSFE1.DAT"
            OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
            READ(5,*) 
            READ(5,*) NPSF(1)
            READ(5,*) 
            DO IEN = 1, NPSF(1)
              READ(5,*) TABENPSF(1,IEN),TABPSF(1,IEN) 
            ENDDO
            CLOSE(5)
          ENDIF
    
          IF ((NOPTM1.GE.11).AND.(NOPTM1.LE.15)) THEN
            NAME = "PSFM1.DAT"
            OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
            READ(5,*) 
            READ(5,*) NPSF(2)
            READ(5,*) 
            DO IEN = 1, NPSF(2)
              READ(5,*) TABENPSF(2,IEN),TABPSF(2,IEN) 
            ENDDO
            CLOSE(5)
          ENDIF
    
          IF ((NOPTE2.EQ.11).AND.(NOPTE2.LE.15)) THEN
            NAME = "PSFE2.DAT"
            OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
            READ(5,*) 
            READ(5,*) NPSF(3)
            READ(5,*) 
            DO IEN = 1, NPSF(3)
              READ(5,*) TABENPSF(3,IEN),TABPSF(3,IEN) 
            ENDDO
            CLOSE(5)
          ENDIF
          WRITE(*,*) ' Inputs have been successfully loaded.'
    
          RETURN
          END
    c
    c***********************************************************************
          SUBROUTINE READ_INT(NAME)
    c***********************************************************************
    c
          INTEGER*4    MAXIS
          REAL         Im2Res
          PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
         *              MAXBIN = 700,  !maximum # of bins in continuum
         *              MAXIS  = 76000000        )   
          CHARACTER*12   NAME
          INTEGER*4      IR1,IR2,IR3,IR4,STEPS,IRCON,IRCONc,NTOTAL
          INTEGER        DEKOD,DENUM,DELEV,DEPARITY
    c
    C     Hereafter, article 'CON' is related to a 'CONtinuum' part of
    C     the level system, while 'DIS' -- to the 'DIScrete' part.
    C     Suffixes 'p' and 's' are related to primary and secondary
    C     transitions, respectively.
    C
          COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
         &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
         &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
         &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
         &               denum(99),delev(99,20),despin(99,20),
         &               deparity(99,20),irconc(2),deltx(99,20)
         &       /PHOTO/ NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
         &               SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
         &               DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
         &               PAR_E1(3),DIPSLP,DIPZER
         &       /readi/ nevents,NSREAL,NRL,ecrit,xrayk,xrayl,numlev,
         &               elowlev(99),elowsp(99),ilowip(99),
         &               STDISa(0:2,0:30,-2:2,0:1),RADW
         &       /intr/  last,ir1,ir2,ir3,ir4,radwid,radwdi,
         &               fnoav,fnodi,rnn(18,14),
         &               rmult(18,126),rovfl(18),avea(18,14),erra(18,14),
         &               puvea(18,14),cova(18,18,14),aves(1:18,0:MAXBIN),
         &               errs(1:18,0:MAXBIN),
         &               diss(1:18,0:MAXBIN),poptlev(99),popslev(99)
         &       /incr/  nfilev,STEPS,enrgfin(18),buff(251),
         &               ncum,cumwidth(14),edo(18),eup(18)
         &       /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,NTOTAL,iregi,
         &               NOPTCS,NLINc,CAPFR(2)
         &       /angle/ spinc(2),ipinc,spintl(18),ipintl(18),
         &               smear,Fk(4,4,0:16,0:16),F4(0:16,0:16)
         &       /SFCE/  SFCEE1,SFCEM1,SFCEE2,GE1,GM1,GE2,
         &               GE1DIS(2),GM1DIS(2),GE2DIS(2)
         &       /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
         &               CONVk(0:1,5,50)
         &       /vars/  factnrm
         &       /GAU/   U,IFLAG
         &       /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
         *               DENLO,DENHI,DENPAR(4),
         *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
         &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
         &               ISWIC,ISWMX,ISWWI,ISWLS
         &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)             ! TNC
    c
          iflag=0
          DO J=0,MAXJC
            DO K=0,1
              NDIS(J,K)=0
            ENDDO
          ENDDO
    c
          OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
    c
    c         Reading of unimportant variables
    c
          DO k=1,15
            READ (5,*)
          ENDDO
          DO 15 I=1, NGIGE
            READ (5,*)
       15 CONTINUE
          READ (5,*)
          READ (5,*)
          DO 16 I=1, NGIGM
            READ (5,*)
       16 CONTINUE
          READ (5,*)
          READ (5,*)
          DO 17 I=1, NGIGE2
            READ (5,*)
       17 CONTINUE
          DO k=1,18
            READ (5,*)
          ENDDO
    c
    c     Conversion coefficients
    c
          read (5,*)
          read (5,*) xrayk,xrayl
          READ (5,*) NMU
          READ (5,*) NENT
          READ (5,*) (ELENT(K),K=1,NENT)
    c      READ (5,*) (((CONVT(ITY,MU,K),MU=1,NMU),ITY=0,1),K=1,NENT)
          READ (5,*) (((CONVT(ITY,MU,K),ITY=0,1),MU=1,NMU),K=1,NENT)
          READ (5,*) NENK
          READ (5,*) (ELENK(K),K=1,NENK)
    c      READ (5,*) (((CONVK(ITY,MU,K),MU=1,NMU),ITY=0,1),K=1,NENK)
          READ (5,*) (((CONVK(ITY,MU,K),ITY=0,1),MU=1,NMU),K=1,NENK)
    c
    c        Intensities of transitions not fixed - error taken into account)
    c
    c      IF (IPRIM.EQ.0) THEN
    c       DO k=1,3
    c        READ (5,*)
    c       ENDDO
    c      ELSE
           DO k=1,5
            READ (5,*)
           ENDDO
    c      ENDIF
    
          READ (5,*) numlev
          DO i=1,numlev
           IF (IPRIM.EQ.0) THEN
             READ (5,*) enrg,spfi,ipfi,denum(i)
           ELSE 
             READ (5,*) enrg,spfi,ipfi,denum(i),prim,errprim
           ENDIF
           ndis(ISUBSC(spfi),ipfi)=ndis(ISUBSC(spfi),ipfi)+1
           endis(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=enrg
           dekod(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=i
           IF (denum(i).GT.0) then
            DO k=1,denum(i)
             READ (5,*) enrgf,sal,errsal,desp,ipar,dlt,alpha
             DO il=1,ndis(ISUBSC(desp),ipar)
               IF (endis(il,ISUBSC(desp),ipar).EQ.enrgf) THEN
                delev(i,k)=il
                despin(i,k)=desp
                deparity(i,k)=ipar
                deltx(i,k)=dlt
                dmixsq=dlt**2
                if (alpha.LE.1e-6) alpha=ALPH_TOT(enrg,spfi,ipfi,enrgf,
         *                            desp,ipar,dmixsq,nent,elent,convt)
                alphak=ALPH_TOT(enrg,spfi,ipfi,enrgf,desp,ipar,dmixsq,
         *                     nenk,elenk,convk)
                if (alphak.GT.alpha) alphak=alpha
                ponv(i,k)=alpha/(1+alpha)
                ponvk(i,k)=alphak/(1+alpha)
       51       sal=sal+errsal*GAUSS(IR4)             !!! not fixed
                if (sal.LE.0) goto 51
                sall(i,k)=sall(i,k-1)+sal*(1+alpha)
               ENDIF
             ENDDO
            ENDDO
           ENDIF
          ENDDO
    c
          CLOSE(5)
          RETURN
          END
    c
    c***********************************************************************
          SUBROUTINE READ_AGAIN(NAME)
    c***********************************************************************
    c
          INTEGER*4    MAXIS
          PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
         *              MAXBIN = 700,  !maximum # of bins in continuum
         *              MAXIS  = 76000000        )   
          CHARACTER*12   NAME
          INTEGER*4      IR1,IR2,IR3,IR4,STEPS,IRCON,IRCONc,NTOTAL
          INTEGER        DEKOD,DENUM,DELEV,DEPARITY
          REAL           Im2Res
          DIMENSION      QG(2),ITYP(2),ISBSPIN(2)
    
          CHARACTER*80  xxx2
    
    c
    C     Hereafter, article 'CON' is related to a 'CONtinuum' part of
    C     the level system, while 'DIS' -- to the 'DIScrete' part.
    C     Suffixes 'p' and 's' are related to primary and secondary
    C     transitions, respectively.
    C
          COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
         &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
         &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
         &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
         &               denum(99),delev(99,20),despin(99,20),
         &               deparity(99,20),irconc(2),deltx(99,20)
         &       /PHOTO/ NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
         &               SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
         &               DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
         &               PAR_E1(3),DIPSLP,DIPZER
         &       /readi/ nevents,NSREAL,NRL,ecrit,xrayk,xrayl,numlev,
         &               elowlev(99),elowsp(99),ilowip(99),
         &               STDISa(0:2,0:30,-2:2,0:1),RADW
         &       /intr/  last,ir1,ir2,ir3,ir4,radwid,radwdi,
         &               fnoav,fnodi,rnn(18,14),
         &               rmult(18,126),rovfl(18),avea(18,14),erra(18,14),
         &               puvea(18,14),cova(18,18,14),aves(1:18,0:MAXBIN),
         &               errs(1:18,0:MAXBIN),
         &               diss(1:18,0:MAXBIN),poptlev(99),popslev(99)
         &       /incr/  nfilev,STEPS,enrgfin(18),buff(251),
         &               ncum,cumwidth(14),edo(18),eup(18)
         &       /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,NTOTAL,iregi,
         &               NOPTCS,NLINc,CAPFR(2)
         &       /angle/ spinc(2),ipinc,spintl(18),ipintl(18),
         &               smear,Fk(4,4,0:16,0:16),F4(0:16,0:16)
         &       /SFCE/  SFCEE1,SFCEM1,SFCEE2,GE1,GM1,GE2,
         &               GE1DIS(2),GM1DIS(2),GE2DIS(2)
         &       /pver/  NENT,ELENt(50),CONVt(0:1,5,50),NENK,ELENk(50),
         &               CONVk(0:1,5,50)
         &       /GAU/   U,IFLAG
         &       /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
         *               DENLO,DENHI,DENPAR(4),
         *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
         &       /vars/  factnrm
         &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
         &               ISWIC,ISWMX,ISWWI,ISWLS
         &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)
    c
    c         Preparing of some variables
    c
          DO K=0,1
            DO L=0,20
              DO M=-2,2
                DO MM=0,2
                  STDISa(MM,L,M,K)=0.
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          DO J=0,MAXJC
            DO K=0,1
              NDIS(J,K)=0
            ENDDO
          ENDDO
          DO i=1,2
            GE1DIS(i)=0.
            GM1DIS(i)=0.
            GE2DIS(i)=0.
          ENDDO
    c
          OPEN (UNIT=5,FILE=NAME,STATUS='OLD')
    c
    c         Reading of unimportant variables
    c
          DO k=1,15
            READ (5,*)
          ENDDO
          DO 15 I=1, NGIGE
            READ (5,*)
       15 CONTINUE
          READ (5,*)
          READ (5,*)
          DO 16 I=1, NGIGM
            READ (5,*)
       16 CONTINUE
          READ (5,*)
          READ (5,*)
          DO 17 I=1, NGIGE2
            READ (5,*)
       17 CONTINUE
          DO k=1,18
            READ (5,*) xxx2
          ENDDO
    C
    C     Data about resonances
    C
    c
    c     Conversion coefficients
    c
          read (5,*)
          read (5,*) xrayk,xrayl
          READ (5,*) NMU
          READ (5,*) NENT
          READ (5,*) (ELENT(K),K=1,NENT)
    c      READ (5,*) (((CONVT(ITY,MU,K),MU=1,NMU),ITY=0,1),K=1,NENT)
          READ (5,*) (((CONVT(ITY,MU,K),ITY=0,1),MU=1,NMU),K=1,NENT)
          READ (5,*) NENK
          READ (5,*) (ELENK(K),K=1,NENK)
    c      READ (5,*) (((CONVK(ITY,MU,K),MU=1,NMU),ITY=0,1),K=1,NENK)
          READ (5,*) (((CONVK(ITY,MU,K),ITY=0,1),MU=1,NMU),K=1,NENK)
    c
    c        Computing of primary intensities corresponding to each capture state
    c
          DO k=1,5
            READ (5,*)
          ENDDO
          READ (5,*) numlev
          DO i=1,numlev
           READ (5,*) enrg,spfi,ipfi,denum(i),prim,errprim
           ndis(ISUBSC(spfi),ipfi)=ndis(ISUBSC(spfi),ipfi)+1
           endis(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=enrg
           dekod(ndis(ISUBSC(spfi),ipfi),ISUBSC(spfi),ipfi)=i
           IF (denum(i).GT.0) then
            DO k=1,denum(i)
             READ (5,*) enrgf,sal,errsal,desp,ipar,dlt,alpha
            ENDDO
           ENDIF
           DO ILINc=1,NLINc
             isbspin(ilinc)=NINT(spfi+.25)-NINT(spinc(ILINc)+.25)
             ityp(ilinc)=ITYPE(spinc(ilinc),ipinc,spfi,ipfi)
           ENDDO
           IF ((ityp(1).EQ.0).AND.(ityp(2).EQ.0)) THEN
             IF (prim.GT.0) THEN
               write(*,*) "No possible transition"
               STOP
             ENDIF
           ELSEIF ((ityp(1).NE.0).AND.(ityp(2).EQ.0)) THEN
            STDISa(1,ndis(ISUBSC(spfi),ipfi),isbspin(1),ipfi)=
         &        STDISa(1,ndis(ISUBSC(spfi),ipfi)-1,isbspin(1),ipfi)
         &        +prim*factnrm/capfr(1)
            IF (ityp(1).EQ.1) ge1dis(1)=ge1dis(1)+prim*factnrm/capfr(1)
            IF ((ityp(1).EQ.3).OR.(ityp(1).EQ.2)) gm1dis(1)=gm1dis(1)
         &                               +prim*factnrm/capfr(1)
            IF (ityp(1).EQ.4) ge2dis(1)=ge2dis(1)+prim*factnrm/capfr(1)
           ELSEIF ((ityp(1).EQ.0).AND.(ityp(2).NE.0)) THEN
             STDISa(2,ndis(ISUBSC(spfi),ipfi),isbspin(2),ipfi)=
         &        STDISa(2,ndis(ISUBSC(spfi),ipfi)-1,isbspin(2),ipfi)
         &        +prim*factnrm/capfr(2)
            IF (ityp(2).EQ.1) ge1dis(2)=ge1dis(2)+prim*factnrm/capfr(2)
            IF ((ityp(2).EQ.3).OR.(ityp(2).EQ.2)) gm1dis(2)=gm1dis(2)
         &                               +prim*factnrm/capfr(2)
            IF (ityp(2).EQ.4) ge2dis(2)=ge2dis(2)+prim*factnrm/capfr(2)
           ELSEIF ((ityp(1).NE.0).AND.(ityp(2).NE.0)) THEN
             QGnfl=0.
             DO ILINc=1,NLINc
              if (prim.GT.0.) then
                QGnfl=QGnfl+capfr(ilinc)*sgamma(bn-enrg,bn,ityp(ilinc))/
         *            density(bn,spinc(ilinc),ipinc)/(prim*factnrm)
              else
                QGnfl=1.
              endif
             ENDDO
      151    CONTINUE
             DO ILINc=1,NLINc
              qg(ilinc)=0.
              iflag=0
              IF (prim.GT.0) THEN
                IF (IVER.LE.1) THEN
                  G1=GAUSS(ir4)
                  QG(ilinc)=sgamma(bn-enrg,bn,ityp(ilinc))*G1*G1/
         *                  density(bn,spinc(ilinc),ipinc)/(prim*factnrm)
                ELSEIF (IVER.EQ.2) THEN
                  G1=GAUSS(ir4)
                  G2=GAUSS(ir4)+CORRI(ilinc)*G1
                  GSQ=(Re2Res(ilinc)*G2*G2/(1+CORRI(ilinc)**2)+
         *            Im2Res(ilinc)*G1*G1)/(Re2Res(ilinc)+Im2Res(ilinc))
                  QG(ilinc)=sgamma(bn-enrg,bn,ityp(ilinc))*GSQ/
         *                  density(bn,spinc(ilinc),ipinc)/(prim*factnrm)
                ENDIF
              ELSE
                QG(ilinc)=0.
              ENDIF
             ENDDO
             IF (((abs(qg(1)*capfr(1)+qg(2)*capfr(2)-QGnfl))/QGnfl.LE.0.01)
         &      .OR.(prim.EQ.0)) THEN
              DO ilinc=1,nlinc
               STDISa(ilinc,ndis(ISUBSC(spfi),ipfi),isbspin(ilinc),ipfi)=
         &       STDISa(ilinc,ndis(ISUBSC(spfi),ipfi)-1,isbspin(ilinc),ipfi)
         &       +prim*factnrm*QG(ilinc)/QGnfl
               IF (ityp(ilinc).EQ.1) ge1dis(ilinc)=ge1dis(ilinc)+
         &                                  prim*factnrm*QG(ilinc)/QGnfl
               IF ((ityp(ilinc).EQ.3).OR.(ityp(ilinc).EQ.2))
         &          gm1dis(ilinc)=gm1dis(ilinc)+prim*factnrm*QG(ilinc)/QGnfl
               IF (ityp(ilinc).EQ.4) ge2dis(ilinc)=ge2dis(ilinc)+
         &                                  prim*factnrm*QG(ilinc)/QGnfl
              ENDDO
              GOTO 152
             ENDIF
             GOTO 151
           ENDIF
      152  CONTINUE
          ENDDO
          CLOSE(5)
          RETURN
          END
    c
    c***********************************************************************
          SUBROUTINE READ_RES_PAR(NAME,En,SPINt)
    c***********************************************************************
    c
          REAL           Im2Res
          CHARACTER*12   NAME
          DIMENSION      ERES(50,2),GRED(50,2),GGAM(50,2),NRES(2),
         &               gf(2),sg(2),rho(2),RI(2)
         
          COMMON /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,NTOTAL,iregi,
         &               NOPTCS,NLINc,CAPFR(2)
         &       /angle/ spinc(2),ipinc,spintl(18),ipintl(18),
         &               smear,Fk(4,4,0:16,0:16),F4(0:16,0:16)
         &       /reson/ Re2Res(2),Im2Res(2),CORRI(2)
    C
    C     Data about resonances
    C
          OPEN(UNIT=6,FILE=NAME,STATUS='OLD')
          READ(6,*)
          READ(6,*) NRES(1)
          DO l=1,NRES(1)
           READ(6,*) ERES(l,1),GRED(l,1),GGAM(l,1)
          ENDDO
          IF (SPINt.NE.0.) THEN
            READ(6,*)
            READ(6,*) NRES(2)
            DO l=1,NRES(2)
             READ(6,*) ERES(l,2),GRED(l,2),GGAM(l,2)
            ENDDO
          ENDIF
          CLOSE(6)
    c
    c     Auxiliary variables connected with resonances
    c
          DO j=1,2
           Re2Res(j)=0.
           Im2Res(j)=0.
           RI(j)=0.
           sg(j)=0.
          ENDDO
          IF (SPINt.NE.0) THEN
           jmax=2
           gf(1)=0.5*(2*(spint-0.5)+1)/(2*spint+1)
           gf(2)=0.5*(2*(spint+0.5)+1)/(2*spint+1)
          ELSE
           jmax=1
           gf(1)=0.5*(2*(spint+0.5)+1)/(2*spint+1)
           gf(2)=gf(1)
          ENDIF
    c
          DO j=1,jmax
           DO l=1,NRES(j)
            gtot=(GGAM(l,j)+sqrt(En)/2/gf(j)*GRED(l,j))
            denom=((En-ERES(l,j))**2+1.e-6*gtot**2./4)
            rnoma=(ERES(l,j)-En)*sqrt(GRED(l,j)/2./gf(j))
            rnomb=1.e-3*gtot/2.*sqrt(GRED(l,j)/2./gf(j))
            Re2Res(j)=Re2Res(j)+(rnoma/denom)**2
            Im2Res(j)=Im2Res(j)+(rnomb/denom)**2
            RI(j)=RI(j)+rnoma/denom*rnomb/denom
            sigma=GGAM(l,j)*GRED(l,j)/2/((En-ERES(l,j))**2+
         *        1e-6/4*(GGAM(l,j)+sqrt(En)/2/gf(j)*GRED(l,j))**2)
            sg(j)=sg(j)+sigma
           ENDDO
           rho(j)=RI(j)/(sqrt(Re2Res(j))*sqrt(Im2Res(j)))
           CORRI(j)=rho(j)/sqrt(1-rho(j)**2)
          ENDDO
          capfr(1)=sg(1)/(sg(1)+sg(2))
          capfr(2)=1.-capfr(1)
    c
    c      write(*,*) 'Capturing spin fractions:',capfr(1),capfr(2)
    c
          IF (jmax.EQ.1) THEN
           NLINc=1
           NOPTCS=1
           spinc(1)=spint+0.5
           spinc(2)=spint+1.5
          ELSEIF (capfr(2).LT.0.001) THEN    ! One per cent contribution is the minimum
           capfr(1)=1.                      ! for consideration of capture spin
           capfr(2)=0.
           NLINc=1
           NOPTCS=1
           spinc(1)=spint-0.5
           spinc(2)=spint+0.5
          ELSEIF (capfr(1).LT.0.001) THEN
           capfr(1)=1.
           capfr(2)=0.
           NLINc=1
           NOPTCS=1
           spinc(1)=spint+0.5
           spinc(2)=spint+1.5
           Re2Res(1)=Re2Res(2)
           Im2Res(1)=Im2Res(2)
           CORRI(1)=CORRI(2)
          ELSE
           NLINc=2
           NOPTCS=2
           spinc(1)=spint-0.5
           spinc(2)=spint+0.5
          ENDIF
    
          RETURN
          END
    c

    c************************************************************************
        SUBROUTINE WRITE_PARAMS(NAME)
c************************************************************************
      INTEGER*4    MAXIS
      PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
     *              MAXBIN = 700,  !maximum # of bins in continuum
     *              MAXIS  = 76000000        )   
      CHARACTER*12   NAME
      CHARACTER*1    CHPAR
      COMMON /PHOTO/ NGIGE,NGIGM,NGIGE2,ER(5),SIG(5),W0(5),ERM(5),
     &               SIGM(5),WM0(5),ERE(5),SIGE(5),WE0(5),FERMC,DEG,
     &               DMG,QEL,NOPTE1,NOPTM1,NOPTE2,EK0,EGZERO,PAR_M1(3),
     &               PAR_E1(3),DIPSLP,DIPZER
     &       /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
     *               DENLO,DENHI,DENPAR(4),
     *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
     &       /MAINE/ eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,NTOTAL,iregi,
     &               NOPTCS,NLINc,CAPFR(2)
     &       /vars/  factnrm
     &       /angle/ spinc(2),ipinc,spintl(18),ipintl(18),
     &               smear,Fk(4,4,0:16,0:16),F4(0:16,0:16)
c
c       Writing down some inoput parameters to output file
c
      OPEN(9,FILE=NAME,STATUS='UNKNOWN')
      WRITE(9,*) 'Parameters of nuclear realizations:'
      WRITE(9,*) 'GDER parameters: Erez[MeV]   Width[MeV]   Sigma[mb]'
      DO i=1,NGIGE
        WRITE(9,111) ER(i),W0(i),SIG(i)
        ENDDO
  111   FORMAT('               ',2F11.2,F13.1)
      WRITE(9,*) 'GDMR parameters: Erez[MeV]   Width[MeV]   Sigma[mb]'
      DO i=1,NGIGM
        WRITE(9,111) ERM(i),WM0(i),SIGM(i)
      ENDDO
      WRITE(9,*) 'GQER parameters: Erez[MeV]   Width[MeV]   Sigma[mb]'
      DO i=1,NGIGE2
       WRITE(9,111) ERE(i),WE0(i),SIGE(i)
      ENDDO
      WRITE(9,112) DEG,DMG,QEL
  112 FORMAT(' SP strength:   E1:',E9.2,'   M1:',E9.2,'   E2:',E9.2)
      WRITE(9,114) TEMPER,EZERO,AMASS
  114 FORMAT(' Density parameters:  CTF:   Temper:',F5.3,'     E0:',
     &F6.3,'    Amass:',F5.0)
      WRITE(9,115) ASHELL,EONE,PAIRING
  115 FORMAT('                     BSFG:   Ashell:',F7.3,'  E1:',
     &F6.3,'  Pairing:',F6.3)
      WRITE(9,113) fermc,ek0,egzero,par_m1(1),par_m1(2),par_m1(3)
  113 FORMAT(' Fermi liq. par.:',F3.1,'            Ek0:',F4.1,'     Eg0:
     &',F4.1,'    PAR_M1(1):',E9.2,'  PAR_M1(2):',E9.2,
     &'  PAR_M1(3):',E9.2)
      IF (IPINc.EQ.0) THEN
       CHPAR='+'
      ELSE
       CHPAR='-'
      ENDIF
      WRITE(9,117) BN,SPINc(1),CHPAR,CAPFR(1)
  117 FORMAT(' Binding energy:',F7.4,'         Jpi:',F4.1,A1,'  Fract:',
     &F6.3)
      WRITE(9,116) factnrm
  116 FORMAT(' Normalization factor: ',E9.3)
      WRITE(9,*)'*******************************************************
     &***************'
      CLOSE(UNIT=9,STATUS='KEEP')
      END
C***********************************************************************
SUBROUTINE DO_IT (NSTEP)
    C***********************************************************************
    C
          COMMON /monit/ ELQQ(0:126),SPQQ(0:126),IPQQ(0:126),ICQQ(126),
         &               DMQQ(126),WIQQ(0:126)
         &       /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
         &               ISWIC,ISWMX,ISWWI,ISWLS
    
    C
          IF (ISWBN.EQ.0) THEN
    c         WRITE (12,*)   ' '
             WRITE (12,100) SPQQ(0),IPQQ(0),NSTEP,ELQQ(0)
             IF (ISWEL.EQ.1) WRITE (12,200) (ELQQ(J),J=1,NSTEP)
             IF (ISWSP.EQ.1) WRITE (12,300) (SPQQ(J),J=1,NSTEP)
             IF (ISWPA.EQ.1) WRITE (12,400) (IPQQ(J),J=1,NSTEP)
             IF (ISWIC.EQ.1) WRITE (12,400) (ICQQ(J),J=1,NSTEP)
             IF (ISWMX.EQ.1) WRITE (12,500) (DMQQ(J),J=1,NSTEP)
             IF (ISWWI.EQ.1) WRITE (12,500) (WIQQ(J),J=1,NSTEP)
          ELSE
             WRITE (12) SPQQ(0),IPQQ(0),NSTEP,ELQQ(0)
             IF (ISWEL.EQ.1) WRITE (12) (ELQQ(J),J=1,NSTEP)
             IF (ISWSP.EQ.1) WRITE (12) (SPQQ(J),J=1,NSTEP)
             IF (ISWPA.EQ.1) WRITE (12) (IPQQ(J),J=1,NSTEP)
             IF (ISWIC.EQ.1) WRITE (12) (ICQQ(J),J=1,NSTEP)
             IF (ISWMX.EQ.1) WRITE (12) (DMQQ(J),J=1,NSTEP)
             IF (ISWWI.EQ.1) WRITE (12) (WIQQ(J),J=1,NSTEP)
          ENDIF
      100 FORMAT (F5.1,I3,I5,F11.5)
      200 FORMAT (2X,8F11.5)
      300 FORMAT (2X,8F11.1)
      400 FORMAT (2X,8I11)
      500 FORMAT (2X,8E11.3)
          RETURN
          END
    C
    C***********************************************************************
          SUBROUTINE CLOSE_IT
    C***********************************************************************
    C
          CLOSE (UNIT=12,STATUS='KEEP')
          RETURN
          END
    C
C***********************************************************************
SUBROUTINE OPEN_IT (INR,IRL)
    C***********************************************************************
    C
          COMMON /opts/  IVER,IPRIM,ISWWR,ISWBN,ISWEL,ISWSP,ISWPA,
         &               ISWIC,ISWMX,ISWWI,ISWLS
    
    C
          CHARACTER*3  EXTR,EXTS
          IERR=0
          IF (ISWBN.LT.0.OR.ISWBN.GT.1) THEN
             IERR=1
             RETURN
          ENDIF
          WRITE (EXTS,100) INR
          WRITE (EXTR,100) IRL
      100 FORMAT (I3.3)
          IF (ISWBN.EQ.1) THEN
             OPEN (UNIT=12,FILE='EVENTS.S'//EXTS//'.R'//EXTR,
         *        FORM='UNFORMATTED',ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
             ELSE
             OPEN (UNIT=12,FILE='EVENTS.S'//EXTS//'.R'//EXTR,
         *        STATUS='UNKNOWN')
          ENDIF
          RETURN
          END
    C
    