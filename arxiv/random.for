c***********************************************************************
SUBROUTINE RANDINIT()
    c***********************************************************************
    c     Matrix of 4xNNR random seeds is generated 
          INTEGER*4  IRINIT,IDUMMY
          REAL       dummy
    
          COMMON /rnd/ IRINIT(4,999)
         &       /readi/ nevents,NSREAL,NRL,ecrit,xrayk,xrayl,numlev,
         &               elowlev(99),elowsp(99),ilowip(99),
         &               STDISa(0:2,0:30,-2:2,0:1),RADW
    
          DO I = 1, 4
            IDUMMY = IRINIT(I,1)
            DO J = 2, NSREAL
              DO K = 1, 100  
                dummy = RAN0(IDUMMY)
              ENDDO
              IRINIT(I,J) = IDUMMY - J
              IF (IRINIT(I,J).LE.0) IRINIT(I,J) = IRINIT(I,J) + 2 * J
            ENDDO
          ENDDO
          RETURN
          END SUBROUTINE
    !***********************************************************************
C***********************************************************************
FUNCTION SEEDS(MODE,SP,IP,IL,IB)
    C
    C     This function yields the seed that is needed for the intrinsic
    C     function RAN inside the subroutine WIDTHS immediately before
    C     this subroutine starts to generate the values of the subscripted
    C     variable STCON.
    c
    C***********************************************************************
    C
          INTEGER*4    MAXIS
          PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
         *              MAXBIN = 700,  !maximum # of bins in continuum
         *              MAXIS  = 76000000        )   
          INTEGER*4      SEEDS,IRCON,IRCONc
          integer        dekod,denum,delev,deparity
          COMMON /WID/   sall(99,0:20),ponv(99,0:20),ponvk(99,0:20),
         &               LEVCON(0:MAXBIN,0:MAXJC,0:1),IRCON(MAXIS),
         &               ENDIS(30,0:MAXJC,0:1),NDIS(0:MAXJC,0:1),
         &               levdis(0:MAXJC,0:1),nddd, dekod(30,0:MAXJC,0:1),
         &               denum(99),delev(99,20),despin(99,20),
         &               deparity(99,20),IRCONc(2),deltx(99,20)
         *      /MAINE/  eall,EIN,EFI,NOPTFL,BN,DELTA,NBIN,ntotal,iregi,
         *               noptcs,nlinc,capfr(2)
    C
          IF (MODE.GT.0) THEN
            SEEDS=IRCONc(MODE)
          ELSE
            SEEDS=IRCON(LEVCON(IB-1,ISUBSC(SP),IP)+IL)
          ENDIF
    C
          RETURN
          END
    C
    C*********************************************************************
          double precision function dran(ir)
    c
    c   this subroutine generates random values between 0.0 and 1.0 using
    c   an integer seed
    c   it is based on the imsl routine ggubs.
    c   double precision version
    c
          implicit double precision (a-h,o-z)
          integer*4 ir
          parameter(da=16807.d0,db=2147483647.d0,dc=2147483648.d0)
          ir=abs(mod(da*ir,db)+0.5d0)
          dran=dfloat(ir)/dc
          return
          end
    c
    c
          function sran(ir)
    c
    c   this subroutine generates random values between 0.0 and 1.0 using a 32-bit integer seed
    c   it is based on the imsl routine ggubs.
    c
          implicit double precision (a-h,o-z)
          integer*4 ir
          real*4 sran
          parameter(da=16807.d0,db=2147483647.d0,dc=2147483648.d0)
          ir=abs(mod(da*ir,db)+0.5d0)
          sran=float(ir)/dc
          return
          end
    c
    c
    c   "Minimal" RND generator of Park and Miller - from Numerical Recipes
    c
          function ran0(ir)
          integer*4 ir,IA,IM,IQ,IS,mask,K
    c      real*8    AM
          parameter (IA=16807,IM=2147483647,IQ=127773,IS=2836,
         *           mask=123459876)
    c
    c  other posible values of constants are:
    c
    c      parameter (IA=48271,IM=2147483647,IQ=44488,IS=3399,
    c     *           mask=123459876)
    c      parameter (IA=69621,IM=2147483647,IQ=30845,IS=23902,
    c     *           mask=123459876)
    c
    c      ir=ieor(ir,mask)
       1  continue                   ! line added by MK
          AM=1./float(IM)
          K=ir/IQ
          ir=IA*(ir-K*IQ)-IS*K
          if (ir.LT.0) ir=ir+IM
    c      ran0=AM*float(ir)
          ran0=AM*dble(ir)
    c      ir=ieor(ir,mask)
           if (ran0.GE.1.0) goto 1   ! line added by MK
          return
          end
    c
    c  *********************************************************
    c
    c   "Minimal" RND generator of Park and Miller with Bays-Durham
    c   shuffle and added safeguards - from Numerical Recipes
    c
          function ran1(ir)
          integer*4 ir,IA,IM,IQ,IS,NTAB,NDIV,K
          parameter (IA=16807,IM=2147483647,IQ=127773,IS=2836,
         *           NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    c
    c  other posible values of constants are:
    c
    c      parameter (IA=48271,IM=2147483647,IQ=44488,IS=3399)
    c      parameter (IA=69621,IM=2147483647,IQ=30845,IS=23902)
    c
          integer iv(NTAB)
          save iv,iy
          data iv /NTAB*0/, iy /0/
    c
          am=1./float(IM)
          if (ir.LE.0.or.iy.EQ.0) then
              ir=max(-ir,1)
              do j=NTAB+8,1,-1
                 K=ir/IQ
                 ir=IA*(ir-K*IQ)-IS*K
                 if (ir.LT.0) ir=ir+IM
                 if (j.LE.NTAB) iv(j)=ir
              enddo
              iy=iv(1)
          endif
          K=ir/IQ
          ir=IA*(ir-K*IQ)-IS*K
          if (ir.LT.0) ir=ir+IM
          j=1+iy/NDIV
          iy=iv(j)
          iv(j)=ir
          ran1=min(AM*iy,RNMX)
          return
          end
    c  ********************************************************
    c
    c   Congruental random number generator
    c
            function ran2(ir)
            integer*4 ir,IA,IC,IM
          parameter (IA=7141,IM=259200,IC=54773)
    c      parameter (IA=8121,IM=134456,IC=28411)
    
        1 continue
          ir=mod(ir*IA+IC,IM)
          if (ir.LT.0) goto 1
          ran2=float(ir)/float(IM)
          end
    
          C***********************************************************************
          FUNCTION NPOISS(IR,AM)
    C
    C     Poisson's random numbers
    c
    C***********************************************************************
    C
          INTEGER*4 IR
          COMMON /GAU/ U,IFLAG
          IF (AM.GT.50.) GOTO 1
          Q=1.
          NPOISS=0
        2 Q=Q*RAN0(IR)
          IF (Q.LT.EXP(-AM)) RETURN
          NPOISS=NPOISS+1
          GOTO 2
        1 NPOISS=INT(AM+SQRT(AM)*GAUSS(IR))
          RETURN
          END
    C
    C
    C
    C***********************************************************************
          FUNCTION GAUSS(IR)
    C
    C     Normally distributed random numbers with
    C     a zero mean and a unit variance
    c
    C***********************************************************************
    C
          INTEGER*4 IR
          COMMON /GAU/ U,IFLAG
          IF(IFLAG) 2,1,2
        1 X1=2.*RAN0(IR)-1.
          X2=2.*RAN0(IR)-1.
          S=X1*X1+X2*X2
          IF (S.GE.1.) GOTO 1
          S=SQRT(-2.*ALOG(S)/S)
          U=X1*S
          IFLAG=1
          GAUSS=X2*S
          RETURN
        2 IFLAG=0
          GAUSS=U
          RETURN
          END
    C