C***********************************************************************
FUNCTION TERM(EEXC)
    c
    c     The effective temperature corresponding to the exc.energy EEXC
    c
    C***********************************************************************
    C attention!! : PAIRING (not EONE) is used for the shift 
    c
          COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
         *               DENLO,DENHI,DENPAR(4),
         *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
    c
          EEFF=EEXC-PAIRING
          IF (EEFF.LT.0.) EEFF=0.
          TERM=SQRT((EEFF)/ASHELL)
          RETURN
          END
    C
    C***********************************************************************
          FUNCTION TERM1(EEXC)
    c
    c     The effective temperature corresponding to the exc.energy EEXC
    c
    C***********************************************************************
    c
          COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
         *               DENLO,DENHI,DENPAR(4),
         *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
    c
          EEFF=EEXC
          IF (EEFF.LT.0.) EEFF=0.
          TERM1=SQRT((EEFF)/ASHELL)
          RETURN
          END
    C
    C***********************************************************************
          FUNCTION TERMOSLO(EEXC)
    c
    c     The effective temperature corresponding to the exc.energy EEXC
    c
    C***********************************************************************
    C attention!! : PAIRING (not EONE) used for shift
    c
          COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
         *               DENLO,DENHI,DENPAR(4),
         *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
    c
          EEFF=EEXC-PAIRING+6.6/(AMASS**0.32)
          IF (EEFF.LT.0.) EEFF=0.
          TERMOSLO=SQRT((EEFF)/ASHELL)
          RETURN
          END
    C
    C***********************************************************************
          FUNCTION EFEXC(EEXC)
    c
    c     The effective temperature corresponding to the exc.energy EEXC
    c
    C***********************************************************************
    C attention!! : PAIRING (not EONE) used for shift
    c
          COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
         *               DENLO,DENHI,DENPAR(4),
         *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
    c
          EFEXC=EEXC-PAIRING
          RETURN
          END
    C
    C***********************************************************************
          FUNCTION TERMDILG(EEXC)
    c
    c     The effective temperature corresponding to the exc.energy EEXC
    c
    C***********************************************************************
    C attention!! : PAIRING (not EONE) used for shift
    c
          COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
         *               DENLO,DENHI,DENPAR(4),
         *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
    c
          EEFF=EEXC-PAIRING
          IF (EEFF.LT.0.) EEFF=0.
          TERMDILG=(1.0 + SQRT(1.0+4.0*ASHELL*EEFF))/2.0/ASHELL
          RETURN
          END
    C
    C***********************************************************************
          FUNCTION SMFREQ()
    c
          COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
         *               DENLO,DENHI,DENPAR(4),
         *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
    c
          SMFREQ=41./AMASS**(0.33333)
          RETURN
          END
    C
    C***********************************************************************
          FUNCTION EREL1()
    c
          COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
         *               DENLO,DENHI,DENPAR(4),
         *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
    c
          AMASS13=AMASS**(0.33333)
          EREL1=31.2/AMASS13+20.6/sqrt(AMASS13)
          RETURN
          END
    C
    C***********************************************************************
          FUNCTION WREL1()
    c
          WREL1=0.026*EREL1()**1.91
          RETURN
          END
    C
    C***********************************************************************
          FUNCTION R0()
    c
          COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
         *               DENLO,DENHI,DENPAR(4),
         *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
    c
          R0=1.27*AMASS**(0.33333)
          RETURN
          END
    C
    C***********************************************************************
          FUNCTION FKS(EX,FNS,FKs0,WWALL,ALPPL)
    c
          EREL=EREL1()
          if (EX.GT.(2*EREL)) then
            FKS=FKs0     
          else
            FACT=WREL1()-EREL**2/ALPPL
            FKS=FACT/WWALL+(FKs0-FACT/WWALL)*ABS((EX-EREL)/EREL)**FNS
          endif
          RETURN
          END
    C
    C***********************************************************************
          FUNCTION WWALL1()
    c
          COMMON /DEN/   NOPTDE,EZERO,EONE,TEMPER,ASHELL,AMASS,PAIRING,FJ,
         *               DENLO,DENHI,DENPAR(4),
         *               ASHELL09,EONE09,TEMPER09,EZERO09,PAIRING09
    c
              WWALL1=36.43*AMASS**(-0.3333333) 
          RETURN
          END
    C
    C***********************************************************************

    C*******************************************************************************
C
      subroutine tred1(nm,n,a,d,e,e2)
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h,scale
c
c     this subroutine is a translation of the algol procedure tred1,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix
c     to a symmetric tridiagonal matrix using
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        a contains information about the orthogonal trans-
c          formations used in the reduction in its strict lower
c          triangle.  the full upper triangle of a is unaltered.
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      do 100 i = 1, n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
  100 continue
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
c
         do 125 j = 1, l
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = 0.0d0
  125    continue
c
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 300
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         e2(i) = scale * scale * h
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
         if (l .eq. 1) go to 285
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            g = e(j) + a(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + a(k,j) * d(k)
               e(k) = e(k) + a(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         h = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - h * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
c
  280    continue
c
  285    do 290 j = 1, l
            f = d(j)
            d(j) = a(l,j)
            a(l,j) = a(i,j)
            a(i,j) = f * scale
  290    continue
c
  300 continue
c
      return
      end
C
C*******************************************************************************
C
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing floating point
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to
c            the accuracy used in floating point variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger floating point number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end
C
C*******************************************************************************
C
      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
C
C*******************************************************************************
C
SUBROUTINE TQLRAT(N,D,E2,IERR)
    C
          INTEGER I,J,L,M,N,II,L1,MML,IERR
          DOUBLE PRECISION D(N),E2(N)
          DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG
    C
    C     This subroutine is a translation of the Algol procedure tqlrat,
    C     Algorithm 464, Comm. ACM 16, 689(1973) by Reinsch.
    C
    C     This subroutine finds the eigenvalues of a symmetric
    C     tridiagonal matrix by the rational QL method.
    C
    C     On input
    C
    C        N is the order of the matrix.
    C
    C        D contains the diagonal elements of the input matrix.
    C
    C        E2 contains the squares of the subdiagonal elements of the
    C          input matrix in its last N-1 positions.  E2(1) is arbitrary.
    C
    C      On output
    C
    C        D contains the eigenvalues in ascending order.  If an
    C          error exit is made, the eigenvalues are correct and
    C          ordered for indices 1,2,...IERR-1, but may not be
    C          the smallest eigenvalues.
    C
    C        E2 has been destroyed.
    C
    C        IERR is set to
    C          zero       for normal return,
    C          J          if the J-th eigenvalue has not been
    C                     determined after 30 iterations.
    C
    C     Calls PYTHAG for  DSQRT(A*A + B*B) .
    C
    C     Questions and comments should be directed to Burton S. Garbow,
    C     Mathematics and Computer Science Div, Argonne National Laboratory
    C
    C     This version dated August 1987.
    C     Modified by C. Moler to fix underflow/overflow difficulties,
    C     especially on the VAX and other machines where epslon(1.0d0)**2
    C     nearly underflows.  See the loop involving statement 102 and
    C     the two statements just before statement 200.
    C
    C     ------------------------------------------------------------------
    C
    
          IERR = 0
          IF (N .EQ. 1) GO TO 1001
    C
          DO 100 I = 2, N
      100 E2(I-1) = E2(I)
    C
          F = 0.0D0
          T = 0.0D0
          E2(N) = 0.0D0
    C
          DO 290 L = 1, N
             J = 0
             H = DABS(D(L)) + DSQRT(E2(L))
             IF (T .GT. H) GO TO 105
             T = H
             B = EPSLON(T)
             C = B * B
             if (c .ne. 0.0d0) go to 105
    C        Spliting tolerance underflowed.  Look for larger value.
             do 102 i = l, n
                h = dabs(d(i)) + dsqrt(e2(i))
                if (h .gt. t) t = h
      102    continue
             b = epslon(t)
             c = b * b
    C     .......... LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ..........
      105    DO 110 M = L, N
                IF (E2(M) .LE. C) GO TO 120
    C     .......... E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
    C                THROUGH THE BOTTOM OF THE LOOP ..........
      110    CONTINUE
    C
      120    IF (M .EQ. L) GO TO 210
      130    IF (J .EQ. 30) GO TO 1000
             J = J + 1
    C     .......... FORM SHIFT ..........
             L1 = L + 1
             S = DSQRT(E2(L))
             G = D(L)
             P = (D(L1) - G) / (2.0D0 * S)
             R = PYTHAG(P,1.0D0)
             D(L) = S / (P + DSIGN(R,P))
             H = G - D(L)
    C
             DO 140 I = L1, N
      140    D(I) = D(I) - H
    C
             F = F + H
    C     .......... RATIONAL QL TRANSFORMATION ..........
             G = D(M)
             IF (G .EQ. 0.0D0) G = B
             H = G
             S = 0.0D0
             MML = M - L
    C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
             DO 200 II = 1, MML
                I = M - II
                P = G * H
                R = P + E2(I)
                E2(I+1) = S * R
                S = E2(I) / R
                D(I+1) = H + S * (H + D(I))
                G = D(I) - E2(I) / G
    C           Avoid division by zero on next pass
                if (g .eq. 0.0d0) g = epslon(d(i))
                h = g * (p / r)
      200    CONTINUE
    C
             E2(L) = S * G
             D(L) = H
    C     .......... GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ..........
             IF (H .EQ. 0.0D0) GO TO 210
             IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
             E2(L) = H * E2(L)
             IF (E2(L) .NE. 0.0D0) GO TO 130
      210    P = D(L) + F
    C     .......... ORDER EIGENVALUES ..........
             IF (L .EQ. 1) GO TO 250
    C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
             DO 230 II = 2, L
                I = L + 2 - II
                IF (P .GE. D(I-1)) GO TO 270
                D(I) = D(I-1)
      230    CONTINUE
    C
      250    I = 1
      270    D(I) = P
      290 CONTINUE
    C
          GO TO 1001
    C     .......... SET ERROR -- NO CONVERGENCE TO AN
    C                EIGENVALUE AFTER 30 ITERATIONS ..........
     1000 IERR = L
     1001 RETURN
          END
    C

    C*******************************************************************************
C
      SUBROUTINE SORT (N,X,XN)
C
C        ALGORITHM AS 304.8 APPL.STATIST. (1996), VOL.45, NO.3
C
C        Sorts the N values stored in array X in ascending order
C
      INTEGER N
      DOUBLE PRECISION X(N),XN(N),XTEMP(N)
C
      INTEGER I, J, INCR
      DOUBLE PRECISION TEMP
C
	DO I = 1, N
       XTEMP(I) = X(I) 
	ENDDO
C
      INCR = 1
C
C        Loop : calculate the increment
C
   10 INCR = 3 * INCR + 1
      IF (INCR .LE. N) GOTO 10

C
C        Loop : Shell-Metzner sort
C
   20 INCR = INCR / 3
      I = INCR + 1
   30 IF (I .GT. N) GOTO 60
      TEMP = X(I)
      J = I
   40 IF (X(J - INCR) .LT. TEMP) GOTO 50
      X(J) = X(J - INCR)
      J = J - INCR
      IF (J .GT. INCR) GOTO 40
   50 X(J) = TEMP
      I = I + 1
      GOTO 30
   60 IF (INCR .GT. 1) GOTO 20
C
	DO I = 1, N
        XN(I) = X(I)
        X(I) = XTEMP(I) 
	ENDDO
C
      RETURN
      END
C
C*******************************************************************************
C
      SUBROUTINE rsm1(nm,n,a,w)
c
      integer n,nm,ierr
      integer k1,k2
c      double precision a(nm,n),w(n),fwork(1)
      double precision a(nm,n),w(n),fwork2(n),fwork1(n)
c
c     this subroutine calls the recommended sequence of
c     subroutines from the eigensystem subroutine package (eispack)
c     to find all of the eigenvalues and some of the eigenvectors
c     of a real symmetric matrix.
c
c     on input
c
c        nm  must be set to the row dimension of the two-dimensional
c        array parameters as declared in the calling program
c        dimension statement.
c
c        n  is the order of the matrix  a.
c
c        a  contains the real symmetric matrix.
c
c        m  the eigenvectors corresponding to the first m eigenvalues
c           are to be computed.
c           if m = 0 then no eigenvectors are computed.
c           if m = n then all of the eigenvectors are computed.
c
c     on output
c
c        w  contains all n eigenvalues in ascending order.
c
c        z  contains the orthonormal eigenvectors associated with
c           the first m eigenvalues.
c
c        ierr  is an integer output variable set equal to an error
c           completion code described in the documentation for tqlrat,
c           imtqlv and tinvit.  the normal completion code is zero.
c
c        fwork  is a temporary storage array of dimension 8*n.
c
c        iwork  is an integer temporary storage array of dimension n.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 10 * n
      if (n .gt. nm) go to 50
      k1 = 1
      k2 = k1 + n

c     .......... find eigenvalues only ..........

c      call  tred1(nm,n,a,w,fwork(k1),fwork(k2))
c      call  tqlrat(n,w,fwork(k2),ierr)

      call  tred1(nm,n,a,w,fwork1,fwork2)
      call  tqlrat(n,w,fwork2,ierr)

   50 return
      end
C
C***********************************************************************
FUNCTION GOE_EIGEN_VAL(N,IS,IP)
    C***********************************************************************
    C
    C     Eigenvalues were generated at the beginning 
    C     by SUBROUTINE GENERATE_GOE_EIGEN_VAL(IR,N)
    C
          INTEGER*4    MAXIS
          PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
         *              MAXBIN = 700,  !maximum # of bins in continuum
         *              MAXIS  = 76000000        )   
          INTEGER*4 N
          COMMON /GOE/   EIGENVAL(0:MAXBIN,0:MAXJC,0:1),
         *               NEIGENVAL(0:MAXJC,0:1),LMODE
    C
          NN=NEIGENVAL(IS,IP)
          K=N-(N/NN)*NN
          GOE_EIGEN_VAL=EIGENVAL(K,IS,IP)+FLOAT((N/NN)*NN)
          RETURN
          END
    C
    C***********************************************************************
          SUBROUTINE GENERATE_GOE_EIGEN_VAL(IR,N)
    C***********************************************************************
    C
    C     Generates eigenvalues of random matrices - they are stored in
    C     the COMMON /GOE/ and are used for generating level in the 
    C     subroutine LEVELS via calling GOE_EIGEN_VAL
    C
          INTEGER*4    MAXIS
          PARAMETER    (MAXJC  = 49,   !maximum spin generated in level scheme
         *              MAXBIN = 700,  !maximum # of bins in continuum
         *              MAXIS  = 76000000        )   
          INTEGER   N,LDA
          PARAMETER (LDA=1001)
          REAL*8    A(LDA,N),RA(N),X,RAORD(N)
          INTEGER*4 IR
    C
          COMMON /GOE/   EIGENVAL(0:MAXBIN,0:MAXJC,0:1),
         *               NEIGENVAL(0:MAXJC,0:1),LMODE
          COMMON  /GAU/ U,IFLAG
    C
          DATA      DPI  /3.141592653589793d0/ 
    C
          IFLAG=0
          DO IP=0,1
           DO IS=0,MAXJC
            DO I=1,N
             A(I,I)=DBLE(GAUSS(IR)/SQRT(2.*FLOAT(N-1)))
             DO J=I+1,N
              A(I,J)=DBLE(GAUSS(IR)/SQRT(4.*FLOAT(N-1)))
              A(J,I)=A(I,J)
             ENDDO !J
            ENDDO !I
    c        CALL DEVLSF(N,A,LDA,RA)
            CALL RSM1(LDA,N,A,RA)     
    c        CALL SORTQQ(LOC(RA),N,SRT$REAL8)
            CALL SORT(N,RA,RAORD)
            DO I = 1,N
              RA(I) = RAORD(I)
            ENDDO
    
            M=-1
            DO L=1,N
             X=RA(L)	      
             IF (DABS(X).lt.0.9d0) THEN         !Only eigenval. in the middle are treated
              M=M+1
              EIGENVAL(M,IS,IP)=FLOAT(N-1)*
         *            (SNGL((X*DSQRT(1.d0-X**2)+DASIN(X))/DPI+.5d0))	    	      
             ENDIF		 			
            ENDDO !L
            NEIGENVAL(IS,IP)=M
            EIG0=EIGENVAL(0,IS,IP)
            DO MM=M,0,-1
             EIGENVAL(MM,IS,IP)=EIGENVAL(MM,IS,IP)-EIG0
            ENDDO !MM
           ENDDO !IS
          ENDDO !IP
          RETURN
          END
    
    C
    
C
c This file contain basic functions used in various programs concerning
c DICEBOX
c Namely: ISUBSC (SPIN)
c         ITYPE (SPIN,IPIN,SPFI,IPFI)
c         NPOISS (IR,AM)
c         GAUSS (IR)
c
C***********************************************************************
FUNCTION ISUBSC(A)
    C
    C     The function ISUBSC is needed in such situations
    C     when values of spin are to be transformed to values
    C     that can serve as subscripts. The function converts
    C     values of a real variable A to integer values
    C     according to the following rule:
    C
    C                A = 0.0  ISUBSC = 0
    C                    0.5           0
    C                    1.0           1
    C                    1.5           1
    C                    2.0           2
    C                    2.5           2
    C                     .            .
    C                     .            .
    C                    etc.         etc.
    C
    C     This function is to be immune against a finite precission in
    C     manipulations with real values.
    C
    C***********************************************************************
    c
          I=INT(A)
          X=A-FLOAT(I)
          IF (ABS(X-.5)-0.25) 1,2,2
    C
    C     The case of half-integer spin:
    C
        1 ISUBSC=I
          RETURN
    C
    C     The case of integer spin:
    C
        2 ISUBSC=NINT(A)
          RETURN
          END
    C
    C
    C
