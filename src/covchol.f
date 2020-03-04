      SUBROUTINE PRXGRD(N,SIGMA,L,LAMBDA,EPS,ALPHA,MAXITR)
      INTEGER N,MAXITR
      DOUBLE PRECISION SIGMA(N,N),L(N,N),LAMBDA,EPS,ALPHA      
c      
c     internal variables
      INTEGER I,J,ITR
      DOUBLE PRECISION F,FNEW,TMP(N,N),GRD(N,N), D(N),
     *                 ONE, ZERO, STEP, DIFF, G, GNW
      ITR = 0
      ONE = 1.0
      ZERO = 0.0
      DO 20 J = 1,N - 1
         DO 10 I = J + 1,N
            TMP(J,I) = SIGMA(I,J)
            TMP(I,J) = SIGMA(I,J)
  10     CONTINUE          
            TMP(J,J) = SIGMA(J,J)
  20  CONTINUE 
      TMP(N,N) = SIGMA(N,N)
c     compute initial objective function
c     compute TMP = L^(-1) * SIGMA * L^(-t)  
      CALL DTRSM("L","L","N","N",N,N,ONE,L,N,TMP,N) 
      CALL DTRSM("R","L","T","N",N,N,ONE,L,N,TMP,N) 
c     compute tr(TMP) + LAMBDA * ||L||_1,off
      F = 0
      G = 0
      DO 30 J=1,N-1
         F = F + TMP(J,J) + 2 * LOG(ABS(L(J,J))) 
         DO 25 I=J+1, N
            G = G + LAMBDA * ABS(L(I,J))
  25     CONTINUE
  30  CONTINUE
      F = F + TMP(N,N) + 2* LOG(ABS(L(N,N)))
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
c     compute GRD = I - TMP
      DO 35 J=1,N
         DO 32 I=1,N
            GRD(I,J) = -TMP(I,J)
  32     CONTINUE   
            GRD(J,J) = GRD(J,J) + 1 
  35  CONTINUE
c     compute GRD = L**(-t) * GRD 
      CALL DTRSM("L","L","T","N",N,N,2*ONE,L,N,GRD,N) 
c     copy old L before starting line search 
      DO 90 J = 1,N - 1
         DO 80 I = J + 1,N
             L(J,I) = L(I,J)
  80     CONTINUE          
             D(J) = L(J,J)
  90  CONTINUE 
             D(N) = L(N,N)
      STEP = 1
c     line search loop here
  600 CONTINUE     
c     gradient step
      DO 110 J = 1,N-1 
         DO 100 I = J + 1,N
            L(I,J) = L(J,I) - STEP * GRD(I,J) 
  100    CONTINUE
            L(J,J) = D(J) - STEP * GRD(J,J)
            IF (L(J,J) .LT. 0) THEN
                   STEP = STEP * ALPHA
                   GOTO 600 
            ENDIF
  110 CONTINUE
            L(N,N) = D(N) - STEP * GRD(N,N)
            IF (L(N,N) .LT. 0) THEN
                   STEP = STEP * ALPHA
                   GOTO 600 
            ENDIF
c     soft thresholding
      DO 130 J =1,N - 1
         DO 120 I=J + 1,N
            L(I,J) = SIGN(ONE,L(I,J))*(ABS(L(I,J))-STEP*LAMBDA) 
            IF (ABS(L(I,J)) .LE. STEP*LAMBDA) THEN
                     L(I,J) = 0
            ENDIF
            TMP(I,J) = SIGMA(I,J)
            TMP(J,I) = SIGMA(I,J)
 120     CONTINUE
         TMP(J,J) = SIGMA(J,J)
 130  CONTINUE
      TMP(N,N) = SIGMA(N,N)
      CALL DTRSM("L","L","N","N",N,N,ONE,L,N,TMP,N) 
      CALL DTRSM("R","L","T","N",N,N,ONE,L,N,TMP,N) 
c     compute FNW, objective function in new L
      DIFF = 0
      FNW = 0 
      GNW = 0
      DO 140 J=1,N - 1
         FNW = FNW + TMP(J,J) + 2 * LOG(ABS(L(J,J))) 
         DO 135 I=J+1, N
            GNW = GNW + LAMBDA * ABS(L(I,J))
c   new - old
            DIFF = DIFF + ((L(I,J) - L(J,I)) ** 2) / (2 * STEP) +  
     *             (L(I,J) - L(J,I)) * GRD(I,J) 
 135     CONTINUE
            DIFF = DIFF + ((L(J,J) - D(J)) ** 2) / (2 * STEP) + 
     *             (L(J,J) - D(J)) * GRD(J,J) 
 140  CONTINUE
            DIFF = DIFF + ((L(N,N) - D(N)) ** 2) / (2 * STEP) + 
     *             (L(N,N) - D(N)) * GRD(N,N) 
      FNW = FNW + TMP(N,N) + 2 * LOG(ABS(L(N,N)))
c     line search with descent condition
      IF (FNW .GT. (F + DIFF)) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF
c     check stopping criteria
      IF (( ((F+G-FNW-GNW)/ABS(F+G)).LE.EPS).OR.(ITR .GE. MAXITR)) THEN
c     terminate, clean L and save additional outputs
         ALPHA = FNW 
         EPS = (F + G - FNW - GNW) / ABS(F+FNW)   
         MAXITR = ITR
         DO 180 J=1,N-1
            DO 170 I =J + 1,N
               L(J,I) = 0
 170        CONTINUE   
 180     CONTINUE   
         GOTO 900 
      ENDIF  
c     update value of objective function and repeat
      F = FNW
      G = GNW
      GOTO 500
 900  CONTINUE
      RETURN
c     last line of PRXGRD
      END

