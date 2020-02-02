      SUBROUTINE PRXGRD(N,SIGMA,L,LAMBDA,EPS,ALPHA,MAXITR)
      INTEGER N,MAXITR
      DOUBLE PRECISION SIGMA(N,N),L(N,N),LAMBDA,EPS,ALPHA      
c      
c     internal variables
      INTEGER I,J,ITR
      DOUBLE PRECISION F,FNEW,TMP(N,N),GRD(N,N), D(N),
     *                 ONE, ZERO, STEP 
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
      DO 30 J=1,N-1
         F = F + TMP(J,J) + 2 * LOG(L(J,J)) 
         DO 25 I=J+1, N
            F = F + LAMBDA * ABS(L(I,J))
  25     CONTINUE
  30  CONTINUE
      F = F + TMP(N,N) + 2* LOG(L(N,N))
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
c     obtain gradient / 2
      DO 35 J=1,N
         DO 32 I=1,N
            GRD(I,J) = TMP(I,J)
  32     CONTINUE   
  35  CONTINUE
      CALL DTRSM("R","L","N","N",N,N,ONE,L,N,GRD,N) 
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
            L(I,J) = L(J,I) + 2 * STEP * GRD(I,J) 
  100    CONTINUE
            L(J,J) = D(J) + 2 * STEP * GRD(J,J)
  110 CONTINUE
            L(N,N) = D(N) + 2 * STEP * GRD(N,N)
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
      FNW = 0 
      DO 140 J=1,N - 1
         FNW = FNW + TMP(J,J) + 2 * LOG(L(J,J)) 
         DO 135 I=J+1, N
            FNW = FNW + LAMBDA * ABS(L(I,J))
 135     CONTINUE
 140  CONTINUE
      FNW = FNW + TMP(N,N) + 2 * LOG(L(N,N))
c     line search with descent condition
      IF (FNW .GT. F) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF
c     check stopping criteria
      IF (((F - FNW) / ABS(F) .LE. EPS) .OR.  (ITR .GE. MAXITR)) THEN
c     terminate, clean L and save additional outputs
         ALPHA = STEP 
         EPS = (F - FNW) / ABS(F)   
         MAXITR = ITR
         DO 160 J=1,N-1
            DO 150 I =J + 1,N
               L(J,I) = 0
 150        CONTINUE   
 160     CONTINUE   
         GOTO 900 
      ENDIF  
c     update value of objective function and repeat
      F = FNW
      GOTO 500
 900  CONTINUE
      RETURN
c     last line of PRXGRD
      END

