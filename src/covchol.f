      SUBROUTINE PRXGRDL(N,SIGMA,L,LAMBDA,EPS,ALPHA,MAXITR)
      INTEGER N,MAXITR,JOB
      DOUBLE PRECISION SIGMA(N,N),L(N,N),LAMBDA,EPS,ALPHA      
c      
c     internal variables
      INTEGER I,J,ITR
      DOUBLE PRECISION F,FNEW,TMP(N,N),
     *                 ONE, ZERO, STEP 
      ITR = 0
      ONE = 1.0
      ZERO = 0.0
      DO 20 J = 1,N - 1
         DO 10 I = J + 1,N
            TMP(J,I) = SIGMA(I,j)
            TMP(I,J) = SIGMA(I,j)
            L(I,J) = L(I,J) / L(J,J) 
  10     CONTINUE          
            TMP(J,J) = SIGMA(J,J)
  20  CONTINUE 
      TMP(N,N) = SIGMA(N,N)
c     compute initial objective function
c     compute TMP = L^(-1) * SIGMA * L^(-t)  
      CALL DTRSM("L","L","N","U",N,N,ONE,L,N,TMP,N) 
      CALL DTRSM("R","L","T","U",N,N,ONE,L,N,TMP,N) 
c     compute tr(D**(-1) * TMP)
      F = 0
      DO 30 J=1,N-1
         F = F + TMP(J,J) / (L(J,J)**2) 
         DO 25 I=J+1, N
            F = F + LAMBDA * ABS(L(I,J))
  25     CONTINUE
  30  CONTINUE
      F = F + TMP(N,N) / (L(N,N)**2) 
c     main loop here, increase iteration counter
 500  CONTINUE      
      ITR = ITR + 1
c     obtain half gradient
c     TMP = L**(-1) * SIGMA * L**(-t) * D**-1 * L**(-1)
      DO 60 J=1,N
         DO 50 I=1,N
            TMP(I,J) = TMP(I,J)  
  50     CONTINUE
  60  CONTINUE       
      CALL DTRSM("R","L","N","U",N,N,ONE,L,N,TMP,N) 
c     copy old L before starting line search 
      DO 90 J = 1,N - 1
         DO 80 I = J + 1,N
            L(J,I) = L(I,J)
  80     CONTINUE          
  90  CONTINUE 
      STEP = 1
c     line search loop here
  600 CONTINUE     
c     gradient step
      DO 110 J = 1,N - 1
         DO 100 I = J + 1,N
            L(I,J) = L(J,I) + 2 * STEP * TMP(I,J) 
  100    CONTINUE
  110 CONTINUE
c     soft thresholding
      DO 130 J =1,N - 1
         DO 120 I=J + 1,N
            L(I,J) = SIGN(ONE,L(I,J))*(ABS(L(I,J))-STEP*LAMBDA) 
            IF (ABS(L(I,J)) .LT. STEP*LAMBDA) THEN
                    L(I,J) = 0
            ENDIF
            TMP(I,J) = SIGMA(I,J)
            TMP(J,I) = SIGMA(I,J)
 120     CONTINUE
         TMP(J,J) = SIGMA(J,J)
 130  CONTINUE
      TMP(N,N) = SIGMA(N,N)
      CALL DTRSM("L","L","N","U",N,N,ONE,L,N,TMP,N) 
      CALL DTRSM("R","L","T","U",N,N,ONE,L,N,TMP,N) 
c     compute FNW, objective function in new L
      FNW = 0 
      DO 140 J=1,N - 1
         FNW = FNW + TMP(J,J) / (L(J,J)**2) 
         DO 135 I=J+1, N
            FNW = FNW + LAMBDA * ABS(L(I,J))
 135     CONTINUE
 140  CONTINUE
      FNW = FNW + TMP(N,N) / (L(N,N)**2) 
c     line search with descent condition
      IF (FNW .GT. F) THEN
         STEP = STEP * ALPHA
         GOTO 600
      ENDIF
c     check stopping criteria
      IF (((F - FNW) / ABS(F) .LE. EPS) .OR.  (ITR .GE. MAXITR)) THEN
c     terminate, clean L and save additional outputs
         ALPHA = FNW 
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
c     last line of PRXGRDLLB
      END

