************************************************************************
*                                                                      *
      SUBROUTINE MAXARR (J)
*                                                                      *
*   This subroutine finds the least self-consistent orbital            *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 10 Dec 1992   *
*                                                                      *
************************************************************************
*
! J initialized to zero
! XHH 1997.02.14
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LFIX
*
      COMMON/FIXD/NFIX,LFIX(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
*
      J = 0
      DLRGST = 0.0D 00
      DO 1 I = 1,NW
         IF (.NOT. LFIX(I)) THEN
            IF (SCNSTY(I) .GT. DLRGST) THEN
               DLRGST = SCNSTY(I)
               J = I
            ENDIF
         ENDIF
    1 CONTINUE
*
      RETURN
      END
