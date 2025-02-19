************************************************************************
*                                                                      *
      SUBROUTINE KEINT (IA,IB,TEGRAL)
*                                                                      *
*   This routine returns kinetic energy integrals.                     *
*      A list of previously computed integrals is maintained. If the   *
*   required  integral is not already  available it is computed  and   *
*   stored in the list.                                                *
*      Note that indices IA and IB may be interchanged on output.      *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, RALLOC, RINTI.                         *
*                                                                      *
*   Written by Farid A Parpia             Last revision: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL FRSTKI,FOUND
*
      POINTER (PINDKE,INDKEI(1))
      POINTER (PVALKE,VALKEI(1))
*
      COMMON/KEILST/NDKEA,NKEI,PINDKE,PVALKE,FRSTKI
*
      include 'parameters.def'
CGG      PARAMETER (KEYORB = 121)
*
*   Ensure that the indices are in `canonical' order
*
      IF (IA .GT. IB) THEN
         ISWAP = IB
         IB = IA
         IA = ISWAP
      ENDIF
*
*   Compute the composite (packed) index
*
      INDEX = IA*KEYORB+IB
*
      IF (.NOT. FRSTKI) THEN
*
*   This branch is executed on all entries except the first
*
         IF (INDEX .GT. INDKEI(NKEI)) THEN
*
*   The index is greater than the largest stored
*
            FOUND = .FALSE.
            LOC = NKEI
*
         ELSEIF (INDEX .LT. INDKEI(1)) THEN
*
*   The index is less than the smallest stored
*
            FOUND = .FALSE.
            LOC = 0
*
         ELSE
*
*   The index is within the range of the indices stored; search
*   for it in the list of indices
*
            JU = NKEI
            JL = 1
*
    1       IF (JU-JL .GT. 1) THEN
*
               JM = (JU+JL)/2
               IF (INDKEI(JM) .GT. INDEX) THEN
                  JU = JM
               ELSE
                  JL = JM
               ENDIF
               GOTO 1
*
            ELSE
*
*   The range is bracketed to the extent possible
*
               IF (INDEX .EQ. INDKEI(JU)) THEN
*
                  FOUND = .TRUE.
                  LOC = JU
*
               ELSEIF (INDEX .EQ. INDKEI(JL)) THEN
*
                  FOUND = .TRUE.
                  LOC = JL
*
               ELSE
*
                  FOUND = .FALSE.
                  LOC = JL
*
               ENDIF
*
            ENDIF
*
         ENDIF
*
         IF (FOUND) THEN
*
*   Found the index in the list; return the value of the integral
*   from storage
*
            TEGRAL = VALKEI(LOC)
*
         ELSE
*
*   Index not found; compute the integral
*
            TEGRAL = RINTI (IA,IB,1)
*
*   Increment the integral counter
*
            NKEI = NKEI+1
*
*   Increase array length by half the present length if the latter
*   is inadequate to store another pair
*
            IF (NKEI .GT. NDKEA) THEN
               NEWSIZ = NDKEA+NDKEA/2
               CALL RALLOC (PINDKE,NDKEA,NEWSIZ,4)
               CALL RALLOC (PVALKE,NDKEA,NEWSIZ,8)
               NDKEA = NEWSIZ
            ENDIF
*
            DO 2 I = NKEI,LOC+2,-1
               INDKEI(I) = INDKEI(I-1)
               VALKEI(I) = VALKEI(I-1)
    2       CONTINUE
*
*   Put the new index and value into storage
*
            INDKEI(LOC+1) = INDEX
            VALKEI(LOC+1) = TEGRAL
*
         ENDIF
*
      ELSE
*
*   This branch is executed only once
*
         FRSTKI = .FALSE.
*
*   Allocate the initial storage for arrays INDOEI and VALOEI;
*   NDCOEA stores the array dimension
*
         NDKEA = 2
         CALL ALLOC (PINDKE,NDKEA,4)
         CALL ALLOC (PVALKE,NDKEA,8)
*
*   Compute the integral's value
*
         TEGRAL = RINTI (IA,IB,1)
*
*   Initialise the integral counter
*
         NKEI = 1
*
*   Store the integral and its value
*
         INDKEI(1) = INDEX
         VALKEI(1) = TEGRAL
*
      ENDIF
*
      RETURN
      END
