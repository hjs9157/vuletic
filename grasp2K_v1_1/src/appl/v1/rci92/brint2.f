************************************************************************
*                                                                      *
      SUBROUTINE BRINT2 (IA,IB,IC,ID,NU,TEGRAL)
*                                                                      *
*   Returns integrals for the transverse photon interaction.           *
*      Integrals are stored in ordered lists. If the integral cannot   *
*   be read from a list, it is computed by calling BRINTF.             *
*                                                                      *
*   Observ that it is not possible to use more than 100 orbitals when  *
*   the Breit interaction is included. If 100 orbitals are used then   *
*   NU <= 19.                                                          *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
*               [RCI92]: BRINTF,                                       *
*                                                                      *
*   Written by Farid A Parpia               Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      LOGICAL FIRST,FOUND
*
Cww   INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)

      POINTER (PINDT1,INDTP1(1))
      POINTER (PVALT1,VALTP1(1))
      POINTER (PINDT2,INDTP2(1))
      POINTER (PVALT2,VALTP2(1))
      POINTER (PINDT3,INDTP3(1))
      POINTER (PVALT3,VALTP3(1))
      POINTER (PINDT4,INDTP4(1))
      POINTER (PVALT4,VALTP4(1))
      POINTER (PINDT5,INDTP5(1))
      POINTER (PVALT5,VALTP5(1))
      POINTER (PINDT6,INDTP6(1))
      POINTER (PVALT6,VALTP6(1))
*
      COMMON/BILST/PINDT1,PINDT2,PINDT3,PINDT4,PINDT5,PINDT6,
     :             PVALT1,PVALT2,PVALT3,PVALT4,PVALT5,PVALT6,
     :             NDTPA(6),NTPI(6),FIRST(6)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /BALLOC/NB(2)
*
      KEY = NW + 1
*
*   Compute the integral label
*
      INDEX = (((NU*KEY+IA)*KEY+IB)*KEY+IC)*KEY+ID
*
      IF (.NOT. FIRST(2)) THEN
*
*   This branch is executed on all entries except the first
*
         IF (INDEX .GT. INDTP2(NTPI(2))) THEN
*
*   The index is greater than the largest stored
*
            FOUND = .FALSE.
            LOC = NTPI(2)
*
         ELSEIF (INDEX .LT. INDTP2(1)) THEN
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
            JU = NTPI(2)
            JL = 1
    1       IF (JU-JL .GT. 1) THEN
               JM = (JU+JL)/2
               IF (INDTP2(JM) .GT. INDEX) THEN
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
               IF (INDEX .EQ. INDTP2(JU)) THEN
*
                  FOUND = .TRUE.
                  LOC = JU
*
               ELSEIF (INDEX .EQ. INDTP2(JL)) THEN
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
            TEGRAL = VALTP2(LOC)
*
         ELSE
*
*   Index not found; compute the integral
*
            TEGRAL = BRINTF (2,IA,IB,IC,ID,NU)
*
*   Increment the integral counter
*
            NTPI(2) = NTPI(2)+1
*
*   Increase array length by half the present length if the latter
*   is inadequate to store another pair
*
*
            IF (NTPI(2) .GT. NDTPA(2)) THEN
               NEWSIZ = NDTPA(2)+NDTPA(1)/2
               CALL RALLOC (PINDT2,NDTPA(2),NEWSIZ,4)
               CALL RALLOC (PVALT2,NDTPA(2),NEWSIZ,8)
               NDTPA(2)= NEWSIZ
            ENDIF
            DO 3 I = NTPI(2),LOC+2,-1
               INDTP2(I) = INDTP2(I-1)
               VALTP2(I) = VALTP2(I-1)
    3       CONTINUE
*
*   Put the new index and value into storage
*
            INDTP2(LOC+1) = INDEX
            VALTP2(LOC+1) = TEGRAL
*
         ENDIF
*
      ELSE
*
*   This branch is executed only once per type of integral
*
         FIRST(2) = .FALSE.
*
*   Designate the initial storage for arrays INDTPx and VALTPx;
*   Array NDTPA stores the array dimensions
*
C         NDTPA(2) = NB(2)
         NDTPA(2) = 100000
*
*   Compute the integral's value
*
         TEGRAL = BRINTF (2,IA,IB,IC,ID,NU)
*
*   Initialise the integral counter
*
         NTPI(2) = 1
*
*   Store the integral and its value
*
         CALL ALLOC (PINDT2,NDTPA(2),4)
         CALL ALLOC (PVALT2,NDTPA(2),8)
         INDTP2(1) = INDEX
         VALTP2(1) = TEGRAL
*
      ENDIF
*
      RETURN
      END
