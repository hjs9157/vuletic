************************************************************************
*                                                                      *
      FUNCTION RINTHF (I,J,K)
*                                                                      *
*   The value of RINTHF is an approximation to:                        *
*                                                                      *
*              K                                                       *
*         I ( r  * ( P (r)*Q (r) + Q (r)*P (r) ; 0 to infinity)        *
*                     I     J       I     J                            *
*                                                                      *
*   where   I ( G(r) ; Range )  denotes  the  integral  of G(r) over   *
*   Range. This is a modification of RINT for the Type B operator.     *
*                                                                      *
*   Call(s) to: [LIB92]: QUAD.                                         *
*                                                                      *
*   Written by Per O. Jonsson             Last revision: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)
*
      POINTER (PNTRPF,PF(NNNP,1)),(PNTRQF,QF(NNNP,1))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Tabulate integrand as required for SUBROUTINE QUAD
*
      MTP = MIN (MF(I),MF(J))
*
*   Value at first tabulation point is arbitrary
*
      TA(1) = 0.0D 00
      DO 1 L = 2,MTP
        TA(L) = (R(L)**K)*(PF(L,I)*QF(L,J)+QF(L,I)*PF(L,J))*RP(L)
    1 CONTINUE
*
*   Perform integration
*
      CALL QUAD (RESULT)
      RINTHF = RESULT
*
      RETURN
      END
