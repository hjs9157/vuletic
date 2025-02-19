************************************************************************
*                                                                      *
      FUNCTION VPINTF (IA,IB)
*                                                                      *
*   Computes nuclear vacuum polarization integrals.                    *
*                                                                      *
*   Call(s) to: [LIB92]: QUAD.                                         *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)      
CGG      PARAMETER (NNN1 = NNNP+10)      
CGG      PARAMETER (NNNW = 120)
      LOGICAL LDBPR
      CHARACTER*2 NH
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NCDIST/ZDIST(NNNP)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      MTP = MIN (MF(IA),MF(IB))
      TA(1) = 0.0D 00
      DO 1 K = 2,MTP
         TA(K) = (PF(K,IA)*PF(K,IB)+QF(K,IA)*QF(K,IB))
     :           *ZDIST(K)
    1 CONTINUE
      CALL QUAD (VPINTF)
*
      IF (LDBPR(9)) WRITE (99,300) NP(IA),NH(IA),NP(IB),NH(IB),VPINTF
*
      RETURN
*
  300 FORMAT (/'VPINTF: V (',1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
*
      END
