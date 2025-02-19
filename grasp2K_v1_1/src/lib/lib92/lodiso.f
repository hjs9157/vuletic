************************************************************************
*                                                                      *
      SUBROUTINE LODISO
*                                                                      *
*   Loads the data from the  .iso  file.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 29 Sep 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z   ! out, ioncty and nelec not used
     :      /DEF11/FMTOAU,AUMAMU       ! in
     :      /NPAR/PARM(2),NPARM        ! out
     :      /NSMDAT/SQN,DMOMNM,QMOMB   ! out
!-----------------------------------------------------------------------
*
*   Read and echo pertinent information from  .iso  file
*
*   Atomic number
*
      READ (22,*) Z
*
*   Nuclear geometry
*
      READ (22,*)
      READ (22,*) A
      READ (22,*)
      READ (22,*) APARM
      READ (22,*)
      READ (22,*) CPARM
*
      IF (A .NE. 0.D0) THEN
         NPARM = 2
         PARM(1) = CPARM*FMTOAU
         PARM(2) = APARM*FMTOAU
      ELSE
         NPARM = 0
      ENDIF
*
*   Nuclear mass
*
      READ (22,*)
      READ (22,*) EMNAMU
*
      IF (EMNAMU .NE. 0.D0) THEN
         EMN = EMNAMU/AUMAMU
      ELSE
         EMN = 0.D0
      ENDIF
*
*   Nuclear spin and moments
*
      READ (22,*)
      READ (22,*) SQN
      READ (22,*)
      READ (22,*) DMOMNM
      READ (22,*)
      READ (22,*) QMOMB

      RETURN
      END
