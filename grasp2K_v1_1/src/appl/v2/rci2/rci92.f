************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***             ******    *****   ****   *****    *****              ***
***             **   **  **   **   **   **   **  **   **             ***
***             **   **  **        **   **   **       **             ***
***             ******   **        **    *****       **              ***
***             **  **   **        **      **       **               ***
***             **   **  **   **   **     **      **                 ***
***             **   **   *****   ****   **      *******             ***
***                                                                  ***
***          Relativistic Configuration-Interaction Program          ***
***                                                                  ***
***   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
***   Grant, and C. F. Fischer, 1990).                               ***
***                                                                  ***
***                            GRASP92                               ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM RCI92
*                                                                      *
*   Entry routine for RCI92. Controls the entire computation.          *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC, SETCON.                                *
*               [RCI92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIX,       *
*                        SETRES, SETSUM, STRSUM.                       *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
*   Updated by Xinghong He                Last revision: 23 Jun 1998   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)

! cpath uses

      CHARACTER*128 NAME, tmpdir, permdir, isofile

      PARAMETER (nblk0 = 20)
      CHARACTER*8 idblk(nblk0)

      LOGICAL GETYN,YES

      COMMON/DEFAULT/NDEF
     :      /BLIM/IPRERUN,NCSFPRE,COEFFCUT1,COEFFCUT2
     :      /WHERE/IMCDF

      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN

! Memories allocated in setmix/lodmix/lodstate/items tree

      ! ...lib92/items
      POINTER (PCCMIN,ICCMIN(1))
      COMMON/DEF7/PCCMIN,NCMIN,NCMAX   ! NCMAX not used throughout

      ! ...lodmix
      POINTER (pncfblk, ncfblkdum  )
      COMMON/hblock/nblock, pncfblk

      ! ...lodmix
      POINTER (pnevblk, nevblk(1))
      POINTER (pncmaxblk, ncmaxblk(1))
      COMMON/hblock2/pnevblk, pncmaxblk

      ! ...lodmix
      POINTER (pidxblk, idxblk(1))
      COMMON/blkidx/pidxblk

! Memories allocated in setres/getcid

      POINTER (piccutblk, iccutblk(1))
      COMMON/iccu/piccutblk

      COMMON/iounit/istdi,istdo,istde

! Things for timing

      INTEGER ncount1, ncount2, ncount_rate, ncount_max

      CHARACTER chdate*8, chtime*10, chzone*5
               !ccyymmdd  hhmmss.sss  Shhmm
      INTEGER  nYMDUHMSM(8)
               !Year Month Day Universal Hour Minute Sesond Millisecond

		CHARACTER str*8, msg*128
!-----------------------------------------------------------------------

      imcdf = 26	! Unit for rci.res file
      IPRERUN = 0
      myid = 0
      nprocs = 1
      open(UNIT=31,STATUS="SCRATCH",FORM="FORMATTED")

!
! Start timing
!
      CALL STARTTIME (ncount1, 'RCI2')
CGG      CALL SYSTEM_CLOCK (ncount1, ncount_rate, ncount_max)

CGG      CALL DATE_AND_TIME (chdate, chtime, chzone, nYMDUHMSM)
CGG      msg = ' Date: ' // chdate //
CGG     &      ' Time: ' // chtime //
CGG     &      ' Zone: ' // chzone
CGG      PRINT *, msg
!
! Get NDEF 
!
CGG         WRITE (istde,*) 'RCI2: Execution begins ...'
CGG         WRITE (istde,*)
         WRITE (istde,'(A)',ADVANCE='NO') 'Default settings? '
         YES = GETYN ()
         IF (YES) THEN
            NDEF = 0
         ELSE
            NDEF = 1
         ENDIF
!
! Get name of the state (used in files like <name>.c, <name>.s)
!
         DO
            WRITE (istde,'(A)',ADVANCE='NO') 'Name of state: '
            READ (*,'(A)') NAME
            K = INDEX (NAME,' ')
            IF (K .GT. 1) EXIT
            WRITE (istde,*) 'Name may not start with a blank. redo...'
         ENDDO

!         ...Form the full name of the files used on node-0

         lenname = LEN_TRIM (NAME)
         isofile = 'isodata'
         print *, 'isofile = ', isofile(1:LEN_TRIM (isofile))
         print *, 'name = ', name(1:LEN_TRIM (name))

   99 CONTINUE
!
! Check compatibility of plant substitutions. 
!
      PRINT *, 'Calling CHKPLT...'
      CALL CHKPLT ('RCI92')

!
! In SETDBG of this version all control logicals are set to 
! false thus no debug output will be made
!
      PRINT *, 'Calling SETDBG...'
      CALL SETDBG
!
! Perform machine- and installation-dependent setup
!
      PRINT *, 'Calling SETMC...'
      CALL SETMC
!
! Set up the physical constants 
!
      PRINT *, 'Calling SETCON...'
      CALL SETCON
!
! Open summary file
!
      PRINT *, 'Calling SETSUM...'
      CALL SETSUM (NAME)

      PRINT *, 'Calling setcsl...'
      CALL setcsl (name(1:lenname) // '.c', ncore, nblk0, idblk)
!
! Set up the  .res  file; determine if this is a restart.
!
      PRINT *, 'Calling SETRES...'
      CALL SETRES (isofile, name(1:lenname) // '.w', idblk)
*
*   Open the  .mix  file; determine the eigenpairs required
*
      PRINT *, 'Calling SETMIX...'
      CALL SETMIX (NAME, idblk)
*
*   Append a summary of the inputs to the  .sum  file
*
      PRINT *, 'Calling STRSUM...'
      CALL STRSUM
*
*   Set up the table of logarithms of factorials
*
      PRINT *, 'Calling FACTT...'
      CALL FACTT
*
*   Calculate all the needed Rk integrals 
*
      PRINT *, 'Calling GENINTRK...'
      CALL GENINTRK (myid, nprocs, ndum, j2max)
*
*   Proceed with the CI calculation
*
      PRINT *, 'Calling MATRIX...'

      CALL MATRIX (ncore, (j2max))

      IF (IPRERUN .EQ. 1) THEN
         IPRERUN = 2
         GOTO 99
      ENDIF

      IF (myid .EQ. 0) THEN
         PRINT *
         PRINT *
         PRINT *, 'Finish time, Statistics'
         PRINT *
      ENDIF
      CALL STOPTIME (ncount1, 'RCI2')
CGG      CALL SYSTEM_CLOCK (ncount2, ncount_rate, ncount_max)
CGG      ncount2 = ncount2 - ncount1
CGG      nseconds = ncount2 / ncount_rate
CGG      WRITE (str, '(I8)') nseconds
CGG      msg = str // ' seconds '
CGG      PRINT *, msg

CGG      CALL DATE_AND_TIME (chdate, chtime, chzone, nYMDUHMSM)

CGG      msg = ' Date: ' // chdate //
CGG     &      ' Time: ' // chtime //
CGG     &      ' Zone: ' // chzone
CGG      PRINT *, msg
*
*   Print completion message
*
CGG      PRINT *, 'RCI2: Execution complete.'
*
      STOP
      END
