      PROGRAM extmix
*
* Extract mixing coefficients and the CSF from files
*   <name>.c, <name>.m / <name>.cm 
*
      IMPLICIT DOUBLE PRECISION (a-h, o-z)
      COMMON/iounit/istdi,istdo,istde
      CHARACTER*100, line(3), g92mix*6, head*500
      CHARACTER*64 StrInput, basnam, from*1, suffix*3, filnam*69, dotc*2
      LOGICAL sort, getyn, first_of_the_block
      DATA nfmix,nfcsf,nfout,nfscratch/20,21,22,23/
      DATA dotc/'.c'/

      POINTER (pntiset, iset(1))
      POINTER (pntjset, jset(1))
      POINTER (pnteval, eval(1))
      POINTER (pntevec, evec(1))

************************************************************************
*  Input conversation on the purpose, filenames, control parameters
************************************************************************

*     ...General description
      WRITE (istde,*) 'extmix - Extract mixing coefficients and CSF'
      WRITE (istde,*)  'Files to be read :'
      WRITE (istde,*)  '    <name>.c,  <name>.m / <name>.cm'
      WRITE (istde,*)  'Files to be written :'
      WRITE (istde,*)  '    rcsl.out, screen output'

*     ...Asking the base name

  123 WRITE (istde,*) 'Enter the file <name>'
      READ (istdi, '(A)') StrInput
      IF (LEN_TRIM (StrInput) .EQ. 0) THEN
         WRITE (istde,*) 'Blank line not acceptable, redo'
         GOTO 123
      ENDIF

      basnam = ADJUSTL (StrInput)

*     ...Determing the suffix

  234 WRITE (istde,*) 'Mixing coefficients from CI or RSCF calc. ?'
      WRITE (istde,*) '   a -- CI ;      b -- RSCF'
      READ (istdi,'(A)') StrInput
      StrInput = ADJUSTL (StrInput)
      from = StrInput(1:1)
      IF (from .EQ. 'a' .OR. from .EQ. 'A') THEN
         suffix = '.cm'
      ELSEIF (from .EQ. 'b' .OR. from .EQ. 'B') THEN
         suffix = '.m'
      ELSE
         WRITE (istde,*) StrInput(1:LEN_TRIM (StrInput)), 
     &                  ' is not a valid choice, redo'
         GOTO 234
      ENDIF

*     ...Form the mixing coefficient filename and check existence

      filnam = basnam(1:LEN_TRIM (basnam)) // suffix
      INQUIRE (FILE = filnam, EXIST = sort) ! Borrow sort
      IF (.NOT. sort) THEN
         WRITE (istde,*) 'File "', filnam(1:LEN_TRIM (filnam)),
     &             '" does not exist, redo'
         GOTO 123
      ENDIF

*     ...The cut-off parameter

      WRITE (istde,*) 'Enter the cut-off value for the coefficients ',
     &                '[0--1]'
      READ (istdi,*) cutoff
      cutoff = ABS (cutoff)

*     ...Sort or not

      WRITE (istde,*) 'Sort the eigen vector components ? (y/n)'
      sort = getyn()

************************************************************************
*  Open mix file, check header
************************************************************************
      OPEN (nfmix, FILE = filnam, FORM = 'UNFORMATTED', STATUS = 'OLD'
     &      , IOSTAT = IOS)
      IF (IOS .NE. 0) THEN
         WRITE (istde,*) 'Failed to to open file "', 
     &                filnam(1:LEN_TRIM (filnam)), '"'
         CLOSE (nfmix)
         STOP
      ENDIF

      READ (nfmix) g92mix
      IF (g92mix .NE. 'G92MIX') 
     &   STOP 'Not a mixing coefficient file'

      READ (nfmix) nelec, ncftot, nw, nvectot, nvecsiz, nblock
      PRINT *
      PRINT '(2X,A,2X,I2,2X,A,2X,I8,2X,A,2X,I2,2X,A,2X,I2)', 
     &       'nblock = ', nblock, '  ncftot = ', ncftot, 
     &                '  nw = ', nw, '  nelec = ', nelec
      PRINT *

************************************************************************
*  Open CSF file
************************************************************************
      filnam = basnam(1:LEN_TRIM (basnam)) // dotc
      OPEN (nfcsf, FILE = filnam, FORM = 'FORMATTED', STATUS = 'OLD'
     &      , IOSTAT = IOS)
      IF (IOS .NE. 0) THEN
         WRITE (istde,*) 'Failed to to open file "', 
     &                filnam(1:LEN_TRIM (filnam)), '"'
         STOP
      ENDIF

      OPEN (nfout, FILE = 'rcsl.out', STATUS = 'UNKNOWN')

      ! The header lines of the CSF
      DO i = 1, 5
         READ (nfcsf,'(A)') head
         WRITE (nfout,'(A)') head(1:LEN_TRIM (head))
      ENDDO

************************************************************************
*  Loop over blocks. For each block
*    Read/write block info
*    Read CSF and write to a direct access file
*    Proceed if nevblk > 0
*       Allocate memory and read mix files
************************************************************************
      DO 432 jblock = 1, nblock

         first_of_the_block = .TRUE.

         READ (nfmix) nb, ncfblk, nevblk, iatjp, iaspa
      PRINT *, ('=',i=1,75)
      PRINT '(2X,A,2X,I2,2X,A,2X,I8,3(2X,A,2X,I2))', 
     &      'nb = ', nb, 'ncfblk = ', ncfblk, 
     &      'nevblk = ', nevblk, '2J+1 = ', iatjp, 'parity = ',iaspa
      WRITE (istde, '(2X,A,2X,I2,2X,A,2X,I8,3(2X,A,2X,I2))') 
     &      'nb = ', nb, 'ncfblk = ', ncfblk, 
     &      'nevblk = ', nevblk, '2J+1 = ', iatjp, 'parity = ',iaspa
      PRINT *, ('=',i=1,75)
         IF (jblock .NE. nb) STOP 'jblock .NE. nb'

         CALL iocsf (nfcsf, nfscratch, jblock, ncfblk, line)

         IF (nevblk .LE. 0) GOTO 432

         CALL alloc (pnteval, nevblk, 8)
         CALL alloc (pntevec, nevblk*ncfblk, 8)
         CALL alloc (pntiset, ncfblk, 4)

         !READ (nfmix) (ivec(i), i = 1, nevblk)
         READ (nfmix) (ivecdum, i = 1, nevblk)
         READ (nfmix) eav, (eval(i), i = 1, nevblk)
         READ (nfmix) ((evec(i + (j-1)*ncfblk ), 
     &                  i = 1, ncfblk), j = 1, nevblk)

*        ...Find the "OR" set
         icount = 0
         DO icf = 1, ncfblk
            DO ip = 1, nevblk
               IF (ABS ( evec(icf+(ip-1)*ncfblk) ) .GT. cutoff) THEN
                  icount = icount + 1
                  iset(icount) = icf
                  EXIT
               ENDIF
            ENDDO
         ENDDO

*        ...Make a copy of the original set which is to be altered if
*           sorted
         IF (sort) THEN
            CALL alloc (pntjset, icount, 4)
            DO i = 1, icount
               jset(i) = iset(i)
            ENDDO
         ENDIF

         PRINT *, 'Average Energy = ', eav, '    ncf_reduced = ', icount

         DO 321 ip = 1, nevblk
            layer = (ip-1)*ncfblk
	    eval(ip) = eval(ip) + eav

            IF (sort) THEN
               DO i = 1, icount
                  icf = iset(i)
                  maxbub = i
                  bubmax = ABS (evec(icf+layer))
                  DO j = i+1, icount
                     jcf = iset(j)
                     IF (ABS ( evec(jcf+layer) ) .GT. bubmax)
     &                                      THEN
                        maxbub = j
                        bubmax = ABS (evec(jcf+layer))
                     ENDIF
                  ENDDO

                  iset(i) = iset(maxbub)
                  iset(maxbub) = icf
               ENDDO
            ENDIF

            PRINT *
            PRINT *, 'Energy = ', eval(ip), '    Coefficients and CSF :'
            PRINT *

            DO i = 1, icount
               icf = iset(i)
               PRINT *, i, evec(icf + layer)
               READ (nfscratch, REC = icf) line
               PRINT *, line(1)(1:LEN_TRIM (line(1)))
               PRINT *, line(2)(1:LEN_TRIM (line(2)))
               PRINT *, line(3)(1:LEN_TRIM (line(3)))

               IF (first_of_the_block) THEN
                  WRITE (nfout,'(A)') line(1)(1:LEN_TRIM (line(1)))
                  WRITE (nfout,'(A)') line(2)(1:LEN_TRIM (line(2)))
                  WRITE (nfout,'(A)') line(3)(1:LEN_TRIM (line(3)))
               ENDIF

*              ...Recover original set
               IF (sort) iset(i) = jset(i)
            ENDDO

            first_of_the_block = .FALSE.

  321    CONTINUE

         IF (jblock .LT. nblock) WRITE (nfout,'(A)') ' *'

         CALL dalloc (pntiset)
         CALL dalloc (pnteval)
         CALL dalloc (pntevec)
         IF (sort) CALL dalloc (pntjset)

	 CLOSE (nfscratch, STATUS = 'DELETE')

  432 CONTINUE

      CLOSE (nfmix)
      CLOSE (nfcsf)
      CLOSE (nfout)

      STOP 'Successful'
      END

      SUBROUTINE iocsf (nfcsf, nfscratch, jblock, ncfblk, line)
      IMPLICIT NONE
      INTEGER nfcsf, nfscratch, jblock, ncfblk, icf, i
      CHARACTER*(*) line(3), star*2
      !CHARACTER star*2

      OPEN (nfscratch, 
!     & FILE = 'scratct', 
     & STATUS = 'SCRATCH',
     & ACCESS = 'DIRECT',
     & RECL = 300)

      DO icf = 1, ncfblk
         READ (nfcsf,'(A)') line(1)
         READ (nfcsf,'(A)') line(2)
         READ (nfcsf,'(A)') line(3)
         WRITE (nfscratch,REC = icf) line
      ENDDO

      READ (nfcsf,'(A)',END=123) star
      IF (star .NE. ' *') STOP ' CSF read wrong'

  123 CONTINUE

      RETURN
      END
