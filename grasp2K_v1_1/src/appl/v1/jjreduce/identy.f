************************************************************************
*                                                                      *
      SUBROUTINE IDENTY
*                                                                      *
*   Determine the position of the multireference CSFs in the csl list  * 
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
************************************************************************
*
      CHARACTER*24 NAME(2)
      CHARACTER*500 CMR1(1000),CMR2(1000),CMR3(1000)
      CHARACTER*500 LINE1C,LINE2C,LINE3C
      CHARACTER*500 LINE1,LINE2,LINE3,LINE4,LINE5
      COMMON/IDENT/NPOS(1000),NMR,NBEGIN2
*
*  Open the multireference and the file to be 
*
*
*   Open the files
*
      OPEN (UNIT = 14,FILE='rcsl.inp',FORM='FORMATTED',STATUS='OLD')

      OPEN (UNIT = 15,FILE='mrlist',FORM='FORMATTED',STATUS='OLD')

    9 FORMAT(A)

      DO I = 1,1000
        NPOS(I) = 0
      ENDDO
*
*  Determine where to start reading the CSFs in the multireference list
*
      REWIND (15)
      NBEGIN1 = 0
   55 READ(15,9) LINE1
      IF (LINE1(6:6).NE.'(') THEN
        NBEGIN1 = NBEGIN1 + 1
      GOTO 55
      ENDIF
*
*  Determine where to start reading the CSFs in the clist
*
      REWIND (14)
      NBEGIN2 = 0
   56 READ(14,9) LINE1
      IF (LINE1(6:6).NE.'(') THEN
        NBEGIN2 = NBEGIN2 + 1
      GOTO 56
      ENDIF
*
*  Read the header information in csl and csl.mr
*
      REWIND (14)
      DO I = 1,NBEGIN2
        READ(14,9) LINE1
      ENDDO
*

      REWIND (15)
      DO K = 1,NBEGIN1
        READ(15,9) LINE1
      ENDDO

      NMR = 0
   10 READ(15,9,END=25) CMR1(NMR+1)
      READ(15,9) CMR2(NMR+1)
      READ(15,9) CMR3(NMR+1)
        NMR = NMR + 1
        IF (NMR.GT.1000) THEN
          WRITE(*,*) ' Not more than 1000 CSFs in the multiref.'
          STOP
        ENDIF
      GOTO 10
   25 CONTINUE
*
      NCF = 0
   50 READ(14,9,END=60) LINE1C
      READ(14,9) LINE2C
      READ(14,9) LINE3C
      NCF = NCF + 1
      DO I = 1,500
        LINE4(I:I) = ' '
      END DO
      ICL3 = LEN_TRIM(LINE3C)
      NCL3 = 0
      DO I = 1,ICL3
        IF (LINE3C(I:I).NE.' ') THEN
          NCL3 = NCL3 + 1
          LINE4(NCL3:NCL3) = LINE3C(I:I)
        END IF
      END DO

      DO I = 1,NMR

C Note that there can be different format for line3 depending on if we use jjgen or csl.
C For this reason we can not simply do a string comparison but we need a more detailed
C treatment. Per J August 2011

        DO J = 1,500
          LINE5(J:J) = ' '
        END DO
        IML3 = LEN_TRIM(CMR3(I))
        NML3 = 0
        DO J = 1,IML3
          IF (CMR3(I)(J:J).NE.' ') THEN
            NML3 = NML3 + 1
            LINE5(NML3:NML3) = CMR3(I)(J:J)
          END IF
        END DO

        IF (CMR1(I).EQ.LINE1C.AND.CMR2(I).EQ.LINE2C) THEN
          IF (NCL3.EQ.NML3) THEN
            NSAME = 0
            DO J = 1,NCL3
              IF (LINE4(J:J).EQ.LINE5(J:J)) NSAME = NSAME + 1
            END DO
          END IF
          IF (NSAME.EQ.NCL3) THEN
            NPOS(I) = NCF
            GOTO 50
          END IF
        END IF
      END DO
      GOTO 50

   60 CONTINUE 
*
*   Close the files
*
      CLOSE (14)
      CLOSE (15)


      RETURN
      END
