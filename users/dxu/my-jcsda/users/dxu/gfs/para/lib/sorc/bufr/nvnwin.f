      FUNCTION NVNWIN(NODE,LUN,INV1,INV2,INVN,NMAX)

C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM:    NVNWIN
C   PRGMMR: WOOLLEN          ORG: NP20       DATE: 1994-01-06
C
C ABSTRACT: THIS FUNCTION LOOKS FOR AND RETURNS ALL OCCURRENCES OF A
C   SPECIFIED NODE WITHIN THE PORTION OF THE CURRENT SUBSET BUFFER
C   BOUNDED BY THE INDICES INV1 AND INV2.  THE RESULTING LIST IS A
C   STACK OF "EVENT" INDICES FOR THE REQUESTED NODE. 
C
C PROGRAM HISTORY LOG:
C 1994-01-06  J. WOOLLEN -- ORIGINAL AUTHOR
C 1998-07-08  J. WOOLLEN -- REPLACED CALL TO CRAY LIBRARY ROUTINE
C                           "ABORT" WITH CALL TO NEW INTERNAL BUFRLIB
C                           ROUTINE "BORT"
C 1999-11-18  J. WOOLLEN -- THE NUMBER OF BUFR FILES WHICH CAN BE
C                           OPENED AT ONE TIME INCREASED FROM 10 TO 32
C                           (NECESSARY IN ORDER TO PROCESS MULTIPLE
C                           BUFR FILES UNDER THE MPI)
C 2003-11-04  S. BENDER  -- ADDED REMARKS/BUFRLIB ROUTINE
C                           INTERDEPENDENCIES
C 2003-11-04  D. KEYSER  -- MAXJL (MAXIMUM NUMBER OF JUMP/LINK ENTRIES)
C                           INCREASED FROM 15000 TO 16000 (WAS IN
C                           VERIFICATION VERSION); UNIFIED/PORTABLE FOR
C                           WRF; ADDED DOCUMENTATION (INCLUDING
C                           HISTORY); OUTPUTS MORE COMPLETE DIAGNOSTIC
C                           INFO WHEN ROUTINE TERMINATES ABNORMALLY OR
C                           UNUSUAL THINGS HAPPEN
C 2009-03-23  J. ATOR    -- USE 1E9 TO PREVENT OVERFLOW WHEN
C                           INITIALIZING INVN; USE ERRWRT
C 2009-03-31  J. WOOLLEN -- ADDED ADDITIONAL DOCUMENTATION
C
C USAGE:    NVNWIN (NODE, LUN, INV1, INV2, INVN, NMAX)
C   INPUT ARGUMENT LIST:
C     NODE     - INTEGER: JUMP/LINK TABLE INDEX TO LOOK FOR
C     LUN      - INTEGER: I/O STREAM INDEX INTO INTERNAL MEMORY ARRAYS
C     INV1     - INTEGER: STARTING INDEX OF THE PORTION OF THE SUBSET
C                BUFFER IN WHICH TO LOOK
C     INV2     - INTEGER: ENDING INDEX OF THE PORTION OF THE SUBSET
C                BUFFER IN WHICH TO LOOK
C     NMAX     - INTEGER: DIMENSIONED SIZE OF INVN; USED BY THE
C                FUNCTION TO ENSURE THAT IT DOES NOT OVERFLOW THE
C                INVN ARRAY
C
C   OUTPUT ARGUMENT LIST:
C     INVN     - INTEGER: ARRAY OF STACK "EVENT" INDICES FOR NODE
C     NVNWIN   - INTEGER: NUMBER OF INDICES RETURNED WITHIN INVN
C
C REMARKS:
C    THIS ROUTINE CALLS:        BORT     ERRWRT
C    THIS ROUTINE IS CALLED BY: UFBEVN
C                               Normally not called by any application
C                               programs.
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 77
C   MACHINE:  PORTABLE TO ALL PLATFORMS
C
C$$$

      INCLUDE 'bufrlib.prm'

      COMMON /USRINT/ NVAL(NFILES),INV(MAXSS,NFILES),VAL(MAXSS,NFILES)
      COMMON /QUIET / IPRT

      CHARACTER*128 BORT_STR
      DIMENSION     INVN(NMAX)
      REAL*8        VAL

C----------------------------------------------------------------------
C----------------------------------------------------------------------

      NVNWIN = 0

      IF(NODE.EQ.0) THEN
         IF(IPRT.GE.1) THEN
      CALL ERRWRT('+++++++++++++++++++++WARNING+++++++++++++++++++++++')
      CALL ERRWRT('BUFRLIB: NVNWIN - NODE=0, IMMEDIATE RETURN')
      CALL ERRWRT('+++++++++++++++++++++WARNING+++++++++++++++++++++++')
      CALL ERRWRT(' ')
         ENDIF
         GOTO 100
      ENDIF

      DO I=1,NMAX
         INVN(I) = 1E9
      ENDDO

C  SEARCH BETWEEN INV1 AND INV2
C  ----------------------------

      DO N=INV1,INV2
         IF(INV(N,LUN).EQ.NODE) THEN
            IF(NVNWIN+1.GT.NMAX) GOTO 900
            NVNWIN = NVNWIN+1
            INVN(NVNWIN) = N
         ENDIF
      ENDDO

C  EXITS
C  -----

100   RETURN
900   WRITE(BORT_STR,'("BUFRLIB: NVNWIN - THE NUMBER OF EVENTS, '//
     . 'NVNWIN (",I5,") EXCEEDS THE LIMIT, NMAX (",I5,")")') NVNWIN,NMAX
      CALL BORT(BORT_STR)
      END
