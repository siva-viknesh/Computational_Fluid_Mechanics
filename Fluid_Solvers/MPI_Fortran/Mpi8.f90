MODULE VARIABLES
	INTEGER :: IP(100),guess
	INTEGER :: nprocs,rank,ierr
END MODULE

PROGRAM MPI7
!This code demostrates mpi in subroutines
	USE VARIABLES
	USE MPI

	DO i=1,100
		IP(i)=i**2
	END DO
	
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	guess=RF(0)	
!	CALL Test1
	WRITE(*,*)guess
	CALL MPI_FINALIZE(ierr)
END PROGRAM

SUBROUTINE Test1
	USE VARIABLES

	guess=IP(rank+10)
	WRITE(*,*)guess
	
END SUBROUTINE

FUNCTION RF(IDUM)
! generates a uniformly distributed random fraction between 0 an 1
! IDUM will generally be 0, but negative values may be used to
! re-initialize the seed
      SAVE MA,INEXT,INEXTP
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DIMENSION MA(55)
      DATA IFF/0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF    = 1
        MJ     = MSEED-IABS(IDUM)
        MJ     = MOD(MJ,MBIG)
        MA(55) = MJ
        MK     = 1
        DO 50 I=1,54
          II     = MOD(21*I,55)
          MA(II) = MK
          MK     = MJ-MK
          IF (MK.LT.MZ) MK = MK+MBIG
          MJ     = MA(II)
50      CONTINUE
        DO 100 K=1,4
          DO 60 I=1,55
            MA(I) = MA(I)-MA(1+MOD(I+30,55))
            IF (MA(I).LT.MZ) MA(I) = MA(I)+MBIG
60        CONTINUE
100     CONTINUE
        INEXT  = 0
        INEXTP = 31
      END IF
200   INEXT = INEXT+1
      IF (INEXT.EQ.56) INEXT = 1
      INEXTP = INEXTP+1
      IF (INEXTP.EQ.56) INEXTP = 1
      MJ = MA(INEXT)-MA(INEXTP)
      IF (MJ.LT.MZ) MJ = MJ+MBIG
      MA(INEXT) = MJ
      RF = MJ*FAC
      IF (RF.GT.1.E-8.AND.RF.LT.0.99999999) RETURN
      GO TO 200
END FUNCTION RF

