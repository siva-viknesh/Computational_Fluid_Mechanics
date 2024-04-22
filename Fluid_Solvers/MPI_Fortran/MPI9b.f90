PROGRAM MPI9b
USE MPI
INTEGER :: nprocs,rank,ierr
INTEGER :: PP(1:2000),PP_local(1:2000),displs(8),rcvcount(8)
INTEGER :: INI,FIN,CPR,NM,N


CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	NM=1003
	N=0
	DO I=1,NM
		PP(I)=I*100
	END DO

	CPR = NM/nprocs
	INI = rank*CPR
	FIN = INI+CPR
	IF(rank.EQ.(nprocs-1))THEN
		FIN = NM
	END IF
	CPR=FIN-INI
	!WRITE(*,*)rank,INI,FIN,CPR
10 N=N+1
IF(N.LE.CPR) THEN
	X=PP(N+INI)
	X=X+1
	PP_local(N)=X
	GO TO 10
END IF

IF (rank.EQ.0) THEN
	DO i=1,nprocs
		rcvcount(i) = CPR
		IF(i.EQ.nprocs) rcvcount(i) = NM-(CPR*(nprocs-1))	
	END DO
END IF
CALL MPI_BCAST(rcvcount,nprocs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	DO i=1,nprocs
		displs(i) = (i-1-rank)*(NM/nprocs)
	END DO
	DO i=1,nprocs
		WRITE(*,*)i,rcvcount(i),displs(i)
	END DO


CALL MPI_GATHERV(PP_local(1:CPR),CPR,MPI_INTEGER,PP(INI+1),rcvcount,displs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!CALL MPI_ALLGATHERV(PP_local(1:CPR),CPR,MPI_INTEGER,PP(INI+1),rcvcount,displs,MPI_INTEGER,MPI_COMM_WORLD,ierr)
IF (rank.EQ.0) THEN
	DO i = 1,NM
		WRITE(*,*)PP(I)
	END DO
END IF
CALL MPI_FINALIZE(ierr)
END PROGRAM
