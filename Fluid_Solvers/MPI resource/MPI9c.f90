PROGRAM MPI9b
USE MPI
INTEGER :: nprocs,rank,ierr
INTEGER :: PP(1:2000,2),PP_local(1:2000,2),displs(8),rcvcount(8)
INTEGER :: INI,FIN,CPR,NM,N,X,Y


CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	NM=13
	N=0
	DO I=1,NM
		PP(I,1)=I*100
		PP(I,2)=I*200
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
	X=PP(N+INI,1)
	X=X+1
	Y=PP(N+INI,2)
	Y=Y+2
	PP_local(N,1)=X
	PP_local(N,2)=Y
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

do i = 1,2
!CALL MPI_GATHERV(PP_local(1:CPR,i),CPR,MPI_INTEGER,PP(INI+1,i),rcvcount,displs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_ALLGATHERV(PP_local(1:CPR,i),CPR,MPI_INTEGER,PP(INI+1,i),rcvcount,displs,MPI_INTEGER,MPI_COMM_WORLD,ierr)
end do
!CALL MPI_GATHERV(PP_local(2,1:CPR),CPR,MPI_INTEGER,PP(2,INI+1),rcvcount,displs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!CALL MPI_ALLGATHERV(PP_local(1:CPR),CPR,MPI_INTEGER,PP(INI+1),rcvcount,displs,MPI_INTEGER,MPI_COMM_WORLD,ierr)
IF (rank.EQ.1) THEN
	DO i = 1,NM
		WRITE(*,*)i,PP(I,1),PP(I,2)
	END DO
END IF
CALL MPI_FINALIZE(ierr)
END PROGRAM
