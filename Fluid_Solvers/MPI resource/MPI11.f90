PROGRAM MPI10
!Demo of a complicated send and recieve
USE MPI
 INTEGER :: status(MPI_STATUS_SIZE)
 INTEGER :: rank,nprocs,ierr
 INTEGER :: a(1000),N,sen,rev
 
		DO 100 I = 1, 1000
			a(I) = I*2
100	CONTINUE
CALL MPI_INIT(ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
N=1000/nprocs

DO 200 I=2,N-1
			a(I) = a(I)*2
200	CONTINUE
sen=a(N)

IF(rank.NE.procs-1)THEN
	CALL MPI_SEND( sen, 1, MPI_INTEGER, rank+1,1, MPI_COMM_WORLD, ierr ) 
END IF

IF(rank.NE.0)THEN
	CALL MPI_RECV(rev, 1, MPI_INTEGER,rank,1, MPI_COMM_WORLD, status, ierr )  
	a(1)=rev
END IF

IF(rank.EQ.0)THEN
 OPEN(10,FILE='test.txt',FORM="formatted")
	DO I=1,N
		WRITE(10,*)a(I)
	END DO
 CLOSE(10)
END IF

CALL MPI_FINALIZE(ierr)
END PROGRAM
