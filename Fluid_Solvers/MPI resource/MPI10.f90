PROGRAM MPI10
!Demo of Simple send and recieve
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
IF(rank.EQ.0)THEN
	sen=12
	CALL MPI_SEND( sen, 1, MPI_INTEGER,1,1, MPI_COMM_WORLD, ierr ) 
END IF

IF(rank.EQ.1)THEN
	CALL MPI_RECV(rev, 1, MPI_INTEGER,0,1, MPI_COMM_WORLD, status, ierr )  
	WRITE(*,*)rev
END IF

CALL MPI_FINALIZE(ierr)
END PROGRAM
