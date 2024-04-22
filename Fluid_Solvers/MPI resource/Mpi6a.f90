PROGRAM MPI6
! This code demostrates mpi reduce
USE MPI
INTEGER :: IP(100),output
INTEGER :: nprocs,rank,ierr



CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

	DO i=1,100
		IP(i)=i!**2                                                             
	END DO
	
CALL MPI_REDUCE(IP,output,100,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF(rank.EQ.0)THEN
	WRITE(*,*)"Output ",output
END IF
CALL MPI_FINALIZE(ierr)

END
