PROGRAM MPI6
! This code demostrates mpi reduce
USE MPI
INTEGER :: IP(100),output(100)
INTEGER :: nprocs,rank,ierr
INTEGER :: inp,outp



CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

	DO i=1,100
		IP(i)=i!**2                                          
		inp=IP(I)                   
		CALL MPI_REDUCE(Inp,outp,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
		OUTPUT(I)=outp
	END DO


CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
IF(rank.EQ.0)THEN
	DO i=1,100
	WRITE(*,*)"Output ",output(I)
END DO
END IF
CALL MPI_FINALIZE(ierr)

END
