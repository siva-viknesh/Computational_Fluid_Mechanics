MODULE VARIABLES
	INTEGER :: IP(100)
	INTEGER :: nprocs,rank,ierr
END MODULE

PROGRAM MPI7
!This code demostrates mpi in subroutines
	USE VARIABLES
	USE MPI

	DO i=1,100
		IP(i)=i
	END DO
	
	CALL MPI_INIT(ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	
	CALL Test1

	CALL MPI_FINALIZE(ierr)
END PROGRAM

SUBROUTINE Test1
	USE VARIABLES
	DO i=1,100
		WRITE(*,*)"My rank is: ",rank," The value is: ",IP(i)
	END DO
	
END SUBROUTINE
