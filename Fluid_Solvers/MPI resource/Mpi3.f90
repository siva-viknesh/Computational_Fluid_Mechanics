PROGRAM MPI3
!Proves that the array delcared b4 MPI Calls is shared by all ranks
USE MPI
INTEGER :: IP(10)
INTEGER :: nprocs,rank,ierr

    DO i=1,10
        IP(i)=i**2
    END DO

CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

    DO i=1,10
        WRITE(*,*)"Array element",i," rank ",rank," Value ", IP(i)
    END DO

CALL MPI_FINALIZE(ierr)

END PROGRAM MPI3
