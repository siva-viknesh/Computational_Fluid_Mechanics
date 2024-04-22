PROGRAM MPI4
!Demo of MPI bcast. Note that the Bcast must be accessed by all ranks
USE MPI
INTEGER :: IP(10)
INTEGER :: nprocs,rank,ierr

Write(*,*)"Aer "

CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
IF(rank.EQ.0)THEN
    DO i=1,10
        IP(i)=i**2
    END DO
END IF
    CALL MPI_BCAST(IP,10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    DO i=1,10
        WRITE(*,*)"Array element",i," rank ",rank," Value ", IP(i)
    END DO

CALL MPI_FINALIZE(ierr)

END
