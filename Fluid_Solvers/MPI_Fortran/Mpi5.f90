PROGRAM MPI5
! This code demostrates the scatter and gather commands
USE MPI
INTEGER :: IP(100)
INTEGER :: nprocs,rank,ierr



CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
IF(rank.EQ.0)THEN
    DO i=1,100
        IP(i)=i!**2                                                             
    END DO
END IF
CALL MPI_SCATTER(IP,10,MPI_INTEGER,IP,10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    DO i=1,100
        IP(i)=IP(i)+100
    END DO
CALL MPI_GATHER(IP,10,MPI_INTEGER,IP,10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

IF(rank.EQ.0)THEN
    DO i=1,100
        WRITE(*,*)"Array element",i," rank ",rank," Value ", IP(i)
    END DO
END IF
CALL MPI_FINALIZE(ierr)

END
