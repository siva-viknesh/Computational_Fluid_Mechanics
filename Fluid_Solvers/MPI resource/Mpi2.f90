PROGRAM MPI2
! Aims at exhibiting assigning ranks to processors
USE MPI
INTEGER :: nprocs,rank,ierr

CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nprocs, ierr )
CALL MPI_COMM_RANK ( MPI_COMM_WORLD, rank, ierr )

    IF (rank.eq.0) THEN
        WRITE(*,*)"Master"
    ELSE
        WRITE(*,*)"Slave", rank 
    END IF

CALL MPI_FINALIZE(ierr)

END PROGRAM MPI2

