PROGRAM MPI12
USE MPI
INTEGER :: rank,nprocs,ierr
INTEGER :: ARRAY(10),SENDER(10),RECIEVER(10)
INTEGER :: destin,source,tag,Ssize,Rsize
INTEGER :: Sdata,Rdata
INTEGER :: status(MPI_STATUS_SIZE)

CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

	DO I=1,10
		ARRAY(I) = I+(rank*10)
	END DO

	DO I = 1,5
		SENDER(I)=ARRAY(I)
	END DO


	DO I=1,5
		Sdata = SENDER(I)
		Ssize = 1
		tag = I
		IF (rank.EQ.0)destin = rank+1
		IF (rank.EQ.1)destin = rank-1
		CALL MPI_SEND( Sdata, Ssize, MPI_INTEGER,destin,tag, MPI_COMM_WORLD, ierr ) 
!		WRITE(*,*)"SENT",rank,I,Sdata
	END DO

	DO I=1,5
		Rsize = 1
		tag = I
		IF (rank.EQ.0)Source = rank+1
		IF (rank.EQ.1)Source = rank-1		
		CALL MPI_RECV( Rdata, Rsize, MPI_INTEGER,source,tag, MPI_COMM_WORLD, status, ierr )  
!		WRITE(*,*)"RECIEVED",rank,I,Rdata
		RECIEVER(I) = Rdata
	END DO

	DO I = 1,5
		ARRAY(I)=RECIEVER(I)
	END DO

IF(rank.EQ.0)THEN
DO I=1,10
WRITE(*,*)rank,ARRAY(I)
END DO
END IF
CALL MPI_FINALIZE(ierr)

END PROGRAM
