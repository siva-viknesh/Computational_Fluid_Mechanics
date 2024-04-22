program main

!*****************************************************************************80
!
!! MAIN is the main program for HELLO_MPI.
!
!  Discussion:
!
!    HELLO is a simple MPI test program.
!
!    Each process prints out a "Hello, world!" message.
!
!    The master process also prints out a short message.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    30 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323,
!    LC: QA76.642.G76.
!
   use mpi

  integer ( kind = 4 ) error
  integer ( kind = 4 ) id
  integer ( kind = 4 ) p
  real ( kind = 8 ) wtime
!
!  Initialize MPI.
!
  call MPI_Init ( error )
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, p, error )
!
!  Get the individual process ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
!
!  Print a message.
!
  if ( id == 0 ) then

    wtime = MPI_Wtime ( )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HELLO_MPI - Master process:'
    write ( *, '(a)' ) '  FORTRAN90/MPI version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  An MPI test program.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  The number of processes is ', p
    write ( *, '(a)' ) ' '

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a)' ) '  Process ', id, ' says "Hello, world!"'

  if ( id == 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'HELLO_MPI - Master process:'
    write ( *, '(a)' ) '  Normal end of execution: "Goodbye, world!".'

    wtime = MPI_Wtime ( ) - wtime
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6,a)' ) &
      '  Elapsed wall clock time = ', wtime, ' seconds.'

  end if
!
!  Shut down MPI.
!
  call MPI_Finalize ( error )

  stop
end program main
