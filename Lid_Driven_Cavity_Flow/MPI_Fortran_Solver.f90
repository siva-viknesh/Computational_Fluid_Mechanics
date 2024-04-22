MODULE VAR
  implicit none
  integer, parameter :: nglo=1025,mglo=nglo
  integer :: iter,itermax,iter_write
  integer :: restart
  double precision :: dt,time,Re,maxtime,time_write 
  double precision , parameter ::pi=4.0d0*datan(1.0d0)

  integer :: rank, nprocs
  integer :: n,m
  integer :: ist_glo,jst_glo
  integer :: ien_glo,jen_glo
  integer :: ist,jst,ien,jen
  integer :: ist_wb,jst_wb,ien_wb,jen_wb !_wb is with physical boundary included
  integer :: all_interfaceids_i_size !=2*(nprocs-1)
  integer,dimension(:), allocatable :: glo_interfaceids_i !dimension will be 2*(nprocs-1)

  double precision, parameter :: dx=1.0d0/dble(nglo-1), dy=1.0d0/dble(mglo-1)

  double precision :: xg(nglo), yg(mglo)
!  double precision :: psig(0:nglo+1,0:mglo+1), vortg(nglo,mglo)
!  double precision :: dvort_dx_g(0:nglo+1,0:mglo+1), d2vort_dx2_g(nglo,mglo)
  double precision :: psi0
  double precision :: max_diff,max_dvort_dt
  double precision, dimension(:,:,:,:), allocatable :: D_inv !For CCD scheme

  double precision, dimension(:,:), allocatable :: psi,vort,vort_base,vort_rk,conv_term,u,v 
  double precision, dimension(:,:), allocatable :: dvort_dx,d2vort_dx2,dvort_dy,d2vort_dy2 
END MODULE VAR

MODULE BICG_VAR
  implicit none
  integer , parameter :: bicg_iter_max=10000
  double precision, dimension(:,:), allocatable :: pk,rk,tk,vk,r00,uk
  double precision :: r00rk_g, r00uk_g,vktk_g,vkvk_g
  double precision :: alpha1_g,zi_g,beta1_g,maxdiff_g
  double precision, parameter :: eps=1.0d-8 
END MODULE BICG_VAR

MODULE RK4_VAR
  implicit none
  integer  :: rk_stage
  double precision, parameter :: f_q(4)=(/ 0.50d0,0.50d0,1.0d0,1.0d0 /)
  double precision, parameter :: fac(4)=(/ (1.0d0/6.0d0),(1.0d0/3.0d0),(1.0d0/3.0d0),(1.0d0/6.0d0) /)
END MODULE RK4_VAR

subroutine Partition_Domain (my_rank, my_nprocs)
use VAR
use BICG_VAR
implicit none
integer :: i_rank,nr
integer :: my_rank, my_nprocs, remainder

!Domain decomposition along x or xi direction only
rank = my_rank
nprocs = my_nprocs

remainder = mod(nglo,nprocs)

nr = nglo/nprocs
if(rank.lt.remainder)then
n = nr+1
else 
n = nr
end if

m = mglo

!Global indexing{
if(rank.eq.0)then
ist_glo = 1
jst_glo = 1
else if(rank.le.remainder)then
ist_glo = (rank*(nr+1))+1
jst_glo = 1
else
ist_glo = (remainder*(nr+1)) + ((rank-remainder)*nr) + 1
jst_glo = 1
end if

if(rank.eq.0)then
ien_glo = ((rank+1)*(nr+1))
jen_glo = mglo
else if(rank.lt.remainder)then
ien_glo = ((rank+1)*(nr+1))
jen_glo = mglo
else
ien_glo = (remainder*(nr+1))+((rank-remainder+1)*nr)
jen_glo = mglo
end if

!if(rank.eq.nprocs-1)then
!ien_glo = nglo
!jen_glo = mglo
!else
!ien_glo = (rank+1)*nr
!jen_glo = mglo
!end if
!}

!Local indexing{
  if(rank.eq.0)then
   ist = 2
  else
   ist = 1
  end if

  if(rank.eq.nprocs-1)then
   ien = n-1
  else
   ien = n
  end if

  jst = 2
  jen = m-1

ist_wb = 1
ien_wb = n

jst_wb = 1
jen_wb = m
!}

!Interface definition for global compact scheme based derivative calculation at interface{
all_interfaceids_i_size = 2*(nprocs-1)
ALLOCATE (glo_interfaceids_i(all_interfaceids_i_size)) !dimension will be 2*(nprocs-1)
do i_rank =0,nprocs-2
if(i_rank.eq.0)then
 glo_interfaceids_i(2*i_rank+1) = ((i_rank+1)*(nr+1))
else if(i_rank.lt.remainder)then
 glo_interfaceids_i(2*i_rank+1) = ((i_rank+1)*(nr+1))
else
 glo_interfaceids_i(2*i_rank+1) = (remainder*(nr+1))+((i_rank-remainder+1)*nr)
end if
glo_interfaceids_i(2*i_rank+2) = glo_interfaceids_i(2*i_rank+1) + 1
end do
!}

!Allocate variables
!flow variables
ALLOCATE(psi(0:n+1,0:m+1),vort(0:n+1,0:m+1),vort_base(0:n+1,0:m+1),vort_rk(0:n+1,0:m+1))
ALLOCATE(conv_term(0:n+1,0:m+1),u(0:n+1,0:m+1),v(0:n+1,0:m+1))
ALLOCATE(dvort_dx(0:n+1,0:m+1),d2vort_dx2(0:n+1,0:m+1))
ALLOCATE(dvort_dy(0:n+1,0:m+1),d2vort_dy2(0:n+1,0:m+1))
!for exact derivative using CCD scheme
ALLOCATE(D_inv(1:all_interfaceids_i_size,1:n,1:2,1:2))
!bicg variables
ALLOCATE(pk(0:n+1,0:m+1),rk(0:n+1,0:m+1),tk(0:n+1,0:m+1),vk(0:n+1,0:m+1),r00(0:n+1,0:m+1), uk(0:n+1,0:m+1))

write(*,*)rank, ist_glo, ien_glo, n
if(rank.eq.0)then
do i_rank =1,all_interfaceids_i_size
write(*,*)glo_interfaceids_i(i_rank)
end do
end if
end subroutine Partition_Domain

subroutine Free_Allocated_Var
use VAR
use BICG_VAR
implicit none

!flow variables
DEALLOCATE(psi,vort,vort_base,vort_rk)
DEALLOCATE(conv_term,u,v)
DEALLOCATE(dvort_dx,d2vort_dx2)
DEALLOCATE(dvort_dy,d2vort_dy2)
DEALLOCATE(D_inv)
!bicg variables
DEALLOCATE(pk,rk,tk,vk,r00, uk)
end subroutine Free_Allocated_Var

subroutine INITIALIZE_VAR
use VAR
implicit none

 psi = 0.0d0
 vort = 0.0d0

end subroutine INITIALIZE_VAR

module Interface_DataTransfer
use VAR,only:rank,nprocs,n,m 
implicit none
include 'mpif.h'
contains
subroutine Datatransfer_at_Interface(mydata)
implicit none
integer :: left,right
integer :: nsendcount_left,nrecvcount_left
integer :: nsendcount_right,nrecvcount_right
integer :: ierr,trans_status(mpi_status_size)
double precision :: send_data_left(m), recv_data_left(m)
double precision :: send_data_right(m), recv_data_right(m)
double precision, dimension(:,:), allocatable,intent(inout) :: mydata

if(rank.eq.0)then
nsendcount_left = -1
nrecvcount_left = -1
nsendcount_right = m
nrecvcount_right = m
left = -1
right = rank+1
send_data_right(1:m) = mydata(n,1:m)
else if(rank.eq.nprocs-1)then
nsendcount_left = m
nrecvcount_left = m
nsendcount_right = -1
nrecvcount_right = -1
left = rank-1
right = -1
send_data_left(1:m) = mydata(1,1:m)
else
nsendcount_left = m
nrecvcount_left = m
nsendcount_right = m
nrecvcount_right = m
left = rank-1
right = rank+1
send_data_left(1:m) = mydata(1,1:m)
send_data_right(1:m) = mydata(n,1:m)
endif

if(rank.eq.0)then
call MPI_Sendrecv(send_data_right, nsendcount_right, MPI_double_precision,right, rank*right,&
                recv_data_right, nrecvcount_right, MPI_double_precision,right, rank*right,&
                MPI_COMM_WORLD, trans_status, ierr)
else if(rank.eq.nprocs-1)then
call MPI_Sendrecv(send_data_left, nsendcount_left, MPI_double_precision,left, rank*left,&
                recv_data_left, nrecvcount_left, MPI_double_precision,left, rank*left,&
                MPI_COMM_WORLD, trans_status, ierr)
else
CALL MPI_SEND(send_data_left, nsendcount_left, MPI_double_precision, left, rank*left, MPI_COMM_WORLD, ierr)
CALL MPI_SEND(send_data_right, nsendcount_right, MPI_double_precision, right, rank*right, MPI_COMM_WORLD, ierr)

CALL MPI_RECV(recv_data_left, nrecvcount_left, MPI_double_precision,left, rank*left,&
                MPI_COMM_WORLD, trans_status, ierr)
CALL MPI_RECV(recv_data_right, nrecvcount_right, MPI_double_precision,right, rank*right,&
                MPI_COMM_WORLD, trans_status, ierr)
endif

if(rank.eq.0)then
mydata(n+1,1:m) = recv_data_right(1:m)
else if(rank.eq.nprocs-1)then
mydata(0,1:m) = recv_data_left(1:m)
else
mydata(0,1:m) = recv_data_left(1:m)
mydata(n+1,1:m) = recv_data_right(1:m)
endif

end subroutine Datatransfer_at_Interface

subroutine Datatransfer2_at_Interface(mydata, mydata1)
implicit none
integer :: left,right
integer :: nsendcount_left,nrecvcount_left
integer :: nsendcount_right,nrecvcount_right
integer :: ierr,trans_status(mpi_status_size)
double precision :: send_data_left(2*m), recv_data_left(2*m)
double precision :: send_data_right(2*m), recv_data_right(2*m)
double precision, dimension(:,:), allocatable,intent(inout) :: mydata, mydata1

if(rank.eq.0)then
nsendcount_left = -1
nrecvcount_left = -1
nsendcount_right = 2*m
nrecvcount_right = 2*m
left = -1
right = rank+1
send_data_right(1:m) = mydata(n,1:m)
send_data_right(m+1:2*m) = mydata1(n,1:m)
else if(rank.eq.nprocs-1)then
nsendcount_left = 2*m
nrecvcount_left = 2*m
nsendcount_right = -1
nrecvcount_right = -1
left = rank-1
right = -1
send_data_left(1:m) = mydata(1,1:m)
send_data_left(m+1:2*m) = mydata1(1,1:m)
else
nsendcount_left = 2*m
nrecvcount_left = 2*m
nsendcount_right = 2*m
nrecvcount_right = 2*m
left = rank-1
right = rank+1
send_data_left(1:m) = mydata(1,1:m)
send_data_left(m+1:2*m) = mydata1(1,1:m)
send_data_right(1:m) = mydata(n,1:m)
send_data_right(m+1:2*m) = mydata1(n,1:m)
endif

if(rank.eq.0)then
!write(*,*)right,nsendcount_right, nrecvcount_right
call MPI_Sendrecv(send_data_right, nsendcount_right, MPI_double_precision,right, rank*right,&
                recv_data_right, nrecvcount_right, MPI_double_precision,right, rank*right,&
                MPI_COMM_WORLD, trans_status, ierr)
else if(rank.eq.nprocs-1)then
!write(*,*)left,nsendcount_left, nrecvcount_left
call MPI_Sendrecv(send_data_left, nsendcount_left, MPI_double_precision,left, rank*left,&
                recv_data_left, nrecvcount_left, MPI_double_precision,left, rank*left,&
                MPI_COMM_WORLD, trans_status, ierr)
else
CALL MPI_SEND(send_data_left, nsendcount_left, MPI_double_precision, left, rank*left, MPI_COMM_WORLD, ierr)
CALL MPI_SEND(send_data_right, nsendcount_right, MPI_double_precision, right, rank*right, MPI_COMM_WORLD, ierr)

CALL MPI_RECV(recv_data_left, nrecvcount_left, MPI_double_precision,left, rank*left,&
                MPI_COMM_WORLD, trans_status, ierr)
CALL MPI_RECV(recv_data_right, nrecvcount_right, MPI_double_precision,right, rank*right,&
                MPI_COMM_WORLD, trans_status, ierr)
endif

if(rank.eq.0)then
mydata(n+1,1:m) = recv_data_right(1:m)
mydata1(n+1,1:m) = recv_data_right(m+1:2*m)
else if(rank.eq.nprocs-1)then
mydata(0,1:m) = recv_data_left(1:m)
mydata1(0,1:m) = recv_data_left(m+1:2*m)
else
mydata(0,1:m) = recv_data_left(1:m)
mydata1(0,1:m) = recv_data_left(m+1:2*m)
mydata(n+1,1:m) = recv_data_right(1:m)
mydata1(n+1,1:m) = recv_data_right(m+1:2*m)
endif

end subroutine Datatransfer2_at_Interface
end module Interface_DataTransfer

!--------------------------------------------------------------------------------------------------------------------!
!                           Main program starts
!--------------------------------------------------------------------------------------------------------------------!
program drivencavity
USE VAR
  USE Interface_DataTransfer
  implicit none
!  include 'mpif.h'
  integer :: i,j,k,l,ierr,my_rank,my_nprocs
  double precision :: uvel,vvel 
  character(len=1024) :: filename
  double precision, parameter :: vort_monitor_at_x=0.950d0,vort_monitor_at_y=0.950d0
  integer, parameter :: i_vort_monitor_at_x=(0.950d0/dx)+1,j_vort_monitor_at_y=(0.950d0/dx)+1
  logical :: file_exists

  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,my_nprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)

  call Partition_Domain(my_rank,my_nprocs)

!{
  dt = 0.0005d0
  Re = 8300.0d0            
  psi0 = 0.0d0
  maxtime = 5.0d3
  itermax = (maxtime/dt)
  time_write = 50.0d0
  iter_write = time_write/dt

!{Initialized above
!  dx = 1.0d0/dble(nglo-1)
!  dy = 1.0d0/dble(mglo-1)
!}

  restart = 0 
!Broadcast these items to all processors for consistency
!}

  do i=1,nglo
     xg(i) = dble(i-1) * dx
  end do

  do j=1,mglo
     yg(j) = dble(j-1) * dy
  end do

  call Datatransfer2_at_Interface (psi,vort)

  call global_ccd_x

  if (restart == 0) then
     CALL INITIALIZE_VAR

     CALL STREAM_BOUNDARY
     call BOVRT
     call Datatransfer2_at_Interface (psi,vort)

!Clearing the file if it exists
        if(rank.eq.0)then
        inquire(file='residual.dat',exist=file_exists)
        if(file_exists)then
        open(12,file = 'residual.dat',status='unknown')
        else
        open(12,file = 'residual.dat',status='replace')
        end if
        close(12)
        open(12,file = 'residual.dat',position='append')

        inquire(file='vort.dat',exist=file_exists)
        if(file_exists)then
        open(13,file = 'vort.dat',status='unknown')
        else
        open(13,file = 'vort.dat',status='replace')
        end if
        close(13)
        end if
  end if

  do iter = 1, itermax

     time = iter*dt

     CALL RUNGE

if(rank.eq.0)then
     write(*,*) "**************************************************"
     write(*,*) "time =",time
     WRITE (*,*) "max dvort_dt = ", max_dvort_dt
     if( mod(iter,50) == 0) then 
        write(12,*) time, max_dvort_dt
        call flush(12)
     end if
     write(*,*) "**************************************************" 
endif

if(i_vort_monitor_at_x.ge.ist_glo.and.i_vort_monitor_at_x.le.ien_glo)then
if(j_vort_monitor_at_y.ge.jst_glo.and.j_vort_monitor_at_y.le.jen_glo)then
     if( mod(iter,50) == 0) then 
        open(13,file = 'vort.dat',position='append')
        write(13,*) time, vort(ist_wb+i_vort_monitor_at_x-ist_glo,jst_wb+j_vort_monitor_at_y-jst_glo)
        !call flush(13)
        close(13)
     end if
end if
end if

     if (mod(iter,iter_write) == 0) then

      write (filename, "(A4,F0.4,A9,I0)") "Time",time,".dat_rank",rank
      open(1,file=trim(filename)) 
      write(1,*)"VARIABLES = x y psi vort"
      write(1,*)"ZONE I=",n,", J=",m
        do j=1,m
        do i=1,n
         write(1,*)xg(ist_glo+i-1),yg(jst_glo+j-1),psi(i,j),vort(i,j)
        end do
        end do
      close(1)  
     endif

  end do
!
!        call write_global_solution_file
 
if(rank.eq.0)then
close(12) 
endif

 CALL MPI_FINALIZE (ierr)

end program drivencavity

!--------------------------------------------------------------------------------------------------------------------!
!				Main program ends
!--------------------------------------------------------------------------------------------------------------------!

SUBROUTINE BOVRT
  USE VAR
  IMPLICIT NONE
  integer :: i,j

!Bottom (except corners)
if(jst_glo.eq.1)then
  do i=ist,ien
     vort(i,1) = -(psi(i,0)-(2.0d0*psi(i,1))+psi(i,2))/(dy**2)
  end do
end if
!Top (except corners)
if(jen_glo.eq.mglo)then
  do i=ist,ien
     vort(i,m)= -(psi(i,m+1)-(2.0d0*psi(i,m))+psi(i,m-1))/(dy**2)
  end do
end if

!Left
if(ist_glo.eq.1)then
  do j=jst_wb,jen
     vort(1,j) = -(psi(0,j)-(2.0d0*psi(1,j))+psi(2,j))/(dx**2)
  end do
end if
!Right
if(ien_glo.eq.nglo)then
  do j=jst_wb,jen
     vort(n,j) = -(psi(n+1,j)-(2.0d0*psi(n,j))+psi(n-1,j))/(dx**2)
  end do
end if

!Top Left Corner
if(jen_glo.eq.mglo.and.ist_glo.eq.1)then
  vort(1,m) = (vort(1,m-1) + vort(2,m))/2.0d0
endif
!Top Right Corner
if(jen_glo.eq.mglo.and.ien_glo.eq.nglo)then
  vort(n,m) = (vort(n,m-1) + vort(n-1,m))/2.0d0
end if

END SUBROUTINE BOVRT
!--------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------!
SUBROUTINE STREAM_BOUNDARY
  USE VAR
  IMPLICIT NONE
  INTEGER :: i,j

!Bottom
if(jst_glo.eq.1)then
     psi(1:n,1) = 0.0d0
end if
!Top
if(jen_glo.eq.mglo)then
     psi(1:n,m)= 0.0d0
end if

!Left
if(ist_glo.eq.1)then
     psi(1,1:m) = 0.0d0
end if
!Right
if(ien_glo.eq.nglo)then
     psi(n,1:m) = 0.0d0
end if

!Ghost BC
!Bottom
if(jst_glo.eq.1)then
  do i=ist,ien
     psi(i,0) = (-(3.0d0*psi(i,1))+(6.0d0*psi(i,2))-psi(i,3))/2.0d0
  end do
end if
!Top
if(jen_glo.eq.mglo)then
  do i=ist,ien
     psi(i,m+1)=((6.0d0*dy)-(3.0d0*psi(i,m))+(6.0d0*psi(i,m-1))-psi(i,m-2))/2.0d0  
  end do
end if

!Left
if(ist_glo.eq.1)then
  do j=jst-1,jen
     psi(0,j) = (-(3.0d0*psi(1,j))+(6.0d0*psi(2,j))-psi(3,j))/2.0d0
  end do
end if
!Right
if(ien_glo.eq.nglo)then
  do j=jst-1,jen
     psi(n+1,j) = (-(3.0d0*psi(n,j))+(6.0d0*psi(n-1,j))-psi(n-2,j))/2.0d0
  end do
end if

END SUBROUTINE STREAM_BOUNDARY

!--------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------!
SUBROUTINE RUNGE
  USE VAR
  USE RK4_VAR
  USE Interface_DataTransfer
  IMPLICIT NONE
  integer :: i,j
  double precision :: differ,dvort_dt
  double precision :: temp_buff(2), temp_buff_glo(2)

  vort_base = vort
  vort_rk = 0.0d0

  DO rk_stage = 1,4
     CALL VTE
     CALL SOLVE_STREAMFN
     CALL BOVRT 
!Transfer psi at interface{
     call Datatransfer2_at_Interface (psi,vort)
!}
  END DO

  dvort_dt = -10.0d0

  do i=ist_wb,ien_wb
     do j=jst_wb,jen_wb
        differ = vort(i,j)-vort_base(i,j)
        if (dabs(differ)>dvort_dt) then
           dvort_dt = dabs(differ)
        end if
     end do
  end do

  dvort_dt = dvort_dt/dt

  temp_buff(1) = dvort_dt 
  CALL GLOBAL_REDUCE(1, temp_buff, temp_buff_glo, 1)
  max_dvort_dt = temp_buff_glo(1)

END SUBROUTINE RUNGE
!--------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------!
SUBROUTINE VTE
  USE VAR
  USE RK4_VAR
  IMPLICIT NONE
  INTEGER :: i,j,k
  double precision :: f1d, f2d, u1, v1, rhsf
  double precision :: theq

  call ccd_x !for vort
  call ccd_y !for vort

  do j=jst,jen
     do i=ist,ien
        u(i,j) = (psi(i,j+1)-psi(i,j-1))/(2.0d0*dy)
        v(i,j) = -(psi(i+1,j)-psi(i-1,j))/(2.0d0*dx)

        conv_term(i,j) = dvort_dx(i,j) * u(i,j)  !rhs3
        conv_term(i,j) = conv_term(i,j) + (dvort_dy(i,j) * v(i,j)) !rhs3-rhs4
     end do
  end do

if(rk_stage.lt.4)then
  do i=ist,ien
     do j=jst,jen

        rhsf= (1.0d0/Re)*(d2vort_dx2(i,j)+d2vort_dy2(i,j))-conv_term(i,j)

        ! ***********RK4 CONVERSION STARTS************** !
        theq = dt * rhsf 
        vort_rk(i,j) = vort_rk(i,j) + fac(rk_stage) * theq

           vort(i,j) = vort_base(i,j) + f_q(rk_stage) * theq
     enddo
  enddo

else
  do i=ist,ien
     do j=jst,jen

        rhsf= (1.0d0/Re)*(d2vort_dx2(i,j)+d2vort_dy2(i,j))-conv_term(i,j)

        theq = dt * rhsf 
        vort_rk(i,j) = vort_rk(i,j) + fac(rk_stage) * theq

           vort(i,j) = vort_base(i,j) + vort_rk(i,j)
        ! ***********RK4 CONVERSION ENDS**************** !
     enddo
  enddo
endif

END SUBROUTINE VTE
!--------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------!
Subroutine SOLVE_STREAMFN
  USE VAR
  USE BICG_VAR
  USE Interface_DataTransfer
  implicit none
  integer :: i,j,k,idiff,jdiff,bicg_iter
  double precision :: diff,coff
  double precision :: d2psi_dx2, d2psi_dy2
  double precision :: d2pk_dx2, d2pk_dy2
  double precision :: d2tk_dx2, d2tk_dy2
  double precision :: temp_buff(2), temp_buff_glo(2)
  double precision :: r00rk, r00uk,vktk,vkvk
  double precision :: alpha1,zi,beta1,maxdiff

  coff = - (2.0d0/(dx**2.0d0)) - (2.0d0/(dy**2.0d0))

  do i=ist,ien
     do j=jst,jen
        d2psi_dx2 = (psi(i+1,j) - 2.0d0*psi(i,j) + psi(i-1,j))/dx**2
        d2psi_dy2 = (psi(i,j+1) - 2.0d0*psi(i,j) + psi(i,j-1))/dy**2
        pk(i,j) =  ((-d2psi_dx2 - d2psi_dy2) - vort(i,j))/coff
        rk(i,j) = pk(i,j)
        r00(i,j) = pk(i,j)
     end do
  end do

  do bicg_iter = 1, bicg_iter_max

     r00rk = 0.0d0
     r00uk = 0.0d0
     vktk = 0.0d0
     vkvk = 0.0d0 

     r00rk_g = 0.0d0
     r00uk_g = 0.0d0
     vktk_g = 0.0d0
     vkvk_g = 0.0d0 

!Transfer pk at interface{
     call Datatransfer_at_Interface (pk)
!}
     do i=ist,ien
        do j=jst,jen
        d2pk_dx2 = (pk(i+1,j) - 2.0d0*pk(i,j) + pk(i-1,j))/dx**2
        d2pk_dy2 = (pk(i,j+1) - 2.0d0*pk(i,j) + pk(i,j-1))/dy**2

        uk(i,j) = (d2pk_dx2 + d2pk_dy2)/coff 
        r00rk = r00rk + r00(i,j) * rk(i,j)
        r00uk = r00uk + r00(i,j) * uk(i,j)
        end do
     end do

!Calculate global sum of r00rk and r00uk{
     temp_buff(1) = r00rk
     temp_buff(2) = r00uk
     temp_buff_glo = 0.0d0
     call GLOBAL_REDUCE (2,temp_buff,temp_buff_glo, 2)
     r00rk_g = temp_buff_glo(1)
     r00uk_g = temp_buff_glo(2)
!}

     alpha1_g = r00rk_g / r00uk_g

     if (r00rk_g == 0.0d0)exit
     if(r00uk_g.eq.0.0)then
        write(*,*)'Prob in BCG, code stops for Mean!'
        stop
     end if

     do i=ist,ien
        do j=jst,jen
           tk(i,j) = rk(i,j) - alpha1_g * uk(i,j)
        end do
     end do

!Transfer tk at interface{
     call Datatransfer_at_Interface (tk)
!}
     do i=ist,ien
        do j=jst,jen
        d2tk_dx2 = (tk(i+1,j) - 2.0d0*tk(i,j) + tk(i-1,j))/dx**2
        d2tk_dy2 = (tk(i,j+1) - 2.0d0*tk(i,j) + tk(i,j-1))/dy**2

           vk(i,j) = ((d2tk_dx2 + d2tk_dy2))/coff    
           vktk = vktk + vk(i,j) * tk(i,j)
           vkvk = vkvk + vk(i,j) * vk(i,j)
        end do
     end do

!Calculate global sum of vktk_g and vkvk_g{
     temp_buff(1) = vktk
     temp_buff(2) = vkvk
     temp_buff_glo = 0.0d0
     call GLOBAL_REDUCE (2,temp_buff,temp_buff_glo,2)
     vktk_g = temp_buff_glo(1)
     vkvk_g = temp_buff_glo(2)
!}

     zi_g = vktk_g / vkvk_g
     maxdiff = -10.0d0
     beta1 = 0.0d0

     do i=ist,ien
        do j=jst,jen
           psi(i,j) = psi(i,j) + alpha1_g * pk(i,j) + zi_g * tk(i,j)
           rk(i,j) = tk(i,j) - zi_g * vk(i,j)
           beta1 = beta1 + r00(i,j) * rk(i,j)
           if (dabs(rk(i,j)) > maxdiff) then
              maxdiff=dabs(rk(i,j))
           end if
        end do
     end do

     temp_buff(1) = maxdiff
     temp_buff_glo = 0.0d0
     call GLOBAL_REDUCE (1,temp_buff,temp_buff_glo,1)
     maxdiff_g = temp_buff_glo(1)

     if(maxdiff_g<eps) exit

!Calculate global sum of beta1_g{
     temp_buff(1) = beta1
     temp_buff_glo = 0.0d0
     call GLOBAL_REDUCE (1,temp_buff,temp_buff_glo,2)
     beta1_g = temp_buff_glo(1)
!}

     beta1_g = (beta1_g / r00rk_g) * ( alpha1_g / zi_g )
     do i = ist,ien
        do j = jst,jen
           pk(i,j) = rk(i,j) + beta1_g * ( pk(i,j) - zi_g * uk(i,j) )
        end do
     end do

  !if(rank.eq.0)write(*,*)alpha1_g,beta1_g,zi_g

  end do!main loop

     if(rank.eq.0) then
        write(*,*) "BICG took", bicg_iter, maxdiff_g!, idiff, jdiff
     end if

     CALL STREAM_BOUNDARY !Usually only normal 

end subroutine SOLVE_STREAMFN
!--------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------!
subroutine ccd_x
  use VAR
  implicit none
  integer :: i,j,k,jloop
  double precision :: ac(n,2,2),bc(n,2,2),cc(n,2,2),dc(n,2),alphac(n,2,2),det
  double precision :: betac(n,2,2),uc(0:n+1),ubar(n,2)
  double precision :: temp11, temp12, temp21, temp22, wc(n,2)
  double precision, parameter :: beta2c=-0.025d0,betanm1c=0.09d0

  call CALC_EX_DER_AT_INTERFACE

  do jloop=1,m

!{Interface derivative is a Dirichlet condition
if(rank.ne.0)then
        ubar(1,1) = dvort_dx(1,jloop) !(vort(2,jloop)-vort(0,jloop))/(2.0d0*dx)!dvort_dx(1,jloop)
        ubar(1,2) = d2vort_dx2(1,jloop)!(vort(2,jloop)-2.0d0*vort(1,jloop)+vort(0,jloop))/(dx**2.0d0)!d2vort_dx2(1,jloop)
endif
if(rank.ne.nprocs-1)then
        ubar(n,1) = dvort_dx(n,jloop)!(vort(n+1,jloop)-vort(n-1,jloop))/(2.0d0*dx)!dvort_dx(n,jloop)
        ubar(n,2) = d2vort_dx2(n,jloop)!(vort(n+1,jloop)-2.0d0*vort(n,jloop)+vort(n-1,jloop))/(dx**2.0d0)!d2vort_dx2(n,jloop)
endif
!}

  uc(0:n+1) = vort(0:n+1,jloop)  !changed here by Suman

  do i=1,n
     bc(i,1,1) = 1.0d0
     bc(i,2,2) = 1.0d0
     bc(i,1,2) = 0.0d0
     bc(i,2,1) = 0.0d0
  end do

if(rank.eq.0)then
ac(1:2,:,:) = 0.0d0
cc(1:2,:,:) = 0.0d0
  cc(2,2,1) = 9.0d0/(8.0d0*dx)
  cc(2,2,2) = -1.0d0/8.0d0
  ac(2,2,1) = -9.0d0/(8.0d0*dx)
  ac(2,2,2) = -1.0d0/8.0d0
else
ac(1,:,:) = 0.0d0
cc(1,:,:) = 0.0d0

     ac(2,1,1) = 7.0d0 / 16.0d0
     ac(2,2,2) = -1.0d0/8.0d0
     ac(2,1,2) = dx/16.0d0
     ac(2,2,1) = -9.0d0/(8.0d0*dx)

     cc(2,1,1) = 7.0d0 / 16.0d0
     cc(2,2,2) = -1.0d0/8.0d0
     cc(2,1,2) = -dx/16.0d0
     cc(2,2,1) = 9.0d0/(8.0d0*dx)
endif

if(rank.eq.nprocs-1)then
ac(n-1:n,:,:) = 0.0d0
cc(n-1:n,:,:) = 0.0d0

  cc(n-1,2,1) = 9.0d0/(8.0d0*dx)
  cc(n-1,2,2) = -1.0d0/8.0d0
  ac(n-1,2,1) = -9.0d0/(8.0d0*dx)
  ac(n-1,2,2) = -1.0d0/8.0d0

else
ac(n,:,:) = 0.0d0
cc(n,:,:) = 0.0d0

     ac(n-1,1,1) = 7.0d0 / 16.0d0
     ac(n-1,2,2) = -1.0d0/8.0d0
     ac(n-1,1,2) = dx/16.0d0
     ac(n-1,2,1) = -9.0d0/(8.0d0*dx)

     cc(n-1,1,1) = 7.0d0 / 16.0d0
     cc(n-1,2,2) = -1.0d0/8.0d0
     cc(n-1,1,2) = -dx/16.0d0
     cc(n-1,2,1) = 9.0d0/(8.0d0*dx)
endif

  do i=3,n-2
     ac(i,1,1) = 7.0d0 / 16.0d0
     ac(i,2,2) = -1.0d0/8.0d0
     ac(i,1,2) = dx/16.0d0
     ac(i,2,1) = -9.0d0/(8.0d0*dx)
     cc(i,1,1) = 7.0d0 / 16.0d0
     cc(i,2,2) = -1.0d0/8.0d0
     cc(i,1,2) = -dx/16.0d0
     cc(i,2,1) = 9.0d0/(8.0d0*dx)
  end do

if(rank.eq.0)then
  dc(1,1) = (0.50d0/dx) * (-3.0d0 * uc(1) + 4.0d0 * uc(2) - uc(3))
  dc(1,2) = (1.0d0/(dx*dx)) * (uc(1) - 2.0d0*uc(2) + uc(3))

!Bug here{
!  beta2c = -0.025d0
!  betanm1c = 0.09d0
!}

  dc(2,1) = (1.0d0/dx) * (((2.0d0*beta2c/3.0d0)-1.0d0/3.0d0)*uc(1)  &
       - ((8.0d0*beta2c/3.0d0)+1.0d0/2.0d0)*uc(2)              &
       + (4.0d0*beta2c + 1.0d0)*uc(3)                          &
       - ((8.0d0*beta2c/3.0d0)+1.0d0/6.0d0)*uc(4)              &
       + (2.0d0*beta2c/3.0d0)*uc(5))
  dc(2,2) = (1.0d0/(dx*dx)) * (uc(1) - 2.0d0*uc(2) + uc(3))
else
  dc(1,1) = ubar(1,1)
  dc(1,2) = ubar(1,2)

  dc(2,1) = (15.0d0/(16.0d0*dx)) * (uc(3)-uc(1))
endif

if(rank.eq.nprocs-1)then
  dc(n,1) = -(0.50d0/dx) * (-3.0d0 * uc(n) + 4.0d0 * uc(n-1) - uc(n-2))
  dc(n,2) = (1.0d0/(dx*dx)) * (uc(n) - 2.0d0*uc(n-1) + uc(n-2))

  dc(n-1,1) = (-1.0d0/dx) * (((2.0d0*betanm1c/3.0d0)-1.0d0/3.0d0)*uc(n)  &
       - ((8.0d0*betanm1c/3.0d0)+1.0d0/2.0d0)*uc(n-1)               &
       + (4.0d0*betanm1c + 1.0d0)*uc(n-2)                           &
       - ((8.0d0*betanm1c/3.0d0)+1.0d0/6.0d0)*uc(n-3)               &
       + (2.0d0*betanm1c/3.0d0)*uc(n-4))

  dc(n-1,2) = (1.0d0/(dx*dx)) * (uc(n) - 2.0d0*uc(n-1) + uc(n-2))
else
  dc(n,1) = ubar(n,1)
  dc(n,2) = ubar(n,2)

     dc(n-1,1) = (15.0d0/(16.0d0*dx)) * (uc(n)-uc(n-2))
endif

  do i=3,n-2
     dc(i,1) = (15.0d0/(16.0d0*dx)) * (uc(i+1)-uc(i-1))
  end do
  do i=2,n-1
     dc(i,2) = (3.0d0/(dx*dx)) * (uc(i+1)- 2.0d0*uc(i) + uc(i-1))
  end do

  do j=1,2
     do k=1,2
        betac(1,j,k) = bc(1,j,k)
     end do
  end do

  do i=2,n
     det = (betac(i-1,2,2)*betac(i-1,1,1) - betac(i-1,2,1)*betac(i-1,1,2))
     temp11 = (betac(i-1,2,2)*cc(i-1,1,1) - betac(i-1,1,2)*cc(i-1,2,1))/det
     temp12 = (betac(i-1,2,2)*cc(i-1,1,2) - betac(i-1,1,2)*cc(i-1,2,2))/det
     temp21 = (-betac(i-1,2,1)*cc(i-1,1,1) + betac(i-1,1,1)*cc(i-1,2,1))/det
     temp22 = (-betac(i-1,2,1)*cc(i-1,1,2) + betac(i-1,1,1)*cc(i-1,2,2))/det
     betac(i,1,1) = bc(i,1,1) - (ac(i,1,1)*temp11 + ac(i,1,2)*temp21)
     betac(i,1,2) = bc(i,1,2) - (ac(i,1,1)*temp12 + ac(i,1,2)*temp22)
     betac(i,2,1) = bc(i,2,1) - (ac(i,2,1)*temp11 + ac(i,2,2)*temp21)
     betac(i,2,2) = bc(i,2,2) - (ac(i,2,1)*temp12 + ac(i,2,2)*temp22)
  end do

  do i=2,n
     det = (betac(i-1,2,2)*betac(i-1,1,1) - betac(i-1,2,1)*betac(i-1,1,2))
     alphac(i,1,1) = (ac(i,1,1)*betac(i-1,2,2) - ac(i,1,2)*betac(i-1,2,1))/det
     alphac(i,1,2) = (-ac(i,1,1)*betac(i-1,1,2) + ac(i,1,2)*betac(i-1,1,1))/det
     alphac(i,2,1) = (ac(i,2,1)*betac(i-1,2,2) - ac(i,2,2)*betac(i-1,2,1))/det
     alphac(i,2,2) = (-ac(i,2,1)*betac(i-1,1,2) + ac(i,2,2)*betac(i-1,1,1))/det
  end do

  wc(1,1) = dc(1,1)
  wc(1,2) = dc(1,2)
  do i=2,n
     wc(i,1) = dc(i,1) - (alphac(i,1,1)*wc(i-1,1)+alphac(i,1,2)*wc(i-1,2))
     wc(i,2) = dc(i,2) - (alphac(i,2,1)*wc(i-1,1)+alphac(i,2,2)*wc(i-1,2))
  end do

  det =(betac(n,2,2)*betac(n,1,1) - betac(n,2,1)*betac(n,1,2))
  ubar(n,1) = (betac(n,2,2)*wc(n,1) - betac(n,1,2)*wc(n,2))/det
  ubar(n,2) = (-betac(n,2,1)*wc(n,1) + betac(n,1,1)*wc(n,2))/det

  do i=n-1,1,-1
     temp11 = wc(i,1) - (cc(i,1,1)*ubar(i+1,1) + cc(i,1,2)*ubar(i+1,2))
     temp21 = wc(i,2) - (cc(i,2,1)*ubar(i+1,1) + cc(i,2,2)*ubar(i+1,2))
     det=(betac(i,2,2)*betac(i,1,1)-betac(i,2,1)*betac(i,1,2))
     ubar(i,1) = (betac(i,2,2)*temp11 - betac(i,1,2)*temp21)/det
     ubar(i,2) = (-betac(i,2,1)*temp11 + betac(i,1,1)*temp21)/det
  end do

     dvort_dx(1:n,jloop) = ubar(1:n,1)
     d2vort_dx2(1:n,jloop) = ubar(1:n,2)

!  do i=ist,ien
!   dvort_dx(i,jloop) = (vort(i+1,jloop)-vort(i-1,jloop))/(2.0d0*dx)!dvort_dx(1,jloop)
!   d2vort_dx2(i,jloop) = (vort(i+1,jloop)-2.0d0*vort(i,jloop)+vort(i-1,jloop))/(dx**2)!d2vort_dx2(1,jloop)
!  end do

!check if 1 and n match with exact der
end do

end subroutine ccd_x

subroutine ccd_y
  use VAR
  implicit none
  integer :: i,j,k,iloop
  double precision :: ac(m,2,2),bc(m,2,2),cc(m,2,2),dc(m,2),alphac(m,2,2),det
  double precision :: betac(m,2,2),uc(m),ubar(m,2),beta2c,betanm1c
  double precision :: temp11, temp12, temp21, temp22, wc(m,2)
  double precision :: max_first=0.0d0, max_second=0.0d0

     do iloop=1,n
        uc(1:m) = vort(iloop,1:m)

  do i=1,m
     bc(i,1,1) = 1.0d0
     bc(i,2,2) = 1.0d0
     bc(i,1,2) = 0.0d0
     bc(i,2,1) = 0.0d0
  end do

  ac(1:2,:,:) = 0.0d0
  cc(1:2,:,:) = 0.0d0
  cc(2,2,1) = 9.0d0/(8.0d0*dy)
  cc(2,2,2) = -1.0d0/8.0d0
  ac(2,2,1) = -9.0d0/(8.0d0*dy)
  ac(2,2,2) = -1.0d0/8.0d0

  ac(m-1:m,:,:) = 0.0d0
  cc(m-1:m,:,:) = 0.0d0
  cc(m-1,2,1) = 9.0d0/(8.0d0*dy)
  cc(m-1,2,2) = -1.0d0/8.0d0
  ac(m-1,2,1) = -9.0d0/(8.0d0*dy)
  ac(m-1,2,2) = -1.0d0/8.0d0

  do i=3,m-2
     ac(i,1,1) = 7.0d0 / 16.0d0
     ac(i,2,2) = -1.0d0/8.0d0
     ac(i,1,2) = dy/16.0d0
     ac(i,2,1) = -9.0d0/(8.0d0*dy)
     cc(i,1,1) = 7.0d0 / 16.0d0
     cc(i,2,2) = -1.0d0/8.0d0
     cc(i,1,2) = -dy/16.0d0
     cc(i,2,1) = 9.0d0/(8.0d0*dy)
  end do

  dc(1,1) = (0.50d0/dy) * (-3.0d0 * uc(1) + 4.0d0 * uc(2) - uc(3))
  dc(1,2) = (1.0d0/(dy*dy)) * (uc(1) - 2.0d0*uc(2) + uc(3))
  dc(m,1) = -(0.50d0/dy) * (-3.0d0 * uc(m) + 4.0d0 * uc(m-1) - uc(m-2))
  dc(m,2) = (1.0d0/(dy*dy)) * (uc(m) - 2.0d0*uc(m-1) + uc(m-2))

  beta2c = -0.025d0
  betanm1c = 0.09d0
  dc(2,1) = (1.0d0/dy) * (((2.0d0*beta2c/3.0d0)-1.0d0/3.0d0)*uc(1)  &
       - ((8.0d0*beta2c/3.0d0)+1.0d0/2.0d0)*uc(2)              &
       + (4.0d0*beta2c + 1.0d0)*uc(3)                          &
       - ((8.0d0*beta2c/3.0d0)+1.0d0/6.0d0)*uc(4)              &
       + (2.0d0*beta2c/3.0d0)*uc(5))
  dc(m-1,1) = (-1.0d0/dy) * (((2.0d0*betanm1c/3.0d0)-1.0d0/3.0d0)*uc(m)  &
       - ((8.0d0*betanm1c/3.0d0)+1.0d0/2.0d0)*uc(m-1)               &
       + (4.0d0*betanm1c + 1.0d0)*uc(m-2)                           &
       - ((8.0d0*betanm1c/3.0d0)+1.0d0/6.0d0)*uc(m-3)               &
       + (2.0d0*betanm1c/3.0d0)*uc(m-4))

  dc(2,2) = (1.0d0/(dy*dy)) * (uc(1) - 2.0d0*uc(2) + uc(3))
  dc(m-1,2) = (1.0d0/(dy*dy)) * (uc(m) - 2.0d0*uc(m-1) + uc(m-2))

  do i=3,m-2
     dc(i,1) = (15.0d0/(16.0d0*dy)) * (uc(i+1)-uc(i-1))
  end do
  do i=2,m-1
     dc(i,2) = (3.0d0/(dy*dy)) * (uc(i+1)- 2.0d0*uc(i) + uc(i-1))
  end do

  do j=1,2
     do k=1,2
        betac(1,j,k) = bc(1,j,k)
     end do
  end do

  do i=2,m
     det = (betac(i-1,2,2)*betac(i-1,1,1) - betac(i-1,2,1)*betac(i-1,1,2))
     temp11 = (betac(i-1,2,2)*cc(i-1,1,1) - betac(i-1,1,2)*cc(i-1,2,1))/det
     temp12 = (betac(i-1,2,2)*cc(i-1,1,2) - betac(i-1,1,2)*cc(i-1,2,2))/det
     temp21 = (-betac(i-1,2,1)*cc(i-1,1,1) + betac(i-1,1,1)*cc(i-1,2,1))/det
     temp22 = (-betac(i-1,2,1)*cc(i-1,1,2) + betac(i-1,1,1)*cc(i-1,2,2))/det
     betac(i,1,1) = bc(i,1,1) - (ac(i,1,1)*temp11 + ac(i,1,2)*temp21)
     betac(i,1,2) = bc(i,1,2) - (ac(i,1,1)*temp12 + ac(i,1,2)*temp22)
     betac(i,2,1) = bc(i,2,1) - (ac(i,2,1)*temp11 + ac(i,2,2)*temp21)
     betac(i,2,2) = bc(i,2,2) - (ac(i,2,1)*temp12 + ac(i,2,2)*temp22)
  end do

  do i=2,m
     det = (betac(i-1,2,2)*betac(i-1,1,1) - betac(i-1,2,1)*betac(i-1,1,2))
     alphac(i,1,1) = (ac(i,1,1)*betac(i-1,2,2) - ac(i,1,2)*betac(i-1,2,1))/det
     alphac(i,1,2) = (-ac(i,1,1)*betac(i-1,1,2) + ac(i,1,2)*betac(i-1,1,1))/det
     alphac(i,2,1) = (ac(i,2,1)*betac(i-1,2,2) - ac(i,2,2)*betac(i-1,2,1))/det
     alphac(i,2,2) = (-ac(i,2,1)*betac(i-1,1,2) + ac(i,2,2)*betac(i-1,1,1))/det
  end do

  wc(1,1) = dc(1,1)
  wc(1,2) = dc(1,2)
  do i=2,m
     wc(i,1) = dc(i,1) - (alphac(i,1,1)*wc(i-1,1)+alphac(i,1,2)*wc(i-1,2))
     wc(i,2) = dc(i,2) - (alphac(i,2,1)*wc(i-1,1)+alphac(i,2,2)*wc(i-1,2))
  end do

  det =(betac(m,2,2)*betac(m,1,1) - betac(m,2,1)*betac(m,1,2))
  ubar(m,1) = (betac(m,2,2)*wc(m,1) - betac(m,1,2)*wc(m,2))/det
  ubar(m,2) = (-betac(m,2,1)*wc(m,1) + betac(m,1,1)*wc(m,2))/det

  do i=m-1,1,-1
     temp11 = wc(i,1) - (cc(i,1,1)*ubar(i+1,1) + cc(i,1,2)*ubar(i+1,2))
     temp21 = wc(i,2) - (cc(i,2,1)*ubar(i+1,1) + cc(i,2,2)*ubar(i+1,2))
     det=(betac(i,2,2)*betac(i,1,1)-betac(i,2,1)*betac(i,1,2))
     ubar(i,1) = (betac(i,2,2)*temp11 - betac(i,1,2)*temp21)/det
     ubar(i,2) = (-betac(i,2,1)*temp11 + betac(i,1,1)*temp21)/det
  end do

     dvort_dy(iloop,1:m) = ubar(1:m,1)
     d2vort_dy2(iloop,1:m) = ubar(1:m,2)
end do

end subroutine ccd_y

subroutine global_ccd_x
  use VAR
  implicit none
  integer :: i,j,k,jloop,block_col,i_interface
  double precision :: ac(nglo,2,2),bc(nglo,2,2),cc(nglo,2,2),dc(nglo,2),alphac(nglo,2,2),det
  double precision :: betac(nglo,2,2),uc(nglo),ubar(nglo,2),beta2c,betanm1c
  double precision :: temp11, temp12, temp21, temp22, wc(nglo,2)
!  double precision :: inv_block_matrix(nglo,nglo,2,2)
  double precision, dimension(:,:,:,:), allocatable ::inv_block_matrix
  double precision :: temp_vec(2)
  double precision :: max_first=0.0d0, max_second=0.0d0
!  double precision :: D1_mat(nglo,nglo), D2_mat(nglo,nglo)

  allocate(inv_block_matrix(nglo,nglo,2,2))

  do i=1,nglo
     bc(i,1,1) = 1.0d0
     bc(i,2,2) = 1.0d0
     bc(i,1,2) = 0.0d0
     bc(i,2,1) = 0.0d0
  end do

  ac(1:2,:,:) = 0.0d0
  cc(1:2,:,:) = 0.0d0
  cc(2,2,1) = 9.0d0/8.0d0 !9.0d0/(8.0d0*dx)
  cc(2,2,2) = -1.0d0/8.0d0
  ac(2,2,1) = -9.0d0/8.0d0!-9.0d0/(8.0d0*dx)
  ac(2,2,2) = -1.0d0/8.0d0

  ac(nglo-1:nglo,:,:) = 0.0d0
  cc(nglo-1:nglo,:,:) = 0.0d0

  cc(nglo-1,2,1) = 9.0d0/8.0d0!9.0d0/(8.0d0*dx)
  cc(nglo-1,2,2) = -1.0d0/8.0d0
  ac(nglo-1,2,1) = -9.0d0/8.0d0!-9.0d0/(8.0d0*dx)
  ac(nglo-1,2,2) = -1.0d0/8.0d0

  do i=3,nglo-2
     ac(i,1,1) = 7.0d0 / 16.0d0
     ac(i,2,2) = -1.0d0/8.0d0
     ac(i,1,2) = 1.0d0/16.0d0!dx/16.0d0
     ac(i,2,1) = -9.0d0/8.0d0!-9.0d0/(8.0d0*dx)
     cc(i,1,1) = 7.0d0 / 16.0d0
     cc(i,2,2) = -1.0d0/8.0d0
     cc(i,1,2) = -1.0d0/16.0d0!-dx/16.0d0
     cc(i,2,1) = 9.0d0/8.0d0!9.0d0/(8.0d0*dx)
  end do

  do jloop=1,nglo

  do block_col = 1,2

  do i = 1,nglo
  dc(i,1) = 0.0d0
  dc(i,2) = 0.0d0
  end do
  dc(jloop,block_col) = 1.0d0

  do j=1,2
     do k=1,2
        betac(1,j,k) = bc(1,j,k)
     end do
  end do

  do i=2,nglo
     det = (betac(i-1,2,2)*betac(i-1,1,1) - betac(i-1,2,1)*betac(i-1,1,2))
     temp11 = (betac(i-1,2,2)*cc(i-1,1,1) - betac(i-1,1,2)*cc(i-1,2,1))/det
     temp12 = (betac(i-1,2,2)*cc(i-1,1,2) - betac(i-1,1,2)*cc(i-1,2,2))/det
     temp21 = (-betac(i-1,2,1)*cc(i-1,1,1) + betac(i-1,1,1)*cc(i-1,2,1))/det
     temp22 = (-betac(i-1,2,1)*cc(i-1,1,2) + betac(i-1,1,1)*cc(i-1,2,2))/det
     betac(i,1,1) = bc(i,1,1) - (ac(i,1,1)*temp11 + ac(i,1,2)*temp21)
     betac(i,1,2) = bc(i,1,2) - (ac(i,1,1)*temp12 + ac(i,1,2)*temp22)
     betac(i,2,1) = bc(i,2,1) - (ac(i,2,1)*temp11 + ac(i,2,2)*temp21)
     betac(i,2,2) = bc(i,2,2) - (ac(i,2,1)*temp12 + ac(i,2,2)*temp22)
  end do

  do i=2,nglo
     det = (betac(i-1,2,2)*betac(i-1,1,1) - betac(i-1,2,1)*betac(i-1,1,2))
     alphac(i,1,1) = (ac(i,1,1)*betac(i-1,2,2) - ac(i,1,2)*betac(i-1,2,1))/det
     alphac(i,1,2) = (-ac(i,1,1)*betac(i-1,1,2) + ac(i,1,2)*betac(i-1,1,1))/det
     alphac(i,2,1) = (ac(i,2,1)*betac(i-1,2,2) - ac(i,2,2)*betac(i-1,2,1))/det
     alphac(i,2,2) = (-ac(i,2,1)*betac(i-1,1,2) + ac(i,2,2)*betac(i-1,1,1))/det
  end do

  wc(1,1) = dc(1,1)
  wc(1,2) = dc(1,2)
  do i=2,nglo
     wc(i,1) = dc(i,1) - (alphac(i,1,1)*wc(i-1,1)+alphac(i,1,2)*wc(i-1,2))
     wc(i,2) = dc(i,2) - (alphac(i,2,1)*wc(i-1,1)+alphac(i,2,2)*wc(i-1,2))
  end do

  det =(betac(nglo,2,2)*betac(nglo,1,1) - betac(nglo,2,1)*betac(nglo,1,2))
  ubar(nglo,1) = (betac(nglo,2,2)*wc(nglo,1) - betac(nglo,1,2)*wc(nglo,2))/det
  ubar(nglo,2) = (-betac(nglo,2,1)*wc(nglo,1) + betac(nglo,1,1)*wc(nglo,2))/det
  inv_block_matrix(nglo,jloop,1,block_col) = ubar(nglo,1)
  inv_block_matrix(nglo,jloop,2,block_col) = ubar(nglo,2)

  do i=nglo-1,1,-1
     temp11 = wc(i,1) - (cc(i,1,1)*ubar(i+1,1) + cc(i,1,2)*ubar(i+1,2))
     temp21 = wc(i,2) - (cc(i,2,1)*ubar(i+1,1) + cc(i,2,2)*ubar(i+1,2))
     det=(betac(i,2,2)*betac(i,1,1)-betac(i,2,1)*betac(i,1,2))
     ubar(i,1) = (betac(i,2,2)*temp11 - betac(i,1,2)*temp21)/det
     ubar(i,2) = (-betac(i,2,1)*temp11 + betac(i,1,1)*temp21)/det

     inv_block_matrix(i,jloop,1,block_col) = ubar(i,1)!/dx
     inv_block_matrix(i,jloop,2,block_col) = ubar(i,2)!/(dx*dx)
  end do

end do
end do

do i_interface = 1,all_interfaceids_i_size
D_inv(i_interface,ist_wb:ien_wb,1,1) = inv_block_matrix(glo_interfaceids_i(i_interface),&
                                           ist_glo:ien_glo,1,1)
D_inv(i_interface,ist_wb:ien_wb,1,2) = inv_block_matrix(glo_interfaceids_i(i_interface),&
                                           ist_glo:ien_glo,1,2)
D_inv(i_interface,ist_wb:ien_wb,2,1) = inv_block_matrix(glo_interfaceids_i(i_interface),&
                                           ist_glo:ien_glo,2,1)
D_inv(i_interface,ist_wb:ien_wb,2,2) = inv_block_matrix(glo_interfaceids_i(i_interface),&
                                           ist_glo:ien_glo,2,2)
end do

  deallocate(inv_block_matrix)

!dvort_dx_g = 0.0d0
!d2vort_dx2_g = 0.0d0
!max_first = 0.0d0
!max_second = 0.0d0
!do jloop = 1,mglo
!
!  uc(1:nglo) = vortg(1:nglo,jloop)
!
!  dc(1,1) = (0.50d0/dx) * (-3.0d0 * uc(1) + 4.0d0 * uc(2) - uc(3))
!  dc(1,2) = (1.0d0/(dx*dx)) * (uc(1) - 2.0d0*uc(2) + uc(3))
!  dc(nglo,1) = -(0.50d0/dx) * (-3.0d0 * uc(nglo) + 4.0d0 * uc(nglo-1) - uc(nglo-2))
!  dc(nglo,2) = (1.0d0/(dx*dx)) * (uc(nglo) - 2.0d0*uc(nglo-1) + uc(nglo-2))
!
!  beta2c = -0.025d0
!  betanm1c = 0.09d0
!  dc(2,1) = (1.0d0/dx) * (((2.0d0*beta2c/3.0d0)-1.0d0/3.0d0)*uc(1)  &
!       - ((8.0d0*beta2c/3.0d0)+1.0d0/2.0d0)*uc(2)              &
!       + (4.0d0*beta2c + 1.0d0)*uc(3)                          &
!       - ((8.0d0*beta2c/3.0d0)+1.0d0/6.0d0)*uc(4)              &
!       + (2.0d0*beta2c/3.0d0)*uc(5))
!  dc(nglo-1,1) = (-1.0d0/dx) * (((2.0d0*betanm1c/3.0d0)-1.0d0/3.0d0)*uc(nglo)  &
!       - ((8.0d0*betanm1c/3.0d0)+1.0d0/2.0d0)*uc(nglo-1)               &
!       + (4.0d0*betanm1c + 1.0d0)*uc(nglo-2)                           &
!       - ((8.0d0*betanm1c/3.0d0)+1.0d0/6.0d0)*uc(nglo-3)               &
!       + (2.0d0*betanm1c/3.0d0)*uc(nglo-4))
!
!  dc(2,2) = (1.0d0/(dx*dx)) * (uc(1) - 2.0d0*uc(2) + uc(3))
!  dc(nglo-1,2) = (1.0d0/(dx*dx)) * (uc(nglo) - 2.0d0*uc(nglo-1) + uc(nglo-2))
!
!  do i=3,nglo-2
!     dc(i,1) = (15.0d0/(16.0d0*dx)) * (uc(i+1)-uc(i-1))
!  end do
!  do i=2,nglo-1
!     dc(i,2) = (3.0d0/(dx*dx)) * (uc(i+1)- 2.0d0*uc(i) + uc(i-1))
!  end do
!
!do i=1,nglo
!do j=1,nglo
!!inv_block_matrix(i,j,1:2,1:2)xdc(j)
!temp_vec(1)=inv_block_matrix(i,j,1,1)*dc(j,1) + inv_block_matrix(i,j,1,2)*dc(j,2)
!temp_vec(2)=inv_block_matrix(i,j,2,1)*dc(j,1) + inv_block_matrix(i,j,2,2)*dc(j,2)
!dvort_dx_g(i,jloop) = dvort_dx_g(i,jloop) + temp_vec(1)
!d2vort_dx2_g(i,jloop) = d2vort_dx2_g(i,jloop) + temp_vec(2)
!end do
!if(max_first .lt. dabs(dvort_dx_g(i,jloop)))max_first = dabs(dvort_dx_g(i,jloop))
!if(max_second .lt. dabs(d2vort_dx2_g(i,jloop)))max_second = dabs(d2vort_dx2_g(i,jloop))
!end do
!end do

!!write(*,*) max_first, max_second

end subroutine global_ccd_x

subroutine CALC_EX_DER_AT_INTERFACE
use VAR
implicit none
integer :: i,j,i_interface,der_count = 1
double precision :: partial_sum_for_der(2*all_interfaceids_i_size*m)
double precision :: glo_sum_for_der(2*all_interfaceids_i_size*m)
double precision :: uc(ist_wb-1:ien_wb+1),dc(ist_wb:ien_wb,2),temp_vec(2)
double precision, parameter :: beta2c=-0.025d0,betanm1c=0.09d0

partial_sum_for_der = 0.0d0
glo_sum_for_der = 0.0d0
der_count = 1

!First calculate partial sum{
do j=jst_wb,jen_wb
  uc(ist_wb-1:ien_wb+1) = vort(ist_wb-1:ien_wb+1,j)

if(rank.eq.0)then
  dc(1,1) = (0.50d0/dx) * (-3.0d0 * uc(1) + 4.0d0 * uc(2) - uc(3))
  dc(1,2) = (1.0d0/(dx)) * (uc(1) - 2.0d0*uc(2) + uc(3))!(1.0d0/(dx*dx)) * (uc(1) - 2.0d0*uc(2) + uc(3))

!{Bug here
!  beta2c = -0.025d0
!  betanm1c = 0.09d0
!}

  dc(2,1) = (1.0d0/dx) * (((2.0d0*beta2c/3.0d0)-1.0d0/3.0d0)*uc(1)  &
       - ((8.0d0*beta2c/3.0d0)+1.0d0/2.0d0)*uc(2)              &
       + (4.0d0*beta2c + 1.0d0)*uc(3)                          &
       - ((8.0d0*beta2c/3.0d0)+1.0d0/6.0d0)*uc(4)              &
       + (2.0d0*beta2c/3.0d0)*uc(5))
  dc(2,2) = (1.0d0/(dx)) * (uc(1) - 2.0d0*uc(2) + uc(3))!(1.0d0/(dx*dx)) * (uc(1) - 2.0d0*uc(2) + uc(3))
else
  dc(1,1) = (15.0d0/(16.0d0*dx)) * (uc(2)-uc(0))
  dc(1,2) = (3.0d0/(dx)) * (uc(2)- 2.0d0*uc(1) + uc(0))!(3.0d0/(dx*dx)) * (uc(2)- 2.0d0*uc(1) + uc(0))

  dc(2,1) = (15.0d0/(16.0d0*dx)) * (uc(3)-uc(1))
  dc(2,2) = (3.0d0/(dx)) * (uc(3)- 2.0d0*uc(2) + uc(1))!(3.0d0/(dx*dx)) * (uc(3)- 2.0d0*uc(2) + uc(1))
endif

if(rank.eq.nprocs-1)then
  dc(n,1) = -(0.50d0/dx) * (-3.0d0 * uc(n) + 4.0d0 * uc(n-1) - uc(n-2))
  dc(n,2) = (1.0d0/(dx)) * (uc(n) - 2.0d0*uc(n-1) + uc(n-2))!(1.0d0/(dx*dx)) * (uc(n) - 2.0d0*uc(n-1) + uc(n-2))

  dc(n-1,1) = (-1.0d0/dx) * (((2.0d0*betanm1c/3.0d0)-1.0d0/3.0d0)*uc(n)  &
       - ((8.0d0*betanm1c/3.0d0)+1.0d0/2.0d0)*uc(n-1)               &
       + (4.0d0*betanm1c + 1.0d0)*uc(n-2)                           &
       - ((8.0d0*betanm1c/3.0d0)+1.0d0/6.0d0)*uc(n-3)               &
       + (2.0d0*betanm1c/3.0d0)*uc(n-4))

  dc(n-1,2) = (1.0d0/(dx)) * (uc(n) - 2.0d0*uc(n-1) + uc(n-2))!(1.0d0/(dx*dx)) * (uc(n) - 2.0d0*uc(n-1) + uc(n-2))
else
  dc(n,1) = (15.0d0/(16.0d0*dx)) * (uc(n+1)-uc(n-1))
  dc(n,2) = (3.0d0/(dx)) * (uc(n+1)- 2.0d0*uc(n) + uc(n-1))!(3.0d0/(dx*dx)) * (uc(n+1)- 2.0d0*uc(n) + uc(n-1))

  dc(n-1,1) = (15.0d0/(16.0d0*dx)) * (uc(n)-uc(n-2))
  dc(n-1,2) = (3.0d0/(dx)) * (uc(n)- 2.0d0*uc(n-1) + uc(n-2))!(3.0d0/(dx*dx)) * (uc(n)- 2.0d0*uc(n-1) + uc(n-2))
endif

  do i=3,n-2
     dc(i,1) = (15.0d0/(16.0d0*dx)) * (uc(i+1)-uc(i-1))
  end do
  do i=2,n-1
     dc(i,2) = (3.0d0/(dx)) * (uc(i+1)- 2.0d0*uc(i) + uc(i-1))!(3.0d0/(dx*dx)) * (uc(i+1)- 2.0d0*uc(i) + uc(i-1))
  end do


do i_interface = 1,all_interfaceids_i_size
do i=ist_wb,ien_wb
temp_vec(1)=D_inv(i_interface,i,1,1)*dc(i,1) + D_inv(i_interface,i,1,2)*dc(i,2)
temp_vec(2)=D_inv(i_interface,i,2,1)*dc(i,1) + D_inv(i_interface,i,2,2)*dc(i,2)

partial_sum_for_der(der_count) = partial_sum_for_der(der_count) + &
                                 temp_vec(1)
partial_sum_for_der(der_count+1) = partial_sum_for_der(der_count+1) + &
                                 temp_vec(2)
end do
der_count = der_count+2
end do

end do
!}

!GLOBAL SUM{
     call GLOBAL_REDUCE (2*all_interfaceids_i_size*m,partial_sum_for_der,glo_sum_for_der,2)
!}

!Store the derivative at the interface{
der_count = 1
do j=jst_wb,jen_wb
do i_interface = 1,all_interfaceids_i_size

if(glo_interfaceids_i(i_interface).eq.ist_glo)then
dvort_dx(1,j) = glo_sum_for_der(der_count)
d2vort_dx2(1,j) = glo_sum_for_der(der_count+1)/dx
end if

if(glo_interfaceids_i(i_interface).eq.ien_glo)then
dvort_dx(n,j) = glo_sum_for_der(der_count)
d2vort_dx2(n,j) = glo_sum_for_der(der_count+1)/dx
end if

der_count = der_count+2
end do
end do
!}

!if(rank.eq.0)then
!do j=1,m
!write(*,*)dvort_dx(n,j),d2vort_dx2(n,j)
!end do
!end if
end subroutine CALC_EX_DER_AT_INTERFACE
!--------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------!

!MPI Routines are called here{
!red_opn = Reduction operation
subroutine GLOBAL_REDUCE(nitems, loc_data, glo_data, red_opn)
use VAR
implicit none
include 'mpif.h'
integer :: ierr
integer,intent(in) :: nitems, red_opn
double precision,intent(in)::loc_data(nitems)
double precision,intent(inout)::glo_data(nitems)

if(red_opn.eq.1)then
call MPI_ALLREDUCE(loc_data,glo_data,nitems,MPI_double_precision,MPI_MAX,MPI_COMM_WORLD,ierr)
else if (red_opn.eq.2)then
call MPI_ALLREDUCE(loc_data,glo_data,nitems,MPI_double_precision,MPI_SUM,MPI_COMM_WORLD,ierr)
else
write(*,*)"Reduction Operation not defined"
end if
end subroutine GLOBAL_REDUCE


!}

!--------------------------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------------------------!

subroutine write_global_solution_file
  use VAR
  implicit none
  include 'mpif.h'
  INTEGER :: ifh
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: istatus
!  double precision :: subarray(n,m)!(12,8) 
  integer :: ifiletype, ierror 
  integer :: sizes(2), subsizes(2), starts(2) 
  integer :: l
  character(11) :: str, string

  l=time
  if(l/10 .eq. 0) then
     write(str,'(F11.6)') time
     string = "000"//str(4:11)
  elseif(l/10 .gt. 0 .and. l/10 .lt. 10) then
     write(str,'(F11.6)') time
     string = "00"//str(3:11)
  elseif(l/10 .gt. 9 .and. l/10 .lt. 100) then
     write(str,'(F11.6)') time
     string = "0"//str(2:11)
  else
     write(str,'(F11.6)') time
     string = str(1:11)
  endif

sizes(1)=mglo 
sizes(2)=nglo 
subsizes(1)=m 
subsizes(2)=n 
starts(1)=0 
starts(2)=rank*subsizes(2)

call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, & 
             MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, & 
            ifiletype, ierror) 

call MPI_TYPE_COMMIT(ifiletype,ierror) 

call MPI_FILE_OPEN(MPI_COMM_WORLD, 'svt_'//trim(string),MPI_MODE_CREATE+MPI_MODE_RDWR, MPI_INFO_NULL,ifh,ierror)
call MPI_FILE_SET_VIEW(ifh,0, MPI_DOUBLE_PRECISION,ifiletype, "native", MPI_INFO_NULL,ierror)
call MPI_FILE_WRITE_ALL(ifh,vort,subsizes(1)*subsizes(2), MPI_DOUBLE_PRECISION,istatus,ierror)
call MPI_FILE_CLOSE(ifh,ierror) 

end subroutine write_global_solution_file
