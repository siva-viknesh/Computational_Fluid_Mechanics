!********************************* HYPERBOLIC GRID GENERATION *************************************!
MODULE INPUT_OUPTUT
    implicit none
    CONTAINS

!********************************** READING THE AIRFOIL DATA **************************************!
    SUBROUTINE READ_FILE (X,Y,N)
        implicit none
        integer                                         ::  i
        integer,                        intent(in)      ::  N 
        real*8,     dimension(1:N),     intent(out)     ::  X, Y                    
        character                                       ::  File_Name*120

     !***************************** READING THE FILE STARTS HERE **********************************!
        File_Name = "CO-ORDINATES_NACA_0015.dat"
        OPEN (18, file = File_Name)
        
        READ(18,*)
        do i = 1,N
          READ (18,*) X(i), Y(i)  
        enddo

        CLOSE (18)

     !****************************** READING THE FILE ENDS HERE ***********************************!
    END SUBROUTINE READ_FILE

    SUBROUTINE WRITE_OUTPUT_AIRFOIL(N,N_ETA,X,Y)
        implicit none
        integer                                         :: N, N_ETA
        real*8, dimension(1: N),     intent(in)         :: X, Y
        integer                                         :: i
        character                                       :: File_Name*120

     !***************************** WRITING THE FILE STARTS HERE **********************************!

        File_Name = "GENERATED_GRID.dat"
        OPEN(18, file = File_Name)
        WRITE(18,*) "Variable = X,Y"
        WRITE(18,*) "ZONE I = ", N, "J = ", N_ETA, "F = POINT"

        do i = 1, N
            WRITE(18,*) X(i), Y(i)  
        enddo             
    
        CLOSE(18)

     !****************************** WRITING THE FILE ENDS HERE ***********************************!
    END SUBROUTINE WRITE_OUTPUT_AIRFOIL

    SUBROUTINE WRITE_OUTPUT(N,X,Y)
        implicit none
        integer                                         :: N
        real*8, dimension(1: N),     intent(in)         :: X, Y
        integer                                         :: i
        character                                       :: File_Name*120

     !***************************** WRITING THE FILE STARTS HERE **********************************!

        File_Name = "GENERATED_GRID.dat"
        OPEN(18, file = File_Name, access = "Append")

        do i = 1, N
            WRITE(18,*) X(i), Y(i)  
        enddo             
    

        CLOSE(18)

     !****************************** WRITING THE FILE ENDS HERE ***********************************!
    END SUBROUTINE WRITE_OUTPUT

    SUBROUTINE ETA_LINE(N_ETA,ETA_MIN,ETA_MAX,ETA,d_eta)
        implicit none
        integer,                    intent(in)          :: N_ETA
        real*8,                     intent(in)          :: ETA_MIN , ETA_MAX 
        real*8, dimension(1:N_ETA), intent(out)         :: ETA
        real*8,                     intent(out)         :: d_eta
        integer                                         :: i        
        ! ETA_MAX -> Radius of Enclosure 
     !******************************** ETA SPACING STARTS HERE ************************************!

        d_eta =  ( ETA_MAX - ETA_MIN ) / real (N_ETA-1)
        ETA(1) = 0.00d0

        do i = 2, N_ETA
            ETA(i) = ETA_MIN + (i-1) * d_eta
        enddo  
     

     !******************************** ETA SPACING ENDS HERE **************************************!

    END SUBROUTINE ETA_LINE

    SUBROUTINE SHI_LINE (N,d_shi,SHI)
        implicit none
        integer,                     intent(in)         ::  N
        real*8,                      intent(out)        ::  d_shi
        real*8, dimension (1:N),     intent(out)        ::  SHI         
        integer                                         ::  i
        real*8                                          ::  SHI_0 = 0.0d0, SHI_MAX = 1.0d0

     !******************************** SHI SPACING STARTS HERE ************************************!
     
     d_shi = ( SHI_MAX - SHI_0 ) / real (N-1)

     SHI (1) = SHI_0
     do i = 2, N
         SHI(i) = SHI_0 + (i-1) * d_shi
     enddo

     !********************************* SHI SPACING ENDS HERE *************************************!

    END SUBROUTINE SHI_LINE

    SUBROUTINE H1_FACTOR (N,d_shi, X, Y,H1)
        implicit none
        integer,                      intent(in)        ::  N
        real*8,                       intent(in)        ::  d_shi
        real*8, dimension(1:N),       intent(in)        ::  X, Y
        real*8, dimension(1:N),       intent(out)       ::  H1
        integer                                         ::  i


        do i = 1, N-1
            if (i == 1) then
             H1 (i) = dsqrt (  ((X(i+1) - X(N-1)) / (2.0d0*d_shi) ) **2 + &
              & ( (Y(i+1) - Y(N-1)) / (2.0d0*d_shi) ) **2)     
            else
             H1 (i) = dsqrt (  ((X(i+1) - X(i-1)) / (2.0d0*d_shi) ) **2 + &
              & ( (Y(i+1) - Y(i-1)) / (2.0d0*d_shi) ) **2)
            endif                
        enddo
        H1 (N) = H1 (1)

    END SUBROUTINE H1_FACTOR

    SUBROUTINE H2_FACTOR(N_ETA,H,ETA, H2)
        implicit none
        integer,                        intent(in)      ::  N_ETA
        real*8,                         intent(in)      ::  H
        real*8, dimension (1:N_ETA),    intent(in)      ::  ETA
        real*8, dimension (1:N_ETA),    intent(out)     ::  H2
        integer                                         ::  i
        real*8                                          ::  b = 2.50d0

     !******************************* ETA CLUSTERING STARTS HERE **********************************!
        do i = 1, N_ETA
            H2(i) =  -H*( 1- tanh(b*(1- ETA(i))) / tanh(b) )             
        enddo 
     !******************************* ETA CLUSTERING STARTS HERE **********************************!

    END SUBROUTINE H2_FACTOR

    SUBROUTINE FUNCTION_F (N,H1,H2,f)
        implicit none
        integer,                        intent(in)      ::  N
        real*8, dimension(1:N),         intent(in)      ::  H1
        real*8,                         intent(in)      ::  H2
        real*8, dimension(1:N),         intent(out)     ::  f
        integer                                         ::  i
        
        do i = 1, N-1
            f(i) = H2 / H1(i)
        enddo  
        f (N) = f (1)

    END SUBROUTINE FUNCTION_F

END MODULE INPUT_OUPTUT

MODULE DERIVATIVE_FINDER
    implicit none
    CONTAINS

    SUBROUTINE I_DERIVATIVE(N,Y,dx,dy_dx)
        implicit none
        integer,                        intent(in)      ::  N
        real*8, dimension(1:N),         intent(in)      ::  Y
        real*8,                         intent(in)      ::  dx
        real*8, dimension(1:N),         intent(out)     ::  dy_dx
        integer                                         ::  i

        do i = 1, N-1  

            if (i==1) then
                dy_dx (i) = ( Y(i+1) - Y(N-1))  / (2.0d0*dx)
            else   
                dy_dx (i) = (Y(i+1) - Y(i-1))   / (2.0d0*dx)     
            endif
        enddo  
       
        dy_dx (N) = dy_dx (1)

    END SUBROUTINE I_DERIVATIVE

END MODULE DERIVATIVE_FINDER

!*********************************** PROGRAMS STARTS HERE *****************************************!

PROGRAM HYPERBOLIC_GRID_GEN
    USE INPUT_OUPTUT
    USE DERIVATIVE_FINDER
    implicit none
    integer, parameter                                  ::  N = 2001, N_ETA = 351
    real*8                                              ::  ETA_MIN = 0.0d0, ETA_MAX = 1.0d0
    real*8,  dimension(1:N)                             ::  X, Y
    real*8,  dimension(1:N)                             ::  SHI, Y_SHI, X_SHI, H2
    real*8,  dimension(1:N_ETA)                         ::  ETA
    real*8,  dimension(1:N)                             ::  X_ETA, Y_ETA
    real*8,  dimension(1:N)                             ::  f, H1
    integer                                             ::  i, j
    real*8                                              ::  d_shi, d_eta, X_ETA_TEMP, H = 30.0d0
    ! N     -> No. of the Points on the  Airfoil 
    ! N_ETA -> No. of ETA Lines
 
 !************************************ COMPILATION STARTS HERE ************************************!

    CALL READ_FILE              (X,Y,N)
    CALL WRITE_OUTPUT_AIRFOIL   (N,N_ETA,X,Y)
    CALL SHI_LINE               (N,d_shi,SHI)
    CALL ETA_LINE               (N_ETA,ETA_MIN,ETA_MAX, ETA, d_eta)
    CALL H2_FACTOR              (N_ETA,H,ETA, H2)
    
    X_ETA = X
    Y_ETA = Y

    do i = 2, N_ETA
        CALL H1_FACTOR          (N,d_shi, X, Y,H1)
        CALL FUNCTION_F         (N, H1, H2(i), f)
        do j = 1, N-1
            CALL I_DERIVATIVE   (N,Y,d_shi,Y_SHI)
            CALL I_DERIVATIVE   (N,X,d_shi,X_SHI)
            X_ETA_TEMP  =   X_ETA(j)
            X_ETA(j)    =   X_ETA(j) + ( f(j) * Y_SHI(j) ) * d_eta
            Y_ETA(j)    =   Y_ETA(j) - X_SHI(j) / Y_SHI(j) * (X_ETA(j) - X_ETA_TEMP)
        enddo       
            Y_ETA (N)   =   Y_ETA (1)
            X_ETA (N)   =   X_ETA (1)

        CALL WRITE_OUTPUT       (N,X_ETA,Y_ETA)
    enddo    
    
END PROGRAM HYPERBOLIC_GRID_GEN
