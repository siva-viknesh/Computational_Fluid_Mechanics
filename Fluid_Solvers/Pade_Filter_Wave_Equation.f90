!***************** FILTER EULER CD2- 1D CONVECTION EQUATION *******************!
MODULE NUMERICAL_SCHEME
    IMPLICIT NONE
    CONTAINS

    SUBROUTINE EULER_SCHEME (Nx, Nc, dx, u)
        IMPLICIT NONE
        INTEGER,                           INTENT (IN)  ::  Nx        
        REAL (KIND = 8),                   INTENT (IN)  ::  Nc, dx
        REAL (KIND = 8), DIMENSION (1:Nx), INTENT (OUT) ::  u
        REAL (KIND = 8), DIMENSION (1:Nx)               ::  du_dx
        INTEGER                                         ::  i

        CALL CD2_SCHEME (Nx, u, dx, du_dx)
            
        u  = u - Nc * du_dx / 2.0d0 

    END SUBROUTINE EULER_SCHEME

    SUBROUTINE CD2_SCHEME (Nx, u, dx, du_dx)
        IMPLICIT NONE
        INTEGER,                           INTENT (IN)  ::  Nx
        REAL (KIND = 8), DIMENSION (1:Nx), INTENT (IN)  ::  u  
        REAL (KIND = 8),                   INTENT (IN)  ::  dx     
        REAL (KIND = 8), DIMENSION (1:Nx), INTENT (OUT) ::  du_dx
        INTEGER                                         ::  i

        DO i = 2, Nx -1
            du_dx(i) = ( u (i+1) - u(i-1) ) / (2.0d0*dx)
        ENDDO
        
        du_dx(1)  = ( u (2) - u(Nx) )  / (2.0d0*dx)
        du_dx(Nx) = ( u (1) - u(Nx-1) ) / (2.0d0*dx)

    END SUBROUTINE CD2_SCHEME

    SUBROUTINE CENTRAL_1D_FILTER (u,Nx,k0x ,Nc,alpha)
        IMPLICIT NONE
        INTEGER,                           INTENT (IN)     ::  Nx
        REAL (KIND = 8), DIMENSION (1:Nx), INTENT (INOUT)  ::  u
        REAL (KIND = 8),                   INTENT (IN)     ::  alpha, k0x, Nc
        REAL (KIND = 8)                                    ::  a0, a1
        REAL (KIND = 8), DIMENSION (1:Nx)                  ::  RHS 
        REAL (KIND = 8), DIMENSION (3,1:Nx)                ::  a
        INTEGER                                            ::  i

        a0 =   ( (1.0d0 + 2.0d0* alpha *COS (k0x)) / SQRT (1.0d0 + Nc **2 * &
            & SIN (k0x)**2) -1.0d0 - 2.0d0*alpha ) / COS (k0x)
        a1 = 1.0d0 + 2.0d0*alpha- a0

        DO i = 2, Nx -1
            RHS(i) = a0 * u (i) + a1 * 0.50d0 * (u(i-1) + u(i+1)) 
        ENDDO 
        
        RHS(1)  = a0 * u (1)  + a1 * 0.50d0 * (  u(Nx)  + u(2) ) 
        RHS(Nx) = a0 * u (Nx) + a1 * 0.50d0 * ( u(Nx-1) + u(1) )

        a(1,:) = alpha
        a(2,:) = 1.0d0
        a(3,:) = alpha 

        CALL PERIODIC_TDMA (a,RHS,u)

    END SUBROUTINE CENTRAL_1D_FILTER

    SUBROUTINE TDMA(a,b,c,d,x,n)
        implicit none
     !    a - sub-diagonal 
     !    b - the main diagonal
     !    c - sup-diagonal
     !    d - right part
     !    x - the answer
     !    n - number of equations

        INTEGER,                            INTENT (IN)    ::  n
        REAL (KIND = 8), DIMENSION (1:n),   INTENT (IN)    ::  a,b,c,d
        REAL (KIND = 8), DIMENSION (1:n),   INTENT (OUT)   ::  x
        REAL (KIND = 8), DIMENSION (1:n)                   ::  cp,dp
        REAL (KIND = 8)                                    ::  m
        INTEGER                                            ::  i
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)

        do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
        enddo
         
        x(n) = dp(n)

        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do
    END SUBROUTINE TDMA

    SUBROUTINE PERIODIC_TDMA (a,b,c)
        IMPLICIT NONE
        REAL    (kind=8), DIMENSION(:,:), INTENT(IN)                   :: a
        REAL    (kind=8), DIMENSION(  :), INTENT(IN)                   :: b
        REAL    (kind=8), DIMENSION(1:SIZE(b, dim=1)),INTENT (OUT)     :: c
        REAL    (kind=8), DIMENSION(1:SIZE(b, dim=1)  )                :: v
        REAL    (kind=8), DIMENSION(1:SIZE(b, dim=1)  )                :: q
        REAL    (kind=8), DIMENSION(2:SIZE(b, dim=1)  )                :: p
        REAL    (kind=8), DIMENSION(1:SIZE(b, dim=1)-1)                :: u
        REAL    (kind=8), DIMENSION(1:SIZE(b, dim=1)-2)                :: r,w
        REAL    (kind=8)           :: sum
        INTEGER            :: n
        INTEGER            :: i
        n = SIZE(b, dim=1)

        q(1) = a(2,1)
        w(1) = a(1,1)/q(1)
        u(1) = a(3,1)/q(1)
        r(1) = a(3,n)

        DO i=2,n-1
            p(i) = a(1,i)
            q(i) = a(2,i) - p(i)*u(i-1)
            IF(i .ne. n-1) THEN
                u(i) = a(3,i)/q(i)
                w(i) = -(p(i  )*w(i-1)/q(i))
                r(i) = -(r(i-1)*u(i-1))
            ELSE
                u(i) = (a(3,i)-p(i)*w(i-1))/q(i)
            END IF
        
        END DO

        p(n) = a(1,n) - (r(n-2)*u(n-2))
        sum = 0.0d0
        DO i=1,n-2
            sum = sum + r(i)*w(i)
        END DO
        q(n) = a(2,n) - p(n)*u(n-1) - sum
        
        !-------solving Lower matrix------------!
        
            v(1) = b(1)/q(1)
        
            DO i=2,n-1
                v(i) = (b(i)-p(i)*v(i-1))/q(i)
            END DO
        
            sum=0.0d0
            DO i=1,n-2
                sum = sum + r(i)*v(i)
            END DO
        
            v(n) = (b(n) - sum - p(n)*v(n-1))/q(n)  
                              
        !-------solving upper matrix------------!
        
            c(n  ) = v(n  )
            c(n-1) = v(n-1) - c(n)*u(n-1)
        
            DO i=n-2,1,-1
                c(i) = v(i) - w(i)*c(n) - u(i)*c(i+1)
            END DO
    END SUBROUTINE PERIODIC_TDMA

    SUBROUTINE WRITE_OUTPUT (Nx,T,u, X)
        IMPLICIT NONE
        CHARACTER                                       ::  FILE_NAME*120
        INTEGER,                           INTENT (IN)  ::  Nx
        REAL (KIND = 8), DIMENSION (1:Nx), INTENT (IN)  ::  u, X
        REAL (KIND = 8),                   INTENT (IN)  ::  T
        INTEGER                                         ::  i
        CHARACTER                                       ::  STR1*10 

        WRITE(STR1(1:4),'(i4.4)') int(T)
        WRITE(STR1(5:10),'(f6.5)') T - int(T)
        FILE_NAME = "TIME_SERIES_"//STR1//".dat"
        OPEN  (25,file = trim (FILE_NAME))
        WRITE (25,*) "VARIABLES=x,u"
        WRITE (25,*) "ZONE I = ", Nx-1
        WRITE (25,*) "SOLUTIONTIME=", T

        DO i = 1, Nx
            WRITE (25,*) X(i), u(i)
        ENDDO

        CLOSE (25)

    END SUBROUTINE WRITE_OUTPUT

    SUBROUTINE INITIAL_SOLUTION (Nx, u, x, k0x, dx)
        IMPLICIT NONE
        INTEGER,                            INTENT (IN)  ::  Nx
        REAL (KIND = 8), DIMENSION (1:Nx),  INTENT (OUT) ::  u
        REAL (KIND = 8), DIMENSION (1:Nx),  INTENT (IN)  ::  x
        REAL (KIND = 8),                    INTENT (IN)  ::  k0x, dx 
        INTEGER                                          ::  i
        REAL (KIND = 8)                                  ::  k0 
        REAL (KIND = 8)                                  ::  T  = 0.0d0

        k0 = k0x / dx

        DO i = 1, Nx
            u(i) = ( 2.0d0 + COS (k0*X(i)) ) * EXP ( -6.9314718060d-3*X(i)**2)  
        ENDDO

        CALL WRITE_OUTPUT (Nx,T,u, X)

    END SUBROUTINE INITIAL_SOLUTION

END MODULE NUMERICAL_SCHEME
!***************************** PROGRAM STARTS HERE ****************************!
PROGRAM FILTER_STABLIZATION

    USE NUMERICAL_SCHEME

    IMPLICIT NONE
    INTEGER,         PARAMETER                      ::  Nx = 513
    REAL (KIND = 8), DIMENSION (1:Nx)               ::  u,  X
    REAL (KIND = 8)                                 ::  X0, Xmin, Xmax, dx
    REAL (KIND = 8)                                 ::  T0, Tmax, dt, T
    REAL (KIND = 8), PARAMETER                      ::  Nc = 0.60d0    ! CFL
    REAL (KIND = 8), PARAMETER                      ::  k0x   = 0.417825d0
    REAL (KIND = 8), PARAMETER                      ::  C = 0.10d0
    REAL (KIND = 8)                                 ::  alpha = 0.450d0
    INTEGER                                         ::  i, Nt = 1000

 !********************** GRID INITILISATION STARTS HERE ***********************!
    Xmin = -20.0d0
    Xmax =  20.0d0
    X0   =  0.0d0                       ! LOCATION OF WAVE PACKET
    T0   =  0.0d0
    Tmax =  20.0d0

    dx = (Xmax - Xmin) / REAL (Nx-1)
    dt = Nc * dx / C
    DO  i = 1, Nx
        X (i) = Xmin + (i-1) * dx
    ENDDO
 !*********************** GRID INITILISATION ENDS HERE ************************!

 !************************* COMPUTATION STARTS HERE ***************************!
    CALL INITIAL_SOLUTION (Nx, u, x, k0x, dx)
    DO i = 1, Nt
        T = i*dt                ! CURRENT SOLUTION TIME 
        CALL EULER_SCHEME         (Nx, Nc, dx, u)
        CALL CENTRAL_1D_FILTER    (u, Nx, k0x,Nc,alpha)
        CALL WRITE_OUTPUT         (Nx,T,u, X)
    ENDDO

!************************** COMPUTATION ENDS HERE *****************************!

END PROGRAM FILTER_STABLIZATION
