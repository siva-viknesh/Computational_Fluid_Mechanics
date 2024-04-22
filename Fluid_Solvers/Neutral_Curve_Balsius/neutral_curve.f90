module global_variables
    implicit none

    ! Blasius profile 
        integer                       :: n = 12000
        double precision              :: etamax = 12.0d0
        double precision, allocatable :: y(:),U(:), U2(:)
        double precision              :: blstol = 1d-10, dely,deta
        logical                       :: write_balsius = .true.

        integer                       :: nproc, mRank, ist, ien , ierr1 

    ! OSE
        integer                       :: nre = 30001, glo_nw, pow1, pow2
        double precision              :: omega0, aMax, wMax, wi, ai 
        double complex                :: va0, van, vw0, vwn
        double complex                :: xi = dcmplx(0.0d0,1.0d0), one = dcmplx(1.0d0,0.0d0), zro = dcmplx(0.0d0,0.0d0)

        double precision              :: re0 = 520.0d0, ren = 15520.0d0 
        double precision, allocatable :: vw(:), vRe(:)
        double precision              :: bicgtol = 1.0d-8
        integer                       :: opt_integration = 2 , j_in   
        character*100                 :: outputfile, inputfile, outputfile2, zname
        character*800                 :: aux_info, outdir 
        character*3                  :: xstr        

        double complex, allocatable, dimension(:)   :: y1,y2,y3,y4,y5,y6
        double complex, allocatable, dimension(:)   :: phi0, phi1
        double precision, allocatable, dimension(:) :: y1vec
    end module global_variables


module mylib
    use global_variables 
    implicit none
    contains
    subroutine  linespace(x,x0,xn,lenx)
        implicit none
        double precision, dimension(:), intent(out) :: x
        double precision :: x0,xn,dx
        integer :: lenx,i 
        dx        = (xn - x0)/(lenx - 1) 
        x(1:lenx) = [( x0 + ((i-1)*dx), i=1,lenx ) ]       
        end subroutine linespace
    end module mylib

module Blasius
    use global_variables  
    use mylib 
    implicit none
    contains
    subroutine  blasius_profile
        implicit none
        integer :: i
        double precision,dimension(n) :: eta,f,f1,f2,f3
        double precision :: du1,du2,du3,du4
        double precision :: f21p,f1np,f21,f1n,delf21
        double precision  :: disp_thick
        
        call linespace(eta,0.0d0,etamax,n)

        !Newton Rapshon
            f21p  = 0.33205733621527d0
            f21   = 0.33205733621527d0 + 1.0d-4
            call blasius_integrator(f1np,f,f1,f2,f3,eta,f21p)
            call blasius_integrator(f1n,f,f1,f2,f3,eta,f21)
        
            do i=1,1000
                delf21 = (1.0d0-f1n)*(f21p - f21)/(f1np - f1n)
                f21p   = f21
                f21    = f21 + delf21
                f1np   = f1n
                call blasius_integrator(f1n,f,f1,f2,f3,eta,f21)
               
                if (abs(delf21) < blstol) exit
            enddo
        !Scaling (eta to displacement thickness conversion)
            disp_thick = eta(n) - f(n)
            y    = eta/disp_thick
            U    = f1
            U2   = f3*(disp_thick)**2
            dely = deta/disp_thick     
            if(write_balsius) then 
            open(31,file='blasius.dat')
            write(31,*) "Variables = y, U, U'' " 
            do i = 1,n 
                write(31,*) y(i), U(i), U2(i)
            enddo
            close(31)
            endif
        end subroutine  blasius_profile
    subroutine  blasius_integrator(f1n,f,f1,f2,f3,eta,f21)
        implicit none
        integer :: i
        double precision,dimension(n),intent(out) :: f,f1,f2,f3
        double precision,dimension(n),intent(in) :: eta
        double precision :: f21,f1n
        double precision,dimension(3) :: fn,k1,k2,k3,k4

        f = 0.0d0; f1 = 0.0d0; f2 = 0.0d0; f3 = 0.0d0;
        fn = (/0.0d0,0.0d0,f21/)
        do i=1,n
            f(i)  = fn(1)
            f1(i) = fn(2)
            f2(i) = fn(3)
            f3(i) = -0.50d0*fn(1)*fn(3)

            k1 = deta*fun_blasius(fn)
            k2 = deta*fun_blasius(fn + 0.50d0*k1)
            k3 = deta*fun_blasius(fn + 0.50d0*k2)
            k4 = deta*fun_blasius(fn +     k3)

            fn = fn + (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0
        enddo
        f1n = f1(n)
        end subroutine blasius_integrator
    function    fun_blasius(temp) result(funf)
        implicit none
        double precision, dimension(3) :: temp,funf
        funf(1) = temp(2)
        funf(2) = temp(3)
        funf(3) = -0.50d0*temp(1)*temp(3)
        end function fun_blasius
    end module Blasius

module CMM
    use global_variables
    use mylib
    implicit none
    contains
    subroutine NewtonRahpson(alpout, Wout, iter,Re,  inalp, inw )
        double precision :: inw, inalp, Re
        double precision :: alpout, Wout, alpha, omega 
        double precision :: df1_da, df2_da, df1_dw, df2_dw, iJacob(2,2), f1(2), f2(2), fun(2)
        double precision :: ap_fun(2), wp_fun(2)
        integer          :: iter 
        
        f2(1:2) = (/ inalp, inw /) 
        do iter = 1, 10000
            f1    = f2 
            alpha = f1(1)       ; omega = f1(2)       ;     fun = CMM_Integrator(Re, omega, alpha)
            alpha = f1(1) + 1d-4; omega = f1(2)       ;  ap_fun = CMM_Integrator(Re, omega, alpha)
            alpha = f1(1)       ; omega = f1(2) + 1d-4;  wp_fun = CMM_Integrator(Re, omega, alpha)
            
            df1_da = 1.0d4*(ap_fun(1) - fun(1) )
            df2_da = 1.0d4*(ap_fun(2) - fun(2) )
            df1_dw = 1.0d4*(wp_fun(1) - fun(1) )
            df2_dw = 1.0d4*(wp_fun(2) - fun(2) )
            
            iJacob(1,1:2) = (/ df2_dw,  -df1_dw /)/(df1_da*df2_dw - df2_da*df1_dw)
            iJacob(2,1:2) = (/-df2_da,   df1_da /)/(df1_da*df2_dw - df2_da*df1_dw)
            
            f2(1) = f1(1) - (iJacob(1,1)*fun(1) + iJacob(1,2)*fun(2))
            f2(2) = f1(2) - (iJacob(2,1)*fun(1) + iJacob(2,2)*fun(2))
            
            if(maxval(abs(f2-f1)) < 1d-10)   exit 
        enddo
        alpout = f1(1)
        wout   = f1(2)
        end subroutine NewtonRahpson

    function CMM_Integrator(Re, omega, alpha ) result (dpr)
        implicit none
        integer                     :: m,j, m1, i1, iRe
        double precision            :: alpha, Re, omega, p, dpr(2)
        double complex              ::  q , rescale
        double complex,dimension(4) :: ky1,ky2,ky3,ky4,ky5,ky6
        double complex,dimension(2) :: kphi1,kphi2, kphi3, kphi4
        double complex,dimension(2) :: utemp 
        double precision            :: dy, ymax, dy2
        character                   :: filename*31

        m    = n/2
        m1   = m/2
        dy   = 2.0d0*dely
        dy2  = 4.0d0*dely 
        ymax = y(n)

        p= dsqrt(alpha**2)
        q=cdsqrt(one*alpha**2+xi*Re*(alpha-omega))
        
        y1(m+1)= (one*p-q)*cdexp(-ymax*(one*p+q))
        y2(m+1)=-(one*p+q)*y1(m+1)
        y3(m+1)=(one*p**2+p*q+q**2)*y1(m+1)
        y4(m+1)=(p*q)*y1(m+1)
            y5(m+1)=-p*q*(one*p+q)*y1(m+1)
            y6(m+1)=((p*q)**2)*y1(m+1)
            rescale = 1.0d0 !exp(-dy*(p*one + q))
            
            y1(m+1)=y1(m+1)*rescale
            y2(m+1)=y2(m+1)*rescale 
            y3(m+1)=y3(m+1)*rescale
            y4(m+1)=y4(m+1)*rescale
            y5(m+1)=y5(m+1)*rescale
            y6(m+1)=y6(m+1)*rescale

             do j=m+1,2,-1
                !  print*, Re, alpha, omega
                ky1(1)=y2(j)
                ky2(1)=y4(j)+y3(j)
              ky3(1)=b1(Re, alpha, omega, U(2*j-1))*y2(j)+y5(j)
              ky4(1)=y5(j)
              ky5(1)=y6(j)+b1(Re, alpha, omega, U(2*j-1))*y4(j)+b2(Re, alpha, omega, U(2*j-1),U2(2*j-1))*y1(j)
              ky6(1)=b2(Re, alpha, omega, U(2*j-1),U2(2*j-1))*y2(j)
              
              ky1(2)=y2(j)-0.5d0*dy*ky2(1)
              ky2(2)=(y4(j)-0.5d0*dy*ky4(1))+(y3(j)-0.5d0*dy*ky3(1))
              ky3(2)=b1(Re, alpha, omega, U(2*j-2))*(y2(j)-0.5d0*dy*ky2(1))+(y5(j)-0.5d0*dy*ky5(1))
              ky4(2)=y5(j)-0.5d0*dy*ky5(1)
              ky5(2)=y6(j)-0.5d0*dy*ky6(1)+b1(Re, alpha, omega, U(2*j-2))*(y4(j)-0.5d0*dy*ky4(1))+&
              b2(Re, alpha, omega, U(2*j-2),U2(2*j-2))*(y1(j)-0.5d0*dy*ky1(1))
              ky6(2)=b2(Re, alpha, omega, U(2*j-2),U2(2*j-2))*(y2(j)-0.5d0*dy*ky2(1))
              
              ky1(3)=y2(j)-0.5d0*dy*ky2(2)
              ky2(3)=(y4(j)-0.5d0*dy*ky4(2))+(y3(j)-0.5d0*dy*ky3(2))
              ky3(3)=b1(Re, alpha, omega, U(2*j-2))*(y2(j)-0.5d0*dy*ky2(2))+(y5(j)-0.5d0*dy*ky5(2))
              ky4(3)=y5(j)-0.5d0*dy*ky5(2)
              ky5(3)=y6(j)-0.5d0*dy*ky6(2)+b1(Re, alpha, omega, U(2*j-2))*(y4(j)-0.5d0*dy*ky4(2))+&
              b2(Re, alpha, omega, U(2*j-2),U2(2*j-2))*(y1(j)-0.5d0*dy*ky1(2))
              ky6(3)=b2(Re, alpha, omega, U(2*j-2),U2(2*j-2))*(y2(j)-0.5d0*dy*ky2(2))
              
              ky1(4)=y2(j)-dy*ky2(3)
              ky2(4)=(y4(j)-dy*ky4(3))+(y3(j)-dy*ky3(3))
              ky3(4)=b1(Re, alpha, omega, U(2*j-3))*(y2(j)-dy*ky2(3))+(y5(j)-dy*ky5(3))
              ky4(4)=y5(j)-dy*ky5(3)
              ky5(4)=y6(j)-dy*ky6(3)+b1(Re, alpha, omega, U(2*j-3))*(y4(j)-dy*ky4(3))+&
                     b2(Re, alpha, omega, U(2*j-3),U2(2*j-3))*(y1(j)-dy*ky1(3))
                     ky6(4)=b2(Re, alpha, omega, U(2*j-3),U2(2*j-3))*(y2(j)-dy*ky2(3)) 
                     
                     y1(j-1)=y1(j)-(dy*(ky1(1)+2.0d0*ky1(2)+2.0d0*ky1(3)+ky1(4))/6.0d0)
              y2(j-1)=y2(j)-(dy*(ky2(1)+2.0d0*ky2(2)+2.0d0*ky2(3)+ky2(4))/6.0d0)
              y3(j-1)=y3(j)-(dy*(ky3(1)+2.0d0*ky3(2)+2.0d0*ky3(3)+ky3(4))/6.0d0)
              y4(j-1)=y4(j)-(dy*(ky4(1)+2.0d0*ky4(2)+2.0d0*ky4(3)+ky4(4))/6.0d0)
              y5(j-1)=y5(j)-(dy*(ky5(1)+2.0d0*ky5(2)+2.0d0*ky5(3)+ky5(4))/6.0d0)
              y6(j-1)=y6(j)-(dy*(ky6(1)+2.0d0*ky6(2)+2.0d0*ky6(3)+ky6(4))/6.0d0)
              
              y1(j-1)=y1(j-1)*rescale
             y2(j-1)=y2(j-1)*rescale 
             y3(j-1)=y3(j-1)*rescale
             y4(j-1)=y4(j-1)*rescale
             y5(j-1)=y5(j-1)*rescale
             y6(j-1)=y6(j-1)*rescale
            enddo
            dpr(1:2) = (/real(y1(1)), aimag(y1(1))/)
            
            ! call system('rm '//trim(outputfile))
            ! open(31,file= trim(outputfile))
            ! write(31,*) "Variables = a,w,v1r,v1i,v2r,v2i"
            ! write(31,*) 'zone T = "',trim(zname), '"'," I = ", na," J = ", nw 
            ! write(31,*) "DT = (DOUBLE, DOUBLE, DOUBLE,DOUBLE, DOUBLE,DOUBLE)"
            ! write(31,*) ' AUXDATA      Re                       = "',  Re, '"' 
            ! write(31,*) ' AUXDATA      omega0                   = "',  omega0, '"' 
            ! write(31,*) ' AUXDATA      na                       = "',  na, '"' 
            ! write(31,*) ' AUXDATA      nw                       = "',  nw, '"' 
            ! write(31,*) ' AUXDATA      aMax                     = "',  aMax, '"' 
            ! write(31,*) ' AUXDATA      wMax                     = "',  wMax, '"' 
            ! write(31,*) ' AUXDATA      ai                       = "',  ai, '"' 
            ! write(31,*) ' AUXDATA      wi                       = "',  wi, '"' 

            ! ! write(31,*) "Zone T = S I = ",na,"J = ",glo_nw
            ! close(31)
        end function CMM_Integrator
    function fun(arg1, var1, var2, var4) result(arg2)
        double complex :: arg1(1:2), arg2(1:2), var1, var2, var4 
        arg2(1) = arg1(2)
        arg2(2) = ( var2*arg1(2) - var4*arg1(1) )/var1
        end function fun
    double complex function b1(Re, alpha, omega,V)
        implicit none
        double precision,intent(in) :: V
        double precision :: Re, alpha, omega
        b1=2.0d0*(alpha**2)+xi*Re*(alpha*V-omega)
        end function b1
    double complex function b2(Re, alpha, omega,V,W)
        double precision,intent(in) :: V,W
        double precision :: Re, alpha, omega
        b2=alpha**4+xi*Re*alpha*W+xi*Re*(alpha*V-omega)*(alpha**2)
        end function b2
    end module CMM

program BCIM_main
    use global_variables
    use mylib
    use Blasius
    use CMM
    implicit none
    integer :: ire, iter 
    double precision :: var1, var2 , g1, g2, Re

    allocate(y(n), u(n), U2(n)) 
    allocate(vre(nre), vw(nre))
    allocate(y1(n/2 +1),y2(n/2 +1),y3(n/2 +1),y4(n/2 +1),y5(n/2 +1),y6(n/2 +1)) 

    deta = etamax/dble(n-1)
    call blasius_profile
    call linespace(vre, re0, ren, nre)

    ! g1 = 0.350d0 
    ! g2 = 0.1250d0 
    ! do ire = 1, nre ,10
    !     Re = vre(ire)
    !     call NewtonRahpson(var1, var2, iter, Re, g1, g2 )
    !     write(31,*) Re, var2, var1
    !     g1 = var1 
    !     g2 = var2
    !     write(*,'(I4, 3F25.16)') iter,Re, var1, var2
    ! enddo
    ! stop 

    open(31,file = 'neutralcurve_up.dat')
    write(31,*) " Variables = Re<sub><greek>d*, <greek>w<sub>0, <greek>a"
    write(31,*) 'zone T = "up" I = ',2*nre 
    write(31,*) 'DT = (DOUBLE DOUBLE DOUBLE)'
    
    g1 = 0.190400d0 
    g2 = 0.044630d0

    do ire = nre, 1,-1 
        Re = vre(ire)
        call NewtonRahpson(var1, var2, iter, Re, g1, g2 )
        write(31,*) Re, var2, var1
        g1 = var1 
        g2 = var2
        write(*,'(I4, 3F25.16)') iter,Re, var1, var2
    enddo

    g1 = 0.30d0 
    g2 = 0.10d0

    do ire = 1, nre 
        Re = vre(ire)
        call NewtonRahpson(var1, var2, iter, Re, g1, g2 )
        write(31,*) Re, var2, var1
        g1 = var1 
        g2 = var2
        write(*,'(I4, 3F25.16)') iter,Re, var1, var2
    enddo




        close(31)

    end program BCIM_main
