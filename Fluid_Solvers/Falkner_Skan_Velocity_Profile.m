%--------------- FALKNER-SKAN VELOCITY PROFILE COMPUTATION ---------------%
clear all ;
close all ;
clc ;
format long
N = 1024 ; 
eta_max = 12 ;
eta_min = 0 ;
deta = ( eta_max - eta_min ) / (N-1) ;
eta = eta_min : deta : eta_max ;
m = 0 ;                           %-----> Parameter
b = 2*m /(m+1);
f = zeros (N,1) ; f (:) = 1 ;
g = zeros (N,1) ; g (:) = 1 ;
h = zeros (N,1) ; h (:) = 1 ;
%---------> Approaching from wall towards freestream
f(1) = 0 ; g(1) = 0 ; 
residual = 1  ;
h_temp = h (1);
g_temp = g (N);
m = 1 ;
while residual > 10d-6
for i = 2:N
    f (i) = f(i-1) + g(i-1)* deta;
    g (i) = g(i-1) + h(i-1)* deta;
    h (i) = h(i-1) - (f(i-1)*h(i-1)+b*(1- g(i-1)^2))* deta ;   
end
    residual = abs(g(N)- 1);
    if m < 2
    %------> Newton Raphson method for better guess of h(1)
    h(1) = h(1) - (g(N) - 1.0) * (h(1)-h_temp) / (g(N) - g_temp) ;  
    else
        h(1) = h (1) + 1 ;
    end
    h_temp = h(1) ;
    g_temp = g(N) ;
    m = m + 1 ;
end
plot (g , eta)
hold on  
plot (h, eta)
