%************************* Property Chart ********************************%
format long;
clc;
clear all;
N     =  301;                %************* No. of Points ****************%
Vkh = linspace(1d-6, pi, N); 
Lkh = N; 
j     =  N-3 ;         %************ Interested Point ***************%
%************************* First Derivative ******************************%

alpha =   1/3   ;
a     =  7/9     ;
b     =  1/36    ;
keq_k_1 =  zeros (1,Lkh);
keq_k_2 =  zeros (1,Lkh);

%********************** Declaration of Matrix A1 *************************%
p32=alpha;      p33=1;      p34=alpha;      
A1 = toeplitz([p33 p32 zeros(1, N-2)], [p33 p34 zeros(1, N-2)]) ;
A1(1,:) = 0 ;
A1(N,:) = 0 ;
A1(2,:)   = 0 ;
A1(N-1,:) = 0 ;

A1(1,1) = 1 ;
A1(N,N) = 1 ;
A1(2,2) = 1 ;
A1(N-1,N-1) = 1 ;
%  A1 (1,N) = p32 ;
%  A1 (N,1) = p34 ; 
%********************** Declaration of Matrix B1 *************************%
p32=-a;      p33=0;      p34=a;  
B1 = toeplitz([p33 p32 zeros(1, N-2)], [p33 p34 zeros(1, N-2)]) ;
B1(1,:)   = 0 ;
B1(N,:)   = 0 ; 
B1(2,:)   = 0 ;
B1(N-1,:) = 0 ;
%  B1 (1,N) = p32 ;
%  B1 (N,1) = p34 ;    
%********************** Declaration of Matrix C1**************************%
q32 = -b;      q33 = 0;        q34 = b;   
C1 = toeplitz([q33 0 q32 zeros(1, N-3)], [q33 0 q34 zeros(1, N-3)]) ;
C1(1,:) = 0 ;
C1(N,:) = 0 ;
C1(2,:)   = 0 ;
C1(N-1,:) = 0 ;
%    C1 (1,N) = q32 ;
%    C1 (N,1) = q34 ;

E1 = B1 + C1 ;   

E1(1,:) = 0 ;
E1(N,:) = 0 ;
E1(1,1) = -1.5 ;        E1(N,N)  = 1.5 ;
E1(1,2) = 2 ;           E1(N,N-1)= -2;
E1(1,3) = -0.5;         E1(N,N-2)= 0.5 ;

% j = 1
beta = -0.025 ; 
E1(2,1) = (2*beta - 1)/3   ;
E1(2,2) = -(8*beta/3 + 0.5);
E1(2,3) = 4*beta +1 ;
E1(2,4) = -(8*beta/3 +1/6);
E1(2,5) = 2*beta/3  ;

% j = N-1
beta = 0.09 ; 
E1(N-1,1) = -(2*beta - 1)/3   ;
E1(N-1,2) = (8*beta/3 + 0.5);
E1(N-1,3) = -(4*beta +1) ;
E1(N-1,4) = (8*beta/3 +1/6);
E1(N-1,5) = -2*beta/3  ;

D1 = A1 \ E1;

%************************* Second Derivative *****************************%
alpha =  2/11     ;
a     =  12/11    ;
b     =  3/44    ;
%********************** Declaration of Matrix A2 ************************
p32=alpha;      p33=1;      p34=alpha;              
A2 = toeplitz([p33 p32 zeros(1, N-2)], [p33 p34 zeros(1, N-2)]) ;
A2(1,:) = 0 ;
A2(N,:) = 0 ;
A2(2,:)   = 0 ;
A2(N-1,:) = 0 ;

A2(1,1) = 1 ;
A2(N,N) = 1 ;
A2(2,2) = 1 ;
A2(N-1,N-1) = 1 ;
%    A2 (1,N) = p32 ;
%    A2 (N,1) = p34 ;
   
%********************** Declaration of Matrix B2 *************************%
   p32=b;      p33=-2*b;      p34=b;      

   B2 = toeplitz([p33 0 p32 zeros(1, N-3)], [p33  0 p34 zeros(1, N-3)]) ;
B2(1,:)   = 0 ;
B2(N,:)   = 0 ; 
B2(2,:)   = 0 ;
B2(N-1,:) = 0 ;
%    B2 (1,N) = p32 ;
%    B2 (N,1) = p34 ;    
   
%********************* Declaration of Matrix C2 *************************%
q32 = a;      q33 =-2*a;        q34 = a;      
C2 = toeplitz([q33 q32 zeros(1, N-2)], [q33 q34 zeros(1, N-2)]) ;
C2(1,:) = 0 ;
C2(N,:) = 0 ;
C2(2,:)   = 0 ;
C2(N-1,:) = 0 ;
%    C2 (1,N) = q32 ;
%    C2 (N,1) = q34 ; 
   
E2 = B2 + C2 ;

E2(1,:) = 0 ;
E2(N,:) = 0 ;

E2(1,1) = 1 ;           E2(N,N)  = 1 ;
E2(1,2) = -2 ;          E2(N,N-1)= -2;
E2(1,3) = 1;            E2(N,N-2)= 1 ;  

% j = 2
beta = -0.025 ; 
E2(2,1) = (2*beta - 1)/3   ;
E2(2,2) = -(8*beta/3 + 0.5);
E2(2,3) = 4*beta +1 ;
E2(2,4) = -(8*beta/3 +1/6);
E2(2,5) = 2*beta/3  ;

% j = N-1
beta = 0.09 ; 
E2(N-1,1) = -(2*beta - 1)/3   ;
E2(N-1,2) = (8*beta/3 + 0.5);
E2(N-1,3) = -(4*beta +1) ;
E2(N-1,4) = (8*beta/3 +1/6);
E2(N-1,5) = -2*beta/3  ;   

D2 = A2\ E2;

for m = 1 : Lkh   
   kh   = Vkh (m);
   sum1 = 0 ;
   sum2 = 0 ;
   for n = 1 : N
      sum1  = sum1 + D1(j,n)*(exp(1i*kh*(n-j))); 
      sum2  = sum2 + D2(j,n)*(exp(1i*kh*(n-j)));  
   end   
    keq_k_1 (m) =  (sum1 / (kh*1i));
    keq_k_2 (m) =  -(sum2 / ((kh)^2));      
end

%**************************** File Writing *******************************%
y  = [Vkh; real(keq_k_1)];
y1 = [Vkh; imag(keq_k_1)];
y2 = [Vkh; real(keq_k_2)];
y3 = [Vkh; imag(keq_k_2)];

fileID = fopen ('Real_keq_k_1_LELE_N_3.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID, '%6.10f %6.14f\n',y);
fclose(fileID);

fileID = fopen ('Imag_keq_k_1_LELE_N_3.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID,'%6.10f %6.14f\n',y1);
fclose(fileID);

fileID = fopen ('Real_keq_k_2_LELE_N_3.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID, '%6.10f %6.14f\n',y2);
fclose(fileID);

fileID = fopen ('Imag_keq_k_2_LELE_N_3.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID,'%6.10f %6.14f\n',y3);
fclose(fileID);
 
 
 
