format long;
clc;
clear all;
N     =  301;                %************* No. of Points ****************%
Vkh = linspace(1d-6, pi, N); 
Lkh = N; 
j     =  N-3 ;         %************ Interested Point ***************%
keq_k_1 =  zeros (1,Lkh);
keq_k_2 =  zeros (1,Lkh);
a = -1;
b = 0;
c = 1;
%********************** Declaration of Matrix C1**************************%
   q32 = a;      q33 = b;        q34 = c;   
C1 = toeplitz([q33 q32 zeros(1, N-2)], [q33 q34 zeros(1, N-2)]) ;
C1 = C1 * 0.5;

C1(1,:)   = 0 ;
C1(N,:)   = 0 ;

C1(1,1) = -1.5 ;        C1(N,N)  = 1.5 ;
C1(1,2) = 2 ;           C1(N,N-1)= -2;
C1(1,3) = -0.5;         C1(N,N-2)= 0.5 ;

%    C (1,N) = q32 ;
%    C (N,1) = q34 ;  

C2 = C1*C1;

for m = 1 : Lkh   
   kh   = Vkh (m);
   sum1 = 0 ;
   sum2 = 0 ;
   for n = 1 : N
      sum1  = sum1 + C1(j,n)*(exp(1i*kh*(n-j))); 
      sum2  = sum2 + C2(j,n)*(exp(1i*kh*(n-j)));  
   end   
    keq_k_1 (m) =  sum1 / (kh*1i);
    keq_k_2 (m) = - sum2 / (kh*kh);      
end

%**************************** File Writing *******************************%
y  = [Vkh; real(keq_k_1)];
y1 = [Vkh; imag(keq_k_1)];
y2 = [Vkh; real(keq_k_2)];
y3 = [Vkh; imag(keq_k_2)];

fileID = fopen ('Real_keq_k_1_CD2_N_3.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID, '%6.10f %6.14f\n',y);
fclose(fileID);

fileID = fopen ('Imag_keq_k_2_CD2_N_3.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID,'%6.10f %6.14f\n',y1);
fclose(fileID);

fileID = fopen ('Real_keq_k_2_CD2_N_3.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID, '%6.10f %6.14f\n',y2);
fclose(fileID);

fileID = fopen ('Imag_keq_k_2_CD2_N_3.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID,'%6.10f %6.14f\n',y3);
fclose(fileID);

figure (3)
plot (Vkh,real(keq_k_1))
xlabel('kh')
ylabel('keq / k')
