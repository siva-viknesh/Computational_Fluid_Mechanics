%************************* Property Chart ********************************%
format short;
clc;
N     =  301;                %************* No. of Points ****************%
Vkh   =  linspace(1d-6, 2.0*pi, N);
dx    =  (2.0*pi)/(N+1);
ye    =  sin (Vkh);
yde   =  cos (Vkh);
ydde  = -sin(Vkh);
Lkh   =  length (Vkh);
j     =  (N-1)/2    ;         
f1    =  7/16       ;      
f2    =  -dx/16      ;      
f3    =  15/(16*dx) ;  
f4    =  9/(8*dx)        ;     
f5    =  -1/8       ;      
f6    =  3/(dx^2)          ;     
keq_k_1 =  zeros (1,Lkh);
keq_k_2 =  zeros (1,Lkh);

%********************** Declaration of Matrix A1 *************************%
p31 =0;        p32=f1;      p33=1;      p34=f1;      p35=0;
 A1 = toeplitz([p33 p32 zeros(1, N-2)], [p33 p34 zeros(1, N-2)]) ;
 A1 (1,N) = p32 ;
 A1 (N,1) = p34 ;

 A1
   
%********************** Declaration of Matrix B1 *************************%
    p31 =0;        p32=-f2;      p33=0;      p34=f2;      p35=0;

 B1 = toeplitz([p33 p32 zeros(1, N-2)], [p33 p34 zeros(1, N-2)]) ;
 B1 (1,N) = p32 ;
 B1 (N,1) = p34 ;
   
%********************** Declaration of Matrix C1**************************%
 q31 = 0;        q32 = -f3;      q33 = 0;        q34 = f3;         q35 = 0;
 
   C1 = toeplitz([q33 q32 zeros(1, N-2)], [q33 q34 zeros(1, N-2)]) ;
   C1 (1,N) = q32 ;
   C1 (N,1) = q34 ;   
%********************** Declaration of Matrix A2 *************************%
        p31 =0;        p32=-f4;      p33=0;      p34=f4;      p35=0;
        
   A2 = toeplitz([p33 p32 zeros(1, N-2)], [p33 p34 zeros(1, N-2)]) ;
   A2 (1,N) = p32 ;
   A2 (N,1) = p34 ;
%********************** Declaration of Matrix B2 *************************%
    p31 =0;        p32=f5;      p33=1;      p34=f5;      p35=0;

   B2 = toeplitz([p33 p32 zeros(1, N-2)], [p33 p34 zeros(1, N-2)]) ;
   B2 (1,N) = p32 ;
   B2 (N,1) = p34 ;    
   
%********************* Declaration of Matrix C2 *************************%
q31 = 0;        q32 = f6;      q33 =-2*f6;        q34 = f6;       q35 = 0;

   C2 = toeplitz([q33 q32 zeros(1, N-2)], [q33 q34 zeros(1, N-2)]) ;
   C2 (1,N) = q32 ;
   C2 (N,1) = q34 ; 
    
%********************** Declaration of Matrix D1 & D2 ********************%
   
D1 = inv(A1 - ((B1*inv(B2))*A2)) * (C1-((B1*inv(B2))*C2)) ; 
D2 = inv(B2 - ((A2*inv(A1))*B1)) * (C2-((A2*inv(A1))*C1)) ;


ypn  = (D1 * ye') ;
yppn = (D2 * ye') ;

figure(1)
plot(Vkh,yde,Vkh,ypn)

figure(2)
plot(Vkh,ydde,Vkh,yppn)
xlim([0 2*pi]) 
ylim([-1 1])

%%

for m = 1 : Lkhtoeplitz
   
   kh  = Vkh (m);
   sum1 = 0 ;
   sum2 = 0 ;
    
   for n = 3 : N-3
      sum1  = sum1 + D1(j,n)*exp(-1i*kh*(j-n)); 
      sum2  = sum2 + D2(j,n)*exp(-1i*kh*(j-n));  
   end
   
    keq_k_1 (m) =  sum1 / (1i*kh);
    keq_k_2 (m) =  -(sum2 / ((kh)^2));
      
end
figure (3)
plot (Vkh,imag(keq_k_2))
xlabel('kh')
ylabel('keq / k')
%%
%**************************** File Writing *******************************%
y  = [Vkh; real(keq_k_1)];
y1 = [Vkh; imag(keq_k_1)];
y2 = [Vkh; real(keq_k_2)];
y3 = [Vkh; imag(keq_k_2)];

fileID = fopen ('Real_keq_k_1_NCCD.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID, '%6.10f %6.14f\n',y);
fclose(fileID);

fileID = fopen ('Imag_keq_k_1_NCCD.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID,'%6.10f %6.14f\n',y1);
fclose(fileID);

fileID = fopen ('Real_keq_k_2_NCCD.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID, '%6.10f %6.14f\n',y2);
fclose(fileID);

fileID = fopen ('Imag_keq_k_2_NCCD.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID,'%6.10f %6.14f\n',y3);
fclose(fileID);

%************************* Plotting the values ***************************%
