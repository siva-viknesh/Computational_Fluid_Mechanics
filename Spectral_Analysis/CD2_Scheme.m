%************************* Property Chart ********************************%
format long;
clc;
N     =  101;                %************* No. of Points ****************%
Vkh   =  0:0.01:pi ;
Lkh   =  length (Vkh);
j     =  (N-1)/2 ;           %*********** Interested Point ***************%
A1    =  zeros (N,N);
A2    =  zeros (N,N);
B1    =  zeros (N,N);
B2    =  zeros (N,N);
C1    =  zeros (N,N);
C2    =  zeros (N,N); 
keq_k =  zeros (1,Lkh);
f     =  1;                %********** Divident of 'C' Matrix **********%

%********************** Declaration of Matrix A1 *************************%
p31 =0;        p32=0;      p33=1;      p34=0;      p35=0;

for i = 3: N-2
   A1(i,i-2) = p31;  A1(i,i-1) = p32; A1(i,i) = p33;   
   A1(i,i+1) = p34;  A1(i,i+2) = p35;    
end
   A1(1,1)= p33; A1(1,2)= p34; A1(1,3)= p35; A1(1,N)= p32; A1(1,N-1)= p31;
   
   A1(2,1)= p32; A1(2,2)= p33; A1(2,3)= p34; A1(2,4)= p35; A1(2,N)= p31;
   
   A1(N-1,1)=p35; A1(N-1,N-3)=p31; A1(N-1,N-2)=p32; A1(N-1,N-1)= p33; 
   A1(N-1,N)= p34;
   
   A1(N,N-2)=p31; A1(N,N-1)=p32; A1(N,N)=p33; A1(N,1) = p34; A1(N,2)=p35;   
   
%********************** Declaration of Matrix B1**************************%
q31 = 0;        q32 = -1;      q33 = 0;        q34 = 1;         q35 = 0;

for i = 3: N-2
   C1(i,i-2) = q31;  C1(i,i-1) = q32; C1(i,i) = q33;   
   C1(i,i+1) = q34;  C1(i,i+2) = q35;    
end
   C1(1,1)= q33; C1(1,2)= q34; C1(1,3)= q35; C1(1,N)= q32; C1(1,N-1)= q31;
   
   C1(2,1)= q32; C1(2,2)= q33; C1(2,3)= q34; C1(2,4)= q35; C1(2,N)= q31;
   
   C1(N-1,1)= q35; C1(N-1,N-3)= q31; C1(N-1,N-2)= q32; C1(N-1,N-1)= q33; 
   C1(N-1,N)= q34;
   
   C1(N,N-2)= q31; C1(N,N-1)= q32; C1(N,N)= q33; C1(N,1) = q34; C1(N,2)= q35;
    
%********************** Declaration of Matrix C **************************%
   
D1 = inv ( A1 - (B1/B2)*A1)* (      )   

for m = 1 : Lkh
   
   kh  = Vkh (m);
   sum = 0 ;
   
   for n = 1 : N
       
      sum     = sum + D1(j,n)*exp(-1i*kh*(j-n)); 
       
   end
    keq_k (m) =   -1i* (sum /kh);    
end

%**************************** File Writing *******************************%
y  = [Vkh; real(keq_k)];
y1 = [Vkh; imag(keq_k)];

fileID = fopen ('Real_keq_k_CD2.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID,'%6.10f %6.14f\n',y);
fclose(fileID);

fileID = fopen ('Imag_keq_k_CD2.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID,'%6.10f %6.14f\n',y1);
fclose(fileID);

%************************* Plotting the values ***************************%
plot (Vkh,real(keq_k))
xlabel('kh')
ylabel('keq / k')
