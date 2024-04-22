%************************* Property Chart ********************************%
format long;
clc;
N     =  101;                %************* No. of Points ****************%
dx = 0.1;
x = 0: dx: 1;
a = 1.0;
b = 4.0;
c = 1.0;
d = 3.0;


%********************** Declaration of Matrix A **************************%
D = 0.3793894912 ;
E = 1.57557379   ;
F = 0.183205192  ;
G = -2.0 ;
p31 =0;        p32=D-G/60;      p33=1;      p34=D+G/60;      p35=0;

for i = 3: N-2
   A(i,i-2) = p31;  A(i,i-1) = p32; A(i,i) = p33;   
   A(i,i+1) = p34;  A(i,i+2) = p35;    
end
   A(1,1)= p33; A(1,2)= p34; A(1,3)= p35; A(1,N)= p32; A(1,N-1)= p31;
   
   A(2,1)= p32; A(2,2)= p33; A(2,3)= p34; A(2,4)= p35; A(2,N)= p31;
   
   A(N-1,1)=p35; A(N-1,N-3)=p31; A(N-1,N-2)=p32; A(N-1,N-1)= p33; 
   A(N-1,N)= p34;
   
   A(N,N-2)=p31; A(N,N-1)=p32; A(N,N)=p33; A(N,1) = p34; A(N,2)=p35;   
   
%********************** Declaration of Matrix B **************************%
q31 = -F/4+G/300;        q32 = -E/2+G/30;      q33 = -11*G/150;
q34 = E/2+G/30;         q35 = F/4+G/300;

for i = 3: N-2
   B(i,i-2) = q31;  B(i,i-1) = q32; B(i,i) = q33;   
   B(i,i+1) = q34;  B(i,i+2) = q35;    
end
   B(1,1)= q33; B(1,2)= q34; B(1,3)= q35; B(1,N)= q32; B(1,N-1)= q31;
   
   B(2,1)= q32; B(2,2)= q33; B(2,3)= q34; B(2,4)= q35; B(2,N)= q31;
   
   B(N-1,1)= q35; B(N-1,N-3)= q31; B(N-1,N-2)= q32; B(N-1,N-1)= q33; 
   B(N-1,N)= q34;
   
   B(N,N-2)= q31; B(N,N-1)= q32; B(N,N)= q33; B(N,1) = q34; B(N,2)= q35;
    
%********************** Declaration of Matrix C **************************%
   
C = (A\B)/f;   

for m = 1 : Lkh
   
   kh  = Vkh (m);
   sum = 0 ;
   
   for n = 1 : N
       
      sum     = sum + C(j,n)*exp(-1i*kh*(j-n)); 
       
   end
    keq_k (m) =   -1i* (sum /kh);    
end

%**************************** File Writing *******************************%
y  = [Vkh; real(keq_k)];
y1 = [Vkh; imag(keq_k)];

fileID = fopen ('Real_keq_k_OUCS3.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID,'%6.10f %6.14f\n',y);
fclose(fileID);

fileID = fopen ('Imag_keq_k_OUCS3.dat','w');
fprintf(fileID, 'variables = kh,keq/k \n');
fprintf(fileID,'%6.10f %6.14f\n',y1);
fclose(fileID);

%************************* Plotting the values ***************************%
plot (Vkh,real(keq_k))
xlabel('kh')
ylabel('keq / k')
