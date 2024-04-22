clc
clear all 
close all 

%% Parameter Spcification
U_inf = 1;
L = 10;
mu = 1.789E-5;
beta = 1/35;
rho = 1.225;
nu = mu/rho;
A = sqrt(nu/U_inf);
h = 0.1;
% Boundary Limits
eta = 0:h:10;
x = 0:h:10;
% Boundary Conditions
y1(1) = 0;
y2(1) = 0;
y3(1) = 0.5055;


%% Falkner Skan Solution with Runge-Kutta Method
f1 = @(x, y1, y2, y3) y2;
f2 = @(x, y1, y2, y3) y3;
f3 = @(x, y1, y2, y3) -y1*y3-beta*(1-y2^2);
for i = 1:(length(eta)-1)
a = h.*[f1(eta(i), y1(i), y2(i), y3(i)), f2(eta(i), y1(i), y2(i), y3(i)), f3(eta(i), y1(i), y2(i), y3(i))];
b = h.*[f1(eta(i), y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2), f2(eta(i)+h/2, y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2), f3(eta(i)+h/2, y1(i)+a(1)/2, y2(i)+a(2)/2, y3(i)+a(3)/2)];
c = h.*[f1(eta(i), y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2), f2(eta(i)+h/2, y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2), f3(eta(i)+h/2, y1(i)+b(1)/2, y2(i)+b(2)/2, y3(i)+b(3)/2)];
d = h.*[f1(eta(i), y1(i)+c(1), y2(i)+c(2), y3(i)+c(3)), f2(eta(i)+h, y1(i)+c(1), y2(i)+c(2), y3(i)+c(3)), f3(eta(i)+h, y1(i)+c(1), y2(i)+c(2), y3(i)+c(3))];
y3(i+1) = y3(i)+ 1/6*(a(3)+2*b(3)+2*c(3)+d(3));
y2(i+1) = y2(i)+ 1/6*(a(2)+2*b(2)+2*c(2)+d(2));
y1(i+1) = y1(i)+ 1/6*(a(1)+2*b(1)+2*c(1)+d(1));
end

%% Plotting
figure(1)
plot(y1, eta, y2, eta, y3, eta, 'LineWidth', 2)
xlim([0 2])
title('Solution of Falkner Skan eqution \beta = 1/35', 'FontSize', 14);
xlabel('f, f'' and f''''', 'FontSize', 20);
ylabel('\eta', 'FontSize', 20);
xlim([-0.2 1.2])
grid on
Legend1 = {'f(\eta)', 'f''(\eta)', 'f''''(\eta)'};
legend(Legend1, 'FontSize', 14);

for i=1:length(eta)-1
y4(i) = (y3(i+1)-y3(i))/(eta(i+1)-eta(i));
end

