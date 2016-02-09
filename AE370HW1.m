close all; clear all; clear figure;
%%Problem 1
f = @(x) sinh(x);
Problem1_Check = integral(f,2,6);
%Rectangle
a = 2; b = 6;
h = (b-a)/3;
Rectangle_Rule_Approximation = h * (f(8/3)+f(12/3)+f(16/3));
%Trapezoid
Trapezoidal_Rule_Approximation = (h/2)*(f(2)+2*f(4)+f(6));
%Simpson
Simpson_Rule_Approximation = (h/2)*(f(2)+4*f(4)+f(6));
%%Problem 2
f2 = @(x) cos((x.^3)+3.*x)./exp(x.^2);
Problem2_Check = integral(f2,0, 2*pi);
TrapezoidalValues = zeros(30,1);
k = 0;
MAX_INTERVAL = 30;
RectangleValues = zeros(2,MAX_INTERVAL);
TrapValues = zeros(2,MAX_INTERVAL);
while k < MAX_INTERVAL
    k = k+1;
    w = (2*pi)/(k);
    Ir = 0;
    It = 0;
    j = 0;
    while j < k
        j = j + 1;
        x = (j-1)*w+(w/2);
        z0 = x - (w/2);
        z1 = x + (w/2);
        Ir = Ir + f2(x)*w;
        It = It + (f2(z0) + f2(z1))*.5*w; 
    end
    RectangleValues(1,k) = k;
    RectangleValues(2,k) = Ir;
    TrapValues(1,k) = k;
    TrapValues(2,k) = It;
end
figure(1)
plot(RectangleValues(1,:),RectangleValues(2,:),'b.-',TrapValues(1,:),TrapValues(2,:),'r.-');
legend('Rectangle', 'Trapezoidal');
ylabel('Integral value');
xlabel('Number of Intervals');
%%Problem 3
f3 = @(x,y) x.*sin(x.^2)+log(2 + y);
Problem3_Check = integral2(f3, -1, 1, -1, 1);
% 1 x 1
GQ1 = 4 * f3(0,0);
% 2 x2 
eta = sqrt(3)^-1;
zeta = sqrt(3)^-1;
GQ2 = 1 * (f3(eta, zeta)) + 1 * (f3(-1*eta, zeta)) + 1 * (f3(eta, -1*zeta)) + 1 * (f3(-1*eta, -1 * zeta));
%3 x 3
w1 = 8/9;
w2 = 5/9;
eta = sqrt(3/5);
zeta = sqrt(3/5);
GQ3 = w1 * w1 * f3(0,0) + w1 * w2*f3(eta, 0) + w1 * w2 * f3(-1 * eta, 0) + w1 * w2 * f3(0, eta) + w1 * w2 * f3(0, -1 * eta) + w2 * w2 * f3(eta, zeta) + w2 * w2 * f3(-1 * eta, zeta) + w2 * w2 * f3(eta, -1 * zeta) + w2 * w2 * f3(-1 * eta, -1 *zeta);

