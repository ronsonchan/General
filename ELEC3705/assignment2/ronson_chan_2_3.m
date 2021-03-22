%% ronson_chan_2_3

close all
clear all
clc

%% definition
sx =  [0,1;1,0];
sy =  [0,1i;-1i,0];
sz =  [1,0;0,-1];
I = eye(2);
sz1 = kron(sz,I);
sz2 = kron(I,sz);
sx1 = kron(sx,I);
sx2 = kron(I,sx);
sy1 = kron(sy,I);
sy2 = kron(I,sy);

Sz1 = kron(sz1,I);
Sz2 = kron(sz2,I);
Sz3 = kron(I,sz2);

Sx1 = kron(sx1,I);
Sx2 = kron(sx2,I);
Sx3 = kron(I,sx2);

Sy1 = kron(sy1,I);
Sy2 = kron(sy2,I);
Sy3 = kron(I,sy2);

up = [1;0];
down = [0;1];

upup = kron(up,up);
updown = kron(up,down);
downup = kron(down,up);
downdown = kron(down,down);

uuu = kron(upup,up);
uud = kron(upup,down);
udu = kron(updown,up);
udd = kron(updown,down);
duu = kron(downup,up);
dud = kron(downup,down);
ddu = kron(downdown,up);
ddd = kron(downdown,down);

%% calculating expectation value for 1
clc
% define time here
resolution = 1000;
time  = linspace(0,1,resolution);
T = length(time);
hbar = 1;
ketX = (1/sqrt(2)) * (up+down);
ketXXX = kron(ketX, kron(ketX, ketX));
psi0 = ketXXX;
psi = zeros(8,T);
J = 1;
% hamiltonian is written in units of angular freq.
Bz = [0,1,100];
exp_x = zeros(3,T,3);   
exp_y = zeros(3,T,3);   
exp_z = zeros(3,T,3);   
% reallocating all the S into a 3D matrix
Sx = zeros(8,8,3);
Sy = zeros(8,8,3);
Sz = zeros(8,8,3);
Sx(:,:,1) = Sx1;
Sx(:,:,2) = Sx2;
Sx(:,:,3) = Sx3;
Sy(:,:,1) = Sy1;
Sy(:,:,2) = Sy2;
Sy(:,:,3) = Sy3;
Sz(:,:,1) = Sz1;
Sz(:,:,2) = Sz2;
Sz(:,:,3) = Sz3;

for i = 1:3
    H = J*(Sz1.*Sz2 + Sz2.*Sz3) - Bz(i)*(Sz1+Sz2+Sz3);
    for t = 1:T
        U = expm(-1i*H*time(t));
        psi(:,t) = U*psi0;
        for j = 1:3
            exp_x(i,t,j) = real(psi(:,t)' * Sx(:,:,j) * psi(:,t));
            exp_y(i,t,j) = real(psi(:,t)' * Sy(:,:,j) * psi(:,t));
            exp_z(i,t,j) = real(psi(:,t)' * Sz(:,:,j) * psi(:,t));
        end
    end
end
%% plotting
close all;

figure(1);
subplot(3,1,1)
hold on
plot(time, exp_x(1,:,1));
plot(time, exp_y(1,:,1));
plot(time, exp_z(1,:,1));
title('Bz = 0');
xlabel('time');
ylabel('Expectation Value');
legend('expectation value of x1','expectation value of y1','expectation value of z1');
hold off
subplot(3,1,2)
hold on
plot(time, exp_x(1,:,2));
plot(time, exp_y(1,:,2));
plot(time, exp_z(1,:,2));
title('Bz = 0');
xlabel('time');
ylabel('Expectation Value');
legend('expectation value of x2','expectation value of y2','expectation value of z2');
hold off
subplot(3,1,3)
hold on
plot(time, exp_x(1,:,3));
plot(time, exp_y(1,:,3));
plot(time, exp_z(1,:,3));
title('Bz = 0');
xlabel('time');
ylabel('Expectation Value');
legend('expectation value of x3','expectation value of y3','expectation value of z3');
hold off


figure(2);
subplot(3,1,1)
hold on
plot(time, exp_x(2,:,1));
plot(time, exp_y(2,:,1));
plot(time, exp_z(2,:,1));
title('Bz = 1');
xlabel('time');
ylabel('Expectation Value');
legend('expectation value of x2','expectation value of y2','expectation value of z2');
hold off
subplot(3,1,2)
hold on
plot(time, exp_x(2,:,2));
plot(time, exp_y(2,:,2));
plot(time, exp_z(2,:,2));
title('Bz = 1');
xlabel('time');
ylabel('Expectation Value');
legend('expectation value of x2','expectation value of y2','expectation value of z2');
hold off
subplot(3,1,3)
hold on
plot(time, exp_x(2,:,3));
plot(time, exp_y(2,:,3));
plot(time, exp_z(2,:,3));
title('Bz = 1');
xlabel('time');
ylabel('Expectation Value');
legend('expectation value of x2','expectation value of y2','expectation value of z2');
hold off

figure(3);
subplot(3,1,1)
hold on
plot(time, exp_x(3,:,1));
plot(time, exp_y(3,:,1));
plot(time, exp_z(3,:,1));
title('Bz = 100');
xlabel('time');
ylabel('Expectation Value');
legend('expectation value of x3','expectation value of y3','expectation value of z3');
hold off
subplot(3,1,2)
hold on
plot(time, exp_x(3,:,2));
plot(time, exp_y(3,:,2));
plot(time, exp_z(3,:,2));
title('Bz = 100');
xlabel('time');
ylabel('Expectation Value');
legend('expectation value of x3','expectation value of y3','expectation value of z3');
hold off
subplot(3,1,3)
hold on
plot(time, exp_x(3,:,3));
plot(time, exp_y(3,:,3));
plot(time, exp_z(3,:,3));
title('Bz = 100');
xlabel('time');
ylabel('Expectation Value');
legend('expectation value of x3','expectation value of y3','expectation value of z3');
hold off