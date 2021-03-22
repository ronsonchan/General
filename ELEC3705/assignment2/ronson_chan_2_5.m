%% ronson_chan_2_5
clc
close all
clear all


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
resolution = 100;
time  = linspace(0,1,resolution);
T = length(time);
hbar = 1;
psi = zeros(8,T);
J = 1;
% hamiltonian is written in units of angular freq.
Bz = [0,1,100];
exp_x = zeros(1,T,3);   
exp_y = zeros(1,T,3);   
exp_z = zeros(1,T,3);   
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
psi0 = kron(up,kron(down,up));
for i = 1:3
    H = J*(Sz1.*Sz2 + Sz2.*Sz3) - Bz(i)*(Sz1+Sz2+Sz3);
    for t = 1:T
        U = expm(-1i*H*time(t));
        psi(:,t) = U*psi0;
        exp_x(1,t,i) = real(psi(:,t)' * Sx(:,:,i) * psi(:,t));
        exp_y(1,t,i) = real(psi(:,t)' * Sy(:,:,i) * psi(:,t));
        exp_z(1,t,i) = real(psi(:,t)' * Sz(:,:,i) * psi(:,t));
        S(i,t) = eee(psi(:,t));
    end
end
%% plotting
close all;
figure(1);
plot(time, S(1,:));
title('Bz = 0');
xlabel('time');
ylabel('Entropy');


figure(2);
plot(time, S(2,:));
title('Bz = 1');
xlabel('time');
ylabel('Entropy');

figure(3);
plot(time, S(3,:));
title('Bz = 100');
xlabel('time');
ylabel('Entropy');
% explanation: cuz 2.5, only one state is occupied, S = 0
% for 2.4, 8 states are occupied, therefore S is greater than 0