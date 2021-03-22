%% ronson_chan_1_2

%% definition
clc;
close all
clear all

up = [1;0];
down = [0;1];
upup = [1;0;0;0];
updown = [0;1;0;0];
downup = [0;0;1;0];
downdown = [0;0;0;1];
h = 6.626e-34;
hbar = 6.626e-34 / (2*pi);
sx =  [0,1;1,0];
sy =  [0,1i;-1i,0];
sz =  [1,0;0,-1];
I = eye(2);
szL = kron(sz,I);
szR = kron(I,sz);
sxL = kron(sx,I);
sxR = kron(I,sx);
V_0 = 1;    % 1V
G = 1e9 * h;
V_TL = 0; % V
V_TR = 0;
V_UL = 1;
V_DL = V_UL/2;
V_UR = 1;
V_DR = V_UR/2;
dVL = 0;
dVR = 0;
%% Calculate Hamiltonian
clc 
H = hamiltonian(V_TL, V_TR, dVL, dVR);
[V,D] = eig(H);
exp_sxL = zeros(4,1);
exp_sxR = zeros(4,1);
exp_szL = zeros(4,1);
exp_szR = zeros(4,1);
for i = 1:4
    exp_sxL(i) = V(:,i)' * sxL * V(:,i);
    exp_sxR(i) = V(:,i)' * sxR * V(:,i);
    exp_szL(i) = V(:,i)' * szL * V(:,i);
    exp_szR(i) = V(:,i)' * szR * V(:,i);
end

%% plotting
clc;
close all

subplot(2,1,1)
hold on
plot(exp_sxL);
plot(exp_szL);
xlabel('state vectors');
ylabel('expectation value');
title('expectation value vs state vector in LHS');
legend('exp sxL', 'exp szL');
hold off

subplot(2,1,2)
hold on
plot(exp_sxR);
plot(exp_szR);
xlabel('state vectors');
ylabel('expectation value');
title('expectation value vs state vector in RHS');
legend('exp sxR', 'exp szR');
hold off

