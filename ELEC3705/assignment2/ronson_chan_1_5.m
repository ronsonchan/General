%% ronson_chan_1_5
clc
% V_TR = V_TL = -5;
% dVL = 5
% dVR = -5;

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
V_0 = 1;    % 1V
G = 1e9 * h;
res = 100;
time = linspace(0,2*pi, res); % time vector
T = length(time);
psi = zeros(4,T);
psi(:,1) = updown;
% sizeV = length(V_TL);
% set up dVL and dVR to oscillate between 1,-1

% dVL = 1
% dVR = -1;
% sizedV = length(dVL);
%% hamiltonian
clc;
close all;

V_TL = -5; % V
V_TR = -5;
dVL = 5*zeros(1,T);
dVR = -5*zeros(1,T);

eVec = zeros(4,T,4);
s_zL_evo = zeros(1,T);
s_zR_evo = zeros(1,T);
S = zeros(1,T);
for t = 1:T
    H = hamiltonian(V_TL, V_TR, dVL(t), dVR(t));
    [V,D] = eig(H); % get eigenvector
    U = expm(-1i*H*time(t)/hbar); % time evlotion operator
    psi(:,t) = U*psi(:,1);
    s_zL_evo(1,t) = psi(:,t)' * szL * psi(:,t);
    s_zR_evo(1,t) = psi(:,t)' * szR * psi(:,t);
    S(1,t) = ee(psi(:,t));
end
figure(1)
hold on
plot(s_zL_evo(1,:))
plot(s_zR_evo(1,:))
legend('expectation value of LHS','expectation value of RHS');
title('expeectation value vs time');
xlabel('time,sec');
ylabel('expectation value');
hold off

figure(2)
plot(time,S)
title('Entropy vs time')
xlabel('time, sec');
ylabel('entropy');
