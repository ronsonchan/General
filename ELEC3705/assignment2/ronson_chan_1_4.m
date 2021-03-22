%% ronson_chan_1_4

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
V_TL = -5:0.01:5; % V
V_TR = -5:0.01:5;
sizeV = length(V_TL);
V_UL = 1;
V_DL = V_UL/2;
V_UR = 1;
V_DR = V_UR/2;
dVL = -5:0.01:5;
dVR = -5:0.01:5;
sizedV = length(dVL);
%% Calculate Hamiltonian
clc 
exp_szL = zeros(4,sizeV,4);
exp_szR = zeros(4,sizeV,4);
S = zeros(4,sizeV,4);
% run a for loop 
% run 4 hamiltonian, each with 1 variable and 3 controls
% to find out the maximum contribution to entropy
for i = 1:sizeV
    H = hamiltonian(V_TL(i), 0, 0, 0);
    [V,D] = eig(H);
    for j = 1:length(V)
       S(1,i,j) = ee(V(:,j)); 
       exp_szL(1,i,j) = real(V(:,j)' * szL * V(:,j));
       exp_szR(1,i,j) = real(V(:,j)' * szR * V(:,j));
    end
    
    H = hamiltonian(0, V_TR(i), 0, 0);
    [V,D] = eig(H);
    for j = 1:length(V)
       S(2,i,j) = ee(V(:,j)); 
       exp_szL(2,i,j) = real(V(:,j)' * szL * V(:,j));
       exp_szR(2,i,j) = real(V(:,j)' * szR * V(:,j));
    end
    
    H = hamiltonian(0, 0, dVL(i), 0);
    [V,D] = eig(H);
    for j = 1:length(V)
       S(3,i,j) = ee(V(:,j)); 
       exp_szL(3,i,j) = real(V(:,j)' * szL * V(:,j));
       exp_szR(3,i,j) = real(V(:,j)' * szR * V(:,j));
    end
    
    H = hamiltonian(0, 0, 0,dVR(i));
    [V,D] = eig(H);
    for j = 1:length(V)
       S(4,i,j) = ee(V(:,j)); 
       exp_szL(4,i,j) = real(V(:,j)' * szL * V(:,j));
       exp_szR(4,i,j) = real(V(:,j)' * szR * V(:,j));
    end
end
%% plotting
figure(1);
hold on
subplot(4,1,1)
plot(V_TR, S(1,:,1))
xlabel('VT');
ylabel('Entropy');
title('Entropy vs VT');
subplot(4,1,2)
plot(V_TR, S(2,:,1))
xlabel('V_T');
ylabel('Entropy');
title('Entropy vs VT');
subplot(4,1,3)
plot(V_TR, S(3,:,1))
xlabel('V_T');
ylabel('Entropy');
title('Entropy vs VT');
subplot(4,1,4)
plot(V_TR, S(4,:,1))
xlabel('V_T');
ylabel('Entropy');
title('Entropy vs VT');
hold off

figure(2)
hold on
subplot(4,1,1)
plot(V_TR, exp_szL(1,:,1))
xlabel('V_T');
ylabel('Exp value, LHS');
title('Expectation value of Z,LHS vs V_T');
subplot(4,1,2)
plot(V_TR, exp_szL(2,:,1))
xlabel('V_T');
ylabel('Exp value, LHS');
title('Expectation value of Z,LHS vs V_T');
subplot(4,1,3)
plot(V_TR, exp_szL(3,:,1))
xlabel('V_T');
ylabel('Exp value, LHS');
title('Expectation value of Z,LHS vs V_T');
subplot(4,1,4)
plot(V_TR, exp_szL(4,:,1))
xlabel('V_T');
ylabel('Exp value, LHS');
title('Expectation value of Z,LHS vs V_T');
hold off

figure(3)
hold on
subplot(4,1,1)
plot(V_TR, exp_szR(1,:,1))
title('Expectation value of Z, RHS vs V_T');
xlabel('V_T');
ylabel('Exp value, RHS');
subplot(4,1,2)
plot(V_TR, exp_szR(2,:,1))
title('Expectation value of Z, RHS vs V_T');
xlabel('V_T');
ylabel('Exp value, RHS');
subplot(4,1,3)
plot(V_TR, exp_szR(3,:,1))
title('Expectation value of Z, RHS vs V_T');
xlabel('V_T');
ylabel('Exp value, RHS');
subplot(4,1,4)
plot(V_TR, exp_szR(4,:,1))
title('Expectation value of Z, RHS vs V_T');
xlabel('V_T');
ylabel('Exp value, RHS');
hold off

%% answer
%% a
% from the above graph we can see that
% for maximum entropy:
% V_TR = V_TL = -5V
% dVL = dVR = 0;

%% b 
% for both electron to be in the upper dot
% V_TR = V_TL = -5
% dVL = 5
% dVR = 5

%% c
% V_TR = V_TL = -5;
% dVL = 5
% dVR = -5;