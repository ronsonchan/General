%% ronson_chan_1_1
 
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
del_0 = 1e9 * h; % J 
G = 1e9 * h;
e = -1.602e-19; % J
% V_TL = 0:1e-3:1;
V_TL = 0.9; % V
sizeofVTL = length(V_TL);
% V_TR = 0:1e-3:1;
V_TR = 0.9;
sizeofVTR = length(V_TR);
V_UL = 1;
V_DL = V_UL/2;
V_UR = 1;
V_DR = V_UR/2;
dVL = V_UL - V_DL;
dVR = V_UR - V_DR;
%% calculate hamiltonian
H = hamiltonian(V_TL,V_TR,dVL,dVR);