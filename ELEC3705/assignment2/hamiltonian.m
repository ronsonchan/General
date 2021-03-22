function H = hamiltonian(V_TL, V_TR, dVL, dVR)

sx = [0,1;1,0];
sz = [1,0;0,-1];
I = eye(2);
szL = kron(sz,I);
szR = kron(I,sz);
sxL = kron(sx,I);
sxR = kron(I,sx);
h = 6.626e-34;
del_0 = 1e9 * h; % J 
G = 1e9 * h;
V_0 = 1;
e = -1.602e-19;

%% Calculate Hamiltonian

del_L = del_0*exp(V_TL/V_0);
del_R = del_0*exp(V_TR/V_0);
eps_L = e*dVL;
eps_R = e*dVR;

H_L = del_L*sxL + eps_L*szL;    % LHS contribution to the Hamiltonian
H_R = del_R*sxR + eps_R*szR;    % RHS contributin to the Hamiltonian
H_C = G* (szL.*szR) ;           % Colomb pootential contribution to Hamiltonian

H = H_L + H_R + H_C;
