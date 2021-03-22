%% coupledspin
clear all;
clc;

%% Pauli matrices in the basis |+>, |->

sigma_x = [0 1; 1 0];
sigma_y = [0 -1i; 1i 0];
sigma_z = [1 0; 0 -1];

%% two-spin matrices in the basis |++>,|+->,|-+>,|-->
sigma_x1 = kron(sigma_x,eye(2));
sigma_x2 = kron(eye(2),sigma_x);

sigma_y1 = kron(sigma_y,eye(2));
sigma_y2 = kron(eye(2),sigma_y);

sigma_z1 = kron(sigma_z,eye(2));
sigma_z2 = kron(eye(2),sigma_z);

%% Hamiltonian and eigenstates

J = 1;   %exchange interaction constant

H = J*(sigma_x1 * sigma_x2 + sigma_y1 * sigma_y2 + sigma_z1 * sigma_z2);

[V,lambda] = eig(H)

S = V(:,1);
Tplus = V(:,2);
T0 = V(:,3);
Tminus = V(:,4);

%% Spin expectation values

sigma_ztot = sigma_z1 + sigma_z2;

S_ztot = S' * sigma_ztot * S
T0_ztot = T0' * sigma_ztot * T0

sigma2 = (sigma_x1 + sigma_x2)^2 + (sigma_y1 + sigma_y2)^2 + (sigma_z1 + sigma_z2)^2

S_tot = S' * sigma2 * S
T0_tot = T0' * sigma2 * T0