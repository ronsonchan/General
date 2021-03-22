%% Entropy
%%
clear all
clc

%% States and density matrices

% Basis vectors
u = [1 0; 0 1];

psi_A = [0.6;0.8];
psi_B = [1;0];

rho_A = psi_A * psi_A'
rho_B = psi_B * psi_B'

rho_C = kron(rho_A,rho_B)

%% Partial trace

rho_partialB = zeros(2);

for i = 1:2
    for j = 1:2
        for n = 1:2
            rho_partialB(i,j) = rho_partialB(i,j) + kron(u(:,i)',u(:,n)') * rho_C * kron(u(:,j),u(:,n));
        end
    end
end

rho_partialB

rho_partialA = zeros(2);

for i = 1:2
    for j = 1:2
        for n = 1:2
            rho_partialA(i,j) = rho_partialA(i,j) + kron(u(:,n)',u(:,i)') * rho_C * kron(u(:,n),u(:,j));
        end
    end
end

rho_partialA

%% entanglement entropy of product state

S = -trace(rho_partialB * logm(rho_partialB))
% note: if I take logm(rho_partialA), Matlab will return an error, because that matrix has a zero eigenvalue

%% entanglement entropy of entangled state

singlet = (1/sqrt(2)) * [0;1;-1;0];

rho_singlet = singlet * singlet'

rho_partialA_singlet = zeros(2);

for i = 1:2
    for j = 1:2
        for n = 1:2
            rho_partialA_singlet(i,j) = rho_partialA_singlet(i,j) + kron(u(:,n)',u(:,i)') * rho_singlet * kron(u(:,n),u(:,j));
        end
    end
end

rho_partialA_singlet

S = -trace(rho_partialA_singlet * logm(rho_partialA_singlet))