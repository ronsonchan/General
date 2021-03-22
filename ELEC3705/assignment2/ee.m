%% This function takes a 4-dimensional ket as input, and outputs the entanglement entropy

function [S] = ee(psi)

% Basis vectors
u = [1 0; 0 1];

% density matrix
rho = psi * psi';

rho_partialB = zeros(2);

for i = 1:2
    for j = 1:2
        for n = 1:2
            rho_partialB(i,j) = rho_partialB(i,j) + kron(u(:,i)',u(:,n)') * rho * kron(u(:,j),u(:,n));
        end
    end
end

% the instruction below ensures that the reduced density matrix never has z
% zero-valued eigenvalue, which would give an error when calculating the
% logarithm
epsilon = 1e-6;
if min(eig(rho_partialB)) < epsilon 
    rho_partialB = rho_partialB + 2*epsilon*eye(2);
end

% entanglement entropy
S = -trace(rho_partialB * logm(rho_partialB));