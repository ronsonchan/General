%% QFT_generalized
function U_qft = QFT_generalized(n)
% Standard gates
X = [0 1; 1 0]; % Pauli X
H = (1/sqrt(2))*[1 1; 1 -1]; % Hadamard
H_n = zeros(2^n,2^n,n); % malloc for Hadamard applied to jth qubit in hilbert space of n qubit
Rot = zeros(2,2,n-1); % the last col of R stores the "magnitude" of how far it needs to be rotated
for i = 1:n
    H_n(:,:,i) = tensor_op(H,n,i); % Hadamard applied to qubit 1 in the combined Hilbert space of 3 qubits
    Rot(:,:,i) = R(i);  % storing nth R-gate
end
qft = 1;
a = 1; % counter for target qubit
L = 0; % counter for control qubit
for  i = 1:n-1  % i counts for number of qubits up to n-1
    Q = 1;
    qft = H_n(:,:,i)*qft;
    for j = 2:n-L % count rotation gate, truncates as you go down the line
        Q = control_op(Rot(:,:,j),n,j+L,a)*Q;
    end
    a = a+1; % increment the ath target qubit
    L = L+1; % truncate number of gates as you go down the circuit line
    qft = Q*qft; 
end
qft = H_n(:,:,end)*qft; % sandwich with the final (last) hadamard gate Hn
% implement swapgate here
SWAP = zeros(2^n,2^n,floor(n/2));
L = 0;
for i = 1:floor(n/2)
    SWAP(:,:,i) = control_op(X,n,i,n-L)*control_op(X,n,n-L,i)*control_op(X,n,i,n-L);
    L = L + 1;
end
[~,~,C] = size(SWAP);
prod = 1;
for i = 1:C
    prod = SWAP(:,:,i)*prod;
end
U_qft = prod*qft;
end