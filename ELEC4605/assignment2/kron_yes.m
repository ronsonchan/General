% function that tensor products for the corresponding N qubits
% where N is the total number of matrices
% A is the operator X Y Z
function Kron = kron_yes(A,N,n)
    I = eye(2);
    dum = 1;
    for i = 1:N
        if(i == n)
            dum = kron(A,dum);
        else
            dum = kron(I,dum);
        end
    end
    Kron = dum;
end