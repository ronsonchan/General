function op = tensor_op(A,n,p)
% This function produces an operator acting in a larger Hilbert space.
% A is a 2x2 unitary, p is the qubit number on which A is to operate
% and n is the number of qubits in the system.
% Example: tensor_op(A,2,3) generates the operator U = I x A x I (where x is
% the tensor product)
I = eye(2,2);
if p == 1
    op = A;
else
    op = I;
end
for i = 1:n-1
    if i+1 == p
        op2 = A;
    else
        op2 = I;
    end
    %op = tensor(op,op2);
    op = kron(op,op2);
end
end