%%
function op = control_op(A,n,c,t)
% this function produces a controlled unitary operation, applying operator
% A to qubit t if qubit c is 1.
% A is a unitary 2x2 operator, n is the number of qubits c is the control
% qubit number and t is the target qubit
op = eye(2^n);
for i = 0:(2^n-1)
    rb = dec2bin(i,n); % convert computational basis index from decimal to binary format
    if strcmp(rb(c),'1') % control qubit is 1
        if strcmp(rb(t),'0') % target qubit is 0
            op(i+1,i+1) = A(1,1);
            cb = rb; cb(t) = '1';
            op(i+1,bin2dec(cb)+1) = A(1,2);
        else % target qubit is 1
            cb = rb; cb(t) = '0';
            op(i+1,bin2dec(cb)+1) = A(2,1);
            op(i+1,i+1) = A(2,2);
        end
    end
end