function ccu = ccu(U,n)
% this function creates a controlled unitary with two controls (hard-coded
% here as qubits 1 and 2). The circuit is described in the lecture notes.
% The unitary acts on qubit 3. Try make this a general function (where you can select 
% two qubits out of n as the controls and apply U to a qubit of your choice).
X = [0 1; 1 0];
V = U^0.5;
ccu = control_op(V,n,1,3)*control_op(X,n,1,2)*control_op(V',n,2,3)*control_op(X,n,1,2)*control_op(V,n,2,3);
end