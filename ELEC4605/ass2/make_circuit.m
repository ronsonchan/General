function U_PE = make_circuit(control_0,control_1,n)
T = n-1; % t = target qubit conter
I = n-1; % i = actual counter
U_PE = 1;
for i = I:-1:0
    for t = T:-1:0
        fprintf("%d%d\n",t,i)
    end
end
end