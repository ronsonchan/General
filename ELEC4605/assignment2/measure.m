% n = number of qubits
function Meas = measure(n,final_state)
    Meas = zeros(n,1);
    I = eye(2^n); % stores computational basis state
    for i = 1:2^n
        Meas(i) = abs(I(:,i)' * final_state).^2;
    end
end