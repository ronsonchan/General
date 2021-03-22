function [S] = eee(ket)
    % Basis vectors
    u = [1 0 ; 0 1];
    psi = ket;
    % density matrix
    rho = psi * psi';

    rho_partial = zeros(2);

    for i = 1:2
        for j = 1:2
            for n = 1:2
                for p = 1:2
                    rho_partial(i,j) = rho_partial(i,j) + kron(u(:,i)',kron(u(:,n)',u(:,p)'))  * rho * kron(kron(u(:,j),u(:,n)),u(:,p));
                end
            end
        end
    end
    
    epsilon = 1e-6;
    if min(eig(rho_partial)) < epsilon 
        rho_partial = rho_partial + 2*epsilon*eye(2);
    end

    % entanglement entropy
    S = -trace(rho_partial * logm(rho_partial));
end