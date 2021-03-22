function [meas_X, meas_Y, meas_Z] = meas(n,A,QQQ)
    %tensor to appropriate dimension
    % malloc for tensors for n qubits
    
    A = upper(A);
    
    if mean(ismember(A,['X','Y','Z'])) < 1
        error('Specified axis does not exist')
    end
    
%     %% Measure X
%     
    if contains(A,'X')
        X = [0 1;1 0];
        meas_X = zeros(1,n);
        %meas_X = zeros(2^n,n);
        %rho_QQQ = QQQ*QQQ'; % density matrix for last state
        for i = 1:n
            tensor = kron_yes(X,n,i); % create the ith tensor of X
            QQQ'*tensor*QQQ
            format long
            meas_X(i) = real(QQQ'*tensor*QQQ);
            %meas_X(:,i) = real(QQQ'*tensor*QQQ);
        end        
    end
%     
%     %% Measure Y
%     
    if contains(A,'Y')
        Y = [0 -1j;1j 0];
        meas_Y = zeros(2^n,n);
        %rho_QQQ = QQQ*QQQ'; % density matrix for last state
        for i = 1:n
            tensor = kron_yes(Y,n,i); % create the ith tensor of X
            meas_Y(:,i) = real(tensor*QQQ);
        end     
    end
%     
% Measure Z
    if contains(A,'Z')
        Z = [1 0;0 -1];
        meas_Z = zeros(2^n,n);
        %rho_QQQ = QQQ*QQQ'; % density matrix for last state
        for i = 1:n
            tensor = kron_yes(Z,n,i); % create the ith tensor of X
            meas_Z(:,i) = real(tensor*QQQ);
        end    
    end
end