function U_KE = KE(k,alpha)
    hbar = (6.626e-34)/(2*pi);
    K = diag(k);
    K2 = K^2;
    P2 = hbar^2 .* K2;
    U_KE = expm(-1i*alpha*P2);
end

