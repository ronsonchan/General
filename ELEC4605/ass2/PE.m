function U_PE = PE(x,beta)
    X = diag(x);
    X2 = X^2;
    U_PE = expm(-1i*beta*X2);
end

