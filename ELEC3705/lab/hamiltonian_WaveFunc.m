function [H,phi_n] = hamiltonian_WaveFunc(N)
% Hamiltonian and wavefunctions

hbar = 6.62606876e-34 / (2*pi); %Planck constant 
m = 9.10938188e-31;  % electron mass
a = 1e-6;  % width of the well
% Calculate the eigenenergies of the bound states in the square well

for n = 1:N
    En(n) = (n*pi*hbar)^2 / (2*m*a^2);
end

% Construct the Hamiltonian, truncated according the the temperature

H = zeros(N,N);

x = 0:5e-9:a;
L = length(x);
phi_n = zeros(N,L);   % array with the space wavefunctions on the rows

% Creates hamiltonian, setting diagonal as eigenvalues
% Creates wave function 'n' in respect to position
for n = 1:N
    H(n,n) = En(n);
    for i = 1:L
    phi_n(n,i) = sqrt(2/a) * sin(n*pi*x(i)/a);
    end
end





end

