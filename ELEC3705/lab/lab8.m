%% lab8

%% Q1

% if V = 0, that means G = 0 for all n, then the hamiltonian 
% <n_H_m> = hbar^2 * k_alpha^2 / 2m0 iff n = m
% then it is just a free electron travelling in freespace
% with discrete velocities.

%% Task
clc
clear all
close all

size  = 25;
a = 400e-12;
beta = 0; % beta can go up to inf
n = floor(-size/2):floor(size/2);
alpha = 1:size+1;  
N = 25;
k0 = 2*pi / (N*a);
k_a = linspace(-pi/a,pi/a,1000);  % crystal k number
G = 2*pi*n/a; % "k" in frequency doamin
hbar = 1; % using natural units
m = 1; % mass of electron, which is 1 in natural units
V = zeros(1,size); % malloc a V vector in FREQUENCY domain
% for counter = 1:size+1
   V = (hbar^2/(2*m)) *(beta*k0)^2 ; 
% end

% making the hamiltonian
% note that this hamiltonian is in fourier domain
H = zeros(size, size);
eigenvalues = zeros(size,length(k_a));
H = V*ones(size,size);
for j = 1:length(k_a)
    for i = 1:size
        H(i,i) = (hbar^2 / (2*m))*(k_a(j) + G(i))^2;
    end
    [eigVec,eVal]= eig(H);
    eigenvalues(:,j)= diag(eVal);
end

plot(eigenvalues')