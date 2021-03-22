%% lab6

%% Q1
% intuitive calculation done on paper

%% Q2
clc
clear all
close all

x_end = 1e-6;
x_step = 1e-9;
x = 0:x_step:x_end;
sigma = 50e-9;
X = length(x);
freq = 1/x_step
dx = 1e-11;  % (k/x=0 = inf)
Fs = 2*freq; % sampling frequency
k = linspace(0,Fs, X);

% resolution = 1001. same as res of x
% range = [10^6, 10^11]

% x(i)*k(i) = 2*pi for all i

psi = (1/sqrt(2*pi*sigma^2))*exp((-(x-x(500)).^2)/(2*sigma^2));
figure(1)
PSI = fft(psi);
PSI_SHIFT = fftshift(PSI);
plot(abs(PSI_SHIFT))
sigma_k = 1/sigma_x;
%% Q3
% solution of type
% psi = A*exp(-alpha*x)
% where alpha = sqrt((2m*(V-E))/hbar^2) is the transition prob

%% Task
clc;
clear all
close all

N = 400;
len = 100;
k0 = 2; % radm-1
sigma_k = 0.2; % radm-1
x = linspace(0,len,N);
k = linspace(-2*pi,2*pi, N);
x0 = 25;
T = 30;
dt = 0.5;
time = 0:dt:T;
% construct energy barrier at [55, 65]m 
height = 1;  % where you can change the height
EB = height*(heaviside(x-55) - heaviside(x-65)); 

% construct psi in frequency domain
psi_k = zeros(1,N);
for i = 1:N
    psi_k(i) = (1/sqrt(2*pi*sigma_k^2))*exp(-((k(i)-k0)^2)/(2*sigma_k^2));
end
% coulomb potential
V = diag(EB); 

% construct Hamiltonian
m = 1; % using natural units
hbar = 1; 
dx = x(2) - x(1);
L = (hbar^2)/(2*m);
H = zeros(N+1, N+1);
e = ones(N,1);
S = spdiags([e -2*e e],-1:1 ,N,N);
H = full(S);
laplacian = (1/dx^2).*H;
H = V + L*laplacian;

% initialize ket vector

psi_x = ifft(psi_k);
psi_x = circshift(psi_x, min(find(x > 25)));
hold on
subplot(3,1,1)
plot(x, abs(psi_x));
subplot(3,1,2)
plot(x,EB);
hold off;

% time evolution

psi = zeros(N,length(time));
psi(:,1) = psi_x;
U = expm(-1i * H * dt);



for i = 2:length(time)
    psi(:,i) = U * psi(:,i-1);
    subplot(3,1,3)
    plot(x, real(abs(psi(:,i))))
    hold on
end
hold off