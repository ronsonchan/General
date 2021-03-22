%% lab3

%% Q1
clc
clear all
close all
h_bar = 6.626e-34 / (2*pi);
m0 = 9.11e-31;
kb = 1.38e-23;
n = [1:500];
a = 1e-6;
k = n*pi/a;
t_length = 10
psi0 = @(x) sqrt(2/a)*sin(k*x);
E = zeros(1, length(n));
for count = 1:length(n)
    E(count) = (k(count)*h_bar)^2 / (2*m0);
end
Z = zeros(1,length(n));
Z(1) = 0;
for count = 1:length(n)
   Z(count) = exp(-E(count)/(kb*t_length)) ;
end
Z = sum(Z, 'all')
P = zeros(1,length(n));
for count = 1:length(n)
   P(count) = (1/Z)*exp(-E(count)/(kb*t_length));
end
Num = 0;
for count  = 1:length(n)
   if((P(count)/P(1)) < 0.01)
      Num = count; 
      break
   end
end

%% define parameters
clear all;
close all;
clc

hbar = 6.626e-34 / (2*pi);
m0 = 9.1094e-31;
kb = 1.3806e-23;
a = 1e-6;
time = 0:1e-11:1e-8;
T = length(time);   % length of time vector
x = 0:5e-9:a;
length = length(x); % length of x
Temp = 4;

%% a
N = 2;
U = zeros(length, length, T);
H = zeros(N, N);
for n = 1:N
    k = n*pi/a;
    E = (k*hbar)^2 / (2*m0);
    H(n,n) = E;
    % define wavefunction here
    for i = 1:length
       phi(n,i) = (2/sqrt(a)) * sin(k*x(i)); 
    end
end

c_n = ones(N, 1);  % malloc for coef vector 

%% part b
N = 65;
p_n = zeros(N, 1);
% find energy
E = zeros(1, N);
H = zeros(N,N);
for n = 1:N
    k = n*pi/a;
    E(n) = (k*hbar)^2 / (2*m0); 
    H(n,n) = E(n);
    for i = 1:length
       phi(n,i) = (2/sqrt(a))* sin(k*x(i));
    end
end
Z = 0;
for n = 1:N
    Z = Z + exp(-E(n)/(kb*Temp));
end
for n = 1:N
    p_n(n) = (1/Z)*exp(-E(n)/(kb*T));
end
c_n = p_n;

%% part c
N = 65;
p_n = zeros(N,1);
E = zeros(N,1);
H = zeros(N,N);
% find energy
c_n = ones(N,1);
phi = zeros(N,1);
for n = 1:N
    k = n*pi/a;
    E(n) = (k*hbar)^2/(2*m0);
    H(n,n) = E(n);
    for i = 1:length
        phi(n,i) = (sqrt(2/a))*sin(k*x(i));
    end
end
n = 1;
for n = 1:N
    if( mod(n,2) == 0)
        c_n(n) = 0;
    end
end
%% task d
N = 65;
p_n = zeros(N,1);
E = zeros(N,1);
H = zeros(N,N);
% find energy
c_n = ones(N,1);
phi = zeros(N,1);
for n = 1:N
    k = n*pi/a;
    E(n) = (k*hbar)^2/(2*m0);
    H(n,n) = E(n);
    for i = 1:length
        phi(n,i) = (sqrt(2/a))*sin(k*x(i));
    end
end
n = 1;
for n = 1:N
    if( mod(n,2) == 1)
        c_n(n) = 0;
    end
end
%% initial state def
clc;
c_n = c_n/(norm(c_n));
psi_0 = zeros(1, length);

for i = 1:N
   psi_0 = psi_0 + c_n(i)*phi(i,:);  
end
figure(1)
subplot(1,2,1);
plot(x,(abs(psi_0)).^2)
xlabel('x coordinate');
ylabel('|psi(0)|^2');

%% Time evolution
clc
ct = zeros(N, length);
psi = zeros(T, length);
for i = 1:T
    % time evolution operator U
   U = expm(-1i*H*time(i)/hbar);
   ct(:,i) = U*c_n;
   psi(i,:) = ct(:,i)' *phi;
end

figure(1)
subplot(1,2,2)
imagesc(x,time,abs(psi).^2);
title(strcat('a=',num2str(a)));
xlabel('position');
ylabel('time');