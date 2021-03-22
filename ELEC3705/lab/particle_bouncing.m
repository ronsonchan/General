clear all
clc

% This code calculates the time evolution of a particle (an electron, in this case) in an infinite potential well
% of width a.
% The wavefunction of the particle is expressed in the basis of the stationary states of the system, i.e. the eigenvectors of the Hamiltonian.

%% Definitions and parameters
clc
close all
clear all
hbar = 6.62606876e-34 / (2*pi); %Planck constant 
kB = 1.3806503e-23;  % Boltzmann constant
m = 9.10938188e-31;  % electron mass
T =4;
a = 1e-6;  % width of the well
x0 = a/2;   % mid-point of the well

time = 0:1e-11:1e-8;   % time range during which we observe the evolution

%% Part a)
N = 2;   % number of stationary states used for the basis

% Part a)
c_n = ones(N,1);     % This instruction makes an array of ones. 
                     % If these are the coefficients of the initial superposition, it means all coefficients are the same. 
                     % However, I still need to normalize the wavefunction.
 

    
%% Part b)
% Question is given as N, lets estimate at 100
N = 65

% For part b
p_n = zeros(1,N);
z = 0;
for i =1:N
     z = z + exp(-energy(i)/(kB*T));
end

for i =1:N
    p(i) = sqrt(exp(-energy(i)/(kB*T))/z);
end

c_n = p';

%% Part c)
% Odd numbers
N = 65
c_n = ones(N,1);
i = 1;

while (2*i+1) <= 65
    c_n(2*i) = 0;
    i = i+1;
end

%% Part d)
% Even values
N = 65
c_n = ones(N,1);
i = 1;

c_n(i) = 0;
while (2*i) <= 65
    c_n(2*i+1) = 0;
    i = i+1;
end

%% Hamiltonian and wave function
clc
[H,phi_n] = hamiltonian_WaveFunc(N);

%% Initial state
x = 0:5e-9:a;
L = length(x);              
c_n = c_n / norm(c_n);  % normalize the coefficients

psi0 = zeros(1,length(x));   % Preallocate a row vector to describe the wavefunction in real space

% Adds all wave components @ the same space together
for n = 1:N
    psi0 = psi0 + c_n(n)*phi_n(n,:);  % create the wave packet with normalized coefficients
end

figure(1)
subplot(1,2,1);
plot(x,(abs(psi0)).^2)
xlabel('x coordinate');
ylabel('|psi(0)|^2');

%% Time evolution

Ti = length(time);

% Preallocate the array that will contain on the columns the wavefunction at
% time t. This array has L columns, as many as the x-points we used to
% describe the wavefunctions, and T rows, as many points in time we are
% going to evaluate

psi_t = zeros(Ti,L);

% preallocate an array of the coefficients of the basis states, which will
% be used to express the actual wavefunctions. Note: we defined the
% coefficient vectors as column vectors, so this array must have N columns
% and T rows

c_t = zeros(N,Ti);

for i = 1:Ti
   U = expm(-1i * H * time(i) / hbar);   
   % time evolution operator
      
   c_t(:,i) = U * c_n ; % time-evolved coefficients
   psi_t(i,:) = c_t(:,i)' * phi_n;
end


%% Plot the results 

% I am using a color-scale plot to show the square of the
% wavefunction as a function of both position and time.
% MATLAB tip: the command "imagesc" does a color-plot of a function of two
% variables

figure(1)
subplot(1,2,2)
imagesc(x,time,abs(psi_t).^2);
title(strcat('a=',num2str(a)));
xlabel('position');
ylabel('time');

%% Save a figure with the results

% The instruction below saves the figure, and gives it a name that contains
% some of the paramaters of the calculation. Check the Matlab help to see
% how this was done

% hgsave(1, strcat('a_',num2str(a),'.fig'));


