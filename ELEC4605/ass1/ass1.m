%% Assignment1
% the provided document is some typical data of some trace of the effective
% magnetic field noise. These traces are sampled with time res of 2.5 us,
% units in T, where individual time traces run along each ROW with 100
% sample traces provided.

% hence there are 100 spins in the ensemble, and each of them are
% sampled over 80,000 samples with 2.5us apart.

%% q1
clc;
close all;
clear all;

load('Assignment1_NoiseSamples_MAT', '-mat'); % samples along y  time along x
g_e = 28e9 ; % Hz/T
resolution = 2.5e-6;
fs = 1/(resolution) ; 
sample_size = length(noise);
trace_number = 100;

% take pwelch -> average over frequency
w_noise = g_e*noise ;
%Sb = zeros(trace_number, 16385); % malloc to store the pwelch
N = size(noise,2);

[Sb,f] = pwelch(w_noise', N, N/2,[],fs);

f(1) = 1e-9;
w = 2*pi*f;
%ave_Sb = zeros(1,length(Sb)); % malloc for average pwelch signals
ave_Sb = mean(Sb,2); 

% Ramsey
% define filter here
N = 100;  % resolution for t
initial_time = 0;
total_time = 5e-4;
tau = linspace(initial_time, total_time, N);
for i = 1:length(tau)
      F_Ramsey(i,:) = tau(i)^2.*(sinc(w*tau(i)/(2*pi))).^2;
%     F_Ramsey(i,:) = tau(i)^2*(sin(w*tau(i)/2)./(w*tau(i)/2)).^2;
end
% 
% figure(1)
% plot(f,F_Ramsey)
% xlabel('\omega (rad/s)'); ylabel('|F(\omega, t)|^2')
% legend('t = 0.5 s', 't = 1 s', 't = 2 s', 't = 3 s')

dw = linspace(10,fs, length(ave_Sb)); % range of w to integrate over (up to the 100rads-1)
figure(2)
loglog(dw, ave_Sb, '-x')
ylabel('S(w) (s-1)');
xlabel('frequency (Hz)')
title('noise spectral density');

% integration here
for i = 1:length(tau)
    x(i)= (1/(2*pi))*trapz(w,F_Ramsey(i,:).*ave_Sb') ; % coherence function
end
exp_y1 = exp(-x);

figure(3)
plot(tau,exp_y1,'--x')
xlabel('Tau (s)')
ylabel('<sigma_y>')
title('Ramsey decay')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T2_ramsey of order 1.31e-4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hahn echo
clc;
close all;
clear all;

load('Assignment1_NoiseSamples_MAT', '-mat'); % samples along y  time along x
g_e = 28e9 ; % rads-1/T
resolution = 2.5e-6;
fs = 1/(resolution) ; 
N = size(noise,2);
% take pwelch -> average over frequency
w_noise = g_e*noise ;
%Sb = zeros(trace_number, 16385); % malloc to store the pwelch
N = size(noise,2);

[Sb,f] = pwelch(w_noise', N, N/2,[],fs);

f(1) = 1e-9;
w = 2*pi*f;
%ave_Sb = zeros(1,length(Sb)); % malloc for average pwelch signals
ave_Sb = mean(Sb,2);

N = 100;  % resolution for t
initial_time = 1e-9;
total_time = 3e-3;
tau = linspace(initial_time, total_time, N);
for i = 1:length(tau)
    F2_Hahn(i,:) = tau(i)^2* ((sin(w*tau(i)/4)./(w*tau(i)/4)).^2) .*(sin(w*tau(i)/4)).^2;
end
% figure(1)
% plot(w,F2_Hahn)
% xlabel('\omega (rad/s)'); ylabel('|F(\omega, tau)|^2')

% integration here
% w = 2*pi*f;
for i = 1:length(tau)
   % ]A = F_Ramsey(i,:).*ave_Sb';
    x(i)= (1/(2*pi))*trapz(w,F2_Hahn(i,:).*ave_Sb') ; % decoherence function
end
exp_y1 = exp(-x);
figure(2)
plot(tau,exp_y1,'-x')
xlabel('Tau (s)')
ylabel('<sigma_y>')
title('Hahn Echo')

% T2 Hahn ~ order of -4;
%% Dynamic decoupling 
clc;
close all;
clear all;

load('Assignment1_NoiseSamples_MAT', '-mat'); % samples along y  time along x
g_e = 28e9 ; % Hz/T
resolution = 2.5e-6;
fs = 1/(resolution) ; 
sample_size = length(noise);
trace_number = 100;

w_noise = g_e*noise ;
%Sb = zeros(trace_number, 16385); % malloc to store the pwelch
N = size(noise,2);
[Sb,f] = pwelch(w_noise', N, N/2,[],fs);

f(1) = 1e-9;
w = 2*pi*f;
%ave_Sb = zeros(1,length(Sb)); % malloc for average pwelch signals
ave_Sb = mean(Sb,2);

%w = 2*pi*linspace(1e-9, fs, length(ave_Sb)); % define range of w
K = 100;  % samples for t
initial_time = 1e-9;
total_time = 1e-2;
tau = linspace(initial_time, total_time, K);
N = 1; % N here refers to number of pulses, N >= 2
% filter design

for i = 1:length(tau)
    start = 1/(2*N);
    finish = 1 - 1/(2*N);
    delta = linspace(start, finish, N);   
    sum_term = 0;   
    for k = 1:N
        sum_term = sum_term + ((-1)^k)*exp(1j*w*delta(k)*tau(i));
    end
    F2_DD(i,:) = (abs(1 + ((-1)^(1+N))*exp(1j*w*tau(i)) + 2*sum_term).^2)./((w).^2);
end
% figure(4)
% plot(w,F2_DD)
% xlabel('\omega (rad/s)'); ylabel('|F(\omega, tau)|^2')

for i = 1:length(tau)
    x(i)= (1/(2*pi))*trapz(w,F2_DD(i,:).*ave_Sb') ; % decoherence function
end
K = 100;  % samples for t
initial_time = 1e-9;
total_time = 1e-2;
tau = linspace(initial_time, total_time, K);
exp_y1 = exp(-x);
figure(3)
plot(tau,exp_y1,'-x')
xlabel('Tau (s)')
ylabel('<sigma_y>')

