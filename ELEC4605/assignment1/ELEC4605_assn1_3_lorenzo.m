%% Setup
clear; clc; close all;

%% Setup and Constants
load('Assignment1_NoiseSamples_MAT');
hbar = 6.62606876e-34 / (2*pi); 
gamma_e = 28e9;
t_res = 2.5e-6;
fs = 1/t_res;

%% Noise
[pxx,f] = pwelch(gamma_e.*noise',length(noise(1,:)),length(noise(1,:))/2,[],fs);

S_b = mean(pxx,2);

f(1) = 1e-9;
w = 2*pi.*f';

%% Pulse Variables

N = 30;
t_res = 2.5e-6;
tau = t_res;

t_max = N*tau;
t = 0:t_res:t_max;

%%
tot_prev = 0;
figure;
hold on;
for i = 1:N
    if i == 1
    tic;
    end
    [F, tot_prev] = filter_N(w,t_max,i, tot_prev);
%     F = filter_N(w,t_max,i);
    
    if i == 1
        dur = toc;
    end
    clc;
    fprintf("Generating Filter...\nProgress: %d%%\nTime Remaining: %f\n",ceil(i/N*100), dur*(N-i));
    
    plot(f,F, 'DisplayName', 'N = ' + string(i))
%     drawnow
%     legend;
    
end
% plot(f,F);
%xlim([0 5e4])
hold off

%% Helper Functions

% Input:
%     w = Centre Frequency
%     tau = Pulse width
%     N = Number  of Pulses
% Output:
%     F = Dynamical decoupling filter function with duration t and N pulses.
function [F,tot] = filter_N(w,tau,N, tot_prev)
    t_max = tau*N;
    
    offset = 1/(2*N);
    delta = linspace(offset, 1-offset,N);
    
%     tot = tot_prev;
    tot = tot_prev +((-1)^N.*(exp(1i*w*t_max.*delta(end))));

    F = abs(1+(-1)^(1+N)*exp(1i.*w.*t_max)+2*tot).^2 ./(w.^2);
end