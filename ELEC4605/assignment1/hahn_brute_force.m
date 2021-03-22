%% hahn echo brute force
clc;
close all;
clear all;

load('Assignment1_NoiseSamples_MAT', '-mat'); % samples along y  time along x
g_e = 28e9 ; % Hz/T
resolution = 2.5e-6;    % units in seconds
fs = 1/(resolution) ;   % units in Hz
hbar = 6.626e-34/(2*pi);
theta = pi/2;
s_x = [0,1;1,0];
s_z = [1,0;0,-1];
Sz = (hbar/2)*s_z;
s_y = [0,-1j;1j,0];

rot_x_half = [cos(theta/2), -1j*sin(theta/2);
              -1j*sin(theta/2), cos(theta/2)]; % rotational matrix along x with pi/2
psi_t0 = rot_x_half*[0;1];     
den_psi_t0 = psi_t0 * (psi_t0)';
rot_x = -1j*s_x;
w_noise = g_e*noise;

rot_x_half = [cos(theta/2), -1j*sin(theta/2);
              -1j*sin(theta/2), cos(theta/2)]; % rotational matrix along x with pi/2
psi_t0 = rot_x_half*[0;1];     
den_psi_t0 = psi_t0 * (psi_t0)';
w_noise = g_e*noise;

S = 100;
T = 50;
t_start = 1e-5;
t_end = 1.5e-3;
tau_half_length = linspace(t_start, t_end, T);
tau_half_steps = fix(tau_half_length/resolution);
for s = 1:S
    for tauNum = 1:length(tau_half_length)
        phase = 0; 
        for t = 1:tau_steps(tauNum)
           w = w_noise(s,t);
           phase = phase + w*resolution;
           U = expm(-1j*Sz*phase/hbar);
           den_rho(:,:,t) = U*den_psi_t0*U';
        end
        exp_y(s,tauNum) = trace(den_rho(:,:,t)*s_y);
    end
end
[rol, col] = size(exp_y);
% for i = 1:col
%    ave_exp_y(i) = mean(exp_y(:,i)); 
% end
ave_exp_y= mean(exp_y);
figure(1)
plot(tau_half_length,ave_exp_y,'-x');
xlabel('Tau (s)');
title('Ramsey, Brute Force');