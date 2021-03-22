%% partc
%% fix time interval between pi
clc;
close all;
clear all;

load('Assignment1_NoiseSamples_MAT', '-mat'); % samples along y  time along x
g_e = 28e9 ; % rads-1/T
resolution = 2.5e-6;
fs = 1/(resolution) ; 
sample_size = length(noise);
trace_number = 100;
hbar = 6.626e-34 / (2*pi);
theta = pi/2;
s_x = [0,1;1,0];
Sx = (hbar/2)*s_x;
s_z = [1,0;0,-1];
Sz = (hbar/2)*s_z;
s_y = [0,-1j;1j,0];
Sy = (hbar/2)*s_y;

rot_x_half = [cos(theta/2), -1j*sin(theta/2);
              -1j*sin(theta/2), cos(theta/2)]; % rotational matrix along x with pi/2
psi_t0 = rot_x_half*[0;1];     
den_psi_t0 = psi_t0 * (psi_t0)';
w_noise = g_e*noise;
rot_x = -1j*s_x;
rot_y = -1j*s_y;

% idea: fix the time between PI-PULSES
% and acording to N, add the sequence after sequence
pi_half_length = 200e-6;
pi_half_steps = pi_half_length/resolution;
pi_length = pi_half_length*2;
pi_steps = pi_length/resolution;
fc = 1/(2*pi_length)
wc = 2*pi*fc;
N = [1:40];
S = 100;
for s = 1:S
    for p = 1:length(N)
        den_rho(:,:,1) = den_psi_t0;
        for i = 2:pi_half_steps
            w = w_noise(s,i);
            U = expm(-1j*s_z*w*resolution);
            den_rho(:,:,i) = U*den_rho(:,:,i-1)*U';
        end
        last = i;
        % pi pulse
        den_rho(:,:,last) = rot_y*den_rho(:,:,last)*rot_y';
        for n = 1:(N(p)-1)
            for i = 1:pi_steps
                w = w_noise(s,i+last);
                U = expm(-1j*s_z*w*resolution);
                den_rho(:,:,i+last) = U*den_rho(:,:,i+last-1)*U';
            end
            last = last + i;
            den_rho(:,:,last) = rot_y*den_rho(:,:,last)*rot_y';
        end
        for i = 1:pi_half_steps
            w = w_noise(s,i+last);
            U = expm(-1j*s_z*w*resolution);
            den_rho(:,:,i+last) = U*den_rho(:,:,i+last-1)*U';
        end
        last = last + i;
        exp_y(s,p) = trace(den_rho(:,:,last)*s_y);
    end
    s
end
[rol,col] = size(exp_y);
for a = 1:col
   ave(a) = mean(exp_y(:,a)); 
end
tot_tau = linspace(pi_half_length, N(end)*pi_length + 2*pi_half_length, col);
figure(1)
plot(tot_tau,ave)
xlabel("Total Time Tau (s)")
ylabel("<sy>");
title("decoherence at fc = 5000 Hz")
%%
tau = tot_tau;
w = linspace(600, 1e5, length(tau));
N = 4;
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

figure(4)
plot(w,F2_DD)
xlabel('\omega (rad/s)'); ylabel('|F(\omega, tau)|^2')
