%% Dynamic Decoupling
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

S = 100; % number of spin entry
N = 1; % number of Pi pulses
pi_half_length = [10e-6:8*resolution:1e-3]; % units in s
pi_half_steps = round(pi_half_length/resolution);
if(N == 1)      % Hahn echo case fo Tau/2
    pi_length = pi_half_length;
    pi_steps = pi_half_steps;
else
    pi_length = pi_half_length*2;
    pi_steps = round(pi_length/resolution);
end
wc = 1./(2*pi_length);   % range of f_centre
for s = 1:S    % loop through spins
   for tau_idx = 1:length(pi_half_steps)  % loop through each period (fixed)
        den_rho(:,:,1) = den_psi_t0;
        for i = 2:pi_half_steps(tau_idx)    % loop through the first T/2
            w = w_noise(s,i);
            U = expm(-1j*s_z*w*resolution);
            den_rho(:,:,i) = U*den_rho(:,:,i-1)*U';
        end
        last = i;
        den_rho(:,:,last) = rot_y*den_rho(:,:,last)*rot_y'; % pi pulse about Y
        for n = 1:N-1     % loops through number of pulses
            %N = 1 does not behave like hahn because 
            %pi_steps is different than pi_half_steps
           %copy_pi = rot_y*den_rho(:,:,last)*rot_y'; % use this to check what rot_y does
            %fprintf("In Here\n");
            for i = 1:pi_steps(tau_idx)    % loops through T
                w = w_noise(s,i+last);
                U = expm(-1j*s_z*w*resolution);
                den_rho(:,:,i+last) = U*den_rho(:,:,i+last-1)*U';
            end
            last = last + i;
            den_rho(:,:,last) = rot_y*den_rho(:,:,last)*rot_y';
        end
        for i = 1:pi_half_steps(tau_idx)
            w = w_noise(s,i+last);
            U = expm(-1j*s_z*w*resolution);
            den_rho(:,:,i+last) = U*den_rho(:,:,i+last-1)*U';
        end
        last = last + i;
       exp_y(s,tau_idx) = trace(den_rho(:,:,last)*s_y);
    end
    s
end
[rol, col] = size(exp_y);
for c = 1:col   
    ave_exp_y(c) = mean(exp_y(:,c));  
end
%time_scale = linspace(10e-6, 2*i*N*resolution, length(pi_half_steps));
time_scale = linspace(10e-6,2*pi_half_length(end)+pi_length(end)*(N-1), col);
plot(time_scale, ave_exp_y, '--x');