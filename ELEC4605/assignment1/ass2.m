%% assignment part b
%% Ramsey
clc;
close all;
clear all;

load('Assignment1_NoiseSamples_MAT', '-mat'); % samples along y  time along x
g_e = 28e9 ; % Hz/T
resolution = 2.5e-6;
fs = 1/(resolution); 
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

w_noise = g_e*noise;

S = 100;
N = 100;
t_start = 1e-5;
t_end = 2e-3;
TAU = linspace(t_start, t_end, N);
for s = 1:S    % loop tthrough spin
    for t = 1:N   %  loop through tau
        phase = 0;
        T = TAU(t)/resolution;
        tt = T;  % step through different tau value 
        for i = 1:tt
            w = w_noise(s,i);
            phase  = phase + w*resolution;  % phase is in radian cuz it's an angle
            U = expm(-1j*Sz*phase/hbar);
            den_rho(:,:,i) = U*den_psi_t0*U';
        end
        exp_y(s,t) = trace(den_rho(:,:,i)*s_y);
    end
end
ave_exp_y = mean(exp_y);
plot(TAU, ave_exp_y, '-x')
xlabel('Tau, (s)')
ylabel('<sigma_y>')
title('Ramsey brute force');
%% Hahn

clc;
close all;
clear all;
hold on;
load('Assignment1_NoiseSamples_MAT', '-mat'); % samples along y  time along x
g_e = 28e9 ; % Hz/T
resolution = 2.5e-6;
fs = 1/(resolution); 
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
w_noise = g_e*noise;
rot_x = -1j*s_x;
S = 100;
t_start = 1e-5;
t_end = 5e-3;
TAU = [4*t_start:40*resolution:t_end];
for s = 1:S    % loop tthrough spin
    for t = 1:length(TAU)   %  loop through tau
        pi_flag = 0;
        T = round(TAU(t)/resolution);
        tt = T;  % step through different tau value 
        den_rho = zeros(2,2,T);
        den_rho(:,:,1) = den_psi_t0;
        for i = 2:tt
            w = w_noise(s,i);       
            if(pi_flag == 0 && i > (floor(TAU(t)/resolution)/2))
                den_rho(:,:,i-1) = rot_x*den_rho(:,:,i-1)*rot_x;
                pi_flag = 1;
            end
            U = expm(-1j*s_z*w*resolution);
            den_rho(:,:,i) = U*den_rho(:,:,i-1)*U';
        end
        exp_y(s,t) = trace(den_rho(:,:,end)*s_y);
    end
    s
end
ave_exp_y = mean(exp_y);
plot(TAU, ave_exp_y, '-x')
xlabel('Tau, (s)')
ylabel('<sigma_y>')
title('Hahn brute force');

%% dynamic
clc;
close all;
clear all;
hold on;
load('Assignment1_NoiseSamples_MAT', '-mat'); % samples along y  time along x
g_e = 28e9 ; % Hz/T
resolution = 2.5e-6;
fs = 1/(resolution); 
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
w_noise = g_e*noise;
rot_y = -1j*s_y;

S = 100;
t_start = 1e-5;
t_end = 1e-2;
TAU = [4*t_start:40*resolution:t_end];
N = 5; % number of pulses
start = 1/(2*N);
finish = 1 - 1/(2*N);

for n = 1:N
    del_pulse = linspace(start,finish,N);
end
for s = 1:S    % loop tthrough spin
    for t = 1:length(TAU)   %  loop through tau
        pi_flag = 0;
        T = round(TAU(t)/resolution);
        tt = T;  % step through different tau value 
        den_rho = zeros(2,2,T);
        den_rho(:,:,1) = den_psi_t0;
        for i = 2:tt
            w = w_noise(s,i);       
            if(ismember(i, round(del_pulse*floor(TAU(t)/resolution))))
               den_rho(:,:,i-1) =  rot_y*den_rho(:,:,i-1)*rot_y';
            end
            U = expm(-1j*s_z*w*resolution);
            den_rho(:,:,i) = U*den_rho(:,:,i-1)*U';
        end
        exp_y(s,t) = trace(den_rho(:,:,end)*s_y);
    end
    s
end

ave_exp_y = mean(exp_y);
plot(TAU, ave_exp_y, '-x')
xlabel('Tau, (s)')
ylabel('<sigma_y>')
title('CPMG brute force');