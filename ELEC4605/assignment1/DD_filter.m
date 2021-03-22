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
f = linspace(1e-9,1e3,length(TAU));
w = 2*pi*f;
N = [1:8];
N = 1; % number of pulses
for n = 1:length(N)
    start = 1/(2*N(n));
    finish = 1 - 1/(2*N(n));
    delta = linspace(start,finish,N(n));
    for s = 1:S    % loop tthrough spin
        for t = 1:length(TAU)   %  loop through tau
            pi_flag = 0;
            T = round(TAU(t)/resolution);
            tt = T;  % step through different tau value 
            den_rho = zeros(2,2,T);
            den_rho(:,:,1) = den_psi_t0;
            for i = 2:tt
                w = w_noise(s,i);       
                if(ismember(i, round(delta*floor(TAU(t)/resolution))))
                   den_rho(:,:,i-1) =  rot_y*den_rho(:,:,i-1)*rot_y';
                end
                U = expm(-1j*s_z*w*resolution);
                den_rho(:,:,i) = U*den_rho(:,:,i-1)*U';
            end
            exp_y(s,t) = trace(den_rho(:,:,end)*s_y);
        end
        s
    end
end
ave_exp_y = mean(exp_y);
plot(TAU, ave_exp_y, '-x')
xlabel('Tau, (s)')
ylabel('<sigma_y>')
title('CPMG brute force');

%%
t_start = 1e-5;
t_end = 1e-2;
TAU = [4*t_start:40*resolution:t_end];
f = linspace(1e-9,1e3,length(TAU));
w = 2*pi*f;
N = [1:8];
for n = 1:length(N)
    parfor i = 1:length(TAU)
        start = 1/(2*N(n));
        finish = 1 - 1/(2*N(n));
        delta = linspace(start, finish, N(n));   
        sum_term = 0;   
        for k = 1:N(n)
            sum_term = sum_term + ((-1)^k)*exp(1j*w*delta(k)*TAU(i));
        end
        F2_DD(i,:) = (abs(1 + ((-1)^(1+N(n)))*exp(1j*w*TAU(i)) + 2*sum_term).^2)./((w).^2);
    end
    hold on
    figure(4)
    plot(w,F2_DD(end,:))
    xlabel('\omega (rad/s)'); ylabel('|F(\omega, tau)|^2')
end
hold off
legend;