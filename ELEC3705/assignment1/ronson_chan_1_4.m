%% ronson_chan_1_4
close all;
clear all; 
clc

t_finish = 7.31e-7;            % Simulate for 1us
% t_finish = 1e-6;
psi_0 = 1/sqrt(2) * [1;1]; % Start with |+X> as initial state
% psi_0 = 1/sqrt(2) * [1;1i];  % Start with |+y> as initial state

pulse1.start_time = 0;      % We apply the pulse from t=0
pulse1.stop_time = t_finish;    % We apply the pulse all the way until the end
pulse1.frequency = 20.95e6;    
pulse1.magnitude = 40e-3;   
% pulse1.phase = pi/2;  % along Y
pulse1.phase = 0; % along X

pulse_sequence = [pulse1];  

t_step = 1e-9;
time = 0:t_step:t_finish;
T = length(time);
parfor rep = 1:200
% This call will run the system with our desired initial state
    [t, meas(rep,:)] = quantum_system(t_finish, pulse_sequence, psi_0);
end

for i = 1:T
    expectationvalue(1,i) = mean(meas(:,i));
end
plot(time, expectationvalue, 'LineWidth', 2);
xlabel("time(ns)");
ylabel("expectational value");
title('expectational value along X vs time with phase = 0');
set(gca, 'Fontsize', 18);

% along the X axis
% with phase = 0, 
% the time it takes to rotate 90: t = 7.31e-7s = duration
% at frequency = 20.95MHz

% along Y 
% phase  = pi/2
% duration = 7.31e-7 s
% at frequency = 20.95MHz