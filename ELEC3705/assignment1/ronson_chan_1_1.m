%% ronson_chan_1_1

close all; 
clear all;

for i = 1:100
% We will define all the parameters of our experiment here
t_finish = 1e-6;            % Simulate for 1us
% psi_0 = 1/sqrt(2) * [1;1];  % Start with |+x> as initial state
psi_0 = [1;0];  % Start with |+z> as initial state
% We next have to define the parameters of our pulse. These are stored
% as a struct in MATLAB in the following way.
pulse1.start_time = 0;      % We apply the pulse from t=0
pulse1.stop_time = 1e-6;    % We apply the pulse all the way until the end
pulse1.frequency = 18e6;    % Pulse frequency is 18MHz
pulse1.magnitude = 20e-3;   % The field strength is 20mT
pulse1.phase = pi/2;        % We have a phase of pi/2 radians

% List all pulses you want to apply in this array
pulse_sequence = [pulse1];  

% This call will run the system with our desired initial state
% [t, meas] = quantum_system(t_finish, pulse_sequence, psi_0);
% This runs the system with the hidden initial state
 [t, meas] = quantum_system(t_finish, pulse_sequence);


meas_over_time(i,:) = meas;
%ave(i) = sum(meas, 'all')/1001;
ave(i) = mean(meas);
end
% figure(1)
% plot(t, meas_over_time);
% xlabel('time(ns)')
% ylabel('expectation value')

figure(2)
plot(t(1:100), ave)