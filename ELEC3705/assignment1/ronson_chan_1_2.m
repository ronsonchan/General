%% ronson_chan_1_2

close all; 
clear all;

% We will define all the parameters of our experiment here
t_finish = 10e-6;            % Simulate for 1us
% psi_0 = 1/sqrt(2) * [1;1];  % Start with |+x> as initial state
psi_0 = [1;0];  % Start with |+z> as initial state
% We next have to define the parameters of our pulse. These are stored
% as a struct in MATLAB in the following way.
pulse1.start_time = 0;      % We apply the pulse from t=0
pulse1.stop_time = 1e-5;    % We apply the pulse all the way until the end

pulse1.magnitude = 40e-3;   % The field strength is 20mT
pulse1.phase = pi/2;        % We have a phase of pi/2 radians
t_step = 1e-9;

% note to self
% I have picked 3 sections of frequency range of size 20 
% which were: 18M to 20M, 20M to 22M and 22M to 24M
% and the expectation value drops to the lowest in the section
% 20M to 22M, therefore I believe the transition frequency 
% lies within this range
% the idea is similar to newton's method where you keep
% halving your section until you get to the closest one.
% I had to do it manuelly since my coding skills are not
% as proficient as others.

pulse1.frequency = 20.95e6;    % Pulse frequency is 18MHz initially
pulse_sequence = [pulse1]; 
 
time = [0:t_step:pulse1.stop_time];
T = length(time);

   
parfor rep = 1:100
     [t, meas(rep,:)] = quantum_system(t_finish, pulse_sequence, psi_0);  
end

for i = 1:T
    expectationvalue(1,i) = mean(meas(:,i));
end

figure(1)
plot(time, expectationvalue, 'LineWidth', 2);  % found 24.56MHz has the lowest
xlabel('time(ns)')
ylabel('average expectation value')
set(gca, 'Fontsize', 18);

% time required for 180 deg rotation: 1.46us
