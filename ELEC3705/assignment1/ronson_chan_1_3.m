%% ronson_chan_1_3

clear all;

t_finish = 50e-6;            

psi_0 = [1;0];  % Start with |+z> as initial state

pulse1.start_time = 0;      % We apply the pulse from t=0
pulse1.stop_time = 50e-6;    % We apply the pulse all the way until the end
pulse1.frequency = 20.95e6;    % Pulse frequency at transition frequency
pulse1.magnitude = 40e-3;   % The field strength is 20mT
pulse1.phase = pi/2;        % We have a phase of pi/2 radians
t_step = 1e-9;
% List all pulses you want to apply in this array
pulse_sequence = [pulse1];  

time = [0:t_step:pulse1.stop_time];
T = length(time);
% This call will run the system with our desired initial state
% parfor i = 1:100
%     [t, meas(i,:)] = quantum_system(t_finish, pulse_sequence, psi_0);
% end
% for i = 1:T
%     expectationvalue(1,i) = mean(meas(:,i));
% end
% plot(time, expectationvalue);
% xlabel("time(ns)");
% ylabel("expectationvalue");

