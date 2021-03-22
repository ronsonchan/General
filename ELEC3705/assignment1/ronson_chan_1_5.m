%% ronson_chan_1_5
close all;
clc; 

t_finish = 7.31e-7;            % Simulate for 1us
%t_finish = 1e-6;
 psi_0 = 1/sqrt(2) * [1;1]; % Start with |+X> as initial state
% psi_0 = 1/sqrt(2) * [1;1i];  % Start with |+y> as initial state
% psi_0 = [1;0];  % Start with |+Z> as initial state
pulse1.start_time = 0;      % We apply the pulse from t=0
pulse1.stop_time = t_finish;    % We apply the pulse all the way until the end
pulse1.frequency = 20.95e6;  
% pulse1.frequency = 0;  % along z
pulse1.magnitude = 40e-3; 
% pulse1.magnitude = 0;   % along z
pulse1.phase = pi/2;  % along Y
% pulse1.phase = 0; % along X, z

pulse_sequence = [pulse1];  

t_step = 1e-9;
time = 0:t_step:t_finish;
T = length(time);

parfor rep = 1:200
% This call will run the system with our desired initial state
    [t, meas(rep,:)] = quantum_system(t_finish, pulse_sequence);
end
for i = 1:T
    expectationvalueY(1,i) = mean(meas(:,i));
end
%% 

I = eye(2);
sx = [0,1;1,0];
sy = [0,-1i;1i,0];
sz = [1,0;0,-1];
p = I/2 + expectationvalueX(732)*sx + expectationvalueY(732)*sy + expectationvalueZ(1)*sz;
[V,D] = eig(p);
% get rid of the global phase of V(:,2)
newV = V(:,2)*exp(-1i*atan(0.4234/0.2914))
hiddenvariable = newV;


%%
pulse1.frequency = 0;  
pulse1.magnitude = 0; 
pulse1.phase = 0;  