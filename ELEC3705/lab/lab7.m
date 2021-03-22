%% lab7

% definition 
clc;
clear all;
close all;

up = [1;0];
down = [0;1];
upup = [1;0;0;0];
updown = [0;1;0;0];
downup = [0;0;1;0];
downdown = [0;0;0;1];
hbar = 6.626e-34 / (2*pi);
sx = 1/2 * [0,1;1,0];
sy = 1/2 * [0,1i;-1i,0];
sz = 1/2 * [1,0;0,-1];
I = eye(2);
sz1 = kron(sz, I);
sz2 = kron(I,sz);
sztot = sz1+sz2;

sx1 = kron(sx, I);
sx2 = kron(I,sx);
sxtot =  sx1+sx2;

sy1 = kron(sy,I);
sy2 = kron(I,sy);
sytot = sy1+sy2;

stot = sxtot^2 + sytot^2 + sztot^2;

%% Q1a

psi = (1/sqrt(2)) * ( upup + updown);

%% Q1b

psi = (1/sqrt(2)) * (updown + downup);
%% Q1

exp_z = psi' * sztot * psi
exp_y = psi' * sytot * psi
exp_x = psi' * sxtot * psi

exp_tot = psi' * stot * psi


%% Q3a

psi = (1/sqrt(2)) * (upup+updown);

%% Q3b

psi = (1/sqrt(2))*(updown + downup);

%% Q3c

psi = (1/sqrt(3)) * (updown + downup + downdown);
%%
u = [1,0;0,1];
rho_C = psi * psi';
rho_partialA = zeros(2);

for i = 1:2
    for j = 1:2
        for n = 1:2
            rho_partialA(i,j) = rho_partialA(i,j) + kron(u(:,i)',u(:,n)') * rho_C * kron(u(:,j),u(:,n));
        end
    end
end

rho_partialB = zeros(2);

for i = 1:2
    for j = 1:2
        for n = 1:2
            rho_partialB(i,j) = rho_partialB(i,j) + kron(u(:,n)',u(:,i)') * rho_C * kron(u(:,n),u(:,j));
        end
    end
end

S = -trace(rho_partialB * logm(rho_partialB))


%% task

%% definition
clear all;
clc
close all;

up = [1;0];
down = [0;1];
upup = [1;0;0;0];
updown = [0;1;0;0];
downup = [0;0;1;0];
downdown = [0;0;0;1];
%hbar = 6.626e-34 / (2*pi);
hbar = 1;
sx = 1/2 * [0,1;1,0];
sy = 1/2 * [0,-1i;1i,0];
sz = 1/2 * [1,0;0,-1];
I = eye(2);
sz1 = kron(sz, I);
sz2 = kron(I,sz);
sztot = sz1+sz2;

sx1 = kron(sx, I);
sx2 = kron(I,sx);
sxtot =  sx1+sx2;

sy1 = kron(sy,I);
sy2 = kron(I,sy);
sytot = sy1+sy2;

stot = sxtot^2 + sytot^2 + sztot^2;

wh = 2*pi*10e6; % rads-1
H = (wh/2)*[1,1,1,1;1,-1,1,-1;1,1,-1,-1;1,-1,-1,1];
period = (pi)/wh;
time = linspace(0,period, 1000);
T = length(time);
delta = 1e-12;

%% task a
clc
psi_0 = downdown;

%% task b

clc;
psi_0 = (1/sqrt(2))* (downdown + upup);
%% time evolation/expectation values
clc
u = [1 0; 0 1]; % basis vector
U = zeros(4,4);
psi = zeros(4,T);
exp_x = zeros(1,T);
exp_y = zeros(1,T);
exp_z = zeros(1,T);
exp_sum = zeros(1,T);
exp_tot = zeros(1,T);
S = zeros(1,T);
rho_partialA = zeros(2);
rho_partialB = zeros(2);
rho_C = zeros(4,4);
for t = 1:T
   U = expm(-1i*H*time(t));
   psi(:,t) = U*psi_0;
   exp_x(:,t) = real(psi(:,t)' * sxtot * psi(:,t));
   exp_y(:,t) = real(psi(:,t)' * sytot * psi(:,t));
   exp_z(:,t) = real(psi(:,t)' * sztot * psi(:,t));
   exp_sum(:,t) = real(exp_x(:,t)+exp_y(:,t)+exp_z(:,t));
   exp_tot(:,t) = real(psi(:,t)' * stot * psi(:,t)); 
   
   rho_C = psi(:,t) * psi(:,t)';
   rho_partialB = zeros(2);
    for i = 1:2
        for j = 1:2
            for n = 1:2
                rho_partialB(i,j) = rho_partialB(i,j) + kron(u(:,i)',u(:,n)') * rho_C * kron(u(:,j),u(:,n));
            end
        end
    end
    rho_partialA = zeros(2);
    for i = 1:2
        for j = 1:2
            for n = 1:2
                rho_partialA(i,j) = rho_partialA(i,j) + kron(u(:,n)',u(:,i)') * rho_C * kron(u(:,n),u(:,j));
            end
        end
    end
   S(t) = -trace(rho_partialB * logm(rho_partialB));
end

%% plot result
subplot(2,1,1)
hold on
plot(time,exp_x);
plot(time,exp_y);
plot(time,exp_z);
plot(time,exp_sum);
plot(time,exp_tot);
xlim([0 50e-9]) 
hold off
subplot(2,1,2)
plot(time, S)



