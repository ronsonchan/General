%% lab2

%% Q1
clc;
clear all
B = [0,0,1]; % magnetic field
s_x = [0,1;1,0];
s_y = [0,-i; i,0];
s_z = [1,0;0,-1];
h_bar = 6.626e-34 / (2*pi);
S = (h_bar/2)*[s_x, s_y, s_z];
gyro = 28e9;  % gyrotio units in Hz/T

Ham_1 = (pi*gyro) * (B(1).*s_x + B(2).*s_y + B(3).*s_z)
Eig_values = eig(Ham_1)

%% Q2
clc;
clear all

Ham_2 =     (2*pi)*[1,0;0,-1] % units in rad/s
[eVector, eValue] = eig(Ham_2);
t = [0:0.1:1];
ket = [1;0];
psi_0 = ket;
c1 = zeros(1,length(t))% coef
c2 = zeros(1,length(t))% coef
psi = zeros(2,length(t));
c1_t0 = ctranspose(eVector(:,1))*ket
c2_t0 = ctranspose(eVector(:,2))*ket
for n = 1:length(t)
   c1(n) = c1_t0 * exp(-i*eValue(1,1)*t(n));
   c2(n) = c2_t0 * exp(-i*eValue(2,2)*t(n));
end
% subplot(2,1,1)
% plot(c1)
% subplot(2,1,2)
% plot(c2)

for n = 1:length(t)
   psi(:,n) = c1(n).*eVector(:,1) + c2(n).*eVector(:,2);
end
teso = zeros(2,length(t));
for n=1:length(t)
   teso(:,n) = (abs(ctranspose(psi(:,n))*ket)).^2
end
%% Q3
clc;
clear all
B = [0,0,1]; % magnetic field
s_x = [0,1;1,0];
s_y = [0,-i; i,0];
s_z = [1,0;0,-1];
h_bar = 6.626e-34 / (2*pi);
S = (h_bar/2)*[s_x, s_y, s_z];
gyro = 28e9;  % gyrotio units in Hz/T

Ham_3 = (pi*gyro) * (B(1).*s_x + B(2).*s_y + B(3).*s_z)
[eValue_3, eVector_3] = eig(Ham_3)
    
up = [1;0];
down = [0;1];
t = [0:0.1:1];
psi_0 = (1/sqrt(2))*(up + down);
c1 = zeros(1,length(t));% coef
c2 = zeros(1,length(t));% coef
psi = zeros(2,length(t));
c1_t0 = ctranspose(eVector_3(:,1))*psi_0;
c2_t0 = ctranspose(eVector_3(:,2))*psi_0;
for n = 1:length(t)
   c1(n) = c1_t0 * exp(-i*eValue_3(1,1)*t(n));
   c2(n) = c2_t0 * exp(-i*eValue_3(2,2)*t(n));
end
for n = 1:length(t)
   psi(:,n) = c1(n).*eVector_3(:,1) + c2(n).*eVector_3(:,2)
end

%% Q4
% expectation value = ctanspose(psi) * S_x * psi
% etc

%% task

clc;
clear all;
up = [1;0];
down = [0;1];
psi_a0 = (1/sqrt(2))*(up + down);
psi_b0 = (1/sqrt(2))*(up + i*down);
psi_c0 = up;
Bz = 57e-6*cos(64*pi/180);
Bx = 57e-6*sin(64*pi/180);
By = 0;
theta = 0;
phi = pi/2 - 64*pi/180; % rad
s_x = [0,1;1,0];
s_y = [0,-i; i,0];
s_z = [1,0;0,-1];
h_bar = 6.626e-34 / (2*pi);
S = (h_bar/2)*[s_x, s_y, s_z]; 
gyro = 28e9;
Ham = (2*pi*gyro*h_bar/2)*(Bx.*s_x + By.*s_y + Bz.*s_z) % units in J
[eVector, eValue] = eig(Ham);

% finding total energy
total_energy_a = eValue(1,1)*ctranspose(psi_a0)*eVector(:,1) + eValue(2,2)*ctranspose(psi_a0)*eVector(:,2)
total_energy_b = abs(eValue(1,1)*ctranspose(psi_b0)*eVector(:,1) + eValue(2,2)*ctranspose(psi_b0)*eVector(:,2))
total_energy_c = eValue(1,1)*ctranspose(psi_c0)*eVector(:,1) + eValue(2,2)*ctranspose(psi_c0)*eVector(:,2)

% finding Lamour frequency
L_a = total_energy_a/h_bar
L_b = total_energy_b/h_bar;
L_c = total_energy_c/h_bar;

% idk
t = linspace(0,2*pi/L_a,50);
psi_at = zeros(2, length(t));
psi_bt = zeros(2, length(t));
psi_ct = zeros(2, length(t));
for n = 1:length(t)
    psi_at(:,n) = expm(-1i * Ham * t(n) / h_bar) * psi_a0;
    psi_bt(:,n) = expm(-1i * Ham * t(n) / h_bar) * psi_b0;
    psi_ct(:,n) = expm(-1i * Ham * t(n) / h_bar) * psi_c0;
end
exp_ax = zeros(2,length(t));
exp_ay = zeros(2,length(t));
exp_az = zeros(2,length(t));
exp_bx = zeros(2,length(t));
exp_by = zeros(2,length(t));
exp_bz = zeros(2,length(t));
exp_cx = zeros(2,length(t));
exp_cy = zeros(2,length(t));
exp_cz = zeros(2,length(t));
for n = 1: length(t)
   exp_ax(:,n) = abs((psi_at(:,n))' * (h_bar/2)*s_x * (psi_at(:,n)));
   exp_ay(:,n) = abs((psi_at(:,n))' * (h_bar/2)*s_y * (psi_at(:,n)));
   exp_az(:,n) = abs((psi_at(:,n))' * (h_bar/2)*s_z * (psi_at(:,n)));
   exp_bx(:,n) = abs((psi_at(:,n))' * (h_bar/2)*s_x * (psi_at(:,n)));
   exp_by(:,n) = abs((psi_at(:,n))' * (h_bar/2)*s_y * (psi_at(:,n)));
   exp_bz(:,n) = abs((psi_at(:,n))' * (h_bar/2)*s_z * (psi_at(:,n)));
   exp_cx(:,n) = abs((psi_at(:,n))' * (h_bar/2)*s_x * (psi_at(:,n)));
   exp_cy(:,n) = abs((psi_at(:,n))' * (h_bar/2)*s_y * (psi_at(:,n)));
   exp_cz(:,n) = abs((psi_at(:,n))' * (h_bar/2)*s_z * (psi_at(:,n)));
end
figure(1)
hold on
plot(t, exp_ax)
plot(t, exp_ay)
plot(t, exp_az)

figure(2)
hold on
plot(t, exp_bx)
plot(t, exp_by)
plot(t, exp_bz)

figure(3)
hold on
plot(t, exp_cx)
plot(t, exp_cy)
plot(t, exp_cz)
hold off
