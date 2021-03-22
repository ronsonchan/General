%% lab 4

%% Q1
gyro = 2;   % electron gyration
sz = [1,0;0,-1];
t = [0:100];
w = 1 % angular frequency
Bx = 1;
By = Bx;
Bz = 1;
B = [Bx, By, Bz];
delta = Bx*exp(i*w*t);
eps = Bz;
H = pi*gyro*[Bz, Bx*cos(w*t) - i*By*sin(w*t);...
             Bx*cos(w*t) - i*By*sin(w*t), -Bz];

% eigenvalues = lamda = +-sqrt(Bz^2 + Bx^2)
% eigenvectors in terms of sigma_z: 
% v1 = sin(theta/2)*up + cos(theta/2)*down
% v2 = cos(tehta/2)*down - sin(theta/2)*up 
% for tan(theta/2 = delta / eps

%% Q2 
% original state  of electron = v1 or  v2 depending on the spin
% now since we're rotating the lab from instead, we can transform our basis vector
% from up and down to v1' = [exp(-iwt/2), 0; 0, exp(iwt/2)]*v1 or
% v2' = [exp(-iwt/2), 0; 0, exp(iwt/2)]*v2
% with a new hamiltonian:
% H' = h_bar/2 * [-delta_w, w1; w1, delta_w]; where w1 = gyro*B_perpen. w0 = gyro*B_parallel, 
% H' = -h_bar/2 * delta_w * sigma_z + h_bar/2 * w1 * sigma_x 
% delta_w = w - w0.
% 
% if w = w0 = gyro*B_parallel, then delta_w = 0
% then H' solely depends on w1 ~ B_perpendicular, ie the electron spin is rotating at the same speed as the lab frame, 
% ie, no relative motion seen from the rotating frame.
% from the lab frame, the electron spin will only precess around the new X basis vector

%% Q3
% clc
% h_bar = 6.626e-34 / (2*pi);
% H
% psi_0 = up
% psi = U*psi_0, where U = exp(-1i*En*t/h_bar) in the diagonal, zero else
% where, E = eignevalue of H

%% task  a

clc;
clear all
close all


h_bar = 6.626e-34/(2*pi);
g = 28e9;
B_para = 1; % Bz
B_perp = 20e-3; % Bx
w0 = 2*pi*g*B_para;
% del_w = 2*pi*10e6;  % for part a
del_w = 0;
w1 = 2*pi*g*B_perp; 
w = w0-del_w;
eps = h_bar/2 * del_w;
delta = (h_bar/2) * w1;
wrabi_lab = (sqrt(delta^2 + eps^2))/(h_bar);
t = linspace(0,2*pi/wrabi_lab,10000);
T = length(t);
H_lab = zeros(2,2,length(t));
for i = 1:T
       H_lab(:,:,i) = (pi*g*h_bar)*[B_para, B_perp*exp(-1i*(w+del_w)*t(i)); B_perp*exp(1i*(w+del_w)*t(i)), -B_para] ;
end
H_rot = (h_bar/2) * [-del_w, w1; w1, del_w];
wrabi_rot = (sqrt((del_w*h_bar/2)^2 + (w1*h_bar/2)^2))/(h_bar);

% now do ket from lab 2

down = [0;1];
up = [1;0];
ket_0 = up;
for i = 1:T
    phi1 = [exp(-1i*w*t(i)/2), 0; 0, exp(1i*w*t(i)/2)] * up;
    phi2 = [exp(-1i*w*t(i)/2), 0; 0, exp(1i*w*t(i)/2)] * down;
end
psi_lab = zeros(2,T);
psi_rot = zeros(2,T);
psi_lab(:,1) = up;
% choose up as ket0 or down for labframe, your choice
time_interval = t(2)-t(1);
for i = 2:T
    U = expm(-1i*H_lab(:,:,i)*time_interval/h_bar);
    psi_lab(:,i) = U*psi_lab(:,i-1);
    
end
% rotating frame
psi_rot(:,1) = phi1;
for  i = 2:T
    U = expm(-1i*H_rot*t(i)/h_bar);
    psi_rot(:,i) = U*psi_rot(:,1);
end

% expectational value
s_x = [0,1;1,0];
s_y = [0,-1i;1i,0];
s_z = [1,0;0,-1];
ex_x1 = zeros(1,T);
ex_y1 = zeros(1,T);
ex_z1 = zeros(1,T);
ex_x2 = zeros(1,T);
ex_y2 = zeros(1,T);
ex_z2 = zeros(1,T);
for count = 1:T
   ex_x1(:, count) = psi_lab(:,count)' * s_x * psi_lab(:, count); 
   ex_y1(:, count) = psi_lab(:,count)' * s_y * psi_lab(:, count); 
   ex_z1(:, count) = psi_lab(:,count)' * s_z * psi_lab(:, count); 
   ex_x2(:, count) = psi_rot(:,count)' * s_x * psi_rot(:, count); 
   ex_y2(:, count) = psi_rot(:,count)' * s_y * psi_rot(:, count);
   ex_z2(:, count) = psi_rot(:,count)' * s_z * psi_rot(:, count);
end

hold on
plot(t, ex_x2);
plot(t, ex_y2);
plot(t, ex_z2);
hold off