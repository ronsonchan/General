%% Lab1

%% Q1
clc
close all

%h_bar = 6.626e-34 / (2*pi);
Sx = [0,1;1,0];
Sy = [0,-i;i,0];
Sz = [1,0;0,-1];
up = [1;0];
down = [0;1];
psi = (1/sqrt(2))*(up + down);

expectation_value_x = ctranspose(psi) * Sx * psi
expectation_value_y = ctranspose(psi) * Sy * psi
expectation_value_z = ctranspose(psi) * Sz * psi
%% Q2

% state vector for x = (1/sqrt(2))*[1:1] or (1/sqrt(2))*[1:-1]

%% Q3 

% also 50%

%% task
%% a
clear all
clc

up = [1;0];
down = [0;1];
Sx = [0,1;1,0];
Sy = [0,-i;i,0];
Sz = [1,0;0,-1];

psi1 = [1;0];
count = 3;
V = zeros(3:3);
S = zeros(1:3);
psi = psi1;
for i = 1:count
    if i ~= 3  
        [V,D] = eig(Sx);
        expectation_value_x = ctranspose(psi) * Sx * psi
        probability_minus1 = abs(ctranspose(V(:,1)) * psi)^2
        probability_plus1 = abs(ctranspose(V(:,2)) * psi)^2
        psi = (1/sqrt(2))*(D(:,1) - D(:,2));
    else
        [V,D] = eig(Sz);
        expectation_value_z = ctranspose(psi) * Sz * psi
        probability_minus1 = abs(ctranspose(V(:,1)) * psi)^2
        probability_plus1 = abs(ctranspose(V(:,2)) * psi)^2
        psi = (1/sqrt(2))*(D(:,1) - D(:,2));
    end
end
%% b
clear all
clc

up = [1;0];
down = [0;1];
Sx = [0,1;1,0];
Sy = [0,-i;i,0];
Sz = [1,0;0,-1];

psi1 = (1/sqrt(2))*up + ((1+i)/2)*down;
count = 3;
V = zeros(3:3);
S = zeros(1:3);
psi = psi1;
for i = 1:count
    if i ~= 3  
        [V,D] = eig(Sx);
        expectation_value_x = ctranspose(psi) * Sx * psi
        probability_minus1 = abs(ctranspose(V(:,1)) * psi)^2
        probability_plus1 = abs(ctranspose(V(:,2)) * psi)^2
        psi = (1/sqrt(2))*(D(:,1) - D(:,2));
    else
        [V,D] = eig(Sz);
        expectation_value_z = ctranspose(psi) * Sz * psi
        probability_minus1 = abs(ctranspose(V(:,1)) * psi)^2
        probability_plus1 = abs(ctranspose(V(:,2)) * psi)^2
        psi = (1/sqrt(2))*(D(:,1) - D(:,2));
    end
end

%% c
clc
clear all

down = [0;1];
up = [1;0];
Sx = [0,1;1,0];
Sy = [0,-i;i,0];
Sz = [1,0;0,-1];
% first B field (S.G.)
[V,D] = eig(Sx);
psi = down;
expectation_value_x = ctranspose(psi) * Sx * psi
probability_minus1 = abs(ctranspose(V(:,1)) * psi)^2
probability_plus1 = abs(ctranspose(V(:,2)) * psi)^2
psi = (1/sqrt(2))*(D(:,1) + D(:,2))

% second B field
theta = 35*pi/180 % radian
phi = 20*pi/180 % radian
Su = [cos(theta),sin(theta)*exp(-i*phi); sin(theta)*exp(i*phi), -cos(theta)]  
% in radian
[V,D] = eig(Su);
expectation_value_u = abs(ctranspose(psi) * Su * psi)
probability_minus1 = abs(ctranspose(V(:,1)) * psi)^2
probability_plus1 = abs(ctranspose(V(:,2)) * psi)^2
psi = V(:,1) 

% third B field
% [V,D] = eig(Sz);
% expectation_value_z = ctranspose(psi) * Sz * psi
% probability_minus1 = abs(ctranspose(V(:,1)) * psi)^2
% probability_plus1 = abs(ctranspose(V(:,2)) * psi)^2
% psi = probability_minus1 * D(:,1) + probability_plus1 * D(:,2)
