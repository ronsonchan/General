%% ronson_chan_2_3

clc
clear all
close all

% parameters definition
% B0 = linspace(0,0.1,32); % T change B strength from 0 to 0.1 with 32 increments
B0 = 0.1;
gnp = -17.235e6 * 2*pi; % in radians
gnsi = 8.458e6 * 2*pi;
ge = 2*pi*28e9;
Ap = 117e6 *2*pi;
sx = [0,1;1,0];
sy = [0,-1i;1i,0];
sz = [1,0;0,-1];
hbar = 1;
I = eye(2);

N = 3;
% prompt = "Enter value for N: ";
% N = input(prompt);  % N is the number of silicon atom 
% for q2.1, choose N = 3
size = N+2; % it's always N silicon atom plus 1 electron 
            % plus 1 phosphours 

spin_x = zeros(2^size, 2^size, size);  
spin_y = zeros(2^size, 2^size, size); 
spin_z = zeros(2^size, 2^size, size); 
% initialize a 3D array where X/Y plane stores 
% the data of your tensor
% and have the order of spins stacking
% on top of each other
           
% for X orientation
for i = 1:size
    if i == 1
        spin_x(1:2,1:2,i) = sx;      
    else
        spin_x(1:2,1:2,i) = I;
   
    end
    for j = 2:size
        if i == j
            spin_x(1:2^j,1:2^j,i) = kron(spin_x(1:2^(j-1),1:2^(j-1),i),sx);
        else
            spin_x(1:2^j,1:2^j,i) = kron(spin_x(1:2^(j-1),1:2^(j-1),i),I);
        end
    end
end
% for y orientation
for i = 1:size
    if i == 1
        spin_y(1:2,1:2,i) = sy;
        
    else
        spin_y(1:2,1:2,i) = I;
   
    end
    for j = 2:size
        if i == j
            spin_y(1:2^j,1:2^j,i) = kron(spin_y(1:2^(j-1),1:2^(j-1),i),sy);
        else
            spin_y(1:2^j,1:2^j,i) = kron(spin_y(1:2^(j-1),1:2^(j-1),i),I);
        end
    end    
end
% for z orientation
for i = 1:size
    if i == 1
        spin_z(1:2,1:2,i) = sz;
        
    else
        spin_z(1:2,1:2,i) = I;
   
    end
    for j = 2:size
        if i == j
            spin_z(1:2^j,1:2^j,i) = kron(spin_z(1:2^(j-1),1:2^(j-1),i),sz);
        else
            spin_z(1:2^j,1:2^j,i) = kron(spin_z(1:2^(j-1),1:2^(j-1),i),I);
        end
    end   
end
% q2.2 materials

% set up initial state here
psi_x = (1/sqrt(2))*[1;1];
psi_z =  [0;1];
psi_0 = kron(psi_x, psi_z);
for i = 1:N
    v = [1; exp(1i*2*pi*rand())].*rand(2,1);
    v = v/norm(v);
    psi_0 = kron(psi_0,v);
end
% time evolution
t_start = 0;
t_ends = 2e-8;
% t_ends = 1e-9
t_step = 1e-11;
time = t_start:t_step:t_ends;
T = length(time);

psi = zeros(2^size, T);
psi(:,1) = psi_0;
A = zeros(N,1);     % malloc for A = rand()
expec_x = zeros(1,T);
rep = 50;
for n = 1
    Hi = zeros(2^size, 2^size) ;
    Hsi = zeros(2^size, 2^size);
    for i = 1:N
        A(i) = rand() * 2*pi*50e6/hbar; 
        Hsi = Hsi + gnsi*B0*spin_z(:,:,i+2);
        Hi = Hi + A(i)*(spin_x(:,:,1)*spin_x(:,:,i+2)...
            + spin_y(:,:,1)*spin_y(:,:,i+2) + spin_z(:,:,1)*spin_z(:,:,i+2));
    end
    H_tot = -ge*B0*spin_z(:,:,1) + gnp*B0*spin_z(:,:,2)...
            + Ap*(spin_x(:,:,1)*spin_x(:,:,2) + spin_y(:,:,1)*spin_y(:,:,2) ...
            + spin_z(:,:,1)*spin_z(:,:,2)) + Hsi + Hi;
    for i = 1:T
        % time evolution operator U
        U = expm(-1i*H_tot*time(i)/hbar);
        psi(:,i) = U*psi_0;
        expec_x(n,i) = psi(:,i)' * spin_x(:,:,1) * psi(:,i);
    end
end
ave_exp = zeros(1,T);
% construct expectation value
for i = 1:T
   ave_exp(i) = mean(expec_x(:,i)); 
end

plot(time, ave_exp, 'LineWidth', 2)
xlabel('time');
ylabel('expectation value');
set(gca, 'Fontsize', 18);

%  comment
% 