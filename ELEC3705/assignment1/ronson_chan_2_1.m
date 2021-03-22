%% ronson_chan_2_1
clc
clear all
close all

% parameters definition
B0 = 0.1; %T
gnp = -17.235e6 * 2*pi; 
gnsi = 8.458e6 * 2*pi;
Ap = 117e6 *2*pi;
sx = [0,1;1,0];
sy = [0,-1i;1i,0];
sz = [1,0;0,-1];
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
