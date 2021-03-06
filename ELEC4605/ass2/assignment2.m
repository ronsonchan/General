%% assignent2
clc; close all; clear all;
%% define constant
n = 3; % number of qubits
hbar = 6.626e-34/(2*pi);
mass = 9.31e-31;   % mass of electron
%w = 2*pi*2.45e9; % res freq of electron
w = 2*pi*1e6;% frequency
%steps = 5e5;
steps = 5e3;
%eps = 1/(w*steps);
eps = 1e-11; % small time step
d = 1e-6; % width of box
beta = (mass*w^2)*eps/(2*hbar);
alpha = eps/(2*hbar*mass);
dx = (2*d)/(2^n);
%% PEdx = 2*d/(2^n);
j = [0:(2^n)-1]; % computational basis state for position
for i = 1:2^n
   x(i) = -d + (j(i) + 1/2)*dx;
end

U_PE = zeros(2^(n-1), 2^(n-1));
U_PE = PE(x,beta);
theta_PE = zeros(1,2^n);
pos_basis = [0:(2^n)-1];
for i = 1:2^n
   theta_PE(i) = angle(U_PE(i,i)); % phase of PE
end
theta_PE
%% KE
dk = pi/d;
m = [0:(2^n)-1]; % basis states in momentum domain
k = (-pi*2^n)/(2*d) + (m).*dk;

U_KE = KE(k,alpha);
theta_KE = zeros(1,2^n);
for i = 1:length(k)
   theta_KE(i) = angle(U_KE(i,i));  % phase of KE
end
theta_KE;
%% time evolution
pauli_x = [0 1;1 0];
psi_s = eye(2^n); % each col stand for a pure state, 
                  % also my computational basis states in position
U_qft = QFT_generalized(n);       % invoke function QFT
Ucs = kron(pauli_x,eye(2^(n-1))); % Cyclic Shift matrix
time = 0:eps:(steps)*eps;         % define time vector 
QQQ = zeros(2^n,length(time));
U_tot = U_qft' * Ucs' * U_KE * Ucs * U_qft * U_PE; 
[V,D] = eig(U_tot); % V stores eigenvector
psi0 = V(:,4);  % initialise state vector of interest
QQQ(:,1) = psi0/norm(psi0); % normalize vector

for t = 1:length(time)-1
    QQQ(:,t+1) = U_qft' * Ucs' * U_KE * Ucs * U_qft * U_PE * QQQ(:,t);
    % plotting in real time
%     plot(j, abs(QQQ(:,t)).^2, 'o')
%     title("Time evolution of Particle in a Box");
%     xlabel("Computational basis in position");
%     ylabel("Expected Probability");
%     ylim([0 1]);
%     xlim([0 (2^n)-1]);
%     colororder({'r'})
%     drawnow;
%     t    
%     pause(0.005)
end
fprintf("Evolution completed!\n");

%
% define puali matrices
X = [0 1; 1 0];
Y = [0 -1j; 1j 0];
Z = [1 0; 0 -1];
final_state = QQQ(:,end);
% measuring in computational basis state  |<k|final>|^2
meas_X = measure(n,final_state);
sum(meas_X) % always add up to 1
range = zeros(1,length(QQQ)); % range stroes norm squared of the wavefunction
% look at time evo of individual qubits
for i = 1:2^n
    for j = 1:length(range)
        range(i,j) = abs(QQQ(i,j)).^2;
    end
%     hold on 
%     plot(time, range(i,:));
%     title("Time Evolution of qubits in Pure State")
%     xlabel("Time (s)")
%     ylabel("Expectation value")
%     legend("qubit 1","qubit 2","qubit 3","qubit 4","qubit 5","qubit 6","qubit 7","qubit 8");
end
hold off