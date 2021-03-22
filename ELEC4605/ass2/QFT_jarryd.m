%% QFT
clear all; close all;

hbar = 6.626e-34/(2*pi); % Planck's reduced constant

%% Standard gates
X = [0 1; 1 0]; % Pauli X
Y = [0 -1i; 1i 0]; % Pauli Y
Z = [1 0; 0 -1]; % Pauli Z
I = eye(2); % Identity

H = (1/sqrt(2))*[1 1; 1 -1]; % Hadamard
S = [1 0; 0 1i]; % Phase
T = [1 0; 0 exp(1i*pi/4)]; % Pi/8

H1 = tensor_op(H,3,1); % Hadamard applied to qubit 1 in the combined Hilbert space of 3 qubits
H2 = tensor_op(H,3,2);
H3 = tensor_op(H,3,3);

X1 = tensor_op(X,3,1);
X2 = tensor_op(X,3,2);
X3 = tensor_op(X,3,3);

%% Conditional gates
S21 = control_op(S,3,2,1); % Conditional phase gate applied to qubit 1 with qubit 2 as the control
S32 = control_op(S,3,3,2);
T31 = control_op(T,3,3,1);
SWAP13 = control_op(X,3,1,3)*control_op(X,3,3,1)*control_op(X,3,1,3); % Unitary that swaps the state of qubits 1 and 3

%% Quantum fourier transform for n = 3 qubits
UQFT = SWAP13*H3*S32*H2*T31*S21*H1; % Unitary for the QFT (note the order in which gates are applied)
%% Initial state in the computational basis
psi0 = [1 0 0 0 0 0 0 0];

%% Transform state to the Fourier basis
psi = UQFT*psi0';

%% plot result
subplot(1,2,1);
plot([0:7],psi0,'o','MarkerSize',6,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
ylim([0 1]);
xlabel('Computational basis state |j>');
ylabel('Magnitude');
title('State in the computational basis','FontSize',12)

colororder({'r','b'})

subplot(1,2,2);
yyaxis left
plot([0:7],abs(psi),'o','MarkerSize',6,...
    'MarkerEdgeColor','r','MarkerFaceColor','r');
ylim([0 1]);
xlabel('Fourier basis state |k>');
ylabel('Magnitude');
hold on;
yyaxis right
plot([0:7],angle(psi)/(2*pi),'o','MarkerSize',6,...
    'MarkerEdgeColor','b','MarkerFaceColor','b')
ylim([-0.5 0.5]);
ylabel('Phase (/2\pi)');
hold off;
title('State in the Fourier basis','FontSize',12)