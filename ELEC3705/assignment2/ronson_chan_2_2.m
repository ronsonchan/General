%% ronson_chan_2_2
close all
clear all
clc

%% definition
sx =  [0,1;1,0];
sy =  [0,1i;-1i,0];
sz =  [1,0;0,-1];
I = eye(2);
sz1 = kron(sz,I);
sz2 = kron(I,sz);
sx1 = kron(sx,I);
sx2 = kron(I,sx);
sy1 = kron(sy,I);
sy2 = kron(I,sy);

Sz1 = kron(sz1,I);
Sz2 = kron(sz2,I);
Sz3 = kron(I,sz2);

Sx1 = kron(sx1,I);
Sx2 = kron(sx2,I);
Sx3 = kron(I,sx2);

Sy1 = kron(sy1,I);
Sy2 = kron(sy2,I);
Sy3 = kron(I,sy2);

up = [1;0];
down = [0;1];

upup = kron(up,up);
updown = kron(up,down);
downup = kron(down,up);
downdown = kron(down,down);

uuu = kron(upup,up);
uud = kron(upup,down);
udu = kron(updown,up);
udd = kron(updown,down);
duu = kron(downup,up);
dud = kron(downup,down);
ddu = kron(downdown,up);
ddd = kron(downdown,down);

ketX = (1/sqrt(2)) * (up+down);
%% a

Sx_sum = Sx1 + Sx2 + Sx3;
[V,D] = eig(Sx_sum);
max_eVal = max(D);
max_eVec = V(:,8);
%% b
ketXXX = kron(ketX, kron(ketX, ketX) );

%% expectation value

psi0 = ketXXX;

exp_XXX1 = psi0' * Sx1 * psi0
exp_XXX2 = psi0' * Sx2 * psi0
exp_XXX3 = psi0' * Sx3 * psi0

% yes they all equal to 1
