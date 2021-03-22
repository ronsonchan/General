%% ronson_chan_2_1
close all;
clear all
clc;

%% definition
sx =  [0,1;1,0];
sy =  [0,1i;-1i,0];
sz =  [1,0;0,-1];
I = eye(2);
up = [1;0];
down = [0;1];

upup = kron(up,up);
updown = kron(up,down);
downup = kron(down,up);
downdown = kron(down,down);

sz1 = kron(sz,I);
sz2 = kron(I,sz);
sx1 = kron(sx,I);
sx2 = kron(I,sx);
sy1 = kron(sy,I);
sy2 = kron(I,sy);
%% contruct 3-coupled spin
clc
% uuu = upupup
% udu = updownup
uuu = kron(upup,up);
uud = kron(upup,down);
udu = kron(updown,up);
udd = kron(updown,down);
duu = kron(downup,up);
dud = kron(downup,down);
ddu = kron(downdown,up);
ddd = kron(downdown,down);

%% construct 3 coupled operator
clc
% capital S donates the operator for 3 coupled state

Sz1 = kron(sz1,I);
Sz2 = kron(sz2,I);
Sz3 = kron(I,sz2);

Sx1 = kron(sx1,I);
Sx2 = kron(sx2,I);
Sx3 = kron(I,sx2);

Sy1 = kron(sy1,I);
Sy2 = kron(sy2,I);
Sy3 = kron(I,sy2);




I = eye(2);
sz =  [1,0;0,-1];
sz1 = kron(sz,I);
Sz1 = kron(sz1,I);