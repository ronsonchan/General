clc;
clear all;
h_bar = 6.626e-34;

Sx = (h_bar/2)*[0,1;1,0];
Sy = (h_bar/2)*[0,-i;i,0];
Sz = (h_bar/2)*[1,0;0,-1];
up = [1;0];
down = [0;1];
sum_x = (up + down);
diff_x = (up - down);
sum_y = up + i*down;
diff_y = up - i*down;

X_plus = Sx*(1/sqrt(2))*sum_x
X_minus = Sx*(1/sqrt(2))*diff_x

Y_plus = Sy*(1/sqrt(2))*sum_y
Y_minus = Sy*(1/sqrt(2))*diff_y


%% time evolution
clc
psi_0 = up;
H = Sz;
t = [0;1];
psi_t = expm(-1i * H .* t/h_bar) * psi_0;
plot(psi_t)