function [energy] = energy(n)
% Calculates the energy in the nth state
a = 1e-6;  % width of the well
m = 9.11*10^-31;
h_bar = 1.05*10^-34;

kn = n*pi/a;



energy = (h_bar^2*kn^2)/(2*m);

end

