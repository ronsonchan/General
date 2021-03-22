%% FID_filter
close all;
clear all;

w = 2*pi*[0.001:0.001:1];
t = [0.5 1 2 3];

for i = 1:length(t)
    F2(i,:) = t(i)^2*(sin(w*t(i)/2)./(w*t(i)/2)).^2;
end

plot(w,F2)
xlabel('\omega (rad/s)'); ylabel('|F(\omega, t)|^2')
legend('t = 0.5 s', 't = 1 s', 't = 2 s', 't = 3 s')

