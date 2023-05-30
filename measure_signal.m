function [y_abs,y_ph,A] = measure_signal(m,x)
n = length(x);
%% signal measurement
A = randn(m,n);
y = A*x; %measurements
y_abs = abs(y);
y_ph = sign(y); %actual phase
end