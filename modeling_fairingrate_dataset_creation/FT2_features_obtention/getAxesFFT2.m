
function [k,f,dx,dt] = getAxesFFT2(EE1,F_new,t)

Nx = size(EE1,2); % Number of samples collected along first dimension
Nt = size(EE1,1); % Number of samples collected along second dimension
dx = t(2)- t(1);  % Distance increment (i.e., Spacing between each column)
dt = hz2erb(F_new(2))-hz2erb(F_new(1)); % Time increment (i.e., Spacing between each row)
x = 0 : dx : (Nx-1)*dx; % distance
t = 0 : dt : (Nt-1)*dt; % time
data_spacetime = randn(Nt,Nx); % random 2D matrix
Nyq_k = 1/(2*dx); % Nyquist of data in first dimension
Nyq_f = 1/(2*dt); % Nyquist of data in second dimension
dk = 1/(Nx*dx);   % Wavenumber increment
df = 1/(Nt*dt);   % Frequency increment
k = -Nyq_k : dk : Nyq_k-dk; % wavenumber
f = -Nyq_f : df : Nyq_f-df; % frequency