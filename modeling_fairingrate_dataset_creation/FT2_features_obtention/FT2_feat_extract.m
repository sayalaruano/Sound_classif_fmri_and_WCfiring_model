function [FT2_feat] = FT2_feat_extract(FR_mat,F_new, timei)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Carry out 2DFT
[k,f,dx,dt] = getAxesFFT2(FR_mat,F_new,timei);
fft2result = fftshift(fft2(FR_mat))*dx*dt;

% Keep spectral modulations larger than 0
indf = find(f>0);

% Keep temporal modulation frequencies below 40 Hz
indk = intersect(find(k>-40),find(k<40));

% Get the output
out = abs(fft2result(indf, indk));

% Calculate mean across spectral and temporal dimensions
FT2_feat = [mean(out), mean(out, 2)'];

end