%% Set frequency variables
Size = 100;  % Spatial (Frequency) size of array - number of units
lim1_frqaxis = 50;
lim2_frqaxis = 8000;
duration = 1; % second
Fs = 16000;
mod_rate = 6; % Hz
mod_depth = 1; % from 0 to 1
carrier_freq = 1000; % Hz

%% Load the data
load('../Stim288.mat')

%% Test for an individual sound 
% Get sound
s = stim(1, :)';

% Calling function of the model
[EE1, EE2, EE3, EE4, F_new, timei] = WC_FiringRate_CortAud_2FT(Size, lim1_frqaxis, lim2_frqaxis, duration, Fs, mod_rate, mod_depth, carrier_freq, s);

% Get the average firing rate across time for the four regions 
avg_EE1 = mean(EE1,2)';
avg_EE2 = mean(EE2,2)';
avg_EE3 = mean(EE3,2)';
avg_EE4 = mean(EE4,2)';

% Carry out 2DFT
FT2_featA1 = FT2_feat_extract(EE1, F_new, timei);
FT2_featR = FT2_feat_extract(EE2, F_new, timei);
FT2_featS = FT2_feat_extract(EE3, F_new, timei);
FT2_featF = FT2_feat_extract(EE4, F_new, timei);

fr_allregions = horzcat(avg_EE1, FT2_featA1, avg_EE2, FT2_featR, avg_EE3, FT2_featS, avg_EE4, FT2_featF);

%% Run the model for the 288 sounds

% Set an empty array to store the fairing rate of all sounds
n_units_allregions = 98*4 + 128*4;
fr_allsounds = zeros(size(stim, 1), n_units_allregions);
n_units_perregion = 98 + 128;
fr_allsoundsA1 = zeros(size(stim, 1), n_units_perregion);
fr_allsoundsR = zeros(size(stim, 1), n_units_perregion);
fr_allsoundsS = zeros(size(stim, 1), n_units_perregion);
fr_allsoundsF = zeros(size(stim, 1), n_units_perregion);


for i = 1:size(stim, 1)
    % Select the sound
    s = stim(i, :)';

    % Calling function of the model
    [EE1, EE2, EE3, EE4, F_new, timei] = WC_FiringRate_CortAud_2FT(Size, lim1_frqaxis, lim2_frqaxis, duration, Fs, mod_rate, mod_depth, carrier_freq, s);
    
    % Get the average firing rate across time for the four regions 
    avg_EE1 = mean(EE1,2)';
    avg_EE2 = mean(EE2,2)';
    avg_EE3 = mean(EE3,2)';
    avg_EE4 = mean(EE4,2)';

    % Carry out 2DFT
    FT2_featA1 = FT2_feat_extract(EE1, F_new, timei);
    FT2_featR = FT2_feat_extract(EE2, F_new, timei);
    FT2_featS = FT2_feat_extract(EE3, F_new, timei);
    FT2_featF = FT2_feat_extract(EE4, F_new, timei);
    
    % Combine the info in a single vector
    fr_allregions = horzcat(avg_EE1, FT2_featA1, avg_EE2, FT2_featR, avg_EE3, FT2_featS, avg_EE4, FT2_featF);

    % Append the results of sound i into the array
    fr_allsounds(i, :) = fr_allregions;
    fr_allsoundsA1(i, :) = horzcat(avg_EE1, FT2_featA1);
    fr_allsoundsR(i, :) = horzcat(avg_EE2, FT2_featR);
    fr_allsoundsS(i, :) = horzcat(avg_EE3, FT2_featS);
    fr_allsoundsF(i, :) = horzcat(avg_EE4, FT2_featF);

end

%% Export datasets
filename_allregions = 'fr_allsounds_2FT.csv';
writematrix(fr_allsounds, filename_allregions);

filenameA1 = 'fr_allsoundsA1_2FT.csv';
writematrix(fr_allsoundsA1, filenameA1);

filenameR = 'fr_allsoundsR_2FT.csv';
writematrix(fr_allsoundsR, filenameR);

filenameS = 'fr_allsoundsS_2FT.csv';
writematrix(fr_allsoundsS, filenameS);

filenameF = 'fr_allsoundsF_2FT.csv';
writematrix(fr_allsoundsF, filenameF);

