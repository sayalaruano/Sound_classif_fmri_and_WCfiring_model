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
load('Stim288.mat')

%% Run the model for the 288 sounds

% Set an empty array to store the fairing rate of all sounds
n_units_allregions = 98*4;
fr_allsounds = zeros(size(stim, 1), n_units_allregions);

for i = 1:size(stim, 1)
    % Select the sound
    s = stim(i, :)';

    % Calling function of the model
    [EE1, EE2, EE3, EE4] = WC_FiringRate_CortAud(Size, lim1_frqaxis, lim2_frqaxis, duration, Fs, mod_rate, mod_depth, carrier_freq, s);
    
    % Get the average firing rate across time for the four regions 
    avg_EE1 = mean(EE1,2)';
    avg_EE2 = mean(EE2,2)';
    avg_EE3 = mean(EE3,2)';
    avg_EE4 = mean(EE4,2)';
    
    % Combine the info in a single vector
    fr_allregions = horzcat(avg_EE1, avg_EE2, avg_EE3, avg_EE4);

    % Append the results of sound i into the array
    fr_allsounds(i, :) = fr_allregions;

end

%% Export dataset
filename = 'fr_allsounds.csv';
writematrix(fr_allsounds, filename);
%% Test for an individual sound 
% Calling function of the model
%[EE1, EE2, EE3, EE4] = WC_FiringRate_CortAud(Size, lim1_frqaxis, lim2_frqaxis, duration, Fs, mod_rate, mod_depth, carrier_freq, s);

% Get the average firing rate across time for the four regions 
%avg_EE1 = mean(EE1,2)';
%avg_EE2 = mean(EE2,2)';
%avg_EE3 = mean(EE3,2)';
%avg_EE4 = mean(EE4,2)';

%fr_allregions = horzcat(avg_EE1, avg_EE2, avg_EE3, avg_EE4);
