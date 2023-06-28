%% Spectro-temporal processing in a two-stream computational model of the Auditory Cortex

% This script embedded functions from :
%       Wilson Cowan cortical model 
%       Wilson HR. Computation by excitatory and inhibitory networks. In: Spikes, Decisions & Actions: Dynamical Foundations of Neuroscience. Oxford UK: Oxford University Press; 1999. p. 88–115.
%
%       4th order gammatone filterbank, Ma N, Green P, Barker J, Coy A. Exploiting correlogram structure for robust speech recognition with multiple speech sources. Speech Commun. 2007;
%       Available at: http://www.dcs.shef.ac.uk/~ning/resources/gammatone/

%       LIN implementation based on NSL toolbox
%       Powen Ru, Multiscale Multirate Spectro-Temporal Auditory Model\
%       Available at: https://isr.umd.edu/Labs/NSL/Software.htm

% The script generates an amplitude modulated tone signal, which is then processed by the 4 simulated areas:
%                         A1 and R (core)
%                         Slow and Fast (belt)
% It computes strongest osciallation frequency across output channel closest to carrier frequency for all 4 areas, the Vector Strength of these oscillations and the mean firing rates.

% Author: Isma Zulfiqar
% Modified by Sebastian Ayala Ruano

function [EE1, EE2, EE3, EE4] = WC_FiringRate_CortAud(size, lim1_frqax, lim2_frqax, time, fs, rate, depth, carr_frq, sound)

% close all; clear; clc; 
tic;
disp(['Modeling spectro-temporal processing in 2 core (A1,R) and 2 belt (Slow, Fast) areas of the AC:' newline])

%% Initializing Frequency axis
Size = size;  % Spatial (Frequency) size of array - number of units
lim1_frqaxis = lim1_frqax;
lim2_frqaxis = lim2_frqax;
F = MakeErbCFs(lim1_frqaxis,lim2_frqaxis,Size);  % ERB frequency axis

disp(['Model spans ' num2str(lim1_frqaxis) ' Hz to ' num2str(lim2_frqaxis) ' Hz on tonotopic axis...'])

%% Generating input sound 
duration = time; % second
Fs = fs;
mod_rate = rate; % Hz
mod_depth = depth; % from 0 to 1
carrier_freq = carr_frq; % Hz

disp(['Input signal is ' num2str(carrier_freq) ' Hz tone, amplitude modulated at ' num2str(mod_rate) ' Hz...'])

s = sound;
t = [0:(duration*Fs)-1]'/Fs;
%s = (1+(mod_depth*sin(2*pi*mod_rate*t -pi/2))) .* sin(2*pi*carrier_freq*t);
figure(1); plot(t,s); xlabel('Time'); ylabel('Amplitude'); title('Sound waveform');

% =============================================================================================
% Peripheral Processing Stage
% =============================================================================================
  
%% Applying gammatone filterbank
disp('Passing the sound through gammatone filterbank ... ')

bm_upsampled = gammatoneFast(s,F,Fs); 
color='k'; figureSpec(1:length(bm_upsampled),F,squeeze(bm_upsampled'),2,'Gammatone filter output',color);

%% Lateral Inhibitory Network
disp('Running tonotopic derivative and half-wave rectification ...')

y2_h = bm_upsampled(:,Size);
y3_h = 0; y4 = []; y5 =[]; v5 = [];
for ch = Size-1:-1:1
    y2 = bm_upsampled(:,ch);
    
    % ===== Lateral inhibition (masked by higher
    % (frequency) spatial response)
    
    y3 = y2 - y2_h;
    y2_h = y2;
    
    % ==== half wave rectification --- > y4
    y4(:,ch) = max(y3,0);
    
    
    % temporal integration window ---> y5
    y5 = filter(1, [1 -0.9845],y4(:,ch));
    v5(:,ch) = decimate(y5,64);
    
end

F_new = F(2:Size-1);
bm = y4(:,2:Size-1);
color='k'; figureSpec(1:64:length(bm_upsampled),F_new,squeeze(v5(:,2:end)'),3,'Temporally integrated Spectrogram (LIN output)',color);

% converting LIN output to high reslution model input
Last = floor(1000 * (length(s)/Fs));  %last time in computation in ms
DelT = 1000/Fs;  % dt - time resolution of the simulations in ms
time = 0:DelT:Last;
step = DelT/1000; % ms
timei = 0:step:Last/1000;  % time axis with DelT for simulations (in s.)
FS = 1/step;
bm = bm(ceil(linspace(1,length(s),length(timei))),:);
bm = 1E3*bm;

Stim = bm';

[~,channel] = min(abs(F-carrier_freq)); % computing channel with CF closes to carrier freq
    
% =============================================================================================
% Cortical Processing Stage
% =============================================================================================
    
%% inizialize the network and define the inputs
disp('Initializing network...')

EE1 = zeros(Size-2, length(time)); IN1 = EE1; EEresp1 = EE1; INresp1 = EE1;
EE2 = zeros(Size-2, length(time)); IN2 = EE2; EEresp2 = EE2; INresp2 = EE2;
EE3 = zeros(Size-2, length(time)); IN3 = EE3; EEresp3 = EE3; INresp3 = EE3;
EE4 = zeros(Size-2, length(time)); IN4 = EE4; EEresp4 = EE4; INresp4 = EE4;
input1 = zeros(1,length(time)); input2 = EE1; input3 = EE1;

%% Stimulus
disp('Initial input is only for excitatory units, inhibitory input is set to zero...')

P = Stim(:,1:length(time));  % input to excitatory network
color='k'; figureSpec(timei,F_new,Stim,4,'Model Input', color); 
Q=zeros(Size-2,length(time));  % input to the IN population is set to zero

%% Model Parameters
disp('Setting model parameters and checking for stability...')

Xsyn1 = 20*(-5:5); % filter width
Xsyn2 = 20*(-5:5);
Xsyn3 = 20*(-2:2);
Xsyn4 = 20*(-2:2);

% naka rushton constants
% AI
thExc1 = 60;
thIn1 = 80;
max1 = 100;

% R
thExc2 = 60;
thIn2 = 80;
max2 = 100;

% slow
thExc3 = 60;
thIn3 = 80;
max3 = 100;

% fast
thExc4 = 60;%20;
thIn4 = 80;%40;
max4 = 100;%12.5;

% maximum synaptic strength for active transient response
EEgain = 1.5; 
EIgain = 1.3;
IIgain = 1.5;
IEgain = EIgain;

% tau - time constant of the network in ms; 
DT1 = 10; % A1
DT2 = 20; % R
DT3 = linspace(300,200,Size-2)'; % Slow, Time constants change across tonotopic axis
DT4 = linspace(3,1,Size-2)'; % Fast, Time constants change across tonotopic axis

% spatial spread, EE = excitatory to excitatory, EI = excitatory to
% inhibitory, IE = inhibitory to excitatory, II = inhibitory to inhibitory

% A1
sigmaEE1 = 40; 
sigmaEI1 = 160; 
sigmaIE1 = sigmaEI1;
sigmaII1 = 10; 
disp('For A1:'); checkSigmas(sigmaEE1, sigmaEI1, sigmaIE1, sigmaII1, thExc1, thIn1, EEgain, EIgain, IEgain, IIgain , max1);

% R
sigmaEE2 = 40;
sigmaEI2 = 160;
sigmaIE2 = sigmaEI2;
sigmaII2 = 10;
disp('For R:'); checkSigmas(sigmaEE2, sigmaEI2, sigmaIE2, sigmaII2, thExc2, thIn2, EEgain, EIgain, IEgain, IIgain , max2);

% Slow
sigmaEE3 = 20;
sigmaEI3 = 80;
sigmaIE3 = sigmaEI3;
sigmaII3 = 10;
disp('For Slow:'); checkSigmas(sigmaEE3, sigmaEI3, sigmaIE3, sigmaII3, thExc3, thIn3, EEgain, EIgain, IEgain, IIgain , max3);

% Fast
sigmaEE4 = 200;
sigmaEI4 = 300;
sigmaIE4 = sigmaEI4;
sigmaII4 = 30; 
disp('For Fast:'); checkSigmas(sigmaEE4, sigmaEI4, sigmaIE4, sigmaII4, thExc4, thIn4, EEgain, EIgain, IEgain, IIgain , max4);


% synaptic weights
synEE1 = EEgain*exp(-abs(Xsyn1)./sigmaEE1); 
synEI1 = EIgain*exp(-abs(Xsyn1)./sigmaEI1);
synII1 = IIgain*exp(-abs(Xsyn1)./sigmaII1); 

synEE2 = EEgain*exp(-abs(Xsyn2)./sigmaEE2); 
synEI2 = EIgain*exp(-abs(Xsyn2)./sigmaEI2);
synII2 = IIgain*exp(-abs(Xsyn2)./sigmaII2); 

synEE3 = EEgain*exp(-abs(Xsyn3)./sigmaEE3); 
synEI3 = EIgain*exp(-abs(Xsyn3)./sigmaEI3);
synII3 = IIgain*exp(-abs(Xsyn3)./sigmaII3); 

synEE4 = EEgain*exp(-abs(Xsyn4)./sigmaEE4); 
synEI4 = EIgain*exp(-abs(Xsyn4)./sigmaEI4);
synII4 = IIgain*exp(-abs(Xsyn4)./sigmaII4); 

%% Simulating response

disp('Simulating responses... ');
%% AI 
smoothing_kernel = [0.5 1 0.5]; % weight of connections from input to single unit
for T = 2:length(time)  %Loop in ms, Euler solution method
    
    input1 = conv(P(:,T),smoothing_kernel,'same')./sum(smoothing_kernel);
    
    EEresp1(:,T) = NeuralConv(synEE1, EE1(:,T-1)) - NeuralConv(synEI1, IN1(:,T-1)) + input1;
    EEresp1(:,T) = (EEresp1(:,T).*(EEresp1(:,T) > 0)).^2;
    INresp1(:,T) = NeuralConv(synEI1, EE1(:,T-1)) - NeuralConv(synII1, IN1(:,T-1)) + Q(:,T);
    INresp1(:,T) = (INresp1(:,T).*(INresp1(:,T) > 0)).^2;
    
    EE1(:,T) = EE1(:,T-1) + (DelT/DT1)*(-EE1(:,T-1) + (max1)*EEresp1(:,T)./(thExc1^2 + EEresp1(:,T)));
    IN1(:,T) = IN1(:,T-1) + (DelT/DT1)*(-IN1(:,T-1) + (max1)*INresp1(:,T)./(thIn1^2 + INresp1(:,T)));
end

color='b'; figureSpec(timei,F_new,squeeze(EE1),5,'A1 - Excitatory Response',color);
color='r'; figureSpec(timei,F_new,squeeze(IN1),6,'A1 - Inhibitory Response',color);

fft_EE1 = abs(fft(squeeze(EE1(channel,:))));
[maxValue,indexMax] = max(fft_EE1(2:ceil(length(fft_EE1)/2))); 
osc_EE1 =  ceil((indexMax) * (FS / length(fft_EE1)));% strongest oscialltion in the signal
vs_EE1 = maxValue/fft_EE1(1); % Vector Strength of strongest oscillation
fr_EE1 = fft_EE1(1);
disp(['  For A1, strongest oscillation at unit' num2str(channel) ' is ' num2str(osc_EE1) ' Hz, ' ...
    'with Vector Strength = ' num2str(vs_EE1)]);

%% R
for T = 2:length(time)  %Loop in ms, Euler solution method
    input2 = P(:,T);%(((100-CC_inRatio)/100) * P(:,T)) + ((CC_inRatio/100) * CC_in(:,T));

    
    EEresp2(:,T) = NeuralConv(synEE2, EE2(:,T-1)) - NeuralConv(synEI2, IN2(:,T-1)) + input2;
    EEresp2(:,T) = (EEresp2(:,T).*(EEresp2(:,T) > 0)).^2;
    INresp2(:,T) = NeuralConv(synEI2, EE2(:,T-1)) - NeuralConv(synII2, IN2(:,T-1)) + Q(:,T);
    INresp2(:,T) = (INresp2(:,T).*(INresp2(:,T) > 0)).^2;
    
    EE2(:,T) = EE2(:,T-1) + (DelT./DT2).*(-EE2(:,T-1) + (max2)*EEresp2(:,T)./(thExc2^2 + EEresp2(:,T)));
    IN2(:,T) = IN2(:,T-1) + (DelT./DT2).*(-IN2(:,T-1) + (max2)*INresp2(:,T)./(thIn2^2 + INresp2(:,T)));
end

color='b'; figureSpec(timei,F_new,squeeze(EE2),7,'R - Excitatory Response',color);
color='r'; figureSpec(timei,F_new,squeeze(IN2),8,'R - Inhibitory Response',color);

fft_EE2 = abs(fft(squeeze(EE2(channel,:))));
[maxValue,indexMax] = max(fft_EE2(2:ceil(length(fft_EE2)/2))); 
osc_EE2 =  ceil((indexMax) * (FS / length(fft_EE2)));% strongest oscialltion in the signal
vs_EE2 = maxValue/fft_EE2(1); % Vector Strength of strongest oscillation
fr_EE2 = fft_EE2(1);

disp(['  For R, strongest oscillation at unit' num2str(channel) ' is ' num2str(osc_EE2) ' Hz, ' ...
    'with Vector Strength = ' num2str(vs_EE2)]);

%% Slow
for T = 2:length(time)
    input3 = (3*(EE2(:,T)));% - minVal) / ( maxVal - minVal ));% .* filter1;
    
    EEresp3(:,T) = NeuralConv(synEE3, EE3(:,T-1)) - NeuralConv(synEI3, IN3(:,T-1)) + input3;
    EEresp3(:,T) = (EEresp3(:,T).*(EEresp3(:,T) > 0)).^2;
    INresp3(:,T) = NeuralConv(synEI3, EE3(:,T-1)) - NeuralConv(synII3, IN3(:,T-1)) + Q(:,T);
    INresp3(:,T) = (INresp3(:,T).*(INresp3(:,T) > 0)).^2;
    
    EE3(:,T) = EE3(:,T-1) + (DelT./DT3).*(-EE3(:,T-1) + (max3)*EEresp3(:,T)./(thExc3^2 + EEresp3(:,T)));
    IN3(:,T) = IN3(:,T-1) + (DelT./DT3).*(-IN3(:,T-1) + (max3)*INresp3(:,T)./(thIn3^2 + INresp3(:,T)));
end

color='b'; figureSpec(timei,F_new,squeeze(EE3),9,'Slow - Excitatory Response',color);
color='r'; figureSpec(timei,F_new,squeeze(IN3),10,'Slow - Inhibitory Response',color);

fft_EE3 = abs(fft(squeeze(EE3(channel,:))));
[maxValue,indexMax] = max(fft_EE3(2:ceil(length(fft_EE3)/2))); 
osc_EE3 =  ceil((indexMax) * (FS / length(fft_EE3)));% strongest oscialltion in the signal
vs_EE3 = maxValue/fft_EE3(1); % Vector Strength of strongest oscillation
fr_EE3 = fft_EE3(1);
disp(['  For Slow, strongest oscillation at unit' num2str(channel) ' is ' num2str(osc_EE3) ' Hz, ' ...
    'with Vector Strength = ' num2str(vs_EE3)]);

%% Fast
smoothing_kernel = [0.25 0.5 0.5 1 1 1 0.5 0.5 0.25];% weight of connections from input to single unit
for T = 2:length(time)

    % adding a smoothing kernel (simple 1D)
    input4 = conv(EE1(:,T),smoothing_kernel,'same')./sum(smoothing_kernel);
                
    EEresp4(:,T) = NeuralConv(synEE4, EE4(:,T-1)) - NeuralConv(synEI4, IN4(:,T-1)) + input4;
    EEresp4(:,T) = (EEresp4(:,T).*(EEresp4(:,T) > 0)).^2;
    INresp4(:,T) = NeuralConv(synEI4, EE4(:,T-1)) - NeuralConv(synII4, IN4(:,T-1)) + Q(:,T);
    INresp4(:,T) = (INresp4(:,T).*(INresp4(:,T) > 0)).^2;
    
    EE4(:,T) = EE4(:,T-1) + (DelT./DT4).*(-EE4(:,T-1) + (max4)*EEresp4(:,T)./(thExc4^2 + EEresp4(:,T)));
    IN4(:,T) = IN4(:,T-1) + (DelT./DT4).*(-IN4(:,T-1) + (max4)*INresp4(:,T)./(thIn4^2 + INresp4(:,T)));   
end

color='b'; figureSpec(timei,F_new,squeeze(EE4),11,'Fast - Excitatory Response',color);
color='r'; figureSpec(timei,F_new,squeeze(IN4),12,'Fast - Inhibitory Response',color);

fft_EE4 = abs(fft(squeeze(EE4(channel,:))));
[maxValue,indexMax] = max(fft_EE4(2:ceil(length(fft_EE4)/2))); 
osc_EE4 =  ceil((indexMax) * (FS / length(fft_EE4)));% strongest oscialltion in the signal
vs_EE4 = maxValue/fft_EE4(1); % Vector Strength of strongest oscillation
fr_EE4 = fft_EE4(1);
disp(['  For Fast, strongest oscillation at unit' num2str(channel) ' is ' num2str(osc_EE4) ' Hz, ' ...
    'with Vector Strength = ' num2str(vs_EE4)]);

disp(['Finished! Time Elapsed = ' num2str(toc) ' sec']);

end

