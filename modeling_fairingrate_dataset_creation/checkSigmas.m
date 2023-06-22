%% Stability conditions for the model
function checkSigmas(sigmaEE, sigmaEI, sigmaIE, sigmaII, thExc, thIn, EEgain, EIgain, IEgain, IIgain, n)

for i=1:length(sigmaEE)
% condition 1: inhibition spread is larger than excitation
if ~((sigmaEI(i) == sigmaIE(i)) && (sigmaIE(i) > sigmaEE(i)))
    error('The sigmaIE/sigmaEI should be greater than sigmaEE.')
end

% condition 2: remove unstable states

if ~((n *( (2 * EEgain * sigmaEE(i)) - (2 * IEgain * sigmaIE(i)))) < thExc && ...
        (n *( (2 * EIgain * sigmaEI(i)) - (2 * IIgain * sigmaII))) > thIn)
    error(['Try changing values of parameter ' num2str(i)]);
end

end
    disp('Condition for activity spread satisfied by sigma values.')
    disp('Condition for stability of response satisfied by set parameters.')
end