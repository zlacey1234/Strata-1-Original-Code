% This function will probe various candidate contact cutoff lengths, and
% produce a plot of the peak in the cos \alpha distribution as a function
% of cutoff length.  It will be up to the user to use this plot to select a
% proper cutoff length given that plot, and extract actual contacts using
% that length.

% INPUTS: tracks, dt (for computing relative velocities), rList

% OUTPUTS: cos \alpha distributions for each candidate length, plot of cos
% \alpha peak as a function of candidate length

function [dists,peaks] = probeContactLength(tracks, rList, dt, tStart, tEnd,plotOn)

dists = cell(length(rList),1);
peaks = zeros(length(rList),1);
cosines_bins = linspace(-1,1,21);

for rCutoff = rList
    disp('Analyzing contacts for r-cutoff length')
    [~,cosines] = getContacts(tracks, rCutoff, dt, tStart, tEnd);
    cosine_list = [];
    for i = 1:size(cosines,1)
        for t = 1:length(cosines{i})
            cosine_list = [cosine_list; cosines{i}{t}];% build up the full list of cosines
        end
    end
    p_cos = hist(cosine_list,cosines_bins);
    p_cos = p_cos / trapz(cosines_bins,p_cos); % normalize
    peaks(rList == rCutoff) = p_cos(cosines_bins == 0);
    if plotOn
        figure;
        plot(cosines_bins,p_cos,'.-');
        title(num2str(rCutoff));
    end
end

figure;
plot(rList,peaks,'.')
end