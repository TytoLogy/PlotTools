binsize_ms = 10;
window_ms = [-100 100];


spiketimes = {	2:2:10, ...
				3:3:9, ...
				10:5:50, ...
				[], ...
				0:10:100, ...
				-24:6:-1, ...
				};
		
% use old PSTH
[P, b] = psth(spiketimes, binsize_ms, window_ms);
figure(1)
subplot(211)
bar(b, P, 1);
ylim([0 15])

%--------------------------------------------------------------------
% figure out max and min time for psth.
%--------------------------------------------------------------------
if length(window_ms) == 2
	% set mintime and maxtime from window
	mintime_ms = window_ms(1);
	maxtime_ms = window_ms(2);
else
	% use 0 for mintime_ms, single value from window_ms as maxtime
	mintime_ms = 0;
	maxtime_ms = window_ms;
end

%--------------------------------------------------------------------
% create bins vector from mintime, maxtime and binsize
%--------------------------------------------------------------------
bins = mintime_ms:binsize_ms:maxtime_ms;
nbins = length(bins);

%--------------------------------------------------------------------
% use histc to get histogram values using bins as the "edges"
% of the histogram bins
%--------------------------------------------------------------------
Htrial = zeros(ntrials, nbins);
for trial = 1:ntrials
	% need to trap empty matrices
	if isempty(spiketimes{trial})
		Htrial(trial, :) = zeros(1, nbins);
	else
		Htrial(trial, :) = histc(spiketimes{trial}, bins);
	end
	
end
% and compute final H value
H = sum(Htrial);

subplot(212)
bar(bins, H, 1)
ylim([0 15])


