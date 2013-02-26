% binsize
binsize_ms = 10;
% window for psth
window_ms = [-100 500];

% sample spiketimes
spiketimes = {	0:100:1000, ...
					0:10:1000, ...
					-99:20:0, ...
					[], ...
					[1100 1200 1300], ...
					-2000:500:2000, ...
				};
		
% test PSTH
[P, b] = psth(spiketimes, binsize_ms, window_ms);
figure(1)
subplot(211)
bar(b, P, 1);


% now test raster
subplot(212)
rasterplot(spiketimes, window_ms)


figure

plotoptions.timelimits = window_ms;
plotoptions.psth_binwidth = binsize_ms;
rasterpsthmatrix(spiketimes, plotoptions)