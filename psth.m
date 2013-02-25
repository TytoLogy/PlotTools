function varargout = psth(spiketimes, binsize_ms, window_ms, Fs)
%------------------------------------------------------------------------
% [H, bins] = psth(spiketimes, binsize_ms, maxtime_ms)
% [H, bins] = psth(spiketimes, binsize_ms, [mintime maxtime])
% [H, bins] = psth(spiketimes, binsize_ms, maxtime_ms, Fs)
%------------------------------------------------------------------------
% Given an input cell array of spiketimes (in units of bins or milliseconds) 
% and desired bin width in milliseconds, compute a peri-stimulus time 
% histogram (PSTH) and return it in vector H.  
% 
% maxtime_ms is used to limit the length of the computed PSTH to something 
% reasonable (e.g., duration of stimulus, or some other user-defined 
% analysis window).
% 
% bins contains vector of time bins; this is useful for plotting the psth using 
% the bar() plot command.
% 
% default units for spiketimes is milliseconds;
% if units are in bins (i.e., samples), call psth with 
% Fs (sampling rate in samples/second) specified
% 
% Examples:
%		generate psth using 10 ms bins, 1000 ms max psth time with spiketimes
%	 	s given in units of milliseconds
%	 	[h, b] = psth(s, 10, 1000);
% 
%	 	% plot the psth using the bar() command and specify no space between
%	 	% histogram bars using '1' for the bar width option
%	 	bar(b, h, 1)
% 
%		generate psth using 20 ms bins, 500 ms max psth time with spiketimes
%	 	s given in units of samples, and sampling rate of 40000 samples/sec
%	 	[h, b] = psth(s, 20, 500, 40000);
% 
%		generate psth using 10 ms bins, with a psth window [mintime maxtime]
%		of [-200 500] ms and spiketimes given in units of msec.
%	 	[h, b] = psth(s, 20, [-200 500]);
		
% 
%------------------------------------------------------------------------
% Input:
%	spiketimes		cell array of spike time data
% 						spike time data in units of milliseconds (default), 
% 						or bins (a.k.a. samples) that are 
% 						referenced to start time of data acquisition
% 						if spiketimes are in units of samples, Fs should be
% 						provided as an input argument (see below), or output
% 						will be strange and you will be unhappy.
% 
%	binsize_ms		psth bin size in milliseconds
%
%	The time window for the psth may be specified in two ways.
%
% 	A single value specifying the maximum time (time of final psth bin) that
% 	will use 0 ms as the minimum or start time for the psth:
%		maxtime_ms		maximum time length of psth window (ms)
% 
% 	or a [1X2] vector containing the minimum and maximum values for the psth
%		window_ms	[mintime maxtime] min and maximum time length 
% 						of psth window (ms)
%
% 	Optional:
% 		Fs					sample rate (samples/second) for spiketimes data
%							should be provided when spiketimes are in units 
%							of bins or samples.  
%							DO NOT provide this if spiketimes are already 
%							in units of milliseconds!!!!!!  also, all other values
%							window, maxtime, binsize must be in milliseconds!
% 
% Output:
% 	H					Nbins X 1  vector containing spike counts / time bin
%	bins				Nbins X 1  vector with time (ms) bins for plotting H using
%	Htrial			Nbins X Ntrials  matrix with individual histograms per
%						trial
%						
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 9 June, 2010 (SJS) from older psth.m used in RFAnalyze package
%
% Revisions:
%	9 June, 2010 (SJS):
% 	 -	updated comments and documentation
% 	 -	changed names of input args by adding _ms to the end to make
% 	 	them more descriptive
%	 -	added Fs input variable
%	7 February, 2012 (SJS): updated email address
%	25 Feb, 2013 (SJS):	
%		MAJOR reworking to allow for spec. of mintime so that	neg. spiketimes 
%		will work
% 	 -	uses MATLAB builtin histc to compute histograms
%------------------------------------------------------------------------

%--------------------------------------------------------------------
% # trials are the # of elements in the spiketimes cell array
% (if not a cell, ntrials is 1, and convert spiketimes into
% cell for compatibility with rest of algorithm)
%--------------------------------------------------------------------
if iscell(spiketimes)
	ntrials = length(spiketimes);
else
	ntrials = 1;
	spiketimes = {spiketimes};
end

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
% preallocate Htrial matrix to store results for each spiketimes vector
Htrial = zeros(ntrials, nbins);
% loop through trials...
for trial = 1:ntrials
	% need to trap empty matrices
	if isempty(spiketimes{trial})
		% if empty, assign zeros...
		Htrial(trial, :) = zeros(1, nbins);
	else
		% if not, compute histogram using bins as "edges"
		Htrial(trial, :) = histc(spiketimes{trial}, bins);
	end
end
% ...and compute final H value
H = sum(Htrial);

%--------------------------------------------------------------------
% assign outputs (or plot)
%--------------------------------------------------------------------
if nargout == 0
	% just plot the psth, since no output args were asked for
	plot(bins, H);
	return
end
% otherwise, assign output vars.
if any(nargout == [1 2 3])
	varargout{1} = H';
end
if any(nargout == [2 3])
	varargout{2} = bins';
end
if nargout == 3
	varargout{3} = Htrial';
end


%% %%%%%%%%%%%%%%%%%%%%%%  OLD (pre 25 Feb 2013) METHOD %%%%%%%%%%%%%%%%%%%%%
%{
%--------------------------------------------------------------------
% figure out lengths of vectors and # of bins
%--------------------------------------------------------------------
% determine number of trials in spiketime data - this will be equal
% to the number of elements in the cell array spiketimes
ntrials = length(spiketimes);

%--------------------------------------------------------------------
% check if Fs was provided as input arg
%--------------------------------------------------------------------
if nargin == 4
	if Fs > 0
		% convert spiketimes from samples to ms
		for trial = 1:ntrials
			spiketimes{trial} = bin2ms(spiketimes{trial}, Fs);
		end
	else
		warning('%s: Fs must be > 0, ignoring Fs value %.2f', mfilename, Fs);
	end
end

%--------------------------------------------------------------------
% check if window was given as a pair of numbers or a singleton
%--------------------------------------------------------------------
if length(window_ms) == 1
	% set maxtime from window
	maxtime_ms = window_ms;
	% pre-allocate the psth array, H
	H = zeros(ceil(maxtime_ms/binsize_ms), 1);
	nbins = length(H);
end


%--------------------------------------------------------------------
% fill the bins for PSTH
%--------------------------------------------------------------------
% loop through trials (or sweeps)
for trial = 1:ntrials
	% get the times and # of elements for the current trial
	sdata = spiketimes{trial};
	sbins = length(sdata);
	% loop through the # of spike times in this trial
	for s = 1:sbins
		% determine the psth bin in which this spike falls
		% by dividing the spike time in milliseconds by the psth binsize in
		% ms... 1 is added in order to account for Matlab indexing starting
		% at 1
		psthbin = floor(sdata(s)/binsize_ms) + 1 ;
		% increment the spike count at the respective psth bin by 1
		if psthbin <= nbins
			H(psthbin) = H(psthbin)+1;
		else
			H(end) = H(end)+1;
		end
	end
end

%--------------------------------------------------------------------
% assign outputs
%--------------------------------------------------------------------
if nargout == 2
	% NEW BINS COMPUTATION (25 Feb 2013, SJS)
	bins = (binsize_ms * (0:(nbins - 1)))';
	% OLD BINS COMPUTATION 
	%	bins = (1:binsize_ms:maxtime_ms)';
end
%}

