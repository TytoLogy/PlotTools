function [H, Hrep] = rasterplot(spiketimes, timeMinMax, ...
											ticksymbol, ticksize, tickcolor, ...
											axesHandle)
%------------------------------------------------------------------------
% H = rasterplot(spiketimes, timeMinMax, 
%						ticksymbol, ticksize, tickcolor, axesHandle)
%------------------------------------------------------------------------
% PlotTools toolbox
%------------------------------------------------------------------------
% 
% Given a cell "vector" (NX1 or 1XN) of spiketimes, e.g.
% 
% 			spiketimes = 
% 				 [1x4 double]
% 				 [1x3 double]
% 							  []
% 				 [1x3 double]
% 				 [1x3 double]
%  
% draw a raster plot of the spiketimes, where each element in the spiketimes
% cell vector corresponds to a vector of spiketimes in milliseconds.
% 
% so, for example, spiketimes{1} = [100 120 200 300]
%	 
% timeMinMax is optional x-axis limits as a [1X2] vector; if not provided
% rasterplot will automatically determine limits
% 
% axes is an optional handle to an axes object - if not provided, 
% rasterplot will create a new one in a new figure
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	spiketimes		cell "vector" (NX1 or 1XN) of vectors containing 
% 						times (NOT interspike intervals!) of spikes
% 	
% 	Optional:
% 		timeMinMax		x-axis limit vector in form [min max]
%		ticksymbol		character symbol for raster ticks (default is '|')
%		ticksize			font size for tick symbol (default is 10)
%		tickcolor		color for raster ticks (default is blue)
% 		axesHandle		handle to axes
% 
% Output Arguments:
% 	H	handle to plot
%	Hrep	handles to reps
%------------------------------------------------------------------------
% See also: psth
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 3 December, 2009 (SJS)
%
% Revisions:
%	7 July, 2011 (SJS)
% 		-	fixed an issue with timeMinMax syntax
% 		-	added ticksymbol input argument to vary symbol
% 		-	added return value of handles to plot
% 	13 July, 2011 (SJS)
% 		-	added ticksize input argument and code to change tick size
%	25 Feb 2013 (SJS): some cleanup
%------------------------------------------------------------------------

%------------------------------------------------------------
% Defaults
%------------------------------------------------------------
TICKASCII = double('|');
TICKSIZE = 10;
TICKCOLOR = 'b';

%------------------------------------------------------------
% some checks on inputs
%------------------------------------------------------------
% check if figure handle was provided at input
if exist('axesHandle', 'var')
	% if so, make sure it is a proper handle, if not, create new figure
	if ishandle(axesHandle)
		H = gca(axesHandle);
	% otherwise create axis		
	elseif isempty(axesHandle)
		H = gca;
	% or throw error
	else
		warning([mfilename ': invalid input axes handle, creating new handle']);
		H = gca;
	end
else
	% otherwise, create axes, save handle to return as output of function
	H = gca;
end

cla

% check if timeMaxMin was specified
if ~exist('timeMinMax', 'var')
	timeMinMax = [0 0];
	maxtimeSearchFlag = 1;
else
	maxtimeSearchFlag = 0;
end
% check if ticksymbol was provided
if ~exist('ticksymbol', 'var')
	ticksymbol = TICKASCII;
end
% check if ticksize was user-specified
if ~exist('ticksize', 'var')
	ticksize = TICKSIZE;
end
% check for tickcolor
if ~exist('tickcolor', 'var')
	tickcolor = TICKCOLOR;
end

%------------------------------------------------------------
% draw plot
%------------------------------------------------------------

% check if spiketimes is a cell
if iscell(spiketimes)
	% if so use length of spiketimes as # of reps
	nReps = length(spiketimes);
else
	% otherwise, nReps = 1 and convert array of spiketimes to cell
	nReps = 1;
	spiketimes = {spiketimes};
end

% create Hrep output if necessary
if nargout > 1
	Hrep = cell(nReps, 1);
end

% first need to find max time if asked by user
if maxtimeSearchFlag
	% convert to matrix
	% debugging sjs 2/23/2015
	% 	if iscolumn(spiketimes)
	% 		tmpval = cell2mat(spiketimes');
	% 	else
	% 		tmpval = cell2mat(spiketimes);
	% 	end
	if iscolumn(spiketimes)
		tmpval = cell2mat(spiketimes);
	else
		tmpval = cell2mat(spiketimes');
	end
	
	% if it's not empty, find max
	if ~isempty(tmpval)
		if isvector(tmpval)
			% if vector, simple max value is ok
			timeMinMax(2) = max(tmpval);
		else
			% otherwise, need overall max
			timeMinMax(2) = max(max(tmpval));
		end
	end
	clear tmpval;
end


% loop through reps
ry = nReps;
for r = 1:nReps
	valid_indices = (spiketimes{r} >= timeMinMax(1)) & ...
							(spiketimes{r} <= timeMinMax(2));
	% find the timestamps in range of timeMinMax(1) and timeMinMax(2)
	ts = spiketimes{r}(valid_indices);	
	% x locations for ticks == spike times
	xlocs = ts;
	% ylocations are set by rep index (r)
	ylocs = (ry)*ones(size(ts));
	% need a row vector of ticks due to a peculiarity of the text() function
	% in Matlab 
	tickchars = char(ticksymbol * ones(length(ts), 1));
	% draw the ticks, return a vector of handles
	h = text(xlocs, ylocs, tickchars, 'Interpreter', 'none', 'HorizontalAlignment', 'center');
	% use the handles vector to set color
	set(h, 'Color', tickcolor);
	set(h, 'FontSize', ticksize);
	% assign to output if needed
	if nargout > 1
		Hrep{r} = h;
	end
	% decrement rep counter
	ry = ry - 1;
end
% set x limit to manual, set limit
xlim('manual')
xlim(timeMinMax);
% set ylimit to manual, set limit
ylim('manual');
ylim([-1*floor(nReps/20) nReps+1]);


	
