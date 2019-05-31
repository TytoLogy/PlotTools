function [H, plotopts] = rasterpsthmatrix(Spikes, varargin) %#ok<*STOUT>
%------------------------------------------------------------------------
% [H, plotopts] = rasterpsthmatrix(Spikes, plotopts)
%------------------------------------------------------------------------
% PlotTools toolbox
%------------------------------------------------------------------------
%	plots matrix of rasters and psth  
%------------------------------------------------------------------------
% Input Arguments:
% 	Spikes		Cell matrix of {Nrows, Ncols}, each element of which contains
% 					a cell array of spike time vectors (NOT interspike intervals!)
% 					in milliseconds 
% 
%	Optional inputs:
%		plotopts		Plot options structure
% 			timelimits: [0 1000]
% 			psth_binwidth: 5
%			psth_color:	[0 0 1]
% 			raster_tickmarker: '|'
% 			raster_ticksize: 12
%			raster_color: [0 0 1]
% 			horizgap: 0.0500
% 			vertgap: 0.0550
% 			plotgap: 0.0125
% 			filelabel: '768_4_1q_Bat_1_output.txt'
% 			columnlabels: {7x1 cell}
% 			rowlabels: {4x1 cell}
% 			idlabel: 'Unit 11'
%			stimulus_times: {nrows, ncols}[nstim, 2]
%			plot_titles: {nrows, ncols} of strings used for titles
%			plot_titles_color: {nrows, ncols} of colors used for titles
%			psthtext: {nrows, ncols} structs with following fields for 
% 							placing text on PSTH plots
% 					.text		strings to place on plot
% 					.xy		[x y] locations for text
% 					.color	string or [r g b] color for text
% 					.size		size of text (default is 10 pt)
% 					.font		font name (string)
% 			rastertext: {nrows, ncols} same as psthtext
% 			stimulus_times_plot:	2 or 3 plots stimulus onset/offset (not
% 										sure what 2 or 3 means....?????)
% 			stimulus_times:	{nrows, ncols} [nstims, 2], where columns
% 									are onset and offset times in seconds
% 			stimulus_on_color {nstims}
%			stimulus_off_color {nstims} colors for stimulus on/off
% 			stimulus_onoff_pct: [value in %] amount from ylimit(2) to start 
% 										onset/offset lines. default is 90;
%			psth_ylimits:	[min max] y-axis limits for psth plots
% 
% Output Arguments:
% 	H				structure of handles to plots figure
% 		figure		figure handle for plots
% 		rasters		axes handles for raster plots
% 		psths			axes handles for psth plots
% 		
%	plotopts		plot options structure, with handles updated
% 		timelimits: [0 1000]
% 		psth_binwidth: 5
% 		raster_tickmarker: '.'
% 		raster_ticksize: 12
% 		horizgap: 0.0500
% 		vertgap: 0.0550
% 		plotgap: 0.0125
% 		filelabel: '768_4_1q_Bat_1_output.txt'
% 		columnlabels: {7x1 cell}
% 		rowlabels: {4x1 cell}
% 		idlabel: 'Unit 11'
% 		plotwidth: 0.0857
% 		plotheight: 0.0837
% 		pos1: {4x7 cell}
% 		pos2: {4x7 cell}
% 
%------------------------------------------------------------------------
% See also: rasterplot, psth
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 	7 July, 2011 (SJS)
%
% Revisions:
%	13 July, 2011 (SJS)
% 	 -	functionalized script
% 	 -	removed Nrows and Ncols as inputs (redundant)
% 	 -	updated comments/documentation
% 	15 July, 2011 (SJS)
% 	 -	made changes to incorporate automatic scaling of psth yaxis to
% 	 	uniform value (global max of all psths)
% 	 -	added documentation
% 	 -	added psthlimits to plotoptions output struct
% 	 -	changed H from cell array to struct with descriptive fields
%	15 Feb, 2013(SJS)
%	 -	fixed default/varargin settings of plotopts
%	 -	cleaned up code, added comments
%	12 Mar 2013 (SJS)
%		added check of max y value for psth - if 0, resets to 1 to avoid
%		ylim() error
%	16 May, 2013 (SJS)
%		- moved rasterplot function into this m file, renamed to raster
% 			this was done since some slightly different capabilities are
% 			required and I didn't want to break anything that relies
% 			on rasterplot it its current form
%		- added raster_color to plotopts
%	20 May, 2013 (SJS)
% 	 -	if element in Spikes is either a character or NaN, that 
% 		plot will be skipped (left empty).  
% 	 -	checks if plotwidth or height are included in plotopts
% 	 -	added widthscale and heightscale to add scaling to all plots
%	 -	added vertoffset and horizoffset to shift all plots
%	31 May, 2013 (SJS)
% 	 -	added stimulus_times
%	3 June 2013 (SJS)
%	 -	added psth_color, tickdir, ticklabelsize, [x,y]labelsize options
%	1 Dec 2014 (SJS)
%	 - added plot title plotting, envelopes
%	23 Mar 2015 (SJS)
%	 - reworked plot title implementation
%	10 Oct 2017 (SJS): 
% 	 - added functionality for different colored on/off stimulus lines
% 		(i.e., for use with sound + opto stimulation
%	 - user can specify psth y-axis limits using plotopts.psth_ylimits
%------------------------------------------------------------------------
% TO DO:
%
% plotting characters for rasters messes up time alignment with psth
% when font size > 12 pt.  use this approach instead:
% plot(a, 1*ones(size(a)), '.')
% hold off
% close all
% plot(a, 1*ones(size(a)), '.', 'MarkerSize', 2)
% hold on, plot(a, 2*ones(size(a)), '.', 'MarkerSize', 2), hold off
% ylim([0 3])
% hold on, plot(1000*rand(size(a)), 3*ones(size(a)), '.', 'MarkerSize', 2), hold off
% ylim([0 4])
% get(gca, 'Children')
%------------------------------------------------------------------------

%-----------------------------------------------------------
% default plotopts
%-----------------------------------------------------------
defaultopts = struct( ...
	'timelimits',				[0 1000]			, ...
	'psth_binwidth',			5					, ...
	'psth_color',				[0 0 1]			, ...
	'raster_tickmarker',		'|'				, ...
	'raster_ticksize',		10					, ...
	'raster_color',			[0 0 1]			, ...
	'horizgap',					0.05				, ...
	'vertgap',					0.055				, ...
	'plotgap',					0.0125			, ...
	'widthscale',				1					, ...
	'heightscale',				1					, ...
	'vertoffset',				0					, ...
	'horizoffset',				0					, ...
	'xlabel',					'msec'			, ...
	'xlabelsize',				8					, ...
	'ylabel',					''					, ...
	'ylabelsize',				8					, ...
	'ticklabelsize',			8					, ...
	'tickdir',					'out'				, ...
	'plot_titles',				{{}}				, ...
	'plot_titles_color',		{{}}				, ...
	'envelopes',				{{}}				, ...
	'envelope_dt',				[]					, ...	
	'psthtext',					{{}}				, ...
	'rastertext',				{{}}				, ...
	'psth_ylimits',			[]					...
); 

ONCOLOR = [0 0.75 0]; 
OFFCOLOR = [1 0 0]; 

%-----------------------------------------------------------
% get dimensions of Spikes cell matrix
%-----------------------------------------------------------
[Nrows, Ncols] = size(Spikes);
fprintf('... %s: creating figure with %d rows, %d columns\n', ...
															mfilename, Nrows, Ncols)
% check size
if ~Nrows || ~Ncols
	error('%s: error in size of Spikes cell matrix', mfilename)
end

%-----------------------------------------------------------
% initialize plotopts struct if not passed in
%-----------------------------------------------------------
if isempty(varargin)
	plotopts = defaultopts;
else
	plotopts = varargin{1};
	pfields =	fieldnames(defaultopts);
	% assign missing options from defaultopts to plotopts
	for n = 1:length(pfields)
		if ~isfield(plotopts, pfields{n})
			plotopts.(pfields{n}) = defaultopts.(pfields{n});
		end
	end
end

%-----------------------------------------------------------
% set plot title colors
%-----------------------------------------------------------
if ~isempty(plotopts.plot_titles)
	if isempty(plotopts.plot_titles_color) || ...
				all(size(plotopts.plot_titles) ~= size(plotopts.plot_titles_color))
		plotopts.plot_titles_color = cell(Nrows, Ncols);
		for r = 1:Nrows
			for c = 1:Ncols
				plotopts.plot_titles_color{r, c} = [0 0 0];
			end
		end
	end
end

%-----------------------------------------------------------
% check settings for stimulus_times
%-----------------------------------------------------------
if isfield(plotopts, 'stimulus_times')
	if ~isfield(plotopts, 'stimulus_times_plot')
		fprintf('%s:\n', mfilename);
		fprintf('\tstimulus_times provided but stimulus_times_plot not set\n');
		fprintf('\tusing default (3 == draw stim onset/offset on all plots)\n');
		plotopts.stimulus_times_plot = 3;
	end
	if ~isfield(plotopts, 'stimulus_on_color')
		plotopts.stimulus_on_color = ONCOLOR;
	end
	if ~isfield(plotopts, 'stimulus_off_color')
		plotopts.stimulus_off_color = OFFCOLOR;
	end
	if ~isfield(plotopts, 'stimulus_onoff_pct')
		plotopts.stimulus_onoff_pct = 90;
	end
end

%----------------------------------------------------------------------------
% compute plot widths and plot heights if values were not provided in
% plotopts input struct
%----------------------------------------------------------------------------
% logic:
%-----------
% any figure window is measured in dimensionless units, with x-axis (horizontal)
% having range of [0 1] and y-axis (vertical) [0 1].  The upper-left corner of
% the figure window has coordinates (0, 1) and lower-left corner is (1, 0)
% 
% Therefore to have all plots the same width and height, and to have them evenly
% spaced with horizontal gaps of size plotopts.horizgap and vertical gaps of
% plotopts.vertgap, compute the width as follows:
% 
% 	plotwidth = ([window width] - [# of columns + 1] * [horizgap]) / [# of
% 	columns]
% 	
% and the height as follows (noting that the # of plots in the vertical
% dimension is 2 * Nrows, since there are 2 plots (raster & psth) per row):
% 
% 	plotheight = ([window height] - [# of rows + 2] * [vertgap]) / 2*[# of rows]
% 
%----------------------------------------------------------------------------
if ~isfield(plotopts, 'plotwidth') || isempty(plotopts.plotwidth)
	plotopts.plotwidth = plotopts.widthscale * ...
									(1 - ((Ncols+1) * plotopts.horizgap)) / Ncols;
end
if ~isfield(plotopts, 'plotheight') || isempty(plotopts.plotheight)
	plotopts.plotheight = plotopts.heightscale * ...
									(1 - ((Nrows+2) * plotopts.vertgap)) / (2*Nrows);
end

%----------------------------------------------------------------------------
% compute positions for rasters (pos1) and psths (pos2)
%----------------------------------------------------------------------------
% logic:
%-----------
%----------------------------------------------------------------------------
pos1 = cell(Nrows, Ncols);
pos2 = cell(Nrows, Ncols);
for r = 1:Nrows
	ypos1(r) = 1 - r*plotopts.plotheight - (r-1)*plotopts.plotheight ...
					- r*plotopts.vertgap - (r-1)*plotopts.plotgap ...
					- plotopts.vertoffset; %#ok<AGROW>
	% only need to include vertoffset once
	ypos2(r) = ypos1(r) - plotopts.plotheight - plotopts.plotgap; %#ok<AGROW>
	for c = 1:Ncols
		xpos(c) = plotopts.horizgap ...
					 + ((c-1) * (plotopts.plotwidth + plotopts.horizgap)) ...
					 + plotopts.horizoffset; %#ok<AGROW>
		pos1{r, c} = [xpos(c) ypos1(r) plotopts.plotwidth plotopts.plotheight];
		pos2{r, c} = [xpos(c) ypos2(r) plotopts.plotwidth plotopts.plotheight];
	end
end

%----------------------------------------------------------------------------
% In order to scale all the psths to same min max y values, need to
% pre-compute the psths and keep track of min and max values.  do this here
%----------------------------------------------------------------------------
% data will be stored in psth struct with histvals and bins cell matrices
%		histvals{Nrows, Ncols} will hold histogram data
%		bins{Nrows, Ncols} will hold the respective x-axis bins for the histvals
%		maxval(Nrows, Ncols) will hold the max values of histvals for each element
%----------------------------------------------------------------------------
% pre-allocate them here.
psthdata.histvals = cell(Nrows, Ncols);
psthdata.bins = cell(Nrows, Ncols);
psthdata.maxval = zeros(Nrows, Ncols);
% loop through rows and cols of Spikes and plot data
for row = 1:Nrows
	for col = 1:Ncols
		% compute psth

		% if Spikes{row, col} is a character, skip it
		if ~ischar(Spikes{row, col})
			% build psth from spike data and plot using bar() function
			% modified call to use full time limits after updating psth function
			% 25 Feb 2013 (SJS)
			[tmpvals, tmpbins] = psth(	Spikes{row, col}, ...
												plotopts.psth_binwidth, ...
												plotopts.timelimits);
	%%%%% OLD
	% 		[tmpvals, tmpbins] = psth(	Spikes{row, col}, ...
	% 											plotopts.psth_binwidth, ...
	% 											plotopts.timelimits(2));
			psthdata.histvals{row, col} = tmpvals;
			psthdata.bins{row, col} = tmpbins;
			psthdata.maxval(row, col) = max(tmpvals);
		end
	end
end

% find the overall maximum value in psthdata.maxval matrix and build ylimit
% vector
maxPSTHval = max(max(psthdata.maxval));
if maxPSTHval == 0
	warning('%s: maxPSTHval == 0... using 1...');
	maxPSTHval = 1;
end

if isempty(plotopts.psth_ylimits)
	psthdata.ylimits = [0 maxPSTHval];
else
	psthdata.ylimits = plotopts.psth_ylimits;
end

%----------------------------------------------------------------------------
% Now, plot the data
%----------------------------------------------------------------------------
% initialize cell matrices for storing handles to raster plots (handles1) and 
% psths (handles2)
handles1 = cell(Nrows, Ncols);
handles2 = cell(Nrows, Ncols);

% loop through rows and cols of Spikes and plot data
for row = 1:Nrows
	for col = 1:Ncols
		%-------------------------------------------------------
		% First, plot raster for this row and column
		%-------------------------------------------------------
		% if Spikes{row, col} is a character, skip it
		if ~ischar(Spikes{row, col})
			% select subplot location for rasters (pos1)
			subplot('Position', pos1{row, col});
			% store the axes handle returned by raster in the handles2 cell array
			handles1{row, col} = raster(	Spikes{row, col}, ...
														plotopts.timelimits, ...
														plotopts.raster_tickmarker, ...
														plotopts.raster_ticksize, ...
														plotopts.raster_color );
			% turn off xtick labels, and turn off yticks
			set(gca, 'XTickLabel', []);
			set(gca, 'ytick', []);
			% if id label provided in plotopts, display it at col 1, row 1
			if isfield(plotopts, 'idlabel')
				if (col == 1) && (row == 1)
					idstr = plotopts.idlabel;
				else
					idstr = '';
				end
			else
				idstr = '';
			end
			% if columnlabel was given, plot on row 1 plots
			if isfield(plotopts, 'columnlabels')
				if (row == 1)
					colstr = plotopts.columnlabels{col};
				else
					colstr = '';
				end
			else
				colstr = '';
			end
			% if row labels given in plotopts, display on y label of column 1
			% plots
			if isfield(plotopts, 'rowlabels')
				if col == 1
					rowstr = plotopts.rowlabels{row};
				else 
					rowstr = '';
				end
			else
				rowstr = '';
			end
			% check if separator (space) is needed between idstr and colstr
			if isempty(idstr) || isempty(colstr)
				sepstr = '';
			else
				sepstr = ' ';
			end
			titlestr = [idstr sepstr colstr ];

			% plot titles
			% check if plot_titles is empty
			if ~isempty(plotopts.plot_titles)
				if ~isempty(plotopts.plot_titles{row, col})
					if ~isempty(titlestr)
						title({titlestr, plotopts.plot_titles{row, col}}, ...
								'Interpreter', 'none', ...
								'Color', plotopts.plot_titles_color{row, col});
					else
						title(plotopts.plot_titles{row, col}, ...
								'Interpreter', 'none', ...
								'Color', plotopts.plot_titles_color{row, col});	
					end
				end
			end
			if ~isempty(titlestr)
				title(titlestr, 'Interpreter', 'none');			
			end
			
			if ~isempty(rowstr)
				ylabel(rowstr, 'Interpreter', 'none', ...
									'FontSize', plotopts.ylabelsize);			
			end
			
			% text for raster plot
			if isfield(plotopts, 'rastertext')
				if ~isempty(plotopts.rastertext)
					if ~isempty(plotopts.rastertext{row, col})
						addtext(plotopts.rastertext{row, col});
					end
				end
			end
			
			%-------------------------------------------------------
			% then, plot psth
			%-------------------------------------------------------
			% select subplot location for psth (pos2)
			subplot('Position', pos2{row, col});
			% store the axes handle returned by bar in the handles2 cell array
			handles2{row, col} = bar(	psthdata.bins{row, col}, ...
												psthdata.histvals{row, col}, ...
												1, ...
												'EdgeColor', plotopts.psth_color, ...
												'FaceColor', plotopts.psth_color);
			% update time limits to match raster
			xlim(plotopts.timelimits)
			% set ylimits to overall value in psthdata
			ylim(psthdata.ylimits);
			% turn off x tick labels for all but the bottom row 
			if row ~= Nrows
				set(gca, 'XTickLabel', []);
			else
				set(gca, 'FontSize', plotopts.ticklabelsize);
				if isfield(plotopts, 'xtick') && isfield(plotopts, 'xticklabel')
					set(gca, 'Xtick', plotopts.xtick);
					set(gca, 'XTickLabel', plotopts.xticklabel);
				end
			end
			
			% turn off y tick labels for all but the left column
			if col ~= 1
				set(gca, 'ytick', []);
			else
				set(gca, 'FontSize', plotopts.ticklabelsize);
				if isfield(plotopts, 'ytick') && isfield(plotopts, 'yticklabel')
					set(gca, 'Ytick', plotopts.ytick);
					set(gca, 'YTickLabel', plotopts.yticklabel);
				end
			end
			
			% set tick direction
			set(gca, 'TickDir', plotopts.tickdir);
			% turn off box
			set(gca, 'Box', 'off');
			% label the x-axis 'msec' on the lower left psth plot
			if (col == 1) && (row == Nrows)
				xlabel(plotopts.xlabel, 'Interpreter', 'none', ...
								'FontSize', plotopts.xlabelsize);
			end
			
			% set Tick length
			if isfield(plotopts, 'ticklength')
				set(gca, 'TickLength', plotopts.ticklength);
			end
			
			% plot stimulus onset/offset lines if stimulus_times provided
			if isfield(plotopts, 'stimulus_times')
				if any(plotopts.stimulus_times_plot == [2 3])
					if iscell(plotopts.stimulus_times)
						try
							[nstim, ~] = size(plotopts.stimulus_times{row, col});
						catch
							keyboard
						end
						if (nstim > 1) && ...
									( length(plotopts.stimulus_onoff_pct) == 1 )
							plotopts.stimulus_onoff_pct = ...
										plotopts.stimulus_onoff_pct * ones(nstim, 1);
						end
						for t = 1:nstim
							% get ylimits
							L = ylim;
							L(1) = 0.01*plotopts.stimulus_onoff_pct(t)*L(2);
							onset = 1000*plotopts.stimulus_times{row, col}(t, 1);
							offset = 1000*plotopts.stimulus_times{row, col}(t, 2);
							% onset line
							line(onset.*[1 1], L, 'Color', ...
														plotopts.stimulus_on_color{t});
							% offset line
							line(offset.*[1 1], L, 'Color', ...
														plotopts.stimulus_off_color{t});	
						end
					else
						% get ylimits
						L = ylim;
						L(1) = 0.01*plotopts.stimulus_onoff_pct*L(2);
						if any(length(plotopts.stimulus_times) == [1 2])
							%draw onset line
							onset = 1000*plotopts.stimulus_times(1);
							line(onset.*[1 1], L, 'Color', ...
													plotopts.stimulus_on_color);
						end
						if length(plotopts.stimulus_times) == 2
							% draw offset line
							offset = 1000*plotopts.stimulus_times(2);
							line(offset.*[1 1], L, 'Color', ...
													plotopts.stimulus_off_color);
						end
					end
				end
			end	% END OF if isfield(plotopts, 'stimulus_times')
			
			% plot stimulus envelopes
			if ~isempty(plotopts.envelopes)
				% generate time vector
				tvec = plotopts.envelope_dt(row, col) * ...
								((1:length(plotopts.envelopes{row, col})) - 1);
				% plot envelope
				hold on
				h1 = plot(tvec, ...
								0.95*maxPSTHval*...
								(1 + 0.05*plotopts.envelopes{row, col}), 'm'); %#ok<NASGU>
				h2 = plot(tvec, ...
								0.95*maxPSTHval*...
								(1 - 0.05*plotopts.envelopes{row, col}), 'm'); %#ok<NASGU>
				% reset ylimits
				ylim(psthdata.ylimits);
				hold off
			end	% END of isempty(plotopts.envelopes)

			% text for psth plot
			if isfield(plotopts, 'psthtext')
				if ~isempty(plotopts.psthtext)
					if ~isempty(plotopts.psthtext{row, col})
						if isempty(plotopts.psthtext{row, col}.xy)
							plotopts.psthtext{row, col}.xy = ...
										[0.70*plotopts.timelimits(2) 0.80*maxPSTHval];
						end
						addtext(plotopts.psthtext{row, col});
					end
				end
			end

		else
			fprintf('%s: Spikes{%d, %d} skipped\n', mfilename, row, col);
			
		end	%	END OF if ~ischar(Spikes{row, col}) || ~isnan(Spikes{row, col})

	end		% END OF col LOOP
end		% END OF row LOOP

drawnow

% plot file title
if isfield(plotopts, 'filelabel')
	plotopts.annH = annotation(	'textbox', ...
											[0 0 0.2 0.05], ...
											'String', plotopts.filelabel, ...
											'Interpreter', 'none', ...
											'FitBoxToText', 'on', ...
											'HorizontalAlignment', 'left', ...
											'VerticalAlignment', 'middle', ...
											'EdgeColor', 'none' 	);
	set(gcf, 'Name', plotopts.filelabel);
	set(gcf, 'FileName', plotopts.filelabel);
end

% adjust figure size
if isfield(plotopts, 'position')
	set(gcf, 'Position', plotopts.position)
end

% adjust paper orientation
if isfield(plotopts, 'PaperOrientation')
	set(gcf, 'PaperOrientation', plotopts.paperorientation);
end

%-------------------------------------------------------
% set output values
%-------------------------------------------------------
H.figure = gcf;
H.rasters = handles1;
H.psths = handles2;
% return plotopts is requested
if nargout == 2
	plotopts.pos1 = pos1;
	plotopts.pos2 = pos2;
	plotopts.psthlimits = psthdata.ylimits;
end
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
function [H, Hrep] = raster(spiketimes, timeMinMax, ...
										ticksymbol, ticksize, tickcolor, ...
										offset) %#ok<INUSD>
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
%------------------------------------------------------------
% Defaults
%------------------------------------------------------------
TICKASCII = double('|');
TICKSIZE = 10;
TICKCOLOR = 'b';
OFFSET = 0; %#ok<NASGU>

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
	tmpval = cell2mat(spiketimes);
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
ry = nReps+1;
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
% ylim([-1*floor(nReps/20) nReps+1]);
ylim([-2 nReps+1]);
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------


%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
function H = addtext(T)
if ~isstruct(T)
	error('%s: not a text struct!', mfilename);
end
% create text object
H = text;
% set properties

% string
if isfield(T, 'text')
	set(H, 'String', T.text);
end
% position
if isfield(T, 'xy')
	set(H, 'Position', [T.xy(1) T.xy(2) 0]);
end
% color
if isfield(T, 'color')
	set(H, 'Color', T.color);
end
% font size
if isfield(T, 'size')
	set(H, 'FontSize', T.size);
end
% font name
if isfield(T, 'font')
	set(H, 'FontName', T.font);
end
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------
			
