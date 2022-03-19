function varargout = plotCurveAndCI(curveStruct, varargin)
%------------------------------------------------------------------------
% TytoLogy:PlotTools
%------------------------------------------------------------------------
% plot frequency tuning curve
%------------------------------------------------------------------------
% assumes curveStruct is in format:
%		curveStruct.fname			filename
% 		curveStruct.spikeCount	cell array of spike counts/trial at each x
% 		curveStruct.xdata			stimulus freqs
%		curveStruct.xlabel		label for x axis
%		curveStruct.window		time window [tstart tend] in ms used for analysis
% 		curveStruct.mean			mean values at each freq
% 		curveStruct.std			std. dev. at each freq
% 		curveStruct.mean_ci		cell array of 95% conf. intervals for mean
% 		curveStruct.median		median spike counts at each freq
% 		curveStruct.median_ci	cell array of 95% conf interval for median
%
%  'MEAN'   plots mean
%  'MEDIAN' plot median
%  'FIG'    plots in given figure
%  'AX'     plots in given axes
%------------------------------------------------------------------------
% See Also: computeFTC, computeRLF, opto
%------------------------------------------------------------------------
 		
%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 28 Mar 2019  (SJS) 
%	- adapted from plotRLF.m from OptoAnalysis
% Revisions:
% 
%------------------------------------------------------------------------

%---------------------------------------------------------------------
% Defaults
%---------------------------------------------------------------------
dataToPlot = 'MEAN';
figH = [];
axH = [];

%---------------------------------------------------------------------
% Parse inputs
%---------------------------------------------------------------------
if nargin > 1
	argIndx = 1;
	while argIndx <= (nargin - 1)
		switch upper(varargin{argIndx})
			case {'MEAN', 'AVERAGE', 'AVG'}
				dataToPlot = 'MEAN';
				argIndx = argIndx + 1;
			case {'MEDIAN', 'MED'}
				dataToPlot = 'MEDIAN';
				argIndx = argIndx + 1;
         case {'FIG', 'FIGH', 'FIGHANDLE'}
            if isgraphics(varargin{argIndx+1}, 'figure')
               figH = varargin{argIndx+1};
            else
               error('argument to %s is not a figure handle!', ...
                           varargin{argIndx});
            end
				argIndx = argIndx + 2;
         case {'AX', 'AXH', 'AXES'}
            if isgraphics(varargin{argIndx+1}, 'axes')
               axH = varargin{argIndx+1};
            else
               error('argument to %s is not an axes handle!', ...
                           varargin{argIndx});
            end
				argIndx = argIndx + 2;
			otherwise
				error('%s: unknown option %s', mfilename, varargin{argIndx});
		end
	end
end

% compute confidence intervals
cimatrix = zeros(length(curveStruct.xdata), 2);
for l = 1:length(curveStruct.xdata)
	if strcmpi(dataToPlot, 'MEAN')
		cimatrix(l, :) = curveStruct.mean_ci{l}';
	else
		cimatrix(l, :) = curveStruct.median_ci{l}';
	end
end

% create or set figure
if isempty(figH) && isempty(axH)
   % if no figure or axes provided, create a figure
   figH = figure;
elseif isempty(figH) && ~isempty(axH)
   % if no figure and axes provided, use parent as figure
   figH = get(axH(1), 'Parent');
   % make current
   figure(figH);
else
   % make provided figure current
   figH = figure(figH);
end

% create or set axes
if isempty(axH)
   % create single axes in current figure
   figure(figH);
   axH = axes(figH);
else
   % set current axes to specified handle
   axes(axH);
end


% plot mean
if strcmpi(dataToPlot, 'MEAN')
	ploterrea(curveStruct.xdata, curveStruct.mean, cimatrix);
	ylabel('Mean Spike Count');
% plot median
else
	ploterrea(curveStruct.xdata, curveStruct.median, cimatrix);
	ylabel('Median Spike Count');	
end
% x label
xlabel(axH, curveStruct.xlabel)
grid on
	
if nargout
	varargout{1} = figH;
	if nargout >= 2
		varargout{2} = cimatrix;
	end
end

