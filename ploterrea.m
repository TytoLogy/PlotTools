function varargout = ploterrea(x, y, errlohi)
%------------------------------------------------------------------------
% H = ploterrea(x, y, errlohi)
%------------------------------------------------------------------------
% PlotTools toolbox
%------------------------------------------------------------------------
% 
% Given vectors x, y and NX2 matrix (or vector) errlohi of error/stddev
% measurements, plot graph of y vs x with area bounded by errlohi (or +/-
% errlohi if errlohi is a vecotr).
%  
% Returns [ha hl] where ha is handle to area plot axis and hb is line plot
% axis
% 
%------------------------------------------------------------------------
% Input Arguments:
%	x			vector of x coordinates to plot 
%	y			vector of y coordinates to plot 
%	errlohi		vector of error coordinate to plot as area 
%					(same length as x and y vectors)
% 								-OR-
% 					N X 2 matrix of - and + bounds of error area
%
% Output Arguments:
%	[ha hl] handles to area and line plots
%------------------------------------------------------------------------
% See also: psth, rasterpsth_matrix
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: ??? (SJS)
%
% Revisions:
%	8 May, 2017 (SJS):	added to PlotTools toolbox from spike spectra 
%								analysis program, cleaned up code, added comments
%	15 Jun 2020 (SJS):	explicitly set color of line plot to blue
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% input checks
%------------------------------------------------------------------------
if nargin ~= 3
	error('%s: needs x, y and error values as inputs!', mfilename);
end


%------------------------------------------------------------------------
% force x and y to column vector
%------------------------------------------------------------------------
x = force_col(x);
y = force_col(y);
nx = length(x);
ny = length(y);
%------------------------------------------------------------------------
% convert errlohi to column vector or matrix if needed
% assume if ncols in errlohi == 2, that it is ok
%------------------------------------------------------------------------
[nrows, ncols] = size(errlohi); %#ok<ASGLU>
% if ncols > 2, force to column vector using transpose
if ncols > 2
	errlohi = errlohi';
elseif ncols == 1
	errlohi = force_col(errlohi);
end
%------------------------------------------------------------------------
% make sure lengths match
%------------------------------------------------------------------------
nerr = length(errlohi);
if nerr ~= nx || nerr ~= ny || nx ~= ny
	error('%s: x, y, errlohi must be same length', mfilename);
end

%------------------------------------------------------------------------
% matrix for error area
%------------------------------------------------------------------------
if isvector(errlohi)
	errcoords = [ (y-errlohi) (y+errlohi)];
elseif ismatrix(errlohi)
	errcoords = [ (y-errlohi(:, 1)) (y+errlohi(:, 2))];
end

% plot area
ha = area(x, errcoords);
set(ha(1), 'LineStyle', 'none');
set(ha(1), 'FaceColor', 'none');
set(ha(2), 'LineStyle', 'none');
set(ha(2), 'FaceColor', 0.75 * [1 1 1]);

% plot line
hold on
hb = plot(x, y, 'Color', 'b');
hold off

% return handles
if nargout
	varargout{1} = [ha hb];
end

