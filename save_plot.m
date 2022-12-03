function varargout = save_plot(figH, outputFormat, outputPath, varargin)
%------------------------------------------------------------------------
% figH = save_plot(figH, outputFormat, outputPath, <<plot options>>)
%------------------------------------------------------------------------
% TytoLogy:PlotTools Toolbox
%------------------------------------------------------------------------
% 
% saves figure pointed to by figH. outputFormat determines format,
% outputPath is location of output file
% 
% name of file will be taken from name of figure 
% so, if 'Name' property of figH os 'aplot' and saveFormat is 'PNG', output
% file will by aplot.png
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	figH	figure handle
%	outputFormat	output format for plot.
% 				'PNG'		      .png, 300 dpi resolution
% 				'PDF'		      .pdf
%           'PDF_EXPORT'   .pdf, using exportgraphics()
% 				'FIG'		      .fig (saves MATLAB figure)
%  outputpath	output directory
%
% Output Arguments:
%	figH	figure handle
%
%------------------------------------------------------------------------
% See also: SpikeData.plotAllData
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created:8 July, 2020 (SJS)
%
% Revisions:
% 9 July 2020 (SJS): moved to PlotTools toolbox
% 2 Dec 2022 (SJS): added pdf_export option
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% default resolution
png_resolution = 600;

% check inputs
if nargin < 3
	error('save_plot: need 3 input args, fig handle, format and path');
end

% convert outputFormat to cell if it is a simple char string
if ~iscell(outputFormat)
   if ischar(outputFormat)
      outputFormat = {outputFormat};
   end
end

% if format is invalid, throw error
for s = 1:length(outputFormat)
   if ~any(strcmpi(outputFormat{s}, {'PNG', 'PDF', 'PDF-EXPORT', 'FIG'}))
   	error('save_plot: unknown save format %s', outputFormat{s});
   end
end

% get figure name from handle and append to output path
pname = fullfile(outputPath, get(figH, 'Name'));
% if Name is empty, try getting name from outputPath
if isempty(pname)
   warning('%s: figure Name is empty, trying name from outputPath', ...
               mfilename);
   [~, pname] = fileparts(outputPath);
   if isempty(pname)
      error('%s: cannot set file name', mfilename);
   end
end

% save plot in desired format(s)
for s = 1:length(outputFormat)
   switch upper(outputFormat{s})
	   case 'FIG'
		   savefig(figH, pname, 'compact');
	   case 'PDF'
		   print(figH, pname, '-dpdf');
      case 'PDF-EXPORT'
         exportgraphics(figH, [pname '.pdf'], 'ContentType', 'vector'); 
	   case 'PNG'
         if isempty(varargin)
            print(figH, '-dpng', sprintf('-r%d', png_resolution), pname);
         else
            print(figH, '-dpng', pname, varargin{:});
         end            
   end
end

if nargout
	varargout{1} = figH;
end



