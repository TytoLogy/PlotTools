function varargout = rasterpsth(varargin)
% [psth_bar_handle, raster_handle] =  
% 		rasterpsth(spiketimes, P, options)



%---------------------------------------------------------------------
% setup
%---------------------------------------------------------------------
% without args, return default opts struct
if nargin == 0
	varargout{1} = default_rasterpsth_opts;
	return
else
	% initialize opts
	opt = default_rasterpsth_opts;
end
if nargin >= 2
	spiketimes = varargin{1};
	P = varargin{2};
end
% if opts were provided, merge with default opts
if nargin == 3
	opt = merge_fields(opt, varargin{3});
end

% plot PSTH on left y axis
yyaxis left
B = bar(P.bins, P.counts, 1, 'EdgeColor', opt.barcolor, ...
	      'FaceColor', opt.barcolor);
set(gca, 'TickDir', 'out');
set(gca, 'YColor', opt.tickColor);
if ~isempty(opt.xrange)
	xlim(opt.xrange);
end
if ~isempty(opt.psth_yrange)
	ylim(opt.psth_yrange);
end
box off
ylabel('Spike Count')
xlabel('ms');

% plot rasters on right y axis
yyaxis right
R = rasterplot(spiketimes, P.range, opt.tickSymbol, ...
	             opt.tickSize, opt.tickColor);
set(R, 'YColor', opt.tickColor);
if ~isempty(opt.raster_yrange)
	ylim(opt.raster_yrange);
end
if ~isempty(opt.xrange)
	xlim(opt.xrange);
end
if ~isempty(opt.title)
	title_noInterp(opt.title);
end
ylabel('trial');

varargout{1} = B;
varargout{2} = R;
varargout{3} = opt;


end


function opt = default_rasterpsth_opts
	% default plot options
	opt.barcolor = 0.7*[1 1 1];
	opt.tickColor = [0 0 0];
	opt.tickSize = 22;
	opt.tickSymbol = '.';
	opt.onsetLineColor = 'g';
	opt.offsetLineColor = 'r';
	opt.title = '';
	opt.xrange = [];
	opt.psth_yrange = [];
	opt.raster_yrange = [];
end