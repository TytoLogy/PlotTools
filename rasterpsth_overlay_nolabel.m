function varargout = rasterpsth_overlay_nolabel(spiketimes, bin_ms, ...
					xrange, yrange, barcolor, tickColor, tickSize, ...
					tickSymbol, tstr)
% [psthvalues, bins, psth_bar_handle, raster_handle] =  
% 		rasterpsth_overlay(spiketimes, bin_ms, ...
% 					xrange, yrange, barcolor, tickColor, tickSize, ...
% 					tickSymbol, tstr)
% no labels drawn on axes

[P, bins] = psth(spiketimes, bin_ms, xrange);
yyaxis left
B = bar(bins, P, 1, 'EdgeColor', barcolor, 'FaceColor', barcolor);
set(gca, 'TickDir', 'out');
set(gca, 'YColor', tickColor);
box off

yyaxis right
R = rasterplot(spiketimes, xrange, tickSymbol, tickSize, tickColor);
set(R, 'YColor', tickColor);
ylim(yrange);
if ~isempty(tstr)
	title_noInterp(tstr);
end

varargout{1} = P;
varargout{2} = bins;
varargout{3} = B;
varargout{4} = R;
