function varargout = yyhistogram(Lvals, Lctr, Lci, Rvals, Rctr, Rci, ...
	                               histbins, varargin)
% [H_L, H_R] = yyhistogram(Lvals, Lctr, Lci, Rvals, Rctr, Rci, ...
% 	                               histbins, <Llabel, Rlabel>)
% 
% plots two histograms, one for each y axis (Left, Right)
% 
% Lvals     left axis data vector [# left values, 1]
% Lctr      left axis histogram bin centers [nbins, 1]
% Lci       left axis histogram confidence intervals [2, 1]


yyaxis left
H_L = histogram(Lvals, histbins);
yrange = ylim;
line(Lctr*[1 1], yrange, 'Color', [0.00,0.45,0.74], 'LineWidth', 3);
line(Lci(1)*[1 1], yrange, 'Color', [0.00,0.45,0.74]);
line(Lci(2)*[1 1], yrange, 'Color', [0.00,0.45,0.74]);
if ~isempty(varargin)
	ylabel(varargin{1}, 'Interpreter', 'none');
end

yyaxis right
H_R = histogram(Rvals, histbins);
yrange = ylim;
line(Rctr*[1 1], yrange, 'Color', [0.85,0.33,0.10], 'LineWidth', 3);
line(Rci(1)*[1 1], yrange, 'Color', [0.85,0.33,0.10]);
line(Rci(2)*[1 1], yrange, 'Color', [0.85,0.33,0.10]);

if length(varargin) > 1
	ylabel(varargin{2}, 'Interpreter', 'none');
end

if nargout
	varargout{1} = H_L;
	varargout{2} = H_R;
end