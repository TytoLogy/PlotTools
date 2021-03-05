function Hpsth = draw_psth(b, h, varargin)
% Hpsth = draw_psth(b, h, barcolor)
% 
% with bins (edges) b and histogram counts h, plot histogram using bar()
% with edge and face color set to (optional) barcolor;
%
% returns handle from  bar() MATLAB function
%
% b and h are usually returned by the histogram() MATLAB function

if isempty(varargin)
	barcolor = 0.5*[1 1 1];
else
	barcolor = varargin{1};
end

Hpsth = bar(b, h, 1, 'EdgeColor', barcolor, 'FaceColor', barcolor);
