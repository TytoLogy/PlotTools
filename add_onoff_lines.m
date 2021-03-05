function Lh = add_onoff_lines(axH, onsetX, offsetX, varargin)

onsetLineColor = 'g';
offsetLineColor = 'r';
lineWidth = [];

if ~isempty(varargin)
	argI = 1;
	while argI <= length(varargin)
		switch upper(varargin{argI})
			case {'ONSET', 'ONSETCOLOR', 'ONSETLINE', 'ONSETLINECOLOR'}
				onsetLineColor = varargin{argI + 1};
				argI = argI + 2;
			case {'OFFSET', 'OFFSETCOLOR', 'OFFSETLINE', 'OFFSETLINECOLOR'}
				offsetLineColor = varargin{argI + 1};
				argI = argI + 2;
			case {'LINEWIDTH', 'LINE_WIDTH', 'LINEWEIGHT', 'LINE_WEIGHT'}
				lineWidth = varargin{argI + 1};
				argI = argI + 2;
			otherwise
				error('%s: unknown option %s', mfilename, varargin{argI});
		end
	end
end

Lh = [0 0];
yl = ylim(axH);
% onset line
Lh(1) = line(axH, onsetX.*[1 1], [0 yl(2)], 'Color', onsetLineColor);
% offset line
Lh(2) = line(axH, offsetX.*[1 1], [0 yl(2)], 'Color', offsetLineColor);
if ~isempty(lineWidth)
	set(Lh(1), 'LineWidth', lineWidth);
	set(Lh(2), 'LineWidth', lineWidth);
end
