function varargout = title_noInterp(tstr, varargin)

H = title(tstr, 'Interpreter', 'none', varargin{:});

if nargout
	varargout{1} = H;
end

