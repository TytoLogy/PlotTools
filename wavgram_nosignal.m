function varargout = wavgram_nosignal(x, n, Fs, window, overlap, maxfreq, lowdBlimit)
% plot spectrograms
%  Original version myspecgram() by Paul Kienzle; 
%  modified by Sean Fulop March 2002
%  modified by Sharad Shanbhag Jan 2021


% if only the window length is given, generate hanning window
if length(window) == 1, window = hanning(window); end

% should be extended to accept a vector of frequencies at which to
% evaluate the fourier transform (via filterbank or chirp
% z-transform)
if length(n)>1
 error('specgram doesn''t handle frequency vectors yet') 
end

% compute window offsets
win_size = length(window);
if (win_size > n)
 n = win_size;
 warning('specgram fft size adjusted---must be at least as long as frame')
end
step = win_size - overlap;

% build matrix of windowed data slices
offset = 1:step:(length(x)-win_size);
S = zeros (n, length(offset));
for i=1:length(offset)
  S(1:win_size, i) = x(offset(i):offset(i)+win_size-1) .* window;
end

% compute fourier transform
STFT = fft(S);

% extract the positive frequency components
if rem(n,2)==1
 ret_n = (n+1)/2;
else
 ret_n = n/2;
end

STFT = STFT(1:ret_n, :);
f = (0:(ret_n-1))*Fs/n;
t = offset/Fs;

spts = floor(2:n*maxfreq/Fs);

% magnitude of STFT
STFTmag = abs(STFT(spts,:));
% normalize so max magnitude will be 0 db
STFTmag = STFTmag/max(max(STFTmag));
% clip everything below lowdBlimit
STFTmag = max(STFTmag, 10^(lowdBlimit/10));

% imagesc(t, f(2:n*maxfreq/Fs), 20*log10(STFTmag)); axis xy; colormap(flipud(gray));
% display as an indexed grayscale image showing the log magnitude of the STFT, 
% i.e. a spectrogram; the colormap is flipped to invert the default setting,
% in which white is most intense and black least---in speech science we want
% the opposite of that.

% Plot Signal

% subplot(211);
% % plot((1000/Fs) * (0:(length(x)-1)), x, 'k');
% % plot(1000*t(1, end), x(1, end), 'k')
% H(1) = gca;
% H(1) = 0;

% Plot Spectrogram
subplot(211);
pcolor(1000*t, 0.001*f(spts), 20*log10(STFTmag));
axis xy;
% Grey colormap
% colormap(gca, flipud(gray));
colormap(gca, jet);
shading interp;
H(2) = gca;

varargout{1} = H;
varargout{2} = STFTmag; 
varargout{3} =  f;
varargout{4} =  t;

end