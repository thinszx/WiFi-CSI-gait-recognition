function ipeakdemo6
% Demonstration script for iPeak function. It generates a test signal
% consisting of two peaks, adds additive and proportional random noise, runs
% ipeak, then prints out the actual and the measured peak positions,
% heights, widths and areas. Each time you run it you get a different
% sample of noise. You can easily evaluate the accuracy of the measurements
% because the actual peak parameter values in this simulation are always
% integers.
%   T. C. O'Haver, September 2011
increment=1;
x=[1:increment:500];

% For each simulated peak, compute the amplitude, position, and width
amp=[1 1];  % Amplitudes of the peaks  (Change if desired)
pos=[100 300];   % Positions of the peaks (Change if desired)
wid=[30 60];   % Widths of the peaks (Change if desired)
Noise=0.2; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positio ns
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*propnoise(z)+Noise.*pinknoise(length(z));  % Optionally adds proportional random noise
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.7*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=33;  % Larger-than-usual values of SmoothWidth and FitWidth 
FitWidth=40; % will help deal with the very poor signal-to-noise ratio.

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW.
% (You can change theses values if desired).
MeasuredPeaks=ipeak(demodata,0,AmpThreshold,SlopeThreshold,SmoothWidth,FitWidth,ActualPeaks(1,2),200);

% Compare MeasuredPeaks to ActualPeaks
disp('-----------------------------------------------------------------')
disp('Detection and measurement of two peaks in a very noisy signal')
disp('You can change the peak heights, positions, and widths in lines 14-16,')
disp('the noise level in line 17, and the noise type in line 33.')
disp(['Signal-to-noise ratio of first peak = ' num2str(amp(1)/Noise) ] )
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
MeasuredPeaks

if size(ActualPeaks)==size(MeasuredPeaks),
    PercentErrors=100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks
    AveragePercentErrors=mean(abs(100.*(ActualPeaks-MeasuredPeaks)./ActualPeaks))
end
disp('The accuracy of peak parameter measurement ')
disp('is poor because of the low signal-to-noise ratio.')
disp(' ')
disp('An even better way to measure this signal is to use the peakfit')
disp('function, now that we know the number of peaks and their')
disp('approximate, positions and widths. Press the M key and ')
disp('enter the peak shape (1) and number of trials. Compare to ')
disp('the iPeak results obtained by pressing the P key. ')
disp(' ')
disp('End of script.')

function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);

function g = lorentzian(x,position,width)
% lorentzian(x,position,width) Lorentzian function.
% where x may be scalar, vector, or matrix
% position and width scalar
% T. C. O'Haver, 1988
% Example: lorentzian([1 2 3],2,2) gives result [0.5 1 0.5]
g=ones(size(x))./(1+((x-position)./(0.5.*width)).^2);

function y=propnoise(x)
% Random noise whose amplitude is proportional to the amplitude of signal
% with white power spectrum with mean zero and unit standard deviation,
% equal in length to x.
% Tom O'Haver, 2012
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*propnoise(model);
% plot(model)
% 
z=x.*randn(size(x));
sz=stdev(z);
y=z./sz;

function stddev=stdev(a)
% Octave and Matlab compatible standard deviation function
sa=size(a);
reshape(a,1,length(a));
if sa(1)>sa(2),
  stddev=std(a);
else
  stddev=(std(a'));
end;

function y=whitenoise(x)
% Random noise with white power spectrum with mean zero 
% and unit standard deviation, equal in length to x
% Tom O'Haver, 2008
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*whitenoise(model);
% plot(model)
% 
y=randn(size(x));

function ry=pinknoise(n)
% Random noise with pink (1/f) power spectrum with mean zero 
% and unit standard deviation. n is number of points.
% Tom O'Haver, 2008
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*pinknoise(length(model));
% plot(model)
x=[1:n];
y=randn(size(x));  % Random normally-distributed white noise
% Fourier filter 
fy=fft(y); % Compute Fourier transform of signal y
% Compute filter shape
lft1=[1:(length(fy)/2)+1];
lft2=[(length(fy)/2):length(fy)];
ffilter1=ones(size(lft1))./(sqrt(lft1));
ffilter2=ones(size(lft2))./(sqrt(lft2));
ffilter=[ffilter1,ffilter2];
if length(fy)>length(ffilter), ffilter=[ffilter ffilter(1)];end
ffy=fy.*ffilter(1:length(fy));  % Multiply filter by Fourier transform of signal
ry=real(ifft(ffy)); % Inverse transform to recover filtered signal 'ry'
ry=((ry-mean(ry))./std(ry)); % Normalize to unit standard deviation

function y=sqrtnoise(x)
% Random noise whose amplitude is proportional to the square root of the
% amplitude of signal with white power spectrum with mean zero and unit
% standard deviation, equal in length to x.
% Tom O'Haver, 2012
% Example:
% model=gaussian([1:1000],500,200);
% model=model+.1.*sqrtnoise(model);
% plot(model)
% 
z=sqrt(abs(x)).*randn(size(x));
sz=std(z);
y=z./sz;