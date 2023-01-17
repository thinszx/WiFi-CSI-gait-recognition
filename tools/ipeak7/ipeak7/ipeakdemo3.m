function ipeakdemo3
% Demonstration script for iPeak function. It generates a test signal
% consisting of several peaks, adds random noise, runs ipeak, then prints out 
% the actual and the measured peak positions, heights, widths and areas. 
% Each time you run it you get a different sample of noise. You can easily
% evaluate the accuracy of the measurements because the actual peak
% parameter values in this simulation are always integers.
%   T. C. O'Haver, September 2011
increment=1;
x=[1:increment:1200];

% For each simulated peak, compute the amplitude, position, and width
amp=[1 2 3 4 5 6 7 8 9 10];  % Amplitudes of the peaks  (Change if desired)
pos=[100 200 300 400 500 600 700 800 900 1000];   % Positions of the peaks (Change if desired)
wid=[10 15 20 25 30 35 40 45 50 55];   % Widths of the peaks (Change if desired)
Noise=.1; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positions
      % A(k,:)=exp(-((x-pos(k))./(0.6006.*wid(k))).^2); % Gaussian peaks
      A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.7*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=round(WidthPoints/2);  % SmoothWidth should be roughly equal to 1/2 the peak width (in points)
FitWidth=round(WidthPoints/2); % FitWidth should be roughly equal to 1/2 the peak widths(in points)

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW.
% (You can change theses values if desired).
MeasuredPeaks=ipeak(demodata,0,AmpThreshold,SlopeThreshold,SmoothWidth,FitWidth,ActualPeaks(1,2),200);


% Compare MeasuredPeaks to ActualPeakst
disp('-----------------------------------------------------------------')
disp(['Signal to noise ratio of smallest peak=' num2str(min(amp)./Noise)])
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
MeasuredPeaks(:,1:5)

SizeResults=size(MeasuredPeaks);
merror=zeros(SizeResults(1),5);
for n=1:SizeResults(1),
    indexn=val2ind(ActualPeaks(1:5,2),MeasuredPeaks(n,2));
    merror(n,:)=100.*abs(ActualPeaks(indexn,:)-MeasuredPeaks(n,1:5))./ActualPeaks(indexn,:);
    merror(n,1)=MeasuredPeaks(n,1);
    % MeasuredPeaks(n,1)=indexn;
end
AveragePercentError=mean(merror)
disp('Demonstration of overlapping Lorentzians peaks, without a background.')
disp('Overlap of peaks causes significant errors in peak height, ')
disp('width, and area. Hint: turn OFF the Autozero mode (T key) and use the Normal ')
disp('curve fit (N key) or Multiple curve fit (M key) with peak shape 2 (Lorentzian).') 
disp('End of demo.')

% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);