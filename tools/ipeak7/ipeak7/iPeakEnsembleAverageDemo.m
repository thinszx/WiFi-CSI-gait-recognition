% A self-contained demonstration of the iPeak function. In this example, 
% the signal contains a repeated pattern of two overlapping Gaussian peaks
% of width 12, with a 2:1 height ratio. These patterns occur a random
% intervals, and the noise level is about 10% of the average peak height.
% Using iPeak's ensemble average function (Shift-E), the patterns can be
% averaged and the signal-to-noise ratio significantly improved. 
format short g
format compact
increment=1;
x=[1:increment:18000];

% For each simulated peak, compute the amplitude, position, and width
pos=[200:50:17800];  % Positions of the peaks (Change if desired)
amp=round(10.*randn(1,length(pos)));  % Amplitudes of the peaks  (Change if desired)
wid=12.*ones(size(pos));   % Widths of the peaks (Change if desired)
Noise=1; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplitude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>9,  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2)+.5.*exp(-((x-pos(k)-wid(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k)*wid(k)]; 
      p=p+1;
  end; 
end 
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Adds constant random noise
% y=z+Noise.*sqrtnoise(z);  % Adds signal-dependent random noise
% y=y+5.*gaussian(x,0,4000); % Optionally adds a broad background signal
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
% demodata=[x' y']; % Assembles x and y vectors into data matrix
% 
% figure(1);plot(x,y,'r')  % Graph the signal in red
% title('Detected peaks are numbered. Peak table is printed in Command Window')

% Initial values of variable parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.5*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=3;
SmoothWidth=15;  % SmoothWidth should be roughly equal the peak width (in points)
FitWidth=15; % FitWidth should be roughly equal to the peak widths (in points)

% % Lavel the x-axis with the parameter values
% xlabel(['SlopeThresh. = ' num2str(SlopeThreshold) '    AmpThresh. = ' num2str(AmpThreshold) '    SmoothWidth = ' num2str(SmoothWidth) '    FitWidth = ' num2str(FitWidth) ])

% Find the peaks
tic;
Measuredpeaks=ipeak([x;y],0,AmpThreshold,0.00015625,SmoothWidth,FitWidth,300,60,0);
ElapsedTime=toc;
% Zoom in on average peak
AveragePeakIndex=val2ind(Measuredpeaks(:,3),15);
AveragePeakPosition=Measuredpeaks(AveragePeakIndex,2);
Measuredpeaks=ipeak([x;y],0,AmpThreshold,0.00015625,SmoothWidth,FitWidth,AveragePeakPosition,73,0);
PeaksPerSecond=length(Measuredpeaks)/ElapsedTime;

% Display results
disp('---------------------------------------------------------')
disp('Demonstration of iPeak''s ensemble average function (Shift-E)')
disp(['SlopeThreshold = ' num2str(SlopeThreshold) ] )
disp(['AmpThreshold = ' num2str(AmpThreshold) ] )
disp(['SmoothWidth = ' num2str(SmoothWidth) ] )
disp(['FitWidth = ' num2str(FitWidth) ] )
disp(['Number of peaks detected = ' num2str(length(Measuredpeaks)) ] )
disp(['Speed = ' num2str(round(PeaksPerSecond)) ' Peaks Per Second' ] )
disp(' ')
disp('Press Shift-E to compute the ensemble average of all peaks, then look at ' )
disp('the plot of the ensemble average in Figure 2 and compare its signal-to-noise ' )
disp('ratio to that of the original signal. The theoretical signal-to-noise improvement ')
disp(['is the square root of the number of peaks averaged = ' num2str(sqrt(length(Measuredpeaks))) ] )
disp(' ')
disp('To curve-fit the ensemble-averaged signal, type this: ')
disp('[parameters,error]=peakfit([1:length(EnsembleAverage);EnsembleAverage],40,60,2)')
disp(' ')
disp('To fit one of the peaks in the original signal, type this: ')
disp('[FitResults,MeanFitError]=peakfit([x;y],AveragePeakPosition,70,2,1,0,10)')
disp(' ')
disp('Compare those two fitting results. The actual repeated signal pattern')
disp('consists of two Gaussians, 12 points apart, width 12, with a 2:1 height ratio.')
disp('The results for the averaged signal pattern should be much closer.')