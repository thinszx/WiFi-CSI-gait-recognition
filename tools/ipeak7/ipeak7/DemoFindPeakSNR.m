% A simple self-contained demonstration of the findpeaksnr function (line 54)
% applied to noisy synthetic data set consisting of a random number of narrow 
% peaks.  Each time you run this, a different set of peaks is generated.
% Calls the fundpeaksnr function, which must be in the Matlab path.
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and 
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% Tom O'Haver (toh@umd.edu). Version 3.1 February 2013
% You can change the signal characteristics in lines 9-16
format short g
format compact
increment=2;
x=[1:increment:8000];

% For each simulated peak, compute the amplitude, position, and width
pos=[200:200:7800];   % Positions of the peaks (Change if desired)
amp=round(5.*randn(1,length(pos)));  % Amplitudes of the peaks  (Change if desired)
wid=50.*ones(size(pos));   % Widths of the peaks (Change if desired)
Noise=.2; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.1,  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      % A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k)*wid(k)]; 
      p=p+1;
  end; 
end 
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
sqrtnoise=sqrt(abs(z)).*randn(size(z));
propnoise=abs(z).*randn(size(z));
sz=stdev(sqrtnoise);
pz=stdev(sqrtnoise);
nz=Noise.*sqrtnoise./sz;
y=z+nz;  % Adds square-root signal-dependent random noise
% y=y+.1.*randn(size(z));  % Adds constant random noise
% y=y+5.*gaussian(x,0,4000); % Optionally adds a broad background signal
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
% demodata=[x' y']; % Assembles x and y vectors into data matrix

figure(1);plot(x,y,'r')  % Graph the signal in red
title('Detected peaks are numbered. Peak table is printed in Command Window')

% Initial values of variable parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.5*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.02*max(y);
SmoothWidth=round(WidthPoints);  % SmoothWidth should be roughly equal the peak width (in points)
FitWidth=round(WidthPoints); % FitWidth should be roughly equal to the peak widths (in points)

% Lavel the x-axis with the parameter values
xlabel(['SlopeThresh. = ' num2str(SlopeThreshold) '    AmpThresh. = ' num2str(AmpThreshold) '    SmoothWidth = ' num2str(SmoothWidth) '    FitWidth = ' num2str(FitWidth) ])

% Find the peaks
tic;
Measuredpeaks=findpeaksnr(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
ElapsedTime=toc;
PeaksPerSecond=length(Measuredpeaks)/ElapsedTime;

% Display results
disp('---------------------------------------------------------')
disp(['SlopeThreshold = ' num2str(SlopeThreshold) ] )
disp(['AmpThreshold = ' num2str(AmpThreshold) ] )
disp(['SmoothWidth = ' num2str(SmoothWidth) ] )
disp(['FitWidth = ' num2str(FitWidth) ] )
disp(['Speed = ' num2str(round(PeaksPerSecond)) ' Peaks Per Second' ] )
disp('          Peak #     Position      Height      Width        SNR')
Measuredpeaks  % Display table of peaks
figure(1);text(Measuredpeaks(:, 2),Measuredpeaks(:, 3),num2str(Measuredpeaks(:,1)))  % Number the peaks found on the graph
clear PercentErrors
if length(ActualPeaks)==length(Measuredpeaks),
    PercentErrors(:,1:4)=100.*(ActualPeaks(:,1:4)-Measuredpeaks(:,1:4))./ActualPeaks(:,1:4);
    PercentErrors(:,1)=Measuredpeaks(:,1)
    AverageAbsolutePercentErrors=mean(abs(100.*(ActualPeaks(:,1:4)-Measuredpeaks(:,1:4))./ActualPeaks(:,1:4)))
end