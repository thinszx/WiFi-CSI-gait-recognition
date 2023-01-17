% A simple self-contained demonstration of the findpeak function (line 54)
% applied to noisy synthetic data set consisting of a random number of narrow 
% peaks.  Each time you run this, a different set of peaks is generated.
% Calls the fundpeaks function, which must be in the Matlab path.
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and 
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
% Tom O'Haver (toh@umd.edu). Version 3 February 2013
% You can change the signal characteristics in lines 9-16
format short g
format compact
increment=1;
x=[1:increment:1800];

% For each simulated peak, compute the amplitude, position, and width
pos=[200:50:1780];   % Positions of the peaks (Change if desired)
amp=round(10.*randn(1,length(pos)));  % Amplitudes of the peaks  (Change if desired)
wid=20.*ones(size(pos));   % Widths of the peaks (Change if desired)
Noise=.4; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>9,  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
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
clf
figure(1);plot(x,y,'r')  % Graph the signal in red
title('Detected peaks are numbered. Peak table is printed in Command Window')

% Initial values of variable parameters

SlopeThreshold=.0001;
AmpThreshold=2;
SmoothWidth=27;
FitWidth=18;
% Label the x-axis with the parameter values
xlabel(['SlopeThresh. = ' num2str(SlopeThreshold) '    AmpThresh. = ' num2str(AmpThreshold) '    SmoothWidth = ' num2str(SmoothWidth) '    FitWidth = ' num2str(FitWidth) ])

% Find the peaks
tic;
Measuredpeaks=findpeaksG(x,y,.0001,2,27,18,3);
ElapsedTime=toc;
PeaksPerSecond=length(Measuredpeaks)/ElapsedTime;

% Display results
disp('---------------------------------------------------------')
disp(['SlopeThreshold = ' num2str(SlopeThreshold) ] )
disp(['AmpThreshold = ' num2str(AmpThreshold) ] )
disp(['SmoothWidth = ' num2str(SmoothWidth) ] )
disp(['FitWidth = ' num2str(FitWidth) ] )
disp(['Speed = ' num2str(round(PeaksPerSecond)) ' Peaks Per Second' ] )
disp('         Peak #     Position      Height      Width       Area')
disp(Measuredpeaks)  % Display table of peaks
figure(1);text(Measuredpeaks(:,2),Measuredpeaks(:,3),num2str(Measuredpeaks(:,1)))  % Number the peaks found on the graph

% if length(ActualPeaks)==length(Measuredpeaks),
%     PercentErrors=100.*(ActualPeaks-Measuredpeaks)./ActualPeaks;
%     PercentErrors(:,1)=Measuredpeaks(:,1);
%     AverageAbsolutePercentErrors=mean(abs(100.*(ActualPeaks-Measuredpeaks)./ActualPeaks));
% end