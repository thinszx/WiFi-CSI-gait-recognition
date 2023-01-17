% Comparison of the accuracy of peak measurement between findpeaks and
% peakfit.
increment=1;
x=1:increment:500;

% For each simulated peak, compute the amplitude, position, and width
amp=[1 2 3 4];  % Amplitudes of the peaks  (Change if desired)
pos=[100 200 300 400];   % Positions of the peaks (Change if desired)
wid=[40 30 20 10];   % Widths of the peaks (Change if desired)
Noise=0.3; % Amount of random noise added to the signal. (Change if desired) 

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
y=z+Noise.*randn(size(z));  % Optionally adds random noise
plot(x,y)
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=0.7*WidthPoints^-2; % Formula for estimating value of SlopeThreshold
AmpThreshold=0.05*max(y);
SmoothWidth=WidthPoints; 
FitWidth=WidthPoints;

% call findpeaks and peakfit 
P=findpeaksG(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
% For Lorentzian peaks, use the folloiwng line instead
% P=findpeaksL(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,3);
[FitResults,MeanFitError]=peakfit(demodata,0,0,4,1,0,0,0,0,0,0);

% Compute errors by subtracting ActualPeaks
NumPeaks= max(P(:,1));
if NumPeaks==4;
  FindpeaksError=100.*(P(1:4,2:5)-ActualPeaks(1:4,2:5))./ActualPeaks(1:4,2:5);
plot(x,y,'.r')
hold on
for peak=1:NumPeaks
    text(P(peak,2),P(peak,3),['Peak ' num2str(P(peak,1))])
end
hold off
PeakfitError=100.*(FitResults(1:4,2:5)-ActualPeaks(1:4,2:5))./ActualPeaks(1:4,2:5);
disp('Average percent errors for all peaks')
disp('      Position       Height      Width        Area')
PercentFindpeaksError=mean(abs(FindpeaksError))
PercentPeakfitError=mean(abs(PeakfitError))
RatioOfFindpeaksToPeakfitErrors=mean(PercentFindpeaksError./PercentPeakfitError)
disp('  ')
end
% To plot actual peak heights vs measured heights by findpeaks:
% plotit(ActualPeaks(:,3),FitResults(:,3),1)
% To plot actual peak heights vs measured heights by peakfit:
% plotit(ActualPeaks(:,3),P(:,3),1)