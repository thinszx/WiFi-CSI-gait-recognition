function ipeakdemo1
% Demonstration script for iPeak function. It generates a test signal
% consisting of several peaks, adds random noise, runs ipeak, then prints out 
% the actual and the measured peak positions, heights, widths and areas. 
% Each time you run it you get a different set of peaks. You can easily
% evaluate the accuracy of the measurements because the actual peak
% parameter values in this simulation are always integers.
%   T. C. O'Haver, September 2011
format short g
increment=1;
x=[1:increment:4000];
 
% For each simulated peak, compute the amplitude, position, and width
amp=round(5.*randn(1,38));  % Amplitudes of the peaks  (Change if desired)
pos=[200:100:3900];   % Positions of the peaks (Change if desired)
wid=40.*ones(size(pos));   % Widths of the peaks (Change if desired)
Noise=0.2; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
  if amp(k)>.2,  % Keep only those peaks above a certain amplitude
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
  end; 
end 
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
y=y+5.*gaussian(x,0,4000); % Optionally adds a broad background signal
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix

% Initial values of variable peak detection parameters
WidthPoints=mean(wid)/increment; % Average number of points in half-width of peaks
SlopeThreshold=.0004; % Formula for estimating value of SlopeThreshold
AmpThreshold=1;
SmoothWidth=56;  % SmoothWidth should be roughly equal to 1/2 the peak width (in points)
FitWidth=40; % FitWidth should be roughly equal to 1/2 the peak widths(in points)

% Now call iPeak, with specified values of AmpT, SlopeT, SmoothW, and FitW,
% with AUTOZERO=1 (ON) and the top window showing the first peak. Once
% iPeak is running, you can change any of these using 
% keystroke commands (press "K" to see a list of commands).
MeasuredPeaks=ipeak(demodata,0,AmpThreshold,SlopeThreshold,SmoothWidth,FitWidth,ActualPeaks(1,2),120,1);

% Compare MeasuredPeaks to ActualPeaks
disp('-----------------------------------------------------------------')
disp(['Signal to noise ratio of smallest peak = ' num2str(min(ActualPeaks(:,3))./Noise)])
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
disp('Demonstration of autozero background correction, for separated, narrow ')
disp('peaks on a large baseline. Hint: Turn on the Autozero mode (T key),')
disp('adjust the zoom setting so that the peaks are shown one at a time in')
disp('the upper window, then press the P key to display the peak table.')
disp('Jump to the next/previous peaks using the Spacebar/Tab keys.')

disp('End of demo.')
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);
% ----------------------------------------------------------------------
function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
