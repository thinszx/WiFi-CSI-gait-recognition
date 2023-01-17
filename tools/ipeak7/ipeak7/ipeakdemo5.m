function ipeakdemo5
% Demonstration script for iPeak function. It generates a test signal
% consisting of four highly-overlapped peaks, adds random noise, runs 
% ipeak, then prints out the actual and the measured peak positions, 
% heights, widths and areas. Each time you run it you get a different
% sample of noise. 
%   T. C. O'Haver, August 2011

increment=1;
Noise=2;
x=1:increment:501;
% For each simulated peak, specify the amplitude, position, and width
% and the overall random noise level.
amp=[100 200 300 400];  % Amplitudes of the peaks  (Change if desired)
pos=[175 225 275 325];   % Positions of the peaks (Change if desired)
wid=[50 50 50 50];   % Widths of the peaks (Change if desired)
% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
      % Create a series of peaks of different x-positions
      A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks
      % alternatively, A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
      % Assembles actual parameters into ActualPeaks matrix: each row = 1
      % peak; columns are Peak #, Position, Height, Width, Area
      ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k).*wid(k)]; 
      p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
y=z+Noise.*randn(size(z));  % Optionally adds random noise
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
demodata=[x' y']; % Assembles x and y vectors into data matrix
disp('-----------------------------------------------------------------')
disp('In this demo the peaks are so highly overlapped that only one or two of')
disp('the highest peaks yield distinct maxima that are detected by iPeak. The') 
disp('height, width, and area estimates are highly inaccurate because of this')
disp('overlap. The normal peak fit function (''N'' key) would be useful in this')
disp('case, but it depends on iPeak for the number of peaks and for the initial')
disp('guesses, and so it would fit only the peaks that were found and numbered.')  
disp(' ')
disp(['Signal to noise ratio of smallest peak = ' num2str(min(ActualPeaks(:,3))./Noise)])
disp('         Peak #    Position      Height       Width        Area')
ActualPeaks
iPeakResults=ipeak(demodata,0,100,8.2069e-007,15,15,248.5,353,0);
SizeResults=size(iPeakResults);
error=zeros(SizeResults(1),6);

disp(' ')
disp('To help in this case, pressing the ''H'' key will activate the peak ')
disp('sharpen function that decreases peak width and increases peak height')
disp('of all the peaks, making it easier for Findpeaks to detect and number ')
disp('them for use by the peakfit function. Note: peakfit fits the original')
disp('unmodified peaks; the sharpening is used only to help locate the peaks.')
disp(' ')
disp('Peak shape (1-8): 1 ')
disp('Number of trials: 1')
figure(2)
PeakfitResults=peakfit(demodata,251,352,4,1,1,1,[178 60 226 56 276 57 324 46],0); 
SizeResults=size(PeakfitResults);
error=zeros(SizeResults(1),5);
for n=1:SizeResults(1),
    indexn=val2ind(ActualPeaks(:,2),PeakfitResults(n,2));
    error(n,:)=100.*abs(ActualPeaks(indexn,:)-PeakfitResults(n,:))./ActualPeaks(indexn,:);
    error(n,1)=indexn;
    PeakfitResults(n,1)=indexn;
end
if SizeResults(1)==1,
    AveragePercentErrors=error;
else
    AveragePercentErrors=mean(error);
end
PeakfitResults
AveragePercentErrors(1)=0;
Peakfit_Average_Percent_Errors=AveragePercentErrors
disp(' ')
disp('The plot of the peakfit results are shown in Figure 2')
disp('The built-in peakfit function (''N'') key gives much more accurate')
disp('values for the height, width, and area of these peaks')
disp(' ')
disp('End of demo. ')
% ----------------------------------------------------------------------
function g = gaussian(x,pos,wid)
%  gaussian(X,pos,wid) = gaussian peak centered on pos, half-width=wid
%  X may be scalar, vector, or matrix, pos and wid both scalar
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x-pos)./(0.6005615.*wid)).^2);