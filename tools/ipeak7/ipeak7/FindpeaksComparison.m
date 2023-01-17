% Comparison of the peak detection functions findpeaksG, findpeaksb,
% findpeaksb3, and findpeaksfit for a computer-generated signal with
% multiple peaks on a variable baseline plus variable random noise. You may
% change the lines here marked by <<< to modify the character of the signal
% peaks, baseline, and noise. For example, you can change the peak shape in
% line 19, peak positions, heights, and widths in lines 14-16, the sampling
% density in line 12, baseline shape in line 32, baseline amplitude in line
% 20, noise type in line 39, noise amplitude in line 18, and the baseline
% correction mode in line 23. Adjust the signal, baseline, and noise to
% look like your data, to see which methods result in the lowest peak
% parameter errors. Tom O'Haver  toh@umd.edu    2015
increment=.4; %  <<< Change data density here (total number of points)
x=1:increment:100;
Heights=[5 4 3 2 1]; % <<< Height of each peak above the baseline
Positions=[15 30 50 70 85]; %  <<< Peak positions (between 1 and 100)
Widths=[6 6 6 6 6]; % <<< Width (FWHM) of each peak
Extras=[0 0 0 0 0]; %  <<< Needed only for variable shapes (Voigt, etc)
noise=.1;  %  <<< Change noise amplitude herepeakshape=1;  
peakshape=1;  %  <<< Change peak shape here: 1=Gaussian, 2=Lorentzaian, etc.
BaselineIntensity=2; %  <<< Change baseline amplitude here (0-10)
extra=0; % Needed only for variable shapes (Voigt, etc)
NumTrials=5; % <<< Normally 1 - 10
AUTOZERO=3; % <<< Change baseline correction: 0=none; 1=linear; 2=quadratic; 3=flat.
width=mean(Widths);
NumPeaks=length(Heights);
% Baseline functions
FlatBaseline=ones(size(x));
LinearBaseline=(1-x/max(x));
ExpBaseline=exp(-x./(10*width));
GaussBaseline=gaussian(x,0,20*width);
% Baseline may be FlatBaseline, LinearBaseline, ExpBaseline, GaussBaseline
Baseline=BaselineIntensity.*ExpBaseline; %  <<< Change baseline type here
y=modelpeaks(x,NumPeaks,peakshape,Heights,Positions,Widths,Extras)+Baseline;
WhiteNoiseArray=randn(size(x));
PinkNoiseArray=pinknoise(length(x));
PropNoiseArray=propnoise(y);
% Add noise to signal
% noise may be white, pink, blue, proportional, or sqrt.
y=y+noise.*PropNoiseArray; % <<< Change noise type here (WhiteNoiseArray,PinkNoiseArray, or PropNoiseArray

% Input arguments needed for peak detection
SlopeThreshold=0.01*(width./increment)^-2; % <<< make this smaller to detect more peaks
AmpThreshold=.4; %  <<< Adjust to less than the smallest desired peak.
smoothwidth=round(.5*width./increment);
fitwidth=round(.5*width./increment);
smoothtype=3;
% Additional input arguments needed for findpeaksb function
window=3*width./increment;

% Peak detection and measurement
clf
% findpeaks function in Matlab's Signal Processing Toolkit
% If you have the Signal Processing Toolkit, you can uncomment the next
% line
% [PeakHeights,PeakPositions]=findpeaks(y,'MINPEAKHEIGHT',AmpThreshold,'MINPEAKDISTANCE',window)

% findpeaksG and findpeaksL functions
if peakshape==2,
    tic
    Pg=findpeaksL(x,y,SlopeThreshold,AmpThreshold,smoothwidth,fitwidth,smoothtype);
    findpeaksGtime=toc;
else
    tic
    Pg=findpeaksG(x,y,SlopeThreshold,AmpThreshold,smoothwidth,fitwidth,smoothtype);
    findpeaksGtime=toc;
end
PeaksFound=length(Pg(:,1));
plot(x,y,'.')
% findpeaksb function
tic
Pb=findpeaksb(x,y,SlopeThreshold,AmpThreshold,smoothwidth,fitwidth,smoothtype,window,peakshape,extra,NumTrials,AUTOZERO);
findpeaksbtime=toc;

% findpeaksb3 function
tic
Pb3=findpeaksb3(x,y,SlopeThreshold,AmpThreshold,smoothwidth,fitwidth,smoothtype,peakshape,extra,NumTrials,AUTOZERO,0);
findpeaksb3time=toc;

% findpeaksfit function
tic
[P,Pf]=findpeaksfit(x,y,SlopeThreshold,AmpThreshold,smoothwidth,fitwidth,smoothtype,peakshape,extra,NumTrials,AUTOZERO,0,0);
findpeaksfittime=toc;

% Calculation and display of peak parameter errors
disp(' ')
disp('Average absolute percent errors of all peaks')
disp('           Position error  Height error  Width error   Elapsed time, sec')
for peak=1:PeaksFound,
    Perror(peak)=-100.*(Positions(peak)-Pg(peak,2))./Positions(peak);
    Herror(peak)=-100.*(Heights(peak)-Pg(peak,3))./Heights(peak);
    Werror(peak)=-100.*(Widths(peak)-Pg(peak,4))./Widths(peak);
    % disp([Positions(peak) Perror(peak) Herror(peak) Werror(peak)])
    findpeaksGErrors=[Perror(:) Herror(:) Werror(:)];
    figure(2);subplot(4,1,1);bar(findpeaksGErrors)
    title('Blue=position error     Green=Height error     Red=Width error')
end
disp(['findpeaksG     ' num2str(mean(abs(findpeaksGErrors))) '      ' num2str(findpeaksGtime) ])
xlabel('findpeaksG/L Errors')

for peak=1:PeaksFound,
    Perror(peak)=-100.*(Positions(peak)-Pb(peak,2))./Positions(peak);
    Herror(peak)=-100.*(Heights(peak)-Pb(peak,3))./Heights(peak);
    Werror(peak)=-100.*(Widths(peak)-Pb(peak,4))./Widths(peak);
    % disp([Positions(peak) Perror(peak) Herror(peak) Werror(peak)])
    findpeaksbErrors=[Perror(:) Herror(:) Werror(:)];
    figure(2);subplot(4,1,2);bar(findpeaksbErrors)
end
disp(['findpeaksb     ' num2str(mean(abs(findpeaksbErrors))) '      ' num2str(findpeaksbtime)])
xlabel('findpeaksb Errors')


for peak=1:PeaksFound,
    Perror(peak)=-100.*(Positions(peak)-Pb3(peak,2))./Positions(peak);
    Herror(peak)=-100.*(Heights(peak)-Pb3(peak,3))./Heights(peak);
    Werror(peak)=-100.*(Widths(peak)-Pb3(peak,4))./Widths(peak);
    % disp([Positions(peak) Perror(peak) Herror(peak) Werror(peak)])
    findpeaksb3Errors=[Perror(:) Herror(:) Werror(:)];
    figure(2);subplot(4,1,3);bar(findpeaksb3Errors)
end
disp(['findpeaksb3     ' num2str(mean(abs(findpeaksb3Errors))) '      ' num2str(findpeaksb3time)])
xlabel('findpeaksb3 Errors')

for peak=1:PeaksFound,
    Perror(peak)=-100.*(Positions(peak)-Pf(peak,2))./Positions(peak);
    Herror(peak)=-100.*(Heights(peak)-Pf(peak,3))./Heights(peak);
    Werror(peak)=-100.*(Widths(peak)-Pf(peak,4))./Widths(peak);
    % disp([Positions(peak) Perror(peak) Herror(peak) Werror(peak)])
    findpeaksfitErrors=[Perror(:) Herror(:) Werror(:)];
    figure(2);subplot(4,1,4);bar(findpeaksfitErrors)
end
disp(['findpeaksfit   ' num2str(mean(abs(findpeaksfitErrors))) '      ' num2str(findpeaksfittime)])
xlabel('findpeaksfit Errors')

figure(1)
for peak=1:PeaksFound,
    hold on
    xrange=(Pb(peak,2)-1.5*Widths(peak)):increment:(Pb(peak,2)+1.5*Widths(peak));
    TopPart=(P(peak,2)-fitwidth/2):increment:(P(peak,2)+fitwidth/2);
    if peakshape==2,
    plot(TopPart,Pg(peak,3).*lorentzian(TopPart,Pg(peak,2),Pg(peak,4)),'m');    
    plot(x,Heights(peak).*lorentzian(x,Positions(peak),Widths(peak)),'k.');     
    plot(xrange,Pb(peak,3).*lorentzian(xrange,Pb(peak,2),Pb(peak,4)),'r');
    plot(x,Pb3(peak,3).*lorentzian(x,Pb3(peak,2),Pb3(peak,4)),'c');
    plot(x,Pf(peak,3).*lorentzian(x,Pf(peak,2),Pf(peak,4)),'g');
    else
    plot(TopPart,Pg(peak,3).*gaussian(TopPart,Pg(peak,2),Pg(peak,4)),'m'); 
    plot(x,Heights(peak).*gaussian(x,Positions(peak),Widths(peak)),'k.'); 
    plot(xrange,Pb(peak,3).*gaussian(xrange,Pb(peak,2),Pb(peak,4)),'r');
    plot(x,Pb3(peak,3).*gaussian(x,Pb3(peak,2),Pb3(peak,4)),'c');
    plot(x,Pf(peak,3).*gaussian(x,Pf(peak,2),Pf(peak,4)),'g');
    end
    hold off
end
 xlabel(['    AmpThreshold: ' num2str(AmpThreshold) '     SlopeThreshold: ' num2str(SlopeThreshold) '    SmoothWidth: ' num2str(smoothwidth) '    FitWidth: ' num2str(fitwidth)  '    Peak shape: ' num2str(peakshape)  '    AUTOZERO: ' num2str(AUTOZERO)])
 title(['Black:True peaks   Blue:Raw signal    Magenta:findpeaksG/L   Red:findpeaksb    Cyan:findpeaksb3    Green:findpeaksfit.  '   ])
 
 figure(3);bar([mean(abs(findpeaksGErrors));mean(abs(findpeaksbErrors));mean(abs(findpeaksb3Errors));mean(abs(findpeaksfitErrors))])
 xlabel('findpeaksG             findpeaksb         findpeaksb3        findpeaksfit')
 if PeaksFound==1;
     disp('1 peak found. ')
 else
     disp([ num2str(PeaksFound ) ' peaks found. '])
 end
 ylabel('Average absolute percent error of all peaks')
 title('Blue=position error     Green=Height error     Red=Width error')
 % Also, you could compare this to the use of the peakfit.m function, 
 % which requires that you know the peak shape and number of peaks.
 % peakfit([x;y],0,0,NumPeaks,peakshape,extra,NumTrials,0,1). That would
 % be better for highly overlapped peaks that do not exhibit distinct
 % maxima for each peak. Try adding an additional peak to the model to
 % account for the baseline.
 