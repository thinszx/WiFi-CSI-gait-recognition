  clf
x=1:.2:100;
  y=((2-(x./50))/2)+modelpeaks(x,3,5,[1 2 3],[20 50 80],[3 3 3],10)+.05*randn(size(x));
  disp('          Peak      Position      Height      Width        Area')
  plot(x,y);NoBackgroundSubtraction=findpeaksplot(x,y,.00005,.5,30,20,3)
  disp(' ')
  disp('          Peak      Position      Height      Width        Area       % error          R2')
  LinearBackgroundSubtraction=findpeaksb(x,y,.00005,.5,30,20,3,150,5,10,1,1)