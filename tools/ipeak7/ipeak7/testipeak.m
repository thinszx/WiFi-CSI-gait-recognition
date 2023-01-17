% Script that tests for the proper installation and operation of iPeak by
% running quickly through all the examples and demos for iPeak on
% http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm#ipeak
% Assumes that ipeakdata.mat has been loaded into the workspace. 
% EXAMPLE 1: 
y=cos(.1:.1:100);ipeak(y);
% EXAMPLE 2: 
x=[0:.01:5]';y=x.*sin(x.^2).^2;M=[x y];ipeak(M);
% EXAMPLE 3: 
x=[0:.1:100];y=(x.*sin(x)).^2;ipeak(x,y);
% EXAMPLE 4: 
x=[0:.1:100];y=5+5.*cos(x)+randn(size(x));ipeak(x,y,10);
ipeak([x;y],10);
ipeak(humps(0:.01:2),3);
x=[0:.1:10];y=exp(-(x-5).^2);ipeak([x' y'],1);
% EXAMPLE 5: 
x=[0:.01:5]';y=x.*sin(x.^2).^2;M=[x y];ipeak(M,0,0,.0001,20,20);
% EXAMPLE 6: 
load ipeakdata.mat;ipeak(Sample1,0,110,0.06,3,4,249.7,0.4);
% EXAMPLE 7: 
ipeak(Sample1,0,110,0.06,3,4,249.7,0.4,1);
% EXAMPLE 8: 
ipeak(Sample1,0,100,0.05,3,6,296,5,0.1,Positions,Names);
% iPeak demos
ipeakdemo
ipeakdemo1
ipeakdemo2
ipeakdemo3
ipeakdemo4
ipeakdemo5
