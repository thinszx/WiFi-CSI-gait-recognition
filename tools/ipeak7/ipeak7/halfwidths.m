function FWHM=halfwidths(x,y,xo)
% function FWHM=halfwidths(x,y) computes the approximate full width at half
% maximum of any shape peak that peaks at xo, falls to at least one-half
% intensity on both sides, and has a zero baseline. If xo is a vector w=
% of peak positions of multiple peaks in the y vector, it will return a
% vector of FWHM values for each peak. Not highly accurate if the data are
% is too noisy. Tom O'Haver (toh@umd.edu) April 2016
%
% Example 1:
% x=-5:.1:5;
% y=sinc(x);
% plot(x,y)
% FWHM=halfwidths(x,y,0)
%
% Example 2:
% x=1:100;
% pos1=30;
% pos2=70;
% w1=20;
% w2=25;
% y=gaussian(x,pos1,w1)+gaussian(x,pos2,w2);
% ny=y+.01.*randn(size(x));
% plot(x,ny)
% FWHM=halfwidths(x,ny,[pos1 pos2])
%
FWHM=zeros(size(xo));
for p=1:length(xo),
    try
        indmax=val2ind(x,xo(p));
        maxy=y(indmax);
        oy=y-maxy/2;
        
        n=indmax;
        while oy(n)>0,
            n=n-1;
        end
        x1=interp1([y(n) y(n+1)],[x(n) x(n+1)],maxy/2);
        
        n=indmax;
        while oy(n)>0,
            n=n+1;
        end   
        x2= interp1([y(n-1) y(n)],[x(n-1) x(n)],maxy/2);
        
        FWHM(p)=x2-x1;
        % [xo indmax maxy n1 x1 n2 x2 FWHM]
    catch
        FWHM(p)=NaN;
    end
end