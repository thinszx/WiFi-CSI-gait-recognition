function data = bandpass_filter(data, upper_order,lower_order,maxf, minf, hsr)
% CSI_DATA=BANDPASS_FILTER(CSI_DATA,UPORDER,lOWORDER,MAXF,MINF,HSR) bandpass filter the
% signal.
%
% CSI_DATA : Power of raw CSI measurements
% MAXF     : Upper cutoff frequency of the filter
% MINF     : Lower cutoff frequency of the filter
% HSR      : Half sample rate
%
% CSI_DATA : Filtered result

[b,a] = butter(upper_order,maxf/hsr,'low');
[d,c] = butter(lower_order,minf/hsr,'high');
for jj = 1:size(data,2)
    data(:,jj) = filter(b,a,data(:,jj));
    data(:,jj) = filter(d,c,data(:,jj));
end
end