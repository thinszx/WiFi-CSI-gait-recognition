function conj_mult = conj_denoise(csi_data,rx_acnt)
% conj_mult=conj_denoise(csi_data,rx_acnt) filter the
% signal with conj.
% 
% csi_data : raw CSI measurements, eg:90*1479
% rx_acnt  : Antenna count for each receiver, eg:3
% 
% conj_mult : Filtered result


csi_data = csi_data(round(1:1:size(csi_data,1)),:);

csi_mean = mean(abs(csi_data));
csi_var = sqrt(var(abs(csi_data)));
csi_mean_var_ratio = csi_mean./csi_var;
[~,idx] = max(mean(reshape(csi_mean_var_ratio,[30 rx_acnt]),1));
csi_data_ref = repmat(csi_data(:,(idx-1)*30+1:idx*30), 1, rx_acnt);

csi_data_adj = zeros(size(csi_data));
csi_data_ref_adj = zeros(size(csi_data_ref));
alpha_sum = 0;
for jj = 1:30*rx_acnt
    amp = abs(csi_data(:,jj));
    alpha = min(amp(amp~=0));
    alpha_sum = alpha_sum + alpha;
    csi_data_adj(:,jj) = abs(abs(csi_data(:,jj))-alpha).*exp(1j*angle(csi_data(:,jj)));
end
beta = 1000*alpha_sum/(30*rx_acnt);
for jj = 1:30*rx_acnt
    csi_data_ref_adj(:,jj) = (abs(csi_data_ref(:,jj))+beta).*exp(1j*angle(csi_data_ref(:,jj)));
end
conj_mult = csi_data_adj .* conj(csi_data_ref_adj);
end