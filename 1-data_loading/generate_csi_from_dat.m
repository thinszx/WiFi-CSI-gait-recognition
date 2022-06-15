function [cfr_array, timestamp] = generate_csi_from_dat(filename,n_antennas,n_subcarriers)
% [csi_data, ~] = generate_csi([spfx_ges, '-r', num2str(1), '.dat'],n_antennas,n_subcarriers);
% load data from *.dat
%
% filename     : file path
%
% cfr_array          : Frequency time profile of CSI data.
% timestamp            : Frequency axis
%
csi_trace = read_bf_file(filename);
timestamp = zeros(length(csi_trace), 1);
cfr_array = zeros(length(csi_trace), n_antennas*n_subcarriers);

for k = 1:length(csi_trace)
    csi_entry = csi_trace{k}; % for the k_{th} packet
    % test
    test = get_scaled_csi(csi_entry);
    %
    csi_all = squeeze(get_scaled_csi(csi_entry)).'; 
    csi = [csi_all(:,1); csi_all(:,2); csi_all(:, 3)].'; 
    timestamp(k) = csi_entry.timestamp_low;
    cfr_array(k,:) = csi;
end