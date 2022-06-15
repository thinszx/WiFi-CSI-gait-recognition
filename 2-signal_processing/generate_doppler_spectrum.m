function [doppler_spectrum, freq_axis_bin,velocity_axis_bin,loading_data_time,signal_processing_time] = generate_doppler_spectrum(data_file, n_receivers, n_antennas,n_subcarriers,...
    method,wave_length,sample_rate,n_pca)

% Set Parameters
window_size = 128;
upper_order = 6;
upper_stop = 100;
lower_order = 3;
lower_stop = 2;

%% 如果要改上下界，改这个地方的参数
freq_bins_unwrap = [0:sample_rate/2-1 -sample_rate/2:-1]'/sample_rate;%[0,...,0.5,-0.5,...,-1/1000]
% freq_lpf_sele = freq_bins_unwrap <= upper_stop / sample_rate & freq_bins_unwrap >= -upper_stop / sample_rate;
freq_lpf_sele = freq_bins_unwrap <= upper_stop / sample_rate & freq_bins_unwrap >= 0;
freq_lpf_positive_max = sum(freq_lpf_sele(2:length(freq_lpf_sele)/2));
% freq_lpf_negative_min = sum(freq_lpf_sele(length(freq_lpf_sele)/2:end));
freq_lpf_negative_min = 0;

% Frequency Bin(Corresponding to FFT Results)
freq_axis_bin = [0:freq_lpf_positive_max -1*freq_lpf_negative_min:-1];
velocity_axis_bin = freq_axis_bin * wave_length/2; % Speed axis.

% Doppler Spectrum For Each Antenna
% [csi_data, ~] = generate_csi_from_dat([data_file, '-r', num2str(1), '.dat'],n_antennas,n_subcarriers);
[csi_data, ~] = generate_csi_from_dat(data_file,n_antennas,n_subcarriers);
test = floor(size(csi_data, 1));
doppler_spectrum = zeros(n_receivers, n_pca,1+freq_lpf_positive_max + freq_lpf_negative_min,...
    floor(size(csi_data, 1)));

loading_data_time=0;
signal_processing_time=0;
% for ii = 1:n_receivers
% %     spth = [data_file, '-r', num2str(ii), '.dat'];
%     spth = data_file;
%     try
%         tic;
%         [csi_data, ~] = generate_csi_from_dat(spth,n_antennas,n_subcarriers);
%     catch err
%         disp(err)
%         continue
%     end
%     loading_data_time=loading_data_time+toc;
%     tic;
%     csi_data = csi_data(round(1:1:size(csi_data,1)),:);
%     conj_mult = conj_denoise(csi_data,n_antennas);
%     signal_denoised = bandpass_filter(conj_mult,upper_order,lower_order, upper_stop, lower_stop, sample_rate / 2);
%     
%     %     % PCA analysis.
%     pca_coef = pca(signal_denoised);
%     signal_pca = signal_denoised * pca_coef(:,1:n_pca);
%     % TFA With CWT or STFT
%     for pca_index=1:n_pca
%         if strcmp(method, 'cwt')
%             freq_time_prof_allfreq = scaled_cwt(signal_pca(:,pca_index),frq2scal(...
%                 [1:sample_rate/2 -1*sample_rate/2:-1],'cmor4-1',1/sample_rate),'cmor4-1');
%         else if strcmp(method, 'stft')
%                 time_instance = 1:length(signal_pca(:,pca_index));
%                 if(~mod(window_size,2))
%                     window_size = window_size + 1;
%                 end
%                 freq_time_prof_allfreq = tfrsp(signal_pca(:,pca_index), time_instance,...
%                     sample_rate, tftb_window(window_size, 'gauss'));
%             end
%         end
%         % Select Concerned Freq
%         freq_time_prof = freq_time_prof_allfreq(freq_lpf_sele, :);
%         % Spectrum Normalization By Sum For Each Snapshot
%         freq_time_prof = abs(freq_time_prof) ./ repmat(sum(abs(freq_time_prof),1),...
%             size(freq_time_prof,1), 1);
%         % Store Doppler Velocity Spectrum
%         if(size(freq_time_prof,2) >= size(doppler_spectrum,4))
%             doppler_spectrum(ii,pca_index,:,:) = freq_time_prof(:,1:size(doppler_spectrum,4));
%         else
%             doppler_spectrum(ii,pca_index,:,:) = [freq_time_prof zeros(size(doppler_spectrum,3),...
%                 size(doppler_spectrum,4) - size(freq_time_prof,2))];
%         end
%     end
%     signal_processing_time=signal_processing_time+toc;
% end

%     spth = [data_file, '-r', num2str(ii), '.dat'];
spth = data_file;
try
   tic;
   [csi_data, ~] = generate_csi_from_dat(spth,n_antennas,n_subcarriers);
catch err
   disp(err)
end
loading_data_time=loading_data_time+toc;
tic;
csi_data = csi_data(round(1:1:size(csi_data,1)),:);
conj_mult = conj_denoise(csi_data,n_antennas);
signal_denoised = bandpass_filter(conj_mult,upper_order,lower_order, upper_stop, lower_stop, sample_rate / 2);
    
% PCA analysis.
pca_coef = pca(signal_denoised);
signal_pca = signal_denoised * pca_coef(:,1:n_pca);
% TFA With CWT or STFT
for pca_index=1:n_pca
    if strcmp(method, 'cwt')
        freq_time_prof_allfreq = scaled_cwt(signal_pca(:,pca_index),frq2scal(...
            [1:sample_rate/2 -1*sample_rate/2:-1],'cmor4-1',1/sample_rate),'cmor4-1');
    else if strcmp(method, 'stft')
            time_instance = 1:length(signal_pca(:,pca_index));
            if(~mod(window_size,2))
                window_size = window_size + 1;
            end
            freq_time_prof_allfreq = tfrsp(signal_pca(:,pca_index), time_instance,...
                1024, tftb_window(window_size, 'gauss'));
        end
    end
    % Select Concerned Freq
    freq_time_prof = freq_time_prof_allfreq(freq_lpf_sele, :);
    % Spectrum Normalization By Sum For Each Snapshot
    freq_time_prof = abs(freq_time_prof) ./ repmat(sum(abs(freq_time_prof),1),...
        size(freq_time_prof,1), 1);
%     % compute noise floor
%     noise_floor = 
    % Store Doppler Velocity Spectrum
    if(size(freq_time_prof,2) >= size(doppler_spectrum,4))
        doppler_spectrum(1,pca_index,:,:) = freq_time_prof(:,1:size(doppler_spectrum,4));
    else
        doppler_spectrum(1,pca_index,:,:) = [freq_time_prof zeros(size(doppler_spectrum,3),...
            size(doppler_spectrum,4) - size(freq_time_prof,2))];
    end
end
signal_processing_time=signal_processing_time+toc;

loading_data_time=loading_data_time/n_receivers;
signal_processing_time=signal_processing_time/n_receivers;
end