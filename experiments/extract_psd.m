% function psd = extract_psd(csi_signal)
    data_file = "D:\Documents\14self\wireless-recognition\gait-WiDar\Gait_Dataset\CSI_Gait\user1-20190627\user1-1-3-r3.dat";
    [csi_data, ~] = generate_csi_from_dat(data_file, 3, 30);
    csi_data = csi_data(round(1:1:size(csi_data,1)),:);
    conj_mult = conj_denoise(csi_data,3);
    
    % Set Parameters
    window_size = 128;
    upper_order = 6;
    upper_stop = 100;
    lower_order = 3;
    lower_stop = 2;
    sample_rate = 1000;
    n_pca = 20;
    method = 'stft';
    
    freq_bins_unwrap = [0:sample_rate/2-1 -sample_rate/2:-1]'/sample_rate;%[0,...,0.5,-0.5,...,-1/1000]
    freq_lpf_sele = freq_bins_unwrap <= upper_stop / sample_rate & freq_bins_unwrap >= -upper_stop / sample_rate;
    % freq_lpf_sele = freq_bins_unwrap <= upper_stop / sample_rate & freq_bins_unwrap >= 0;
    freq_lpf_positive_max = sum(freq_lpf_sele(2:length(freq_lpf_sele)/2));
    freq_lpf_negative_min = sum(freq_lpf_sele(length(freq_lpf_sele)/2:end));
    
    doppler_spectrum = zeros(3, n_pca,1+freq_lpf_positive_max + freq_lpf_negative_min,...
        floor(size(csi_data, 1)));
    
    signal_denoised = bandpass_filter(conj_mult,upper_order,lower_order, upper_stop, lower_stop, sample_rate / 2);

    % PCA analysis.
    pca_coef = pca(signal_denoised);
    signal_pca = signal_denoised * pca_coef(:,1:n_pca);
%     segment = signal_pca(1:256, 1);
    signal = signal_pca(:, 1);
    stepsize = 100; % 每次滑动窗口的前进步长
    windowsize = 256; % 滑动窗口长度
    noverlap = windowsize - stepsize; % overlap长度
    
    slice_cnt = fix((length(signal)-noverlap)/stepsize); % 应有的切片个数
    psd_sum = zeros(slice_cnt, windowsize);
    kurto_test = zeros(slice_cnt, 1);
    for iSlice=1:slice_cnt
        signal_slice = signal(1+stepsize*(iSlice-1):windowsize+stepsize*(iSlice-1),:);
%         psd_sum(iSlice, :) = abs(fftshift(fft(signal_slice))).^2;
        psd_sum(iSlice, :) = periodogram(signal_slice, 'centered', 'power');
        disp(kurtosis(abs(signal_slice)));
    end
    plot(abs(signal));
%     tttt = kurtosis(abs(signal));
    
    psd_sum = psd_sum.';
    Fs = 1000;
    index=0:windowsize/2-1;
%     fx=index*Fs/(windowsize-1);
%     psd = 10 * log(psd_sum ./ windowsize);
    
%     fx = [0:windowsize-1]*(1/(windowsize*(1/Fs)));    
%     plot(f, psd_sum);
%     x = signal;
kurto = kurtosis(psd_sum);
% end

