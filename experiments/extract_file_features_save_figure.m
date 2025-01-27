function features = extract_file_features_save_figure(data_file, save_dir, filename)
% data_file = "D:\Documents\14self\wireless-recognition\gait-WiDar\Gait_Dataset\CSI_Gait\user1-20190719\user1-1-6-r3.dat";
    %% generate spectrogram
    wave_length = 299792458 / 5.825e9;
    [doppler_spectrum, freq_axis_bin,velocity_axis_bin,~,~] = generate_doppler_spectrum(data_file, 6, 3, 30, 'stft',wave_length,1000,20);

    [~,idx] = max(freq_axis_bin);
    circle_length = length(freq_axis_bin) - idx;
%     doppler_spectrum_noflip = circshift(doppler_spectrum, [0 0 circle_length 0]);
    doppler_spectrum = circshift(doppler_spectrum, [0 0 circle_length 0]);
    freq_axis_bin=circshift(freq_axis_bin, [0 circle_length]);
    velocity_axis_bin=circshift(velocity_axis_bin, [0 circle_length]);
    
    %% flip negative values to positive
%     doppler_spectrum = zeros(size(doppler_spectrum, 1), size(doppler_spectrum, 2), ...
%                              (size(doppler_spectrum, 3)-1)/2+1, ... % 除去0以外频率的二分之一+频率0
%                              size(doppler_spectrum, 4));
%     zero_freq_bin_i = median(1:length(freq_axis_bin)); % 0频率在freq_axis_bin中的index
%     for c=1:20 % pca compnt
%         for i=1:zero_freq_bin_i
%             doppler_spectrum(1,c,i,:) = (doppler_spectrum_noflip(1,c,zero_freq_bin_i+i-1,:) ...
%                                        + doppler_spectrum_noflip(1,c,zero_freq_bin_i-i+1,:)) / 2;
%         end
%     end

    upper_stop = 100;
%     colormap(jet);
    size_spectrum = size(doppler_spectrum);
    doppler_sum = zeros(size(doppler_spectrum));

    for i=1:1 % pca compnt
        % filter
        mdlSpecLpf = fspecial('gaussian', 5, 0.8); % Build filter
        outSpec = imfilter(doppler_spectrum(1,i,:,:),mdlSpecLpf); % Smoothing
        doppler_sum(1,1,:,:) = doppler_sum(1,1,:,:) + outSpec;
    end
    squeezed = squeeze(doppler_sum(1,1,:,:));
%     mesh(1:size(doppler_spectrum,4),0:upper_stop,squeezed);
%     view([0,90]);

    %% extract torso speeds
    spectrumbins = squeeze(doppler_sum(1,1,:,:)); % size: 101*2678
    % spectrumbins_trans = transpose(spectrumbins); % 将

    freq_axis = 0:upper_stop;

    fc_contour = derive_torso_contour_frequency(freq_axis, spectrumbins, 0.03);
%     % plot(fc_contour);

    % lowpass(fc_contour,150,1000)
    fc = 300;
    fs = 1000;
    [b,a] = butter(10, fc/(fs/2));
    fc_filter = filter(b,a,fc_contour);
    fc_contour_freq_smooth = smoothdata(fc_filter,'gaussian',250); % 250大概是500ms的一半，也就是近似步态周期的一半
%     % plot(fc_contour_freq_smooth);
    movement_speed = freq2speed(fc_contour_freq_smooth, wave_length);
%     % plot(movement_speed);

    torso_freq = derive_frequency_by_cumulated_percentage(freq_axis, spectrumbins, 0.5); % 提取躯干频率
%     % plot(torso_freq);
    torso_filter = filter(b,a,torso_freq);
    torso_freq_smooth = smoothdata(torso_filter,'gaussian',250); % 250大概是500ms的一半，也就是近似步态周期的一半
    torso_speed = freq2speed(torso_freq_smooth, wave_length);
%     % plot(torso_speed);

    leg_freq = derive_frequency_by_cumulated_percentage(freq_axis, spectrumbins, 0.95); % 提取躯干频率
    leg_filter = filter(b,a,leg_freq);
    leg_freq_smooth = smoothdata(leg_filter,'gaussian',250); % 250大概是500ms的一半，也就是近似步态周期的一半
    leg_speed = freq2speed(leg_freq_smooth, wave_length);
    
%     % plot(leg_speed);
    x = squeeze(1:size(doppler_spectrum, 4));
%     plot(x,movement_speed, x,torso_speed, x,leg_speed);
%     legend('movement speed', 'torso speed', 'leg speed');

    %% autocorrelation
    % extract period where speed is no less than 80% maxmium speed
    max_speed_80 = max(movement_speed) * 0.8;
    max_speed_75 = max(movement_speed) * 0.75;
    left_index = 0;
    right_index = 0;
    for i=1:length(movement_speed)
       if movement_speed(1, i) > max_speed_80
           left_index = i;
           break;
       end
    end
    for i=length(movement_speed):-1:1
        if movement_speed(1, i) > max_speed_75
           right_index = i;
           break;
        end
    end
    detect_period = movement_speed(left_index:right_index);
%     period = autocorrelation(detect_period);
    remove_average = detect_period - mean(detect_period);
    [period, locs, diffs] = autocorrelation(remove_average);
    
    average_speed = mean(movement_speed);
    steady_average_speed = mean(detect_period);
    
    %% save axes figure
    gca = axes;
    clf(gca, 'reset');
    hold(gca, 'on');
    imagesc(squeezed);
    set(gca,'YDir','normal') 
    yyaxis(gca, 'left');
    gca.YLim = [0, upper_stop];
    gca.XLim = [1, length(squeezed)];
    
    yyaxis(gca, 'right');
    gca.YLim = [0, upper_stop*wave_length/2];
    xline(gca, [left_index, right_index], 'red', 'LineWidth', 2);
    plot(gca, x,movement_speed, x,torso_speed, x,leg_speed, 'LineWidth', 2);
    legend(gca, 'movement speed', 'torso speed', 'leg speed');
    
    locs = locs + left_index;
    for i=1:length(locs)
        if locs(i) <= length(squeezed)
            plot(locs(i), movement_speed(locs(i)), '.', 'MarkerSize', 15, 'Color', 'red');
        end
    end
    if length(diffs) == 1
        difference = 1;
    else
        difference = max(diffs)-min(diffs);
    end
    save_path = fullfile(save_dir, append(...
                        num2str(difference), "_", ...
                        num2str(round(period)), "_", ...
                        num2str(round(average_speed*1000)), "_", ...
                        num2str(round(steady_average_speed*1000)), "_", filename));
    saveas(gca, save_path, 'bmp');
    hold(gca, 'off');
    clf(gca, 'reset');
    

    %% extract spectrogram signature
    % extract walking period
%     plot(detect_period);
    % 获取detect_period平稳周期的第一个峰值和最后一个，分别作为第一个步态开始和最后一个步态结束的半周期点，然后依据这两个锚点确定步态特征提取的起点和终点
%     startpoint = round(locs(1) - period/2);
%     endpoint = round(locs(end) + period/2);
%     % 防止越界
%     if startpoint < 1
%         startpoint = 1;
%     end
%     if endpoint > length(movement_speed)
%         endpoint = length(movement_speed);
%     end
    % 遍历每个周期，进行步态频谱特征的提取
%     % 计算周期数，用作循环遍历
%     period_cnt = ceil((endpoint-startpoint)/period);
    % 存储四个周期的均值
    fisrt_magnitude_array = [];
    second_magnitude_array = [];
    third_magnitude_array = [];
    fourth_magnitude_array = [];
    period_start = round(locs(1) - period/2);
    for i=1:length(locs)
%         period_start = round(locs(i) - period/2);
        if period_start < 1
            period_start = 1;
        end
        period_end = locs(i);
        slice_len = period_end - period_start;
        one_period = spectrumbins(:, period_start:period_end);
        fisrt_stage = one_period(:, 1:round(slice_len/4)); % 半周期的1/4
        second_stage = one_period(:, round(slice_len/4):round(slice_len/2));
        third_stage = one_period(:, round(slice_len/2):round(3*slice_len/4));
        fourth_stage = one_period(:, round(3*slice_len/4):end);

        fisrt_magnitude_array = [fisrt_magnitude_array; extract_spec_signature(fisrt_stage)];
        second_magnitude_array = [second_magnitude_array; extract_spec_signature(second_stage)];
        third_magnitude_array = [third_magnitude_array; extract_spec_signature(third_stage)];
        fourth_magnitude_array = [fourth_magnitude_array; extract_spec_signature(fourth_stage)];

        % 更新周期
        % 到达边界，停止检测
        if i < length(locs) && locs(i+1) < length(squeezed)
            period_between = movement_speed(locs(i):locs(i+1));
            valley_loc = findvalleys([1:length(period_between)], period_between, 0, min(period_between)-1, 1, 1, 1);
            if size(valley_loc, 1) > 1
                period_start = round(locs(i) + period / 2);
            else
                period_start = locs(i) + valley_loc(1, 2);
            end
        else
            break;
        end
%         period_start = round(period_start + period/2);
%         period_end = round(period_end + period/2);
%         if (period_end > endpoint) && (period_end < length(movement_speed))
%             break;
%         end
    end
    %% extract other features
    fisrt_magnitude_average = sum(fisrt_magnitude_array, 1) / size(fisrt_magnitude_array, 1); % 沿第一维度和的平均值
    second_magnitude_average = sum(second_magnitude_array, 1) / size(fisrt_magnitude_array, 1);
    third_magnitude_average = sum(third_magnitude_array, 1) / size(fisrt_magnitude_array, 1);
    fourth_magnitude_average = sum(fourth_magnitude_array, 1) / size(fisrt_magnitude_array, 1);
    gait_cycle = period / 1000; % 步态周期
    estimated_footstep = average_speed * period / 1000; % 模拟步长
    estimated_steady_footstep = steady_average_speed * period / 1000; % 模拟步长
    % slice the steady period of leg speed and torso dpeed accordingt to previous index of steady movement speed
    % extract leg speed feature
    leg_steady_period = leg_speed(left_index:right_index);
    max_leg_speed = max(leg_steady_period);
    min_leg_speed = min(leg_steady_period);
    average_leg_speed = mean(leg_steady_period);
    var_leg_speed = var(leg_steady_period, 1);
    % extract torso speed feature
    torso_steady_speed = torso_speed(left_index:right_index);
    max_turso_speed = max(torso_steady_speed);
    min_turso_speed = min(torso_steady_speed);
    average_turso_speed = mean(torso_steady_speed);
    var_turso_speed = var(torso_steady_speed, 1);

    %% save features
    % 1. gait_cycle
    % 2. estimated_footstep
    % 3. max_leg_speed
    % 4. min_leg_speed
    % 5. average_leg_speed
    % 6. var_leg_speed
    % 7. max_turso_speed
    % 8. min_turso_speed
    % 9. average_turso_speed
    % 10. var_turso_speed
    % 11-111  fisrt_magnitude  - length 101
    % 112-212 second_magnitude - length 101
    % 213-313 third_magnitude_average - length 101
    % 314-414 fourth_magnitude_average - length 101
    features = [gait_cycle, estimated_footstep, estimated_steady_footstep, ...
                max_leg_speed, min_leg_speed, average_leg_speed, var_leg_speed, ...
                max_turso_speed, min_turso_speed, average_turso_speed, var_turso_speed, ...
                fisrt_magnitude_average, second_magnitude_average, third_magnitude_average, fourth_magnitude_average];
end

