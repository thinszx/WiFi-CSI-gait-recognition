function fc_contour = derive_torso_contour_frequency(freq_axis, spectrumbins, threshold)
    freqbin_len = size(spectrumbins, 1);
    timebin_len = size(spectrumbins, 2);

    fc_contour = [];

    for i=1:timebin_len
        contour_f_points = []; % 满足要求的一个时间点的freq points集合
        power_sum = sum(spectrumbins(:, i));
        for k=1:freqbin_len
            power_point = spectrumbins(k, i);
            if power_point / power_sum > threshold
                contour_f_points = [contour_f_points, freq_axis(k)];
            end
            % 提取符合条件的频率最大值
        end
        if size(contour_f_points) == 0
            fc_contour = [fc_contour, 0];
        else
            fc_contour = [fc_contour, max(contour_f_points)];
        end
    end
end

