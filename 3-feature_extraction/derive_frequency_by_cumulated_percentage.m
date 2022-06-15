function contour_f_points= derive_frequency_by_cumulated_percentage(freq_axis, spectrumbins, threshold)
    freqbin_len = size(spectrumbins, 1);
    timebin_len = size(spectrumbins, 2);
    
    contour_f_points = []; % 轮廓的freq points集合
    
    for i=1:timebin_len
        power_sum = sum(spectrumbins(:, i));
        for k=1:freqbin_len
            cumulated_sum = sum(spectrumbins(1:k, i));
            if cumulated_sum / power_sum > threshold
                contour_f_points = [contour_f_points, freq_axis(k)];
                break;
            end
            % 提取符合条件的频率最大值
        end
    end
end

