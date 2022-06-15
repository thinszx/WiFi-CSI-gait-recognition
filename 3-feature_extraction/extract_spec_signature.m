function magnitude_array = extract_spec_signature(one_period_spec_bins)
    magnitude_array = zeros(1, size(one_period_spec_bins, 1));
    for freq_num=1:size(one_period_spec_bins, 1) % 遍历所有的频率，横向加和
        magnitude_array(freq_num) = mean(one_period_spec_bins(freq_num, :));
    end
end

