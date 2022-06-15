function [R] = acf(detect_period)

    % calculate autocorrelation
    remove_average = detect_period - mean(detect_period);
    time_len = length(remove_average);
    % element-wise multiplication
    R = zeros(1, time_len);
    for pau=1:time_len
        R(1, pau) = sum(remove_average(1:time_len-pau+1) .* remove_average(pau:time_len));
    end
%     plot(R)
   
    % [~, locs] = min(R); %     % [~, locs] = min(R); % 其实我觉得按论文里，这个应该是locs后面那个再上去的峰值
    % 按照论文中，半周期应该是第二个峰值点
%     [~, locs] = findpeaks(R);
end

