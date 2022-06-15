%% parse filename
parent_dir = "D:\Documents\14self\wireless-recognition\gait-WiDar\Gait_Dataset\CSI_Gait\";
user_dir = 'user2-20190627';
data_dir = fullfile(parent_dir, user_dir);
data_array  = dir(append(data_dir, '\*.dat'));
test = 1:415;
% test = extract_file_features(append(data_dir, data_array(3).name)); % 先提取一次特征得到特征长度，列表前两个是.和..
feature_len = length(test);
feature_mat = zeros(length(data_array), feature_len);
label_array = zeros(length(data_array), 1); % 标签
filename_array = repmat("empty", length(data_array), 1); % 文件名

in_use_index = 1;
for i = 1:length(data_array) % 列表前两个是.和..，从第三个开始遍历
    % 排除空文件
    if data_array(i).bytes == 0
        continue;
    end
    [~, baseFileNameNoExt, ~] = fileparts(data_array(i).name);
    disp(i);
    % 获取标签
    filename = data_array(i).name;
    filedscrp = split(filename, '-');
    username = char(filedscrp(1));
    userlabel = str2double(username(5:end));
    path = str2double(filedscrp(2));
    reptime = str2double(filedscrp(3));
    rx_char = char(filedscrp(4));
    rx = str2double(rx_char(2));
    % 排除不需要的路径和接收器及方向数据
    try
%         if rx == 6 && path == 3 && mod(reptime,2) == 0 % 用了翻转以后即便是另外一个方向的也可以检测到，所以不用非得是偶数次
        if rx == 6 && path == 3 && mod(reptime, 2) == 0 % 还是得偶数次...不然特征提取有问题
            feature_mat(in_use_index, :) = extract_file_features_save_figure(fullfile(data_dir, data_array(i).name), data_dir, baseFileNameNoExt);
            label_array(in_use_index) = userlabel;
            filename_array(in_use_index) = fullfile(data_dir, data_array(i).name);
            in_use_index = in_use_index+1;
%         elseif rx == 3 && path == 1 && mod(reptime,2) == 0
        elseif rx == 3 && path == 1 && mod(reptime, 2) == 0
            feature_mat(in_use_index, :) = extract_file_features_save_figure(fullfile(data_dir, data_array(i).name), data_dir, baseFileNameNoExt);
            label_array(in_use_index) = userlabel;
            filename_array(in_use_index) = fullfile(data_dir, data_array(i).name);
            in_use_index = in_use_index+1;
        else
            continue;
        end
    catch
        continue;
    end
end
% 拼接特征和标签
in_use_index = in_use_index-1;
result_mat = horzcat(feature_mat(1:in_use_index, :), label_array(1:in_use_index, :));
filename_array = filename_array(1:in_use_index, :);
filename = cellstr(filename_array); % 方便python读取
save(user_dir, 'result_mat', 'filename');


