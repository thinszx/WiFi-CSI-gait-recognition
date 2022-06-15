function [xpeaks,locs] = findallpeaksdistance(x,y,threshold,peakdistance)
% ---------------------------------------------------------------
% x            - 数据横坐标
% y            - 数据纵坐标
% threshold    - 峰值阈值
% peakdistance - 峰值之间的最小间距
% locs         - 返回峰值的真实下标
% xpeaks       - 返回峰值的真实横坐标
% ---------------------------------------------------------------

%% 最大值个数待定
locs =[];

%% 标记众多极大值点
markPeaks = dif(sign(dif(y)));

%% 挖除低于阈值的数据
N = length(x);
for i = 1:N
    if y(i) <= threshold
       y(i) = NaN;
    end
end

%% 估计最多可能出现的峰值数目P
P = 1;
for i = 2:N
    if (y(i)>= threshold) && (markPeaks(i-1)==-2)
        P = P+1;
    end
end

%% 寻找第一个峰值的位置
Peak = threshold;
locs(1) = 1;
for i = 2:N
    if (y(i) > Peak) && (markPeaks(i-1)==-2)
       Peak = y(i);
       locs(1) = i;
    end
end

%% 把峰值附近半径peakdistance以内的数据干掉
M = floor(peakdistance/abs((x(2)- x(1))));
for i = 1:P
    
    if (locs(i)-M>=1) && (locs(i)+M<=N) 
        y(locs(i)-M:locs(i)+M) = NaN;
    elseif (locs(i)-M < 1) && (locs(i)+M<=N)
        y(1:locs(i)+M) = NaN;
    elseif (locs(i)+M > N) && (locs(i)-M>=1)
        y(locs(i)-M:N) = NaN;
    else
        y(1:N) = NaN;
    end
    
    Peak = threshold;
    locs(i+1) = NaN;
    for j = 2:N
        if (y(j) > Peak) && (markPeaks(j-1)==-2)
           Peak = y(j);
           locs(i+1) = j;
        end
    end 
end

%% 真实的峰值个数Q
Q = 0;
for i = 1:P
    if ~isnan(locs(i))
      Q = Q+1;
    end
end
locs = sort(locs);
locs = locs(1:Q);

%% 输出下标对应的真实 x 轴坐标
xpeaks = x(locs);
end
    
%% 一阶导数函数
function dy = dif(y)
dy = zeros(size(y));
dy(1:end-1) = y(2:end) - y(1:end-1);
dy(end) = (y(end-2) + 2*y(end-1)- 3*y(end))/6;
end