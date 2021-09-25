function R = Wavelet(ecg)

ecg_first = ecg(:,1)';
ecg_second = ecg(:,2)';

% 处理一下ecg信号反向的情况
max_first = max(ecg_first);
min_first = min(ecg_first);
mean_first = mean(ecg_first);

if abs(max_first - mean_first) < abs(min_first - mean_first)
    ecg_first = -ecg_first;
end

max_second = max(ecg_second);
min_second = min(ecg_second);
mean_second = mean(ecg_second);

if abs(max_second - mean_second) < abs(min_second - mean_second)
    ecg_second = -ecg_second;
end


fs = 360; %采样频率
range = 25; %连续小波变换尺度系数范围

%% 通道一的结果

%去除信号均值

mean_first_1 = mean(ecg_first);
ecg_first = ecg_first - mean_first_1;
len_first = length(ecg_first);


%% 使用15-25Hz带通滤波器，只留下R波部分，方便寻找R峰

%设计40阶15-25Hz带通FIR滤波器
FIR_1=[0.0041,0.0053,0.0068,0.0080,0.0081,0.0058,-0.0000,-0.0097,-0.0226,...   
   -0.0370,-0.0498,-0.0577,-0.0576,-0.0477,-0.0278,0,0.0318,0.0625,0.0867,...    
    0.1000,0.1000,0.0867,0.0625,0.0318,0,-0.0278,-0.0477,-0.0576,-0.0577,...   
    -0.0498,-0.0370,-0.0226,-0.0097,-0.0000,0.0058,0.0081,0.0080,0.0068,...
    0.0053,0.0041]; % 使用fdatool设计并导出的滤波器系数,带通FIR,15~25Hz,详情使用fdatool打开DS1.fda查看

%为了防止滤波器的边缘效应，需要对原始数据进行延拓
FIR_1_len = length(FIR_1);


%% 滤波

% 通道1 信号滤波
ecg_first_filter = ecg_first;

ecg_first_filter = [ones(1,FIR_1_len) * ecg_first_filter(1) ecg_first_filter ones(1,FIR_1_len) * ecg_first_filter(end)];
ecg_first_filter = filter(FIR_1,1,ecg_first_filter);
ecg_first_filter = ecg_first_filter(FIR_1_len + 1 : len_first + FIR_1_len);

% 通道2 信号滤波
ecg_second_filter = ecg_second;

ecg_second_filter = [ones(1,FIR_1_len) * ecg_second_filter(1) ecg_second_filter ones(1,FIR_1_len) * ecg_second_filter(end)];
ecg_second_filter = filter(FIR_1,1,ecg_second_filter);
ecg_second_filter = ecg_second_filter(FIR_1_len + 1 : len_first + FIR_1_len);


%寻找通道1滤波信号最佳的尺度系数
Cab_first_filter = cwt(ecg_first_filter,1:range,'db3');
[max_val, ~] = max(Cab_first_filter,[],2);
[~,first_a_filter] = max(max_val);


%用最大的因子first_a_filter进行小波变换，结果为wave_first_filter
wave_first_filter = cwt(ecg_first_filter,first_a_filter,'db3');

R_first_max = max(wave_first_filter);

RR_first_interval = 0;

R_first = [];
R_first_val = [];

%% 寻找R峰位置

R_pre = -floor(0.285 * fs);

Peak_mean = 0;

for i = 3 : len_first - 2
    %R峰幅度应大于左右两侧的信号幅度，且大于平均R峰幅度的0.5倍
    if wave_first_filter(i) > wave_first_filter(i-1) && wave_first_filter(i) > wave_first_filter(i-2) && wave_first_filter(i) > wave_first_filter(i+1) && wave_first_filter(i) > wave_first_filter(i+2) && wave_first_filter(i) > 0.1 * R_first_max && wave_first_filter(i) > 0.5 * Peak_mean
        
        if i - R_pre > floor(0.285 * fs) %当前峰值与前一R峰距离超过0.285s，判定当前峰值为R峰，存储
            R_first = [R_first i];
            R_first_val = [R_first_val wave_first_filter(i)];
            R_pre = i;
            Peak_mean = mean(R_first_val);
            
        else %当前峰值与前一R峰距离小于0.285s，更新R峰信息
            if wave_first_filter(i) > wave_first_filter(R_pre);
                R_first(end) = i;
                R_first_val(end) = wave_first_filter(i);
                R_pre = i;
                Peak_mean = mean(R_first_val);
            end
        end
    end
end

%计算R峰平均间隔
RR_first_interval = round(mean(diff(R_first)));

%% 第一个和最后一个R峰可能误判，处理一下
if length(R_first) >= 3
    if R_first(2) - R_first(1) < 0.75 * RR_first_interval
        R_first = R_first(2:end);
        R_first_val = R_first_val(2:end);
    end

    if R_first(end) - R_first(end - 1) < 0.75 * RR_first_interval
        R_first = R_first(1:end - 1);
        R_first_val = R_first_val(1:end - 1);
    end

end


%% 寻找通道1未滤波信号最大的Cab
Cab_first_1 = cwt(ecg_first,1:range,'db3');
[max_val, ~] = max(Cab_first_1,[],2);
[~,first_a1] = max(max_val);

%用最大的因子a进行小波变换，结果为wave_1
wave_first_1 = cwt(ecg_first,first_a1,'db3');

%重新定位R峰位置
for i = 1 : length(R_first)
    if R_first(i) - floor(0.08 * fs) >= 1
        if R_first(i) + floor(0.08 * fs) <= len_first
            [~,ind] = max(wave_first_1(R_first(i) - floor(0.08 * fs):R_first(i) + floor(0.08 * fs)));
            ind = R_first(i) - floor(0.08 * fs) + ind - 1;
        else
            [~,ind] = max(wave_first_1(R_first(i) - floor(0.08 * fs):end));
            ind = R_first(i) - floor(0.08 * fs) + ind - 1;
        end
    else
        [~,ind] = max(wave_first_1(1:R_first(i) + floor(0.08 * fs)));
    end
    
    R_first(i) = ind;
end
    

T_first = [];
T_first_val = [];

P_first = [];
P_first_val = [];

%% 通道1，寻找T峰的粗略位置

%在R峰右侧0.15-0.5s范围内的最大值判定为T峰
for i = R_first
    if i + floor(0.5 * RR_first_interval) < len_first
        
        [~,ind] = max(wave_first_1(i + floor(0.15 * fs):i + floor(0.5 * RR_first_interval)));
    else
        [~,ind] = max(wave_first_1(i + floor(0.15 * fs):end));
    end
    
    ind = i + floor(0.15 * fs) + ind - 1;
    T_first = [T_first ind];
    T_first_val = [T_first_val wave_first_1(ind)];
end

num = length(R_first) + length(T_first);
first_h_1 = (sum(R_first_val) + sum(T_first_val)) / num;


%% 求解 QRS波群的区域
first_m1 = zeros(1,len_first);

for i = 1 : length(R_first)
    left = R_first(i) - 1;
    right = R_first(i) + 1;
    first_m1(R_first(i)) = 1;

    while left > 0 && wave_first_1(left) < wave_first_1(left + 1) 
        first_m1(left) = 1;
        left = left - 1;
    end

    while left > 0 && abs(wave_first_1(left)) > first_h_1
        first_m1(left) = 1;
        left = left - 1;
    end

    while right < len_first && wave_first_1(right) < wave_first_1(right - 1) 
        first_m1(right) = 1;
        right = right + 1;
    end

    while right < len_first && abs(wave_first_1(right)) > first_h_1
        first_m1(right) = 1;
        right = right + 1;
    end
        
end  

%Vqrs_first即通道1的QRS波群范围的信号
Vqrs_first = ecg_first .* first_m1;

ecg_first_2 = ecg_first;
    
%从原始信号中去除QRS波群
for i = 2:len_first
    if first_m1(i) == 1
        ecg_first_2(i) = ecg_first_2(i-1);
    end
end

%寻找最佳的尺度系数，对剩余信号进行CWT
Cab_first_2 = cwt(ecg_first_2,1:range,'db3');
[max_val2, ~] = max(Cab_first_2,[],2);
[~,first_a2] = max(max_val2);

wave_first_2 = cwt(ecg_first_2,first_a2,'db3');

T_first = [];
T_first_val = [];

%% 寻找P峰和T峰的位置
%T峰判定为R峰右侧0.15-0.5s范围内的最大值，P峰判定为R峰左侧 0.1s-0.5RR间隔 范围内的最大值

for i = R_first
    
    if i + floor(0.5 * RR_first_interval) < len_first
        
        [~,ind] = max(wave_first_2(i + floor(0.15 * fs):i + floor(0.5 * RR_first_interval)));
    else
        [~,ind] = max(wave_first_2(i + floor(0.15 * fs):end));
    end
    
    ind = i + floor(0.15 * fs) + ind - 1;
    T_first = [T_first ind];
    T_first_val = [T_first_val wave_first_2(ind)];
    
    
    if i - floor(0.5 * RR_first_interval) >= 1
        
        [~,ind] = max(wave_first_2(i - floor(0.5 * RR_first_interval):i - floor(0.1 * fs) ));
        ind = i - floor(0.5 * RR_first_interval) + ind - 1;
    else
        [~,ind] = max(wave_first_2(1:i - floor(0.1 * fs)));
    end
    
    P_first = [P_first ind];
    P_first_val = [P_first_val wave_first_2(ind)];
end

T_first_mean = mean(ecg_first_2(T_first));
P_first_mean = mean(ecg_first_2(P_first));


num = length(T_first) + length(P_first);
first_h_2 = (sum(T_first_val) + sum(P_first_val)) / num;

first_m2 = zeros(1,len_first);

%% 寻找T峰的所在区域
if T_first_mean > P_first_mean

    for i = 1 : length(T_first)
        left = T_first(i) - 1;
        right = T_first(i) + 1;
        first_m2(T_first(i)) = 1;
        
        while left > 0 && wave_first_2(left) < wave_first_2(left + 1) 
            first_m2(left) = 1;
            left = left - 1;
        end
        
        while left > 0 && abs(wave_first_2(left)) > first_h_2
            first_m2(left) = 1;
            left = left - 1;
        end
        
        while right < len_first && wave_first_2(right) < wave_first_2(right - 1) 
            first_m2(right) = 1;
            right = right + 1;
        end
        
        while right < len_first && abs(wave_first_2(right)) > first_h_2
            first_m2(right) = 1;
            right = right + 1;
        end
        
    end    
else
    for i = 1 : length(P_first)
        left = P_first(i) - 1;
        right = P_first(i) + 1;
        
        first_m2(P_first(i)) = 1;
        
        while left > 0 && wave_first_2(left) < wave_first_2(left + 1) 
            first_m2(left) = 1;
            left = left - 1;
        end
        
        while left > 0 && abs(wave_first_2(left)) > first_h_2
            first_m2(left) = 1;
            left = left - 1;
        end
        
        while right < len_first && wave_first_2(right) < wave_first_2(right - 1) 
            first_m2(right) = 1;
            right = right + 1;
        end
        
        while right < len_first && abs(wave_first_2(right)) > first_h_2
            first_m2(right) = 1;
            right = right + 1;
        end
        
    end  
end


first_Vt = ecg_first_2 .* first_m2;

ecg_first_3 = ecg_first_2;

%将T波从Ecg信号中移除
for i = 2:len_first
    if first_m2(i) == 1
        ecg_first_3(i) = ecg_first_3(i-1);
    end
end


%寻找最佳的尺度系数
Cab_first_3 = cwt(ecg_first_3,1:range,'db3');
[max_val3, ~] = max(Cab_first_3,[],2);
[~,first_a3] = max(max_val3);

%使用该尺度系数对剩余信号进行CWT
wave_first_3 = cwt(ecg_first_3,first_a3,'db3');

first_h_3 = mean(wave_first_3);

first_m3 = zeros(1,len_first);

% 寻找P峰的所在区域
if T_first_mean > P_first_mean
    
    for i = 1 : length(P_first)
        left = P_first(i) - 1;
        right = P_first(i) + 1;
        
        first_m3(P_first(i)) = 1;
        
        while left > 0 && wave_first_3(left) < wave_first_3(left + 1) 
            first_m3(left) = 1;
            left = left - 1;
        end
        
        while left > 0 && abs(wave_first_3(left)) > first_h_2
            first_m3(left) = 1;
            left = left - 1;
        end
        
        while right < len_first && wave_first_3(right) < wave_first_3(right - 1) 
            first_m3(right) = 1;
            right = right + 1;
        end
        
        while right < len_first && abs(wave_first_3(right)) > first_h_2
            first_m3(right) = 1;
            right = right + 1;
        end
        
    end  
     
else
    
    for i = 1 : length(T_first)
        left = T_first(i) - 1;
        right = T_first(i) + 1;
        first_m3(T_first(i)) = 1;
        
        while left > 0 && wave_first_3(left) < wave_first_3(left + 1) 
            first_m3(left) = 1;
            left = left - 1;
        end
        
        while left > 0 && abs(wave_first_3(left)) > first_h_2
            first_m3(left) = 1;
            left = left - 1;
        end
        
        while right < len_first && wave_first_3(right) < wave_first_3(right - 1) 
            first_m3(right) = 1;
            right = right + 1;
        end
        
        while right < len_first && abs(wave_first_3(right)) > first_h_2
            first_m3(right) = 1;
            right = right + 1;
        end
        
    end   
end



%%  通道二的计算  %%%%%%%%%%%%

%去除信号均值
mean_second_1 = mean(ecg_second);
ecg_second = ecg_second - mean_second_1;
len_second = length(ecg_second);


%寻找通道1滤波信号最佳的尺度系数

Cab_second_filter = cwt(ecg_second_filter,1:range,'db3');
[max_val, ~] = max(Cab_second_filter,[],2);
[~,second_a_filter] = max(max_val);


%用最大的因子second_a_filter进行小波变换，结果为wave_second_filter
wave_second_filter = cwt(ecg_second_filter,second_a_filter,'db3');

R_second_max = max(wave_second_filter);

RR_second_interval = 0;
Peak_mean = 0;

R_second = [];
R_second_val = [];
R_pre = -floor(0.285 * fs);

%% 寻找R峰位置
for i = 3 : len_second - 2
    %R峰幅度应大于左右两侧的信号幅度，且大于平均R峰幅度的0.5倍
    if wave_second_filter(i) > wave_second_filter(i-1) && wave_second_filter(i) > wave_second_filter(i-2) && wave_second_filter(i) > wave_second_filter(i+1) && wave_second_filter(i) > wave_second_filter(i+2) && wave_second_filter(i) > 0.1 * R_second_max && wave_second_filter(i) > 0.5 * Peak_mean
        if i - R_pre > floor(0.285 * fs)%当前峰值与前一R峰距离超过0.285s，判定当前峰值为R峰，存储
            R_second = [R_second i];
            R_second_val = [R_second_val wave_second_filter(i)];
            R_pre = i;
            Peak_mean = mean(R_second_val);
        else %当前峰值与前一R峰距离小于0.285s，更新R峰信息
            if wave_second_filter(i) > wave_second_filter(R_pre);
                R_second(end) = i;
                R_second_val(end) = wave_second_filter(i);
                R_pre = i;
                Peak_mean = mean(R_second_val);
            end
        end
    end
end

%计算R峰平均间隔
RR_second_interval = round(mean(diff(R_second)));

%% 第一个和最后一个R峰可能误判，处理一下
if length(R_second) >= 3
    if R_second(2) - R_second(1) < 0.75 * RR_second_interval
        R_second = R_second(2:end);
        R_second_val = R_second_val(2:end);
    end

    if R_second(end) - R_second(end - 1) < 0.75 * RR_second_interval
        R_second = R_second(1:end - 1);
        R_second_val = R_second_val(1:end - 1);
    end
end

%% 寻找通道2未滤波信号的最佳尺度系数
Cab_second_1 = cwt(ecg_second,1:range,'db3');
[max_val, ~] = max(Cab_second_1,[],2);
[~,second_a1] = max(max_val);


%用最大的因子进行小波变换，结果为wave_second_1
wave_second_1 = cwt(ecg_second,second_a1,'db3');


%重定位R峰位置
for i = 1 : length(R_second)
    if R_second(i) - floor(0.08 * fs) >= 1
        if R_second(i) + floor(0.08 * fs) <= len_second
            [~,ind] = max(wave_second_1(R_second(i) - floor(0.08 * fs):R_second(i) + floor(0.08 * fs)));
            ind = R_second(i) - floor(0.08 * fs) + ind - 1;
        else
            [~,ind] = max(wave_second_1(R_second(i) - floor(0.08 * fs):end));
            ind = R_second(i) - floor(0.08 * fs) + ind - 1;
        end
    else
        [~,ind] = max(wave_second_1(1:R_second(i) + floor(0.08 * fs)));
    end
    
    R_second(i) = ind;
end



T_second = [];
T_second_val = [];

P_second = [];
P_second_val = [];

%% 通道2，寻找T峰的粗略位置

%在R峰右侧0.15-0.5s范围内的最大值判定为T峰
for i = R_second
    if i + floor(0.5 * RR_second_interval) < len_second
        
        [~,ind] = max(wave_second_1(i + floor(0.15 * fs):i + floor(0.5 * RR_second_interval)));
    else
        [~,ind] = max(wave_second_1(i + floor(0.15 * fs):end));
    end
    
    ind = i + floor(0.15 * fs) + ind - 1;
    T_second = [T_second ind];
    T_second_val = [T_second_val wave_second_1(ind)];
end

num = length(R_second) + length(T_second);
second_h_1 = (sum(R_second_val) + sum(T_second_val)) / num;


%% 求解 QRS波群的区域
second_m1 = zeros(1,len_second);

for i = 1 : length(R_second)
    left = R_second(i) - 1;
    right = R_second(i) + 1;
    second_m1(R_second(i)) = 1;

    while left > 0 && wave_second_1(left) < wave_second_1(left + 1) 
        second_m1(left) = 1;
        left = left - 1;
    end

    while left > 0 && abs(wave_second_1(left)) > second_h_1
        second_m1(left) = 1;
        left = left - 1;
    end

    while right < len_first && wave_second_1(right) < wave_second_1(right - 1) 
        second_m1(right) = 1;
        right = right + 1;
    end

    while right < len_first && abs(wave_second_1(right)) > second_h_1
        second_m1(right) = 1;
        right = right + 1;
    end
        
end  
 
%Vqrs_second即通道2的QRS波群范围的信号
Vqrs_second = ecg_second .* second_m1;

ecg_second_2 = ecg_second;

%从原始信号中去除QRS波群
for i = 2:len_second
    if second_m1(i) == 1
        ecg_second_2(i) = ecg_second_2(i-1);
    end
end

%寻找最佳的尺度系数，对剩余信号进行CWT
Cab_second_2 = cwt(ecg_second_2,1:range,'db3');
[max_val2, ind2] = max(Cab_second_2,[],2);
[~,second_a2] = max(max_val2);

wave_second_2 = cwt(ecg_second_2,second_a2,'db3');

T_second = [];
T_second_val = [];

%% 寻找P峰和T峰的位置
%T峰判定为R峰右侧0.15-0.5s范围内的最大值，P峰判定为R峰左侧 0.1s-0.5RR间隔 范围内的最大值
for i = R_second
    
    if i + floor(0.5 * RR_second_interval) < len_second
        
        [~,ind] = max(wave_second_2(i + floor(0.15 * fs):i + floor(0.5 * RR_second_interval)));
    else
        [~,ind] = max(wave_second_2(i + floor(0.15 * fs):end));
    end
    
    ind = i + floor(0.15 * fs) + ind - 1;
    T_second = [T_second ind];
    T_second_val = [T_second_val wave_second_2(ind)];
    
    
    if i - floor(0.5 * RR_second_interval) >= 1
        
        [~,ind] = max(wave_second_2(i - floor(0.5 * RR_second_interval):i - floor(0.1 * fs) ));
        ind = i - floor(0.5 * RR_second_interval) + ind - 1;
    else
        [~,ind] = max(wave_second_2(1:i - floor(0.1 * fs)));
    end
    
    P_second = [P_second ind];
    P_second_val = [P_second_val wave_second_2(ind)];
end

T_second_mean = mean(ecg_second_2(T_second));
P_second_mean = mean(ecg_second_2(P_second));

num = length(T_second) + length(P_second);
second_h_2 = (sum(T_second_val) + sum(P_second_val)) / num;

second_m2 = zeros(1,len_second);

%% 寻找T峰的所在区域
if T_second_mean > P_second_mean

    for i = 1 : length(T_second)
        left = T_second(i) - 1;
        right = T_second(i) + 1;
        second_m2(T_second(i)) = 1;
        
        while left > 0 && wave_second_2(left) < wave_second_2(left + 1) 
            second_m2(left) = 1;
            left = left - 1;
        end
        
        while left > 0 && abs(wave_second_2(left)) > second_h_2
            second_m2(left) = 1;
            left = left - 1;
        end
        
        while right < len_second && wave_second_2(right) < wave_second_2(right - 1) 
            second_m2(right) = 1;
            right = right + 1;
        end
        
        while right < len_second && abs(wave_second_2(right)) > second_h_2
            second_m2(right) = 1;
            right = right + 1;
        end
        
    end    
else
    for i = 1 : length(P_second)
        left = P_second(i) - 1;
        right = P_second(i) + 1;
        
        second_m2(P_second(i)) = 1;
        
        while left > 0 && wave_second_2(left) < wave_second_2(left + 1) 
            second_m2(left) = 1;
            left = left - 1;
        end
        
        while left > 0 && abs(wave_second_2(left)) > second_h_2
            second_m2(left) = 1;
            left = left - 1;
        end
        
        while right < len_second && wave_second_2(right) < wave_second_2(right - 1) 
            second_m2(right) = 1;
            right = right + 1;
        end
        
        while right < len_second && abs(wave_second_2(right)) > second_h_2
            second_m2(right) = 1;
            right = right + 1;
        end
        
    end  
end

second_Vt = ecg_second_2 .* second_m2;

ecg_second_3 = ecg_second_2;

%将T波从Ecg信号中移除
for i = 2:len_second
    if second_m2(i) == 1
        ecg_second_3(i) = ecg_second_3(i-1);
    end
end

%寻找最佳的尺度系数
Cab_second_3 = cwt(ecg_second_3,1:range,'db3');
[max_val3, ind3] = max(Cab_second_3,[],2);
[~,second_a3] = max(max_val3);

%使用该尺度系数对剩余信号进行CWT
wave_second_3 = cwt(ecg_second_3,second_a3,'db3');

second_h_3 = mean(abs(wave_second_3));

second_m3 = zeros(1,len_second);

% 寻找P峰的所在区域
if T_second_mean > P_second_mean
    
    for i = 1 : length(P_second)
        left = P_second(i) - 1;
        right = P_second(i) + 1;
        
        second_m3(P_second(i)) = 1;
        
        while left > 0 && wave_second_3(left) < wave_second_3(left + 1) 
            second_m3(left) = 1;
            left = left - 1;
        end
        
        while left > 0 && abs(wave_second_3(left)) > second_h_2
            second_m3(left) = 1;
            left = left - 1;
        end
        
        while right < len_second && wave_second_3(right) < wave_second_3(right - 1) 
            second_m3(right) = 1;
            right = right + 1;
        end
        
        while right < len_second && abs(wave_second_3(right)) > second_h_2
            second_m3(right) = 1;
            right = right + 1;
        end
        
    end  
     
else
    
    for i = 1 : length(T_second)
        left = T_second(i) - 1;
        right = T_second(i) + 1;
        second_m3(T_second(i)) = 1;
        
        while left > 0 && wave_second_3(left) < wave_second_3(left + 1) 
            second_m3(left) = 1;
            left = left - 1;
        end
        
        while left > 0 && abs(wave_second_3(left)) > second_h_2
            second_m3(left) = 1;
            left = left - 1;
        end
        
        while right < len_second && wave_second_3(right) < wave_second_3(right - 1) 
            second_m3(right) = 1;
            right = right + 1;
        end
        
        while right < len_second && abs(wave_second_3(right)) > second_h_2
            second_m3(right) = 1;
            right = right + 1;
        end
        
    end   
end


%% 将两个通道的处理结果综合起来
m1 = first_m1 + second_m1;%通道1的QRS区域和通道2的QRS区域相加
m2 = first_m2 + second_m2;%通道1的T波区域和通道2的T波区域相加
m3 = first_m3 + second_m3;%通道1的P波区域和通道2的P波区域相加

%分别求m1、m2、m3的均值，作为阈值
mean_1 = mean(m1);
mean_2 = mean(m2);
mean_3 = mean(m3);

%大于阈值的范围最为最终的QRS区域、P波区域、T波区域
m1 = m1 > mean_1;
m2 = m2 > mean_2;
m3 = m3 > mean_3;

%处理一下，防止区域混叠
for i = 1 : len_first
    if m1(i) == 1 
        if m2(i) == 1
            m2(i) = 0;
        end
        
        if m3(i) == 1
            m3(i) = 0;
        end
    end
end


%% 在QRS区域内，寻找最终的R峰
i = 1;
m1 = [m1 0];
R = [];

while i <= len_first
    
    if m1(i) == 1
        left = i;
        while i <= len_first && m1(i) == 1
            i = i + 1;
        end
        right = i - 1;
        
        if right - left > 0.05 * fs
            [~,ind] = max(ecg_first(left:right));
            ind = left + ind - 1;
            R = [R ind];
        end
        
    else
        i = i + 1;
    end
end  

    

%% 寻找最终的T峰
T = [];
for i = 1 : length(R)
    
    if R(i) + floor(0.05 * fs) <= len_first
        if R(i) + floor(0.5 * fs) <= len_first
            [~,ind] = max(ecg_first(R(i) + floor(0.05 * fs) : R(i) + floor(0.5 * fs)));
            ind = R(i) + floor(0.05 * fs) + ind - 1;
        else
            [~,ind] = max(ecg_first(R(i) + floor(0.05 * fs) : end));
            ind = R(i) + floor(0.05 * fs) + ind - 1;
        end
    end
    
    T = [T ind];
end


%% 寻找最终的P峰

P = [];
for i = 1 : length(R)
    
    if R(i) - floor(0.05 * fs) >= 1
    
        if R(i) - floor(0.25 * fs) >= 1
            [~,ind] = max(ecg_first(R(i) - floor(0.25 * fs) : R(i) - floor(0.05 * fs)));
            ind = R(i) - floor(0.25 * fs) + ind - 1;
        else
            [~,ind] = max(ecg_first(1 : R(i) - floor(0.05 * fs)));
            
        end
    end
    
    P = [P ind];
end

end

