function R = Wavelet(ecg)

ecg_first = ecg(:,1)';
ecg_second = ecg(:,2)';

% ����һ��ecg�źŷ�������
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


fs = 360; %����Ƶ��
range = 25; %����С���任�߶�ϵ����Χ

%% ͨ��һ�Ľ��

%ȥ���źž�ֵ

mean_first_1 = mean(ecg_first);
ecg_first = ecg_first - mean_first_1;
len_first = length(ecg_first);


%% ʹ��15-25Hz��ͨ�˲�����ֻ����R�����֣�����Ѱ��R��

%���40��15-25Hz��ͨFIR�˲���
FIR_1=[0.0041,0.0053,0.0068,0.0080,0.0081,0.0058,-0.0000,-0.0097,-0.0226,...   
   -0.0370,-0.0498,-0.0577,-0.0576,-0.0477,-0.0278,0,0.0318,0.0625,0.0867,...    
    0.1000,0.1000,0.0867,0.0625,0.0318,0,-0.0278,-0.0477,-0.0576,-0.0577,...   
    -0.0498,-0.0370,-0.0226,-0.0097,-0.0000,0.0058,0.0081,0.0080,0.0068,...
    0.0053,0.0041]; % ʹ��fdatool��Ʋ��������˲���ϵ��,��ͨFIR,15~25Hz,����ʹ��fdatool��DS1.fda�鿴

%Ϊ�˷�ֹ�˲����ı�ԵЧӦ����Ҫ��ԭʼ���ݽ�������
FIR_1_len = length(FIR_1);


%% �˲�

% ͨ��1 �ź��˲�
ecg_first_filter = ecg_first;

ecg_first_filter = [ones(1,FIR_1_len) * ecg_first_filter(1) ecg_first_filter ones(1,FIR_1_len) * ecg_first_filter(end)];
ecg_first_filter = filter(FIR_1,1,ecg_first_filter);
ecg_first_filter = ecg_first_filter(FIR_1_len + 1 : len_first + FIR_1_len);

% ͨ��2 �ź��˲�
ecg_second_filter = ecg_second;

ecg_second_filter = [ones(1,FIR_1_len) * ecg_second_filter(1) ecg_second_filter ones(1,FIR_1_len) * ecg_second_filter(end)];
ecg_second_filter = filter(FIR_1,1,ecg_second_filter);
ecg_second_filter = ecg_second_filter(FIR_1_len + 1 : len_first + FIR_1_len);


%Ѱ��ͨ��1�˲��ź���ѵĳ߶�ϵ��
Cab_first_filter = cwt(ecg_first_filter,1:range,'db3');
[max_val, ~] = max(Cab_first_filter,[],2);
[~,first_a_filter] = max(max_val);


%����������first_a_filter����С���任�����Ϊwave_first_filter
wave_first_filter = cwt(ecg_first_filter,first_a_filter,'db3');

R_first_max = max(wave_first_filter);

RR_first_interval = 0;

R_first = [];
R_first_val = [];

%% Ѱ��R��λ��

R_pre = -floor(0.285 * fs);

Peak_mean = 0;

for i = 3 : len_first - 2
    %R�����Ӧ��������������źŷ��ȣ��Ҵ���ƽ��R����ȵ�0.5��
    if wave_first_filter(i) > wave_first_filter(i-1) && wave_first_filter(i) > wave_first_filter(i-2) && wave_first_filter(i) > wave_first_filter(i+1) && wave_first_filter(i) > wave_first_filter(i+2) && wave_first_filter(i) > 0.1 * R_first_max && wave_first_filter(i) > 0.5 * Peak_mean
        
        if i - R_pre > floor(0.285 * fs) %��ǰ��ֵ��ǰһR����볬��0.285s���ж���ǰ��ֵΪR�壬�洢
            R_first = [R_first i];
            R_first_val = [R_first_val wave_first_filter(i)];
            R_pre = i;
            Peak_mean = mean(R_first_val);
            
        else %��ǰ��ֵ��ǰһR�����С��0.285s������R����Ϣ
            if wave_first_filter(i) > wave_first_filter(R_pre);
                R_first(end) = i;
                R_first_val(end) = wave_first_filter(i);
                R_pre = i;
                Peak_mean = mean(R_first_val);
            end
        end
    end
end

%����R��ƽ�����
RR_first_interval = round(mean(diff(R_first)));

%% ��һ�������һ��R��������У�����һ��
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


%% Ѱ��ͨ��1δ�˲��ź�����Cab
Cab_first_1 = cwt(ecg_first,1:range,'db3');
[max_val, ~] = max(Cab_first_1,[],2);
[~,first_a1] = max(max_val);

%����������a����С���任�����Ϊwave_1
wave_first_1 = cwt(ecg_first,first_a1,'db3');

%���¶�λR��λ��
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

%% ͨ��1��Ѱ��T��Ĵ���λ��

%��R���Ҳ�0.15-0.5s��Χ�ڵ����ֵ�ж�ΪT��
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


%% ��� QRS��Ⱥ������
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

%Vqrs_first��ͨ��1��QRS��Ⱥ��Χ���ź�
Vqrs_first = ecg_first .* first_m1;

ecg_first_2 = ecg_first;
    
%��ԭʼ�ź���ȥ��QRS��Ⱥ
for i = 2:len_first
    if first_m1(i) == 1
        ecg_first_2(i) = ecg_first_2(i-1);
    end
end

%Ѱ����ѵĳ߶�ϵ������ʣ���źŽ���CWT
Cab_first_2 = cwt(ecg_first_2,1:range,'db3');
[max_val2, ~] = max(Cab_first_2,[],2);
[~,first_a2] = max(max_val2);

wave_first_2 = cwt(ecg_first_2,first_a2,'db3');

T_first = [];
T_first_val = [];

%% Ѱ��P���T���λ��
%T���ж�ΪR���Ҳ�0.15-0.5s��Χ�ڵ����ֵ��P���ж�ΪR����� 0.1s-0.5RR��� ��Χ�ڵ����ֵ

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

%% Ѱ��T�����������
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

%��T����Ecg�ź����Ƴ�
for i = 2:len_first
    if first_m2(i) == 1
        ecg_first_3(i) = ecg_first_3(i-1);
    end
end


%Ѱ����ѵĳ߶�ϵ��
Cab_first_3 = cwt(ecg_first_3,1:range,'db3');
[max_val3, ~] = max(Cab_first_3,[],2);
[~,first_a3] = max(max_val3);

%ʹ�øó߶�ϵ����ʣ���źŽ���CWT
wave_first_3 = cwt(ecg_first_3,first_a3,'db3');

first_h_3 = mean(wave_first_3);

first_m3 = zeros(1,len_first);

% Ѱ��P�����������
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



%%  ͨ�����ļ���  %%%%%%%%%%%%

%ȥ���źž�ֵ
mean_second_1 = mean(ecg_second);
ecg_second = ecg_second - mean_second_1;
len_second = length(ecg_second);


%Ѱ��ͨ��1�˲��ź���ѵĳ߶�ϵ��

Cab_second_filter = cwt(ecg_second_filter,1:range,'db3');
[max_val, ~] = max(Cab_second_filter,[],2);
[~,second_a_filter] = max(max_val);


%����������second_a_filter����С���任�����Ϊwave_second_filter
wave_second_filter = cwt(ecg_second_filter,second_a_filter,'db3');

R_second_max = max(wave_second_filter);

RR_second_interval = 0;
Peak_mean = 0;

R_second = [];
R_second_val = [];
R_pre = -floor(0.285 * fs);

%% Ѱ��R��λ��
for i = 3 : len_second - 2
    %R�����Ӧ��������������źŷ��ȣ��Ҵ���ƽ��R����ȵ�0.5��
    if wave_second_filter(i) > wave_second_filter(i-1) && wave_second_filter(i) > wave_second_filter(i-2) && wave_second_filter(i) > wave_second_filter(i+1) && wave_second_filter(i) > wave_second_filter(i+2) && wave_second_filter(i) > 0.1 * R_second_max && wave_second_filter(i) > 0.5 * Peak_mean
        if i - R_pre > floor(0.285 * fs)%��ǰ��ֵ��ǰһR����볬��0.285s���ж���ǰ��ֵΪR�壬�洢
            R_second = [R_second i];
            R_second_val = [R_second_val wave_second_filter(i)];
            R_pre = i;
            Peak_mean = mean(R_second_val);
        else %��ǰ��ֵ��ǰһR�����С��0.285s������R����Ϣ
            if wave_second_filter(i) > wave_second_filter(R_pre);
                R_second(end) = i;
                R_second_val(end) = wave_second_filter(i);
                R_pre = i;
                Peak_mean = mean(R_second_val);
            end
        end
    end
end

%����R��ƽ�����
RR_second_interval = round(mean(diff(R_second)));

%% ��һ�������һ��R��������У�����һ��
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

%% Ѱ��ͨ��2δ�˲��źŵ���ѳ߶�ϵ��
Cab_second_1 = cwt(ecg_second,1:range,'db3');
[max_val, ~] = max(Cab_second_1,[],2);
[~,second_a1] = max(max_val);


%���������ӽ���С���任�����Ϊwave_second_1
wave_second_1 = cwt(ecg_second,second_a1,'db3');


%�ض�λR��λ��
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

%% ͨ��2��Ѱ��T��Ĵ���λ��

%��R���Ҳ�0.15-0.5s��Χ�ڵ����ֵ�ж�ΪT��
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


%% ��� QRS��Ⱥ������
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
 
%Vqrs_second��ͨ��2��QRS��Ⱥ��Χ���ź�
Vqrs_second = ecg_second .* second_m1;

ecg_second_2 = ecg_second;

%��ԭʼ�ź���ȥ��QRS��Ⱥ
for i = 2:len_second
    if second_m1(i) == 1
        ecg_second_2(i) = ecg_second_2(i-1);
    end
end

%Ѱ����ѵĳ߶�ϵ������ʣ���źŽ���CWT
Cab_second_2 = cwt(ecg_second_2,1:range,'db3');
[max_val2, ind2] = max(Cab_second_2,[],2);
[~,second_a2] = max(max_val2);

wave_second_2 = cwt(ecg_second_2,second_a2,'db3');

T_second = [];
T_second_val = [];

%% Ѱ��P���T���λ��
%T���ж�ΪR���Ҳ�0.15-0.5s��Χ�ڵ����ֵ��P���ж�ΪR����� 0.1s-0.5RR��� ��Χ�ڵ����ֵ
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

%% Ѱ��T�����������
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

%��T����Ecg�ź����Ƴ�
for i = 2:len_second
    if second_m2(i) == 1
        ecg_second_3(i) = ecg_second_3(i-1);
    end
end

%Ѱ����ѵĳ߶�ϵ��
Cab_second_3 = cwt(ecg_second_3,1:range,'db3');
[max_val3, ind3] = max(Cab_second_3,[],2);
[~,second_a3] = max(max_val3);

%ʹ�øó߶�ϵ����ʣ���źŽ���CWT
wave_second_3 = cwt(ecg_second_3,second_a3,'db3');

second_h_3 = mean(abs(wave_second_3));

second_m3 = zeros(1,len_second);

% Ѱ��P�����������
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


%% ������ͨ���Ĵ������ۺ�����
m1 = first_m1 + second_m1;%ͨ��1��QRS�����ͨ��2��QRS�������
m2 = first_m2 + second_m2;%ͨ��1��T�������ͨ��2��T���������
m3 = first_m3 + second_m3;%ͨ��1��P�������ͨ��2��P���������

%�ֱ���m1��m2��m3�ľ�ֵ����Ϊ��ֵ
mean_1 = mean(m1);
mean_2 = mean(m2);
mean_3 = mean(m3);

%������ֵ�ķ�Χ��Ϊ���յ�QRS����P������T������
m1 = m1 > mean_1;
m2 = m2 > mean_2;
m3 = m3 > mean_3;

%����һ�£���ֹ������
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


%% ��QRS�����ڣ�Ѱ�����յ�R��
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

    

%% Ѱ�����յ�T��
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


%% Ѱ�����յ�P��

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

