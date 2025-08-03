function J = generate_AM_noise_interference(fc, deltaFn, fs, numSamples, Pj)
% 生成噪声调幅干扰信号
% 输入:
%   fc - 中心频率 (Hz)
%   deltaFn - 基带噪声带宽 (Hz)
%   fs - 采样频率 (Hz)
%   numSamples - 信号长度 (采样点数)
%   Pj - 干扰功率 (W)
%
% 输出:
%   J - 噪声调幅干扰信号

% 计算时间向量
t = (0:numSamples-1)/fs;

% 计算噪声参数：根据功率平衡和约束条件
% 总功率Pj = U0^2/2 + sigma_n^2/2
% 取U0 = 3*sigma_n，则Pj = (9*sigma_n^2)/2 + sigma_n^2/2 = 5*sigma_n^2
sigma_n_sq = Pj / 5;   % 噪声方差
U0 = 3 * sqrt(sigma_n_sq);  % 载波电平

% 生成调制噪声Un(t)：零均值，方差为sigma_n_sq，带宽限制在[0, deltaFn] Hz
% 设计低通滤波器（FIR）以产生带限噪声
nyquist = fs/2;
cutoff = deltaFn / nyquist; % 归一化截止频率
filter_order = 1000; % 滤波器阶数
b = fir1(filter_order, cutoff, 'low'); % FIR低通滤波器

% 生成高斯白噪声（多生成一些以补偿滤波器延迟）
x = randn(1, numSamples + filter_order);
% 滤波
Un = filter(b, 1, x);
% 去除滤波器延迟
Un = Un(filter_order/2 + 1: filter_order/2 + numSamples);
Un = Un - mean(Un); % 确保均值为0

% 调整方差
current_var = var(Un);
Un = Un * sqrt(sigma_n_sq / current_var);

% 随机相位
phi = 2*pi*rand();

% 生成载波
carrier = cos(2*pi*fc*t + phi);

% 合成干扰信号
J = (U0 + Un) .* carrier;
end