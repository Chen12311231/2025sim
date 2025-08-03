function J = generate_FM_noise_interference(fc, deltaFn, f_de, fs, numSamples, Pj)
% 生成噪声调频干扰信号
% 输入:
%   fc - 中心频率 (Hz)
%   deltaFn - 基带噪声带宽 (Hz)
%   f_de - 有效调频带宽 (Hz)
%   fs - 采样频率 (Hz)
%   numSamples - 采样点数
%   Pj - 干扰功率 (W)
% 输出:
%   J - 噪声调频干扰信号

% 时间向量
t = (0:numSamples-1)/fs;

% 生成高斯白噪声
noise_white = randn(1, numSamples);

% 设计低通滤波器（限制噪声带宽）
order = 6;  % 滤波器阶数
[b, a] = butter(order, deltaFn/(fs/2), 'low');

% 滤波得到带宽受限的调制噪声
u = filter(b, a, noise_white);
u = u - mean(u);  % 确保零均值
sigma_n = std(u);  % 计算标准差

% 计算调频斜率
K_FM = f_de / sigma_n;

% 数值积分
integral_u = cumsum(u) * (1/fs);  % ∫u(t')dt'

% 生成随机初始相位
phi0 = 2*pi*rand();

% 生成噪声调频信号
theta = 2*pi*fc*t + 2*pi*K_FM*integral_u + phi0;
J = cos(theta);  % 幅度为1的干扰信号

% 调整干扰信号功率
current_power = mean(J.^2);
scaling_factor = sqrt(Pj / current_power);
J = scaling_factor * J;
end