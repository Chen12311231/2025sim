clear; clc; close all;


% 参数设置
fj = 20e6;          % 中心频率 (20 MHz)
deltaFn = 10e6;     % 基带噪声带宽 (10 MHz)
fs = 80e6;          % 采样频率 (80 MHz)
Pj = 10;            % 干扰信号功率 (10 W)
duration = 200e-6;  % 信号持续时间 (200 μs)

% 计算采样点数和时间向量
N = round(fs * duration);  % 采样点数
t = (0:N-1)/fs;           % 时间向量 (s)

% 计算噪声参数
sigma_n_sq = 2;           % 噪声方差 σn²
U0 = sqrt(18);            % 载波电平 

% 生成带限高斯白噪声 (Un(t))
% 生成高斯白噪声
x = randn(1, N);          

% 设计低通滤波器 
nyquist = fs/2;
cutoff = deltaFn / nyquist;
filter_order = 1000;       % 滤波器阶数
b = fir1(filter_order, cutoff, 'low');  % FIR低通滤波器

%  滤波并调整功率
Un = filter(b, 1, x);      % 滤波
Un = Un - mean(Un);        % 确保均值为0
current_power = var(Un);   % 当前方差
Un = Un * sqrt(sigma_n_sq / current_power); % 调整到目标方差

% 生成随机相位
phi = 2*pi*rand();         % [0, 2π]均匀分布的随机相位

% 生成载波信号
carrier = cos(2*pi*fj*t + phi);

% 合成噪声调幅干扰信号
J = (U0 + Un) .* carrier;

% 绘制时域波形 (前200 μs)
figure;
plot(t*1e6, J);            
xlabel('时间 (\mus)');
ylabel('幅度');
title('噪声调幅干扰信号时域波形 (前200 \mus)');
xlim([0, 200]);
grid on;

% 计算并绘制频谱
[pxx, f] = pwelch(J, hann(4096), 2048, 4096, fs, 'centered', 'power');

figure;
plot(f/1e6, 10*log10(pxx)); % 频率转换为MHz，功率转换为dB
xlabel('频率 (MHz)');
ylabel('功率谱密度 (dB/Hz)');
title('噪声调幅干扰信号功率谱');
xlim([-40, 40]);
grid on;

