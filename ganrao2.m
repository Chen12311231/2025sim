clc;
clear;
close all;

% 参数设置
fs = 80e6;           % 采样频率 80 MHz
T_total = 100e-6;    % 总时长 100 μs
t = 0:1/fs:T_total;  % 时间向量
N = length(t);       % 采样点数

f_j = 20e6;          % 中心频率 20 MHz
DeltaF_n = 10e6;     % 基带噪声带宽 10 MHz
f_de = 1e6;          % 有效调频带宽 1 MHz
U_j = 1;             % 信号幅度

% 生成调制噪声u(t)
rng(0);  % 设置随机种子保证结果可重现
noise_white = randn(1, N);  % 高斯白噪声

% 设计10MHz低通滤波器
order = 6;  % 滤波器阶数
fc = DeltaF_n;  % 截止频率
[b, a] = butter(order, fc/(fs/2), 'low');

% 滤波得到带宽受限的调制噪声
u = filter(b, a, noise_white);
u = u - mean(u);  % 确保零均值
sigma_n = std(u);  % 计算标准差

% 计算调频斜率
K_FM = f_de / sigma_n;

% 数值积分
integral_u = cumsum(u) * (1/fs); 

% 生成随机初始相位
phi0 = 2*pi*rand();

% 生成噪声调频信号
theta = 2*pi*f_j*t + 2*pi*K_FM*integral_u + phi0;
s = U_j * cos(theta);

%% 修改部分：绘制前20μs的时域波形
num_points = round(20e-6 * fs); % 计算20μs对应的采样点数
if num_points > N
    num_points = N; 
end

figure;
plot(t(1:num_points)*1e6, s(1:num_points), 'b');
xlabel('时间 (\mus)');
ylabel('幅度');
title('噪声调频信号时域波形 (前20\mu s)');
grid on;
xlim([0, t(num_points)*1e6]);

%% 计算并绘制功率谱密度
% 使用Welch方法计算功率谱
nfft = 2^nextpow2(N);  % FFT长度
[Pxx, f] = pwelch(s, hamming(1024), 512, nfft, fs, 'centered', 'psd');

% 转换为dBm/Hz (假设阻抗为50欧姆)
Pxx_dBm = 10*log10(Pxx/1e-3 * 50);  % dBm/Hz

% 绘制功率谱
figure;
plot(f/1e6, Pxx_dBm, 'b', 'LineWidth', 1.5);
hold on;

xlabel('频率 (MHz)');
ylabel('功率谱密度 (dBm/Hz)');
title('噪声调频信号功率谱密度');
legend('仿真结果');
grid on;
xlim([-40, 80]);  
ylim([-150, -50]); 

