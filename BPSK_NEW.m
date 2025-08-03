clc;
clear all;
close all;

%% 参数配置 
bitRate = 1000;         % 比特率 
symbol_rate = bitRate;  % 符号率
sps = 16;               % 每个符号的采样点数
fs = bitRate * sps;     % 采样频率 
fc = 1000;              % 载波频率 
rolloff = 0.5;          % 滚降因子
numBits = 1e6;          % 仿真总比特数

% 转发式干扰参数
JSR_dB = 10;            % 干扰信号比 (dB)
T_pr = 0.5e-3;          % 干扰机处理时间 = 0.5 ms (半符号周期)
delay_samples = round(T_pr * fs); % 转换为采样点数
interf_gain = sqrt(10^(JSR_dB/10)); % 干扰幅度增益

fprintf('转发干扰参数:\n');
fprintf('延迟时间: %.2f ms (%d 采样点)\n', T_pr*1000, delay_samples);
fprintf('干信比(JSR): %.1f dB\n', JSR_dB);

%% 2. 发射端 (Transmitter)
bitsTx = randi([0, 1], 1, numBits);
signalBipolar = 2 * bitsTx - 1; % 双极性转换
signalUp = upsample(signalBipolar, sps); % 上采样

% 根升余弦滤波器 
rcosFilter = rcosdesign(rolloff, 6, sps);
txBaseband = filter(rcosFilter, 1, signalUp); % 基带滤波

% 载波调制 
t = (0:length(txBaseband)-1) * (1/fs); % 时间向量从0开始
txSignal = txBaseband .* cos(2 * pi * fc .* t); % 载波调制

%% 生成干扰信号
delayed_signal = [zeros(1, delay_samples), txSignal(1:end-delay_samples)];
delayed_signal = delayed_signal(1:length(txSignal)); % 统一长度
interf_signal = interf_gain * delayed_signal; % 调整干扰功率

%% 误码率计算准备
ebn0_dB = -6:10; % Eb/N0 范围 (dB)
snr_dB = ebn0_dB - 10*log10(0.5 * (1/symbol_rate)/ (1/fs)); % SNR计算

% 初始化结果存储
ber_noInterf = zeros(size(ebn0_dB)); % 无干扰BER
ber_withInterf = zeros(size(ebn0_dB)); % 有干扰BER



%% 主仿真循环
for i = 1:length(ebn0_dB)
    % === 无干扰情况 (与原始代码一致) ===
    rx_noInterf = awgn(txSignal, snr_dB(i), 'measured');
    [ber_noInterf(i), ~] = process_receiver(rx_noInterf, rcosFilter, sps, fs, fc, bitsTx);
    
    % === 有干扰情况 ===
    % 叠加干扰信号
    rx_withInterf = txSignal + interf_signal;
    % 添加AWGN噪声
    rx_withInterf = awgn(rx_withInterf, snr_dB(i), 'measured');
    % 接收处理 (使用相同处理函数)
    [ber_withInterf(i), ~] = process_receiver(rx_withInterf, rcosFilter, sps, fs, fc, bitsTx);

end

%% 7. 结果可视化
figure;
semilogy(ebn0_dB, ber_noInterf, '-ob', 'LineWidth', 1.5);
hold on;
xlabel('Eb/N_0 (dB)');
ylabel('误码率 (BER)');
title('BPSK系统性能: 转发式干扰影响');
legend('无干扰', ['转发干扰 (JSR=', num2str(JSR_dB), 'dB, 延时=', num2str(T_pr*1000), 'ms)'], 'Location', 'best');

% 添加理论BER曲线作为参考
berTheory = berawgn(ebn0_dB, 'psk', 2, 'nondiff');
semilogy(ebn0_dB, berTheory, 'k--', 'LineWidth', 1.5);
legend('有干扰(仿真)', '无干扰(仿真)');

%% 信号可视化 
figure;
segment = 5000:5200; % 选择一段信号
t_segment = t(segment);

subplot(3,1,1);
plot(t_segment, txSignal(segment));
title('原始BPSK信号');
ylabel('幅度');

subplot(3,1,2);
plot(t_segment, interf_signal(segment));
title('转发干扰信号');
ylabel('幅度');

subplot(3,1,3);
plot(t_segment, txSignal(segment) + interf_signal(segment));
title('叠加干扰后的信号');
xlabel('时间 (s)');
ylabel('幅度');

%% 接收处理函数 (与原始代码一致)
function [ber, bitsDecision] = process_receiver(rxSignal, rcosFilter, sps, fs, fc, bitsTx)
    t = (0:length(rxSignal)-1) * (1/fs);
    
    % 相干解调 
    rxDemod = rxSignal .* cos(2 * pi * fc .* t);
    
    % 低通滤波器 
    LPF = fir1(128, 0.2); % 截止频率为0.2*(fs/2)
    rxDemod_lp = filter(LPF, 1, rxDemod);
    
    % 匹配滤波 
    rxMatched = filter(rcosFilter, 1, rxDemod_lp);
    
    % 动态计算滤波器延迟 
    delay_rcos = (length(rcosFilter)-1)/2; % 根升余弦滤波器群延迟
    delay_lpf = (length(LPF)-1)/2; % 低通滤波器群延迟
    filtDelay = round(delay_rcos + delay_lpf); % 总延迟
    
    % 抽样判决 
    rxSampled = rxMatched(filtDelay + 1 : sps : end);
    bitsDecision = sign(rxSampled); % 双极性判决
    
    % 误码率计算 
    valid_bits = min(length(bitsDecision), length(bitsTx));
    [~, ber] = biterr(bitsTx(1:valid_bits), (bitsDecision+1)/2);
end