clc;
clear all;
close all;

%% 1. 参数配置 (Parameters)
bitRate = 1000;         % 比特率 (Hz)
symbol_rate = bitRate;  % 符号率。BPSK调制，一个bit对应一个符号
sps = 16;               % 每个符号的采样点数 (Samples per symbol)
fs = bitRate * sps;     % 采样频率 (Hz)>=2*最高频率
fc = 1000;              % 载波频率 (Hz)
rolloff = 0.5;          % 滚降因子
numBits = 1e6;          % 仿真总比特数

%% 2. 发射端 (Transmitter)
% 生成原始比特流
bitsTx = randi([0, 1], 1, numBits);

% 双极性转换 (0 -> -1, 1 -> 1)
% 在双极性变换中，bit从二进制形式（0和1）被映射到双极性值。数据比特通过这种映射与载波信号相结合以进行传输。
signalBipolar = 2 * bitsTx - 1;

% 上采样，每个符号间插入0点,为脉冲成型做准备
signalUp = upsample(signalBipolar, sps);

% 根升余弦滤波器 (脉冲成型)
%阶数为96=6*sps。滤波器截断为6个符号
rcosFilter = rcosdesign(rolloff, 6, sps);

% 基带信号滤波
txBaseband = filter(rcosFilter, 1, signalUp);

% 载波调制 (上变频)
% 0:length(txBaseband)-1是采样序列
t = (1:length(txBaseband)) .* (1/fs); 
txSignal = txBaseband .* cos(2 * pi * fc.* t);

%% 3. 信道 (Channel)
channel = 1; % 理想信道
ebn0_dB = -6:10; % Eb/N0 范围 (dB)

% 根据 Eb/N0 计算对应的信噪比 SNR
% BPSK中, Es/N0 = Eb/N0
% Es/N0 = Eb/N0+10log10(k)=0.5(Tsym/Tsamp)+SNR(dB)，其中符号周期Tsym=1/symbol_rate，采样周期Tsamp=1/fre_sample。k为1，表示每个符号的信息bit数。
snr_dB = ebn0_dB - 10*log10(0.5 * (1/symbol_rate)/ (1/fs));


%% 4. 接收与误码率计算 (Receiver & BER Calculation)
for i = 1:length(snr_dB)
    % --- 通过AWGN信道 ---
    % 'measured' 选项会先测量信号功率，再根据SNR添加正确功率的噪声
    rxSignal = awgn(channel * txSignal, snr_dB(i), 'measured');

    % --- 接收端处理 ---
    % 相干解调 (下变频)
    rxDemod = rxSignal .* cos(2 * pi * fc .* t);
    
    % 低通滤波器
    LPF = fir1(128,0.2); % 截止频率为0.2*（fs/2）
    rxDemod_lp = filter(LPF,1,rxDemod);
    
    % 匹配滤波 (使用与发射端相同的滤波器)
    rxMatched = filter(rcosFilter, 1, rxDemod_lp);

    % 计算总的滤波器延迟 (发射端+接收端)
    filtDelay = 160;
    
    % 抽样判决
    % 从最佳抽样点提取数据
    rxSampled = rxMatched(filtDelay + 1 : sps : end);
    
    % 判决 (大于0为1, 小于0为0)
    bitsDecision = sign(rxMatched(filtDelay+1:sps:end));
    
    % --- 误码率计算 ---
    % 确保比较的比特数一致
    
    [numErrors(i), berSim(i)] = biterr(bitsTx(1:length(bitsDecision)), (bitsDecision+1)/2);
end

%% 5. 绘图 (Plotting Results)
figure;
semilogy(ebn0_dB, berSim, '-*b', 'LineWidth', 1.5);
xlabel('Eb/N_0 (dB)');
ylabel('BER');
title('BPSK不同Eb/N_0下的误比特率仿真曲线');
grid on;
% legend('仿真误码率 (Simulated BER)');