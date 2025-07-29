clc;
clear;
close all;

%% AFDM系统参数
fc = 4e9;               % 载波频率 (Hz)
delt_f = 15e3;          % 子载波间隔 (Hz)
N = 64;                 % 系统子载波数
BW = N * delt_f;        % 系统带宽 (Hz)
ts = 1/BW;              % 采样间隔 (s)
N_simulation = 100;     % 蒙特卡洛仿真次数

data_positon = [];      % 数据符号位置

Nd = 6;                 % 每帧包含的AFDM符号数
N_frm = 10;             % 每种信噪比下的仿真帧数
M = 4;                  % QPSK调制
SNR = 0:2:20;           % 信噪比范围(dB)

%% AFDM核心参数
alpha_max = 3;          % 最大归一化多普勒频移 (为覆盖信道文件中的3.5，设为4)
l_max = 2;              % 最大时延抽头 (为覆盖信道文件中的3，设为3)

% --- 新增功能开关 ---
is_fractional_doppler = true; % true: 执行分数多普勒估计; false: 执行原有的整数多普勒估计

% 根据Bemani论文(57)式，Q的定义需要考虑分数多普勒扩展
if is_fractional_doppler
    xi_nu = 0; % 分数多普勒扩展因子, 对应论文中的 k_v 或 ξ_ν
    c1 = (2*(alpha_max + xi_nu) + 1)/(2*N);
    Q = (l_max + 1)*(2*(alpha_max + xi_nu) + 1) - 1;
else
    c1 = (2*alpha_max + 1)/(2*N); % 根据文献式(47)
    Q = (l_max + 1)*(2*alpha_max + 1) - 1; % 保护间隔长度
end
c2 = 0;                 % 任意小数
L_cpp = l_max + 1;      % CPP长度


%% 生成DAFT/IDAFT矩阵
n = 0:N-1;
Lambda_c1 = diag(exp(-1j*2*pi*c1*n.^2));
Lambda_c2 = diag(exp(-1j*2*pi*c2*n.^2));
F = dftmtx(N)/sqrt(N);
A = Lambda_c2 * F * Lambda_c1;  % DAFT矩阵
A_H = A';                       % IDAFT矩阵

%% 生成基带信号
P_data = randi([0 1], 1, N * Nd * N_frm); % 产生基带bit流数据

data_temp1 = reshape(P_data, log2(M), [])'; % 以每组2比特进行分组
data_temp2 = bi2de(data_temp1);             % 二进制转化为十进制
qpsk_data = pskmod(data_temp2, M, pi/M);   % QPSK调制

%% 插入导频（单导频+保护符号）
X_pilot = 3+3j;            % 导频值

X_pilot_station = 1;       % 导频符号位置（文献中0位置，这里从1开始）
pilot_num = 1;             % 导频数量

% 数据符号位置（确保与导频和保护间隔不重叠）
% Bemani论文图7将导频置于0，后跟Q个保护符号 [cite: 957, 959]
% 这里为了兼容原有结构，将数据放在导频和保护区之外
for img = 1:N
    if img - X_pilot_station >= Q + 2 && img - X_pilot_station <= N - Q
        data_positon = [data_positon, img];
    end
end


data_row = length(data_positon);
data_col = ceil(length(qpsk_data)/data_row);

% 初始化符号矩阵
pilot_seq = ones(pilot_num, data_col) * X_pilot; % 导频序列
data = zeros(N, data_col); % 预设整个矩阵
data(X_pilot_station, :) = pilot_seq; % 插入导频

% 补零处理数据
if data_row * data_col > length(qpsk_data)
    data2 = [qpsk_data; zeros(data_row*data_col-length(qpsk_data), 1)]; % 补0
else
    data2 = qpsk_data;
end

%% 串并转换
data_seq = reshape(data2, data_row, data_col);
data(data_positon, :) = data_seq; % 将导频与数据合并

%% IDAFT
idaft_data = A_H * data;

%% 插入循环前缀
n_cp = (-L_cpp:-1)';
cpp_phase = exp(-1j*2*pi*c1*(N^2 + 2*N*n_cp));
cpp = idaft_data(end-L_cpp+1:end,:) .* cpp_phase;
Tx_cd = [cpp; idaft_data];

%% 并串转换
Tx_data = reshape(Tx_cd, [], 1);
x = Tx_data;

%% 上变频
N_total = length(x);
tt = ts*(0:N_total-1);
x = x .* exp(1j*2*pi*fc*tt).';

%% 经过信道
multi_signal = multi_channel(x);

%% 模拟
Ber_record = zeros(1, length(SNR));
for i_simulation = 1:N_simulation
    tic
    for jj = 1:length(SNR)
        %% 添加AWGN噪声
        rx_channel = multi_signal;
        rx_channel = awgn(multi_signal, SNR(jj), 'measured');
       
        %% 下变频
        rx_channel = rx_channel .* exp(-1j*2*pi*fc*tt).';

        %% 串并转换
        Rx_data1 = reshape(rx_channel, N+L_cpp, []);

        %% 去掉循环前缀
        Rx_data2 = Rx_data1(L_cpp+1:end, :);

        %% DAFT
        daft_data = A * Rx_data2;
        y = daft_data;

        %% 信道估计
        % 根据Bemani论文图8，提取与导频相关的接收符号 [cite: 962]
        if is_fractional_doppler
             ind_r_E = [1:(alpha_max + xi_nu + 1), (N - Q + 1):(N)];
        else
             ind_r_E = [1:alpha_max+1, N-Q+alpha_max+1:N];
        end

        y_E = y(ind_r_E, :);

        P = 3; % 信道文件中的路径数
        
        l_est = zeros(P, 1);
        alpha_est = zeros(P, 1);
        a_est = zeros(P, 1);
        h_est = zeros(P, data_col);
        
        if ~is_fractional_doppler
            % ================================================================
            % 您原有的整数多普勒估计算法 (基于Yin & Tang论文)
            % ================================================================
            best_params = zeros(P, 2); % [l_i, alpha_i]
            hi = zeros(P, data_col); % 信道增益

            % 找到能量最大的 P 个路径
            [~, sort_idx] = sort(abs(y_E).^2, 'descend');
            top_P_idx = sort_idx(1:P, :);

            for p = 1:P
                for col = 1:data_col
                    idx = top_P_idx(p, col);
                    m = ind_r_E(idx) - 1;

                    if m >= 0 && m <= alpha_max
                        l_i = 0;
                        alpha_i = -m;
                    elseif m >= N - alpha_max && m <= N - 1
                        l_i = 0;
                        alpha_i = N - m;
                    else
                        for d = 1:l_max
                            if m >= N - alpha_max - 2*N*c1*d && m <= N - alpha_max - 2*N*c1*(d-1) - 1
                                l_i = d;
                                alpha_i = N - m - 2*N*c1*d;
                                break;
                            end
                        end
                    end

                    E = exp(1j * (2*pi/N) * (N*c1*l_i^2 - N*c2*m^2));
                    hi(p, col) = y_E(idx, col) / (X_pilot * E) ;
                    if col == 1 % 仅记录第一帧的参数
                        best_params(p, :) = [l_i, alpha_i];
                    end
                end
            end
            l_est = best_params(:, 1);
            alpha_est = best_params(:, 2);
            h_est = hi;

        else
            % ================================================================
            % 新增：分数多普勒信道估计算法 (遵循Bemani论文VI-B节)
            % ================================================================

            for col = 1:data_col
                y_E_col = y_E(:,col);
                % --- 阶段1: 粗估计(l_i, α_i) -- 使用Yin&Tang的快速映射法 ---
                [~, sort_idx] = sort(abs(y_E_col).^2, 'descend');
                top_P_idx = sort_idx(1:P);
                l_est_col = zeros(P,1);
                alpha_est_col = zeros(P,1);

                for p_idx = 1:P
                    idx = top_P_idx(p_idx);
                    m = ind_r_E(idx) - 1;
                    if m >= 0 && m <= alpha_max, l_i = 0; alpha_i = -m;
                    elseif m >= N - alpha_max && m <= N - 1, l_i = 0; alpha_i = N - m;
                    else
                        for d = 1:l_max
                            if m >= N - alpha_max - 2*N*c1*d && m <= N - alpha_max - 2*N*c1*(d-1) - 1
                                l_i = d; alpha_i = N - m - 2*N*c1*d; break;
                            end
                        end
                    end
                    l_est_col(p_idx) = l_i;
                    alpha_est_col(p_idx) = alpha_i;
                end
            end
            % --- 阶段2: 精调(a_i) -- 使用Bemani的网格搜索法 ---
            a_search_grid = -0.5:0.01:0.5;
            a_est_col = zeros(P,1);
            for p_idx = 1:P
                l_p = l_est_col(p_idx); alpha_p = alpha_est_col(p_idx);
                frac_corr = zeros(length(a_search_grid), 1);
                for a_idx = 1:length(a_search_grid)
                    nu = alpha_p + a_search_grid(a_idx);
                    h_pilot_resp = calculate_pilot_response(l_p, nu, N, c1, c2, ind_r_E);
                    frac_corr(a_idx) = abs(h_pilot_resp' * y_E_col)^2;
                end
                [~, max_idx] = max(frac_corr);
                a_est_col(p_idx) = a_search_grid(max_idx);
            end

            % --- 阶段3: 计算信道增益(h_i) ---
            nu_est_col = alpha_est_col + a_est_col;
            G = zeros(P, P); z = zeros(P, 1);
            h_pilot_responses = cell(P, 1);
            for p_idx = 1:P, h_pilot_responses{p_idx} = calculate_pilot_response(l_est_col(p_idx), nu_est_col(p_idx), N, c1, c2, ind_r_E); end
            for i_p = 1:P
                z(i_p) = (h_pilot_responses{i_p}' * y_E_col) / X_pilot;
                for j_p = 1:P, G(i_p, j_p) = h_pilot_responses{i_p}' * h_pilot_responses{j_p}; end
            end
            h_est_col = (G + 1e-6*eye(P)) \ z;

            if col == 1, l_est = l_est_col; alpha_est = alpha_est_col; a_est = a_est_col; end
            h_est(:, col) = h_est_col;
        end
    

        


        %% 构建完整信道矩阵
        data_aftereq = zeros(length(data_positon), data_col);
        nu_est = alpha_est + a_est; % 完整的归一化多普勒

        if ~is_fractional_doppler
        % ===============================================
        % 方案A: 针对整数多普勒的高效实现 O(PN)
        % ===============================================
        for col = 1:data_col
            Heff = zeros(N, N);
            for i = 1:P
                li = l_est(i);
                alpha_i = alpha_est(i);
                loc_i = mod(alpha_i + 2*N*c1*li, N);
                for p = 0:N-1
                    q = mod(p + loc_i, N);
                    phase = 2*pi/N * (N*c1*li^2 - q*li + N*c2*(q^2 - p^2));
                    Hi = exp(1j * phase);
                    Heff(p+1, q+1) = Heff(p+1, q+1) + h_est(i, col) * Hi;
                end
            end
            
            % --- 均衡 ---
            H_eff = Heff(data_positon, data_positon);
            snr_linear = 10^(SNR(jj)/10); 
            W_LMMSE = (H_eff' * H_eff + (1/snr_linear) * eye(length(data_positon))) \ H_eff';
            y_sym = daft_data(data_positon, col);
            x_hat = W_LMMSE * y_sym;
            x_hat(~isfinite(x_hat)) = 0;
            data_aftereq(:, col) = x_hat;
        end
    else
        % ===============================================
        % 方案B: 针对分数多普勒的通用实现 O(PN^2)
        % ===============================================
        nu_est = alpha_est + a_est;
        k_v = xi_nu; % 使用我们之前定义的扩展因子作为带宽参数
        
        for col = 1:data_col
            Heff = zeros(N, N);
            for i = 1:P
                li = l_est(i);
                nu_i = nu_est(i);
                
                % 计算每个路径的中心位置 loc_i
                % 注意：这里的中心位置应基于整数多普勒 alpha_i
                alpha_i_est = alpha_est(i);
                loc_i = mod(alpha_i_est + 2*N*c1*li, N);
                
                % 遍历每一行 p
                for p = 0:N-1
                    % 找到能量带的中心 q_center
                    q_center = mod(p + loc_i, N);
                    
                    % 定义能量带的偏移量，从 -k_v 到 +k_v
                    q_offsets = -k_v:k_v;
                    
                    % 计算需要计算的 q 的索引，并处理环绕效应
                    q_indices = mod(q_center + q_offsets, N);
                    
                    % 只在能量带内计算 Heff 的值
                    for q_idx = 1:length(q_indices)
                        q = q_indices(q_idx);
                        
                        % --- 以下是与之前相同的计算逻辑 ---
                        phase_term = exp(1j * 2*pi/N * (N*c1*li^2 - q*li + N*c2*(q^2 - p^2)));
                        
                        F_i_pq_arg = (p - q + nu_i + 2*N*c1*li);
                        
                        % 环绕校正，确保 F_i_pq_arg 在 [-N/2, N/2] 附近以获得精确的闭式解
                        % 或者直接求和，更稳健
                        n_sum = (0:N-1)';
                        F_i_pq = sum(exp(-1j * 2*pi/N * F_i_pq_arg * n_sum));
                        
                        Hi_pq = (1/N) * phase_term * F_i_pq;
                        
                        Heff(p+1, q+1) = Heff(p+1, q+1) + h_est(i, col) * Hi_pq;
                    end
                end
            end
            
            % --- 均衡 ---
            H_eff = Heff(data_positon, data_positon);
            snr_linear = 10^(SNR(jj)/10); 
            W_LMMSE = (H_eff' * H_eff + (1/snr_linear) * eye(length(data_positon))) \ H_eff';
            y_sym = daft_data(data_positon, col);
            x_hat = W_LMMSE * y_sym;
            x_hat(~isfinite(x_hat)) = 0;
            data_aftereq(:, col) = x_hat;
        end
    end

        %% 并串转换
        data_aftereq = reshape(data_aftereq, [], 1);
        data_aftereq = data_aftereq(1:length(qpsk_data)); 

        %% QPSK解调
        data_aftereq = data_aftereq(1:length(qpsk_data));
        demodulation_data = pskdemod(data_aftereq, M, pi/M);
        De_data1 = reshape(demodulation_data, [], 1);
        De_data2 = de2bi(De_data1);

        De_Bit = reshape(De_data2', 1, []);

        %% 计算误码率
        [err, Ber] = biterr(De_Bit, P_data);
        Ber_record(jj) = Ber_record(jj) + Ber;
    end
         toc
    




end
   


%% 计算平均BER并绘图
Ber_record = Ber_record/N_simulation;
semilogy(SNR, Ber_record,'O-', 'DisplayName', 'AFDM Estimated CSI', 'LineWidth', 1.5);
grid on;
title('AFDM系统误码率 (分数多普勒)');
xlabel('SNR(dB)');
ylabel('BER');
xlim([min(SNR) max(SNR)]);
set(gca, 'YScale', 'log');
legend;

%绘制信道矩阵结构
figure;
imagesc(abs(Heff'));
colorbar;
title('估计的信道矩阵 |Heff^T|');
xlabel('p');
ylabel('q');
axis square;

function [h_pilot_resp] = calculate_pilot_response(l, nu, N, c1, c2, ind_r_E)
% CALCULATE_PILOT_RESPONSE - 计算给定信道参数下导频的理论响应向量
% 该函数遵循Bemani论文中的信道模型，用于信道估计
%
% 输入:
%   l          - 时延抽头 (integer)
%   nu         - 归一化多普勒频移 (alpha + a)
%   N          - 子载波数
%   c1, c2     - DAFT参数
%   ind_r_E    - 接收端用于估计的导频符号索引 (1-based)
%
% 输出:
%   h_pilot_resp - 导频响应向量 (列向量)

% 根据论文(74)式，导频响应是H_i,E矩阵的第一列 [cite: 973]
% H_i,E是H_i的第一列，并由ind_r_E选择行
% H_i的第一列表明 q=0

p = (0:N-1)'; % 接收符号索引 (0-based)
q = 0;        % 发射导频符号索引 (0-based)

% 根据(30)式计算相位项 [cite: 660]
phase_term = exp(1j * 2*pi/N * (N*c1*l^2 - q*l + N*c2*(q^2 - p.^2)));

% 根据(31)式计算F_i(p,q)项 [cite: 661]
n_sum = (0:N-1)';
F_i_pq_arg = p - q + nu + 2*N*c1*l;
F_i_pq = sum(exp(-1j * 2*pi/N * F_i_pq_arg .* n_sum), 1)'; % 逐元素计算后求和

% 根据(30)式计算H_i的第一列 [cite: 660]
h_i_col1 = (1/N) * phase_term .* F_i_pq;

% 根据ind_r_E提取用于估计的行
h_pilot_resp = h_i_col1(ind_r_E);

end