% 主程序文件：main.m
% 功能：参数初始化和算法调用

clear;
close all;
clc;

%% ================== 基本参数初始化 ==================

% 基本物理参数
fc = 26e9;              % 载频 26 GHz
c0 = 3e8;               % 光速 m/s
lambda = c0 / fc;       % 波长

% 天线阵列参数
Nt = 16;                % 发射天线数
Nr = 16;                % 接收天线数
d = lambda / 2;         % 天线间距

% 波束成形参数
Omega0 = -0.5;          % 波束起始角度 (归一化)
B = 1;                  % 波束宽度 (归一化)
Q = 512;                % 角度量化数
M0 = 1000;              % 最大迭代次数

% OFDM参数
Nc = 3276;              % 子载波数量 (共Nc个)
M = Nc-1;               % 子载波索引
fdelta = 120e3;         % 带宽 120 kHz (frequency interval)
Ts = 35.6e-06;          % OFDM符号周期 35.6μs, 包含Tcp
Ns = 70;                % OFDM符号数 ->1?
Ncp = 1;

% 发射功率参数
Ptx_dBm = 45;           % 发射功率 45 dBm
Ptx = 10^((Ptx_dBm - 30) / 10);  % 转换为瓦特

% 基站几何配置 (等边三角形，边长500m)
L = 500;                % 基站间距离 500m
pos_tx = [0, 0];        % 发射基站位置
pos_rx1 = [L*sqrt(3)/2,L/2];   % 接收基站1位置
pos_rx2 = [L*sqrt(3)/2,-L/2];  % 接收基站2位置

% 目标参数
target_pos = [250, 50];      % 目标位置 [x, y] (m)
target_vel = [10, 5];         % 目标速度 [vx, vy] (m/s)  
rcs = 0.01;                   % 雷达散射截面 0.01 m²

% 噪声参数
NF_dB = 6;              % 噪声系数 6 dB
T0 = 290;               % 标准温度 290K
kB = 1.38e-23;          % 玻尔兹曼常数


%% ================== 信号生成 ==================

% 发射信号生成

phi = pi/4 + randi([0 3], Nc, Ns) * pi/2;
s = exp(1i * phi); % QPSK symbols: s(Nc, Ns)
% --- (a) 对每个 OFDM 符号进行 IFFT ---
% ifft(s, Nc, 1) 表示沿着第一个维度 (子载波维度) 进行 Nc 点 IFFT
s_time_parallel = ifft(s, Nc, 1); 
% --- (b) 添加循环前缀 (CP) ---
% 从每个符号末尾取 Ncp 个样本加到开头
% --- (c) 并行转串行 (P/S) ---
% 将所有符号按时间顺序连接成一个长的列向量: (Nc*Ns,1)
s_time = s_time_parallel(:);

% 调用 AMCF_ZCI 函数，根据指定参数计算波束成形权重。
v = AMCF_ZCI(Nt, Omega0, B, Q, M0).';

x = s_time * v;



% 生成接收信号
% [y_rx1, y_rx2] = generate_received_signals();

%% ================== 信道建模 ==================

% 子载波频率
freq_sub = fc + (0:Nc) * fdelta / Nc;

% 信道系数计算
% alpha = sqrt(lambda^2 * rcs / ((4*pi)^3 * range_tx^4));

% 生成信道矩阵
H_channel = generate_ISAC_Hecho_channel(Nt, Nc, fc, lambda, d, Ts, target_pos, target_vel, rcs, fdelta);

%% ================== YOLO==================


%% ================== 性能评估 ==================

% 计算RMSE
% pos_error = norm(pos_est - target_pos);
% vel_error = norm(vel_est - target_vel);

%% ================== 函数定义 ==================

function [y_rx1, y_rx2] = generate_received_signals(s_symbols, H_channel, w_rx, ...
                                                   v_tx, M, Ns, noise_power)
    % 生成接收信号
    y_rx1 = zeros(Ns, M + 1);
    y_rx2 = zeros(Ns, M + 1);
    
    for m = 0:M
        for ns = 1:Ns
            H_cells = H_channel{m + 1, ns};
            H1 = H_cells{1};
            H2 = H_cells{2};
            
            % 接收信号计算
            signal1 = w_rx{m + 1}' * H1 * v_tx * s_symbols(m + 1);
            signal2 = w_rx{m + 1}' * H2 * v_tx * s_symbols(m + 1);
            
            % 添加噪声
            noise1 = sqrt(noise_power/2) * (randn + 1j*randn);
            noise2 = sqrt(noise_power/2) * (randn + 1j*randn);
            
            y_rx1(ns, m + 1) = signal1 + noise1;
            y_rx2(ns, m + 1) = signal2 + noise2;
        end
    end
end