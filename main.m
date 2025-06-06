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
Nc = 3276 / 12;              % 子载波数量 (共Nc个)
M = Nc-1;               % 子载波索引
fdelta = 120e3 * 12;         % 带宽 120 kHz (frequency interval)
Ts = 35.6e-06;          % OFDM符号周期 35.6μs, 包含Tcp
Ns = 70;                % OFDM符号数 ->1?

% 发射功率参数
Ptx_dBm = 44;           % 发射功率 44 dBm
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

% Rainbow-Beamforming 相关参数 (用于计算随子载波变化的接收波束成形向量 wm)
theta_start_deg = 60;       % rainbow-beam波束起始扫描角度 (度)
theta_end_deg = -60;      % rainbow-beam波束终止扫描角度 (度)

% 计算Rainbow-Beamforming所需的频率参数
W_bandwidth = (Nc-1) * fdelta;              % 系统总带宽 (Hz)
f0_lowest_freq = fc - W_bandwidth/2;        % 系统最低频率 (Hz)


%% ================== 信号生成 ==================

% 发射信号生成
% phi = pi/4 + randi([0 3], Nc, Ns) * pi/2;
% s = exp(1i * phi); % QPSK symbols: s(Nc, Ns)
s_data = ones(Nc, Ns); % QPSK symbols: s(Nc, Ns)
% --- (a) 对每个 OFDM 符号进行 IFFT ---
% ifft(s, Nc, 1) 表示沿着第一个维度 (子载波维度) 进行 Nc 点 IFFT
% s_time_parallel = ifft(s, Nc, 1); 
% --- (b) 添加循环前缀 (CP) ---
% 从每个符号末尾取 Ncp 个样本加到开头
% --- (c) 并行转串行 (P/S) ---

% 将所有符号按时间顺序连接成一个长的列向量: (Nc*Ns,1)
s = s_data(:);

% 调用 AMCF_ZCI 函数，根据指定参数计算波束成形权重。
v = AMCF_ZCI(Nt, Omega0, B, Q, M0).';
% x = s * v;

% 所有子载波的接收波束成形向量 (wm)
w = cell(Nc, 1);

for m_idx = 1:Nc
    w{m_idx} = generate_rainbow_beam(Nr, fc, fdelta, Nc, d, theta_start_deg, theta_end_deg, m_idx);
end

%% ================== 信道建模 ==================
% 生成信道矩阵
H_channel_rx1 = generate_ISAC_Hecho_channel(Nt, Nr, Nc, Ns, fc, lambda, d, d, Ts, target_pos, target_vel, rcs, fdelta, pos_rx1);
H_channel_rx2 = generate_ISAC_Hecho_channel(Nt, Nr, Nc, Ns, fc, lambda, d, d, Ts, target_pos, target_vel, rcs, fdelta, pos_rx2);

% 生成接收信号
Y_received_rx1 = generate_signal_Yreceived(Nc, Ns, Nt, Nr, fc, fdelta, d, s_data, v.', H_channel_rx1, theta_start_deg, theta_end_deg);
Y_received_rx2 = generate_signal_Yreceived(Nc, Ns, Nt, Nr, fc, fdelta, d, s_data, v.', H_channel_rx2, theta_start_deg, theta_end_deg);


%% ================== 感知算法 ==================

% 将所有需要的参数打包到一个结构体中
params.fc = fc;
params.c0 = c0;
params.lambda = lambda;
params.Nt = Nt;
params.Nr = Nr;
params.d = d;
params.Nc = Nc;
params.Ns = Ns;
params.fdelta = fdelta;
params.Ts = Ts;
params.v_tx = v.';
params.theta_start_deg = theta_start_deg;
params.theta_end_deg = theta_end_deg;

% 设置要寻找的目标数量
K = 1;

% OMP估计
[est_dist1, est_vel1, est_angle1] = estimate_target_parameters_omp(Y_received_rx1, params, K);
[est_dist2, est_vel2, est_angle2] = estimate_target_parameters_omp(Y_received_rx2, params, K);

%% ================== 性能评估 ==================

% % 计算RMSE (均方根误差)
% pos_error = abs(dist_est - norm(target_pos)); 
% vel_error = abs(vel_est - norm(target_vel * (target_pos' / norm(target_pos)))); 
% 
% fprintf('\n--- 性能评估 ---\n');
% fprintf('真实距离: %.2f m\n', norm(target_pos));
% fprintf('估计距离误差: %.2f m\n', pos_error);
% fprintf('真实径向速度: %.2f m/s\n', norm(target_vel * (target_pos' / norm(target_pos))));
% fprintf('估计速度误差: %.2f m/s\n', vel_error);