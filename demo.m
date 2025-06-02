%% --- 测试用例 ---
clear;
clc;


% 基本物理参数
fc = 26e9;              % 载频 26 GHz
c0 = 3e8;               % 光速 m/s
lambda = c0 / fc;       % 波长

% 天线阵列参数
Nt = 16;                % 发射天线数
Nr = 16;                % 接收天线数
d = lambda / 2;         % 天线间距

N = 16;              % 天线数
Omega0 = -0.5;       % 波束起始角度
B = 1;             % 波束宽度
Q = 512;             % 角度量化数
M = 1000;              % 最大迭代次数
v = AMCF_ZCI(N, Omega0, B, Q, M);


% % 可视化波束图
% Omega_plot = linspace(-1, 1, 1000);
% beam_pattern = zeros(size(Omega_plot));
% for i = 1:length(Omega_plot)
%     a = exp(1j * pi * (0:N-1)' * Omega_plot(i)) / sqrt(N);
%     beam_pattern(i) = abs(a' * v);
% end
% figure;
% plot(Omega_plot, 20*log10(beam_pattern / max(beam_pattern)))
% xlabel('\Omega'); ylabel('Beam Gain (dB)');
% title('AMCF-ZCI Wide Beam Pattern'); grid on;

% ================== 计算和绘制发射波束方向图 ==================
% 定义角度扫描范围 (方位角，假设阵列位于x轴)
theta = linspace(-pi/2, pi/2, 360); % 从 -90 度到 90 度，360个采样点

% 计算波数 k
k = 2 * pi / lambda;

% 初始化阵列因子
AF = zeros(size(theta));

% 计算阵列因子
for i = 1:length(theta)
    steering_vector = exp(1j * k * d * (0:Nt-1)' * sin(theta(i))); % 计算当前角度的导向矢量
    AF(i) = abs(v.' * steering_vector); % v 是 1xNt, steering_vector 是 Ntx1 
end

% 计算波束能量 (功率方向图)，并转换为 dB
power_dB = 10 * log10(abs(AF).^2 / max(abs(AF).^2)); % 归一化到最大值

% 绘制波束方向图
figure;
plot(theta * 180 / pi, power_dB, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Angle (degrees)');
ylabel('Normalized Beam Power (dB)');
title('Transmit Beam Pattern');
ylim([-40 0]); % 设置y轴范围以便观察旁瓣