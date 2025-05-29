function H_echo = generate_ISAC_Hecho_channel(Nt, Nc, Ns, fc, lambda, d, c0, Ts, target_pos_xy, target_vel_xy, rcs)
% 功能: 根据YOLO论文模型生成动态无人机目标的ISAC回波信道矩阵
%       主要参考论文中的公式 (4) 和 (5) (动态目标部分)
%
% 输入参数:
%   Nt            : 天线数量 (对应论文中的 N)
%   Nc            : 子载波数量 (对应论文中的 M+1)
%   Ns            : OFDM符号数量
%   fc            : 中心载频 (Hz) (论文中 f_c, 多普勒计算中的 f_0 近似为此)
%   lambda        : 波长 (m)
%   d             : 天线间距 (m)
%   c0            : 光速 (m/s)
%   Ts            : OFDM符号周期 (s)
%   target_pos_xy : 动态无人机位置 [x, y] (m), 行向量或列向量
%   target_vel_xy : 动态无人机速度 [vx, vy] (m/s), 行向量或列向量
%   rcs           : 无人机的雷达散射截面 (m^2) (论文中 sigma_c,k)
%
% 输出参数:
%   H_echo        : Ns x Nc 的 cell 数组, 每个 cell 是 Nt x Nt 的回波信道矩阵 H_k,ns,m

% 初始化输出
H_echo = cell(Ns, Nc); % MATLAB中cell索引从1开始

% 子载波频率计算
subcarrier_spacing = 120e3; % 从 main.m 获取 (fdelta)
f_start = fc - (Nc-1)/2 * subcarrier_spacing; % 假设的最低频率

% 目标参数转换
target_r = norm(target_pos_xy); % 距离 r_k
target_theta = atan2(target_pos_xy(2), target_pos_xy(1)); % 角度 theta_k (弧度)

% 计算径向速度 v_k
if target_r < 1e-6 % 避免除以零
    radial_vel = 0;
else
    radial_vel = dot(target_vel_xy, target_pos_xy) / target_r;
end

% 计算路径损耗和RCS相关的系数 alpha_k
alpha_k_base = sqrt(lambda^2 * rcs / ((4*pi)^3 * target_r^4));

% 迭代生成每个OFDM符号和每个子载波的信道矩阵
for ns_idx = 1:Ns % OFDM 符号索引 (对应论文 n_s, 从 0 到 Ns-1)
    n_s_paper = ns_idx - 1; % 论文中的 n_s 从0开始

    for m_idx = 1:Nc % 子载波索引 (对应论文 m, 从 0 到 M)
        m_paper = m_idx - 1; % 论文中的 m 从0开始

        fm = f_start + m_paper * subcarrier_spacing; % 当前子载波频率 f_m

        % 1. 计算导向矢量 a_k,m
        a_km = zeros(Nt, 1);
        for ant_idx = 1:Nt
            n_antenna_paper_index = ant_idx - (Nt+1)/2; % 天线索引 n，使其中心为0
            phase_val = (2*pi*fm/c0) * n_antenna_paper_index * d * sin(target_theta);
            a_km(ant_idx) = exp(1j * phase_val);
        end

        % 2. 计算动态目标的回波信道 H_k,ns,m (论文公式(4))
        doppler_phase_term = exp(1j * 2*pi * fc * (2*radial_vel/c0) * n_s_paper * Ts);
        range_phase_term = exp(-1j * 2*pi * fm * (2*target_r/c0));
        
        H_k_ns_m = alpha_k_base * doppler_phase_term * range_phase_term * (a_km * a_km.');
        
        H_echo{ns_idx, m_idx} = H_k_ns_m;
    end
end

end