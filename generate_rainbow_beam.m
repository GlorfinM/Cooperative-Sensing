function w_rainbow = generate_rainbow_beam(N, fc, fdelta, Nc, d_ant, theta_start_deg, theta_end_deg, m_idx)
% 功能: 根据可控波束偏斜(CBS-Beamforming)原理生成单个子载波的“彩虹波束”波束成形向量。
%       该向量随子载波索引 m_idx 变化，使不同子载波的波束指向不同方向。
%
% 输入参数:
%   N               : 天线阵列中的天线数量 (例如 Nt 或 Nr)
%   fc              : 系统中心载频 (Hz)
%   fdelta          : 子载波间隔 (Hz)
%   Nc              : 子载波总数
%   d_ant           : 天线阵列的天线间距 (m) (例如 d_tx 或 d_rx)
%   theta_start_deg : CBS波束起始扫描角度 (度)
%   theta_end_deg   : CBS波束终止扫描角度 (度)
%   m_idx           : 目标子载波的索引 (1-based, 即从1到Nc)
%
% 输出参数:
%   w_rainbow       : (N_ant x 1) 的列向量, 计算得到的波束成形向量

% 定义物理常数
c0 = 3e8; % 光速 m/s

% 将角度从度转换为弧度
theta_start_rad = deg2rad(theta_start_deg);
theta_end_rad = deg2rad(theta_end_deg);

% 计算rainboe-Beamforming所需的频率参数
W_bandwidth = (Nc-1) * fdelta;      % 系统总带宽 (Hz)
f0_lowest_freq = fc - W_bandwidth/2;    % 系统最低频率 (Hz)

% --- 计算指定子载波 m_idx 对应的波束成形向量 w ---

% 获取当前子载波的基带频率 (相对于带宽起始)
m_paper_idx = m_idx - 1; % 子载波的0-indexed索引 (0 to Nc-1)
tilde_fm = m_paper_idx * fdelta; 

% 初始化波束成形向量
w_rainbow = zeros(N, 1);

% 遍历天线计算向量的每个元素
for p_ant = 1:N
    % 中心化天线索引，与论文 n = -(N-1)/2, ..., (N-1)/2 的形式对应
    n_antenna_index = p_ant - (N+1)/2; 

    % 根据论文公式 (11) 及后续推导计算 PS 和 TTD 的等效作用
    % PS (移相器) 的设置使得最低频 f0 指向 theta_start
    phi_p_val = - (f0_lowest_freq * n_antenna_index * d_ant * sin(theta_start_rad)) / c0;
    
    % TTD (真延时线) 的设置使得最高频 fM (f0+W) 指向 theta_end
    t_p_val   = - (phi_p_val / W_bandwidth) - ...
                (((f0_lowest_freq + W_bandwidth) * n_antenna_index * d_ant * sin(theta_end_rad)) / (W_bandwidth * c0));
    
    % 根据论文公式 (9) 计算 w 的第 p_ant 个元素
    % [\tilde{w}]_n = (1/sqrt(N)) * exp(-j*2*pi*phi_n) * exp(-j*2*pi*tilde_f*t_n)
    w_rainbow(p_ant) = (1/sqrt(N)) * exp(-1j * 2 * pi * phi_p_val) * exp(-1j * 2 * pi * tilde_fm * t_p_val);
end

end