function Y_received = generate_signal_Yreceived(Nc, Ns, Nt, Nr, fc, fdelta, d_rx, s, v_tx, H_echo, theta_start_deg, theta_end_deg)
% 功能: 根据指定的发射和接收波束成形策略以及信道，计算接收信号矩阵。
%       - 发射波束成形向量 v_tx: 固定，与子载波和OFDM符号无关。
%       - 接收波束成形向量 wm: 随子载波 m 变化 (实现CBS-Beamforming)。
%
% 输入参数:
%   Nc              : 子载波数量
%   Ns              : OFDM符号数量
%   Nt              : 发射天线数
%   Nr              : 接收天线数
%   fc              : 系统中心载频 (Hz)
%   fdelta          : 子载波间隔 (Hz)
%   d_rx            : 接收天线阵列的天线间距 (m)
%   s               : (Nc x Ns) 矩阵, 原始符号
%   v_tx            : (Nt x 1) 列向量, 固定的发射波束成形向量
%   H_echo          : (Ns x Nc) cell数组, H_echo{ns_idx, m_idx} 是 (Nr x Nt) 信道矩阵
%   theta_start_deg : rainbow-beam波束起始扫描角度 (度)
%   theta_end_deg   : rainbow-beam波束终止扫描角度 (度)
%
% 输出参数:
%   Y_received      : (Nc x Ns) 矩阵, 接收到的复数符号

% 光速 (m/s)
c0 = 3e8;
% rainbow-Beamforming 相关参数转换为弧度
theta_start_rad = deg2rad(theta_start_deg);
theta_end_rad = deg2rad(theta_end_deg);

% rainbow-Beamforming所需的频率参数
W_bandwidth = (Nc-1) * fdelta;              % 系统总带宽 (Hz)
f0_lowest_freq = fc - W_bandwidth/2;        % 系统最低频率 (Hz)

% --- 初始化接收信号矩阵 ---
Y_received = zeros(Nc, Ns); % 存储接收到的复数符号 (维度与 s_qpsk 一致)

% --- 循环计算接收信号 ---
for ns_idx = 1:Ns  % 遍历OFDM符号
    for m_idx = 1:Nc % 遍历子载波
        
        % (a) 获取当前子载波和符号的原始QPSK符号 (标量)
        current_symbol_mn = s(m_idx, ns_idx); 
        
        % (b) 获取对应的信道矩阵 (Nr x Nt)
        current_channel_H_mn = H_echo{ns_idx, m_idx}; 
        
        % (c) 计算当前子载波 m_idx 对应的接收波束成形向量 wm (Nr x 1)
        m_paper_idx = m_idx - 1; % 子载波的0-indexed索引 (0 to Nc-1)
        tilde_fm = m_paper_idx * fdelta; % 当前子载波的基带频率 (相对于带宽起始)
        
        wm = zeros(Nr, 1); % 初始化当前子载波的接收波束成形向量
        for p_ant = 1:Nr % 遍历接收天线
            n_rx_antenna_index = p_ant - (Nr+1)/2; % 中心化天线索引

            % 根据论文公式 (11) 及后续推导计算 PS 和 TTD 的等效作用
            % PS 设置使 f0 指向 theta_start
            phi_p_val = - (f0_lowest_freq * n_rx_antenna_index * d_rx * sin(theta_start_rad)) / c0;
            % TTD 设置使 fM (最高频 f0+W) 指向 theta_end
            t_p_val   = - (phi_p_val / W_bandwidth) - ...
                        (((f0_lowest_freq + W_bandwidth) * n_rx_antenna_index * d_rx * sin(theta_end_rad)) / (W_bandwidth * c0));
            
            % 根据论文公式 (9) 计算 wm 的第p个元素
            wm(p_ant) = (1/sqrt(Nr)) * exp(-1j * 2 * pi * phi_p_val) * exp(-1j * 2 * pi * tilde_fm * t_p_val);
        end
        
        % (d) 计算经过发射波束成形后的信号向量 (Nt x 1)
        x_transmit_antenna_signal_mn = v_tx * current_symbol_mn; 
        
        % (e) 计算信号经过信道后在接收天线端的信号向量 (Nr x 1)
        r_receive_antenna_signal_mn = current_channel_H_mn * x_transmit_antenna_signal_mn;
        
        % (f) 计算经过随子载波变化的接收波束成形 wm 后的最终接收符号 (标量)
        % wm' 是 wm 的共轭转置
        Y_received(m_idx, ns_idx) = wm' * r_receive_antenna_signal_mn;
    end
end

end