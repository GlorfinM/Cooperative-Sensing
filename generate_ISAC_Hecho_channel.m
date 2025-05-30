function H_echo = generate_ISAC_Hecho_channel(Nt, Nr, Nc, Ns, fc, lambda, d_tx, d_rx, Ts, target_pos_xy, target_vel_xy, target_rcs, fdelta, rx_pos_xy, clutter_info)
% 功能: 根据YOLO论文模型及双基站配置生成动态无人机目标的ISAC回波信道矩阵，并包含静态杂波。
%       主要参考论文中的公式 (4) 和 (5) (动态目标部分)，并适配双基站模型。
%       发射机固定在原点 [0,0]。
%
% 输入参数:
%   Nt            : 发射天线数量
%   Nr            : 接收天线数量
%   Nc            : 子载波数量 (对应论文中的 M+1)
%   Ns            : OFDM符号数量
%   fc            : 中心载频 (Hz) (论文中 f_c, 多普勒计算中的 f_0 近似为此)
%   lambda        : 波长 (m)
%   d_tx          : 发射天线间距 (m)
%   d_rx          : 接收天线间距 (m)
%   Ts            : OFDM符号周期 (s)
%   target_pos_xy : 动态无人机位置 [x, y] (m), 行向量或列向量
%   target_vel_xy : 动态无人机速度 [vx, vy] (m/s), 行向量或列向量
%   target_rcs    : 动态无人机的雷达散射截面 (m^2) (论文中 sigma_c,k)
%   fdelta        : 子载波间隔 (Hz)
%   rx_pos_xy     : 接收基站位置 [x, y] (m)
%   clutter_info  : 静态杂波信息, cell数组，每个cell是一个struct，包含:
%                   clutter_info{i}.pos_xy = [x,y] % 杂波位置
%                   clutter_info{i}.rcs    = value % 杂波RCS
%
% 输出参数:
%   H_echo        : Ns x Nc 的 cell 数组, 每个 cell 是 Nr x Nt 的回波信道矩阵 
%                   (H_dynamic_target + Sum_H_clutter)

% 定义物理常数
c0 = 3e8;               % 光速 m/s
pos_tx = [0,0];         % 发射机位置固定在原点

% 初始化输出
H_echo = cell(Ns, Nc); % cell索引从1开始, 每个元素是 Nr x Nt

% 子载波频率计算
f_start = fc - (Nc-1)/2 * fdelta;

% -------------------- 动态目标参数计算 --------------------
% 1. Tx到目标的几何参数
vec_tx_to_target = target_pos_xy - pos_tx;
r_tx = norm(vec_tx_to_target);
if r_tx < 1e-6 % 避免除以零
    r_tx = 1e-6;
end
theta_tx = atan2(vec_tx_to_target(2), vec_tx_to_target(1)); % Tx到目标的角度 (弧度)

% 2. Rx到目标的几何参数
vec_rx_to_target = target_pos_xy - rx_pos_xy;
r_rx = norm(vec_rx_to_target);
if r_rx < 1e-6 % 避免除以零
    r_rx = 1e-6;
end
% theta_rx 是指信号从目标到达接收机的角度 (AoA at Rx)
% 导向矢量通常定义为信号到达方向，所以是从目标看向接收机
% 接收机本地坐标系与全局坐标系平行，为从目标指向接收机的向量与x轴的夹角
vec_target_to_rx = rx_pos_xy - target_pos_xy;
theta_rx = atan2(vec_target_to_rx(2), vec_target_to_rx(1)); % 目标到Rx的角度 (用于构造Rx导向矢量)

% 3. 计算路径损耗和RCS相关的系数 alpha_k (双基站)
% 论文中 alpha_k proportional to 1/r_k^2 (monostatic)
% 双基站中 alpha_k proportional to 1/(r_tx * r_rx)
alpha_k_dynamic = sqrt(lambda^2 * target_rcs / ((4*pi)^3)) / (r_tx * r_rx);

% 4. 计算双基站多普勒频移相关参数
% v_comp_tx: 目标速度在 Tx-目标 连线上的投影分量 (远离Tx为正)
v_comp_tx = dot(target_vel_xy, vec_tx_to_target / r_tx);
% v_comp_rx: 目标速度在 目标-Rx 连线上的投影分量 (远离目标朝向Rx为正)
v_comp_rx = dot(target_vel_xy, vec_target_to_rx / r_rx); 
% 有效速度分量之和，影响多普勒频移。论文中为 2*v_k (v_k为单基地径向速度)
% 这里是两个路径上速度投影的和 (v_path_tx + v_path_target_rx)
effective_doppler_vel = v_comp_tx + v_comp_rx;


% --- 迭代生成每个OFDM符号和每个子载波的信道矩阵 ---
for ns_idx = 1:Ns % OFDM 符号索引
    n_s_paper = ns_idx - 1; % 论文中的 n_s 从0开始

    for m_idx = 1:Nc % 子载波索引
        m_paper = m_idx - 1; % 论文中的 m 从0开始
        fm = f_start + m_paper * fdelta; % 当前子载波频率 f_m

        % --- (A) 计算动态目标的回波信道 H_k,ns,m ---
        % A1. 计算发射导向矢量 a_k,m_tx (Nt x 1)
        a_km_tx = zeros(Nt, 1);
        for ant_idx = 1:Nt
            % 天线索引 n,使其中心为0 (ULA沿y轴)
            % 天线索引方式围绕中心对称
            n_ant_tx = ant_idx - (Nt+1)/2; 
            phase_val_tx = (2*pi*fm/c0) * n_ant_tx * d_tx * sin(theta_tx);
            a_km_tx(ant_idx) = exp(1j * phase_val_tx);
        end

        % A2. 计算接收导向矢量 a_k,m_rx (Nr x 1)
        a_km_rx = zeros(Nr, 1);
        for ant_idx = 1:Nr
            n_ant_rx = ant_idx - (Nr+1)/2;
            phase_val_rx = (2*pi*fm/c0) * n_ant_rx * d_rx * sin(theta_rx); % 使用 theta_rx
            a_km_rx(ant_idx) = exp(1j * phase_val_rx);
        end
        
        % A3. 计算动态目标的相位项
        % 多普勒相位 (使用有效多普勒速度)
        doppler_phase_dynamic = exp(1j * 2*pi * fc * (effective_doppler_vel/c0) * n_s_paper * Ts);
        % 距离相位 (使用双基站路径长度)
        range_phase_dynamic = exp(-1j * 2*pi * fm * (r_tx + r_rx)/c0);
        
        % A4. 计算动态目标信道矩阵 (Nr x Nt)
        H_dynamic_target_ns_m = alpha_k_dynamic * doppler_phase_dynamic * range_phase_dynamic * (a_km_rx * a_km_tx.');
        
        % --- (B) 计算静态环境杂波的回波信道 Sum_i H_i,ns,m ---
        H_static_clutter_ns_m = zeros(Nr, Nt); % 初始化当前符号和子载波的总杂波信道
        
        if ~isempty(clutter_info) && iscell(clutter_info)
            for i_clutter = 1:length(clutter_info)
                clutter_pos = clutter_info{i_clutter}.pos_xy;
                clutter_rcs_val = clutter_info{i_clutter}.rcs;

                % B1. Tx到杂波的几何参数
                vec_tx_to_clutter = clutter_pos - pos_tx;
                r_tx_c = norm(vec_tx_to_clutter);
                if r_tx_c < 1e-6, r_tx_c = 1e-6; end
                theta_tx_c = atan2(vec_tx_to_clutter(2), vec_tx_to_clutter(1));

                % B2. Rx到杂波的几何参数
                vec_rx_to_clutter = clutter_pos - rx_pos_xy;
                r_rx_c = norm(vec_rx_to_clutter);
                if r_rx_c < 1e-6, r_rx_c = 1e-6; end
                vec_clutter_to_rx = rx_pos_xy - clutter_pos;
                theta_rx_c = atan2(vec_clutter_to_rx(2), vec_clutter_to_rx(1));
                
                % B3. 杂波路径损耗系数 beta_i
                beta_i_clutter = sqrt(lambda^2 * clutter_rcs_val / ((4*pi)^3)) / (r_tx_c * r_rx_c);
                
                % B4. 计算杂波的发射和接收导向矢量
                a_cm_tx = zeros(Nt, 1); % Clutter transmit steering vector
                for ant_idx = 1:Nt
                    n_ant_tx = ant_idx - (Nt+1)/2;
                    phase_val_tx_c = (2*pi*fm/c0) * n_ant_tx * d_tx * sin(theta_tx_c);
                    a_cm_tx(ant_idx) = exp(1j * phase_val_tx_c);
                end

                a_cm_rx = zeros(Nr, 1); % Clutter receive steering vector
                for ant_idx = 1:Nr
                    n_ant_rx = ant_idx - (Nr+1)/2;
                    phase_val_rx_c = (2*pi*fm/c0) * n_ant_rx * d_rx * sin(theta_rx_c);
                    a_cm_rx(ant_idx) = exp(1j * phase_val_rx_c);
                end
                
                % B5. 静态杂波的距离相位 (多普勒相位为1，因为是静态的)
                range_phase_clutter = exp(-1j * 2*pi * fm * (r_tx_c + r_rx_c)/c0);
                
                % B6. 单个杂波单元的信道矩阵 (Nr x Nt)
                H_clutter_i_ns_m = beta_i_clutter * range_phase_clutter * (a_cm_rx * a_cm_tx.');
                
                H_static_clutter_ns_m = H_static_clutter_ns_m + H_clutter_i_ns_m;
            end
        end
        
        % --- (C) 总回波信道 ---
        H_echo{ns_idx, m_idx} = H_dynamic_target_ns_m + H_static_clutter_ns_m;
    end
end

end