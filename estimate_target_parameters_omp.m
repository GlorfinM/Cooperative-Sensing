function [estimated_distance, estimated_velocity, estimated_angle_deg] = estimate_target_parameters_omp(Y_received, params, K)
% 功能: 使用OMP算法根据接收信号估计目标的距离、速度和角度。
%
% 输入参数:
%   Y_received : (Nc x Ns) 矩阵, 接收到的复数符号。
%   params     : 包含所有系统参数的结构体，必须包含以下字段:
%                fc, c0, lambda, Nt, Nr, d, Nc, Ns, fdelta, Ts, v_tx,
%                theta_start_deg, theta_end_deg
%   K          : 要检测的目标数量 (通常设为1)。
%
% 输出参数:
%   estimated_distance   : 估计的目标距离 (m)。
%   estimated_velocity   : 估计的目标径向速度 (m/s)。
%   estimated_angle_deg  : 估计的目标角度 (度)。


% 从结构体中提取参数
fc = params.fc; 
c0 = params.c0; 
lambda = params.lambda;
Nt = params.Nt; 
Nr = params.Nr; 
d = params.d;
Nc = params.Nc; 
Ns = params.Ns; 
fdelta = params.fdelta;
Ts = params.Ts; 
v_tx = params.v_tx;
theta_start_deg = params.theta_start_deg;
theta_end_deg = params.theta_end_deg;

% --- 定义参数的搜索网格 ---
num_dist_points = 50;
num_vel_points = 31;
num_angle_points = 61;

dist_grid = linspace(200, 300, num_dist_points);
vel_grid = linspace(-15, 15, num_vel_points);
angle_grid_deg = linspace(-30, 30, num_angle_points);
angle_grid_rad = deg2rad(angle_grid_deg);


y_omp = Y_received(:);

% --- 预计算 ---
W_rainbow = zeros(Nr, Nc);
for m_idx = 1:Nc
    W_rainbow(:, m_idx) = generate_rainbow_beam(Nr, fc, fdelta, Nc, d, theta_start_deg, theta_end_deg, m_idx);
end

A_tx_all = zeros(Nt, num_angle_points);
A_rx_all = zeros(Nr, num_angle_points);
for i = 1:num_angle_points
    theta = angle_grid_rad(i);
    A_tx_all(:, i) = exp(1j * 2 * pi * d * (0:Nt-1)' * sin(theta) / lambda);
    A_rx_all(:, i) = exp(1j * 2 * pi * d * (0:Nr-1)' * sin(theta) / lambda);
end

% 将预计算数据打包，传递给原子生成函数
precomputed.W_rainbow = W_rainbow;
precomputed.A_tx_all = A_tx_all;
precomputed.A_rx_all = A_rx_all;
precomputed.beam_gain_cache = cell(num_angle_points, 1); 


tic;
residual = y_omp;
index_set = [];
A_sub = zeros(Nc * Ns, K); % 预分配小的子字典内存

for iter = 1:K
    % --- 在线寻找最佳原子 ---
    max_projection = -1;
    best_atom_index = -1;
    
    atom_idx_1d = 1;
    for i_angle = 1:num_angle_points
        for i_dist = 1:num_dist_points
            for i_vel = 1:num_vel_points
                % 在线生成当前原子
                atom = generate_atom(i_angle, i_dist, i_vel, params, precomputed, ...
                                    angle_grid_rad, dist_grid, vel_grid);
                
                % 计算投影（相关性）
                current_projection = abs(atom' * residual);
                
                % 更新最大值
                if current_projection > max_projection
                    max_projection = current_projection;
                    best_atom_index = atom_idx_1d;
                end
                atom_idx_1d = atom_idx_1d + 1;
            end
        end
    end
    
    % --- 更新残差 ---
    if best_atom_index == -1
        break; % 如果找不到任何相关的原子，则退出
    end
    
    index_set = [index_set, best_atom_index];
    
    % 重新生成本次迭代找到的最佳原子，并加入子字典
    [i_a, i_d, i_v] = ind2sub([num_angle_points, num_dist_points, num_vel_points], best_atom_index);
    best_atom = generate_atom(i_a, i_d, i_v, params, precomputed, ...
                                angle_grid_rad, dist_grid, vel_grid);
    A_sub(:, iter) = best_atom;
    
    % 最小二乘法求解并更新残差
    residual = y_omp - A_sub(:, 1:iter) * (A_sub(:, 1:iter) \ y_omp);
end
toc;

if isempty(index_set)
    estimated_distance = NaN;
    estimated_velocity = NaN;
    estimated_angle_deg = NaN;
    return;
end

best_atom_index = index_set(1); % 以第一个找到的目标为例
[i_angle, i_dist, i_vel] = ind2sub([num_angle_points, num_dist_points, num_vel_points], best_atom_index);

estimated_angle_deg = angle_grid_deg(i_angle);
estimated_distance = dist_grid(i_dist);
estimated_velocity = vel_grid(i_vel);

end


% ================== 辅助函数：在线生成单个原子 ==================
function atom = generate_atom(i_angle, i_dist, i_vel, params, precomputed, angle_grid_rad, dist_grid, vel_grid)
    % 提取参数
    fc = params.fc; c0 = params.c0; lambda = params.lambda;
    Nc = params.Nc; Ns = params.Ns; fdelta = params.fdelta;
    Ts = params.Ts; v_tx = params.v_tx;

    % 获取当前参数值
    current_angle_rad = angle_grid_rad(i_angle);
    current_dist = dist_grid(i_dist);
    current_vel = vel_grid(i_vel);
    
    % 计算时延和多普勒
    tau = 2 * current_dist / c0;
    f_doppler = 2 * current_vel / lambda;
    
    % --- 计算角度相关的波束增益 ---
    % 检查缓存，如果已计算则直接使用，否则计算并存入缓存
    if isempty(precomputed.beam_gain_cache{i_angle})
        a_tx = precomputed.A_tx_all(:, i_angle);
        a_rx = precomputed.A_rx_all(:, i_angle);
        
        gain_per_carrier = zeros(Nc, 1);
        for m_idx = 1:Nc
            wm = precomputed.W_rainbow(:, m_idx);
            gain_per_carrier(m_idx) = (wm' * a_rx) * (a_tx' * v_tx);
        end
        precomputed.beam_gain_cache{i_angle} = gain_per_carrier;
    else
        gain_per_carrier = precomputed.beam_gain_cache{i_angle};
    end
    
    % --- 构建原子向量 ---
    atom = zeros(Nc * Ns, 1);
    for ns_idx = 1:Ns
        t_ns = (ns_idx - 1) * Ts;
        phase_vel = exp(1j * 2 * pi * f_doppler * t_ns);
        
        for m_idx = 1:Nc
            f_m = fc + (m_idx - 1) * fdelta;
            phase_dist = exp(-1j * 2 * pi * f_m * tau);
            
            linear_idx = (ns_idx - 1) * Nc + m_idx;
            atom(linear_idx) = gain_per_carrier(m_idx) * phase_dist * phase_vel;
        end
    end
end