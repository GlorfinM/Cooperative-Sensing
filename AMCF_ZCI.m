function v = AMCF_ZCI(N, Omega0, B, Q, M)
% 输入参数:
%   N      : 天线数量 (对应论文中的 N_UE)
%   Omega0 : 波束覆盖的起始角度 (例如 -0.5, 对应归一化空间频率)
%   B      : 波束宽度 (例如 0.5, 对应归一化空间频率范围)
%   Q      : AoA (到达角) 的量化点数 (例如 512)
%   M      : 迭代次数 (例如 20)
% 输出参数:
%   v      : 设计得到的 N x 1 维波束赋形向量


% Step 1: 使用类 Zadoff-Chu (ZC) 序列进行初始化
% 对应算法1中的 "Obtain v(0) via (28)"
v = zeros(N, 1); 
for n_init = 1:N % 遍历每个天线单元进行初始化
    if mod(N, 2) == 0
        % ZC序列相位计算公式 (偶数情况)
        phase = pi * (B * n_init^2 / (2 * N) + n_init * Omega0);
        v(n_init) = (1/sqrt(N)) * exp(1j * phase);
    else 
        % ZC序列相位计算公式 (奇数情况)
        phase = pi * (B * n_init * (n_init + 1) / (2 * N) + n_init * Omega0);
        v(n_init) = (1/sqrt(N)) * exp(1j * phase);
    end
end
% 初始化后的 v 即为 v(0)


% Step 2: 构建导向矩阵 A 和期望波束增益向量 g
% Omega 是归一化的空间频率 (-1 到 1 之间: -pi -> pi)
Omega = linspace(-1, 1, Q); % 生成 Q 个均匀分布的 AoA 采样点
A = zeros(N, Q);            % 初始化导向矩阵 A (N x Q)
g = zeros(Q, 1);            % 初始化期望波束增益向量 g (Q x 1)

for q_idx = 1:Q % 遍历每个 AoA 量化点
    % 计算对应 Omega(q_idx) 的导向向量 (A 的列向量)
    % 导向向量的每个元素表示对应天线接收到的信号相位
    A(:, q_idx) = exp(1j * pi * (0:N-1)' * Omega(q_idx)) / sqrt(N);

    % 设置期望增益 g
    % 如果当前 AoA 在期望的波束覆盖范围内 [Omega0, Omega0 + B]
    if Omega(q_idx) >= Omega0 && Omega(q_idx) <= Omega0 + B
        g(q_idx) = sqrt(2 / B); % 目标波束内的期望增益
    else
        g(q_idx) = 0;           % 目标波束外的期望增益为 0 (抑制旁瓣)
    end
end

% % Step 3: Iterative optimization (AMCF)
% for m = 1:M
%     r = g .* exp(1j * angle(A' * v));
%     p = A * r;
% 
%     % Closed-form solution for constant modulus vector v
%     for n = 1:N
%         v(n) = p(n) / abs(p(n));
%     end
% end

% Step 3: 迭代优化 (AMCF - Alternating Minimization Constant modulus Codeword Framework)
% 此处 v 是步骤1中得到的 v(0)
% 循环迭代 M 次
for m = 1:M
    % 对应算法1的第5行: 计算 Theta_m = angle(A^H * v(m-1))
    % A 是 N x Q 矩阵, v 是 N x 1 向量. A' 是 Q x N 矩阵. A'*v 是 Q x 1 向量.
    Theta_m = angle(A' * v); % Theta_m 是 Q x 1 的相位向量

    % 对应算法1的第6行: 计算 r_m = g .* exp(1j * Theta_m)
    % g 是 Q x 1 向量. '.*' 表示元素级乘法.
    r_m = g .* exp(1j * Theta_m); % r_m 是 Q x 1 向量

    % 隐式步骤: 计算 p_m = A * r_m
    % 这个 p_m 用于构造公式 (26) 中的 t 向量分量
    % A 是 N x Q 矩阵, r_m 是 Q x 1 向量. p_m 是 N x 1 复数向量.
    p_m = A * r_m;

    % 对应算法1的第7行: 通过公式 (27) 和 (26) 计算 v(m)
    % 初始化 v_next, 它将成为当前的 v(m)
    v_next = zeros(N, 1); % N x 1 的零向量

    % 遍历每个天线单元 n (对应论文中的 n_idx 从 1 到 N_UE)
    for n_idx = 1:N
        % 获取 p_m(n_idx) 的实部和虚部
        % 根据定义 t = [Re{p}; Im{p}],
        % 对于 v 的第 n 个元素, [t]_n 是 Re{p_m(n_idx)}, [t]_{n+N_UE} 是 Im{p_m(n_idx)}
        t_n_val = real(p_m(n_idx));        % p_m(n_idx) 的实部
        t_n_plus_N_val = imag(p_m(n_idx)); % p_m(n_idx) 的虚部

        % 计算公式 (26) 中的各项
        % 分母项: sqrt(N_UE * ([t]_n^2 + [t]_{n+N_UE}^2))
        % 此处 N_UE 即为 N
        denominator_squared_term = N * (t_n_val^2 + t_n_plus_N_val^2);

        u_n_val = 0;        % 如果 p_m(n_idx) 为零, 则 u_n_val 默认为 0
        u_n_plus_N_val = 0; % 如果 p_m(n_idx) 为零, 则 u_n_plus_N_val 默认为 0

        % 检查 p_m(n_idx) 是否非零, 以避免除以零的错误
        % 使用一个小的 epsilon 进行浮点数比较
        if denominator_squared_term > 1e-12 % 判断分母项是否显著大于零的阈值
            denominator_val = sqrt(denominator_squared_term); % 计算分母
            % 根据公式 (26) 计算 u_n 和 u_{n+N_UE}
            u_n_val = t_n_val / denominator_val;
            u_n_plus_N_val = t_n_plus_N_val / denominator_val;
        end
        % 如果 denominator_squared_term 为零 (或数值上接近零),
        % 则 t_n_val 和 t_n_plus_N_val 也为零, 因此 u_n_val 和 u_n_plus_N_val 保持为 0.

        % 根据公式 (27): [v(m)]_n = [u]_n + j * [u]_{n+N_UE}
        % 构造 v(m) 的第 n_idx 个元素
        v_next(n_idx) = u_n_val + 1j * u_n_plus_N_val;
    end

    % 更新 v 为当前迭代计算得到的 v_next, 用于下一次迭代
    v = v_next;
end
% 输出: v 是经过 M 次迭代后的最终码字 v(M)
% 如果 p_m(n_idx) 非零, v 的每个元素的模长将是 1/sqrt(N);
% 如果 p_m(n_idx) 为零, 则 v 的对应元素为 0.
end