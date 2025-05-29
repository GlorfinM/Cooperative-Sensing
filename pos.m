clear;
clc;
close all;

% --- 参数定义 ---
% 基站几何配置 (等边三角形，边长500m)
L = 500;                % 基站间距离 500m
pos_tx = [0, 0];        % 发射基站位置 [x, y]
pos_rx1 = [L*sqrt(3)/2, L/2];   % 接收基站1位置 [x, y]
pos_rx2 = [L*sqrt(3)/2, -L/2];  % 接收基站2位置 [x, y]

% 目标参数
target_pos = [250, 50];      % 目标位置 [x, y] (m)

% --- 绘图 ---
% 创建一个新的图形窗口
figure;
hold on; % 允许在同一图形上叠加绘制多个对象

% 绘制发射基站 (Tx) - 使用蓝色方块表示
plot(pos_tx(1), pos_tx(2), 'bs', 'MarkerSize', 12, 'MarkerFaceColor', 'b', 'DisplayName', 'Transmitter (Tx)');
text(pos_tx(1) + 10, pos_tx(2) + 10, 'Tx', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'b');

% 绘制接收基站 (Rx1 和 Rx2) - 使用红色上三角表示
plot(pos_rx1(1), pos_rx1(2), 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'Receiver (Rx1)');
text(pos_rx1(1) + 10, pos_rx1(2), 'Rx1', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');
plot(pos_rx2(1), pos_rx2(2), 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'Receiver (Rx2)'); % 不重复显示图例
text(pos_rx2(1) + 10, pos_rx2(2), 'Rx2', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'r');

% 绘制无人机目标 - 使用绿色圆圈表示
plot(target_pos(1), target_pos(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'UAV Target');
text(target_pos(1) + 10, target_pos(2), 'Target', 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'g');

% 绘制基站之间的连线 (可选，用于显示三角形结构)
plot([pos_tx(1), pos_rx1(1)], [pos_tx(2), pos_rx1(2)], 'k--'); % Tx to Rx1
plot([pos_tx(1), pos_rx2(1)], [pos_tx(2), pos_rx2(2)], 'k--'); % Tx to Rx2
plot([pos_rx1(1), pos_rx2(1)], [pos_rx1(2), pos_rx2(2)], 'k--'); % Rx1 to Rx2

% --- 图形设置 ---
title('Base Station and UAV Target Position Diagram'); % 设置标题
xlabel('X Coordinate (m)');                    % 设置 X 轴标签
ylabel('Y Coordinate (m)');                    % 设置 Y 轴标签
grid on;                               % 显示网格
% 设置坐标轴范围
% xlim([-100, 500]); % 设置 X 轴范围，例如从 -50 到 500

axis equal;                            % 设置 X 和 Y 轴具有相同的比例，确保几何形状不失真
ylim([-350, 350]); % 设置 Y 轴范围，使其更大，例如从 -350 到 350
legend('show', 'Location', 'northwest');    % 显示图例
hold off;                              % 结束叠加绘制