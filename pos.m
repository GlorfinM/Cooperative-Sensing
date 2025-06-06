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
target_vel = [10, 5];

% --- 绘图 ---
% 创建一个新的图形窗口
figure;
hold on; % 允许在同一图形上叠加绘制多个对象

% 绘制发射基站 (Tx) - 使用红色方块表示
plot(pos_tx(1), pos_tx(2), 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'Tx');


% 绘制接收基站 (Rx1 和 Rx2) - 使用红色上三角表示
plot(pos_rx1(1), pos_rx1(2), 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'Rx1');
plot(pos_rx2(1), pos_rx2(2), 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'DisplayName', 'Rx2'); 


% 绘制无人机目标 - 使用绿色圆圈表示
plot(target_pos(1), target_pos(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'UAV Target');
% 绘制目标速度矢量
quiver(target_pos(1), target_pos(2), target_vel(1)*5, target_vel(2)*5, 'm', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'AutoScale', 'off', 'DisplayName', 'Velocity Vector'); 
text(target_pos(1) + target_vel(1)*5 + 5, target_pos(2) + target_vel(2)*5 + 5, sprintf('(%.2f, %.2f)', target_vel(1), target_vel(2)), 'FontSize', 9, 'Color', 'm');

% 绘制基站之间的连线
h_line1 = plot([pos_tx(1), pos_rx1(1)], [pos_tx(2), pos_rx1(2)], 'k--'); % Tx 到 Rx1
set(get(get(h_line1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h_line2 = plot([pos_tx(1), pos_rx2(1)], [pos_tx(2), pos_rx2(2)], 'k--'); % Tx 到 Rx2
set(get(get(h_line2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h_line3 = plot([pos_rx1(1), pos_rx2(1)], [pos_rx1(2), pos_rx2(2)], 'k--'); % Rx1 到 Rx2
set(get(get(h_line3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% --- 绘制信号路径 ---
signal_path_color = 'b'; % 信号路径线条颜色
signal_path_linestyle = '--'; % 信号路径线型
signal_path_linewidth = 2;   % 信号路径线宽

% 箭头参数
arrow_color = 'b'; % 箭头颜色 (黑色)
arrow_scale = 0.7; % 箭头整体缩放因子
arrow_linewidth = 2; % 箭头线宽
arrow_headsize_multiplier = 8; % 箭头头部大小乘数 (实际MaxHeadSize = multiplier * scale)
arrow_body_length_factor = 0.25; % 箭头箭身长度因子 (相对于线段向量的部分长度)
arrow_tail_offset_factor = 0.125; % 箭头尾部从线段中点向后的偏移因子

% 1. Tx -> UAV 信号路径
line_path_tx_uav = [pos_tx; target_pos];
plot(line_path_tx_uav(:,1), line_path_tx_uav(:,2), ...
     'Color', signal_path_color, 'LineStyle', signal_path_linestyle, 'LineWidth', signal_path_linewidth, ...
     'HandleVisibility', 'off');

mid_point_tx_uav = mean(line_path_tx_uav);
vec_tx_uav = diff(line_path_tx_uav);
quiver(mid_point_tx_uav(1) - vec_tx_uav(1)*arrow_tail_offset_factor*arrow_scale, ...
       mid_point_tx_uav(2) - vec_tx_uav(2)*arrow_tail_offset_factor*arrow_scale, ...
       vec_tx_uav(1)*arrow_body_length_factor*arrow_scale, ...
       vec_tx_uav(2)*arrow_body_length_factor*arrow_scale, 0, ...
       'Color', arrow_color, 'LineWidth', arrow_linewidth, ... 
       'MaxHeadSize', arrow_headsize_multiplier*arrow_scale, 'HandleVisibility', 'off');

% 2. UAV -> Rx1 信号路径
line_path_uav_rx1 = [target_pos; pos_rx1];
plot(line_path_uav_rx1(:,1), line_path_uav_rx1(:,2), ...
     'Color', signal_path_color, 'LineStyle', signal_path_linestyle, 'LineWidth', signal_path_linewidth, ...
     'HandleVisibility', 'off');

mid_point_uav_rx1 = mean(line_path_uav_rx1);
vec_uav_rx1 = diff(line_path_uav_rx1);
quiver(mid_point_uav_rx1(1) - vec_uav_rx1(1)*arrow_tail_offset_factor*arrow_scale, ...
       mid_point_uav_rx1(2) - vec_uav_rx1(2)*arrow_tail_offset_factor*arrow_scale, ...
       vec_uav_rx1(1)*arrow_body_length_factor*arrow_scale, ...
       vec_uav_rx1(2)*arrow_body_length_factor*arrow_scale, 0, ...
       'Color', arrow_color, 'LineWidth', arrow_linewidth, ... 
       'MaxHeadSize', arrow_headsize_multiplier*arrow_scale, 'HandleVisibility', 'off');

% 3. UAV -> Rx2 信号路径
line_path_uav_rx2 = [target_pos; pos_rx2];
plot(line_path_uav_rx2(:,1), line_path_uav_rx2(:,2), ...
     'Color', signal_path_color, 'LineStyle', signal_path_linestyle, 'LineWidth', signal_path_linewidth, ...
     'DisplayName', 'Signal Path');

mid_point_uav_rx2 = mean(line_path_uav_rx2);
vec_uav_rx2 = diff(line_path_uav_rx2);
quiver(mid_point_uav_rx2(1) - vec_uav_rx2(1)*arrow_tail_offset_factor*arrow_scale, ...
       mid_point_uav_rx2(2) - vec_uav_rx2(2)*arrow_tail_offset_factor*arrow_scale, ...
       vec_uav_rx2(1)*arrow_body_length_factor*arrow_scale, ...
       vec_uav_rx2(2)*arrow_body_length_factor*arrow_scale, 0, ...
       'Color', arrow_color, 'LineWidth', arrow_linewidth, ... 
       'MaxHeadSize', arrow_headsize_multiplier*arrow_scale, 'HandleVisibility', 'off');

% --- 图形设置 ---
title('Base Station and UAV Target Position Diagram'); 
xlabel('X Coordinate (m)');                    
ylabel('Y Coordinate (m)');                    
grid on;                               
% 设置坐标轴范围
axis equal;                            
ylim([-350, 350]); %
legend('show', 'Location', 'northwest');   
hold off;                              