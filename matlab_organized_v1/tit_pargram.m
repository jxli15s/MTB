% 加载数据
data = readmatrix('/Volumes/T9/work/test/sparse/15x1/10nmd/pargram/u-epsilon/onsite-diff/data-matlab-2.dat'); % 确保文件与脚本在同一目录下
x = data(:, 1);
y = data(:, 2);
z = data(:, 3);

% 创建唯一的 x 和 y 网格
[x_unique, ~, x_idx] = unique(x);
[y_unique, ~, y_idx] = unique(y);

% 初始化 Z 矩阵
Z = NaN(length(y_unique), length(x_unique));

% 填充 Z 矩阵
for i = 1:length(z)
    Z(y_idx(i), x_idx(i)) = z(i);
end

% 绘制相图
figure;
imagesc(x_unique, y_unique, Z); % 使用 imagesc 显示网格化数据

% 设置自定义柔和颜色（橙色、灰色、蓝色）
colormap([1 0.6 0.2; 0.7 0.7 0.7; 0.2 0.4 0.8]); % 橙色、灰色、蓝色

% 设置颜色范围和颜色条
caxis([1, 3]); % 限制 z 值范围
colorbar('Ticks', [1, 2, 3], 'TickLabels', {'1', '2', '3'}); % 添加颜色条
xlabel('X');
ylabel('Y');
title('Phase Diagram with Custom Colors (Orange, Gray, Blue)');

% 调整图像显示
axis tight; % 去掉多余空白


%%

% 绘制相图
figure;
imagesc(x_unique, y_unique, Z); % 使用 imagesc 显示网格化数据

% 设置自定义柔和颜色（橙色、灰色、蓝色）
colormap([1 0.6 0.2; 0.7 0.7 0.7; 0.2 0.4 0.8]); % 橙色、灰色、蓝色

% 设置颜色范围和颜色条
caxis([1, 3]); % 限制 z 值范围
colorbar('Ticks', [1, 2, 3], 'TickLabels', {'1', '2', '3'}); % 添加颜色条
xlabel('X');
ylabel('Y');
title('Phase Diagram with Custom Colors (Orange, Gray, Blue)');

% 调整 y 轴方向为从小到大
set(gca, 'YDir', 'normal'); % 修正 y 轴方向
axis tight; % 去掉多余空白
