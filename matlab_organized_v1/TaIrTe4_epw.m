% -------------------------
%   compute_a2F_qv.m
% -------------------------
clc;
clear;
% ==== 1. 读取 linewidth 文件 ====
linewidth_file = 'data/TaIrTe4/sc/ef-129/linewidth.phself.0.000K';
fid = fopen(linewidth_file);
if fid == -1
    error(' 无法打开文件: %s', linewidth_file);
end

% 自动跳过以 # 开头的注释行，合并多个空格
linewidth_data = textscan(fid, '%d %d %f %f', ...
    'MultipleDelimsAsOne', true, ...
    'CommentStyle', '#');
fclose(fid);

iq_list     = linewidth_data{1};  % q 点索引（从 1 开始）
imode_list  = linewidth_data{2};  % mode 索引（从 1 开始）
omega_list  = linewidth_data{3};  % 单位 cm^-1
ph_spectrum=load('data/TaIrTe4/sc/ef-129/freq.plot');
omega_list  = ph_spectrum(:,2);  % 单位 cm-1
% gamma_list  = linewidth_data{4};  % 如有需要可解注
nq = max(iq_list);
nmode = max(imode_list);
npoint = length(iq_list);

%
% ==== 2. 读取 lambda 文件 ====
lambda_file = 'data/TaIrTe4/sc/ef-129/lambda_0k.dat';
lambda_data = readmatrix(lambda_file);  % [nq × nmode]

% omega_matrix=reshape(omega_list,[nmode,nq])';
omega_matrix=reshape(omega_list,[nq,nmode]);
lambda_matrix=lambda_data(:,2:end);
a2f_matrix=omega_matrix.*lambda_matrix./2.0;

% ==== 4. 保存结果 ====
outfile = 'data/TaIrTe4/sc/ef-129/a2F_qv.dat';
fid = fopen(outfile, 'w');
fprintf(fid, '# iq  omega[cm^-1]     lambda_qv        a2F_qv\n');
for v_idx = 1:nmode
    for q_idx = 1:nq
    fprintf(fid, '%4d  %14.6e  %14.6e  %14.6e\n', ...
        q_idx, omega_matrix(q_idx,v_idx), lambda_matrix(q_idx,v_idx),a2f_matrix(q_idx,v_idx));
    end
    fprintf(fid,'\n');
end
fclose(fid);

disp('✅ 成功完成！结果保存在 a2F_qv.dat');

%%
data = readmatrix('data/TaIrTe4/sc/ef-129/kpoints.dat'); 
kpath = data(:,1);
labels={'R','Y','\Gamma','X'}; % labels for k
kk=[0.000000,0.823225,1.073469,1.896693];
linesize=2;
xrange=[kpath(1),kpath(end)];
figure();
% % Energy=omega_matrix*0.2418; %mev to Thz
Energy=omega_matrix/8.1; %cm-1 to meV

yrange=[0,max(Energy,'','all')+1];
for i=1:length(Energy(1,:))
    plot(kpath(:),Energy(:,i),'Color','black','LineWidth',linesize);
    hold on
end

% % % plot(kpath,zeros(1,length(kpath))+1.6,'--red','LineWidth',1)

for i=1:length(kk)-2
     plot([kk(i+1) kk(i+1)],[min(min(Energy))-1 max(max(Energy))+1],'--k','LineWidth',linesize)
end
grid off;
box on;
xlim(xrange)
xticks(kk)
xticklabels(labels)
ylim(yrange)
ylabel('\omega (meV)','FontSize',24)

% 2. scatter 所有点，加上 λ 的大小
[q_idx, ~] = meshgrid(kpath, 1:nmode);
q_idx=q_idx';
q_flat     = q_idx(:);
omega_flat = omega_matrix(:)/8.1;
lambda_flat = lambda_matrix(:);
a2f_flat=a2f_matrix(:);
circle_size = 800 * abs(lambda_flat+0.000001);  % 调节可视大小

scatter(q_flat, omega_flat, circle_size, ...
    'filled', 'MarkerFaceAlpha', 0.5, ...
    'MarkerFaceColor', [0.2 0.4 0.9]);


% 美化
xlabel('q-point index');
% ylabel('\omega (cm^{-1})');
% ylabel('\omega (THz)');
ylabel('\omega (meV)');
title('Phonon dispersion with \lambda_{qv} size overlay');
grid on;

%%

figure;
hold on;
for imode = 1:nmode
    plot(1:nq, omega_matrix(:, imode), 'LineWidth', 1.2);
end
xlabel('q-point index');
ylabel('\omega (cm^{-1})');
title('Phonon Dispersion');
% grid on;
box on;


%%
figure;
hold on;

% 1. plot 所有 mode 的声子谱线
for imode = 1:nmode
    plot(1:nq, omega_matrix(:, imode)*0.2418, 'k-', 'LineWidth', 0.8);  % 黑色细线
end

yline(1.6, 'r--', 'y = 1.6', 'LineWidth', 1.2, ...
    'LabelHorizontalAlignment', 'left', ...
    'LabelVerticalAlignment', 'bottom');

% 2. scatter 所有点，加上 λ 的大小
[q_idx, mode_idx] = meshgrid(1:nq, 1:nmode);
q_idx=q_idx';
q_flat     = q_idx(:);
omega_flat = omega_matrix(:);
lambda_flat = lambda_matrix(:);
a2f_flat=a2f_matrix(:);
circle_size = 500 * abs(lambda_flat+0.000001);  % 调节可视大小

scatter(q_flat, omega_flat*0.2418, circle_size, ...
    'filled', 'MarkerFaceAlpha', 0.5, ...
    'MarkerFaceColor', [0.2 0.4 0.9]);


% 美化
xlabel('q-point index');
% ylabel('\omega (cm^{-1})');
ylabel('\omega (THz)');
title('Phonon dispersion with \lambda_{qv} size overlay');
grid on;

%%
% 设置k空间范围
N = 500;
kx = linspace(-pi, pi, N);
ky = linspace(-pi, pi, N);
[KX, KY] = meshgrid(kx, ky);

% 计算等值函数 f(kx, ky) = cos(kx) + cos(ky)
F = cos(KX) + cos(KY);

% 设置等值线的目标值（可调）
c = 1;  % 你可以改为 -1, 0.5, 1.5 等等

% 绘制等值线
figure;
contour(KX, KY, F, [c c], 'LineWidth', 2);
xlabel('k_x');
ylabel('k_y');
title(['\cos(k_x) + \cos(k_y) = ', num2str(c)]);
axis equal;
xlim([-pi, pi]);
ylim([-pi, pi]);
set(gca, 'XTick', [-pi 0 pi], 'XTickLabel', {'-\pi','0','\pi'});
set(gca, 'YTick', [-pi 0 pi], 'YTickLabel', {'-\pi','0','\pi'});
grid on;
%%
data = readmatrix('data/BaCoAl/a2f_w_0k.dat');   % 第一列: omega, 第二列: a2F(omega)
omega = data(:,1); a2f = data(:,2);

mask = omega > 0;              % 避免除零，从正频开始
omega = omega(mask); a2f = a2f(mask);

lambda_cum = 2 * cumtrapz(omega, a2f./omega);
lambda_tot = lambda_cum(end);
fprintf('lambda_total = %.6f\n', lambda_tot);

% (可选) 计算 ω_log
wlog = exp( (2/lambda_tot) * trapz(omega, (a2f./omega).*log(omega)) );
fprintf('omega_log = %.6f (same unit as omega)\n', wlog);
% 计算并绘制 λ 的累积值
figure;
plot(omega*0.2418, lambda_cum, 'r-', 'LineWidth', 1.5);
hold on;
plot(omega*0.2418,a2f,'b-','LineWidth',1.5)
xlabel('\omega (THz)');
ylabel('\lambda_{cum}');
title('Cumulative \lambda vs \omega');
grid on;


%%
% data = readmatrix('data/BaCoAl/a2f_w_0k.dat');   % 第一列: omega, 第二列: a2F(omega)
data = readmatrix('data/TaIrTe4/sc/ef-129/a2f_w_10k_mesh.dat');   % 第一列: omega, 第二列: a2F(omega)
omega = data(:,1); a2f = data(:,2);

mask = omega > 0;              % 避免除零，从正频开始
omega = omega(mask); a2f = a2f(mask);

f = a2f ./ omega;
% 高斯核
sigma = 0.15;               % 核宽度，单位跟 omega 一致
dx = mean(diff(omega));
w = round(5*sigma/dx);
g = exp(-(-w:w).^2/(2*(sigma/dx)^2));
g = g / sum(g);


% 卷积平滑
f_smooth = conv(f, g, 'same');
f_smooth_ori = conv(a2f, g, 'same');

% 积分
lambda_cum = 2 * cumtrapz(omega, f_smooth);


% (可选) 计算 ω_log
wlog = exp( (2/lambda_tot) * trapz(omega, (a2f./omega).*log(omega)) );
fprintf('omega_log = %.6f (same unit as omega)\n', wlog);
% 计算并绘制 λ 的累积值
figure;
plot(omega, lambda_cum, 'b-', 'LineWidth', 1.5);
hold on;
plot(omega,f_smooth,'b--','LineWidth',1.5)
plot(omega,f_smooth_ori,'r-','LineWidth',1.5)
xlabel('\omega (THz)');
ylabel('\lambda_{cum}');
title('Cumulative \lambda vs \omega');
grid on;
% integrand
% f = a2f ./ omega;

%%
figure('Color','w');
ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);
%
% figure;
% subplot(1,2,1);
% data = readmatrix('data/BaCoAl/meshes/kpoints.dat'); 
data = readmatrix('data/TaIrTe4/sc/ef-129/kpoints.dat'); 
kpath = data(:,1);
% labels={'\Gamma','M','K','\Gamma','A','L','H','A'}; % labels for k
% kk=[0.000000,1.047198,1.855162,3.2199,3.534074,4.581272,5.389236,6.753990];
labels={'R','Y','\Gamma','X'}; % labels for k
kk=[0.000000,0.823225,1.073469,1.896693];
linesize=2;
xrange=[kpath(1),kpath(end)];
yrange=[0,max(Energy,'','all')+1];
hold(ax1,'on')
for i=1:length(Energy(1,:))
    plot(ax1, kpath(:),Energy(:,i),'Color','black','LineWidth',linesize);
end
%
plot(ax1, kpath,zeros(1,length(kpath))+1.6,'--red','LineWidth',1)

for i=1:length(kk)-2
     plot(ax1,[kk(i+1) kk(i+1)],[min(min(Energy))-2 max(max(Energy))+2],'--k','LineWidth',1)
end
% grid off
box(ax1,'on')
xlim(ax1, xrange)
xticks(ax1,kk)
xticklabels(ax1,labels)
ylim(ax1, yrange)
ylabel(ax1,'\omega (THz)','FontSize',24)
ax1.FontSize = 20;
%
% 2. scatter 所有点，加上 λ 的大小
[q_idx, ~] = meshgrid(kpath, 1:nmode);
q_idx=q_idx';
q_flat     = q_idx(:);
omega_flat = omega_matrix(:);
lambda_flat = lambda_matrix(:);
a2f_flat=a2f_matrix(:);
circle_size = 800 * abs(lambda_flat+0.000001);  % 调节可视大小

h0=scatter(ax1, q_flat, omega_flat, circle_size, ...
    'filled', 'MarkerFaceAlpha', 0.5, ...
    'MarkerFaceColor', [0.2 0.4 0.9]);
lgd=legend(ax1, h0, { '$\lambda_{qv}$'}, ...
       'Interpreter','latex', 'Location','best');
lgd.FontSize = 16;

% h1=plot(ax2,lambda_cum,omega*0.2418, 'w-', 'LineWidth', 1.5);
% hold(ax2,"on");
% plot(omega*0.2418,a2f,'b--','LineWidth',1.5)
% h2=plot(ax2, f_smooth*2,omega*0.2418,'r-','LineWidth',1.5)
h2=plot(ax2, f_smooth_ori,omega*8.1,'r-','LineWidth',1.5)
lgd2=legend(ax2, h2, { '$\sum_{qv}\lambda_{qv}\delta(\omega-\omega_{qv})$'}, ...
       'Interpreter','latex', 'Location','best');
lgd2.FontSize = 8
% 
% legend(ax2, [h1 h2], {'$\lambda(\omega)$', '$\alpha^2F(\omega)$'}, ...
%        'Interpreter','latex', 'Location','best');
% xlabel('\omega (THz)');
yticks(ax2,'')
ylim(ax2,yrange)
xlabel(ax2,'\lambda');
xlim(ax2,[0, 2])
% title(ax2,'Cumulative \lambda vs \omega');
ax2.FontSize = 16;
set(ax1,'Units','normalized','Position',[0.08 0.12 0.62 0.82]); % 左：更宽
set(ax2,'Units','normalized','Position',[0.72 0.12 0.1 0.82]); % 右：更窄

%%
% -------------------------
%   compute_a2F_qv.m
% -------------------------
clc;
clear;
% ==== 1. 读取 linewidth 文件 ====
linewidth_file = 'data/TaIrTe4/sc/ef-129/linewidth.phself.0.000K';
fid = fopen(linewidth_file);
if fid == -1
    error(' 无法打开文件: %s', linewidth_file);
end

% 自动跳过以 # 开头的注释行，合并多个空格
linewidth_data = textscan(fid, '%d %d %f %f', ...
    'MultipleDelimsAsOne', true, ...
    'CommentStyle', '#');
fclose(fid);

iq_list     = linewidth_data{1};  % q 点索引（从 1 开始）
imode_list  = linewidth_data{2};  % mode 索引（从 1 开始）
omega_list  = linewidth_data{3};  % 单位 cm^-1
gamma_list  = linewidth_data{4};  % 如有需要可解注
nq = max(iq_list);
nmode = max(imode_list);
npoint = length(iq_list);

%
% ==== 2. 读取 lambda 文件 ====
lambda_file = 'data/TaIrTe4/sc/ef-129/lambda_0k.dat';
lambda_data = readmatrix(lambda_file);  % [nq × nmode]

omega_matrix=reshape(omega_list,[nmode,nq])';

lambda_matrix=lambda_data(:,2:end);
a2f_matrix=omega_matrix.*lambda_matrix./2.0;

% ==== 4. 保存结果 ====
outfile = 'data/TaIrTe4/sc/ef-129/a2F_qv.dat';
fid = fopen(outfile, 'w');
fprintf(fid, '# iq  omega[meV]     lambda_qv        a2F_qv\n');
for v_idx = 1:nmode
    for q_idx = 1:nq
    fprintf(fid, '%4d  %14.6e  %14.6e  %14.6e\n', ...
        q_idx, omega_matrix(q_idx,v_idx), lambda_matrix(q_idx,v_idx),a2f_matrix(q_idx,v_idx));
    end
    fprintf(fid,'\n');
end
fclose(fid);

disp('✅ 成功完成！结果保存在 a2F_qv.dat');
%%

data = readmatrix('data/TaIrTe4/sc/ef-129/kpoints.dat'); 
kpath = data(:,1);
labels={'R','Y','\Gamma','X'}; % labels for k
kk=[0.000000,0.823225,1.073469,1.896693];
linesize=2;
xrange=[kpath(1),kpath(end)];
figure();
Energy=omega_matrix;
%Energy=omega_matrix*0.2418; %mev to Thz


yrange=[0,max(Energy,'','all')+1];
for i=1:length(Energy(1,:))
    plot(kpath(:),Energy(:,i),'Color','black','LineWidth',linesize);
    hold on
end

% % % plot(kpath,zeros(1,length(kpath))+1.6,'--red','LineWidth',1)

for i=1:length(kk)-2
     plot([kk(i+1) kk(i+1)],[min(min(Energy))-1 max(max(Energy))+1],'--k','LineWidth',linesize)
end
grid off;
box on;
xlim(xrange)
xticks(kk)
xticklabels(labels)
ylim(yrange)
ylabel('\omega (THz)','FontSize',24)
%
% 2. scatter 所有点，加上 λ 的大小
[q_idx, ~] = meshgrid(kpath, 1:nmode);
q_idx=q_idx';
q_flat     = q_idx(:);
% omega_flat = omega_matrix(:)*0.2418;% meV to THz
omega_flat = omega_matrix(:);
lambda_flat = lambda_matrix(:);
a2f_flat=a2f_matrix(:);
circle_size = 800 * abs(lambda_flat+0.000001);  % 调节可视大小

% % scatter(q_flat, omega_flat, circle_size, ...
% %     'filled', 'MarkerFaceAlpha', 0.5, ...
% %     'MarkerFaceColor', [0.2 0.4 0.9]);


% 美化
xlabel('q-point index');
ylabel('\omega (THz)');
ylabel('\omega (meV)');
title('Phonon dispersion with \lambda_{qv} size overlay');
grid on;
%%
data = readmatrix('data/BaCoAl/a2f_w_0k.dat');   % 第一列: omega, 第二列: a2F(omega)
omega = data(:,1); a2f = data(:,2);

mask = omega > 0;              % 避免除零，从正频开始
omega = omega(mask); a2f = a2f(mask);

lambda_cum = 2 * cumtrapz(omega, a2f./omega);
lambda_tot = lambda_cum(end);
fprintf('lambda_total = %.6f\n', lambda_tot);

% (可选) 计算 ω_log
wlog = exp( (2/lambda_tot) * trapz(omega, (a2f./omega).*log(omega)) );
fprintf('omega_log = %.6f (same unit as omega)\n', wlog);
% 计算并绘制 λ 的累积值
figure;
plot(omega*0.2418, lambda_cum, 'r-', 'LineWidth', 1.5);
hold on;
plot(omega*0.2418,a2f,'b-','LineWidth',1.5)
xlabel('\omega (THz)');
ylabel('\lambda_{cum}');
title('Cumulative \lambda vs \omega');
grid on;
%%
% data = readmatrix('data/BaCoAl/a2f_w_0k.dat');   % 第一列: omega, 第二列: a2F(omega)
data = readmatrix('data/TaIrTe4/sc/ef-129/a2f_w_10k_mesh.dat');   % 第一列: omega, 第二列: a2F(omega)
omega = data(:,1); a2f = data(:,2);

mask = omega > 0;              % 避免除零，从正频开始
omega = omega(mask); a2f = a2f(mask);

f = a2f ./ omega;
% 高斯核
sigma = 0.15;               % 核宽度，单位跟 omega 一致
dx = mean(diff(omega));
w = round(5*sigma/dx);
g = exp(-(-w:w).^2/(2*(sigma/dx)^2));
g = g / sum(g);


% 卷积平滑
f_smooth = conv(f, g, 'same');
f_smooth_ori = conv(a2f, g, 'same');

% 积分
lambda_cum = 2 * cumtrapz(omega, f_smooth);


% (可选) 计算 ω_log
wlog = exp( (2/lambda_tot) * trapz(omega, (a2f./omega).*log(omega)) );
fprintf('omega_log = %.6f (same unit as omega)\n', wlog);
% 计算并绘制 λ 的累积值
figure;
plot(omega, lambda_cum, 'b-', 'LineWidth', 1.5);
hold on;
plot(omega,f_smooth,'b--','LineWidth',1.5)
plot(omega,f_smooth_ori,'r-','LineWidth',1.5)
xlabel('\omega (THz)');
ylabel('\lambda_{cum}');
title('Cumulative \lambda vs \omega');
grid on;
% integrand
% f = a2f ./ omega;

%%
figure('Color','w');
ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);
%
% figure;
% subplot(1,2,1);
% data = readmatrix('data/BaCoAl/meshes/kpoints.dat'); 
data = readmatrix('data/TaIrTe4/sc/ef-129/kpoints.dat'); 
kpath = data(:,1);
% labels={'\Gamma','M','K','\Gamma','A','L','H','A'}; % labels for k
% kk=[0.000000,1.047198,1.855162,3.2199,3.534074,4.581272,5.389236,6.753990];
labels={'R','Y','\Gamma','X'}; % labels for k
kk=[0.000000,0.823225,1.073469,1.896693];
linesize=2;
xrange=[kpath(1),kpath(end)];
yrange=[0,max(Energy,'','all')+1];
hold(ax1,'on')
for i=1:length(Energy(1,:))
    plot(ax1, kpath(:),Energy(:,i),'Color','black','LineWidth',linesize);
end
%
% plot(ax1, kpath,zeros(1,length(kpath))+1.6,'--red','LineWidth',1)

for i=1:length(kk)-2
     plot(ax1,[kk(i+1) kk(i+1)],[min(min(Energy))-2 max(max(Energy))+2],'--k','LineWidth',1)
end
% grid off
box(ax1,'on')
xlim(ax1, xrange)
xticks(ax1,kk)
xticklabels(ax1,labels)
ylim(ax1, yrange)
ylabel(ax1,'\omega (THz)','FontSize',24)
ax1.FontSize = 20;
%
% 2. scatter 所有点，加上 λ 的大小
[q_idx, ~] = meshgrid(kpath, 1:nmode);
q_idx=q_idx';
q_flat     = q_idx(:);
omega_flat = omega_matrix(:)*0.2418;
lambda_flat = lambda_matrix(:);
a2f_flat=a2f_matrix(:);
circle_size = 800 * abs(lambda_flat+0.000001);  % 调节可视大小

h0=scatter(ax1, q_flat, omega_flat, circle_size, ...
    'filled', 'MarkerFaceAlpha', 0.5, ...
    'MarkerFaceColor', [0.2 0.4 0.9]);
lgd=legend(ax1, h0, { '$\lambda_{qv}$'}, ...
       'Interpreter','latex', 'Location','best');
lgd.FontSize = 16;

% h1=plot(ax2,lambda_cum,omega*0.2418, 'w-', 'LineWidth', 1.5);
% hold(ax2,"on");
% plot(omega*0.2418,a2f,'b--','LineWidth',1.5)
% h2=plot(ax2, f_smooth*2,omega*0.2418,'r-','LineWidth',1.5)
h2=plot(ax2, f_smooth_ori,omega,'r-','LineWidth',1.5)
lgd2=legend(ax2, h2, { '$\sum_{qv}\lambda_{qv}\delta(\omega-\omega_{qv})$'}, ...
       'Interpreter','latex', 'Location','best');
lgd2.FontSize = 8
% 
% legend(ax2, [h1 h2], {'$\lambda(\omega)$', '$\alpha^2F(\omega)$'}, ...
%        'Interpreter','latex', 'Location','best');
% xlabel('\omega (THz)');
yticks(ax2,'')
ylim(ax2,yrange)
xlabel(ax2,'\lambda');
xlim(ax2,[0, 2])
% title(ax2,'Cumulative \lambda vs \omega');
ax2.FontSize = 16;
set(ax1,'Units','normalized','Position',[0.08 0.12 0.62 0.82]); % 左：更宽
set(ax2,'Units','normalized','Position',[0.72 0.12 0.1 0.82]); % 右：更窄