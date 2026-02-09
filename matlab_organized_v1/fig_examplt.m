clc
clear
% 生成一些示例数据
x = rand(100, 1) * 10;
y = rand(100, 1) * 10;
z = 3 * x.^2 + 2 * y.^2 + x .* y + 5 * x + 4 * y + randn(100, 1) * 10; % 生成一些随机z数据

% 外部参数
a = 3; % 将a设置为3
b = 2; % 将b设置为2
a=[3,2]
% 调用函数拟合和绘制数据
fitCustomXYDataWithExternalParams(x, y, z, a);


function fitCustomXYDataWithExternalParams(x, y, z, a)
    % FITCUSTOMXYDATAWITHEXTERNALPARAMS Fits (x, y, z) data using a custom fit type with external parameters
    % Usage: fitCustomXYDataWithExternalParams(x, y, z, a, b)
    
    % Validate the input
    if length(x) ~= length(y) || length(y) ~= length(z)
        error('x, y, and z must have the same length');
    end

    % Define a custom fit type with external parameters using anonymous function
    customFitType = fittype(@(c, d, e, f, xdata, ydata) a(1)*xdata.^2 + a(2)*ydata.^2 + c*xdata.*ydata + d*xdata + e*ydata + f, ...
                            'independent', {'xdata', 'ydata'}, ...
                            'coefficients', {'c', 'd', 'e', 'f'});

    % Perform the fit
    [fitresult, gof] = fit([x, y], z, customFitType);

    % Display the fit result
    disp('Fit result:');
    disp(fitresult);
    disp('Goodness of fit:');
    disp(gof);

    % Plot the original data
    figure;
    scatter3(x, y, z, 'filled');
    hold on;

    % Plot the fitted surface
    [X, Y] = meshgrid(linspace(min(x), max(x), 100), linspace(min(y), max(y), 100));
    Z = feval(fitresult, X, Y);
    mesh(X, Y, Z);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['3D Data Fit with Custom Fit Type and External Parameters a=', num2str(a(1)), ', b=', num2str(a(2))]);
    legend('Data', 'Fitted Surface');
    hold off;
end
