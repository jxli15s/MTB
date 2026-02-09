a=[[1,0];[0.5,sqrt(3)/2]];
b=inv(a')*2*pi;
x=[0,0.333333,0.666667]
y=[0,0.333333,0.666667]
X=x*a(1,1)+y*a(2,1);
Y=x*a(1,2)+y*a(2,2);
% Z=[5.882,-5.706,-0.000177] %MoSe2 AA norelax
% Z=[3.84,-8.96,5.12] %MoSe2 AA relax
% Z = [-2.48,-2.07,4.56] %MoSe2 AB relax
%Z = [1.76,-1.805,0] %MoS2 norelax AA
% Z = [2.45,-6.25,3.80] %MoS2 relax AA
% Z = [1.087,-0.54,-0.54] %MoS2 norelax AB
Z = [-0.784,-4.147,4.931] % MoS2 relax AB
g1=b(2,:);
g2=-b(1,:);
g3=g2-g1;
g4=-g1;
g5=-g2;
g6=-g3;
g=[g1;g2;g3;g4;g5;g6]
d1=X*g1(1)+Y*g1(2);
d2=X*g2(1)+Y*g2(2);
d3=X*g3(1)+Y*g3(2);
d4=X*g4(1)+Y*g4(2);
d5=X*g5(1)+Y*g5(2);
d6=X*g6(1)+Y*g6(2);

X=X'
Y=Y'
Z=Z'
fitCustomXYDataWithExternalParams(X, Y, Z, g);


function fitCustomXYDataWithExternalParams(x, y, z, g)
    % FITCUSTOMXYDATAWITHEXTERNALPARAMS Fits (x, y, z) data using a custom fit type with external parameters
    % Usage: fitCustomXYDataWithExternalParams(x, y, z, a, b)
    
    % Validate the input
    if length(x) ~= length(y) || length(y) ~= length(z)
        error('x, y, and z must have the same length');
    end



    % Define a custom fit type with external parameters using anonymous function

    customFitType=fittype(@(V,phi,X,Y)  real(V*exp(1j*deg2rad(phi))*exp(1j*(X*g(1,1)+Y*g(1,2)))+...
                                        V*exp(1j*deg2rad(phi))*exp(1j*(X*g(3,1)+Y*g(3,2)))+...
                                        V*exp(1j*deg2rad(phi))*exp(1j*(X*g(5,1)+Y*g(5,2)))+...
                                        V*exp(-1j*deg2rad(phi))*exp(1j*(X*g(2,1)+Y*g(2,2)))+...
                                        V*exp(-1j*deg2rad(phi))*exp(1j*(X*g(4,1)+Y*g(4,2)))+...
                                        V*exp(-1j*deg2rad(phi))*exp(1j*(X*g(6,1)+Y*g(6,2)))),...
                               'independent', {'X','Y'})

    % Define lower and upper bounds for the coefficients
    lowerBounds = [-4, -180];
    upperBounds = [4, 180];

    % Perform the fit
    options = fitoptions('Method', 'NonlinearLeastSquares', ...
                         'Lower', lowerBounds, 'Upper', upperBounds);

    % Perform the fit
    [fitresult, gof] = fit([x, y], z, customFitType,options);

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
    % title(['3D Data Fit with Custom Fit Type and External Parameters a=', num2str(a), ', b=', num2str(b)]);
    legend('Data', 'Fitted Surface');
    hold off;
end
