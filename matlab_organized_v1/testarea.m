clc;
clear;
x=linspace(-2,2)
y=linspace(-2,2)
[X,Y]=meshgrid(x,y)
Z=X.^2+Y.^2
[C, h] = contour(X, Y, Z, [2,2]);
sz = size(h.ContourMatrix,2);
nn(1) = h.ContourMatrix(2,1);
xx = h.ContourMatrix(1,2:nn(1)+1);
yy = h.ContourMatrix(2,2:nn(1)+1);
area(1) = polyarea(xx,yy)