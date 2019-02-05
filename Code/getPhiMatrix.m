function [x, y, phiPlot] = getPhiMatrix(xnodes,phi)

x = xnodes(:,1);
y = xnodes(:,2);
xu = unique(x);
yu = unique(y);
nx = length(xu);
ny = length(yu);
[x,y] = meshgrid(xu,yu);
phiPlot = reshape(phi,nx,ny)';


