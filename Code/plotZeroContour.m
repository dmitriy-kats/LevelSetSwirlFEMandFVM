function h = plotZeroContour(xnodes,phi,og)

x = xnodes(:,1);
y = xnodes(:,2);
xu = unique(x);
yu = unique(y);
nx = length(xu);
ny = length(yu);
[x,y] = meshgrid(xu,yu);
phiPlot = reshape(phi,nx,ny)';
if og==1
    h = contour(x,y,phiPlot,[0 0],':k','LineWidth',2);
else
    h = contour(x,y,phiPlot,[0 0],'-','LineWidth',0.5);
end

