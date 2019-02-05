clear all; clc;
close all;

subfldr={'/FEM/testc/'};

%% Problem Parameters
Lx = 1.0;
Ly = 1.0;
mu  = 1;
dt=1e-3; timesteps=8/dt;
alpha=1;

%% Mesh Parameters
nx = 250;      % number of elements in x direction
ny = 250;      % number of elements in y direction
ne = nx * ny; % total number of elements

nnx = nx + 1;    % number of nodes in x direction
nny = ny + 1;    % number of nodes in y direction
nn  = nnx * nny; % total number of nodes (unknowns)

h=Lx/nx;

%% Create mesh
[xnodes, nconn, surfconnD, surfconnR, surfconnU, surfconnL]=createFem2dMesh(Lx, Ly, nx, ny);

%% Initialize sparse matrix-vector problem
phi=sqrt((xnodes(:,1)-ones(nn,1).*0.5).^2+(xnodes(:,2)-ones(nn,1).*0.75).^2)-ones(nn,1).*0.15;
phiO=phi;

M = computeFemMMatrix(xnodes, nconn);

sf_dr=subfldr{1};

[xm, ym, phiP] = getPhiMatrix(xnodes,phi);
ff=figure('Visible','off');
contour(xm,ym,phiP);
hold on
contour(xm,ym,phiP,[0,0],'LineWidth',2);
hold off
axis([0 1 0 1])
axis equal
xlabel('X','Interpreter','LaTex');
ylabel('Y','Interpreter','LaTex');
title(['t= ' num2str(0) '\hspace{0.1in} N=' num2str(nx) 'x' num2str(ny) ...
' (FEM)'],'Interpreter','LaTex');
c=colorbar;
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',16,'FontName','AvantGarde');
set(c,'fontsize',12);
saveas(ff,[pwd sf_dr sprintf('FIG%d.png',1)]);
close all
scount=1;
SS(1).t=0;
SS(1).phi=phi;

for tt=1:timesteps
 sf_dr=subfldr{1};   
M_d = computeFemMdMatrix(xnodes, nconn,h,alpha,tt*dt);
C= computeFemCMatrix(xnodes, nconn,tt*dt);
C_d = computeFemCdMatrix(xnodes, nconn,h,alpha,tt*dt);

M_l=M+M_d;
M_l_inv=sparse(nn,nn);

for jj=1:nn
        tempsum=sum(M_l(jj,:));
        M_l(jj,:)=0;
        M_l(jj,jj)=tempsum;
        M_l_inv(jj,jj)=1/tempsum;
end    
    
    
phi =-M_l_inv*(C+C_d)*phi.*dt+phi;

if mod(tt,100)==0
    [xm, ym, phiP] = getPhiMatrix(xnodes,phi);
    scount=scount+1;
    SS(scount).t=dt*tt;
    SS(scount).phi=phi;
    ff=figure('Visible','off');
    contour(xm,ym,phiP);
    hold on
    contour(xm,ym,phiP,[0,0],'LineWidth',2);
    hold off
    axis([0 1 0 1])
    axis equal
    xlabel('X','Interpreter','LaTex');
    ylabel('Y','Interpreter','LaTex');
    title(['t= ' num2str(tt*dt) '\hspace{0.1in} N=' num2str(nx) 'x' num2str(ny) ...
    ' (FEM)'],'Interpreter','LaTex');
    c=colorbar;
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',16,'FontName','AvantGarde'); 
    set(c,'fontsize',12);
    saveas(ff,[pwd sf_dr sprintf('FIG%d.png',scount)]);
    close all
end

end
matfile = fullfile([pwd fileparts(sf_dr)], sprintf('TEST%d.mat',3));
save(matfile);
