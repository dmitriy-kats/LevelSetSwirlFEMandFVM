clc; clear; close all;

%test cases
MGk=[50 100 150 300]; %number of volumes in one direction
subfldr={'/FVM/testa/','/FVM/testb/','/FVM/testc/','/FVM/testd/'};
 

for k=1:length(MGk)
MG=MGk(k); %number of volumes in one direction
sf_dr=subfldr{k};
dt=0.5e-3;
%Initialize the values of phi where phi=0 is the interface
%Initial interface is circle of radius 1
si=linspace(0,2*pi,2000);
xx=cos(si)*0.15+0.5; %x location of parametrized curve
yy=sin(si)*0.15+0.75; %y location of parametrized curve
Wx=1; %Width of space
Hy=Wx; %Height of space
%MG=30; %number of grids to divide width/height
gs=(Wx)/(MG-1); %grid size
[xm, ym]=meshgrid(gs/2:gs:Wx-gs/2); %phi locations
phi=ones(MG-1,MG-1)*inf; %phi values
phiN=phi; %this will store the new phi value
%sets the phi values so phi=0 is at the interface
k_eq=22;

%intialize the value of phi to be a signed distance function
for j=1:length(ym) 
 for i=1:length(xm)
  phi(j,i)=norm([(xm(1,i)-0.5) (ym(j,1)-0.75)])-0.15;
end
end

%velocity field
uK=@(x,y,t)2*cos(pi*t/8)*(-sin(pi*x)^2*sin(pi*y)*cos(pi*y));
vK=@(x,y,t)2*cos(pi*t/8)*(sin(pi*y)^2*sin(pi*x)*cos(pi*x));

%this is for record keeping of the evolution of phi
scount=1;
SS(1).t=0;
SS(1).phi=phi;

h=figure('Visible','off');
contour(xm,ym,phi);
hold on
contour(xm,ym,phi,[0,0],'LineWidth',2);
hold off
axis([0 1 0 1])
axis equal
xlabel('X','Interpreter','LaTex');
ylabel('Y','Interpreter','LaTex');
title(['t= ' num2str(0) '\hspace{0.1in} N=' num2str(MG) 'x' num2str(MG) ...
' (FVM)'],'Interpreter','LaTex');
c=colorbar;
fig=gcf;
set(findall(fig,'-property','FontSize'),'FontSize',16,'FontName','AvantGarde');
set(c,'fontsize',12);
saveas(h,[pwd sf_dr sprintf('FIG%d.png',scount)]);
close all

for tt=1:8/dt 
for j=2:length(ym)-1
   for i=2:length(xm)-1
     %MAC grid naming
     xp=xm(1,i); yp=ym(j,1); xe=xp+gs/2; yn=yp+gs/2; xw=xp-gs/2; ys=yp-gs/2;
     xpe=xp+gs; xpw=xp-gs; ypn=yp+gs; yps=yp-gs;
     ct=tt*dt;
     ue=uK(xe,yp,ct); uw=uK(xw,yp,ct);
     vn=vK(xe,yn,ct); vs=vK(xe,ys,ct);
     up=uK(xp,yp,ct); vp=vK(xp,yp,ct); 
     upe=uK(xpe,yp,ct); vpe=vK(xpe,yp,ct);
     upw=uK(xpw,yp,ct); vpw=vK(xpw,yp,ct);
     upn=uK(xp,ypn,ct); vpn=vK(xp,ypn,ct);
     ups=uK(xp,yps,ct); vps=vK(xp,yps,ct);
          
     %gradients in x and y
     gpx=getkappa(up,k_eq)*(phi(j,i+1)-phi(j,i))/gs + ...
        (1-getkappa(up,k_eq))*(phi(j,i)-phi(j,i-1))/gs;
     gpy=getkappa(vp,k_eq)*(phi(j+1,i)-phi(j,i))/gs + ...
                (1-getkappa(vp,k_eq))*(phi(j,i)-phi(j-1,i))/gs;
     
     %gradients in other volumes
     if i+2<MG
        gpex=getkappa(upe,k_eq)*(phi(j,i+2)-phi(j,i+1))/gs + ...
        (1-getkappa(upe,k_eq))*(phi(j,i+1)-phi(j,i))/gs;
     else
         gpex=(phi(j,i+1)-phi(j,i))/gs;
     end
     gpey=getkappa(vpe,k_eq)*(phi(j+1,i+1)-phi(j,i+1))/gs + ...
                (1-getkappa(vpe,k_eq))*(phi(j,i+1)-phi(j-1,i+1))/gs;
     if i-2>0
     gpwx=getkappa(upw,k_eq)*(phi(j,i)-phi(j,i-1))/gs + ...
        (1-getkappa(upw,k_eq))*(phi(j,i-1)-phi(j,i-2))/gs;
     else
     gpwx=(phi(j,i)-phi(j,i-1))/gs;
     end
     gpwy=getkappa(vpw,k_eq)*(phi(j+1,i-1)-phi(j,i-1))/gs + ...
                (1-getkappa(vpw,k_eq))*(phi(j,i-1)-phi(j-1,i-1))/gs;
     gpnx=getkappa(upn,k_eq)*(phi(j+1,i+1)-phi(j+1,i))/gs + ...
        (1-getkappa(upn,k_eq))*(phi(j+1,i)-phi(j+1,i-1))/gs;
     if j+2<MG
     gpny=getkappa(vpn,k_eq)*(phi(j+2,i)-phi(j+1,i))/gs + ...
                (1-getkappa(vpn,k_eq))*(phi(j+1,i)-phi(j,i))/gs;
     else
     gpny=(phi(j+1,i)-phi(j,i))/gs;  
     end
     gpsx=getkappa(ups,k_eq)*(phi(j-1,i+1)-phi(j-1,i))/gs + ...
        (1-getkappa(ups,k_eq))*(phi(j-1,i)-phi(j-1,i-1))/gs;
    if j-2>0
     gpsy=getkappa(vps,k_eq)*(phi(j,i)-phi(j-1,i))/gs + ...
                (1-getkappa(vps,k_eq))*(phi(j-1,i)-phi(j-2,i))/gs;
    else
     gpsy=(phi(j,i)-phi(j-1,i))/gs;
    end
            
     if ue>0
         phie=phi(j,i)+gs/2*gpx-dt/2*(gpx*up+gpy*vp);
     else
         phie=phi(j,i+1)-gs/2*gpex-dt/2*(gpex*upe+gpey*vpe);
     end
     
     if uw>0
         phiw=phi(j,i-1)+gs/2*gpwx-dt/2*(gpwx*upw+gpwy*vpw);
     else
         phiw=phi(j,i)-gs/2*gpx-dt/2*(gpx*up+gpy*vp);
     end
     
     if vn>0
         phin=phi(j,i)+gs/2*gpy-dt/2*(gpx*up+gpy*vp);
     else
         phin=phi(j+1,i)-gs/2*gpny-dt/2*(gpnx*upn+gpny*vpn);
     end
     
     if vs>0
         phis=phi(j-1,i)+gs/2*gpsy-dt/2*(gpsx*ups+gpsy*vps);
     else
         phis=phi(j,i)-gs/2*gpy-dt/2*(gpx*up+gpy*vp);
     end
     
     phi_tay=phi(j,i)-dt/2*(gpx*up+gpy*vp);
            
    
     phiN(j,i)=phi(j,i)-dt/gs*(ue*(phie-phi_tay)+vn*(phin-phi_tay)) + ...
         dt/gs*(uw*(phiw-phi_tay)+vs*(phis-phi_tay));
   end
end
%Boundary conditions
phiN(:,1)=phiN(:,2);
phiN(:,end)=phiN(:,end-1);
phiN(1,:)=phiN(2,:);
phiN(end,:)=phiN(end-1,:);
phi=phiN;

if mod(tt,200)==0 %records the phi at a certain time
    scount=scount+1;
    SS(scount).t=dt*tt;
    SS(scount).phi=phi;
    h=figure('Visible','off');
    contour(xm,ym,phi);
    hold on
    contour(xm,ym,phi,[0,0],'LineWidth',2);
    hold off
    axis([0 1 0 1])
    axis equal
    xlabel('X','Interpreter','LaTex');
    ylabel('Y','Interpreter','LaTex');
    title(['t= ' num2str(tt*dt) '\hspace{0.1in} N=' num2str(MG) 'x' num2str(MG) ...
    ' (FVM)'],'Interpreter','LaTex');
    c=colorbar;
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',16,'FontName','AvantGarde'); 
    set(c,'fontsize',12);
    saveas(h,[pwd sf_dr sprintf('FIG%d.png',scount)]);
    close all
end

end
matfile = fullfile([pwd fileparts(sf_dr)], sprintf('TEST%d.mat',k));
save(matfile);
end
