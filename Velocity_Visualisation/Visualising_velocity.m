clear       %% Clears the memory
clc         %% Clear the command window of the matlab
close all   %% Close all the exisitng files opened in matlab

%% Physical parameters for cilia

lc=12.75;               %..... Ciliary length in micrometer (um)
Tc=1/16.9;              %..... Beating period (s)
Ut=(2*pi*lc)/Tc;        %..... Tip velocity (um/s)
ref_len=1.2595*(lc/2);  %..... Position of rotlet (i.e., point torque) above the no-slip wall (reference scaling for evaluating the flow field)

%% Building the hexagonal grid

R=11*ref_len;                 %..... Hexagonal outer radius (in um)
r=(sqrt(3)*R)/2;              %..... Hexagonal inner radius (in um)
wall=1;                       %..... Switching on the no-slip wall contribution on velocity (1 for on/ 0 for off)
Ommag=10*18.384*ref_len^3;    %..... Super rotlet magnitude obtained from Appendix A
super_radius=1*ref_len;       %..... Optimum radius obtained from Appendix A

N=3;        %.... Total number of hexagonal centre in x-direction (in both domain and image layer)
M=3;        %.... Total number of hexagonal centre in y-direction (in both domain and image layer) 

%% Building the centres in terms of inner and outer hex. radius

xrow=(-(N-1):1:N-1)*(1.5*R);
yrow=(-(M-1):1:M-1)*r;

%% Introducing orientational noise

sigma=0;   %.....Possible values {0, 20, 30}

%% Obtaining point torque positions in xy planes, hexagonal centres, torque angles, difference between the mean angle and disturbed angle of torque

[centres_x,centres_y,xhex_cen,yhex_cen,mod_angles, ang_diff] = Hexagon_circle_packing(super_radius,R,xrow,yrow,r,sigma);

%% Evaluating velocity profiles

thres=250;  %..... Extra points in domain 

xvec=linspace(xrow(2)-thres,xrow(end-1)+thres,100);   %.... Query points in x-direction
yvec=linspace(yrow(2)-thres,yrow(end-1)+thres,100);   %.... Query points in y-direction
zvec=linspace(0,10*ref_len,150);                      %.... Query points in z-direction

[X, Z] = meshgrid(xvec, zvec);                        %.... Obtaining meshgrid in xz-direction

%%... This rotation is majorly done for side view, where we fix the
%%y-aixs point and rotate the xz-plane with angle 'theta'

Y=0*ones(size(X));                                    %......Evaluation plane in y-direction, which could be rotated .......
theta=30;                                             %...... Rotation angle

Xd=cosd(theta)*X-sind(theta)*Y;                       %...... x-plane (meshpoints) after rotation by angle 'theta'
Yd=sind(theta)*X+cosd(theta)*Y;                       %...... y-plane (meshpoints) after rotation by angle 'theta'
Zd=Z;                                                 %...... Rotation in 2D, hence z-axis is fixed

%...... Rotation with non-mesh points in x and y direction

xvecd=cosd(theta)*xvec-sind(theta)*Y(1,:);            
yvecd=sind(theta)*xvec+cosd(theta)*Y(1,:);

[mm,nn]=size(Xd);   %...... obtaining the number of points along x and y direction

%...... Measuring the lengths of domain along x and y direction

Lx=abs(xrow(1)-xrow(end));
Ly=abs(yrow(1)-yrow(end));

%...... Initializing the velocity profiles in xy-direction

sumx_xy=zeros(length(yvec),length(xvec));
sumy_xy=zeros(length(yvec),length(xvec));
sumz_xy=zeros(length(yvec),length(xvec));

%...... Initializing the velocity profiles in xz-direction

sumx_xz=zeros(mm,nn);
sumy_xz=zeros(mm,nn);
sumz_xz=zeros(mm,nn);

%%.... Evaluating on each points 

for iii=1:length(centres_x)   %..... Running over each points                                
  
    theta=mod_angles(iii,1);     %....... Getting the reference angle with respect to each circle centre

    Rotat=[cosd(theta) -sind(theta) 0;sind(theta) cosd(theta) 0;0 0 0]; %..... Rotation matrix to orient the torque radially outward
    vec=[0;1;0];               %..... Vector orienting torque such that flow will be in +x-direction
    res_vec=Rotat*vec;         %..... Multiplying with Rotational matrix gives the radially outward oriented torque
    Omega=Ommag*res_vec';     %..... Resulting torque with magnitude

    Xs=[centres_x(iii,1) centres_y(iii,1) 1*ref_len];   %..... Position of each torque, which is the centre of circles packing the hexagon

    [u1_xy,u2_xy,u3_xy]=Periodic_Rotlet_noslip_xyplane(xvec,yvec,2*ref_len,Xs,Omega,wall,Lx,Ly);   %.... Function file for evaluating enface view

    [u1_xz,u2_xz,u3_xz,Zm]=Periodic_Rotlet_noslip_anyplane(Xd,Yd,Zd,Xs,Omega,wall);    %...... Function file for evaluating the side view in any direction

    %%.. Superposition of velocities in the given points

    sumx_xy=sumx_xy+u1_xy;
    sumy_xy=sumy_xy+u2_xy;
    sumz_xy=sumz_xy+u3_xy;

    sumx_xz=sumx_xz+u1_xz;
    sumy_xz=sumy_xz+u2_xz;
    sumz_xz=sumz_xz+u3_xz;
end

%% Adding Background flow

ep=0;  %... Shear rate

sumx_xy1=sumx_xy+ep*(2*ref_len);   %..... Adding the shear rate in the x-component in enface view
sumx_xz1=sumx_xz+ep*Zm;            %..... Adding the shear rate in the x-component in side view

%%... Resulting magnitude

Um_xz=sqrt(sumx_xz1.^2+sumy_xz.^2+sumz_xz.^2);
Um_xy=sqrt(sumx_xy1.^2+sumy_xy.^2+sumz_xy.^2);

%% Setting up the boundaries for limiting using xlim, ylim, zlim commands 

xtic_min=xrow(2)/lc;
xtic_max=xrow(end-1)/lc;

ytic_min=yrow(2)/lc;
ytic_max=yrow(end-1)/lc;

ztic_min=zvec(1)/lc;
ztic_max=zvec(end)/lc;

%% Plotting the enface view

figure
imagesc(xvec/lc,yvec/lc,Um_xy/Ut)
ax = gca;
ax.YDir = 'normal';
hold on
colormap hot
hcb = colorbar;
hcb.Title.String = "$|\textbf{\textit{u}}|/U_t$";
hcb.Title.Interpreter = "latex";
hcb.Location="eastoutside";
oldmap=colormap;
colormap(flipud(oldmap));
h1=streamslice(xvec/lc,yvec/lc,sumx_xy1/Ut,sumy_xy/Ut,10);
set(h1,'color','k','LineStyle','-','linewidth',1.5)
hold on
xlabel('$x/\langle l_c \rangle$','fontsize',25,'interpreter','latex')
ylabel('$y/\langle l_c \rangle$','fontsize',25,'interpreter','latex')
clim([0,1])

xlim([xrow(2)/lc xrow(end-1)/lc])
ylim([yrow(2)/lc yrow(end-1)/lc])

xticks([xtic_min 0.5*(xtic_min+xtic_max) xtic_max])
yticks([ytic_min 0.5*(ytic_min+ytic_max) ytic_max])

xticklabels({0 9 18})
yticklabels({0 5 10})

set(gca,'linewidth',2,'fontsize',25)
hold on
plot(Xd(1,:)/ref_len,Yd(1,:)/ref_len,'k--','LineWidth',4)
hold on

exportgraphics(gca,'enface_sigma_5.eps','ContentType','vector') %..... Saving the file in vector format

%% Plotting the side view

figure
imagesc(xvecd/lc,zvec/lc,Um_xz/Ut)
ax = gca;
ax.YDir = 'normal';
hold on
colormap hot
% hcb = colorbar;
% hcb.Title.String = "$|\textbf{\textit{u}}|~(\mu m/s)$";
% hcb.Title.Interpreter = "latex";
oldmap=colormap;
colormap(flipud(oldmap));
h1=streamslice(xvecd(1:1:100)/lc,zvec(1:1:150)/lc,sumx_xz1(1:1:150,1:1:100)/Ut,sumz_xz(1:1:150,1:1:100)/Ut,5);
set(h1,'color','k','LineStyle','-','linewidth',1.6)
hold on
xlabel('$x/\langle l_c \rangle$','fontsize',25,'interpreter','latex')
ylabel('$z/\langle l_c \rangle$','fontsize',25,'interpreter','latex')
clim([0,1])

ylim([ztic_min ztic_max])
yticks([ztic_min 0.5*(ztic_min+ztic_max) ztic_max])
yticklabels({0 4 8})

xlim([xtic_min xtic_max])
xticks([xtic_min 0.5*(xtic_min+0.5*(xtic_min+xtic_max)) 0.5*(xtic_min+xtic_max) 0.5*(xtic_max+0.5*(xtic_min+xtic_max)) xtic_max])
xticklabels({0 5 10 15})

set(gca,'linewidth',2,'fontsize',25)
daspect([1 1 1])

exportgraphics(gca,'side_view_sigma_5.eps','ContentType','vector')  

%%
