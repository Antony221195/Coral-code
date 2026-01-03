function [centres_x,centres_y,xhex_cen,yhex_cen,mod_angles,ang_diff] = Hexagon_circle_packing(super_radius,Patch_radius,xrow,yrow,r,sigma)

%% General description
%.... This files build the circle packing in hexagon and the centres of
%each circle is rotlet. This also gives you the angles of each centres with
%respect to each origin of the hexagon in the given domain

%% The following will trim down the hexagon centres from incoming array

[xbuild,ybuild]=meshgrid(xrow,yrow);

%...... Filtering hexagonal centres..........

xbuild_tran=xbuild';
ybuild_tran=ybuild';

xbuild_col=xbuild_tran(:);
ybuild_col=ybuild_tran(:);

xhex_cen=NaN(10000,1);   %..... x position of centres
yhex_cen=NaN(10000,1);   %..... y position of centres
count=1;
for iii=1:length(xbuild_col)
    if mod(iii,2)==1
        xhex_cen(count,1)=xbuild_col(iii,1);
        yhex_cen(count,1)=ybuild_col(iii,1);
        count=count+1;
    end
end
xhex_cen(isnan(xhex_cen))=[];
yhex_cen(isnan(yhex_cen))=[];

xhex_cen=xhex_cen'; yhex_cen=yhex_cen';

%% The following will generate the centres of each rotlet in single hexagon patch of radius 'Patch_radius'

ysvec=0:super_radius*sqrt(3):(sqrt(3)*Patch_radius)/2;
xs=cell(length(ysvec),1);
ys=cell(length(ysvec),1);

for iii=1:length(ysvec)-1
    xsvec=-(iii-1)*super_radius:2*super_radius:(iii-1)*super_radius;
    xs{iii,1}=xsvec';
    yv=ysvec(iii)*ones(1,length(xsvec));
    ys{iii,1}=yv';
end

xs=cell2mat(xs);
ys=sqrt(3)*super_radius+cell2mat(ys);

rot_vec=-60*(1:1:5); %.... different angles of each face of the hexagon. Here 5 is remaining faces

%..... Centres in other faces orienting in clockwise direction .........
xs_ph=cosd(rot_vec).*xs-sind(rot_vec).*ys;
ys_ph=sind(rot_vec).*xs+cosd(rot_vec).*ys;

%...... Gives the centres of rotlet in hexagon positioned at origin
xcentres_rot=[xs;xs_ph(:)];
ycentres_rot=[ys;ys_ph(:)];

%% Evaluating the angle extended by the radial vector with respect to positive x-axis for each super rotlets in origin hexagone

origin_angles=atan2(ycentres_rot,xcentres_rot)*(180/pi);
origin_angles(origin_angles<0)=origin_angles(origin_angles<0)+360;

origin_angles_het=zeros(size(origin_angles));
ang_diff=zeros(size(origin_angles));

%.... Initialise the random variables
rng(1)
ran=randn(1,length(origin_angles_het));

for www=1:length(origin_angles_het)
    ranvar=ran(www);
    origin_angles_het(www,1)=mod(origin_angles(www,1) + sigma * ranvar, 360);  %.... establishing heterogeneity in alignment
    angle_wrapped_het = mod(origin_angles_het(www,1) + 180, 360) - 180;               %..... Wrap to [-180, 180) for het
    angle_wrapped = mod(origin_angles(www,1) + 180, 360) - 180;               %..... Wrap to [-180, 180) for ref angle
    ang_diff(www,1)=angle_wrapped-angle_wrapped_het;
end


%% The following will translate the rotlet in origin hexagon to other centres' hexagon
centres_x=xcentres_rot+xhex_cen;
centres_y=ycentres_rot+yhex_cen;

%% Generating angle extended by each centres with respective origin
mod_angles=zeros(size(centres_x));
for ttt=1:length(xhex_cen)
    mod_angles(:,ttt)=origin_angles_het;
end

%% Building the rectangular domain to eliminate the rotlet outside the domain and image layer

%%... Condition on image layer as given Appendix C of manuscript

xstart=xrow(2)-0.5*Patch_radius; xend=xrow(end-1)+0.5*Patch_radius;
ystart=yrow(2)-0.5*r; yend=yrow(end-1)+0.5*r;

x_rect=[xstart xend xend xstart xstart];
y_rect=[ystart ystart yend yend ystart];

%..... columnize the centres and angles.....
centres_x=centres_x(:);
centres_y=centres_y(:);
mod_angles=mod_angles(:);

[in,~] = inpolygon(centres_x,centres_y,x_rect,y_rect);

centres_x=centres_x(in);
centres_y=centres_y(in);
mod_angles=mod_angles(in);


end





