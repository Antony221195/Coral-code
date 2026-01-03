clear       %% Clears the memory
clc         %% Clear the command window of the matlab
close all   %% Close all the exisitng files opened in matlab


%% Physical parameters for cilia

lc=12.75;               %..... Ciliary length in micrometer (um)
Tc=1/16.9;              %..... Beating period (s)
Ut=(2*pi*lc)/Tc;        %..... Tip velocity (um/s)
ref_len=1.2595*(lc/2);  %..... Position of rotlet (i.e., point torque) above the no-slip wall (reference scaling for evaluating the flow field)

%% Building the hexagonal grid (Dimension in um)

R=11*ref_len;                  %..... Outer radius of hexagon
r=(sqrt(3)*R)/2;               %..... inner radius of hexagon
wall=1;                        %..... Determining condition for wall
Ommag=10*18.384*ref_len^3;     %..... Magnitude of super rotlet
super_radius=1*ref_len;        %..... radius covered by super-rotlet
H=0.25*ref_len;                %..... height of bottom source
num_particle =1000;            %..... Number of Brownian particles in source

ep=0;  %.... shear rate of background flow

%...... Building the centres ............
N=3; M=3;
xrow=(-(N-1):1:N-1)*1.5*R;
yrow=(-(M-1):1:M-1)*r;

%...... Defining the boundary ......
xmin=xrow(2); xmax=xrow(end-1);
ymin=yrow(2); ymax=yrow(end-1);
zmin=0; 

%% Bulding the periodic hexagonal unit with circle packing, whose centres denote the position of point torques
sigma=20;   %..... Introducing noise in the orientation of torques
[centres_x,centres_y,~,~,mod_angles] = Hexagon_circle_packing(super_radius,R,xrow,yrow,r,sigma);

%% Randomly distributing the particles in the source positioned above the no-slip wall
[xpos,ypos,zpos] = random_particle_box(num_particle,xrow,yrow,H);  

D=6.5;  %.... Diffusivity in um^2/s

dt=0.01;  %..... Time step 
Tend=200; %..... Total time taken

T=0:dt:Tend;   %..... Discretising the total time
N=numel(T);    
cof_mat=sqrt((2*D*dt));   %...... factor multiplying the noise part of the particle

count1=0;  ccc=1;

for ttt=1:N-1   %%... Running over each time step

    trun=T(ttt);

    %.... Initialising the position of all particles for each time step
    xup=zeros(length(xpos),1); yup=zeros(length(xpos),1); zup=zeros(length(xpos),1);

    for uuu=1:length(xpos)     %..... Running over each particle

        y_initial=[xpos(uuu,1);ypos(uuu,1);zpos(uuu,1)];  %.... Initial position of the particle
        Ct=[randn;randn;randn];   %.... Noise part

        %%.... In general, the position evaluation splits into the
        %%deterministic part and noise part

        %%... Classic Runge-Kutte method for evaluating the position (for deterministic part)

        %..... Evaluating the k1, k2, k3 and k4 of Classic Runge-Kutte
        %Method (Refer the wikipedia)

        k1=dt*velfunc(y_initial,centres_x,mod_angles,centres_y,wall,Ommag,ref_len,ep);
        k2=dt*velfunc(y_initial+0.5*k1,centres_x,mod_angles,centres_y,wall,Ommag,ref_len,ep);
        k3=dt*velfunc(y_initial+0.5*k2,centres_x,mod_angles,centres_y,wall,Ommag,ref_len,ep);
        k4=dt*velfunc(y_initial+k3,centres_x,mod_angles,centres_y,wall,Ommag,ref_len,ep);

        y_rf=y_initial+(1/6)*(k1+2*k2+2*k3+k4); %.... Update the position of particle for deterministic part 

        y_rot=y_rf;   

        %.........Final solution after adding noise term .............
        y_out=y_rot+cof_mat*Ct;
        y_final=y_out;
        y_final=Boundary_check(y_final,zmin,xmin,xmax,ymin,ymax); %.... Checking the boundary condition outlined in the manucript for mass transport
        xup(uuu,1)=y_final(1,1); yup(uuu,1)=y_final(2,1); zup(uuu,1)=y_final(3,1);  %... Update the resultant position at given time-step for each particle
    end

    %%... Checking whether the particles positioned inside the source since
    %%we have to maintain the constant concentration inside the source

    count=0; out=0;

    for ppp=1:length(xup)
        zval=zup(ppp,1);
        if  ((zval<=H) && (zval>zmin))   %%.... Condition for checking the particle position inside the source
            count=count+1;
        else
            out=out+1;
        end
    end

    %%... Update the particles if the number of particles drops below the
    %%threshold number to maintain the constant concentration

    if count==num_particle
        xpos1=[]; ypos1=[]; zpos1=[];
    elseif count>num_particle
        xpos1=[]; ypos1=[]; zpos1=[];
    else
        [xpos1ex, ypos1ex, zpos1ex] = random_particle_box(num_particle-count,xrow,yrow,H);
        xpos1=xpos1ex(:); ypos1=ypos1ex(:); zpos1=zpos1ex(:);
    end

    xpos=[xup;xpos1];  ypos=[yup;ypos1];  zpos=[zup;zpos1];

    count1=count1+1;

    %%.... Updating the mass transport results for every 10 time-steps

    if count1==10
        filename1 = strcat('RESULTS1_',num2str(ccc),'.mat');
        time_taken=trun;
        save([saving_directory,filename1],'xpos','ypos','zpos','time_taken');
        count1=0;
        ccc=ccc+1;
    end

end



