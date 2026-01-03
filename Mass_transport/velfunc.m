function [out]=velfunc(y1,centres_x,mod_angles,centres_y,wall,Ommag,ref_len,ep)

x=y1(1);  y=y1(2);   z=y1(3);

usum=0;

for iii=1:length(centres_x)

    theta=mod_angles(iii,1);

    Rotat=[cosd(theta) -sind(theta) 0;sind(theta) cosd(theta) 0;0 0 0];
    vec=[0;1;0];
    res_vec=Rotat*vec;
    Omega=Ommag*res_vec';

    Xs=[centres_x(iii,1) centres_y(iii,1) 1*ref_len];

    A=[0 -Omega(3) Omega(2);Omega(3) 0 -Omega(1);-Omega(2) Omega(1) 0];
    Bk3=[Omega(2);-Omega(1);0];

    up_xpos=Xs(1); up_ypos=Xs(2);
    
    X1=[x,y,z];

    %......... unbounded rotlet..................
    r1=X1(1)-up_xpos; r2=X1(2)-up_ypos; r3=X1(3)-Xs(3);
    r=sqrt(r1^2+r2^2+r3^2);
    rk=(1/r^3)*[r1;r2;r3];
    U1=A*rk;
    %........ Image rotlet..............
    R1=r1; R2=r2; R3=X1(3)+Xs(3);
    R=sqrt(R1^2+R2^2+R3^2);
    Rk=(1/R^3)*[R1;R2;R3];
    U1im=A*Rk;
    %....................................
    C11=(1/R^3)-((3*R1*R1)/R^5); C12=-(3*R1*R2)/R^5; C13=-(3*R1*R3)/R^5;
    C21=-(3*R2*R1)/R^5; C22=(1/R^3)-((3*R2*R2)/R^5); C23=-(3*R2*R3)/R^5;
    C31=-(3*R3*R1)/R^5; C32=-(3*R3*R2)/R^5; C33=(1/R^3)-((3*R3*R3)/R^5);
    Cik=[C11 C12 C13;C21 C22 C23;C31 C32 C33];
    U2im=2*Xs(3)*(Cik*Bk3);
    %.......................................
    D11=(R1*R1*R3)/R^5;D12=(R1*R2*R3)/R^5;D13=(R1*R3*R3)/R^5;
    D21=(R2*R1*R3)/R^5;D22=(R2*R2*R3)/R^5;D23=(R2*R3*R3)/R^5;
    D31=(R3*R1*R3)/R^5;D32=(R3*R2*R3)/R^5;D33=(R3*R3*R3)/R^5;
    Dik=[D11 D12 D13;D21 D22 D23;D31 D32 D33];
    U3im=6*Dik*Bk3;

    %.........Final solution for advection term .............
    sum1=((1/(8*pi))*(U1-wall*U1im+wall*U2im+wall*U3im));
    usum=usum+sum1;
end

%... Background flow ........
Back_flow=[ep*z;0;0];

out=usum+Back_flow;

end


