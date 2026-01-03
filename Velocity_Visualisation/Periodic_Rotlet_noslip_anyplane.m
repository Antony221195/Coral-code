function [u1,u2,u3,Z]=Periodic_Rotlet_noslip_anyplane(xvec,yvec,zvec,Xs,Omega,wall)

%%.... This function file evaluates the side view velocity profiles u1, u2,
%%u3 along any y-plane. The orientation of plane is decided outside the
%%function file

%%... Refer to Blake and Chwang (1974) and Selvan et al. (2023) for the
%%expression for rotlet

%%.... Getting the size of x and y via meshpoints generated outside

[m,n]=size(xvec);

%%... Initialising velocities
u1=zeros(m,n);
u2=zeros(m,n);
u3=zeros(m,n);
Z=zeros(m,n);

for jjj=1:n
    for iii=1:m
        x=xvec(iii,jjj);  y=yvec(iii,jjj);   z=zvec(iii,jjj);

        A=[0 -Omega(3) Omega(2);Omega(3) 0 -Omega(1);-Omega(2) Omega(1) 0];
        Bk3=[Omega(2);-Omega(1);0];

        up_xpos=Xs(1); up_ypos=Xs(2);

        X1=[x,y,z];  %.... Evalation point

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
        %......... Higher order singularities 
        C11=(1/R^3)-((3*R1*R1)/R^5); C12=-(3*R1*R2)/R^5; C13=-(3*R1*R3)/R^5;
        C21=-(3*R2*R1)/R^5; C22=(1/R^3)-((3*R2*R2)/R^5); C23=-(3*R2*R3)/R^5;
        C31=-(3*R3*R1)/R^5; C32=-(3*R3*R2)/R^5; C33=(1/R^3)-((3*R3*R3)/R^5);
        Cik=[C11 C12 C13;C21 C22 C23;C31 C32 C33];
        U2im=2*Xs(3)*(Cik*Bk3);

        D11=(R1*R1*R3)/R^5;D12=(R1*R2*R3)/R^5;D13=(R1*R3*R3)/R^5;
        D21=(R2*R1*R3)/R^5;D22=(R2*R2*R3)/R^5;D23=(R2*R3*R3)/R^5;
        D31=(R3*R1*R3)/R^5;D32=(R3*R2*R3)/R^5;D33=(R3*R3*R3)/R^5;
        Dik=[D11 D12 D13;D21 D22 D23;D31 D32 D33];
        U3im=6*Dik*Bk3;

        %.........Final solution for advection term .............
        sum=((1/(8*pi))*(U1-wall*U1im+wall*U2im+wall*U3im));

        u1(iii,jjj)=sum(1,1);
        u2(iii,jjj)=sum(2,1);
        u3(iii,jjj)=sum(3,1);
        Z(iii,jjj)=z;
    end
end

end
