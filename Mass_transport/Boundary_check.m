function [out] = Boundary_check(y_final,zmin,xmin,xmax,ymin,ymax)

%.... Bottom boundary check ....... 
if y_final(3,1)<zmin
    y_final(3,1)=zmin-y_final(3,1);
end

%.... first ends of lateral boundary check
if y_final(1,1)<xmin
    clmin=abs(y_final(1,1)-xmin);
    y_final(1,1)=xmax-clmin;
end

if y_final(2,1)<ymin
    dlmin=abs(y_final(2,1)-ymin);
    y_final(2,1)=ymax-dlmin;
end

%.... last ends of lateral boundary check

if y_final(1,1)>xmax
    clmax=abs(y_final(1,1)-xmax);
    y_final(1,1)=xmin+clmax;
end

if y_final(2,1)>ymax
    dlmax=abs(y_final(2,1)-ymax);
    y_final(2,1)=ymin+dlmax;
end

out=y_final;

end