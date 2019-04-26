function dangle = state(t,angle,omega)

a=angle(1); b=angle(2); c=angle(3);
Wx=omega(1); Wy=omega(2); Wz=omega(3);


if abs(cos(b))<0.01
    cosb=sign(cos(b))*0.01;
else
    cosb=cos(b);
end

% R=Rz*Ry*Rx;
dangle1=(Wx*cos(c) + Wy*sin(c))/cosb;
dangle2= Wy*cos(c) - Wx*sin(c);
dangle3= Wz + sin(b)*dangle1;

% R=Rx*Ry*Rz;
% dangle3=(Wz*cos(a) - Wy*sin(a))/cosb;
% dangle1=Wx - sin(b)*dangle3;
% dangle2=Wy*cos(a) + Wz*sin(a);


dangle=[dangle1; dangle2; dangle3];