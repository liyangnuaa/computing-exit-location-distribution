function y=fun(~,x)

global X Y Z Ux Uy Uz xlin ylin zlin

[~,Ix]=min(abs(xlin-x(1)));
[~,Iy]=min(abs(ylin-x(2)));
[~,Iz]=min(abs(zlin-x(3)));
J=3;
X0=X(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);
Y0=Y(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);
Z0=Z(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);
Ux0=Ux(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);
Uy0=Uy(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);
Uz0=Uz(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);

bx=[-2*(x(1)^3-x(1))-(x(2)+x(3));
    -x(2)+2*(x(1)^3-x(1));
    -x(3)+2*(x(1)^3-x(1))];
y=zeros(size(x));
y(1)=-(bx(1)+interp3(X0,Y0,Z0,Ux0,x(1),x(2),x(3),'linear'))/norm(bx);
y(2)=-(bx(2)+interp3(X0,Y0,Z0,Uy0,x(1),x(2),x(3),'linear'))/norm(bx);
y(3)=-(bx(3)+interp3(X0,Y0,Z0,Uz0,x(1),x(2),x(3),'linear'))/norm(bx);

