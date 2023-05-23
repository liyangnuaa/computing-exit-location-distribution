function y=fun(~,x)

global miu gamma X Y Ux Uy

xlin=X(1,:);
ylin=Y(:,1);
[~,Ix]=min(abs(xlin-x(1)));
[~,Iy]=min(abs(ylin-x(2)));
J=5;
X0=X(Iy-J:Iy+J,Ix-J:Ix+J);
Y0=Y(Iy-J:Iy+J,Ix-J:Ix+J);
Ux0=Ux(Iy-J:Iy+J,Ix-J:Ix+J);
Uy0=Uy(Iy-J:Iy+J,Ix-J:Ix+J);

bx=[x(1)-x(1)^3-gamma*x(1)*x(2)^2;-miu*(1+x(1)^2)*x(2)];
y=zeros(size(x));
y(1)=-(bx(1)+interp2(X0,Y0,Ux0,x(1),x(2),'linear'))/norm(bx);
y(2)=-(bx(2)+interp2(X0,Y0,Uy0,x(1),x(2),'linear'))/norm(bx);

