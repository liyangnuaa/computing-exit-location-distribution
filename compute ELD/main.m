clear;
clc;

global miu gamma X Y Ux Uy

xmin=-1.5;
xmax=0;
ymin=-0.8;
ymax=0.8;

miu=1.5;
gamma=1;

%% vector field
n0=20;
x1=linspace(xmin,xmax,n0);
y1=linspace(ymin,ymax,n0);
[x,y]=meshgrid(x1,y1);
u=x-x.^3-gamma*x.*y.^2;
v=-miu*(1+x.^2).*y;
figure;
streamslice(x,y,u,v);
axis([xmin xmax ymin ymax])
figure;
plot([-0.5,-0.5],[-0.8,0.8]);
figure;
plot(-1,0,'*');

%% data processing
N=500;
xlin=linspace(xmin,xmax,N);
ylin=linspace(ymin,ymax,N);
[xmesh,ymesh]=meshgrid(xlin,ylin);

xtest=reshape(xmesh,1,N^2);
ytest=reshape(ymesh,1,N^2);
path = sprintf('xtest.mat');
save(path,'xtest');
path = sprintf('ytest.mat');
save(path,'ytest');

%% NN results
load('Stest.mat');
xmesh1=reshape(xtest,N,N);
ymesh1=reshape(ytest,N,N);
S=Stest(:,1)';
S=S+(xtest+1).^2+ytest.^2;
Smesh1=reshape(S,N,N);
lx=Stest(:,2)';
lxmesh1=reshape(lx,N,N);
ly=Stest(:,3)';
lymesh1=reshape(ly,N,N);

figure;
mesh(xmesh1,ymesh1,Smesh1);
figure;
mesh(xmesh1,ymesh1,lxmesh1);
figure;
mesh(xmesh1,ymesh1,lymesh1);


%% boundary points
Nbp=200;
xb=[-0.5*ones(1,Nbp);linspace(ymin+0.02,ymax-0.02,Nbp)];
INT=zeros(1,Nbp);
V_star=zeros(1,Nbp);
V1=zeros(1,Nbp);

%% Computing prefactor via NN
xnode=[-1;0];
H_bar=[4,0;0,4*miu];
h=0.001;
T=10;
Nstep=floor(T/h);
dx=xlin(2)-xlin(1);
dy=ylin(2)-ylin(1);
% xb=[-0.5;0];
% Nbp=1;

X=xmesh1;
Y=ymesh1;
[Ux,Uy]=gradient(Smesh1,dx,dy);
[L1_1,L1_2]=gradient(lxmesh1,dx,dy);
[L2_1,L2_2]=gradient(lymesh1,dx,dy);
for k=1:Nbp
    delta=0.05;
    t=-h;
    x0=xb(:,k);
    MPEP=x0;
    
    for i=1:Nstep
        x1=rk4(0,h,x0);
        if norm(x1-xnode)<=0.2*delta
            break;
        end
        x0=x1;
        MPEP=[x0 MPEP];
        t=[-(i+1)*h,t];
    end
    %MPEP=[xnode MPEP];
%     figure;
%     plot(MPEP(1,:),MPEP(2,:),'m-');
    
    Np=length(t);
    divL=[];
    normbx=[];
    for i=1:Np
        x=MPEP(:,i);
        [~,Ix]=min(abs(xlin-x(1)));
        [~,Iy]=min(abs(ylin-x(2)));
        J=5;
        X0=X(Iy-J:Iy+J,Ix-J:Ix+J);
        Y0=Y(Iy-J:Iy+J,Ix-J:Ix+J);
        L1_10=L1_1(Iy-J:Iy+J,Ix-J:Ix+J);
        L2_20=L2_2(Iy-J:Iy+J,Ix-J:Ix+J);
        divL=[interp2(X0,Y0,L1_10,x(1),x(2),'linear')+interp2(X0,Y0,L2_20,x(1),x(2),'linear'),divL];
        normbx=[norm([x(1)-x(1)^3-gamma*x(1)*x(2)^2;-miu*(1+x(1)^2)*x(2)]),normbx];
    end
    integ=(sum(divL./normbx)-0.5*(divL(1)/normbx(1)+divL(end)/normbx(end)))*h;
    INT(k)=integ;
    
    [~,Ix]=min(abs(xlin-xb(1,k)));
    [~,Iy]=min(abs(ylin-xb(2,k)));
    V_star(k)=Smesh1(Iy,Ix);
    V1(k)=Ux(Iy,Ix);
end

%% computing exit location distribution
epsilon=1/20;
prefactor=sqrt(abs(det(H_bar)))/(2*pi*epsilon)*exp(-INT);
ELD_NN=V1.*prefactor.*exp(-V_star/epsilon);
ELD=ELD_NN/(sum(ELD_NN)-0.5*(ELD_NN(1)+ELD_NN(end)))/(xb(2,2)-xb(2,1));

figure;
plot(xb(2,:),ELD);


