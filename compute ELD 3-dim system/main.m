clear;
clc;

global X Y Z Ux Uy Uz xlin ylin zlin

xmin=-2;
xmax=0;
ymin=-1.5;
ymax=1.5;
zmin=-1.5;
zmax=1.5;

%% vector field
% n0=20;
% x1=linspace(xmin,xmax,n0);
% y1=linspace(ymin,ymax,n0);
% [x,y]=meshgrid(x1,y1);
% u=x-x.^3-gamma*x.*y.^2;
% v=-miu*(1+x.^2).*y;
% figure;
% streamslice(x,y,u,v);
% axis([xmin xmax ymin ymax])
% figure;
% plot([-0.5,-0.5],[-0.8,0.8]);
% figure;
% plot(-1,0,'*');

%% data processing
N=100;
xlin=linspace(xmin,xmax,N);
ylin=linspace(ymin,ymax,N);
zlin=linspace(zmin,zmax,N);
[xmesh,ymesh,zmesh]=meshgrid(xlin,ylin,zlin);

xtest=reshape(xmesh,1,N^3);
ytest=reshape(ymesh,1,N^3);
ztest=reshape(zmesh,1,N^3);
path = sprintf('xtest.mat');
save(path,'xtest');
path = sprintf('ytest.mat');
save(path,'ytest');
path = sprintf('ztest.mat');
save(path,'ztest');

%% NN results
load('Stest.mat');
xmesh1=reshape(xtest,N,N,N);
ymesh1=reshape(ytest,N,N,N);
zmesh1=reshape(ztest,N,N,N);
% xmesh1=permute(xmesh1,[2 1 3]);
% ymesh1=permute(ymesh1,[2 1 3]);
% zmesh1=permute(zmesh1,[2 1 3]);
S=Stest(:,1)';
S=S+(xtest+1).^2+ytest.^2+ztest.^2;
lx=Stest(:,2)';
ly=Stest(:,3)';
lz=Stest(:,4)';
Smesh1=reshape(S,N,N,N);
lxmesh1=reshape(lx,N,N,N);
lymesh1=reshape(ly,N,N,N);
lzmesh1=reshape(lz,N,N,N);
% Smesh1=permute(Smesh1,[2 1 3]);
% lxmesh1=permute(lxmesh1,[2 1 3]);
% lymesh1=permute(lymesh1,[2 1 3]);
% lzmesh1=permute(lzmesh1,[2 1 3]);

ymesh2=reshape(ymesh1(:,75,:),N,N);
zmesh2=reshape(zmesh1(:,75,:),N,N);
Smesh2=reshape(Smesh1(:,75,:),N,N);
lxmesh2=reshape(lxmesh1(:,75,:),N,N);
lymesh2=reshape(lymesh1(:,75,:),N,N);
lzmesh2=reshape(lzmesh1(:,75,:),N,N);

figure;
mesh(ymesh2,zmesh2,Smesh2);
Strue=9/16+ymesh2.^2+zmesh2.^2;
figure;
mesh(ymesh2,zmesh2,Strue);


%% boundary points
Nbp=100;
dellt00=0.7;
y1=linspace(ymin+dellt00,ymax-dellt00,Nbp);
z1=linspace(zmin+dellt00,zmax-dellt00,Nbp);
[yy,zz]=meshgrid(y1,z1);
xx=-0.5*ones(Nbp,Nbp);
xb=[reshape(xx,1,Nbp^2);reshape(yy,1,Nbp^2);reshape(zz,1,Nbp^2)];
INT=zeros(1,Nbp^2);
V_star=zeros(1,Nbp^2);
V1=zeros(1,Nbp^2);

%% Computing prefactor via NN
xnode=[-1;0;0];
M1=[8,2,2,0,0,0;
    -4,5,0,1,1,0;
    -4,0,5,0,1,1;
    0,-8,0,2,0,0;
    0,-4,-4,0,2,0;
    0,0,-8,0,0,2];
M2=[1;0;0;1;0;1];
A=M1\M2;
H_bar=inv([A(1),A(2),A(3);A(2),A(4),A(5);A(3),A(5),A(6)]);
h=0.005;
T=100;
Nstep=floor(T/h);
dx=xlin(2)-xlin(1);
dy=ylin(2)-ylin(1);
dz=zlin(2)-zlin(1);
% xb=[-0.5;0];
% Nbp=1;

X=xmesh1;
Y=ymesh1;
Z=zmesh1;
[Ux,Uy,Uz]=gradient(Smesh1,dx,dy,dz);
[L1_1,L1_2,L1_3]=gradient(lxmesh1,dx,dy,dz);
[L2_1,L2_2,L2_3]=gradient(lymesh1,dx,dy,dz);
[L3_1,L3_2,L3_3]=gradient(lzmesh1,dx,dy,dz);
for k=1:Nbp^2
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
        [~,Iz]=min(abs(zlin-x(3)));
        J=3;
        X0=X(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);
        Y0=Y(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);
        Z0=Z(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);
        L1_10=L1_1(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);
        L2_20=L2_2(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);
        L3_30=L3_3(Iy-J:Iy+J,Ix-J:Ix+J,Iz-J:Iz+J);
        divL=[interp3(X0,Y0,Z0,L1_10,x(1),x(2),x(3),'linear')+...
            interp3(X0,Y0,Z0,L2_20,x(1),x(2),x(3),'linear')+...
            interp3(X0,Y0,Z0,L3_30,x(1),x(2),x(3),'linear'),divL];
        normbx=[norm([-2*(x(1)^3-x(1))-(x(2)+x(3));-x(2)+2*(x(1)^3-x(1));-x(3)+2*(x(1)^3-x(1))]),normbx];
    end
    integ=(sum(divL./normbx)-0.5*(divL(1)/normbx(1)+divL(end)/normbx(end)))*h;
    INT(k)=integ;
    
    [~,Ix]=min(abs(xlin-xb(1,k)));
    [~,Iy]=min(abs(ylin-xb(2,k)));
    [~,Iz]=min(abs(zlin-xb(3,k)));
    V_star(k)=Smesh1(Iy,Ix,Iz);
    V1(k)=Ux(Iy,Ix,Iz);
end

%% computing exit location distribution
epsilon=1/20;
prefactor=sqrt(abs(det(H_bar)))/(2*pi*epsilon)*exp(-INT);
ELD_NN=V1.*prefactor.*exp(-V_star/epsilon);
y0=reshape(xb(2,:),Nbp,Nbp);
z0=reshape(xb(3,:),Nbp,Nbp);
ELD_NN0=reshape(ELD_NN,Nbp,Nbp);
elds=0;
for iz=1:Nbp-1
    for iy=1:Nbp-1
        elds=elds+0.25*(y1(2)-y1(1))*(z1(2)-z1(1))*(ELD_NN0(iz,iy)+ELD_NN0(iz+1,iy)+ELD_NN0(iz,iy+1)+ELD_NN0(iz+1,iy+1));
    end
end
ELD=ELD_NN0/elds;

figure;
mesh(y0,z0,ELD);


