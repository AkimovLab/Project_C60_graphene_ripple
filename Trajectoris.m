% This code plots the trajectories relating to motion of C60 at different temperatures
clc, clear
load allT2s.mat

% Data column order:
%1Step 2CPU 3PotEng 4KinEng 5Temp 6Lx 7Ly 8Press 
%9v_xc_x 10v_xc_y 11v_xc_z 12c_pe_c60 13c_lennard 14c_ke_c60 
%15v_vc_x 16v_vc_y 17v_vc_z
%18v_x1_x 19v_x1_y 20v_x1_z 21v_x2_x 22v_x2_y 23v_x2_z 24c_pe_sub 25c_ke_sub 
%26v_wc_x 27v_wc_y 28v_wc_z 29v_w12_x 30v_w12_y 31v_w12_z 32c_temp_c60 33c_temp_sub

%T=[1,2,3, 4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21  ];
T= [1,5,10,20,30,35,50,60,75,100,150,200,250,300,400,500,600,700,800,900,1000];

Ti=14;
% P=100e3;                      % lag length (time) should be up to N/4 (less than 50 K or =47.5 K)
dt = 1e-3;                      % 0.001 ps time between trajectory points
thermo=200;
% tk=[1:P]'*dt*thermo/1000;            % ns
x(:,:)=imdata(:,9,:);
y(:,:)=imdata(:,10,:);
z(:,:)=imdata(:,11,:);
len(:,:)=imdata(10000:end,13,:);
pe(:,:)=imdata(:,12,:);
potE(:,:)=imdata(:,3,:);
NT = length(x); 
t = imdata(1:NT,1);
time=dt*t;

%% Trajectory
trans=100;
xt=x;yt=y;
xt(:,1:3)=x(:,1:3)-500;
%  xt(:,4)=x(:,4)-500;
%  yt(:,6)=y(:,6)+50;
 yt(:,1)=y(:,1)+500;
 yt(:,2)=y(:,2)+400;
 yt(:,3)=y(:,3)+300;
 yt(:,4)=y(:,4)-200;

figure(1)
for j=1:6
% Ekave(j)=mean(imdata(10:end,6,j));
% Etot=imdata(:,3,j)+imdata(:,6,j);
% Etot2(j)=mean(imdata(10:end,9,j));
% plot(imdata(:,1,j),Etot);
plot(xt(:,j),yt(:,j),'LineWidth',3); 
hold on
end
xlabel('x (\AA)','Interpreter','latex')
ylabel('y (\AA)','Interpreter','latex')
legend('T=1 K','T=5 K','T=10 K','T=20 K','T=30 K','T=35 K')
set(gca,'FontName','Cambria','FontSize',16);

xt=x;yt=y;
 xt(:,6:8)=x(:,6:8)-1000;
% xt(:,9)=x(:,9)+700;
% 
 yt(:,6)=y(:,6)+800;
 yt(:,7)=y(:,7)+700;
yt(:,8)=y(:,8)+400;
 yt(:,9)=y(:,9)-1300;
 yt(:,10)=y(:,10)-1400;
%  xt(:,10)=x(:,10)+1000;

figure(2)
for j=6:11
plot(xt(:,j),yt(:,j),'LineWidth',3);
hold on
end
xlabel('x (\AA)','Interpreter','latex')
ylabel('y (\AA)','Interpreter','latex')
legend('T=35 K','T=50 K','T=60 K','T=75 K','T=100 K','T=150 K')
set(gca,'FontName','Cambria','FontSize',16);

xt=x;yt=y;
yt(:,11:13)=y(:,11:13)-5000;
% yt(:,12:13)=y(:,12:13)-6000;
 yt(:,15)=y(:,15)+1000;
%  xt(:,15)=x(:,15)+3000;
%  xt(:,11)=x(:,11)-1000;
 xt(:,12)=x(:,12)+2000;
 xt(:,13)=x(:,13)+2000;
% xt(:,14)=x(:,14)-2000;
yt(:,16)=y(:,16)+8000;

figure(3)
for j=11:16
plot(xt(:,j),yt(:,j),'LineWidth',3); 
hold on
end
xlabel('x (\AA)','Interpreter','latex')
ylabel('y (\AA)','Interpreter','latex')
legend('T=150 K','T=200 K','T=250 K','T=300 K','T=400 K','T=500 K')
set(gca,'FontName','Cambria','FontSize',16);

%% x-t levy flight
figure(16)
hold on
plot(time/1000,x(:,Ti),'.'); 
xlabel('t (ns)','Interpreter','latex')
ylabel('x (\AA)','Interpreter','latex')
% legend('T=30 K','T=300 K')
% legend('C60/mono-Gr','C60/2layer-Gr','C60/fixed-Gr','Location','Best')
set(gca,'FontName','Cambria','FontSize',16);

% z-t
figure(15)
hold on
plot(time-20000,z(:,Ti),'.'); 
xlabel('t (ps)','Interpreter','latex')
ylabel('z (\AA)','Interpreter','latex')
% legend('T=30 K','T=300 K')
% legend('C60/SLG','C60/DLG','C60/FLG','Location','Best')
set(gca,'FontName','Cambria','FontSize',20);
axis([0,500,4.5,8])

