% This code integrate (Runge-Kutta method) the angular velocity of C60, omegax omegay,
% omegaz around center of mass (body centered coordination) and gives the
% angle of rotation of the C60 in each dimension
clc
clear all

load 'allT1s.mat'

%% Integration
Ignore=1;
%T= [1,5,10,20,50,75,100,150,200,300,400,500,600,700,800];
%T=[1,2,3 ,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21  ];
T= [1,5,10,20,30,35,50,60,75,100,150,200,250,300,400,500,600,700,800,900];
i=14;                                   % 5K
% TimeI=imdata{i,1}(Ignore:end)'*1e-3;    % ps
% omegax=imdata{i,26}(Ignore:end)';       % rad/ps
% omegay=imdata{i,27}(Ignore:end)';       % rad/ps
% omegaz=imdata{i,28}(Ignore:end)';       % rad/ps
TimeI(:,:)=imdata(:,1,i)'*1e-3;          % ps
% omegax(:,:) = imdata(10000:end,26,i)';   % throw off 5% of initial data points
% omegay(:,:) = imdata(10000:end,27,i)';
% omegaz(:,:) = imdata(10000:end,28,i)';

omegax(:,:) = imdata(:,26,i)';   
omegay(:,:) = imdata(:,27,i)';
omegaz(:,:) = imdata(:,28,i)';

Omega=[omegax; omegay; omegaz];


N=length(TimeI);
h=(TimeI(end)-TimeI(1))/(N-1);
Angle=zeros(3,N);
Angle(:,1)=[0 ; 0 ; 0];

for n=1:N-1
    k1=h*state(TimeI(n),Angle(:,n),Omega(:,n));
    k2=h*state(TimeI(n)+h/2,Angle(:,n)+k1/2,(Omega(:,n)+Omega(:,n+1))/2);
    k3=h*state(TimeI(n)+h/2,Angle(:,n)+k2/2,(Omega(:,n)+Omega(:,n+1))/2);
    k4=h*state(TimeI(n)+h,Angle(:,n)+k3,Omega(:,n+1)); 
    
    Angle(:,n+1)=Angle(:,n)+(k1+2*k2+2*k3+k4)/6;
end


figure(10)
hold on
box on
plot(TimeI(1:500:end)/1000,omegax(1:500:end)*1000,'LineWidth',3,'DisplayName','Rot_x')
%  plot(TimeI(1:500:end)/1000,omegay(1:500:end)*1000,'LineWidth',3,'DisplayName','Rot_y')
plot(TimeI(1:500:end)/1000,omegaz(1:500:end)*1000,'LineWidth',3,'DisplayName','Rot_z')
xlabel('Time (ns)','FontName','Cambria')
ylabel('\omega (rad/ns)','FontName','Cambria')
set(gca,'FontName','Cambria','FontSize',22);
%legend('show','Location','Best');


axes('Position',[.7 .7 .2 .2])
box on
plot(TimeI(1:1000:end)/1000,Angle(1,1:1000:end),'LineWidth',3,'DisplayName','Rot_x')
hold on
%  plot(TimeI(1:1000:end)/1000,Angle(2,1:1000:end),'LineWidth',3,'DisplayName','Rot_y')
hold on
plot(TimeI(1:1000:end)/1000,Angle(3,1:1000:end),'LineWidth',3,'DisplayName','Rot_z')
xlabel('Time (ns)','FontName','Cambria')
ylabel('C_{60} Rotation (rad)','FontName','Cambria')
legend('show','Location','Best');
grid off
set(gca,'FontName','Cambria','FontSize',20);
% plot(TimeI(1:10:end)/1000,Angle(1,1:10:end),'Color',[rand,rand,rand],'LineWidth',2,'DisplayName','Rot_x')
% xlabel('Time (ns)','FontName','Cambria','FontSize',40)
% ylabel('p-carborane rotation (rad)','FontName','Cambria','FontSize',40)
% axis([0 8 -10  70]);
% axis([0 8 -200 400]);