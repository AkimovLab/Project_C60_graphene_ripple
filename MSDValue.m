% This code calculate the MSD value for each temperature by ensemble average 
clc, clear
load allT1.mat 

x(:,:)=imdata(10000:end,9,:);    % throw off 5% of initial data points
y(:,:)=imdata(10000:end,10,:);
NT = length(x);                  % number of data points
t = imdata(1:NT,1);       
Tc=imdata(:,32,:);
%T= [1,5,10,20,35,50,60,75,100,150,200,250,300,400,500,600,700];
%T= [1,2,3,4 ,5 ,6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21  ];
T= [1,5,10,20,30,35,50,60,75,100,150,200,250,300,400,500,600,700]; %,800,900
TT=[1,5,10,20,30,35,50,60,75,100,150,200,250,300,400,500];
%sq = sqrt(x.^2 + y.^2);
%NT = floor(nData); %# for MSD, dt should be up to 1/4 of number of data points

Npar=60;                        % datapoints in each section  %1500 % 5000 so np =40
Nseg=floor(NT/Npar);            %numver of particles   %33
dt = 1e-3;                      % 0.001 ps time between trajectory points
thermo=200;
nt=1:Nseg;
time=dt*thermo*nt;              % time in (ps)

nn1=5;
nn2=Npar/nn1;

msdk = zeros(Nseg,nn1,length(T)); msd = zeros(Nseg,length(T));
X=zeros(Npar,Nseg); Y=zeros(Npar,Nseg);
start=1000;

%% msd

for j=1:length(T)               % calculate msd for all T's
    for ts = 1:Nseg;
        for k=1:nn1
            for i = (k-1)*nn2:k*nn2-1
                msdk(ts,k,j)= msdk(ts,k,j)+((x(1+ts+i*Nseg,j)-x(1+i*Nseg,j))^2+(y(1+ts+i*Nseg,j)-y(1+i*Nseg,j))^2)/nn2;
%               msdk2(ts,k,j)=((X(i+1,ts)-X(i+1,1))^2+(Y(i+1,ts)-Y(i+1,1))^2)/nn2+msdk2(ts,k,j);
            end
        end
        msd(ts,j)=mean(msdk(ts,:,j));
    end  
end



ss=time(start:end)-time(start);  %[1:Nseg];
for j=1:length(TT)
    for k=1:nn1
        mm=msdk(start:end,k,j)-msdk(start,k,j);
        [a]=polyfit(ss',mm,1);
        ft = fittype( 'a*x', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( ft );
        opts.StartPoint = 450; %0.75
        [fitresult, gof] = fit( ss',mm, ft, opts );
        F(k,j)=fitresult.a/4;
        D(k,j)=a(1)/4;
    end
    Fmean(j)=mean(F(:,j));
    Fstd(j)=std(F(:,j));
    Dmean(j)=mean(D(:,j));
    Dstd(j)=std(D(:,j));
end



%% ploting MSD Curves

figure(1)
for p=[1:16]
plot(time,msd(:,p),'.-','LineWidth',3,'MarkerSize',12)
%M(p,1) = time(:)\msd(:,p);     % linear fitting with zero intercept
% y_est = tt*M(p);
% plot(tt, y_est, '-r')
Tnow=round(imdata(end,32:33,p))
hold on
end
xlabel('t (ps)','Interpreter','latex')
ylabel('MSD ($\AA^{2}$)','Interpreter','latex')
set(gca,'FontName','Cambria','FontSize',16);
legend('T=1 K','T=5 K','T=10 K','T=20 K','T=30 K','T=35 K','Location','Best')

figure(2)
for p=6:11
plot(time,msd(:,p),'.-','LineWidth',3,'MarkerSize',12)
hold on
end
xlabel('t (ps)','Interpreter','latex')
ylabel('MSD ($\AA^{2}$)','Interpreter','latex')
legend('T=35 K','T=50 K','T=60 K','T=75 K','T=100 K','T=150 K')
set(gca,'FontName','Cambria','FontSize',16);

figure(3)
for p=11:16
plot(time,msd(:,p),'.-','LineWidth',3,'MarkerSize',12)
hold on
end
xlabel('t (ps)','Interpreter','latex')
ylabel('MSD ($\AA^{2}$)','Interpreter','latex')
legend('T=150 K','T=200 K','T=250 K','T=300 K','T=400 K','T=500 K')
set(gca,'FontName','Cambria','FontSize',16);


%% log log msd
a=1;
b=3166;
lmsd=log10(msd);
ltime=log10(time);

for p=[2:16]
figure(6)
loglog(time(a:end),msd(a:end,p),'.-','color',[rand,rand,rand],'MarkerSize',25)
xlabel('t (ps)','Interpreter','latex')
ylabel('MSD ($\AA^{2}$)','Interpreter','latex')
hold on

figure(7)
% plot(ltime,lmsd(:,p),'-^')
plot(ltime(a:end),lmsd(a:end,p),'.-','color',[rand,rand,rand],'MarkerSize',25)

xlabel('log(t)(ps) ','Interpreter','latex')
ylabel('log(MSD)($\AA^{2}$) ','Interpreter','latex')
hold on
b=polyfit(ltime(a:end),lmsd(a:end,p)',1);
Blog(p,:)=b;
end
legend('T=5 K','T=10 K','T=20 K','T=30 K','T=35 K','T=50 K',...
       'T=60 K','T=75 K','T=100 K','T=150 K','T=200 K','T=250 K','T=300 K',...
       'T=400 K','T=500 K','T=600 K','T=700 K','T=800 K')
set(gca,'FontName','Cambria','FontSize',14);
figure(6)
legend('T=5 K','T=10 K','T=20 K','T=30 K','T=35 K','T=50 K',...
       'T=60 K','T=75 K','T=100 K','T=150 K','T=200 K','T=250 K','T=300 K',...
       'T=400 K','T=500 K','T=600 K','T=700 K','T=800 K')
set(gca,'FontName','Cambria','FontSize',12);
% D=(10.^(Blog(:,2)))/4;
figure(8)
for p=[2 16]
loglog(time,msd(:,p),'.-','color',[rand,rand,rand],'MarkerSize',25)
xlabel('t (ps)','Interpreter','latex')
ylabel('MSD ($\AA^{2}$)','Interpreter','latex')
hold on
end
legend('T=30 K','T=300 K')
set(gca,'FontName','Cambria','FontSize',16);

line([1.5 3],[3 6],'Color','black','LineStyle','--')
line([1.5 3],[0.5 2],'Color','black','LineStyle','--')

%% calculation of alfa power

alfa=smooth(Blog(:,1),'sgolay');
figure(9)
hold on
plot(T(1:16), alfa,'-','LineWidth',3,'MarkerSize',12)
xlabel('Temperature (K)','Interpreter','latex')
ylabel('$$\alpha$$','Interpreter','latex')
set(gca,'FontName','Cambria','FontSize',16);
line([30 30],[0 2],'Color','black','LineStyle','--')
line([150 150],[0 2],'Color','black','LineStyle','--')

% axis([5 500 0.5 2])
% line([150 150],[0 2],'Color','black','LineStyle','--')
% line([35 35],[0 2],'Color','black','LineStyle','--')
% legend('C60/mono-Gr','C60/2layer-Gr','C60/fixed-Gr','Location','Best')

figure(30)
hold on
errorbar(TT,Fmean,Fstd,'linewidth',3,'MarkerSize',10)
% errorbar(TT,Dmean,Dstd,'linewidth',3,'MarkerSize',10)
% plot(TT,Dmean,'-bs','linewidth',3,'MarkerSize',12)
box on
grid off
xlabel('Temperature (K)','fontsize',20,'FontName','Cambria')
ylabel('D (Å^{2}/ps)','fontsize',20,'FontName','Cambria')
set(gca,'FontName','Cambria','FontSize',20);

figure(31)
T_1= 1./TT(2:end);
lnD= log(Fmean(2:end));
lnDstd= log(Fstd(2:end));
plot(T_1, lnD,'.-','LineWidth',4,'MarkerSize',16)
% errorbar(T_1,lnD,lnDstd,'linewidth',3,'MarkerSize',10)
hold on
% yy = smooth(lnD,'sgolay');
% plot(T_1, yy,'.-','LineWidth',4,'MarkerSize',16)
xlabel('$ 1/T (K^{-1}) $','Interpreter','latex')
ylabel('lnD ($\AA^{2}/ps$)','Interpreter','latex')
set(gca,'FontName','Cambria','FontSize',20);
% legend('C60/SLG','C60/DLG','C60/FLG','Location','Best')
at=10;
b1=polyfit(T_1(1:3),lnD(1:3),1);
b2=polyfit(T_1(at:end),lnD(at:end),1);
yhat1=polyval(b1,T_1(1:7));
yhat2=polyval(b2,T_1(at-3:end));
plot(T_1(1:7),yhat1,'k--','MarkerSize',55)
plot(T_1(at-3:end),yhat2,'k--','MarkerSize',55)

% msdp=msd(start:end,p)'; %700
% b=polyfit(time(start:end),msdp,1);
% xlim=time(:);
% yhat=polyval(b,xlim);           % evaluate at origin and points
% figure(2)
% plot(xlim,yhat,'k--','MarkerSize',55)
% B(p,2:3)=b;

