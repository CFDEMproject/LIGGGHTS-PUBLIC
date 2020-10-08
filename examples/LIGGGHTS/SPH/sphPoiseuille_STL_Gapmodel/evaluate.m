% (c) by Andreas Eitzlmayr

clear all

folder='post/';
ResultFolder='';
fileName='dump';
fileExt='.sph';

NoHeadlines = 9;

% case parameters
rhoDef = 1000; % kg/m³
MuDef = 0.5; % Pas
FDef = 10; % gradP / rho [m/s2]
bhalf = 0.0001; % half gap width [mm]
bcenter = 0;
Nu = MuDef/rhoDef;

% time
t0=0;
step=1; % saved steps
tmax=100; % steps total
dt = 1e-6; % timestep [s]

t = [100 10 2]; % plot times [steps]
ts = t * dt; % plot times [s]
Nt = length(t);

N = (tmax-t0)/step + 1;

for i = 1:Nt
    %Idx = t0 + (i-1)*step;
    Idx = t(i);
    dump=importdata([folder,fileName,num2str(Idx),fileExt],' ', NoHeadlines);

    % ITEM: ATOMS id type x y z ix iy iz vx vy vz fx fy fz p rho f_int
    % f_wallCount f_usedGapmodel f_dvdx[1] f_dvdx[2] f_dvdx[3] f_dvdy[1] f_dvdy[2] f_dvdy[3] f_dvdz[1] f_dvdz[2] f_dvdz[3] f_gamma f_omega f_mixidx

    x(:,i) = dump.data(:,3);
    vy(:,i) = -dump.data(:,10);
end

%% Calculate velocity profile for Couette-Poiseuille flow

H = 2*bhalf;

Nx = 201;
xCalc = zeros(Nt,Nx);
VCalc = xCalc;

for it = 1:Nt
    xCalc(it,:) = linspace(0,H,Nx);
    
    for i = 0:100
        VCalc(it,:) = VCalc(it,:) + (2*i+1)^-3 * sin(pi*xCalc(it,:)/H*(2*i+1)) * exp(-(2*i+1)^2*pi^2*Nu*ts(it)/H^2) ;
    end
    
    VCalc(it,:) = - VCalc(it,:) * 4*FDef*H^2/(Nu * pi^3);
    VCalc(it,:) = VCalc(it,:) - FDef/(2*Nu) * xCalc(it,:).*(xCalc(it,:)-ones(1,Nx)*H);
end

xCalc = -xCalc + ones(Nt,Nx) * (bhalf+bcenter);



%% Evolution of v average

for i = 1:(tmax-t0)/step+1
    Idx = t0 + (i-1)*step;
    dump=importdata([folder,fileName,num2str(Idx),fileExt],' ', NoHeadlines);
    vyav(i) = mean(-dump.data(:,10));
    tTotal(i) = Idx;
end
integrity = mean(dump.data(:,17));
wallCount = mean(dump.data(:,18));
usedGapm = mean(dump.data(:,19));
dvdx = mean(-dump.data(:,21));
mixIdx = mean(dump.data(:,31));

disp(['vy Average = ',num2str(vyav(end)),'m/s']);
disp(['integrity = ',num2str(integrity)]);
disp(['wall count = ',num2str(wallCount)]);
disp(['used gapmodel = ',num2str(usedGapm)]);
disp(['dv/dx = ',num2str(dvdx)]);
disp(['mixing index = ',num2str(mixIdx)]);

plot(tTotal,vyav,'-b','LineWidth',2);

set(gca,'FontSize',16);
ylabel('vy average [m/s]');
xlabel('timestep');

% Save graphic:
fileSaveName=[ResultFolder,'figEvolution'];
set(gcf,'PaperpositionMode','auto');
print(gcf,'-depsc','-r200',[fileSaveName,'.eps']);
%print(gcf,'-dtiff','-r200',[fileSaveName,'.tiff']);
disp(['Saved as: ', fileSaveName]);



%% Plot result
lineformat = char('ko','kd','ks','k+','kx','k*');

vmax = 0.00015;
AvVel1 = 2/3*max(VCalc(1,:));
AvVel2 = mean(VCalc(2,:));
AvVel3 = mean(VCalc(3,:));

x = x*1000;
xCalc = xCalc*1000;
ts = ts*1000;

plot(x(:,1),vy(:,1),lineformat(1,:),'LineWidth',2);
hold on
for i=2:Nt
    plot(x(:,i),vy(:,i),lineformat(i,:),'LineWidth',2);
end
plot(xCalc(1,:),VCalc(1,:),'-r','LineWidth',1);
plot([(bcenter-bhalf)*1000 (bhalf+bcenter)*1000],[AvVel1 AvVel1],'--k')
plot([(bcenter-bhalf)*1000 (bhalf+bcenter)*1000],[AvVel2 AvVel2],'--k')
plot([(bcenter-bhalf)*1000 (bhalf+bcenter)*1000],[AvVel3 AvVel3],'--k')
for i=2:Nt
    plot(xCalc(i,:),VCalc(i,:),'-r','LineWidth',1);
end
hold off

set(gca,'FontSize',16);
ylabel('vy [m/s]');
xlabel('x [mm]');

legend(['t = ',num2str(ts(1)),' ms'],['t = ',num2str(ts(2)),' ms'], ...
    ['t = ',num2str(ts(3)),' ms'],'analytical','average','Location','EastOutside')

axis([(bcenter-bhalf)*1000 (bhalf+bcenter)*1000 0 vmax])

% Save graphic:
fileSaveName=[ResultFolder,'figVelProfile'];
set(gcf,'PaperpositionMode','auto');
print(gcf,'-depsc','-r200',[fileSaveName,'.eps']);
%print(gcf,'-dtiff','-r200',[fileSaveName,'.tiff']);
disp(['Saved as: ', fileSaveName]);
close