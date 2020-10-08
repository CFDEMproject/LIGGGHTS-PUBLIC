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
H = 5; % [mm] height

% time
t0=0;
step=100; % saved steps
tmax=2000; % steps total
dt = 1e-6; % timestep [s]

t = [2000]; % plot times [steps]
ts = t * dt; % plot times [s]
Nt = length(t);

N = (tmax-t0)/step + 1;

for i = 1:Nt
    %Idx = t0 + (i-1)*step;
    Idx = t(i);
    dump=importdata([folder,fileName,num2str(Idx),fileExt],' ', NoHeadlines);

    % ITEM: ATOMS id type x y z ix iy iz vx vy vz fx fy fz p rho f_int
    % f_fgradP[1] f_fgradP[2] f_fgradP[3] f_wallCount f_usedGapmodel ...
    % f_dvdx[1] f_dvdx[2] f_dvdx[3] f_dvdy[1] f_dvdy[2] f_dvdy[3] ...
    % f_dvdz[1] f_dvdz[2] f_dvdz[3] f_gamma f_omega f_mixidx 

    y(:,i) = dump.data(:,4);
    p(:,i) = dump.data(:,15);
    integrity(:,i) = dump.data(:,17);
    fgradPy(:,i) = dump.data(:,19);
end

%% Calculate velocity profile for Couette-Poiseuille flow

% H = 2*bhalf;
% 
% Nx = 201;
% xCalc = zeros(Nt,Nx);
% VCalc = xCalc;
% 
% for it = 1:Nt
%     xCalc(it,:) = linspace(0,H,Nx);
%     
%     for i = 0:100
%         VCalc(it,:) = VCalc(it,:) + (2*i+1)^-3 * sin(pi*xCalc(it,:)/H*(2*i+1)) * exp(-(2*i+1)^2*pi^2*Nu*ts(it)/H^2) ;
%     end
%     
%     VCalc(it,:) = - VCalc(it,:) * 4*FDef*H^2/(Nu * pi^3);
%     VCalc(it,:) = VCalc(it,:) - FDef/(2*Nu) * xCalc(it,:).*(xCalc(it,:)-ones(1,Nx)*H);
% end
% 
% xCalc = -xCalc + ones(Nt,Nx) * (bhalf+bcenter);



%% Evolution of p max

for i = 1:(tmax-t0)/step+1
    Idx = t0 + (i-1)*step;
    dump=importdata([folder,fileName,num2str(Idx),fileExt],' ', NoHeadlines);
    pmax(i) = max(dump.data(:,15));
    tTotal(i) = Idx;
end

plot(tTotal,pmax,'-b','LineWidth',2);

set(gca,'FontSize',16);
ylabel('p max [Pa]');
xlabel('timestep');

% Save graphic:
fileSaveName=[ResultFolder,'figEvolution'];
set(gcf,'PaperpositionMode','auto');
print(gcf,'-depsc','-r200',[fileSaveName,'.eps']);
%print(gcf,'-dtiff','-r200',[fileSaveName,'.tiff']);
disp(['Saved as: ', fileSaveName]);



%% Plot result
lineformat = char('ko','kd','ks','k+','kx','k*');


y = y*1000;
ts = ts*1000;

poly = polyfit(y,p,1);
gradP = poly(1); % [Pa/mm]

plot(y(:,1),p(:,1),lineformat(1,:),'LineWidth',2);
hold on

plot([0 H],[poly(2) poly(2)+gradP*H],'-r','LineWidth',1);

% for i=2:Nt
%     plot(x(:,i),vy(:,i),lineformat(i,:),'LineWidth',2);
% end
% plot(xCalc(1,:),VCalc(1,:),'-r','LineWidth',1);
% plot([(bcenter-bhalf)*1000 (bhalf+bcenter)*1000],[AvVel1 AvVel1],'--k')
% plot([(bcenter-bhalf)*1000 (bhalf+bcenter)*1000],[AvVel2 AvVel2],'--k')
% plot([(bcenter-bhalf)*1000 (bhalf+bcenter)*1000],[AvVel3 AvVel3],'--k')
% for i=2:Nt
%     plot(xCalc(i,:),VCalc(i,:),'-r','LineWidth',1);
% end
hold off

set(gca,'FontSize',16);
xlabel('y [mm]');
ylabel('p [Pa]');

legend('SPH-gapmodel',['fit: ',num2str(gradP*1000),' Pa/m'],'Location','NorthEast')

axis([0 H 0 50]);


% Save graphic:
fileSaveName=[ResultFolder,'figPressProfile'];
set(gcf,'PaperpositionMode','auto');
print(gcf,'-depsc','-r200',[fileSaveName,'.eps']);
%print(gcf,'-dtiff','-r200',[fileSaveName,'.tiff']);
disp(['Saved as: ', fileSaveName]);

disp(['gradP = ',num2str(gradP*1000),'Pa/m']);
disp(['integrity av. = ',num2str(mean(integrity))]);

close