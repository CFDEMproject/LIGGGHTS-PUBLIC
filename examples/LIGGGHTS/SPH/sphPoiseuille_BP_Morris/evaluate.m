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
b = 0.005; % gap width [mm]
Nu = MuDef/rhoDef;

% time
t0=0;
step=200; % saved steps
tmax=30000; % steps total
dt = 1e-6; % timestep [s]

t = [30000 6000 1200]; % plot times [steps]
ts = t * dt; % plot times [s]
Nt = length(t);

N = (tmax-t0)/step + 1; % number of saved steps


for i = 1:Nt
    %Idx = t0 + (i-1)*step;
    Idx = t(i);
    dump=importdata([folder,fileName,num2str(Idx),fileExt],' ', NoHeadlines);

    % ITEM: ATOMS id type x y z ix iy iz vx vy vz fx fy fz p rho f_int
    % f_dvdx[1] f_dvdx[2] f_dvdx[3] f_dvdy[1] f_dvdy[2] f_dvdy[3] f_dvdz[1] f_dvdz[2] f_dvdz[3] f_gamma f_omega f_mixidx
    
    x(:,i) = dump.data(:,3);
    vy(:,i) = -dump.data(:,10);
    dvdx(:,i) = -dump.data(:,19);
    mixIdx(:,i) = dump.data(:,29);
end


%% Calculate velocity profile for Poiseuille flow

Nx = 41;
xCalc = zeros(Nt,Nx);
VCalc = xCalc;
dvdxCalc = xCalc(:,1:end-1);

for it = 1:Nt
    xCalc(it,:) = linspace(0,b,Nx);
    
    for i = 0:100
        VCalc(it,:) = VCalc(it,:) + (2*i+1)^-3 * sin(pi*xCalc(it,:)/b*(2*i+1)) * exp(-(2*i+1)^2*pi^2*Nu*ts(it)/b^2) ;
    end
    
    VCalc(it,:) = - VCalc(it,:) * 4*FDef*b^2/(Nu * pi^3);
    VCalc(it,:) = VCalc(it,:) - FDef/(2*Nu) * xCalc(it,:).*(xCalc(it,:)-ones(1,Nx)*b);
    
    dvdxCalc(it,:) = -(VCalc(it,2:end) - VCalc(it,1:end-1))./(xCalc(it,2:end)-xCalc(it,1:end-1));
end

xCalc = -xCalc + ones(Nt,Nx) * b;
xCalc2 = (xCalc(:,2:end)+xCalc(:,1:end-1))/2;

%% Plot result
lineformat = char('ko','kd','ks','k+','kx','k*');

vmax = 0.1;

x = x*1000;
xCalc = xCalc*1000;
xCalc2 = xCalc2*1000;
ts = ts*1000;

plot(x(:,1),vy(:,1),lineformat(1,:),'LineWidth',2);
hold on
for i=2:Nt
    plot(x(:,i),vy(:,i),lineformat(i,:),'LineWidth',2);
end
for i=1:Nt
    plot(xCalc(i,:),VCalc(i,:),'-r','LineWidth',1);
end
hold off

set(gca,'FontSize',16);
ylabel('vy [m/s]');
xlabel('x [mm]');

legend(['SPH/BP, t = ',num2str(ts(1)),' ms'],['SPH/BP, t = ',num2str(ts(2)),' ms'], ...
    ['SPH/BP, t = ',num2str(ts(3)),' ms'],'analytical','Location','NorthEast')

axis([-0.5 5.5 0 vmax])

% Save graphic:
fileSaveName=[ResultFolder,'figVelProfile'];
set(gcf,'PaperpositionMode','auto');
print(gcf,'-depsc','-r200',[fileSaveName,'.eps']);
%print(gcf,'-dtiff','-r200',[fileSaveName,'.tiff']);
disp(['Saved as: ', fileSaveName]);


plot(x(:,1),dvdx(:,1),lineformat(1,:),'LineWidth',2);
hold on
for i=2:Nt
    plot(x(:,i),dvdx(:,i),lineformat(i,:),'LineWidth',2);
end
for i=1:Nt
    plot(xCalc2(i,:),dvdxCalc(i,:),'-r','LineWidth',1);
end
hold off

set(gca,'FontSize',16);
ylabel('dv/dx [1/s]');
xlabel('x [mm]');

legend(['SPH/BP, t = ',num2str(ts(1)),' ms'],['SPH/BP, t = ',num2str(ts(2)),' ms'], ...
    ['SPH/BP, t = ',num2str(ts(3)),' ms'],'analytical','Location','NorthEast')

axis([-0.5 5.5 -60 60])

% Save graphic:
fileSaveName=[ResultFolder,'figVelGradient'];
set(gcf,'PaperpositionMode','auto');
print(gcf,'-depsc','-r200',[fileSaveName,'.eps']);
%print(gcf,'-dtiff','-r200',[fileSaveName,'.tiff']);
disp(['Saved as: ', fileSaveName]);
close;


plot(x(:,1),mixIdx(:,1),lineformat(1,:),'LineWidth',2);
hold on
for i=2:Nt
    plot(x(:,i),mixIdx(:,i),lineformat(i,:),'LineWidth',2);
end
plot([0 5],[0.5 0.5 ],'-r','LineWidth',1);
hold off

set(gca,'FontSize',16);
ylabel('mixing index');
xlabel('x [mm]');

legend(['SPH/BP, t = ',num2str(ts(1)),' ms'],['SPH/BP, t = ',num2str(ts(2)),' ms'], ...
    ['SPH/BP, t = ',num2str(ts(3)),' ms'],'analytical','Location','SouthEast')

axis([-0.5 5.5 0 1])

% Save graphic:
fileSaveName=[ResultFolder,'figMixIdx'];
set(gcf,'PaperpositionMode','auto');
print(gcf,'-depsc','-r200',[fileSaveName,'.eps']);
%print(gcf,'-dtiff','-r200',[fileSaveName,'.tiff']);
disp(['Saved as: ', fileSaveName]);


%% Evolution of v average

for i = 1:(tmax-t0)/step+1
    Idx = t0 + (i-1)*step;
    dump=importdata([folder,fileName,num2str(Idx),fileExt],' ', NoHeadlines);
    vyav(i) = mean(-dump.data(:,10));
    tTotal(i) = Idx;
end

disp(['vy Average = ',num2str(vyav(end)),'m/s']);

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
