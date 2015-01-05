clear
clc
close all
more off
% ---------------------------------------------------------------------
%   USER INPUT
k=1;
data2PlotList{k} = 'heat.json';             k=k+1;

n=1;
particleIDs2Plot(n) = 135;   n=n+1;
particleIDs2Plot(n) = 677; n=n+1;

j=1;
gridPoint2Plot{j} = 1; symbol{j}='ro'; faceC{j}='r'; j=j+1;     %Center
gridPoint2Plot{j} = 5; symbol{j}='k<'; faceC{j}='k'; j=j+1;     %Middle
gridPoint2Plot{j} = 99; symbol{j}='bv'; faceC{j}='none';j=j+1;  %Surface

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
particles2Plot = size(particleIDs2Plot,2);

for iDat=1:size(data2PlotList,2)
data2Plot=data2PlotList{iDat}
%Scan dir for data
allDirs = dir('.'); isub = [allDirs(:).isdir]; allDirs=allDirs(isub);
allDirs(1)='';allDirs(1)=''; %cut-off current and parent dit
myFile = dir(['*/',data2Plot]);

for particleID=1:particles2Plot
disp(['Processing particle: ', num2str(particleIDs2Plot(particleID))]);
figure
raw = struct;
validTime = 0;
for iFile=1:size(allDirs,1)
    myFile = dir([allDirs(iFile).name,'/',data2Plot]);
    if(size(myFile,1)>0)
        currentFile = [allDirs(iFile).name,'/',myFile(1).name];
        disp(['processing file:', currentFile]),
        dat = loadjson(currentFile);
        cell_data = struct2cell(getfield(dat,'data'));
        number_particles = size(cell_data,1);
        i=particleIDs2Plot(particleID); %load current particle!
        particle_data = str2double(cell_data{i});
        if(isnan(particle_data)(1)==0)
            validTime = validTime + 1;
            raw.data{validTime}       = particle_data;
            raw.time(validTime)       = str2num(allDirs(iFile).name);
          
            %Extract data to plot
            lastPoint = size(raw.data{validTime},2);
            raw.fluid(validTime)      = raw.data{validTime}(lastPoint);
            for j=1:size(gridPoint2Plot,2)
                if(gridPoint2Plot{j}>=lastPoint)
                    gridPoint2Plot{j} = lastPoint-1;
                end
                raw.plot{j}.dat(validTime) = raw.data{validTime}(gridPoint2Plot{j});
            end
            raw.gridpoints(validTime) = size(particle_data,2);
        end
    end
end

plot(raw.time, raw.fluid,'kd-');
hold on
for j=1:size(gridPoint2Plot,2)
    plot(raw.time, raw.plot{j}.dat,symbol{j},'MarkerFaceColor', faceC{j});
end
  
%axis([1 28 600 650])  
leg=legend('fluid','particle_center','particle_middle','particle_surf');
legend boxoff
set(leg,'Location','NorthOutside');
xlabel('time (s)','Fontsize',11)
ylabel(data2Plot,'Fontsize',11)
set(gca,'Fontsize',10);
xlhand = get(gca,'xlabel');ylhand = get(gca,'ylabel');
set(ylhand,'Position',get(ylhand,'Position') - [1. 0 0])
set(gcf, 'paperunits', 'centimeters', 'paperposition', [0 0 21 18])

print('-depsc2',['plot_',data2Plot,'_',num2str(particleIDs2Plot(particleID)),'.eps'])

end
end

close all
