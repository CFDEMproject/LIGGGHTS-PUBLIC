% *.m file
clear
clc
more off

numberOfParticles = 100

%================================================================
text = '[293,293,293,293,293,293,293,293,293,293,  323]';

fid=fopen('0/heat.json','w');

% header
fprintf(fid, ['{\n']);
fprintf(fid, ['    "name": "heat",\n']);
fprintf(fid, ['    "data":\n']);
fprintf(fid, ['    {\n']);

for iP=1:numberOfParticles
    index=iP;
    fprintf(fid, ['        "',num2str(index), '":  ', text ]);
    if(iP==numberOfParticles)
        fprintf(fid, ['\n']);
    else
        fprintf(fid, [', \n']);
    end        
end

% footer
fprintf(fid, ['    }\n']);
fprintf(fid, ['}']);
fclose(fid);
%================================================================
text = '[0.001]';

fid=fopen('0/radius.json','w');

% header
fprintf(fid, ['{\n']);
fprintf(fid, ['    "name": "radius",\n']);
fprintf(fid, ['    "data":\n']);
fprintf(fid, ['    {\n']);

for iP=1:numberOfParticles
    index=iP;
    fprintf(fid, ['        "',num2str(index), '":  ', text ]);
    if(iP==numberOfParticles)
        fprintf(fid, ['\n']);
    else
        fprintf(fid, [', \n']);
    end        
end

% footer
fprintf(fid, ['    }\n']);
fprintf(fid, ['}']);
fclose(fid);
%================================================================
