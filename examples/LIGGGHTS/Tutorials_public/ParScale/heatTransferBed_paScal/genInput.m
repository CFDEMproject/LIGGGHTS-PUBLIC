clear
clc
more off

text = '[400,400,400,400,400,400,400,400,400,400,  350]'
%text = '[1]'

numberOfParticles = 800

fid=fopen('output.dat','w')

for iP=1:numberOfParticles
    index=iP-1;
    fprintf(fid, ['        "',num2str(index), '":  ', text ]);
    if(iP==numberOfParticles)
        fprintf(fid, ['\n']);
    else
        fprintf(fid, [', \n']);
    end        
end
    
fclose(fid)
