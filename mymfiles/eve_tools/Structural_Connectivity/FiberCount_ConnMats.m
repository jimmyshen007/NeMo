function FiberCount_ConnMats(PbelongDir,atlassize,OutPath)

pbdir = Dir2Arr(PbelongDir,'Pb*.mat');
MapTFC = zeros(atlassize);
total_fcount = 0;
for j = 1:size(pbdir,1);
    if mod(j,500); 
        disp(['Calculating pair-wise connectivity regions ' num2str(j) ' of ' num2str(size(pbdir,1))])
    end
    str = deblank(pbdir(j,:));
    load(str)
    us = strfind(str,'_');
    dots = strfind(str,'.');
    gs = strfind(str,'g');
    r1 = str2num(str(gs(end)+1:us(end)-1));
    r2 = str2num(str(us(end)+1:dots-1));
    eval(['fcount = size(Pbelong' num2str(r1) '_' num2str(r2) ',1);'])
    MapTFC(r1,r2) = fcount;
    total_fcount = total_fcount + fcount;
end
MapTFC = MapTFC + MapTFC';
nMapTFC = MapTFC./total_fcount;
save([OutPath filesep 'Matrix_TFC' num2str(atlassize) '_' OutPath(strfind(OutPath,'e0'):strfind(OutPath,'e0')+6)],'MapTFC','nMapTFC','total_fcount')
    