function returnarr=system_find(searchdir,ext,maxdepth,tofilenm)
% returns <Nx500> char array with results from system call using 'find'
% command 

MSL=500;

if nargin < 3 || isempty(maxdepth)
	mxdpthexp='';
else
    mxdpthexp=['-maxdepth ' int2str(maxdepth)];
end
    
if size(searchdir,1) > 1
    tmp=[];
    for i=1:size(searchdir,1)
        tmp=[tmp ' ' deblank(searchdir(i,:))];
    end
    searchdir=tmp;
    clearvars tmp;
else
    searchdir=deblank(searchdir);
end

if nargin < 4
    [~,res]=system(['find ' searchdir ' ' mxdpthexp ' -name ''' ext '''']);
    if isempty(res)
        returnarr=[];
    elseif ~strcmp(res(1:5),'find:')
        C=textscan(res,'%[^\n]');
        returnarr=[];
        for i=1:size(C{1},1)
            tmp=C{1}{i};
            returnarr=[returnarr; tmp blanks(MSL-length(tmp))];
        end
    else
        returnarr=[];
    end
else
    system(['find ' searchdir ' ' mxdpthexp ' -name ''' ext ''' > ' tofilenm]);
    returnarr=ReadInTxt(tofilenm);
end

return;

end