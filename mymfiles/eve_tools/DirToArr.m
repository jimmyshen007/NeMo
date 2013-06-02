function [arr_name,sz]=DirToArr(varargin)
%DEPRECATED: DirToArr sends dir search output to matlab array variable

MSL=500;

read_dir=varargin{1};

arr_name=[];

if nargin<2
    g=dir(read_dir);
    S=size(g);
    for i=1:S(1)
        C=g(i).name;
        if ~strcmp(g(i).name(1),'.')
            arr_name=[arr_name; [read_dir filesep C blanks(MSL-(length(C)+length(read_dir)+1))]];
        end
    end    
    sz=S(1);
    return;
else
    for v=2:nargin
        eval(['modif' int2str(v-1) '=varargin{v}; g=dir([read_dir filesep modif' int2str(v-1) ']);']);
        S=size(g);
        for i=1:S(1)
            C=g(i).name;
            arr_name=[arr_name; [read_dir filesep C blanks(MSL-(length(C)+length(read_dir)+1))]];
        end    
    end
    S=size(arr_name); sz=S(1);
end

% S=size(g);
% arr_name=[];

disp(['DirToArr: ' int2str(sz) ' entries.']);

return;

end

