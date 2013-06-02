function findfibre3(Zones, zones0, NumbFiles,Path,PoscDir)
% Usage : findfibre(Zones, zones0, NumbFiles)
%
% Description: Given a set of files that possess the trajectories of
% calculated paths, the begining and arrival zones of these and the
% probability of each path, this function unifies all the paths that unite
% each couple of zones of interest, as well as it unifies the probabilities
% corresponding to these paths.
%
% Input Parameters:
%       Zones       : struct file of fiels: 'name', 'voxels', 'number' and
%                     'roi' corresponding to gray matter defined zones.
%       zones0      : number of the zones of interest.
%       NumbFiles   : number of files kept with information on the paths.
%
% Output Parameters: This function saves two different types of files. For
% example: The file Pbelong1_2 have the position of the original files that
% contains the paths that unite zones 1 and 2, qbelong1_2 is a variable struct
% that in the field 'voxels' possesses the coordinates of each path that unite
% the two zones and also, for the probabilistic methods, in the field 'PT'
% possesses the probability or validity index of each of these paths.
%--------------------------------------------------------------------------
% Authors: Yasser Iturria Medina & Pedro Vald�s-Hern�ndez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 15th 2005
% Version $1.0

cwd = cd;
if ~exist([Path filesep 'Pbelong'],'dir')
    mkdir([Path filesep 'Pbelong'])
end
if ~exist([Path filesep 'qbelong'],'dir')
    mkdir([Path filesep 'qbelong'])
end
%Path = uigetdir('','Dir to load files Poscs_ and Paths_, and to save results ...');
cd(Path)
if nargin < 3, NumbFiles = size(spm_get('files',Path,'Paths_*.mat'),1); end
counter = 0;
zones0 = sort(unique(zones0));
for j = 1:size(zones0,2) - 1
    for i = j+1:size(zones0,2)
        counter = counter + 1;
        zones(counter,:) = [zones0(j) zones0(i)];
        eval(strcat('Pbelong',num2str(zones0(j)),'_',num2str(zones0(i)),'= [];'));
    end
end
for j = 1:size(zones0,2)
    numbers(j) = Zones(zones0(j)).number;
end
for file = 1:NumbFiles
    disp([num2str(file) ' -> ' num2str(NumbFiles)]);
    eval(['load ' PoscDir filesep 'Poscs_' num2str(file)]);
    if ~isempty(posc)
        for i = 1:size(posc,1)
            Z1 = find(numbers == posc(i,1));
            Z2 = find(numbers == posc(i,2));
            if Z1 > Z2, temp = Z2; Z2 = Z1; Z1 = temp; end
            if (~isempty(Z1)) && (~isempty(Z2))
                eval(strcat('Pbelong',num2str(zones0(Z1)),'_',num2str(zones0(Z2)),...
                    ' = [Pbelong',num2str(zones0(Z1)),'_',num2str(zones0(Z2)),'; file i];'));
            end
        end
    end
end
for j = 1:size(zones,1)
    eval(strcat('save',' Pbelong',filesep,'Pbelong',num2str(zones(j,1)),'_',num2str(zones(j,2)),...
        ' Pbelong', num2str(zones(j,1)),'_',num2str(zones(j,2))));
    eval(strcat('clear Pbelong',num2str(zones(j,1)),'_',num2str(zones(j,2))));
end

fileAnt = 0;
for i = 1:size(zones,1)
    disp([num2str(i) ' -> ' num2str(size(zones,1))]);
    eval(strcat('load',' Pbelong',filesep,'Pbelong',num2str(zones(i,1)),'_',num2str(zones(i,2))));
    eval(strcat('Pbelong = Pbelong',num2str(zones(i,1)),'_',num2str(zones(i,2)),';'));
    qbelong = struct('voxels',[],'PT',[]);
    for j = 1:size(Pbelong,1)
        if Pbelong(j,1) ~= fileAnt
            eval(strcat('load',' Paths',filesep,'Paths_',num2str(Pbelong(j,1))));
            fileAnt = Pbelong(j,1);
        end
        eval(strcat('qbelong(',num2str(j),') = q(Pbelong(',num2str(j),',2));'));
    end
    eval(strcat('save',' qbelong',filesep,'qbelong',num2str(zones(i,1)),'_',num2str(zones(i,2)),...
        ' qbelong'));
    eval(strcat('clear qbelong',num2str(zones(i,1)),'_',num2str(zones(i,2)),' qbelong'));
end
cd(cwd)