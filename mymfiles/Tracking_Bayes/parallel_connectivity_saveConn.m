function parallel_connectivity_saveConn(zones0,Zones,Path)
% Usage : [Connectivity_Map] = connectivity(zones0,Zones)
%
% Description: This function is to quantify anatomical connection strength
%              (ACS), anatomical connection density (ACD) and anatomical 
%              connection probability (ACP) between n gray matter zones using
%              information of computed routes that simulate nervous fibers 
%              trajectories.
%
% Input Parameters:
%        zones0   : number of the zones among which will be considered the
%                   conectivity.
%        Zones    : struct file of fiels: 'name', 'voxels', 'number' and
%                   'roi' corresponding to gray matter defined zones.
%
% Output Parameters:
%  Connectivity_Map: two-dimensional matriz that possesses the conectivity
%  values calculated for all the pairs of zones. Example: the element i,j
%  is the connectivity value between zones i and j, this value will be saved
%  as Ci_j in the chosen directory.
%--------------------------------------------------------------------------
% Useful References:
% 1.- Y. Iturria-Medina, P. Valdes-Hernandez, E. Canales-Rodriguez. 2005.
% Measures of anatomical connectivity obtained from neuro diffusion images.
% Presented at the 11th Annual Meeting of the Organization for Human Brain
% Mapping, June 12-16, 2005, Toronto, Ontario, Canada. Available on CD-Rom
% in NeuroImage, Vol. 26, No.1.
%--------------------------------------------------------------------------
% Author: Yasser Iturria Medina
% Neuroimaging Department
% Cuban Neuroscience Center
% November 15th 2005
% Version $1.0

zones0 = sort(zones0);
counter = 0;
for j = 1:size(zones0,2) - 1
    for i = j+1:size(zones0,2)
        counter = counter + 1;
        zones(counter,1) = zones0(j);
        zones(counter,2) = zones0(i);
    end
end
%cwd = cd;
% if nargin>2 && ~isempty(Path)
%     Path = uigetdir('','Dir to load files qbelong, and to save Connectivity_Map ...');
% end
%cd(Path)
conn1 = NaN(size(zones,1),1);
conn2 = NaN(size(zones,1),1);
conn3 = NaN(size(zones,1),1);
OutConPath = [Path filesep 'ConnInfo'];
if ~isdir(OutConPath)
    mkdir(OutConPath)
end
parfor i = 1:size(zones,1)
    disp([num2str(i) ' -> ' num2str(size(zones,1))]);
    Q = load([Path filesep 'qbelong' filesep 'qbelong' num2str(zones(i,1)) '_' num2str(zones(i,2))]);
    ROI1 = Zones(zones(i,1)).roi;
    ROI2 = Zones(zones(i,2)).roi;
    [ACS, ACD, ACP] = Areas_saveConn(ROI1,ROI2,Q.qbelong,zones(i,1),zones(i,2),OutConPath);
    Q = [];
    %     eval(strcat('load Conn',num2str(zones(i,1)),'_',num2str(zones(i,2)),' ACS ACD ACP'));
    conn1(i) = ACS; conn2(i) = ACD; conn3(i) = ACP;
    %save(['Conn' num2str(zones(i,1)) '_' num2str(zones(i,2))],'ACS','ACD','ACP');
end
for i = 1:size(zones,1)
    ACS = conn1(i);
    ACD = conn2(i);
    ACP = conn3(i);    
   % save(['Conn' num2str(zones(i,1)) '_' num2str(zones(i,2))],'ACS','ACD','ACP');
end
counter = 0;
MapACS = zeros(size(zones0,2),size(zones0,2));
MapACD = zeros(size(zones0,2),size(zones0,2));
MapACP = zeros(size(zones0,2),size(zones0,2));
for i = 1:size(zones0,2) - 1
    for j = i+1 : size(zones0,2)
        counter = counter + 1;
        MapACS(i,j) = conn1(counter);
        MapACD(i,j) = conn2(counter);
        MapACP(i,j) = conn3(counter);
    end
end
MapACS = MapACS + MapACS';

Connectivity_Map = {};

figure;
% subplot(131);
imagesc(MapACS); title('ACS'); hold on
plot(1:size(zones,2),1:size(zones,2),'x','Markersize',10,'color','w');
colorbar; Map = MapACS; save([Path filesep 'Matrix_ACS'], 'Map');
Connectivity_Map{1} = [Path filesep 'Matrix_ACS'];

MapACD = MapACD + MapACD';
figure;
% subplot(132);
imagesc(MapACD); title('ACD'); hold on
plot(1:size(zones,2),1:size(zones,2),'x','Markersize',10,'color','w');
colorbar; Map = MapACD;save([Path filesep 'Matrix_ACD'], 'Map');
Connectivity_Map{2} = [Path filesep 'Matrix_ACD'];

MapACP = MapACP + MapACP';
figure;
% subplot(133);
imagesc(MapACP); title('ACP'); hold on
plot(1:size(zones,2),1:size(zones,2),'x','Markersize',10,'color','w');
colorbar; Map = MapACP; save([Path filesep 'Matrix_ACP'], 'Map');
Connectivity_Map{3} = [Path filesep 'Matrix_ACP'];