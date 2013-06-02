function Organize_PathsByRegion(ROI_numFileName,qPpaths,InDir,OutDir)
% This function assigns each voxel to a larger region of interest to perform
% the core search of the white matter tracts. 
%
% Input:   PatientString  [String]  This is the unique Patient ID
%          CleanWMFileName [String]  File name of the white matter mask that
%                                   has been cleaned of the isolated
%                                   pixels.
%          szbox          [Integer] The size (in pixels) of a side of the
%                                   larger region of interest (cube-shaped)
%                                   to be used in the core search of WM
%                                   tracts.
%
% Output: ROI_numFileName [String]  File name of the mask (of the same size
%                                   as the original white matter mask) that
%                                   contains the number of the larger ROI 
%                                   to which each pixel belongs.

Rn = load(ROI_numFileName);
ROI_num = Rn.ROI_num;

tic
Paths = [];
Regions = [];
PathsFiles = Dir2Arr([InDir filesep qPpaths{1}],'Paths*.mat');
PoscsFiles = Dir2Arr([InDir filesep qPpaths{2}],'Poscs*.mat');
numpaths = size(PathsFiles,1);

if ~isdir(OutDir)
    mkdir(OutDir)
end

for i = 1:max(Rn.ROI_num(:))
    eval(['PathsROI_' num2str(i) ' = [];']);
end

for i = 1:numpaths
    if mod(i,100) == 0; disp(['Loading Paths_' num2str(i)]); end
    load(deblank(PathsFiles(i,:)));
    load(deblank(PoscsFiles(i,:)));
    if ~isfield(q,'location')
        for j = 1:length(q)
            if ~isempty(q(j).voxels)
                pix = unique(floor(q(j).voxels)','rows');
                q(j).ROI = nonzeros(unique(ROI_num(sub2ind(size(ROI_num),pix(:,1),pix(:,2),pix(:,3)))));
                q(j).location = [i j];
                q(j).regions = posc(j,:);
            else
                q(j).ROI = [];
                q(j).location = [i j];
                q(j).regions = [];
            end
        end
    end
    if isempty(q) == 0
        Paths = [Paths q];
        Regions = [Regions; posc];
    end
end
toc
save([OutDir filesep 'AllRegions.mat'],'Regions')

for i = 1:length(Paths)
    if mod(i,50000) == 0; disp(['Sorting Path ' num2str(i)]); end
    Pathsj = Paths(i);
    dROI = Paths(i).ROI;
    for k = 1:length(dROI)
        eval(['PathsROI_' num2str(dROI(k)) ' = [PathsROI_' num2str(dROI(k)) ' Pathsj];'])
    end
end

for i = 1:max(ROI_num(:))
    if mod(i,100) == 0; disp(['Saving Paths in Region ' num2str(i)]); end
    PathsROI = [];
    eval(['PathsROI = PathsROI_' num2str(i) ';']);
    eval(['save (''' OutDir filesep 'PathsInROI_' num2str(i) ''', ''PathsROI'');']);
end
PathsByRegion = OutDir;