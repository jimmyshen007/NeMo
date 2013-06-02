function amy_Tracking3(tracking,seed,seed_zone,Tzones,Zones,Diff_name,inf_name,Mask_name,Atlas_name,step,longMax,NumXfiles,repetitions,ang,numVectors,OutFilePath)
%amy updated it to use the package spm5 as opposed to spm2
% Usage : Tracking1(tracking,seed,Tzones,Zones,Diff_name,inf_name,Mask_name,
%             Atlas_name,step,longMax,NumXfiles,repetitions,ang,numVectors)
%
% Description: This function compute nervous fibers trajectories front a set
% of seed points to a set of target gray matter zones using three different
% methods ("FACT", "DTensor", "BayesianF"). Also, for the two last methods,
% it compute the probability or validity index of each path.
%
% Input Parameters:
%        tracking   : variable that specifies the method that will be used
%                     ("FACT", "DTensor" or "BayesianF") to the layout of
%                     nervous fibers trajectories.
%        seed       : matrix [n,3] that contains the coordinates of the n
%                     seed points.
%        Tzones     : vector that contains the numbers of the target zones.
%        Zones      : variable struct that contains information about the defined
%                     zones.
%        Diff_name  :
%        inf_name   :
%        Mask_name  :
%        Atlas_name :
%        step       : step size.
%        longMax    : maximum number of steps that can have each path.
%        NumXfiles  : number of seed points in each file to save.
%        repetitions: for the probabilistic methods, it indicates the number
%                     of paths that begin in each seed point.
%        ang        : maximum angle of bend that can have each path
%                     between two serial steps.
%        numVectors : number of random vectors to be generate in each point
%                     (only for "Dtensor" method).
%
%        It is also necessary: the atlas image, the mask for tracking and
%        diffusion data whose structure depends from the tracking method to use
%        accordind to:
%                     - "FACT", the eigenvectors image E_�subject�.img, , and,
%                       if you want to interpolate, the diffusion tensor image
%                       D_'subject'.img.
%                     - "Dtensor", the diffusion tensor image D_'subject'.img and
%                       the eigenvectors image E_'subject'.img.
%                     - "BayesianF", the ODF image ODF_'subject'.img and a txt
%                       or cod file whit a [N 3] matrix that contain the N
%                       directions of the ODF.
%
% Output Parameters: This function saves two different types of files. For
% example: The file Paths_35 contains the paths calculated for the group 35
% of seeds points, and also, for the probabilistic methods, it contains the
% probability or validity index of each path. Poscs_35 contain the gray matter
% zones of which leave and arrive these paths.
%--------------------------------------------------------------------------
% Useful References:
%  1.-Mori S., Crain J. B., Van Zijl P. C. M. y Chackov V. P. Three Dimensional
% Tracking of Axonal Projections in the Brain by Magnetic Resonance Imaging.
% Ann. Neurol, 45, 263-269, 1999.
%  2.-Y. Iturria-Medina, E. Canales-Rodr�guez, L. Meli�-Garc�a,
% P. Vald�s-Hern�ndez. Bayesian formulation for fiber tracking. Presented
% at the 11th Annual Meeting of the Organization for Human Brain Mapping,
% June 12-16, 2005, Toronto, Ontario, Canada. Available on CD-Rom in
% NeuroImage, Vol. 26, No.1.
%--------------------------------------------------------------------------
% Authors: Yasser Iturria Medina, Erick Canales Rodr�guez & Pedro Vald�s-Hern�ndez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 15th 2005
% Version $1.0

global Diff_Data inf_Data Mask Atlas

for i = 1:length(Tzones)
    rA_zones(i) = Zones(Tzones(i)).number;
end
if nargin < 7
    [MaskFileName,MaskFilePath] = uigetfile({'*.img';'*.nii'},'Reading Mask_for_tracking image ...');
    Mask_name = [MaskFilePath MaskFileName];
end
MaskVol = spm_vol(Mask_name);
%MaskVol.pinfo = [1 0 0]';
Mask = spm_read_vols(spm_vol(MaskVol));
if nargin < 8
    [AtlasFileName,AtlasFilePath] = uigetfile({'*.img';'*.nii'},'Reading Atlas image ...');
    Atlas_name = [AtlasFilePath AtlasFileName];
end
AtlasVol = spm_vol(Atlas_name);
if spm('ver')=='SPM2'
    size_voxel = AtlasVol.private.hdr.dime.pixdim(2:4);
elseif (strcmp(spm('ver'),'SPM5') || strcmp(spm('ver'),'SPM8'))
    mat = abs(diag(AtlasVol.mat));
    size_voxel = sqrt(sum(AtlasVol.mat(1:3,1:3).^2));
end
%AtlasVol.pinfo = [1 0 0]';
Atlas = round(spm_read_vols(AtlasVol));

if nargin < 9, step = 0.5; end
if nargin < 10, longMax = pi*size(Mask,1); end
if nargin < 11, NumXfiles = 10; end

factor = [min(size_voxel)/size_voxel(1) min(size_voxel)/size_voxel(2) min(size_voxel)/size_voxel(3)]';
q = struct('voxels',[],'PT',[]); posc = []; f = 0; k = 0; initial = 1; j0 = 1; 
path = OutFilePath;
% directory = 1;
% while directory
%     path = uigetdir('Select directory to save paths...');
%     NumbFiles = 0;%size(spm_select('files',path,'Paths_*.mat'),1);
%     directory = 0;
%     if NumbFiles ~= 0
%         quest=questdlg('Paths calculated in the suitable directory already exist. Does you want to continue saving in this directory and to kept the previous calculation?');
%         if strcmp(quest,'Yes'),
%             k = NumbFiles;
%             j0 = NumbFiles*NumXfiles + 1;
%             directory = 0;
%         end
%         if strcmp(quest,'No'), directory = 1; end
%         if strcmp(quest,'Cancel'), return; end
%     end
% end
num_seed = size(seed,1);
cwd = cd;

D0Vol = spm_vol(Diff_name);
Diff_Data  = spm_read_vols(D0Vol);
Ns = size(Diff_name,1);
inf_Data = zeros([Ns 3]);
fid = fopen(inf_name);
for i=1:Ns
    inf_Data(i,:) = str2num(fgetl(fid));
end
cd(path)
if ~exist('Paths','dir')
    mkdir('Paths')
end
if ~exist('Poscs','dir')
    mkdir('Poscs')
end

for j = 1:num_seed
    if mod(j,1000)==0; disp([num2str(j) ' -> ' num2str(num_seed)]); end
    seed_r = seed(j,:);
    seed_r_zone = seed_zone(j);
    for i = 1:repetitions
        [r,a,Ptm_T] = amy_pathB(seed_r,seed_r_zone,longMax,factor,ang);
        if (~isempty(r)) && (ismember(a(2),rA_zones))
            f = f + 1;
            L = size(r,2);
            q(f).voxels = r;
            q(f).PT = single(Ptm_T);
            posc(f,:) = int16(a);
        end
    end
    if (~mod(j,NumXfiles)) || (j == size(seed,1))
        k = k + 1;
        eval(['save Paths' filesep 'Paths_' num2str(k) ' q']);
        eval(['save Poscs' filesep 'Poscs_' num2str(k) ' posc']);
        q = struct('voxels',[],'PT',[]); posc = []; f = 0;
    end
end
cd(cwd)