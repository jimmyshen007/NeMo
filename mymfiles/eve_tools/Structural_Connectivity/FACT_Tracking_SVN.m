function FACT_Tracking_SVN(tracking,seed,seed_zone,Tzones,Zones,Diff_name,inf_name,Mask_name,Atlas_name,step,longMax,NumXfiles,repetitions,ang,numVectors,OutFilePath,Ppaths,seed_begin)
%
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
MaskVol.pinfo = [1 0 0]';
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
AtlasVol.pinfo = [1 0 0]';
Atlas = round(spm_read_vols(AtlasVol));
if nargin < 9, step = 0.5; end
if nargin < 10, longMax = pi*size(Mask,1); end
if nargin < 11, NumXfiles = 10; end

factor = [min(size_voxel)/size_voxel(1) min(size_voxel)/size_voxel(2) min(size_voxel)/size_voxel(3)]';
q = struct('voxels',[],'PT',[]); posc = []; f = 0; j0 = 1; 
path = OutFilePath;

num_seed = size(seed,1);
cwd = cd;

load(Diff_name);
%             for v = 1:9
%                 E0Vol(v).pinfo = [1 0 0]';
%             end
Diff_Data = VectorF;
S = size(Diff_Data);
if nargin < 6
    [filename, pathname] = uigetfile({'D*.img';'D*.nii'},'If you want to interpolate, select the diffusion tensor image ...');
    if filename ~= 0
        name = [pathname, filename];
        inf_name = [name(ones(6,1),:) repmat(',',6,1) num2str((1:6)')];
    end
end
if ~isempty(inf_name)
    D0Vol = spm_vol(inf_name);
    %                 for v = 1:6
    %                     D0Vol(v).pinfo = [1 0 0]';
    %                 end
    D0 = spm_read_vols(D0Vol);
    S = size(D0);
    inf_Data = zeros([3 3 S(1:3)]);
    inf_Data(1,1,:,:,:) = D0(:,:,:,1); inf_Data(1,2,:,:,:) = D0(:,:,:,2);
    inf_Data(1,3,:,:,:) = D0(:,:,:,3); inf_Data(2,2,:,:,:) = D0(:,:,:,4);
    inf_Data(2,3,:,:,:) = D0(:,:,:,5); inf_Data(3,3,:,:,:) = D0(:,:,:,6);
    inf_Data(2,1,:,:,:) = inf_Data(1,2,:,:,:); inf_Data(3,1,:,:,:) = inf_Data(1,3,:,:,:);
    inf_Data(3,2,:,:,:) = inf_Data(2,3,:,:,:);
    clear D0 name
    interp = 1;
else
    interp = 0;
end
cwd = cd;
cd(path)
if ~exist(Ppaths{1},'dir') %Paths dir
    mkdir(Ppaths{1})
end
if ~exist(Ppaths{2},'dir') %Poscs dir
    mkdir(Ppaths{2})
end

savestr = (floor(seed_begin/NumXfiles)+1):1:(floor(num_seed/NumXfiles)+1);
k = 0;

for j = j0:num_seed
    if mod(j,1000)==0; disp(['Seed ' num2str(j) ' out of ' num2str(num_seed)]); end
    seed_r = seed(j,:);
    seed_r_zone = seed_zone(j);
    for dir = [1 -1]
        [r,a] = pathM(seed_r,seed_r_zone,step,longMax,factor,dir,interp);
        if (~isempty(r)) && (ismember(a(2),rA_zones))
            f = f + 1;
            L = size(r,2);
            q(f).voxels = r;
            q(f).PT = [];
            posc(f,:) = int16(a);
        end
    end
    if (~mod(j,NumXfiles)) || (j == size(seed,1))
        k = k + 1;
        eval(['save ' Ppaths{1} filesep 'Paths_' num2str(savestr(k)) ' q']);
        eval(['save ' Ppaths{2} filesep 'Poscs_' num2str(savestr(k)) ' posc']);
        q = struct('voxels',[],'PT',[]); posc = []; f = 0;
    end
end
cd(cwd)