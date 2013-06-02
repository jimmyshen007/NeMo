function [GMSegFile, WMSegFile, CSSegFile, AtlasedFile] = Atlasing(Images, Output_dir, AtlasFile, TempFile, Thresh)
%
% Syntax :
% atlasing(Images, Output_dir, AtlasFile, TempFile, Thresh);
%
% This script was developed over SPM2 toolbox and it computes automatically
% individual atlases based on Magnetic Resonance Images.
% The first step is based on the MRI image normalization to a stereotaxic
% space, MNI space(Montreal Neurological Institute). Here a transformations
% matrix is obtained. Then the MRI individual files are segmented in three
% different brain tissues (cerebral spinal fluid, gray and white matter) at
% this stage. During the second step each gray matter voxel is labeled with
% one structure label using the transformation matrix obtained in the
% normalization process and an anatomical atlas constructed by manual segmentation
% for a group of subjects.
%
% Input Parameters:
%   Images     : Individual MRI files.
%   Output_dir : Output directory for segmented, normalized and atlased
%                files. If the user doesn't change the output directory,
%                the resulting files are saved in the same address than the
%                individual MRI files.
%   AtlasFile  : Reference Atlas File used in the automatic labelling step.
%   TempFile   : Template file for the normalisation step. The user can use
%                any of the templates included at the SPM package or select
%                another, always taking care that this template file is in
%                the same space than the atlas that will be used in the
%                labelling step.
%   Thresh     : Threshold for gray matter segmentation.
%                Just the voxels with  higher probability than the threshold
%                are taken into acount in the automatic labelling step. If
%                the threshold isn't specified then an automatic one is taken.
%                All voxels with higher gray matter probabillity than 1-(GM+WM+CSF)
%                are taken into account in the automatic labelling step.
%                Being:
%                GM(A voxel V belongs to Gray Matter tissue with a probability GM)
%                WM(A voxel V belongs to White Matter tissue with a probability WM)
%                CSF(A voxel V belongs to Cerebral Spinal Fluid with a probability CSF)
%
% Related references:
% 1.- Ashburner J, Friston K. Multimodal image coregistration and partitioning--
%     a unified framework. Neuroimage. 1997 Oct;6(3):209-17.
% 2.- Voxel-based morphometry--the methods. Neuroimage. 2000 Jun;11(6 Pt
% 1):805-21.
% 3.- Evans AC, Collins DL, Milner B (1992). An MRI-based Stereotactic Brain
%     Atlas from 300 Young Normal Subjects, in: Proceedings of the 22nd Symposium
%     of the Society for Neuroscience, Anaheim, 408.
%
% See also: spm_normalise  spm_segment  Auto_Labelling
%__________________________________________________
% Authors: Lester Melie Garcia & Yasser Alemán Gómez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 15th 2005
% Version $1.0

spm_defaults
warning off;
global defaults
dseg                   = defaults.segment;
dseg.write.wrt_cor     = 0;
dnrm                   = defaults.normalise;
dnrm.write.vox         = [1 1 1];
dnrm.write.bb          = reshape([-90 90 -126 90 -72 108],2,3); %bounding box template
dnrm.write.preserve    = 0;
dnrm.estimate.graphics = 0;   % Crashes otherwise
dnrm.estimate.smoref  = 8;
dnrm.estimate.smosrc  = 8;
cd =which('IBASPM_Atlasing');
[cd nm ext] = fileparts(cd);
%=====================Checking Input Parameters===========================%
if nargin==0
    [Images,sts] = spm_select([1 Inf],'image','Selecting UnNormalized Images','',cd);
    [NFiles,m] = size(Images);
    [AtlasFileName,AtlasFilePath] = uigetfile({'*.img'},'Reading Reference Atlas File ...');
    AtlasFile = [AtlasFilePath AtlasFileName];
    [TempFile] = Choosing_TempFile;
    Thresh = input('Please select a threshold for gray matter segmentation:   ');
else
    if isempty(Images)
        [Images,sts] = spm_select([1 Inf],'image','Selecting UnNormalized Images','',cd);
    end
    if isempty(AtlasFile)
        [AtlasFileName,AtlasFilePath] = uigetfile('*.img','Reading Reference Atlas File ...');
        AtlasFile = [AtlasFilePath AtlasFileName];
    end
    if isempty(TempFile)
        [TempFile] = Choosing_TempFile;
    end
    if ~exist('Thresh', 'var')
        Thresh = input('Please select a threshold for gray matter segmentation:   ');
    end
end
%=========================================================================%
%
%=========================Main program=====================================
VG0 = spm_vol(TempFile);
V = spm_vol(Images);
Output_dir = char(Output_dir);
Ns = length(V); % Number of subjects (Analyze files)
H = waitbar(0,['Subjects  ' num2str(Ns)],'Resize','on','Position',[233.25 237.75 273 50.25],'Resize','off');
for i=1:Ns
    waitbar(i/Ns,H,['Subject  ' num2str(i) ' of ' num2str(Ns)]);
    disp(['Case ---> ' num2str(i)]);
    [pth, nm, xt] = fileparts(V(i).fname);
    ExOut = exist('Output_dir');
    if (nargin<2)|(ExOut == 0)|(isempty(Output_dir))
        Output_dir(i,:) = [pth filesep 'Atlased'];
        mkdir(pth,'Atlased');
    end
    %%%%%%%%%--------------- Normalization ----------------%%%%%%%%%%%%%%%%
    %
    disp(['Normalising ...']);
    mkdir(pth,'Normalized');
    matname = [spm_str_manip([pth filesep 'Normalized' filesep nm],'sd') '_sn.mat'];
    spm_normalise(VG0, V(i), matname,dnrm.estimate.weight,'',dnrm.estimate);
    spm_write_sn(V(i).fname, matname, dnrm.write);
    if xt == '.img'
        movefile([pth filesep 'w' nm xt],[pth filesep 'Normalized']);
        movefile([pth filesep 'w' nm '.hdr'],[pth filesep 'Normalized']);
    elseif xt == '.nii'
        movefile([pth filesep 'w' nm xt],[pth filesep 'Normalized']);
    end

    %%%%%%%%%-------- End of Normalization Step -----------%%%%%%%%%%%%%%%%
    %
    %%%%%%%%%--------------- Segmentation -----------------%%%%%%%%%%%%%%%%
    %
    disp(['Segmenting ...']);
    mkdir(pth,'Segmented');
    V0 =spm_segment(V(i),VG0,dseg);
    Vol = struct('fname','','dim','','mat',V(i).mat,'pinfo',[1/255 0 0]',...
        'descrip','Atlas image','dt',[2 0]);
    if strcmp(spm('ver'), 'SPM2')
        dim = [V(i).dim(1:3) Vol.dt(1)];
     elseif strcmp(spm('ver'), 'SPM5')
        dim = [V(i).dim(1:3)];
    end
    Vol.dim = dim;
    for j=1:3,
        tmp   = fullfile([pth filesep 'Segmented'],[nm,'_seg',num2str(j), xt]);
        VO(j) = Vol;
        VO(j).fname = tmp;
    end;
    VO = spm_create_vol(VO);
    for pp=1:V(i).dim(3),
        VO(1) = spm_write_plane(VO(1),double(V0(1).dat(:,:,pp))/255,pp);
        VO(2) = spm_write_plane(VO(2),double(V0(2).dat(:,:,pp))/255,pp);
        VO(3) = spm_write_plane(VO(3),double(V0(3).dat(:,:,pp))/255,pp);
    end;
    clear VO; fclose('all');
    %%%%%%%%%--------- End of Segmentation Step -----------%%%%%%%%%%%%%%%%
    %
    %%%%%%%%%---------------  Atlasing --------------------%%%%%%%%%%%%%%%%
    disp('Atlasing ...');
    GMSegFile{i} = [pth filesep 'Segmented' filesep nm '_seg1' xt];
    WMSegFile{i} = [pth filesep 'Segmented' filesep nm '_seg2' xt];
    CSSegFile{i} = [pth filesep 'Segmented' filesep nm '_seg3' xt];
    Transf_matname = matname;
    AtlasedFile{i} = Auto_Labelling(GMSegFile{i}, WMSegFile{i}, CSSegFile{i}, AtlasFile, Transf_matname, Output_dir(i,:),Thresh);
    %%%%%%%%%--------- End of the Atlasing Step -----------%%%%%%%%%%%%%%%%
end;
GMSegFile = char(GMSegFile);
WMSegFile = char(WMSegFile);
CSSegFile = char(CSSegFile);
AtlasedFile = char(AtlasedFile);
close(H);
%========================End of main program==============================%
return;

% replaced internal fn below choosing... by m-file... - AR
