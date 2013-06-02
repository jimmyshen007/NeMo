function [GMSegFile, WMSegFile, CSSegFile, AtlasedFile] = Atlasing_spm8(Images, Output_dir, AtlasFile, TempFile, Thresh,SPMdir)
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
% Authors: Lester Melie Garcia & Yasser Alem�n G�mez
% Neuroimaging Department
% Cuban Neuroscience Center
% November 15th 2005
% Version $1.0

spm_get_defaults;
warning off;
global defaults;
dseg                   = defaults.preproc;
dseg.output.GM = [0 0 1];
dseg.output.WM = [0 0 1];
dseg.output.CSF = [0 0 1];
dseg.tpm   = char(...
               fullfile(SPMdir,'tpm','r1mm_grey.nii'),...
               fullfile(SPMdir,'tpm','r1mm_white.nii'),...
               fullfile(SPMdir,'tpm','r1mm_csf.nii'));
 
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
%H = waitbar(0,['Subjects  ' num2str(Ns)],'Resize','on','Position',[233.25 237.75 273 50.25],'Resize','off');
for i=1:Ns
   % waitbar(i/Ns,H,['Subject  ' num2str(i) ' of ' num2str(Ns)]);
    disp(['Case ---> ' num2str(i)]);
    [pth, nm, xt] = fileparts(V(i).fname);
    ExOut = exist('Output_dir');
    if (nargin<2)|(ExOut == 0)|(isempty(Output_dir))
        Output_dir(i,:) = [pth filesep 'Atlased'];
        mkdir(pth,'Atlased');
    end
    %%%%%%%%%--------------- Normalization ----------------%%%%%%%%%%%%%%%%
    %
   % disp(['Normalising ...']);
    mkdir(pth,'Normalized');
   % matname = [spm_str_manip([pth filesep 'Normalized' filesep nm],'sd') '_sn.mat'];
  %  spm_normalise(VG0, V(i), matname,dnrm.estimate.weight,'',dnrm.estimate);
  %  spm_write_sn(V(i).fname, matname, dnrm.write);


    %%%%%%%%%-------- End of Normalization Step -----------%%%%%%%%%%%%%%%%
    %
    %%%%%%%%%--------------- Segmentation -----------------%%%%%%%%%%%%%%%%
    %
    disp(['Normalizing and Segmenting ...']);
    mkdir(pth,'Segmented');
    V0 = spm_preproc(V(i),dseg);
    % convert segmentation information to sn-files
    [po,pin] = spm_prep2sn(V0);
    
    opts2.biascor   = 0;         %write bias corrected images    
    opts2.GM        = [0 0 1];   %write modulated, unmodulated normalized, and native GM
    opts2.WM        = [0 0 1];   %write modulated, unmodulated normalized, and native WM
    opts2.CSF       = [0 0 1];   %write modulated, unmodulated normalized, and native CSF
    opts2.cleanup   = 0;         %Do not cleanup

    % write segmented data
    spm_preproc_write(po,opts2);
    
    %save the segmentation information
    [pth fn ext] = fileparts(V(i).fname);
    save(fullfile([pth filesep 'Normalized'], [fn '_seg_sn.mat']), '-struct', 'po');
    save(fullfile([pth filesep 'Normalized'], [fn '_seg_inv_sn.mat']), '-struct', 'pin');
    
    %Write the normalized T1 image to the subdirectory. 
    spm_write_sn(V(i).fname,po);
    if ext == '.img'
        movefile([pth filesep 'w' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'w' fn '.hdr'],[pth filesep 'Normalized']);
    elseif xt == '.nii'
        movefile([pth filesep 'w' fn ext],[pth filesep 'Normalized']);
    end
    
    %Save the segmented images to the subdirectory.
    if ext == '.img'
        movefile([pth filesep 'c1' fn ext],[pth filesep 'Segmented']);
        movefile([pth filesep 'c1' fn '.hdr'],[pth filesep 'Segmented']);
        movefile([pth filesep 'c2' fn ext],[pth filesep 'Segmented']);
        movefile([pth filesep 'c2' fn '.hdr'],[pth filesep 'Segmented']);
        movefile([pth filesep 'c3' fn ext],[pth filesep 'Segmented']);
        movefile([pth filesep 'c3' fn '.hdr'],[pth filesep 'Segmented']);
    elseif xt == '.nii'
        movefile([pth filesep 'c1' fn ext],[pth filesep 'Segmented']);
        movefile([pth filesep 'c2' fn ext],[pth filesep 'Segmented']);
        movefile([pth filesep 'c3' fn ext],[pth filesep 'Segmented']);
    end
    %define the images to normalize
    grey = fullfile([pth filesep 'Segmented'], ['c1' fn ext]);
    white = fullfile([pth filesep 'Segmented'], ['c2' fn ext]);
    csf = fullfile([pth filesep 'Segmented'], ['c3' fn ext]);
    Vgrey = spm_vol(grey);
    Vwhite = spm_vol(white);
    Vcsf = spm_vol(csf);
    
    %write unmodulated normalized
    spm_write_sn(Vgrey,po);
    spm_write_sn(Vwhite,po);
    spm_write_sn(Vcsf,po);
    if ext == '.img'
        movefile([pth filesep 'wc1' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'wc1' fn '.hdr'],[pth filesep 'Normalized']);
        movefile([pth filesep 'wc2' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'wc2' fn '.hdr'],[pth filesep 'Normalized']);
        movefile([pth filesep 'wc3' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'wc3' fn '.hdr'],[pth filesep 'Normalized']);
        movefile([pth filesep 'mwc1' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'mwc1' fn '.hdr'],[pth filesep 'Normalized']);
        movefile([pth filesep 'mwc2' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'mwc2' fn '.hdr'],[pth filesep 'Normalized']);
        movefile([pth filesep 'mwc3' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'mwc3' fn '.hdr'],[pth filesep 'Normalized']);
    elseif xt == '.nii'
        movefile([pth filesep 'wc1' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'wc2' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'wc3' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'mwc1' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'mwc2' fn ext],[pth filesep 'Normalized']);
        movefile([pth filesep 'mwc3' fn ext],[pth filesep 'Normalized']);
    end


    %%%%%%%%%--------- End of Segmentation Step -----------%%%%%%%%%%%%%%%%
    %
    %%%%%%%%%---------------  Atlasing --------------------%%%%%%%%%%%%%%%%
    disp('Atlasing ...');
    GMSegFile{i} = [pth filesep 'Segmented' filesep 'c1' fn ext];
    WMSegFile{i} = [pth filesep 'Segmented' filesep 'c2' fn ext];
    CSSegFile{i} = [pth filesep 'Segmented' filesep 'c3' fn ext];
    Transf_matname = fullfile([pth filesep 'Normalized'], [fn '_seg_sn.mat']);
    AtlasedFile{i} = Auto_Labelling(GMSegFile{i}, WMSegFile{i}, CSSegFile{i}, AtlasFile, Transf_matname, Output_dir(i,:),Thresh);
    %%%%%%%%%--------- End of the Atlasing Step -----------%%%%%%%%%%%%%%%%
end;
GMSegFile = char(GMSegFile);
WMSegFile = char(WMSegFile);
CSSegFile = char(CSSegFile);
AtlasedFile = char(AtlasedFile);
%close(H);
%========================End of main program==============================%
return;

% replaced internal fn below choosing... by m-file... - AR
