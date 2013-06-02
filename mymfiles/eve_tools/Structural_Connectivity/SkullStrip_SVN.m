function [nskT1scan] = SkullStrip_SVN(Images, GMSegFile, WMSegFile, CSFSegFile,MRIOPT)
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


%=========================Main program=====================================
V = spm_vol(Images);
Ns = length(V); % Number of subjects (Analyze files)
for i=1:Ns
    [pth, nm, xt] = fileparts(V(i).fname);
    nskT1scan = [pth filesep 'nsk' nm xt];
    Vnsk = V(i);
    Vnsk.fname = nskT1scan;
    Vnsk.private.dat.fname = nskT1scan;
    nskT1 = spm_read_vols(V(i)).*((spm_read_vols(spm_vol(GMSegFile))+spm_read_vols(spm_vol(WMSegFile))+spm_read_vols(spm_vol(CSFSegFile)))>0.5);
    spm_write_vol(Vnsk,nskT1);
    %Update the WM, GM, and CSF segmentations to exclue skull voxels
    wminfo = spm_vol(WMSegFile);
    wm = spm_read_vols(wminfo);
    wm(nskT1==0) = 0;
    spm_write_vol(wminfo,wm);
    gminfo = spm_vol(GMSegFile);
    gm = spm_read_vols(gminfo);
    gm(nskT1==0) = 0;
    spm_write_vol(gminfo,gm);
    csfinfo = spm_vol(CSFSegFile);
    csf = spm_read_vols(csfinfo);
    csf(nskT1==0) = 0;
    spm_write_vol(csfinfo,csf);
end

if MRIOPT == 1;
    openMRIcron([Images ' -b 20 -o ' GMSegFile ' -o ' WMSegFile]);
end
%========================End of main program==============================%
return;

