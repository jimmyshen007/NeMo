function [newDiffFileNames, newGMMaskFileName, newWMMaskFileName, newContMaskFileName, newAtlasFileName] = resampling2DWMRI_MNI_SVN(DiffFileNames, GMMaskFileName, WMMaskFileName, ContMaskFileName, AtlasFileName,T1FileName,MRIOPT)
% added code to reslice DTI volumes to isotropic volumes before calling Yasser's resampling2DWI
startup_varsonly;
%eve_tools='/home/eve/Documents/MATLAB/eve_tools'; addpath(eve_tools);

if nargin==0
    [DiffFileNames,DiffFilePath] = myuigetfiles('*.img; *.nii','Select ALL DTI images, make sure 1st one is b0...');
    [GMMaskFileName,GMMaskFilePath] = uigetfile({'*.img'; '*.nii'},'Reading Mask_Gray ...');
    [WMMaskFileName,WMMaskFilePath] = uigetfile({'*.img'; '*.nii'},'Reading Mask_White ...');
    [ContMaskFileName,ContMaskFilePath] = uigetfile({'*.img'; '*.nii'},'Reading Mask_Contrast ...');
    [AtlasFileName,AtlasFilePath] = uigetfile({'*.img'; '*.nii'},'Reading Atlas image ...');
    GMMaskFileName = [GMMaskFilePath GMMaskFileName];
    WMMaskFileName = [WMMaskFilePath WMMaskFileName];
    ContMaskFileName = [ContMaskFilePath ContMaskFileName];
    AtlasFileName = [AtlasFilePath AtlasFileName];
    for i = 1:length(DiffFileNames)
        DiffFileNames{i} = [DiffFilePath DiffFileNames{i}];
    end
end

%% reslice diffusion images to make them isotropic
[p,f,e] = fileparts(DiffFileNames{1});
minisovoxeldim=2.5; % will be resampled to this vox size if larger than this

if  strcmp(f(1:3),'iso') 
    for i = 1:length(DiffFileNames)
        [DiffFilePt, DiffFileNm, ext] = fileparts(DiffFileNames{i});
        newDiffFileNames{i} = fullfile(DiffFilePt, ['iso-' DiffFileNm ext]);
        tmphdr = load_nii_hdr(DiffFileNames{i});
        %reslice_nii(DiffFileNames{i}, newDiffFileNames{i}, (1/50)*round(50*max(tmphdr.dime.pixdim(2:4)))*ones(1,3));
        reslice_nii(DiffFileNames{i}, newDiffFileNames{i}, (1/50)*round(50*min(max(tmphdr.dime.pixdim(2:4)),minisovoxeldim))*ones(1,3));
    end
else
    newDiffFileNames = DiffFileNames;
end
%img_auto_translate(newDiffFileNames','backup_',0);
%%-------------------------------------------------------------------------- 
% Authors: Yasser Iturria Medina & Pedro Vald�s-Hern�ndez 
% Neuroimaging Department
% Cuban Neuroscience Center
% November 27th 2007
% Version $1.0

% [b0FileName,b0FilePath] = uigetfile({'*.img'; '*.nii'},'Select b0 image...');
% b0Image = [b0FilePath b0FileName];
b0Image = newDiffFileNames{1};
%img_auto_translate(b0Image,'backup_',0,0);
warning off

%Align origins of files for easier SPM processing
%img_auto_translate({GMMaskFileName; WMMaskFileName; ContMaskFileName; AtlasFileName},'backup_',0);
% change_space(GMMaskFileName,b0Image,0);
% change_space(WMMaskFileName,b0Image,0);
% change_space(ContMaskFileName,b0Image,0);
% change_space(AtlasFileName,b0Image,0);

%Make the new files to pass into the coreg.
[pp, nn, xx] = fileparts(GMMaskFileName);
newGMMaskFileName = fullfile(pp, ['s' nn xx]);
copyfile(GMMaskFileName,newGMMaskFileName)
if strcmp(xx,'.img')
    copyfile([pp filesep nn '.hdr'],[pp filesep 's' nn '.hdr'])
end
[pp, nn, xx] = fileparts(WMMaskFileName);
newWMMaskFileName = fullfile(pp, ['s' nn xx]);
copyfile(WMMaskFileName,newWMMaskFileName)
if strcmp(xx,'.img')
    copyfile([pp filesep nn '.hdr'],[pp filesep 's' nn '.hdr'])
end
[pp, nn, xx] = fileparts(ContMaskFileName);
newContMaskFileName = fullfile(pp, ['s' nn xx]);
copyfile(ContMaskFileName,newContMaskFileName)
if strcmp(xx,'.img')
    copyfile([pp filesep nn '.hdr'],[pp filesep 's' nn '.hdr'])
end
[pp, nn, xx] = fileparts(AtlasFileName);
newAtlasFileName = fullfile(pp, ['s' nn xx]);
copyfile(AtlasFileName,newAtlasFileName)
if strcmp(xx,'.img')
    copyfile([pp filesep nn '.hdr'],[pp filesep 's' nn '.hdr'])
end
[pp, nn, xx] = fileparts(T1FileName);
newT1FileName = fullfile(pp, ['s' nn xx]);
copyfile(T1FileName,newT1FileName)
if strcmp(xx,'.img')
    copyfile([pp filesep nn '.hdr'],[pp filesep 's' nn '.hdr'])
end

%Coreg start...
matlabbatch=create_coreg_job_struct(b0Image, {newGMMaskFileName; newWMMaskFileName; newContMaskFileName; newAtlasFileName; newT1FileName});
%Initialize SPM for the command line
spm('defaults','PET');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
%... coreg end

%Overwrite the newly created files that have 'ss' prefix with the new file 
%names with the 's' prefix.
[pp, nn, xx] = fileparts(GMMaskFileName);
movefile([pp filesep 'ss' nn xx],[pp filesep 's' nn xx])
if strcmp(xx,'.img')
    movefile([pp filesep 'ss' nn '.hdr'],[pp filesep 's' nn '.hdr'])
end
[pp, nn, xx] = fileparts(WMMaskFileName);
movefile([pp filesep 'ss' nn xx],[pp filesep 's' nn xx])
if strcmp(xx,'.img')
    movefile([pp filesep 'ss' nn '.hdr'],[pp filesep 's' nn '.hdr'])
end
[pp, nn, xx] = fileparts(ContMaskFileName);
newContMaskFileName = fullfile(pp, ['s' nn xx]);
movefile([pp filesep 'ss' nn xx],[pp filesep 's' nn xx])
if strcmp(xx,'.img')
    movefile([pp filesep 'ss' nn '.hdr'],[pp filesep 's' nn '.hdr'])
end
[pp, nn, xx] = fileparts(AtlasFileName);
movefile([pp filesep 'ss' nn xx],[pp filesep 's' nn xx])
if strcmp(xx,'.img')
    movefile([pp filesep 'ss' nn '.hdr'],[pp filesep 's' nn '.hdr'])
end
[pp, nn, xx] = fileparts(T1FileName);
movefile([pp filesep 'ss' nn xx],[pp filesep 's' nn xx])
if strcmp(xx,'.img')
    movefile([pp filesep 'ss' nn '.hdr'],[pp filesep 's' nn '.hdr'])
end

if MRIOPT == 1;
    openMRIcron([b0Image ' -b 20 -o ' newGMMaskFileName ' -o ' newWMMaskFileName]);
end

