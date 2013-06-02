function at_vol_mats=RunAtlasing_SVN(T1FileName,atlassize,MRIOPT,translateOPT,temptype)

atspm8=tic;

if nargin<4
    translateOPT=0;
end

startup_varsonly;
% load([MS_img filesep 'T1Files.mat']);

if translateOPT
     img_auto_translate(T1FileName,'backup_',atlassize);
end


%This is the directory of the toolbox.
%SCBdir =  'C:\Documents and Settings\Radiology Guest\Desktop\Atlasing\StructrualConnectivity_Bayes';  
SPMdir = SPM_dir; 

%Procedure to load the T1FileName, can be a list of filenames.
disp('Loading files..............................');

disp(T1FileName);

%The atlas to be used in the parcellation of the GM.
%AtlasTextFileName = [Atlasing_dir filesep 'Subject1_Atlas_WedgeLobe.cod'];
[AtlasFileName,~] = ChooseAtlas(atlassize,Atlasing_dir);

%Leave this variable empty so that the subdirectories and output files are 
%creaed in the patient's T1 directory. 
%StrOutputDir = ['Atlasing' num2str(atlassize)];
StrOutputDir = '';
%This is the treshold for the GM segmentation (usually leave empty).
Thresh = [];
%This is the type of structural image, normally T1, but could also be T2,
%etc.

if nargin<5
temptype = 't1';
end
[TempFileName] = Choosing_TempFile(temptype);

%04/23/2012: Added Amy's spm_coreg function as preprocessing before
%Atlasing_spm8_SVN
%approx_coreg2MNI(deblank(TempFileName),deblank(T1FileName),'backup_');

%Input the T1FileName (with directory) and produce the normalized T1, as
%well as the WM, GM, CSF segmentation and parcellated GM atlas. Note: 
%possible to pass cell array of T1filedies at once
if nargin<3
        MRIOPT=0;
end

[~, ~, ~, AtlasedFileName] = Atlasing_spm8_SVN(T1FileName, StrOutputDir, AtlasFileName, TempFileName, Thresh,SPMdir,MRIOPT,atlassize,1);
toc(atspm8)

%Take each of the atlased T1's and count the number of voxels that have
%been assigned to each of the 116 regions.

at_vol_mats=char(zeros(size(AtlasedFileName,1),MSL));
 for i = 1:size(AtlasedFileName,1)  
    [d,f,~] = fileparts(deblank(AtlasedFileName(i,:)));
    at = spm_read_vols(spm_vol(deblank(AtlasedFileName(i,:))));
    at_vol = zeros(atlassize,1);
    for j = 1:atlassize
        at_vol(j) = length(find(at==j));
    end
    save([d filesep f(1:end-6) '_Vol.mat'],'at_vol');
    at_vol_mats(i,:)=[d filesep f(1:end-6) '_Vol.mat' blanks(MSL-length([d filesep f(1:end-6) '_Vol.mat']))];
 end
 return;
end
