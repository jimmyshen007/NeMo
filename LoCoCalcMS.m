function [All_LoCoMNI,All_LoCoT1] = LoCoCalcMS(DamageFileName,Coreg2MNI,CalcSpace,atlassize,StrSave,NumWorkers,dispMask)
% This function takes in a lesion mask and projects the damage onto the 
% areas of the cortex that have losses in connectivity resulting from that
% pattern of injury. It also outputs the network produced by taking intact
% networks and removing fiber streamlines through those areas of injury.
%--------------------------------------------------------------------------
% Input 
% DamageFileName        A string with the filename of the Damage Mask
% Coreg2MNI             A structure provided if the mask is in native 
%                       (individual) space and needs to be coregistered to 
%                       MNI space. If the mask is already in MNI space,  
%                       provide an empty structure. This struct has the 
%                       following fields: 
%                       ImageFileName: A string with the filename of the 
%                       stuctural image.
%                       StructImageType: A string indicating which image 
%                       type the structural image is - this will be used 
%                       in the normalization routine.
% CalcSpace             A string denoting the space in which to calculate 
%                       the LoCo, either 'MNI','T1', or 'both'.   
% atlassize             A number indicating the cortical atlas you would 
%                       like to use, either 86 (FreeSurfer) or 116 (AAL).
% StrSave               A string identifying the LoCo results file name.
% NumWorkers            The number of parallel workers you'd like to use.
% dispMask              A flag indicating if you'd like to display the
%                       injury mask in MNI/T1 space.
% -------------------------------------------------------------------------
% Output
% OutStrMNI             The filename of the LoCo results (MNI space).
% OutStrT1              The filename of the LoCo results (T1 space).
%
% Written by 
% Amy Kuceyeski
% IDEAL
% Weill Cornell Medical College 
% August 24, 2012  
%--------------------------------------------------------------------------
%Set up the defaults
if nargin < 1
    error('You must at least input a damage mask...')
end
if nargin < 2 || isempty(Coreg2MNI)
    Coreg2MNI = [];
end
if nargin < 3 || isempty(CalcSpace)
    CalcSpace = 'MNI';
end
if nargin < 4 || isempty(atlassize)
    atlassize = 116;
end
if nargin < 5 || isempty(StrSave)
    StrSave = '';
end
if nargin < 6 || isempty(NumWorkers)
    NumWorkers = 2;
end
if nargin < 7 || isempty(dispMask)
    dispMask = 1;
end
%Start the program....
if NumWorkers > 1;
    if matlabpool('size')==0
        eval(['matlabpool open ' num2str(NumWorkers)]);
    elseif ~(matlabpool('size')==NumWorkers);
        matlabpool close
        eval(['matlabpool open ' num2str(NumWorkers)]);
    end
end

%Make sure the files are correctly mapped...
SVNEveTools = '/home/amy/Desktop/Matlab2009a/trunk';
AtlasingDir = '/home/amy/Desktop/Matlab2009a/StructrualConnectivity_Bayes/Atlasing';
SPMdir8 = '/home/amy/Desktop/Matlab2009a/spm8';
niitools = '/home/amy/Desktop/Matlab2009a/StructrualConnectivity_Bayes/nifti_toolbox';
ODFtools = '/home/amy/Desktop/Matlab2009a/StructrualConnectivity_Bayes/ODF_Qball';
SCBtrack = '/home/amy/Desktop/Matlab2009a/StructrualConnectivity_Bayes/Tracking_Bayes';
mymfiles = '/home/amy/Desktop/Matlab2009a/mymfiles';

addpath('/home/amy/Desktop/Matlab2009a/BCT');
addpath(SCBtrack)
addpath(mymfiles)
p = genpath(SVNEveTools);
addpath(p);
p3 = genpath(SPMdir8);
addpath(p3);
addpath(AtlasingDir)
addpath(niitools)
addpath(ODFtools)

[pth,fn,ext] = fileparts(DamageFileName);

%First, coregister the native image and corresponding mask to MNI space, if
%needed.
if ~isempty(Coreg2MNI)
    affinereg = 1;
    [AtlasFileName,AtlasTextFileName]=ChooseAtlas(atlassize,AtlasingDir);
    GMWMcontrast = 1;
    %Leave this variable empty so that the subdirectories and output files are
    %created in the patient's T1 directory.
    StrOutputDir ='';
    %This is the treshold for the GM segmentation (usually leave empty).
    Thresh = [];
    %This is the type of structural image, normally T1, but could also be T2,
    %etc.
    Type = Coreg2MNI.StructImageType;
    [TempFileName] = Choosing_TempFile(Type);
    T1FileName = Coreg2MNI.ImageFileName;
    %Input the T1FileName (with directory) and produce the segmented volumes
    %and the normalized T1 image to MNI space.
    if ~isempty(T1FileName)
        [pth2,fn2,ext2] = fileparts(T1FileName);
        if  ~exist([pth2 filesep 'Normalized' filesep 'wa' fn2 ext2],'file')
            disp('Approximate Affine co-registration...')
            [T1FileName,DamageFileName] = approx_coreg2MNI(TempFileName,T1FileName,'',DamageFileName);
            [~,~,~,~] = Atlasing_spm8_SVN(T1FileName, StrOutputDir, AtlasFileName, TempFileName, Thresh,SPMdir8,1,atlassize,affinereg);
            %Coregister the mask to MNI space using the transformation just calculated.
            nflags = struct('interp',0,'vox',NaN,'bb',NaN,'wrap',[0 0 0],'preserve',0,...
                'prefix','w');
            spm_write_sn(DamageFileName,deblank(Dir2Arr([pth2 filesep 'Normalized'],'*seg_sn.mat')),nflags)
        end
        DamageFileName = [pth filesep 'wa' fn ext];
        openMRIcron([pth2 filesep 'Normalized' filesep 'wa' fn2 ext2 ' -b 20 -o ' DamageFileName]);
    else
        disp(['Please provide a structural image for normalization purposes, or if your image is already '...
            'in MNI space, please provide an empty struct input for the variable Coreg2MNI.'])
    end
else
    [TempFileName] = Choosing_TempFile('t1');
end

results_dir = pth;
if ~isdir(results_dir); mkdir(results_dir); end
results_dir2 = [pth filesep 'DetailedLoCoResults'];
if ~isdir(results_dir2); mkdir(results_dir2); end
main_dir = '/mnt/extra/BTFdata/';
FT_MNI = [main_dir 'FiberTracts' num2str(atlassize) '_MNI'];

PatStr = Dir2Arr(FT_MNI,'e0*');

%% Calculate the LoCo in MNI space
%results_dir = [main_dir 'LoCoCalc/MNIresults'];
ROI_num_filename = [FT_MNI filesep 'w1mm_T1_nzROInum.mat'];
Damage = spm_read_vols(spm_vol(DamageFileName));
%Seed_Damage = spm_read_vols(spm_vol('/home/amy/work/BTFdata/Normals/LoCoCalc/w1444_WML.img')); 
%results_dir = '/mnt/extra/ADdata_Zhang/GroupWiseStats/MNI_LoCo';
%Damage = spm_read_vols(spm_vol([results_dir filesep 'MNIrzwbAD_fdr0p05_clust0.img']));
%StrSave = 'GWAD86';
OutStrMNI = [];
OutStrT1 = [];

All_LoCoMNI = {};
if dispMask
    openMRIcron([TempFileName ' -b 20 -o ' DamageFileName])
end
if strcmp(CalcSpace,'MNI') || strcmp(CalcSpace,'both')
    tic
    LoCoResults = [];
    for i = 1:size(PatStr,1)
        disp(['Calculating LoCo for patient ' num2str(i) ' of ' num2str(size(PatStr,1))])
        PathsByRegion = deblank(PatStr(i,:));
        patstr = PathsByRegion(end-6:end);
        OutPathFileName = [results_dir2 filesep 'Dam_Zones' StrSave '_MNI_' patstr '.mat'];
        %Get the paths that go through the damage mask.
        if ~exist(OutPathFileName,'file')
            [AZ_file,~] = PathDisruption(Damage,ROI_num_filename,OutPathFileName,PathsByRegion);
        else
            AZ_file = OutPathFileName;
        end
        load(AZ_file)
        [b,unique_location,j1] = unique(Affect_location,'rows','first');
        Affzon = Affect_zones(unique_location,:);
        Affloc = Affect_location(unique_location,:);
        load([FT_MNI filesep patstr filesep 'AllRegions.mat'])
        [hh,xx] = hist(Regions(:),1:1:atlassize);
        [hha,xxa] = hist(Affzon(:),1:1:atlassize);
        perc_Affect = hha./hh;
        All_PA = perc_Affect;
        OrigMat = load([FT_MNI filesep 'AllConnMatrices' num2str(atlassize) '_FT' filesep 'Matrix_TFC' num2str(atlassize) '_' patstr '.mat']);
        sAffzon = sort(Affzon,2,'ascend');
        Urows = unique(sAffzon,'rows');
        ConMat = OrigMat.MapTFC;
        for j = 1:size(Urows,1)
            if length(find(Urows(j,:)))==2;
                pc = sum(ismember(sAffzon,Urows(j,:),'rows'));
                ConMat(Urows(j,1),Urows(j,2)) = ConMat(Urows(j,1),Urows(j,2))-pc;
                ConMat(Urows(j,2),Urows(j,1)) = ConMat(Urows(j,1),Urows(j,2));
            end
        end
        nConMat = ConMat./OrigMat.total_fcount;
        LoCoResults(i).NormPatID = patstr;
        LoCoResults(i).Regions = All_PA;
        LoCoResults(i).ConMat = ConMat;
        LoCoResults(i).nConMat = nConMat;
        LoCoResults(i).OrigMat = OrigMat;
        netmet = []; netmet2 = []; netmet3 = []; netmet4 = [];
        D = distance_wei(1./LoCoResults(i).nConMat);
        [netmet(i,1),netmet(i,2),~,netmet(i,3),netmet(i,4)] = charpath(D);
        LoCoResults(i).nConMatMets = [netmet(i,1),netmet(i,2),netmet(i,3),netmet(i,4)];
        D = distance_wei(1./LoCoResults(i).ConMat);
        [netmet2(i,1),netmet2(i,2),~,netmet2(i,3),netmet2(i,4)] = charpath(D);
        LoCoResults(i).ConMatMets = [netmet2(i,1),netmet2(i,2),netmet2(i,3),netmet2(i,4)];
        D = distance_wei(1./LoCoResults(i).OrigMat.MapTFC);
        [netmet3(i,1),netmet3(i,2),~,netmet3(i,3),netmet3(i,4)] = charpath(D);
        LoCoResults(i).OrigMat.MetsMapTFC = [netmet3(i,1),netmet3(i,2),netmet3(i,3),netmet3(i,4)];
        D = distance_wei(1./LoCoResults(i).OrigMat.nMapTFC);
        [netmet4(i,1),netmet4(i,2),~,netmet4(i,3),netmet4(i,4)] = charpath(D);
        LoCoResults(i).OrigMat.nMetsMapTFC = [netmet4(i,1),netmet4(i,2),netmet4(i,3),netmet4(i,4)];
    end
    mLoCo = zeros(size(LoCoResults(1).Regions));
    CM = zeros(size(LoCoResults(1).ConMat));
    nCM = zeros(size(LoCoResults(1).nConMat));
    OM = zeros(size(LoCoResults(1).ConMat));
    nOM = zeros(size(LoCoResults(1).nConMat));
    for i = 1:length(LoCoResults)
        mLoCo = mLoCo + LoCoResults(i).Regions;
        CM = CM + LoCoResults(i).ConMat;
        nCM = nCM + LoCoResults(i).nConMat;
        OM = OM + LoCoResults(i).OrigMat.MapTFC;
        nOM = nOM + LoCoResults(i).OrigMat.nMapTFC;
    end
    sLoCo = [];
    sLoCo.NormPatID = 'Mean Results';
    sLoCo.Regions = mLoCo./length(LoCoResults);
    sLoCo.ConMat = CM./length(LoCoResults);
    sLoCo.nConMat = nCM./length(LoCoResults);
    sLoCo.OrigMat.MapTFC = OM./length(LoCoResults);
    sLoCo.OrigMat.nMapTFC = nOM./length(LoCoResults);
    netmet = []; netmet2 = []; netmet3 = []; netmet4 = [];
    D = distance_wei(1./sLoCo.nConMat);
    [netmet(i,1),netmet(i,2),~,netmet(i,3),netmet(i,4)] = charpath(D);
    sLoCo.nConMatMets = [netmet(i,1),netmet(i,2),netmet(i,3),netmet(i,4)];
    D = distance_wei(1./sLoCo.ConMat);
    [netmet2(i,1),netmet2(i,2),~,netmet2(i,3),netmet2(i,4)] = charpath(D);
    sLoCo.ConMatMets = [netmet2(i,1),netmet2(i,2),netmet2(i,3),netmet2(i,4)];
    D = distance_wei(1./sLoCo.OrigMat.MapTFC);
    [netmet3(i,1),netmet3(i,2),~,netmet3(i,3),netmet3(i,4)] = charpath(D);
    sLoCo.OrigMat.MetsMapTFC = [netmet3(i,1),netmet3(i,2),netmet3(i,3),netmet3(i,4)];
    D = distance_wei(1./sLoCo.OrigMat.nMapTFC);
    [netmet4(i,1),netmet4(i,2),~,netmet4(i,3),netmet4(i,4)] = charpath(D);
    sLoCo.OrigMat.nMetsMapTFC = [netmet4(i,1),netmet4(i,2),netmet4(i,3),netmet4(i,4)];
    LoCoResults = [sLoCo LoCoResults];
    All_LoCoMNI = [results_dir filesep 'LoCo' num2str(atlassize) StrSave '_MNI'];
    save(All_LoCoMNI,'LoCoResults')
    time_for_one_mask = toc
end

%% Calculate the LoCo in T1 space
All_LoCoT1 = {};
if strcmp(CalcSpace,'T1') || strcmp(CalcSpace,'both')
    norm_dir = '/home/amy/work/BTFdata/Normals';
    t1results_dir = [results_dir filesep 'T1_LoCo'];
    if ~isdir(t1results_dir); mkdir(t1results_dir); end
    disp('Starting the T1 LoCo calculation...')
    FT_T1 = [main_dir 'FiberTracts' num2str(atlassize)];
    PatStr = Dir2Arr(FT_T1,'e0*');
    LoCoResults = [];
    tic
    for i = 1:size(PatStr,1)
        disp(['Calculating LoCo for patient ' num2str(i) ' of ' num2str(size(PatStr,1))])
        PathsByRegion = deblank(PatStr(i,:));
        patstr = PathsByRegion(end-6:end);
        OutPathFileName = [t1results_dir filesep 'Dam_Zones_' StrSave 'T1_' patstr '.mat'];
        if ~exist(OutPathFileName)
            ROI_num_filename = [FT_T1 filesep patstr filesep 'ROI_num_' patstr '.mat'];
            NLtransmatMNI2T1 = deblank(Dir2Arr([norm_dir filesep patstr filesep 'T1scan' filesep 'Normalized'],'acos*_seg_inv_sn.mat'));
            LtransmatT12DTI = deblank(Dir2Arr([norm_dir filesep patstr filesep 'T1scan'],'*T12DTI.mat'));
            sT1_NC = deblank(Dir2Arr([norm_dir filesep patstr filesep 'T1scan'],'snsk*.img'));
            [dp,dfn,dext] = fileparts(DamageFileName);
            t1DamageFileName = [dp filesep 'T1_LoCo' filesep dfn dext];
            copyfile(DamageFileName,t1DamageFileName)
            if strcmp(dext,'.img')
                copyfile([dp filesep dfn dext],[dp filesep 'T1_LoCo' filesep dfn '.hdr'])
            end
            %Coregister the patient's lesion mask in MNI space to the normal control's T1.
            w_flags = struct('interp',0,'vox',NaN,'bb',NaN,'wrap',[0 0 0],'preserve',0,...
                'prefix',['a' patstr]);
            spm_write_sn(t1DamageFileName,NLtransmatMNI2T1,w_flags);
            
            [dp,dfn,dext] = fileparts(t1DamageFileName);
            aDamage = [dp filesep 'a' patstr dfn dext];
            [pp, nn, xx] = fileparts(aDamage);
            saDamage = fullfile(pp, ['s' nn xx]);
            copyfile(aDamage,saDamage)
            if strcmp(xx,'.img')
                copyfile([pp filesep nn '.hdr'],[pp filesep 's' nn '.hdr'])
            end
            load(LtransmatT12DTI)
            spm_get_space(saDamage,M*spm_get_space(saDamage))
            rflags.mask = 0; rflags.mean = 0; rflags.interp = 0; rflags.which = 1; rflags.prefix = 's';
            
            spm_reslice({sT1_NC,saDamage},rflags)
            [pp, nn, xx] = fileparts(aDamage);
            movefile([pp filesep 'ss' nn xx],[pp filesep 's' nn xx])
            if strcmp(xx,'.img')
                movefile([pp filesep 'ss' nn '.hdr'],[pp filesep 's' nn '.hdr'])
            end
            if dispMask
                if i == 1
                    openMRIcron([T1 ' -b 20 -o ' t1DamageFileName])
                    pause(1);
                end
                openMRIcron([sT1_NC ' -b 20 -o ' saDamage])
            end
            Load_Damage = spm_read_vols(spm_vol(saDamage));
            
            %Get the paths that go through the damage mask.
            
            [AZ_file,Affected_WMTracts] = PathDisruption(Load_Damage,ROI_num_filename,OutPathFileName,PathsByRegion);
        else
            AZ_file = OutPathFileName;
        end
        load(AZ_file)
        [b,unique_location,j1] = unique(Affect_location,'rows','first');
        Affzon = Affect_zones(unique_location,:);
        Affloc = Affect_location(unique_location,:);
        load([FT_T1 filesep patstr filesep 'AllRegions.mat'])
        [hh,xx] = hist(Regions(:),1:1:116);
        [hha,xxa] = hist(Affzon(:),1:1:116);
        perc_Affect = hha./hh;
        All_PA = perc_Affect;
        OrigMat = load([FT_T1 filesep 'AllConnMatrices' num2str(atlassize) '_FT' filesep 'Matrix_TFC' num2str(atlassize) '_' patstr '.mat']);
        sAffzon = sort(Affzon,2,'ascend');
        Urows = unique(sAffzon,'rows');
        ConMat = OrigMat.MapTFC;
        for j = 1:size(Urows,1)
            if length(find(Urows(j,:)))==2;
                pc = sum(ismember(sAffzon,Urows(j,:),'rows'));
                ConMat(Urows(j,1),Urows(j,2)) = ConMat(Urows(j,1),Urows(j,2))-pc;
                ConMat(Urows(j,2),Urows(j,1)) = ConMat(Urows(j,1),Urows(j,2));
            end
        end
        nConMat = ConMat./OrigMat.total_fcount;
        LoCoResults(i).NormPatID = patstr;
        LoCoResults(i).Regions = All_PA;
        LoCoResults(i).ConMat = ConMat;
        LoCoResults(i).nConMat = nConMat;
        LoCoResults(i).OrigMat = OrigMat;
    end
    mLoCo = zeros(size(LoCoResults(1).Regions));
    CM = zeros(size(LoCoResults(1).ConMat));
    nCM = zeros(size(LoCoResults(1).nConMat));
    for i = 1:length(LoCoResults)
        mLoCo = mLoCo + LoCoResults(i).Regions;
        CM = CM + LoCoResults(i).ConMat;
        nCM = nCM + LoCoResults(i).nConMat;
    end
    sLoCo = [];
    sLoCo.NormPatID = 'Mean Results';
    sLoCo.Regions = mLoCo./length(LoCoResults);
    sLoCo.ConMat = CM./length(LoCoResults);
    sLoCo.nConMat = nCM./length(LoCoResults);
    sLoCo.OrigMat = [];
    LoCoResults = [sLoCo LoCoResults];
    All_LoCoT1 = [results_dir filesep 'LoCo' num2str(atlassize) StrSave 'T1'];
    save(All_LoCoT1,'LoCoResults')
    time_for_one_mask = toc
end

function [AZ_file,Affect_location] = PathDisruption(WMseed_Weights,ROI_numFileName,OutPathFileName,PathsByRegion)

% WMmask = load(WMmaskFileName); WMf = fieldnames(WMmask); eval(['WMmask = WMmask.' WMf{1} ';']); WMmask = find(WMmask);
% 
% if ischar(WMseedFileName)
%     WMseed = load(WMseedFileName); WMf = fieldnames(WMseed); eval(['WMseed = WMseed.' WMf{1} ';']); WMseed = find(WMseed);
% else
%     WMseed = find(WMseedFileName);
% end
% 
% WMseed = intersect(WMseed,WMmask);
WMseed = find(WMseed_Weights);

Rmap = load(ROI_numFileName);
ROI_num = Rmap.ROI_num;

ROI_num_WMseeds = ROI_num(WMseed);
runindex = setdiff(unique(ROI_num_WMseeds'),0);
runl = length(runindex);

Affect_zones = [];
Affect_location = [];
Affect_weight = [];
Affect_subject = [];
parfor ii = 1:runl
    i = runindex(ii);
    %First load all the paths within the ROI
    PR = load(strcat(PathsByRegion, [filesep 'PathsInROI_'], num2str(i),'.mat'));
    PathsROI = PR.PathsROI;     
    pix_in_ROI = WMseed(ROI_num_WMseeds==i);    
    if ~isempty(pix_in_ROI)
        [WMi,WMj,WMk] = ind2sub(size(ROI_num),pix_in_ROI);        
        for k = 1:length(PathsROI)
            path_floor = floor(PathsROI(k).voxels);
            intvox = intersect(path_floor',[WMi,WMj,WMk],'rows');
            if  ~isempty(intvox)
                Affect_zones = [Affect_zones; PathsROI(k).regions];
                Affect_location = [Affect_location; PathsROI(k).location];
                Affect_weight = [Affect_weight; max(WMseed_Weights(sub2ind(size(ROI_num),intvox(:,1),intvox(:,2),intvox(:,3))))];
               % Affect_subject = [Affect_subject; PathsROI(k).PatientID];
            end
        end
    end
end
save(OutPathFileName,'Affect_zones','Affect_location','WMseed','Affect_weight','Affect_subject')
AZ_file = OutPathFileName;