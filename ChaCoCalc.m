function All_LoCoMNI = ChaCoCalc(DamageFileName,Coreg2MNI,CalcSpace,atlassize,StrSave,NumWorkers,dispMask,coregOnly)
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
% coregOnly             A flag indicating if you'd only like to run the
%                       coregistration, not the ChaCo analysis.
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
if nargin < 8 || isempty(coregOnly)
    coregOnly = 0;
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
SVNEveTools = ['mymfiles' filesep 'eve_tools']; %change into relative path
AtlasingDir = ['resource' filesep 'Atlasing'];
SPMdir8 = ['mymfiles' filesep 'spm8'];
niitools = ['mymfiles' filesep 'nifti_toolbox'];
ODFtools = ['mymfiles' filesep 'ODF_Qball'];
SCBtrack = ['mymfiles' filesep 'Tracking_Bayes'];
mymfiles = 'mymfiles';

if ~(ismcc() || isdeployed())
addpath(['mymfiles' filesep 'BCT']);
addpath(SCBtrack)
addpath(mymfiles)
p = genpath(SVNEveTools);
addpath(p);
p3 = genpath(SPMdir8);
addpath(p3);
addpath(AtlasingDir)
addpath(niitools)
addpath(ODFtools)
end

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

All_LoCoMNI = {};
if coregOnly
    return
end

results_dir = pth;
if ~isdir(results_dir); mkdir(results_dir); end
results_dir2 = [pth filesep 'DetailedChaCoResults'];
if ~isdir(results_dir2); mkdir(results_dir2); end
main_dir = ['Tractograms' filesep];
FT_MNI = [main_dir 'FiberTracts' num2str(atlassize) '_MNI_BIN'];

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


if dispMask
    openMRIcron([TempFileName ' -b 20 -o ' DamageFileName])
end
if strcmp(CalcSpace,'MNI') || strcmp(CalcSpace,'both')
    tic
    LoCoResults = [];
    for i = 1:size(PatStr,1)
        disp(['Calculating ChaCo for patient ' num2str(i) ' of ' num2str(size(PatStr,1))])
        PathsByRegion = deblank(PatStr(i,:));
        patstr = PathsByRegion(end-9:end);
        OutPathFileName = [StrSave filesep 'Dam_Zones'  '_MNI_' patstr '.mat'];
        %Get the paths that go through the damage mask.
        if ~exist(OutPathFileName,'file')
            [AZ_file,~] = PathDisruption(Damage,ROI_num_filename,OutPathFileName,PathsByRegion);
        else
            AZ_file = OutPathFileName;
        end
        load(AZ_file)
        [b,unique_location,j1] = unique(Affect_location,'rows','first');
        Affzon = Affect_zones(unique_location,:);
        Affloc = b;
        
        AfW = NaN(length(Affloc),1);
        for jj = 1:size(Affloc,1);
            AfW(jj) = min(Affect_weight(j1==jj));
        end
        AfW = [AfW AfW];
        
        load([FT_MNI filesep patstr filesep 'AllRegions.mat'])
        [hh,xx] = hist(Regions(find(Regions)),1:1:atlassize);
        
        WAffect = NaN(atlassize,1);
        for k = 1:atlassize;
            WAffect(k) = sum(AfW(find(Affzon==k)));
        end
        perc_WAffect = WAffect'./hh;
        All_PA = perc_WAffect;
        
        
        OrigMat = load([FT_MNI filesep 'AllConnMatrices' num2str(atlassize) '_FT_MNI' filesep 'Matrix_TFC' num2str(atlassize) '_MNI_' patstr(1:end-3) '.mat']);
        sAffzon = sort(Affzon,2,'ascend');
        Urows = unique(sAffzon,'rows');
        ConMat = OrigMat.MapTFC;
        for j = 1:size(Urows,1)
            if length(find(Urows(j,:)))==2;
                pc = sum(AfW(ismember(sAffzon,Urows(j,:),'rows'),1));
                ConMat(Urows(j,1),Urows(j,2)) = ConMat(Urows(j,1),Urows(j,2))+pc;
                ConMat(Urows(j,2),Urows(j,1)) = ConMat(Urows(j,1),Urows(j,2));
            end
        end
        nConMat = ConMat./OrigMat.total_fcount;
        LoCoResults(i).NormPatID = patstr(1:end-3);
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
    for ii = 1:length(LoCoResults)
        mLoCo = mLoCo + LoCoResults(ii).Regions;
        CM = CM + LoCoResults(ii).ConMat;
        nCM = nCM + LoCoResults(ii).nConMat;
        OM = OM + LoCoResults(ii).OrigMat.MapTFC;
        nOM = nOM + LoCoResults(ii).OrigMat.nMapTFC;
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
    [netmet(1,1),netmet(1,2),~,netmet(1,3),netmet(1,4)] = charpath(D);
    sLoCo.nConMatMets = [netmet(1,1),netmet(1,2),netmet(1,3),netmet(1,4)];
    D = distance_wei(1./sLoCo.ConMat);
    [netmet2(1,1),netmet2(1,2),~,netmet2(1,3),netmet2(1,4)] = charpath(D);
    sLoCo.ConMatMets = [netmet2(1,1),netmet2(1,2),netmet2(1,3),netmet2(1,4)];
    D = distance_wei(1./sLoCo.OrigMat.MapTFC);
    [netmet3(1,1),netmet3(1,2),~,netmet3(1,3),netmet3(1,4)] = charpath(D);
    sLoCo.OrigMat.MetsMapTFC = [netmet3(1,1),netmet3(1,2),netmet3(1,3),netmet3(1,4)];
    D = distance_wei(1./sLoCo.OrigMat.nMapTFC);
    [netmet4(1,1),netmet4(1,2),~,netmet4(1,3),netmet4(1,4)] = charpath(D);
    sLoCo.OrigMat.nMetsMapTFC = [netmet4(1,1),netmet4(1,2),netmet4(1,3),netmet4(1,4)];
    ChaCoResults = [sLoCo LoCoResults];
    
    All_LoCoMNI = [StrSave filesep 'ChaCo' num2str(atlassize) '_MNI'];
    save(All_LoCoMNI,'ChaCoResults');
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
for ii = 1:runl
    i = runindex(ii);
    %First load all the paths within the ROI
    PathsROI = nemo_read(strcat(PathsByRegion, [filesep 'PathsInROI_'], num2str(i),'.bin'), [181 217 181]); % To be replaced by nemo_read function
    pix_in_ROI = WMseed(ROI_num_WMseeds==i);    
    if ~isempty(pix_in_ROI)
        [WMi,WMj,WMk] = ind2sub(size(ROI_num),pix_in_ROI);        
        for k = 1:length(PathsROI)
            path_floor = floor(PathsROI(k).voxels);
            intvox = intersect(path_floor',[WMi,WMj,WMk],'rows');
            if  ~isempty(intvox)
                Affect_zones = [Affect_zones; PathsROI(k).regions];
                Affect_location = [Affect_location; PathsROI(k).location];
                Affect_weight = [Affect_weight; min(WMseed_Weights(sub2ind(size(ROI_num),intvox(:,1),intvox(:,2),intvox(:,3))))];
%                 if ~(max(WMseed_Weights(sub2ind(size(ROI_num),intvox(:,1),intvox(:,2),intvox(:,3))))==min(WMseed_Weights(sub2ind(size(ROI_num),intvox(:,1),intvox(:,2),intvox(:,3)))))
%                     disp('shit')
%                 end
                    % Affect_subject = [Affect_subject; PathsROI(k).PatientID];
            end
        end
    end
end
save(OutPathFileName,'Affect_zones','Affect_location','WMseed','Affect_weight','Affect_subject')
AZ_file = OutPathFileName;