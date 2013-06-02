function [time_for_streamlines, time_for_findfibre, time_for_conn]=CalcPaths_Bayesian_SVN(PatDir,T1FileName,atlassize,diff_dir,diff_ext,diff_omit,bvecext,preatlasedopt,affinereg,overwrite_fiber,redo_ODFonly,FreeSurfer)
%A port of the standalone CalcPaths_Bayesian.m

   
%Load system variables for function access
startup_varsonly;

%This is the directory of the toolbox.
%SCBdir =  SCB_dir;
SPMdir = SPM_dir;

%Initialize the global variables used in this function
global NumZones Zones;
NumZones = atlassize;
global obj tracking seed Tzones step longMax NumXfiles repetitions ang numVectors Diff_name inf_name Mask_name Atlas_name;

%Default parameter settings if not part of the input.
if nargin < 12; 
    FreeSurfer = 0;
end
if nargin < 11; 
    redo_ODFonly = 0;
end
if nargin < 10;
    overwrite_fiber = 1;
end
if nargin < 9;
    affinereg = 1;
end

%Create array of DTI filenames located in diff_dir
if ~isempty(diff_omit)
    [dirsearch,~]=Dir2Arr(diff_dir,diff_ext,diff_omit);
else
    [dirsearch,~]=Dir2Arr(diff_dir,diff_ext);
end

DiffFileNames = cellstr(dirsearch);
    
GradTextFileName = deblank(Dir2Arr(diff_dir,bvecext));

%The atlas to be used in the parcellation of the GM.
if ~FreeSurfer
    [AtlasFileName,AtlasTextFileName]=ChooseAtlas(atlassize,Atlasing_dir);
else
    AtlasTextFileName = '/home/amy/work/BTFdata/Normals/FreeSurfer_NewROIs.txt';
end

qPpaths={['Paths' num2str(atlassize)],['Poscs' num2str(atlassize)],['Pbelong' num2str(atlassize)],['qbelong' num2str(atlassize)]};

GMWMcontrast = 1;

%Leave this variable empty so that the subdirectories and output files are 
%created in the patient's T1 directory. 
[t1pth,t1fn,t1ext]=fileparts(T1FileName);
StrOutputDir ='';

%This is the treshold for the GM segmentation (usually leave empty).
Thresh = [];

%This is the type of structural image, normally T1, but could also be T2,
%etc.
Type = 't1';
[TempFileName] = Choosing_TempFile(Type);

step = 0.7;
longMax = 100;
NumXfiles = 100;
repetitions = 50;
ang = pi/3;
numVectors = [];

%If all of the atlas and segmented files exist, and some fiber tracking has
%already been done, skip directly to the fiber tracking and pick up where
%you left off.
if (overwrite_fiber == 1 || ~isdir([diff_dir filesep 'FiberTracts' num2str(atlassize)])) && redo_ODFonly == 0 && FreeSurfer == 0    
    %Input the T1FileName (with directory) and produce the normalized T1, as
    %well as the WM, GM, CSF segmentation and parcellated GM atlas. Note:
    %possible to pass cell array of T1files at once
    
    %Generate the WM, GM, CSF segmentationsnormalized T1 image.
    if ~preatlasedopt
        disp('Approximate Affine co-registration...')
        T1FileName = approx_coreg2MNI(TempFileName,T1FileName);
        [GMSegFileName, WMSegFileName, CSFSegFileName, AtlasedFileName] = Atlasing_spm8_SVN(T1FileName, StrOutputDir, AtlasFileName, TempFileName, Thresh,SPMdir,1,atlassize,affinereg);
        % [GMSegFileName, WMSegFileName, CSFSegFileName, AtlasedFileName] = Atlasing(T1FileName, StrOutputDir, AtlasFileName, TempFileName, Thresh);
    else
        [GMSegFileName, WMSegFileName, CSFSegFileName, AtlasedFileName]=AutoLabel_Only_SVN(T1FileName,'',atlassize);
    end
    
    %Make Masks
    disp('Skull strip the T1 and impose brain mask')
    nskT1FileName = deblank(Dir2Arr(t1pth,['nsk*' t1ext]));
    if isempty(nskT1FileName)
        [nskT1FileName] = SkullStrip_SVN(T1FileName, GMSegFileName, WMSegFileName, CSFSegFileName,1);
    end
    GMMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['Mask_Gray' t1ext]));
    if isempty(GMMaskFileName)
        disp('Make Masks');
        [GMMaskFileName, WMMaskFileName, ContMaskFileName] = makeMasks_prob_SVN(GMWMcontrast, GMSegFileName, WMSegFileName, CSFSegFileName);
    else
        WMMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['Mask_White' t1ext]));
        ContMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['Mask_contrast' t1ext]));
    end
    rsGMMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['sMask_Gray' t1ext]));
    [pp, nn, xx] = fileparts(AtlasedFileName);
    rsAtlasedFileName = fullfile(pp, ['s' nn xx]);
    if isempty(rsGMMaskFileName)
        disp('Resample 2dWMRI');
        [rsDiffFileNames, rsGMMaskFileName, rsWMMaskFileName, rsContMaskFileName, rsAtlasedFileName] = resampling2DWMRI_SVN2(DiffFileNames, GMMaskFileName, WMMaskFileName, ContMaskFileName, AtlasedFileName,nskT1FileName,1);
    elseif ~exist(rsAtlasedFileName,'file')
        copyfile(AtlasedFileName,rsAtlasedFileName)
        if strcmp(xx,'.img')
            copyfile([pp filesep nn '.hdr'],[pp filesep 's' nn '.hdr'])
        end
        load([t1pth filesep 'LinearTrans_T12DTI'])
        spm_get_space(rsAtlasedFileName,M*spm_get_space(rsAtlasedFileName))
        rflags.mask = 0; rflags.mean = 0; rflags.interp = 0; rflags.which = 1; rflags.prefix = 's';
        rsDiffFileNames = DiffFileNames;
        spm_reslice({rsDiffFileNames{1},rsAtlasedFileName},rflags)                
        [pp, nn, xx] = fileparts(AtlasedFileName);
        movefile([pp filesep 'ss' nn xx],[pp filesep 's' nn xx])
        if strcmp(xx,'.img')
            movefile([pp filesep 'ss' nn '.hdr'],[pp filesep 's' nn '.hdr'])
        end
        rsContMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['sMask_contrast' t1ext]));
        rsWMMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['sMask_White' t1ext]));
        openMRIcron([rsDiffFileNames{1} ' -b 20 -o ' rsAtlasedFileName ' -o ' rsWMMaskFileName]);
    else
        rsDiffFileNames = DiffFileNames;
        rsContMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['sMask_contrast' t1ext]));
        rsWMMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['sMask_White' t1ext]));
        openMRIcron([rsDiffFileNames{1} ' -b 20 -o ' rsAtlasedFileName ' -o ' rsWMMaskFileName]);
    end
    
    disp('Make regions');
    %% Make Regions
    %[Zones_OutDir, nn, xx] = fileparts(rsGMMaskFileName);
    [Zones_OutDir, ~, ~] = fileparts(rsGMMaskFileName);
    %[pp,nn,xx] = fileparts(AtlasTextFileName);
    
    %[ZonesFileName, MaskSurfFileName, MaskTrackingFileName, Zones] = makeZones(NumZones,AtlasTextFileName,rsAtlasedFileName,rsContMaskFileName,rsWMMaskFileName,Zones_OutDir);
    [ZonesFileName, ~, MaskTrackingFileName, Zones] = makeZones_SVN(NumZones,AtlasTextFileName,rsAtlasedFileName,rsContMaskFileName,rsWMMaskFileName,Zones_OutDir);
    %% Compute_QBI
    ODFdirsTextFileName = deblank(Dir2Arr(diff_dir,'ODF_directions.txt'));
    if isempty(ODFdirsTextFileName)
        [~, ODFdirsTextFileName] = compute_QBI_AR_SVN(rsContMaskFileName, GradTextFileName, rsDiffFileNames);
    else
        disp('Orientation distrubution function (ODF) already created, continuing...')
    end
    %[~, ODFdirsTextFileName] = compute_QBI_AR(rsContMaskFileName, GradTextFileName, rsDiffFileNames);
    seed_begin = 1;
elseif redo_ODFonly
    [~,f,ext] = fileparts(DiffFileNames{1});
    minisovoxeldim=2.5; % will be resampled to this vox size if larger than this
    if ~(strcmp(f(1:3),'iso') || strcmp(f(1:2),'rf'))
        for ii = 1:length(DiffFileNames)
            [DiffFilePt, DiffFileNm, ext] = fileparts(DiffFileNames{ii});
            rsDiffFileNames{ii} = fullfile(DiffFilePt, ['iso-' DiffFileNm ext]);
            tmphdr = load_nii_hdr(DiffFileNames{ii});
            %reslice_nii(DiffFileNames{i}, newDiffFileNames{i}, (1/50)*round(50*max(tmphdr.dime.pixdim(2:4)))*ones(1,3));
            reslice_nii(DiffFileNames{ii}, rsDiffFileNames{ii}, (1/50)*round(50*min(max(tmphdr.dime.pixdim(2:4)),minisovoxeldim))*ones(1,3));
        end
    else
        rsDiffFileNames = DiffFileNames;
    end
    rsContMaskFileName = [t1pth filesep 'Segmented' filesep 'sMask_contrast' num2str(atlassize) ext];
    ODFdirsTextFileName = deblank(Dir2Arr(diff_dir,'ODF_directions.txt'));
    if isempty(ODFdirsTextFileName)
        [~, ODFdirsTextFileName] = compute_QBI_AR_SVN(rsContMaskFileName, GradTextFileName, rsDiffFileNames);  
    else
        disp('Orientation distrubution function (ODF) already created, continuing...')
    end
    ZonesFileName = [t1pth filesep 'Segmented' filesep 'Zones' num2str(atlassize) '.mat'];
    MaskTrackingFileName = [t1pth filesep 'Segmented' filesep 'Mask_forTracking' num2str(atlassize) t1ext];
    rsAtlasedFileName = [t1pth filesep 'Atlased' num2str(atlassize) filesep 's' t1fn '_Atlas' t1ext];  
    seed_begin = 1;
elseif FreeSurfer
    disp('Skull strip the T1 and impose brain mask')
    nskT1FileName = deblank(Dir2Arr(t1pth,'abrainmask.nii'));
    GMMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['Mask_Gray' t1ext]));
    WMMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['Mask_White' t1ext]));
    ContMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['Mask_contrast' t1ext]));
    AtlasedFileName = deblank(Dir2Arr([t1pth filesep 'Atlased' num2str(atlassize)],[t1fn '_atlas' num2str(atlassize) t1ext]));
    rsGMMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['sMask_Gray' t1ext]));
    [pp, nn, xx] = fileparts(AtlasedFileName);
    rsAtlasedFileName = fullfile(pp, ['s' nn xx]);
    if isempty(rsGMMaskFileName)
        disp('Resample 2dWMRI');
        [rsDiffFileNames, rsGMMaskFileName, rsWMMaskFileName, rsContMaskFileName, rsAtlasedFileName] = resampling2DWMRI_SVN2(DiffFileNames, GMMaskFileName, WMMaskFileName, ContMaskFileName, AtlasedFileName,nskT1FileName,1);
    else
        rsDiffFileNames = DiffFileNames;
        rsContMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['sMask_contrast' t1ext]));
        rsWMMaskFileName = deblank(Dir2Arr([t1pth filesep 'Segmented'],['sMask_White' t1ext]));
        openMRIcron([rsDiffFileNames{1} ' -b 20 -o ' rsAtlasedFileName ' -o ' rsWMMaskFileName]);
    end
    
    disp('Make regions');
    %% Make Regions
    %[Zones_OutDir, nn, xx] = fileparts(rsGMMaskFileName);
    [Zones_OutDir, ~, ~] = fileparts(rsGMMaskFileName);
    %[pp,nn,xx] = fileparts(AtlasTextFileName);
    ZonesFileName = [t1pth filesep 'Segmented' filesep 'Zones' num2str(atlassize) '.mat'];
    if ~exist(ZonesFileName,'file')
        [ZonesFileName, ~, MaskTrackingFileName, Zones] = makeZones_SVN(NumZones,AtlasTextFileName,rsAtlasedFileName,rsContMaskFileName,rsWMMaskFileName,Zones_OutDir);
    else
        load(ZonesFileName)
       MaskTrackingFileName = [t1pth filesep 'Segmented' filesep 'Mask_forTracking' num2str(atlassize) t1ext];
    end
    %% Compute_QBI
    ODFdirsTextFileName = deblank(Dir2Arr(diff_dir,'ODF_directions.txt'));
    if isempty(ODFdirsTextFileName)
        [~, ODFdirsTextFileName] = compute_QBI_AR_SVN(rsContMaskFileName, GradTextFileName, rsDiffFileNames);
    else
        disp('Orientation distrubution function (ODF) already created, continuing...')
    end
    %[~, ODFdirsTextFileName] = compute_QBI_AR(rsContMaskFileName, GradTextFileName, rsDiffFileNames);
    seed_begin = 1;
else
    ZonesFileName = [t1pth filesep 'Segmented' filesep 'Zones' num2str(atlassize) '.mat'];
    MaskTrackingFileName = [t1pth filesep 'Segmented' filesep 'Mask_forTracking' num2str(atlassize) t1ext];
    rsAtlasedFileName = [t1pth filesep 'Atlased' num2str(atlassize) filesep 's' t1fn '_Atlas' t1ext];  
    ftfiles = dir([diff_dir filesep 'FiberTracts' num2str(atlassize) filesep 'Paths' num2str(atlassize)]);
    numfiles = size(ftfiles,1)-2;
    seed_begin = numfiles*NumXfiles + 1;
    ODFdirsTextFileName = [diff_dir filesep 'ODF_directions.txt'];
    rsDiffFileNames = cellstr(Dir2Arr(diff_dir,['iso-*' t1ext]));
end
    %eval(['load ''' ZonesFileName '''']);
    load(ZonesFileName);
    
    %Get numdir from ODF direction file lines
    odffid=fopen([diff_dir filesep 'ODF_directions.txt']);
    tline=fgets(odffid); numdir=0;
    while ischar(tline)
        numdir=numdir+1;
        tline=fgets(odffid);
    end
    fclose(odffid);
    
    Diff_name = {};
    
    diff_odf=Dir2Arr(diff_dir,'ODF_iso-*',['*.hdr'; '*.mat']);
    %diff_odf=ODFFileName;
    DF=size(diff_odf);
    for jj = 1:DF(1)
        %[dfd,q,dfext]=fileparts(deblank(diff_odf(jj,:)));
        [~,q,dfext]=fileparts(deblank(diff_odf(jj,:)));
        kk = 0;
        for j = 1:numdir
            kk = kk +1;
            if kk < 10
                Diff_name{kk} = [diff_dir filesep q dfext ', ' num2str(kk)];
            else
                Diff_name{kk} = [diff_dir filesep q dfext ',' num2str(kk)];
            end
        end
    end

    Diff_name = char(Diff_name);
    inf_name = ODFdirsTextFileName;
    Mask_name = MaskTrackingFileName;
    Atlas_name = rsAtlasedFileName;
 
    %Get the seed points for the fiber tracking algorithm from the external
    %points to the gray matter regions (contained in Zones.roi).  Also, save
    %the regions that the Zones.roi pixels correspond to as input to the
    %tracking algorithm so that it is not confused about the origin of the
    %tract.
    Tzones = 1:NumZones;
    obj = 2;
    tracking = 'BayesianF';
    seed = [];
    seed_zone = [];
    for zone = Tzones
        seed = [seed; Zones(zone).roi];
        seed_zone = [seed_zone; zone*ones(length(Zones(zone).roi),1)];
    end
    
   
    %Create and save the fiber tracts.
    OutPath = [diff_dir filesep 'FiberTracts' num2str(atlassize)];
    if ~isdir(OutPath)
        mkdir(OutPath)
    end
    if ~isdir([OutPath filesep qPpaths{1}])
        mkdir([OutPath filesep qPpaths{1}])
    end
    if ~isdir([OutPath filesep qPpaths{2}])
        mkdir([OutPath filesep qPpaths{2}])
    end
    tic
    amy_Tracking_SVN(tracking,seed,seed_zone,Tzones,Zones,Diff_name,inf_name,Mask_name,Atlas_name,step,longMax,NumXfiles,repetitions,ang,numVectors,OutPath,qPpaths(1:2),seed_begin);
    time_for_streamlines = toc;
    %Summarize the fiber tracts and put them into files for each pair of gray
    %matter regions.
    tic
    findfibre_SVN(Zones, Tzones, ceil((length(seed_zone)/NumXfiles)),OutPath,[OutPath filesep qPpaths{2}],qPpaths);
    time_for_findfibre = toc;
    
    tic
    parallel_connectivity_saveConn_SVN(Tzones,Zones,OutPath,qPpaths{4},atlassize);
    time_for_conn = toc;

    disp(['Time for streamlines (amy_Tracking_generic): ' int2str(time_for_streamlines/3600.) ' hours']);
    disp(['Time for fiber tracking (findfibre_generic): ' int2str(time_for_findfibre/3600.) ' hours']);
    disp(['Time for connectivity (parallel_connectivity_saveConn_generic): ' int2str(time_for_conn/60.) ' minutes']);

    save([OutPath filesep 'CPB_completion_time_' atlassize '.mat'],'time_for_streamlines','time_for_findfibre','time_for_conn');
    
    %Compute average path length between regions
    avgmatrix=qbelong_trklength_avg_SVN([OutPath filesep qPpaths{4}],atlassize);
    
%     if isempty(find(any(avgmatrix)==0))
%         %Generate TrackVis file of entire brain network for visualization
        trk_amass([1 atlassize],[OutPath filesep qPpaths{4}],deblank(rsDiffFileNames{1}),[OutPath filesep 'all_brain' num2str(atlassize) '.trk'],atlassize,1);
%     end
       
    
    return;
end
