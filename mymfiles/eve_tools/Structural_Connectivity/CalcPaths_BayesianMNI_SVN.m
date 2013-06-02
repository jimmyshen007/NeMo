function [time_for_streamlines, time_for_findfibre, time_for_conn]=CalcPaths_BayesianMNI_SVN(PatDir,T1FileName,atlassize,diff_dir,diff_ext,diff_omit,bvecext,overwrite_fiber,redo_ODFonly)
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

%Create array of DTI filenames located in diff_dir
if ~isempty(diff_omit)
    [dirsearch,~]=Dir2Arr(diff_dir,diff_ext,diff_omit);
else
    [dirsearch,~]=Dir2Arr(diff_dir,diff_ext);
end

DiffFileNames = cellstr(dirsearch);
    
GradTextFileName = deblank(Dir2Arr(diff_dir,bvecext));

%The atlas to be used in the parcellation of the GM.
[AtlasFileName,AtlasTextFileName]=ChooseAtlas(atlassize,Atlasing_dir);

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
if (overwrite_fiber == 1 || ~isdir([diff_dir filesep 'FiberTracts' num2str(atlassize)])) && redo_ODFonly == 0    
    %Input the T1FileName (with directory) and produce the normalized T1, as
    %well as the WM, GM, CSF segmentation and parcellated GM atlas. Note:
    %possible to pass cell array of T1files at once

    %Rename the segmentation files so that you choose the ones in MNI space
    [gms1p,gms1fn,gms1ext] = fileparts(T1FileName);
    GMSegFileName = [gms1p filesep 'Segmented' filesep 'wc1' gms1fn gms1ext];
    WMSegFileName = [gms1p filesep 'Segmented' filesep  'wc2' gms1fn gms1ext];
    CSFSegFileName = [gms1p filesep 'Segmented' filesep 'wc3' gms1fn gms1ext];
    [GMMaskFileName, WMMaskFileName, ContMaskFileName] = makeMasks_prob_SVN(GMWMcontrast, GMSegFileName, WMSegFileName, CSFSegFileName);
    disp('Make regions');
    GMvol = spm_read_vols(spm_vol(GMMaskFileName));
    Atvol = spm_read_vols(spm_vol(AtlasFileName));
    wAtlas = zeros(size(GMvol));
    wAtlas(find(GMvol)) = Atvol(find(GMvol));
    ATinfo = spm_vol([gms1p filesep 'Atlased' num2str(atlassize) filesep gms1fn '_Atlas' gms1ext]);
    ATinfo.fname = [gms1p filesep 'Atlased' num2str(atlassize) filesep 'w' gms1fn '_Atlas' gms1ext];
    ATinfo.private.dat.fname = ATinfo.fname;
    AtlasedFileName = ATinfo.fname;
    spm_create_vol(ATinfo)
    spm_write_vol(ATinfo,wAtlas)
    
    %% Make Regions
    %[Zones_OutDir, nn, xx] = fileparts(rsGMMaskFileName);
    [Zones_OutDir, ~, ~] = fileparts(GMMaskFileName);
    %[pp,nn,xx] = fileparts(AtlasTextFileName);
    
    %[ZonesFileName, MaskSurfFileName, MaskTrackingFileName, Zones] = makeZones(NumZones,AtlasTextFileName,rsAtlasedFileName,rsContMaskFileName,rsWMMaskFileName,Zones_OutDir);
    [ZonesFileName, ~, MaskTrackingFileName, Zones] = makeZones_SVN(NumZones,AtlasTextFileName,AtlasedFileName,ContMaskFileName,WMMaskFileName,Zones_OutDir);
    %% Compute_QBI
    [ODFFileName, ODFdirsTextFileName] = compute_QBI_AR_SVN(ContMaskFileName, GradTextFileName, DiffFileNames);
    %[~, ODFdirsTextFileName] = compute_QBI_AR(rsContMaskFileName, GradTextFileName, rsDiffFileNames);
    seed_begin = 1;
else
    ZonesFileName = [t1pth filesep 'Segmented' filesep 'Zones.mat'];
    MaskTrackingFileName = [t1pth filesep 'Segmented' filesep 'Mask_forTracking' t1ext];
    AtlasedFileName = [t1pth filesep 'Atlased' num2str(atlassize) filesep t1fn '_Atlas' t1ext];  
    ftfiles = dir([diff_dir filesep 'FiberTracts' num2str(atlassize) filesep 'Paths' num2str(atlassize)]);
    numfiles = size(ftfiles,1)-2;
    seed_begin = numfiles*NumXfiles + 1;
    ODFdirsTextFileName = [diff_dir filesep 'ODF_directions.txt'];
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
    
    diff_odf=Dir2Arr(diff_dir,['ODF_war*' t1ext]);
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
    Atlas_name = AtlasedFileName;
 
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
        trk_amass([1 atlassize],[OutPath filesep qPpaths{4}],deblank(DiffFileNames{1}),[OutPath filesep 'all_brain' num2str(atlassize) '.trk'],atlassize,1);
%     end
       
    
    return;
end
