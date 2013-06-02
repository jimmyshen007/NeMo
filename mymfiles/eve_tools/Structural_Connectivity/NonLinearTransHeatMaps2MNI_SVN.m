function [time_to_convert] = NonLinearTransHeatMaps2MNI_SVN(qbelongoriginal_dir,sT1scan,sT1scan_Atlas,wT1scan,ZonesFileName,MRIOPT)
%% Inputs: type, description
% qbelongoriginal_dir: (char) directory that holds the diffusion space 
%                      tacts (qbelong) 
% sT1scan            : (char) The full path to the T1 scan that is in
%                      diffusion space (sacos*.img) or the skull-stripped  
%                      version of the diffusion space T1 (snskacos*.img)
% sT1scan_Atlas      : (char) The full path to the atlased T1 image in
%                      diffusion space (sacos*_Atlas.img).
% wT1scan            : (char) The full path to the normalized T1 scan in
%                      MNI space (wacos*.img), or the skull-stripped 
%                      version of it (nskwacos*.img). 
% atlassize          : (int) The number of regions in the atlas used for 
%                      fiber tracking
% MRIOPT             : (binary) Flag for plotting the normalized T1 and
%                       atlased T1 in MNI space (0 = no, 1 = yes)
% Outputs: type, description
% mbeg               : (int) The number of fibers that had the same
%                       starting region before and after transformation.
% mend               : (int) The number of fibers that had the same
%                       ending region before and after transformation.
% numf               : (int) The number of fibers that were transformed.
% time_to_convert    : (double) Number of seconds for the transformation.
%--------------------------------------------------------------------------
%Description:
% This function takes in the fibers produced by amy_tracking_SVN and
% organized into files that denote region pairs by findfibre_SVN (in the 
% qbelong subdirectory of the FiberTracts directory) and normalizes them 
% into MNI space. It first computes the non-linear transformation from 
% diffusion space to MNI space and applies the transformation to each fiber.
% NOTE: If your transformation occurred without error, mbeg, mend and numf 
% should all be equal.
%
% IMPORTANT: wT1scan and sT1scan must either BOTH be skull-stripped or not-skull
% stripped, or your normalization will fail!

%Author: Amy Kuceyeski
%        Weill Cornell Medical College
%        5/1/2012

%%

%First, create the heatmaps from the original qbelong directory, if not
%already created.
heatmap_dir = [qbelongoriginal_dir filesep 'HeatMaps'];
if ~isdir(heatmap_dir) || (length(dir(heatmap_dir))==2)
    Create_HeatMaps(qbelongoriginal_dir,sT1scan);
end

% % Get the roi boundary information.
% load(ZonesFileName);
% for i = 1:length(Zones);
%     skin = Zones(i).roi;
%     for j = 1:length(skin,1);
%         mask(skin(j,1),skin(j,2),skin(j,3)) = 1;
%         spm_write_sn()
% end

% Identify the flags for the non-linear transformation into MNI space for
%the atlas and the fiber masks.
flags = struct('interp',0,'vox',[1 1 1],'bb',NaN,'wrap',[0 0 0],'preserve',0,...
                   'prefix','w');
               
%First, load in the atlas that is in DTI space.
[at_dir,at_fn,at_ext] = fileparts(sT1scan_Atlas);
atinfo = spm_vol(sT1scan_Atlas);
subj_at = spm_read_vols(atinfo);

%Load the DTI space T1 scan.
[st1_dir,st1_fn,st1_ext] = fileparts(sT1scan);

%Normalise the DTI space T1 scan to the same patient's MNI space T1 scan, 
%and apply that transformation to the DTI space atlas and T1 scan. 
transmat = [st1_dir filesep st1_fn '_sn.mat'];
if ~exist(transmat)
    spm_normalise(wT1scan,sT1scan,transmat,[],[],flags)
    spm_write_sn(sT1scan,transmat,flags)
    spm_write_sn(atinfo,transmat,flags)
end

%Check the normalisation to MNI space was ok for both the T1 and atlased
%T1.
if MRIOPT
    openMRIcron([wT1scan ' -b 20 -o ' st1_dir filesep 'w' st1_fn st1_ext])
end

MNI_atfn = [at_dir filesep 'w' at_fn at_ext];
if MRIOPT
    openMRIcron([wT1scan ' -b 20 -o ' MNI_atfn])
end

path_dir = qbelongoriginal_dir;
[pd_dir,pd_file,~] = fileparts(path_dir);

qbelongnew_dir = [pd_dir filesep pd_file '_NLMNI']; 
if ~isdir(qbelongnew_dir)
    mkdir(qbelongnew_dir)
end

heatmapnew_dir = [qbelongnew_dir filesep 'HeatMaps'];
if ~isdir(heatmapnew_dir)
    mkdir(heatmapnew_dir)
end

ft = dir(heatmap_dir);

tic
parfor k = 3:length(ft);
    hm = ft(k).name;    
    if mod(k,100)==0; disp(['Loading file -> ' num2str(k)]); end
    spm_write_sn([heatmap_dir filesep hm],transmat,flags);
end
time_to_convert = toc;