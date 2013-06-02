function [mbeg,mend,numf,time_to_convert] = NonLinearLookupTracts2MNI_SVN(qbelongoriginal_dir,sT1scan,sT1scan_Atlas,wT1scan,atlassize,MRIOPT,LookupListFileName)
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
    openMRIcron([wT1scan ' -b 20 -o ' st1_dir filesep flags.prefix st1_fn st1_ext])
end
pause(0.5);
MNI_atfn = [at_dir filesep flags.prefix at_fn at_ext];
if MRIOPT
    openMRIcron([wT1scan ' -b 20 -o ' MNI_atfn])
end

MNI_at = spm_read_vols(spm_vol(MNI_atfn));
MNI_at(isnan(MNI_at)) = 0;
mni_vox_sz = size(MNI_at);

path_dir = qbelongoriginal_dir;
[pd_dir,pd_file,~] = fileparts(path_dir);
ft = dir(path_dir);

qbelongnew_dir = [pd_dir filesep pd_file '_NLMNI_LUT']; 
if ~isdir(qbelongnew_dir)
    mkdir(qbelongnew_dir)
end

mbeg = 0;
mend = 0;
numf = 0;
mbeg2 = 0;
mend2 = 0;
%Load the transformed index image so that we can lookup voxel by voxel the 
%corresponding voxel in MNI space. 
if isempty(LookupListFileName)
    CreateLookupList(sT1scan,transmat,flags);
    LookupListFileName = deblank(Dir2Arr(st1_dir,'LookUp_MNIvox_1mm.mat'));
end
load(LookupListFileName)

% Calculate the deformation field corresponding to the DTI to MNI
% transformation matrix.
Def = spm_calc_def_SVN(transmat,flags);

tic
for k = 3:length(ft);
    if mod(k,100)==0; disp(['Loading file -> ' num2str(k)]); end
    %Load each of the qbelong files
    if ~isdir([path_dir filesep ft(k).name])
    Q = load([path_dir filesep ft(k).name]);
    q = Q.qbelong;
    qnew = struct('voxels',[],'PT',[]);
    posc_orig =  NaN(size(q,2),2);
    posc_new = NaN(size(q,2),2);    
    if ~(size(q,2)==1)
        for j = 1:size(q,2)
            if ~isempty(q(j).voxels)
                % Find the voxel indices where the fiber traverses, and
                % find the corresponding MNI space voxel.
                p = sub2ind(size(subj_at),floor(q(1).voxels(1,:))',floor(q(1).voxels(2,:))',floor(q(1).voxels(3,:))');
                % and then apply the transformation to MNI space. 
                [txyz,tval] = lookup_coord(p,DTI_vox,MNI_vox,mni_vox_sz);
                t1 = [tval txyz];
                
                qnew(j).PT = q(j).PT;
                qnew(j).voxels = txyz;
                qnew(j).bme = tval;
                posc_orig(j,:) = [subj_at(floor(q(j).voxels(1,1)),floor(q(j).voxels(2,1)),floor(q(j).voxels(3,1))) ...
                    subj_at(floor(q(j).voxels(1,end)),floor(q(j).voxels(2,end)),floor(q(j).voxels(3,end)))];
                tind =sub2ind(size(MNI_at),txyz(:,1),txyz(:,2),txyz(:,3));
                if isempty(find(tval==1,1));
                    posc_new(j,:) = [unique(MNI_at(tind(tval==3))) unique(MNI_at(tind(tval==3)))];
                else
                    posc_new(j,:) = [unique(MNI_at(tind(tval==1))) unique(MNI_at(tind(tval==3)))];
                end
                
                %Compare to the original way
                mask = zeros(size(subj_at));
                p1 = floor(q(j).voxels(:,1));
                p2 = unique(floor(q(j).voxels(:,2:end-1))','rows');
                p3 = floor(q(j).voxels(:,end));
                mask(sub2ind(size(subj_at),p2(:,1),p2(:,2),p2(:,3))) = 2;
                mask(p1(1),p1(2),p1(3)) = 1;
                mask(p3(1),p3(2),p3(3)) = 3;
                % and then apply the transformation to MNI space. 
                [txyz2,tval2] = nonlin_transform_coord(mask,Def);
                t2 = [tval2 txyz2];
                
               % qnew2(j).PT = q(j).PT;
               % qnew2(j).voxels = txyz;
                tind2 =sub2ind(size(Def.X),txyz2(:,1),txyz2(:,2),txyz2(:,3));
                if isempty(find(tval2==1,1));
                    posc_new2(j,:) = [unique(MNI_at(tind2(tval2==3))) unique(MNI_at(tind2(tval2==3)))];
                else
                    posc_new2(j,:) = [unique(MNI_at(tind2(tval2==1))) unique(MNI_at(tind2(tval2==3)))];
                end
            if ~(size(intersect(t1,t2,'rows'),1)==size(t1,1))
                display('Did not match...FUCK!')
            end
                
            end
        end
        % Check that the fibers start and end in the same regions as they
        % did before the transformation.
        match_orig = ones(size(q,2),2);
        match_orig(find(posc_orig(:,1)-posc_new(:,1)),1) = 0;
        match_orig(find(posc_orig(:,2)-posc_new(:,2)),2) = 0;
        mbeg = mbeg + length(find(match_orig(:,1)));
        mend = mend +  length(find(match_orig(:,2)));
        
        match_orig2 = ones(size(q,2),2);
        match_orig2(find(posc_orig(:,1)-posc_new2(:,1)),1) = 0;
        match_orig2(find(posc_orig(:,2)-posc_new2(:,2)),2) = 0;
        mbeg2 = mbeg2 + length(find(match_orig2(:,1)));
        mend2 = mend2 +  length(find(match_orig2(:,2)));
        numf = numf + size(q,2);
    else
        qnew = q;
    end
    qbelong = qnew;
    %parsave([qbelongnew_dir filesep ft(k).name],qbelong)
    end
end
time_to_convert = toc;



function [txyz,tval] = lookup_coord(p,DTI_vox,MNI_vox,mni_vox_sz)
% Input Parameters:
% C         : The mask that contains the voxels with the fibers going
%           through it, with a 1 in the beginning voxel and a 3 in the end voxel.
% Def       : A 1 x 3 structure, with fields X, Y and Z that give the
%           deformation field
%
% Output Parameter:
% txyz      : A 3 x N matrix with the x,y,z coordinates (in voxels) 
%           of the corresponding points in the MNI template image     
% tval      : A 3 x 1 vector with values 1,2, and 3 denoting the beginning
%           voxel, the middle voxels and the end voxel, respectively.
%--------------------------------------------------------------------------
tind2 = [];
mm1 = MNI_vox(DTI_vox == p(1),:)';
tind1 = mm1(~isnan(mm1));
tval1 = ones(size(tind1));
pmid = unique(p(find(p~=p(1),1,'first'):find(p~=p(end),1,'last')));
for i = 1:length(pmid);
    mm = MNI_vox(DTI_vox == pmid(i),:)';
    tind2 = [tind2; mm(~isnan(mm))];
end
tval2 = 2*ones(size(tind2));
mm3 = MNI_vox(DTI_vox == p(end),:)';
tind3 = mm3(~isnan(mm3));
tval3 = 3*ones(size(tind3));

[txyz(:,1),txyz(:,2),txyz(:,3)] = ind2sub(mni_vox_sz,[tind1; tind2; tind3]);
tval = [tval1; tval2; tval3];
%VO = spm_write_vol(VO,Dat);
return;

function CreateLookupList(sT1scan,transmat,flags)
disp('Creating voxel-by-voxel lookup table from diffusion to MNI space')
[st1_dir,~,st1_ext] = fileparts(sT1scan);
vinfo = spm_vol(sT1scan);
sT1 = spm_read_vols(vinfo);
sIndex = zeros(vinfo.dim);
DTI_vox = find(sT1);
sIndex(DTI_vox) = 1:1:length(DTI_vox);
vinfo.fname = [st1_dir filesep 'sT1_IndexImage' st1_ext];
vinfo.private.dat.fname = vinfo.fname;
vinfo.dt = [8 0];
spm_create_vol(vinfo);
spm_write_vol(vinfo,sIndex);
spm_write_sn(vinfo,transmat,flags);
wIndex = spm_read_vols(spm_vol([st1_dir filesep flags.prefix 'sT1_IndexImage' st1_ext]));
n = hist(wIndex(find(wIndex)),length(unique(wIndex))-1);
MNI_vox = NaN(length(DTI_vox),max(n));
tic
for i = 1:length(DTI_vox)
    if mod(i,10000)==0; disp(['Voxel ' num2str(i) ' of ' num2str(length(DTI_vox))]); end
    MNI_vox(i,1:length(find(wIndex==sIndex(DTI_vox(i))))) = find(wIndex==sIndex(DTI_vox(i)));
end
save([st1_dir filesep 'LookUp_MNIvox_1mm.mat'],'DTI_vox','MNI_vox');
tt = toc;
disp(['Lookup Table Creation time = ' num2str(tt/3600) ' hours'])
return;

function parsave(filename,qbelong)

save(filename,'qbelong')
return;

function [txyz,tval] = nonlin_transform_coord(C,Def)
% Input Parameters:
% C         : The mask that contains the voxels with the fibers going
%           through it, with a 1 in the beginning voxel and a 3 in the end voxel.
% Def       : A 1 x 3 structure, with fields X, Y and Z that give the
%           deformation field
%
% Output Parameter:
% txyz      : A 3 x N matrix with the x,y,z coordinates (in voxels) 
%           of the corresponding points in the MNI template image     
% tval      : A 3 x 1 vector with values 1,2, and 3 denoting the beginning
%           voxel, the middle voxels and the end voxel, respectively.
%--------------------------------------------------------------------------
d  = [0 0 0; 0 0 0];

% Accumulate data
%Dat= zeros(VO.dim(1:3));
[xx,yy,zz] = ind2sub(size(C),find(C));
Dat = zeros(size(Def.X));
slices = max(1,floor(0.5*min(zz))):1:min(floor(3*max(zz)),size(Def.X,3));
for i = slices;
    dat = spm_bsplins(C,Def.X(:,:,i),Def.Y(:,:,i),Def.Z(:,:,i),d);
    dat(isnan(dat)) = 0;
    Dat(:,:,i) = dat;
end
%Dat2 = spm_bsplins(C,Def.X,Def.Y,Def.Z,d);
%Dat2(isnan(Dat2)) = 0;
[txyz(:,1),txyz(:,2),txyz(:,3)] = ind2sub(size(Dat),find(Dat));
tval = Dat(find(Dat));
%VO = spm_write_vol(VO,Dat);
return;