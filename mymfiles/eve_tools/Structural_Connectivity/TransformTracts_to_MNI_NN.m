function [numT,mbeg,mend,timeforconversion,time_for_findfibre,time_for_conn] = TransformTracts_to_MNI_NN(LinearTrans_DTI2T1,NonLinearTrans_seg_inv_sn_T12MNI,Source_AtlasFileName,Template_AtlasFileName,OutPath,Zones,Tzones,in_qPpaths,out_qPpaths,atlassize,ROI_numFileName,redo_MNItrans)

Source_Atlas = spm_read_vols(spm_vol(Source_AtlasFileName));
Template_Atlas = spm_read_vols(spm_vol(Template_AtlasFileName));
files = dir([OutPath filesep in_qPpaths{1}]);
file_poscs = dir([OutPath filesep in_qPpaths{2}]);

if ~isdir([OutPath filesep out_qPpaths{1}])
    mkdir([OutPath filesep out_qPpaths{1}]);
end
if ~isdir([OutPath filesep out_qPpaths{2}])
    mkdir([OutPath filesep out_qPpaths{2}])
end

nbhood1 = [0 0 0; unique(perms([1 0 0]),'rows'); unique(perms([1 1 0]),'rows'); unique(perms([-1 0 0]),'rows');
    unique(perms([-1 -1 0]),'rows'); unique(perms([-1 -1 1]),'rows'); unique(perms([-1 1 1]),'rows');
    unique(perms([-1 1 0]),'rows'); 1 1 1; -1 -1 -1];
[corners,j] = find(sum(double(~(nbhood1==0)),2)==3);

nbhood2 = [];
for ji = 1:length(corners)
    nbhood2 = [nbhood2; repmat(nbhood1(corners(ji),:),[27 1]) + nbhood1];
end
nbhood = unique(nbhood2,'rows');

numT = 0;
mbeg = 0;
mend = 0;

load(LinearTrans_DTI2T1)
A = inv(M_T1)*inv(M)*M_DTI;

Rn = load(ROI_numFileName);
ROI_num = Rn.ROI_num;
%Get the X Y Z elements of the deformation field and the invM matrix that
%defines the real-world coordinate to image coordinate transform for the
%Template image.
%[X,Y,Z,invM2] = get_transform_params(TransDef,Template_Image);
if ~(matlabpool('size')==4); matlabpool open 4; end
tic
if redo_MNItrans == 1
    parfor j = 3:length(files);
        qf = files(j).name;
        unind = findstr(qf,'_');
        poscf = ['Poscs_' qf(unind+1:end)];
        if mod(j-2,10)==0
            disp(['Transforming Paths in file ' num2str(j-2) ' out of ' num2str(length(files))])
        end
        q = load([OutPath filesep in_qPpaths{1} filesep qf]);
        q = q.q;
        if length(q)==1; flag = ~isempty(q.voxels); else flag = 1; end
        posc = load([OutPath filesep in_qPpaths{2} filesep poscf]);
        posc = posc.posc;
        qmni = [];
        poscmni = NaN(length(q),2);
        if  flag %&& ~exist([out_FT filesep qf],'file')
            %disp(['Transforming Path file ' qf])
            for i = 1:length(q)
                % vcoord = A*[double([q(i).voxels(:,1:end-1) floor(q(i).voxels(:,end))+0.5]); ones(1,size(q(i).voxels,2))];
                vcoord = A*[double(q(i).voxels); ones(1,size(q(i).voxels,2))];
                vcoord = vcoord(1:3,:);
                if ~isempty(vcoord)
                    [~,mni_coord] = nii_map_coords(vcoord, '',NonLinearTrans_seg_inv_sn_T12MNI);
                    qmni(i).voxels = mni_coord;
                    qmni(i).PT = q(i).PT;
                    
                    checkind = repmat(floor(mni_coord(:,1))',[size(nbhood,1) 1]) + nbhood;
                    check_begin = Template_Atlas(sub2ind(size(Template_Atlas),checkind(:,1),checkind(:,2),checkind(:,3)));
                    if ismember(posc(i,1),check_begin)
                        poscmni(i,1) = posc(i,1);
                    else
                        gmroi = find(check_begin);
                        mid_nbhood = checkind(gmroi,:) + 0.5;
                        [ii,minii] = min(sum((repmat(mni_coord(:,1)',[size(mid_nbhood,1),1])-mid_nbhood).^2,2));
                        if ~isempty(ii)
                            poscmni(i,1) = Template_Atlas(checkind(gmroi(minii),1),checkind(gmroi(minii),2),checkind(gmroi(minii),3));
                        else
                            poscmni(i,1) = 0;
                        end
                    end
                    checkind = repmat(floor(mni_coord(:,end))',[size(nbhood,1) 1]) + nbhood;
                    check_begin = Template_Atlas(sub2ind(size(Template_Atlas),checkind(:,1),checkind(:,2),checkind(:,3)));
                    gmroi = intersect(find(check_begin),find(check_begin~=poscmni(i,1)));
                    mid_nbhood = checkind(gmroi,:) + 0.5;
                    [ii,minii] = min(sum((repmat(mni_coord(:,end)',[size(mid_nbhood,1),1])-mid_nbhood).^2,2));
                    if ~isempty(ii)
                        poscmni(i,2) = Template_Atlas(checkind(gmroi(minii),1),checkind(gmroi(minii),2),checkind(gmroi(minii),3));
                    else
                        poscmni(i,2) = 0;
                    end
                    qmni(i).regions = poscmni(i,:);
                    qmni(i).location = [j-2 i];
                    pix = unique(floor(qmni(i).voxels)','rows');
                    qmni(i).ROI = nonzeros(unique(ROI_num(sub2ind(size(ROI_num),pix(:,1),pix(:,2),pix(:,3)))));
                else
                    qmni(i).voxels = [];
                    qmni(i).PT = [];
                end
            end
            numT = numT + length(q);
            mbeg = mbeg + (length(q)-length(find(poscmni(:,1)-posc(:,1))));
            mend = mend + (length(q)-length(find(poscmni(:,2)-posc(:,2))));
        end
        q = qmni;
        posc = poscmni;
        savefilesq([OutPath filesep out_qPpaths{1} filesep qf],q)
        savefilesp([OutPath filesep out_qPpaths{2} filesep poscf],posc)
    end
end
timeforconversion = toc

tic
findfibre_SVN(Zones, Tzones,length(dir([OutPath filesep in_qPpaths{1}]))-2,OutPath,[OutPath filesep out_qPpaths{2}],out_qPpaths);
time_for_findfibre = toc;
tic
parallel_connectivity_saveConn_SVN(Tzones,Zones,OutPath,out_qPpaths{4},atlassize);
time_for_conn = toc;

trk_amass([1 atlassize],[OutPath filesep out_qPpaths{4}],Template_AtlasFileName,[OutPath filesep 'all_brain' num2str(atlassize) '_MNI.trk'],atlassize,1);
