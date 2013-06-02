clear all
close all
%matlabpool open 3

SVNEveTools = '/home/amy/Desktop/Matlab2009a/trunk';
AtlasingDir = '/home/amy/Desktop/Matlab2009a/StructrualConnectivity_Bayes/Atlasing';
SPMdir8 = '/home/amy/Desktop/Matlab2009a/spm8';
niitools = '/home/amy/Desktop/Matlab2009a/StructrualConnectivity_Bayes/nifti_toolbox';
ODFtools = '/home/amy/Desktop/Matlab2009a/StructrualConnectivity_Bayes/ODF_Qball';
SCBtrack = '/home/amy/Desktop/Matlab2009a/StructrualConnectivity_Bayes/Tracking_Bayes';
mymfiles = '/home/amy/Desktop/Matlab2009a/mymfiles';

addpath(SCBtrack)
addpath(mymfiles)
p = genpath(SVNEveTools);
addpath(p);
p3 = genpath(SPMdir8);
addpath(p3);
addpath(AtlasingDir)
addpath(niitools)
addpath(ODFtools)
addpath('/home/amy/Desktop/Matlab2009a/BCT');

yeo_dir = '/home/amy/Desktop/Matlab2009a/StructrualConnectivity_Bayes/Atlasing/Yeo_JNeurophysiol11_MNI152';
yeo_labels = spm_read_vols(spm_vol([yeo_dir filesep 'rYeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii']));
fs_labels = spm_read_vols(spm_vol([yeo_dir filesep 'FSatlas86.nii']));

nl = NaN(86,1);
for i =1:86;
    yl = yeo_labels(find(fs_labels==i));
    hist(yl,0:1:7)
    [h,xx] = hist(yl,0:1:7);
    [~,iw] = max(h);
    %nl = input(['Suggested label is ' num2str(xx(iw)) ', which label should this region get?']);
    nl_v(i) = xx(iw);
    pause;
end
% pats = Dir2Arr('/home/amy/work/BTFdata/Normals','e00*');
% all_ati = [];
% for i = 1:size(pats,1);
%     nat = Dir2Arr([deblank(pats(i,:)) filesep 'T1scan_FS' filesep 'Atlased86'],'wa*.nii');
%     if ~isempty(nat);
%         ati = spm_read_vols(spm_vol(nat));
%         all_ati = [all_ati ati(:)];
% %         for j = 1:(length(unique(ati))-1);
% %             yl = yeo_labels(find(ati==j));
% %             hist(yl)
% %         end
%     end
% end
% 
% win_at = zeros(size(all_ati,1),1);
% for i = 1:size(all_ati,1);
%     ulab = unique(all_ati(i,:));
%     if length(ulab)==1
%         win_at(i) = ulab;
%     else
%         kk = [];
%         for j = 1:length(ulab);
%             kk = [kk length(find(ulab(j)==all_ati(i,:)))];
%         end
%         [~,iw] = max(kk);
%         win_at(i) = ulab(iw);
%     end
% end