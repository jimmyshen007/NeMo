function matlabbatch=create_coreg_job_struct(refimg,srcimgs)
%Function created for CalcPaths_Bayesian307 function resampling2DMWRI307 (or
%something like that); will replace change_space
%Note: all fields MUST be in cell format, or else SPM will choke

%load template batch job struct
load('rs_seg2dti_job.mat');

matlabbatch{1,1}.spm.spatial.coreg.estwrite.ref=cellstr([refimg ',1']);
matlabbatch{1,1}.spm.spatial.coreg.estwrite.source=cellstr([srcimgs{5} ',1']);
matlabbatch{1,1}.spm.spatial.coreg.estwrite.other{1,1}=[srcimgs{1} ',1'];
matlabbatch{1,1}.spm.spatial.coreg.estwrite.other{2,1}=[srcimgs{2} ',1'];
matlabbatch{1,1}.spm.spatial.coreg.estwrite.other{3,1}=[srcimgs{3} ',1'];
matlabbatch{1,1}.spm.spatial.coreg.estwrite.other{4,1}=[srcimgs{4} ',1'];

return;
end
