function matlabbatch=create_coreg_job(refimg,srcimgs)
%Function created for CalcPaths_Bayesian307 function resampling2DMWRI307 (or
%something like that) on the JAKE computer; will replace change_space

%load template batch job struct

if ischar(srcimgs)
    srcimgs=cellstr(srcimgs);
end

load('coreg_job.mat');

S=length(srcimgs);


matlabbatch{1,1}.spm.spatial.coreg.estwrite.ref=cellstr([refimg{1} ',1']);
matlabbatch{1,1}.spm.spatial.coreg.estwrite.source=cellstr([srcimgs{1} ',1']);

if length(srcimgs)>1
	for s=2:S
    		matlabbatch{1,1}.spm.spatial.coreg.estwrite.other{s-1,1}=[srcimgs{s} ',1'];
	end
end
return;
end
