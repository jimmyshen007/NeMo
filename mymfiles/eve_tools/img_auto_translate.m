function img_auto_translate(coT1,prefix,RunAtSize)

%The array coT1 points to the images that require shifting dcm2nii converted images to match T1
%template. RunAtSize gives option to normalize and atlas the image by specifying the size of the atlas

if nargin<3
    RunAtSize=0;
    %emlopt=0;
    if nargin<2
        prefix='backup_';
    end
end


%Set RunAtlasing command

RunAtCmd='RunAtlasing_SVN(coT1{t},RunAtSize,1,0,''t1'')';

if ischar(coT1), coT1=cellstr(coT1); end
if size(coT1,2)~=1, coT1=coT1'; end

tsp=[-0.5; -0.6872; -0.4066];

for t=1:size(coT1,1)
    [d,f,e]=fileparts(coT1{t});    
    V=spm_vol(coT1{t});
    for i=1:3
        V.mat(i,4)=tsp(i)*V.dim(i)*V.mat(i,i);
    end
    if strcmp(e,'.img')
        if ~exist([d filesep prefix f '.hdr'],'file')
            copyfile([d filesep f '.hdr'],[d filesep prefix f '.hdr']);
        end
    elseif strcmp(e,'.nii')
        if ~exist([d filesep prefix f '.nii'],'file')
            copyfile([d filesep f '.nii'],[d filesep prefix f '.nii']);
        end
    end
    spm_create_vol(V);
    disp(f);
    if RunAtSize~=0
        eval(RunAtCmd);
    end
end

end
