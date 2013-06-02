function [NewSourceImage,NewAddImage] = approx_coreg2MNI(RefImage,SourceImage,fileprefix,AddImages)

%Approximately align (with a linear transformation) the source and the
%template image. SourceImage can be a list of files.

if nargin<4 || isempty(fileprefix)
    fileprefix = 'a';
end

if size(SourceImage,1)==1
    [t1pth,t1fn,ext]=fileparts(SourceImage);
    NewSourceImage = [t1pth filesep fileprefix t1fn ext];
    %First make a copy of the original T1 header and image
    if strcmp(ext,'.img')
        copyfile([t1pth filesep t1fn '.hdr'],[t1pth filesep fileprefix t1fn '.hdr']);
    end
    copyfile([t1pth filesep t1fn ext],[t1pth filesep fileprefix t1fn ext])
    %Coregister the source to the template.
    x = spm_coreg(RefImage,NewSourceImage);
    M = inv(spm_matrix(x)); % RefImage coords = M*SourceImage coords
    save([t1pth filesep 'ApproxTrans_' t1fn '_to_MNI.mat'],'M')
    spm_get_space(NewSourceImage,M*spm_get_space(NewSourceImage));
    if nargin>3
        for j = 1:size(AddImages,1)
            [t1pth,t1fn,ext]=fileparts(deblank(AddImages(j,:)));
            NewAddImage = [t1pth filesep fileprefix t1fn ext];
            %First make a copy of the original T1 header and image
            if strcmp(ext,'.img')
                copyfile([t1pth filesep t1fn '.hdr'],[t1pth filesep fileprefix t1fn '.hdr']);
            end
            copyfile([t1pth filesep t1fn ext],[t1pth filesep fileprefix t1fn ext])
            spm_get_space(NewAddImage,M*spm_get_space(NewAddImage));
        end
    else
        NewAddImage = {};
    end
else
    NewSourceImage = {};
    for i = 1:length(SourceImage)
        ssi = SourceImage{i};
        [t1pth,t1fn,ext]=fileparts(ssi);
        NewSourceImage{i} = [t1pth filesep fileprefix t1fn ext];
        %First make a copy of the original T1 header and image
        if strcmp(ext,'img')
            copyfile([t1pth filesep t1fn '.hdr'],[t1pth filesep fileprefix t1fn '.hdr']);
        end
        copyfile([t1pth filesep t1fn ext],NewSourceImage{i})
        %Coregister the source to the template.
        x = spm_coreg(RefImage,NewSourceImage{i});
        M = inv(spm_matrix(x)); % RefImage coords = M*SourceImage coords
        spm_get_space(NewSourceImage,M*spm_get_space(NewSourceImage{i}));        
    end
end