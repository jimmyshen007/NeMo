function [GMSegFile, WMSegFile, CSSegFile, AtlasedFile]=AutoLabel_Only_SVN(T1FileName,pth,NewAtlasSize,Atlas_Output_dir,OVERRIDE)

%For images that have already been normalized and segmented, and atlased
%using a given atlas, performs the atlasing portion with a different atlas.
%pth should specify the directory that contains the 'Normalized' and
%'Segmented' folders

startup_varsonly;

[d,fn,ext]=fileparts(deblank(T1FileName));

%Try to account for pth as relative directory
if isempty(pth)
    pth=d;
elseif ~exist(pth,'dir')
    pth=[d filesep pth];
end

[AtlasFile,~]=ChooseAtlas(NewAtlasSize,Atlasing_dir);

if nargin < 4 || isempty(Atlas_Output_dir)
    Atlas_Output_dir=[pth filesep 'Atlased' num2str(NewAtlasSize)];
    OVERRIDE=0;
end


Thresh=[];

%%%%%%%%%---------------  Atlasing --------------------%%%%%%%%%%%%%%%%

disp('Atlasing ...');

GMSegFile = [pth filesep 'Segmented' filesep 'c1' fn ext];
WMSegFile = [pth filesep 'Segmented' filesep 'c2' fn ext];
CSSegFile = [pth filesep 'Segmented' filesep 'c3' fn ext];


if exist(Atlas_Output_dir,'dir')
    if ~OVERRIDE
        if exist([Atlas_Output_dir filesep fn '_Atlas' ext],'file') %&& exist([Atlas_Output_dir filesep fn '_Vol.mat'],'file')
            if ~exist([Atlas_Output_dir filesep fn '_Vol.mat'],'file') %Check to see that Volumetrics have been done
                AtlasedFile=[Atlas_Output_dir filesep fn '_Atlas' ext];
                [d,f,~] = fileparts(deblank(AtlasedFile));
                at = spm_read_vols(spm_vol(deblank(AtlasedFile)));
                at_vol = zeros(NewAtlasSize,1);
                for j = 1:NewAtlasSize
                    at_vol(j) = length(find(at==j));
                end
%                 save([d filesep f(1:end-6) '_Vol.mat'],'at_vol');
                save([Atlas_Output_dir filesep fn '_Vol.mat'],'at_vol');
            end
            AtlasedFile=[Atlas_Output_dir filesep fn '_Atlas' ext];
            disp('File already exists.');
            return;
        end
    end
else
    mkdir(Atlas_Output_dir);
end
    

Transf_matname = fullfile([pth filesep 'Normalized'], [fn '_seg_sn.mat']);
AtlasedFile = Auto_Labelling(GMSegFile, WMSegFile, CSSegFile, AtlasFile, Transf_matname, Atlas_Output_dir,Thresh);

%%%%%%%%%--------- End of the Atlasing Step -----------%%%%%%%%%%%%%%%%


%Volumetrics
[d,f,~] = fileparts(deblank(AtlasedFile));
at = spm_read_vols(spm_vol(deblank(AtlasedFile)));
at_vol = zeros(NewAtlasSize,1);
for j = 1:NewAtlasSize
    at_vol(j) = length(find(at==j));
end
save([d filesep f(1:end-6) '_Vol.mat'],'at_vol');
 
return;
end
