function [AtlasFileName,AtlasTextFileName]=ChooseAtlas(atlassize,Atlasing_dir)

if atlassize==116
    AtlasTextFileName = [Atlasing_dir filesep 'atlas116.cod'];
    AtlasFileName = [Atlasing_dir filesep 'atlas116.img'];
elseif atlassize==307
    AtlasTextFileName = [Atlasing_dir filesep 'Subject1_Atlas_WedgeLobe.cod'];
    AtlasFileName = [Atlasing_dir filesep 'Subject1_Atlas_WedgeLobe.img'];
elseif atlassize==69
    AtlasTextFileName = [Atlasing_dir filesep 'atlas69.cod'];
    AtlasFileName = [Atlasing_dir filesep 'Atlas69.img'];
elseif atlassize==71
    AtlasTextFileName = [Atlasing_dir filesep 'atlas71.cod'];
    AtlasFileName = [Atlasing_dir filesep 'atlas71.img'];
elseif atlassize==84
    AtlasTextFileName = [Atlasing_dir filesep 'atlas84.cod'];
    AtlasFileName = [Atlasing_dir filesep 'atlas84.img'];
elseif atlassize==90
    AtlasTextFileName = [Atlasing_dir filesep 'atlas90.cod'];
    AtlasFileName = [Atlasing_dir filesep 'atlas90.img'];
elseif atlassize==281
    AtlasTextFileName = [Atlasing_dir filesep 'atlas281.cod'];
    AtlasFileName = [Atlasing_dir filesep 'atlas281.nii'];
elseif atlassize==86
    AtlasTextFileName = [Atlasing_dir filesep 'atlas86.cod'];
    AtlasFileName = [Atlasing_dir filesep 'atlas86.img'];
else
    error('CalcPaths_Bayesian_generic:atlassizechk','Yo yo yo, you gonna specify an atlas, homey?');
end

return;
end