function atFileName = getAtlas(atlassize)

if atlassize == 116
    atFileName = ['resource' filesep 'at116' filesep 'sacos007a1001_Atlas.img'];
elseif atlassize == 86
    atFileName = ['resource' filesep 'at86' filesep 'saT1_atlas86.nii'];
else
    disp('Invalid atlas size selection');
    atFileName='';
end