% function [Mask_WM_GM, Mask_GM] = makeMasks_prob_SVN(contrast)
function [outGMFileName, outWMFileName, outContFileName] = makeMasks_prob_SVN(contrast, GMSegFileName, WMSegFileName, CSFSegFileName)

% AR modes: aded filenames to be passed as input params, and returned as output params
% Usage : [Mask_WM_aGM, Mask_GM] = makeMasks_prob(contrast)
%
% Description: Using probabilistic segmentations of matter, this function 
%              calculed a three dimensional mask that only possesses nonzero
%              values for white matter and gray matter points.
%
% Input Parameters:
%        contrast : parameter that tunes relative importance between white
%        and gray matter components.
%        
% Output Parameters:
%        Mask_WM_GM: a three dimensional mask that only possesses nonzero
%                   values for white matter and gray matter points. It will 
%                   be saved as image 'Mask_contrast.img' in the chosen 
%                   directory.
%        Mask_GM   : a three dimensional mask that only possesses nonzero
%                   values for gray matter points. It will be saved as image
%                   'Mask_Gray.img' in the chosen directory.
%-------------------------------------------------------------------------- 
% Authors: Yasser Iturria Medina & Pedro Vald�s-Hern�ndez 
% Neuroimaging Department
% Cuban Neuroscience Center
% November 15th 2005
% Version $1.0

if nargin < 1
    contrast = 1;
end
warning off
if nargin < 2
    [GMSegFileName,GMSegFilePath] = uigetfile({'*.img'; '*.nii'},'Reading Gray Matter Segmentation File ...');
    GMSegFileName =  [GMSegFilePath GMSegFileName];
    [WMSegFileName,WMSegFilePath] = uigetfile({'*.img'; '*.nii'},'Reading White Matter Segmentation File ...');
    WMSegFileName =  [WMSegFilePath WMSegFileName];
    [CSFSegFileName,CSFSegFilePath] = uigetfile({'*.img'; '*.nii'},'Reading CSF Matter Segmentation File ...');
    CSFSegFileName =  [CSFSegFilePath CSFSegFileName];
    [pth,nam,ext] = fileparts(GMSegFileName);
    Path = uigetdir('','Directory to save Mask_Gray.img, Mask_White.img and Mask_contrast.img');
else
    [pth,nam,ext] = fileparts(GMSegFileName);
    Path = pth;
end

GMget = spm_vol(GMSegFileName);
WMget = spm_vol(WMSegFileName);
CSFget = spm_vol(CSFSegFileName);
    

GM = single(spm_read_vols(spm_vol(GMget)));
WM = single(spm_read_vols(spm_vol(WMget)));
CSF = single(spm_read_vols(spm_vol(CSFget)));
None = 1 - (GM + WM + CSF);
indGM = find((GM > WM) & (GM > CSF) & (GM > None));
indWM = find((WM > GM) & (WM > CSF) & (WM > None));
Mask_GM = zeros(size(GM));
Mask_WM = zeros(size(GM));
Mask_GM(indGM) = 1;
Mask_WM(indWM) = 1;
% [Mask_GM] = Iso_Remove(Mask_GM, 8);
% [Mask_WM] = Iso_Remove(Mask_WM, 8);
indGM = find(Mask_GM);
indWM = find(Mask_WM);
Mask_huecos = zeros(size(GM));
Mask_huecos(unique([indGM; indWM])) = 1;
Mask_WM_GM = ((contrast*WM + GM)/max(max(max((contrast*WM + GM))))).*Mask_huecos;
Mask_WM_GM = Mask_WM_GM/max(max(max((Mask_WM_GM))));
outContFileName = [Path filesep 'Mask_contrast' ext];
GMget.fname = outContFileName;
spm_write_vol(GMget,Mask_WM_GM);
outGMFileName = [Path filesep 'Mask_Gray' ext];
GMget.fname = outGMFileName;
spm_write_vol(GMget,Mask_GM);
outWMFileName = [Path filesep 'Mask_White' ext];
GMget.fname = outWMFileName;
spm_write_vol(GMget,Mask_WM);
