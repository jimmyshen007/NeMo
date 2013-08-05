function [at,filenameout,NumRegions]=colorcode_regions(voxorfilename)

% Takes binary/logical image volume or NIFTI filename and codes regions
% numerically by descending size
% @ LoCastro

if ischar(voxorfilename)
    V=spm_vol(voxorfilename);
    vox=spm_read_vols(V);
else
    vox=voxorfilename;
end

CC=bwconncomp(vox);

region_sizes=zeros(CC.NumObjects,1);
NumRegions=CC.NumObjects;
for i=1:CC.NumObjects
    region_sizes(i)=size(CC.PixelIdxList{i},1);
end

[~,idx]=sort(region_sizes,'descend');

at=zeros(size(vox));

for i=1:CC.NumObjects
    lesion=CC.PixelIdxList{idx(i)};
    at(lesion)=i;
end

if ischar(voxorfilename)
    [d,f,e]=fileparts(voxorfilename);
    filenameout=[d filesep f '_numbered' e];
    V.fname=filenameout;
    spm_write_vol(V,at);
else
    filenameout='';
end

return