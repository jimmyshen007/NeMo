function [] = processFreeSurfer(FSdir,FSroi,NEWroi)

wm_fn = deblank(Dir2Arr(FSdir,'wm.nii'));
ribbon_fn = deblank(Dir2Arr(FSdir,'ribbon.nii'));
seg_fn = deblank(Dir2Arr(FSdir,'aparc+aseg.nii'));

wm_hdr = spm_vol(wm_fn);
wm = spm_read_vols(wm_hdr);
ribbon = spm_read_vols(spm_vol(ribbon_fn));
seg = spm_read_vols(spm_vol(seg_fn));

wm_mask = zeros(size(wm));
wm_mask(find(wm>0&wm<250)) = 1;
[pth,fn,ext] = fileparts(wm_hdr.fname);

if ~isdir([pth filesep 'Segmented']); mkdir([pth filesep 'Segmented']); end
wm_hdr.fname = [pth filesep 'Segmented' filesep 'Mask_White.nii']
wm_hdr.private.dat.fname = wm_hdr.fname;
spm_create_vol(wm_hdr)
spm_write_vol(wm_hdr,wm_mask)

old_roi = load(FSroi);
new_roi = load(NEWroi);
atlas = zeros(size(seg));
for i = 1:new_roi(end);
    atlas(find(seg==old_roi(i))) = new_roi(i);
end
if ~isdir([pth filesep 'Atlased' num2str(new_roi(end))]); mkdir([pth filesep 'Atlased' num2str(new_roi(end))]); end
wm_hdr.fname = [pth filesep 'Atlased' num2str(new_roi(end)) filesep 'T1_atlas' num2str(new_roi(end)) '.nii'];
wm_hdr.private.dat.fname = wm_hdr.fname;
spm_create_vol(wm_hdr)
spm_write_vol(wm_hdr,atlas)

gm_mask = zeros(size(atlas));
gm_mask(find(atlas)) = 1;
wm_hdr.fname = [pth filesep 'Segmented' filesep 'Mask_Gray.nii'];
wm_hdr.private.dat.fname = wm_hdr.fname;
spm_create_vol(wm_hdr)
spm_write_vol(wm_hdr,gm_mask)