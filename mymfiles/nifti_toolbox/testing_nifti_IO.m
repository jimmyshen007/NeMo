%function testing_nifti_IO
% script will read in a float 32 nifti volume, create a complex64 nifti
% volume and write it out to disk
% uses nifti_toolbox written in MATLAB by Jimmy Shen

% clear all;
flname = 'L:\CIND_Frank_data\T1.nii';

nii = load_nii(flname, 1);
nii.hdr.dime,
q = double(nii.img); 
clear nii;
[nrows, ncols, nslices] = size(q),
max(max(max(abs(q)))),
min(min(min(abs(q)))),

figure; 
subplot(2,2,1); colormap(gray); imagesc(abs(squeeze(q(40,:,:)))); 
subplot(2,2,2); colormap(gray); imagesc(abs(squeeze(q(80,:,:)))); 
subplot(2,2,3); colormap(gray); imagesc(abs(squeeze(q(120,:,:)))); 
subplot(2,2,4); colormap(gray); imagesc(abs(squeeze(q(160,:,:)))); 


% % If need to display the 3D FFT of this volume, uncomment below
% qq = q;
% for j = 1:ncols
%     for k = 1:nslices
%         qq(:,j,k) = single(myifft(q(:,j,k)));
%     end
% end
% 
% for i = 1:nrows
%     tmpslice = squeeze(qq(i,:,:));
%     qq(i,:,:) = single(myifft2(tmpslice));
% end
% 
% max(max(max(abs(qq)))),
% min(min(min(abs(qq)))),
% 
% figure; 
% subplot(2,2,1); colormap(gray); imagesc(abs(squeeze(qq(40,:,:)))); 
% subplot(2,2,2); colormap(gray); imagesc(abs(squeeze(qq(80,:,:)))); 
% subplot(2,2,3); colormap(gray); imagesc(abs(squeeze(qq(120,:,:)))); 
% subplot(2,2,4); colormap(gray); imagesc(abs(squeeze(qq(160,:,:)))); 

 
%% Create a Complex volume and write it out to disk
[pathstr, name, ext, versn] = fileparts(flname);
flname_cplx = fullfile(pathstr, [name '_cplx' ext]);
datatype = 1792; 
% From Nifti standard: datatype (optional):	Storage data type:
%  		2 - uint8,  4 - int16,  8 - int32,  16 - float32,
%  		32 - complex64,  64 - float64,  128 - RGB24,
%  		256 - int8,  512 - uint16,  768 - uint32, 
%  		1792 - complex128

qq = q + sqrt(-1)*q;        % creates a new complex volume whose real and imag are equal to the original real volume
nii = make_nii(qq, [],[],datatype);      % creates the nifti structure with header info and image data
save_nii(nii, flname_cplx);
