function [VODF_fname, ODFdirsTextFileName] = compute_QBI_AR(MaskFileName, GradFileName, DiffusionDataNames)

% apparently DiffusionDataNames refers to diffusion weighted data - we will assume b0 is also included in this list (is 1st entry)

% global MaskFileName MaskFilePath GradFileName GradFilePath DiffusionDataNames;
global ODF_datatype mxODF;
ODF_datatype = 'uint8';
%mxODF = 0.3;

%% Reading Mask and Diffusion gradient directions:
warning off
if nargin==0
    [MaskFileName,MaskFilePath] = uigetfile({'*.img'; '*.nii'},'Reading Mask_Contrast ...');
    MaskFileName = [MaskFilePath MaskFileName];
    [GradFileName, GradFilePath] = uigetfile({'*.txt';'*.cod'},'Select gradient directions');
    GradFileName = [GradFilePath GradFileName];
    [DiffusionDataNames,path] = myuigetfiles('*.img; *.nii', 'Select diffusion weighted images (in the correct order; select b0 as 1st image)');
    for i=1:length(DiffusionDataNames)
        DiffusionDataNames{i} = [path DiffusionDataNames{i}];
    end
end
[pth,nam,ext] = fileparts(MaskFileName);
MaskVol = load_nii(MaskFileName);
Mask = MaskVol.img;
% MaskVol = spm_vol( MaskFileName);
% Mask = spm_read_vols(MaskVol);
grad = load(GradFileName);
if size(grad,1)==3;
    grad = grad';
end
% below looks redundant...
% filename1 = GradFileName;
% pathname1 = GradFilePath;
% if filename1 ~= 0
%     grad = load([pathname1 filename1]);
% end

pos = ismember(grad, [0 0 0], 'rows');
ind_zeros = find(pos);
if ~isempty(ind_zeros), grad(ind_zeros,:) = []; end
Ns = length(grad);

%% Selecting Diffusion weighted data:
thisdir = pwd;
cd(thisdir);

names = char(DiffusionDataNames{2:end});
if size(names,1) ~= Ns
    errordlg('The number of diffusion gradients and the number of images are different');
    return
end
Vdata = load_nii(names(1,:));

ODF_voxel_size = Vdata.hdr.dime.pixdim(2:4);
szVdata = size(Vdata.img);
[pt,nm,ex] = fileparts(names(1,:));
VODF_fname = fullfile(pt,['ODF_' nm ext]);
%% Computing the design matrix:
[V, F] = getFV(2, 'oct');
Dir = V;
ODFdirsTextFileName = [pt filesep 'ODF_directions.txt'];
eval(['save ''' pt filesep 'ODF_directions.txt'' Dir -ascii'])
gradn = V;

%-----------------------------------------------%

[Lmax Nmin] = calcula_Lmax(grad);
if Lmax >= 6
    Lmax = 6;
end

[ph, th, r] = cart2sph(grad(:,1), grad(:,2), grad(:,3));
th = pi/2 - th;
ind = find(ph(:) < 0);
ph(ind) = 2*pi - abs(ph(ind));

cc = .3;
Ymatrix = zeros(length(th),1); Laplacian = []; Sharp = [];
for i=0:2:Lmax
    for j=-i:i
        Ymatrix = [Ymatrix  Ylm(i, j, th, ph)];
        Laplacian = [Laplacian (i^2)*(i + 1)^2];
        Sharp = [Sharp (1 + cc*i*(i+1) )];
    end
end
Laplacian = diag(Laplacian);
Ymatrix(:,1) = [];
%-----------------------%
[ph, th, r] = cart2sph(gradn(:,1), gradn(:,2), gradn(:,3));
th = pi/2 - th;
ind = find(ph(:) < 0);
ph(ind) = 2*pi - abs(ph(ind));

YmatrixR = [];K  = [];
for L=0:2:Lmax
    for m=-L:L
        YmatrixR = [YmatrixR  Ylm(L, m, th, ph)];
        factor1 = ((-1)^(L/2))*doublefactorial(L+1)./( (L+1)*doublefactorial(L) );
        K = [K; factor1];
    end
end
%--------------------------%
% cond = 0; % without regularize
cond = 1; % regularizing
%-----------------------------------------------%
MX = pinv(Ymatrix'*Ymatrix)*Ymatrix';
%MX = pinv(Ymatrix);
if cond == 1
    MX = diag(Sharp)*MX;
end
MX = diag(K)*MX;
MX = YmatrixR*MX;

%% Computing ODF and writing image:
nrows = szVdata(1); ncols = szVdata(2); nslices = szVdata(3); nDirs = size(V,1);
L1imin = 1; L1imax = nrows;
L2imin = 1; L2imax = ncols;
ODF_mat = zeros([szVdata size(V,1)], ODF_datatype);
mxODF = 0;
for v=1:Ns
    Vdata = load_nii(names(v,:));
    tmp = Vdata.img;
    mxODF = max(mxODF, max(tmp(:)));
end
for v=1:Ns
    Vdata = load_nii(names(v,:));
    tmp = Vdata.img;
    if strcmp(ODF_datatype, 'int16') || strcmp(ODF_datatype, 'uint8')
        %tmp = min(mxODF, tmp);
        tmp = single(tmp)/single(mxODF)*single(intmax(ODF_datatype));
    end
    ODF_mat(:,:,:,v) = eval([ODF_datatype '(tmp)']);
end
clear Vdata tmp;
% for slice = 1:3
for slice = 1:nslices
    disp(['slice -> ' num2str(slice)]);
    ODF_slice = zeros(nDirs, nrows*ncols, ODF_datatype);
    Mascara = squeeze(Mask(:,:,slice));
    data_slice = squeeze(ODF_mat(:,:,slice,1:Ns));
    %size(data_slice), nrows, ncols, Ns,
    data_slice = (reshape(data_slice, [nrows*ncols, Ns])).';
    I = find(Mascara(:));
    dd = single(data_slice(:,I));
    sY = abs(MX*dd);
    Like = sY; 
    mny = min(sY,[], 1);
    for rr = 1:nDirs
        q = Like(rr,:)-mny;
        Like(rr,:) = q;
    end
    smy = sum(sY,1)+eps;
    for rr = 1:nDirs
        Like(rr,:) = Like(rr,:)./smy;
    end
    Like = Like*single(intmax(ODF_datatype));
    if strcmp(ODF_datatype, 'int16') || strcmp(ODF_datatype, 'uint8')
        Like = min(single(intmax(ODF_datatype)), Like);
    end
    ODF_slice(:,I) = eval([ODF_datatype '(Like)']);
    ODF_mat(:,:,slice,:) = reshape(ODF_slice.', [nrows,ncols,nDirs]);
end
ODFVol = make_nii(ODF_mat, ODF_voxel_size, [], spm_type(ODF_datatype));
% ODFVol = Vdata;
% ODFVol.fileprefix = sprintf('%s%s', 'ODF-', Vdata.fileprefix);
% ODFVol.img = ODF_mat;
save_nii(ODFVol, VODF_fname);
cd(thisdir);

