function compute_QBI

global MaskFileName MaskFilePath GradFileName GradFilePath DiffusionDataNames;

%% Reading Mask:
warning off
if isempty(MaskFileName)
    [MaskFileName,MaskFilePath] = uigetfile({'*.img'; '*.nii'},'Reading Mask_Contrast ...');
end
[pth,nam,ext] = fileparts([MaskFilePath MaskFileName]);
MaskVol = spm_vol([MaskFilePath MaskFileName]);
Mask = spm_read_vols(MaskVol);

%% Selecting Diffusion gradient directions:
if isempty(GradFileName)
    [GradFileName, GradFilePath] = uigetfile({'*.txt';'*.cod'},'Select gradient directions');
end
filename1 = GradFileName;
pathname1 = GradFilePath;
if filename1 ~= 0
    grad = load([pathname1 filename1]);
end
pos = ismember(grad, [0 0 0], 'rows');
ind_zeros = find(pos);
if ~isempty(ind_zeros), grad(ind_zeros,:) = []; end
Ns = length(grad);

%% Selecting Diffusion weighted data:
if isempty(DiffusionDataNames)
    if strcmp(spm('ver'),'SPM2')
        DiffusionDataNames = spm_get(Inf,'Select diffusion weighted images (in the correct order)');
    elseif strcmp(spm('ver'),'SPM5')
        DiffusionDataNames = spm_select(Inf,'image', 'Select diffusion weighted images (in the correct order)');
    end
end
names = DiffusionDataNames;
if size(names,1) ~= Ns
    errordlg('The number of diffusion gradients and the number of images are different');
    return
end
for v=1:Ns
    Vdata(v) = spm_vol(names(v,:));
end
[pt,nm,ex] = fileparts(names(1,:));

%% Creating ODF image volume:
[V, F] = getFV(2);
for i = 1:size(V,1)
    VODF(i) = MaskVol(1);
    VODF(i).fname = fullfile(pt,['ODF_' nm ext]);
    VODF(i).pinfo = [1 0 0]';
    if strcmp(spm('ver'),'SPM2')
        VODF(i).dim(4) = spm_type('float');
        VODF(i).n = i;
    elseif strcmp(spm('ver'),'SPM5')
        VODF(i).dt = [16 0];
        VODF(i).n = [i 1];
    end    
end
VODF = spm_create_vol(VODF); Dir = V;
global VODF;
eval(['save ' pt filesep 'ODF_directions.txt Dir -ascii']);
gradn = V;

%-----------------------------------------------%

[Lmax Nmin] = calcula_Lmax(grad);
if Lmax >= 6
    Lmax = 6;
end

%% Computing the design matrix:
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

%% Computing ODF and writing image:
% for slice = 1:3
for slice = 1:Vdata(1).dim(3)
    disp(['slice -> ' num2str(slice)])
    Mxy = [1     0     0     0;
        0     1     0     0;
        0     0     1  slice;
        0     0     0     1];
    for i = 1:Ns
        Data(:,:,i) = spm_slice_vol(Vdata(i),Mxy,Vdata(i).dim(1:2),0);
    end
    if ~isempty(ind_zeros), Data(:,:,ind_zeros) = []; end    
    Mascara = squeeze(Mask(:,:,slice));
    if slice == 1
        [f, c, Ndata] = size(Data);
        L1imin = 1; L1imax = f;
        L2imin = 1; L2imax = c;
    end
    ODF_slice = zeros(f,c,size(V,1), 'single');
    for j=1:f
        for k=1:c
            if Mascara(j,k) > 0
                %-----------------------------------------------------%
                y = squeeze(Data(j,k,:));
                MX = Ymatrix'*Ymatrix;
                sY = pinv(MX)*Ymatrix'*y;

                if cond == 1
                    aSharp = Sharp'.*sY;
                elseif cond == 0
                    aSharp = sY;
                end

                ss = aSharp.*K;
                ODF = abs(YmatrixR*ss);
                %-----------------------------------------------------%
                Like = ODF - min(ODF);
                Like = Like/(sum(Like)+eps);
                ODF_slice(j,k,:) = single(Like);
            end
        end
    end
    for i = 1:size(V,1),
        spm_write_plane(VODF(i),ODF_slice(:,:,i),slice);
    end
end
if strcmp(spm('ver'),'SPM2')
    spm_close_vol(VODF);
end
return
