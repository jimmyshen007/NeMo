function [AtlasedFileName, VAtlas] = Auto_Labelling(GMSegFile, WMSegFile, CSSegFile, AtlasFile, Transf_matname, Output_dir, Thresh)
%
% Syntax : 
% Atlas =   Auto_Labelling(GMSegFile, WMSegFile, CSSegFile, AtlasFile,
% transf_matname)
%
% This function returns an individual Atlas in Analyze format. 
% The first step is to find the voxels that belongs to gray matter tissue 
% using the tissues probabilities maps previously obatained by the segmentation
% process. Due to the thresholding process, some holes, as well as isolated 
% points are present at gray matter volume. To solve this problem is used the 
% matlab function "imfill" to refills the internal holes and an internal
% function to reduce the isolated points(Uncoment if you want to employ it).
% Then each gray matter voxel is labeled with one structure label using the
% space transformation matrix (obtained in the normalization step), and an 
% anatomical atlas (constructed by manual segmentation for a group of
% subjects).
%
% Input Parameters:
%  GMSegFile      : Gray Matter Segmentation File.
%  WMSegFile      : White Matter Segmentation File.
%  CSSegFile      : Cerebral Spinal Fluid Segmentation File.
%  AtlasFile      : Reference Atlas File.
%  Transf_matname : Normalization Transform File.
%  Output_dir     : Output directory for segmented, normalized and atlased
%                    files. If the user doesn't change the output directory,
%                    the resulting files are saved in the same address than the
%                    Gray Matter Segmentation File.
%  Thresh         : Threshold for gray matter segmentation. 
%                   Just the voxels with  higher probability than the threshold
%                   are taken into acount in the automatic labelling step. If 
%                   the threshold isn't specified then an automatic one is taken. 
%                   All voxels with higher gray matter probabillity than 1-(GM+WM+CSF)
%                   are taken into account in the automatic labelling step.
%                   Being:
%                   GM(A voxel V belongs to gray matter tissue with a probability GM).
%                   WM(A voxel V belongs to white matter tissue with a probability WM). 
%                   CSF(A voxel V belongs to cerebral spinal fluid with a probability CSF).
% Output Parameter:
%   Atlas : Individual Gray Matter Atlas.
%
% Related References:
% 1.- Ashburner J, Friston K. Multimodal image coregistration and partitioning-- 
%     a unified framework. Neuroimage. 1997 Oct;6(3):209-17. 
% 3.- Evans AC, Collins DL, Milner B (1992). An MRI-based Stereotactic Brain
%     Atlas from 300 Young Normal Subjects, in: Proceedings of the 22nd Symposium 
%     of the Society for Neuroscience, Anaheim, 408.
%
% Note: The morphological treatment(erosion and dilation) are comented due 
%       to its dependence of the results with the structure element used to 
%       erode and dilate.      
%
% See also: imfill  imerode  imdilate  spm_segment atlasing spm_normalise
%__________________________________________________________________________
% Authors:  Yasser Alem�n G�mez & Lester Melie Garc�a  
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2005
% Version $1.0

%=====================Checking Input Parameters===========================%
if nargin==0
    [GMSegFile,sts] = spm_select([1],'image','Selecting Gray Matter Images','',cd);
    [WMSegFile,sts] = spm_select([1],'image','Selecting White Matter Images','',cd);
    [CSSegFile,sts] = spm_select([1],'image','Selecting Cerebral Spinal Fluid Images','',cd);
    [AtlasFileName,AtlasFilePath] = uigetfile({'*.img'},'Reading Reference Atlas File ...');
    AtlasFile = [AtlasFilePath AtlasFileName];
    [Transf_matname,sts] = spm_select([1],'mat','Selecting Normalization Transform Files','',cd);
    Output_dir{1,1} = uigetdir(cd ,'Selecting Output Directory');
else
    if isempty(GMSegFile)
        [GMSegFile,sts] = spm_select([1],'image','Selecting Gray Matter Images','',cd);
    end;
    if isempty(WMSegFile)
        [WMSegFile,sts] = spm_select([1],'image','Selecting White Matter Images','',cd);
    end;
    if isempty(CSSegFile)
        [CSSegFile,sts] = spm_select([1],'image','Selecting Cerebral Spinal Fluid Images','',cd);
    end;
    if isempty(AtlasFile)
        [AtlasFileName,AtlasFilePath] = uigetfile({'*.img'},'Reading Reference Atlas File ...');
        AtlasFile = [AtlasFilePath AtlasFileName];
    end;
    if isempty(Transf_matname)
        [Transf_matname,sts] = spm_select([1],'mat','Selecting Normalization Transform Files','',cd);
    end;
    if isempty(Output_dir)
        Output_dir{1,1} = uigetdir(cd ,'Selecting Output Directory');
    end
end
fclose('all');
%=========================================================================%
%
%=========================Main program=====================================
spm_defaults;
global defaults
defaults.analyze.flip = 0; % Not flipped Images

VSeg =  spm_vol(GMSegFile);
VWM = spm_vol(WMSegFile);
VCS = spm_vol(CSSegFile);
[Apth,Aname,Aext] = fileparts(AtlasFile);  % Here the matrix Iatlas is loaded.
Vatlas = spm_vol(AtlasFile);Iatlas = spm_read_vols(Vatlas); Atlastype = Aext(2:end);

%%%%%%%%%%%-----------------Thresholding-------------------%%%%%%%%%%%%%%%%
 Iseg_mask = zeros(VSeg.dim(1),VSeg.dim(2),VSeg.dim(3));
 AllBrainMask = Iseg_mask;
 eThresh = exist('Thresh', 'var');
 for z = 1:VSeg.dim(3)
     %% tissues
     GM  = spm_slice_vol(VSeg,spm_matrix([0 0 z]),VSeg.dim(1:2),0);
     WM  = spm_slice_vol(VWM,spm_matrix([0 0 z]),VWM.dim(1:2),0);
     CSF = spm_slice_vol(VCS,spm_matrix([0 0 z]),VCS.dim(1:2),0);
     None = 1 - (GM + WM + CSF);
     %% inds
     Thresh2 = Thresh;
     size(GM); size(WM); size(CSF); size(None); size(Thresh);
     if (eThresh==0) || (isempty(Thresh)) || (Thresh == 0)
         indGM = find((GM > WM) & (GM > CSF) & (GM > None));
     else
         indGM  = find((GM>CSF)&(GM>WM)&(GM>=Thresh));
     end
     %%%% --------------- Creating Brain Mask  ------------------%%%%
     %%%%  Lester Melie  April 11th 2006
     indWM = find((WM > GM) & (WM > CSF) & (WM > None));  %
     %indCSF = find((CSF > GM) & (CSF > WM) & (CSF > None));
     if (~isempty(indGM))&&(~isempty(indWM)) %&(~isempty(indCSF))
         %WM(unique([indGM; indWM; indCSF])) = 1;
         WM(unique([indGM; indWM])) = 1;
         AllBrainMask(:,:,z) = WM;
     end;
     %% brain
     Brain = sparse(VSeg.dim(1),VSeg.dim(2));
     Brain(indGM) = 1;
     Iseg_mask(:,:,z) = full(Brain);
end
%%%%%%%%%-------------- End of Thresholding ---------------%%%%%%%%%%%%%%%%
disp(['Refilling...']);  
Iseg_mask = imfill(logical(Iseg_mask),'holes');
clear IWM ICS Iseg ind WM CSF None;
%Iseg_mask = Iso_Rem(Iseg_mask,7);
tic;VAtlas = Get_Norm_Coord(Atlastype,AtlasFile,Transf_matname,Iseg_mask, Iatlas, Output_dir);toc;
AtlasedFileName = VAtlas.fname;
%========================End of main program==============================%
return;

%=====================Internal functions==================================%
function Vol = Get_Norm_Coord(Atlastype, AtlasFile, Matname, Iseg_mask, Iatlas, Output_dir)
%
%
% Input Parameters:
%   Atlastype   : Gray Matter Segmentation File
%   AtlasFile   : Reference Atlas File
%   matname     : Normalisation Transform File
%   Iseg_mask   : White Matter Segmentation File
%   Iatlas      : Cerebral Spinal Fluid Segmentation File 
%   Output_dir  : Output Directory
%
% Note: This is based on the function 'spm_write_sn' developed by PhD.John Ashburner(FIL,UCL).
%
%__________________________________________________________________________
% Authors:  Yasser Alem�n G�mez & Pedro Vald�s Hern�ndez  
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2005
% Version $1.0
warning off
global defaults
defaults.analyze.flip = 0; % Not flipped Images
Vatlas = spm_vol(AtlasFile);
if length(Vatlas.dim)==4
    dt = [Vatlas.dim(4) 0];
elseif length(Vatlas.dim)==3
    dt = Vatlas.dt;
end
% VAtlas_mat = Vatlas.mat;
load('-mat',Matname); 
[Vpth,Vname,Vext] = fileparts(VF.fname); mat = VF.mat; 
if strcmp(spm('ver'),'SPM2')
dim = [VF.dim(1:3) dt(1)];
elseif strcmp(spm('ver'), 'SPM5') || strcmp(spm('ver'), 'SPM8')
    dim = [VF.dim(1:3)];
end
Vol = struct('fname','','dim',dim,'mat',mat,'pinfo',[1 0 0]',...
    'descrip','Atlas image','dt',dt);
Vol.fname =[char(Output_dir) filesep Vname '_Atlas' Vext];
Vol = spm_create_vol(Vol);
x = 1:VG(1).dim(1); y = 1:VG(1).dim(2); z = 1:VG(1).dim(3);
if ~isempty(Tr)
    BX = spm_dctmtx(VG(1).dim(1),size(Tr,1),x-1);
    BY = spm_dctmtx(VG(1).dim(2),size(Tr,2),y-1);
    BZ = spm_dctmtx(VG(1).dim(3),size(Tr,3),z-1);
end
[X,Y] = ndgrid(x,y); clear x y
y1 = single(0); y1(VG(1).dim(1),VG(1).dim(2),VG(1).dim(3)) = 0;
y2 = single(0); y2(VG(1).dim(1),VG(1).dim(2),VG(1).dim(3)) = 0;
y3 = single(0); y3(VG(1).dim(1),VG(1).dim(2),VG(1).dim(3)) = 0;
M = VG(1).mat;
for j=1:length(z);
    if ~isempty(Tr)
        X1 = X    + BX*get_2Dtrans(Tr(:,:,:,1),BZ,j)*BY';
        Y1 = Y    + BX*get_2Dtrans(Tr(:,:,:,2),BZ,j)*BY';
        Z1 = z(j) + BX*get_2Dtrans(Tr(:,:,:,3),BZ,j)*BY';
    else
        X1 = X; Y1 = Y; Z1 = z(j);
    end
    y1(:,:,j) = single(M(1,1)*X1 + M(1,2)*Y1 + M(1,3)*Z1 + M(1,4));
    y2(:,:,j) = single(M(2,1)*X1 + M(2,2)*Y1 + M(2,3)*Z1 + M(2,4));
    y3(:,:,j) = single(M(3,1)*X1 + M(3,2)*Y1 + M(3,3)*Z1 + M(3,4));
end
clear X1 Y1 Z1 X Y z
disp(['Inverting the deformations field...']); 
M = Affine/VG(1).mat; M(4,:) = [0 0 0 1];
[iy1,iy2,iy3] = spm_invdef(y1,y2,y3,VF.dim(1:3),M,VG(1).mat); 
clear y1 y2 y3
M = inv(Vatlas.mat);
for j = 1:VF.dim(3)
    A = zeros(VF.dim(1),VF.dim(2));
    disp(['Slice ----> ' num2str(j)]);
    if sum(sum(Iseg_mask(:,:,j  )))~=0
        X2 = M(1,1)*double(iy1(:,:,j)) + M(1,2)*double(iy2(:,:,j)) + M(1,3)*double(iy3(:,:,j)) + M(1,4);
        Y2 = M(2,1)*double(iy1(:,:,j)) + M(2,2)*double(iy2(:,:,j)) + M(2,3)*double(iy3(:,:,j)) + M(2,4);
        Z2 = M(3,1)*double(iy1(:,:,j)) + M(3,2)*double(iy2(:,:,j)) + M(3,3)*double(iy3(:,:,j)) + M(3,4);
        A = spm_sample_vol(Iatlas,X2,Y2,Z2,0);
    end
    A(isnan(A)) = 0;
    spm_write_plane(Vol,squeeze(A.*Iseg_mask(:,:,j)),j);
end
fclose all;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T2 = get_2Dtrans(T3,B,j)
d   = [size(T3) 1 1 1];
tmp = reshape(T3,d(1)*d(2),d(3));
T2  = reshape(tmp*B(j,:)',d(1),d(2));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = Iso_Rem(T,Nhood)
%
%This function removes isolated points from white matter mask. 
%
% Input Parameters:
%   T            : White Matter  Mask
%   Nhood        : Minimun number of neighbors. 
% Output Parameters:
%   I            : White Matter Mask without isolated points  
%__________________________________________________________________________
% Authors:  Yasser Alem�n G�mez 
% Neuroimaging Department
% Cuban Neuroscience Center
% Last update: November 15th 2005
% Version $1.0

warning off
I = zeros(size(T)+2);
I(2:end-1,2:end-1,2:end-1) = T;
clear T
ind = find(I>0);
[x,y,z] = ind2sub(size(I), ind);
s = size(x,1);
sROI = zeros(size(I));
% figure;spy(I(:,:,160));
[X, Y, Z] = meshgrid(-1:1,-1:1,-1:1); 
X = X(:);Y = Y(:);Z = Z(:);
Neib = [X Y Z];clear X Y Z;
pos = find((Neib(:,1)==0)&(Neib(:,2)==0)&(Neib(:,3)==0));
Neib(pos,:) = [];
for i =1:26
M = Neib(i,:);
S = [x y z] + M(ones(s,1),:);
ind2 = sub2ind(size(I),S(:,1),S(:,2),S(:,3));
sROI(ind) = sROI(ind) + I(ind2);
end
ind = find(sROI<Nhood);
I(ind) =0;
I = I(2:end-1,2:end-1,2:end-1);
return;